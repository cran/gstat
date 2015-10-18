/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2011 (C) Edzer Pebesma

    Edzer Pebesma, edzer.pebesma@uni-muenster.de
	Institute for Geoinformatics (ifgi), University of Münster 
	Weseler Straße 253, 48151 Münster, Germany. Phone: +49 251 
	8333081, Fax: +49 251 8339763  http://ifgi.uni-muenster.de 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version. As a special exception, linking 
    this program with the Qt library is permitted.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    (read also the files COPYING and Copyright)
*/

/*! \file mapio.c 
	\brief generic grid map I/O library functions 

mapio.c offers a generic grid (raster) map i/o interface to gstat.
It offers functions for reading and writing files, and access to cell
values (get/put). Raster map format is auto-detected.

Current implementation supports the following formats:
	arcinfo gridascii
	arcinfo gridfloat
	idrisi ascii
	idrisi binary real
	PCRaster/csf (all formats in, REAL4 out)
	ER-Mapper (all formats; 4-byte real output)
	GMT (netcdf)
	T2 (mike-she grid maps) (when #ifdef HAVE_T2_GRIDFORMAT)
	Surfer (DSAA/ascii)

TODO:
	? implement optional arcinfo/idrisi map read/writing on a per cel basis
	- integrate map topology information in local searches, esp. for cs.
	[[[ for simulation, do use a storage on a grid basis, not in a
	DATA structure: efficient because of grid ordering -> easy searching]]]

	GSLIB grid definition (is it useful? read as a DATA structure with a tag?):
	xmn ymn zmn # coordinates centre of first block
	nx ny nz    # number of blocks;
	xsiz ysiz zsiz # block size in x,y and z direction
	value1      # values following, x cycles fastest, then y, then z
	value2      # loc = (iz-1)nx.ny + (iy-1)nx + ix
	...	 # iz = 1 + int(loc/(nx.ny))
		    # iy = 1 + int((loc-(iz-1)nx.ny)/nx)
		    # ix = loc - (iz-1)nx.ny-(iy-1)nx
*/

/* 
 * if (gl_rowwise != 0), no complete map is not kept in memory.
 * take care that the complete map should be blanked (RputAllMV())
 * before writing takes place, when applying this method to other
 * formats (e.g. GDAL), as map_set_row() is only applied to WRITE_ONLY
 * maps for rows that contain non-missing valued cells!
 * SECOND: WRITE_ONLY maps should be open in read-write (r+) state!!
 * */

#include <stdio.h>
#include <math.h>				/* floor() */
#include <float.h>				/* FLT_MAX */
#include <string.h>				/* strtok() */
#include <ctype.h>				/* isalpha() */
#include <stdlib.h>				/* exit */

#include "defs.h"

#ifdef HAVE_LIBCSF
#include "csf.h"
#endif

#ifdef HAVE_LIBGIS
# include "gis.h"
#endif

/*! if defined, a stand-alone mapio library can be compiled */
#ifndef MAPIO_LIB
# include "glvars.h"
#else							/* standalone version: */
# include "defaults.h"
int debug_level = 1, gl_secure = 0, gl_rowwise = 1;
double gl_zero = DEF_zero;
char *gl_mv_string = "NA";
#endif

#include "utils.h"
#include "debug.h"
#include "read.h"
#include "userio.h"
#include "mapio.h"

static GRIDMAP *write_error(GRIDMAP * m);

#define SWAP_N(a,n) swap_floats((unsigned char *)a,n)
#define SWAP_M_N(a,m,n) swap_multiformat((unsigned char *)a,m,n)


#define CHECK_ROWS     1
#define CHECK_COLS     2
#define CHECK_CELLSIZE 4
#define CHECK_X_UL     8
#define CHECK_Y_UL    16
#define CHECK_SUM     31		/* sum of all checks */

#define BINARY_NATIVE      1
#define BINARY_NON_NATIVE  2
#define DEFAULT_MISVAL -9999.0
#define SURFER_MISVAL 1.70141E+38

/*
 * create a new GRIDMAP structure
 * allocates memory and initializes all fields for a GRIDMAP structure
 * returns: pointer to GRIDMAP structure
 */
GRIDMAP *new_map(MAP_READ_STATUS status)
{
	GRIDMAP *map;

	map = (GRIDMAP *) emalloc(sizeof(GRIDMAP));
	map->status = status;
	map->type = MT_UNKNOWN;
	map->history = NULL;
	map->description = NULL;
	map->filename = NULL;
	map->rows = 0;
	map->cols = 0;
	map->base_size = 0;
	map->grid = NULL;
	map->base = NULL;
	map->first_time_row = NULL;
	map->is_binary = 0;
	map->celltype = CT_UNKNOWN;
	map->misval = DEFAULT_MISVAL;	/* only for arcgrid */
	map->cellmin = map->cellmax = FLT_MAX;
	map->CSF_MAP = NULL;
#ifdef HAVE_LIBGDAL
	map->GeoTransform = (double *) emalloc(6 * sizeof(double));
#endif
	map->write = write_error;
	map->read_row = map->write_row = NULL;
	map->current_row = 0;
	return map;
}

static GRIDMAP *write_error(GRIDMAP * m)
{
	pr_warning("%s: writing this map format is not supported", m->filename);
	assert(0);
	return NULL;
}
/*
 * give x,y coordinate of cell center for cell [row, col]
 * libcsf has it's own function; other formats assume increasing x
 * for increasing cols and decreasing y for increasing rows
 * returns: non-zero if row or col are outside map limits
 */
int map_rowcol2xy(GRIDMAP * m,	/* pointer to gridmap */
				  unsigned int row,	/* current row number */
				  unsigned int col,	/* current column number */
				  double *x,	/* return value: pointer to x-coordinate */
				  double *y /* return value: pointer to y-coordinate */ )
{
	assert(m);
	assert(x);
	assert(y);

	if (row >= m->rows || col >= m->cols)
		return 1;
#if defined(HAVE_LIBCSF)
	if (m->CSF_MAP)
		RgetCoords((MAP *) m->CSF_MAP, 1, row, col, x, y);
	else
#endif
	{
		*x = m->x_ul + (col + 0.5) * m->cellsizex;
		*y = m->y_ul - (row + 0.5) * m->cellsizey;
	}
	return 0;
}
/*
 * converts x and y coordinate to (row,col) pair.
 *
 * see comment for map_rowcol2xy()
 *
 * returns: non-zero if x or y are outside map limits
 */
int map_xy2rowcol(GRIDMAP * m /* pointer to map */ ,
				  double x,		/* x-coordinate */
				  double y,		/* y-coordinate */
				  unsigned int *row,	/* output value: pointer to row number */
				  unsigned int *col
				  /* output value: pointer to column number */ )
{
	assert(m);
	assert(row);
	assert(col);

#if defined(HAVE_LIBCSF)
	if (m->CSF_MAP) {			/* handle possible non-zero map angle, CSF 2+: */
		if (RgetRowCol
			((MAP *) m->CSF_MAP, x, y, (size_t *) row,
			 (size_t *) col) != 1) return 1;
	} else
#endif
	{
		if (x < m->x_ul || x > m->x_ul + m->cols * m->cellsizex ||
			y > m->y_ul || y < m->y_ul - m->rows * m->cellsizey)
			return 1;
		*row = (unsigned int) floor((m->y_ul - y) / m->cellsizey);
		*col = (unsigned int) floor((x - m->x_ul) / m->cellsizex);
		if (*row == m->rows) /* on the bottom edge */
			*row = *row - 1;
		if (*col == m->cols) /* on the right edge */
			*col = *col - 1;
	}
	return 0;
}

