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

/*
 Module for prediction or simulation:
 - loop over all prediction/simulation locations (regular or random path);
 - make selection at that location;
 - make an estimate (predict/simulate) at that location
   (if mask map has MV at the location, write MV)
 - write estimate to file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defs.h"
#include "mapio.h"
#include "userio.h"
#include "data.h"
#include "utils.h"
#include "debug.h"
#include "block.h"
#include "report.h"
#include "glvars.h"
#include "random.h"
#include "version.h"
#include "getest.h"
#include "msim.h"
#include "sim.h"
#include "select.h"
#include "predict.h"

typedef enum {
	AT_POINTS,
	AT_GRIDMAP
} PRED_AT;

static GRIDMAP **masks = NULL;
#ifdef WITH_SPIRAL
static DATA_GRIDMAP *mask_topology = NULL;
#endif

/* global variables: */
int strata_min;
unsigned int n_pred_locs = 0;
#define STRATUM(val) (floor(val - strata_min + 0.5))

void map_sign(GRIDMAP *m, const char *what);

unsigned int *get_n_sim_locs_table(unsigned int *size) {
	unsigned int i, *table;

	*size = (int) get_n_vars();
	table = (unsigned int *) emalloc(*size * sizeof(int));
	for (i = 0; i < *size; i++)
		table[i] = n_pred_locs;
	return table;
}

/* 
 * procedure map_sign() for putting history and description to maps (csf)
 * history is fixed for a session; description depends on the map contents
 */
void map_sign(GRIDMAP *m, const char *what) {
	char *pwd = NULL, *user = NULL, *timestring = NULL, *fname = NULL;
	time_t tm;

	if ((user = getenv("LOGNAME")) == NULL)
		user = ""; 
	if ((pwd = getenv("PWD")) == NULL)
		pwd = "";
	if ((fname = command_file_name) == NULL)
		fname = "";
	tm = time(&tm);
	if ((timestring = asctime(localtime(&tm))) == NULL)
		timestring = "";
	m->history = (char *) emalloc(1024 * sizeof(char));
#ifdef HAVE_SNPRINTF
	snprintf
#else 
	sprintf
#endif
		(m->history,
#ifdef HAVE_SNPRINTF
		1024,
#endif
		"%s%s\n%s %s %s\n%s%s\n%s %s\n%s %s\n%s %s (seed %lu)\n",
		*user ? 
		"creator:       " : "", user, 
		"program:      ", argv0, VERSION,
		*pwd ? 
		"path:          " : "", pwd,
		"command file: ", fname,
		"method used:  ", method_string(get_method()),
		"date:         ", timestring, get_seed());
	if (what == NULL) { /* paranoia */
		m->description = (char *) emalloc(sizeof(char));
		m->description[0] = '\0';
	} else {
		m->description = (char *) emalloc((strlen(what) + 1) * sizeof(char));
#ifdef HAVE_SNPRINTF
		snprintf(m->description, strlen(what) + 1, "%s\n", what);
#else
		sprintf(m->description, "%s\n", what);
#endif

	}
}

const void *get_mask0(void) {
	if (masks)
		return masks[0];
	return NULL;
}
