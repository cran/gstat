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
 * TODOs:
 * - deal with dummy data (unc.sim); set up quadtree blocks properly
 * - catch errors if no response is present in prediction
 */

#include <time.h> /* for s_gstat_progress function */

#include "config.h" /* may define USING_R */
# include <S.h> /* defines seed_in, also for R */

#ifdef USING_R
# include <R.h>
# include <Rinternals.h>
/* # include <R_ext/Utils.h> */
/* # include <Rinternals.h> */
# define R_UNIFORM unif_rand()
# define R_NORMAL  norm_rand()
# define RANDIN seed_in((long *) NULL)
# define RANDOUT seed_out((long *) NULL)
# define S_EVALUATOR
#else /* some S-Plus version; assuming >= 6 for now: */
# if (!defined(SPLUS_VERSION) || SPLUS_VERSION < 6000)
#  error("no SPLUS_VERSION >= 6.0")
# endif
# define SEXP s_object *
# define PROTECT(x) x
# define UNPROTECT(x)
# define R_UNIFORM unif_rand(S_evaluator)
# define R_NORMAL  norm_rand(S_evaluator)
# define RANDIN seed_in((long *) NULL, S_evaluator)
# define RANDOUT seed_out((long *) NULL, S_evaluator)
# define Rprintf printf
#endif

#include "defs.h"
#include "data.h"
#include "select.h"
#include "utils.h"
#include "userio.h"
#include "vario.h"
#include "fit.h"
#include "lex.h"
#include "sem.h"
#include "glvars.h"
#include "debug.h"
#include "mapio.h"
#include "msim.h"
#include "random.h"
#include "getest.h"
#include "polygon.h"
#include "s.h"

void s_gstat_printlog(const char *mess);
void s_gstat_progress(unsigned int current, unsigned int total);
double s_r_uniform(void);
double s_r_normal(void);
static int seed_is_in = 0;
static DATA_GRIDMAP *gstat_S_fillgrid(SEXP gridparams);
static void gstat_set_block(long i, SEXP block, SEXP block_cols, DPOINT *current);
static void S_no_progress(unsigned int current, unsigned int total);
static const char VarName[] = "(R Data)";

SEXP gstat_init(SEXP s_debug_level) {

	S_EVALUATOR

	remove_all();
	init_userio(1);  
	/* 1: set up for stdio */
	set_gstat_progress_handler(S_no_progress);
	set_gstat_error_handler(s_gstat_error);
	set_gstat_warning_handler(s_gstat_warning);
	set_gstat_log_handler(s_gstat_printlog);
	setup_meschach_error_handler(1);
	init_global_variables();
	init_data_minmax();
	RANDIN; /* load R/S seed into rng */
	seed_is_in = 1;
	set_rng_functions(s_r_uniform, s_r_normal, "S/R random number generator");
	debug_level = INTEGER(s_debug_level)[0];
	if (debug_level < 0) {
		debug_level = -debug_level;
		set_gstat_progress_handler(s_gstat_progress);
	}
	gl_secure = 1;
	return(s_debug_level);
}

SEXP gstat_exit(SEXP x) {

	S_EVALUATOR

	RANDOUT; /* write seed back to R/S engine */
	seed_is_in = 0;
	remove_all();
	/* init_global_variables(); */
	return(x);
}

SEXP gstat_new_data(SEXP sy, SEXP slocs, SEXP sX, SEXP has_intercept, 
			SEXP beta, SEXP nmax, SEXP nmin, SEXP maxdist, 
			SEXP vfn, SEXP sw, SEXP grid, SEXP degree, SEXP is_projected,
			SEXP vdist, SEXP lambda, SEXP omax) {
	double *y, *locs, *X, *w = NULL;
	long i, j, id, n, dim, n_X, has_int;
	DPOINT current;
	DATA **d;
	char name[20];

	S_EVALUATOR

	/* sy = AS_NUMERIC(sy); */
	sy = coerceVector(sy,REALSXP);
	n = LENGTH(sy);
	y = REAL(sy);
	if (n == 0)
		ErrMsg(ER_IMPOSVAL, "no data read");

	if (LENGTH(slocs) % n != 0)
		PROBLEM "dimensions do not match: locations %d and data %ld",
			LENGTH(slocs), n ERROR;
	dim = LENGTH(slocs) / n;
	if (dim <= 0)
		PROBLEM "too few spatial dimensions: %ld", dim ERROR;
	if (dim > 3)
		PROBLEM "too many spatial dimensions: %ld", dim ERROR;
	locs = REAL(slocs);

	if (LENGTH(sw) == n)
		w = REAL(sw);

	if (LENGTH(sX) % n != 0)
		PROBLEM "dimensions do not match: X %d and data %ld: missing values in data?",
			LENGTH(sX), n ERROR;
	n_X = LENGTH(sX) / n;
	X = REAL(sX);

	assert(n_X > 0);
	current.z = 0.0;
	current.bitfield = 0;

	id = get_n_vars();
	sprintf(name, "var%ld", id);
	which_identifier(name);
	d = get_gstat_data();
	d[id]->id = id;

	d[id]->n_list = d[id]->n_max = 0;
	d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
	d[id]->x_coord = "x";
	d[id]->y_coord = "y";
	d[id]->z_coord = "z";
	d[id]->variable = "R data";
	d[id]->fname = "R data";
	d[id]->lambda = REAL(lambda)[0];
	has_int = INTEGER(has_intercept)[0];
	/* increase d[id]->n_X and set d[id]->colX[i]: */
	for (i = d[id]->n_X = 0; i < n_X; i++) 
		data_add_X(d[id], i + (has_int ? 0 : 1)); 
	assert(d[id]->n_X == n_X);
	for (i = 0; i < LENGTH(beta); i++) /* do nothing if beta is numeric(0) */
		d[id]->beta = push_d_vector(REAL(beta)[i], d[id]->beta);
	if (INTEGER(nmax)[0] > 0) /* leave default (large) if < 0 */
		d[id]->sel_max = INTEGER(nmax)[0];
	if (INTEGER(omax)[0] > 0) /* leave default (0) if <= 0 */
		d[id]->oct_max = INTEGER(nmax)[0];
	if (INTEGER(nmin)[0] > 0) /* leave default (0) if <= 0 */
		d[id]->sel_min = INTEGER(nmin)[0];
	if (REAL(maxdist)[0] > 0.0)
		d[id]->sel_rad = REAL(maxdist)[0];
	switch(INTEGER(vfn)[0]) {
		case 1: /* d[id]->variance_fn = v_identity; == leave NULL */ break;
		case 2: d[id]->variance_fn = v_mu; break;
		case 3: d[id]->variance_fn = v_bin; break;
		case 4: d[id]->variance_fn = v_mu2; break;
		case 5: d[id]->variance_fn = v_mu3; break;
		default: PROBLEM "unknown variance function %d", 
				 	INTEGER(vfn)[0] ERROR;
	}
	gl_longlat = (INTEGER(is_projected)[0] == 0);
	d[id]->mode = X_BIT_SET | V_BIT_SET;
	if (dim > 1)
		d[id]->mode = d[id]->mode | Y_BIT_SET;
	if (dim > 2)
		d[id]->mode = d[id]->mode | Z_BIT_SET;
	set_norm_fns(d[id]); /* so gstat can calculate distances */
	if (w != NULL)
		d[id]->colnvariance = 1; /* it is set */
	switch(LENGTH(grid)) {
		case 0: case 1: break; /* empty, i.e., numeric(0) */
		case 6: d[id]->grid = gstat_S_fillgrid(grid); break;
		default: PROBLEM 
			"length of grid topology %d unrecognized", LENGTH(grid) ERROR;
	}
	d[id]->polynomial_degree = INTEGER(degree)[0];
	if (d[id]->polynomial_degree < 0 || d[id]->polynomial_degree > 3) {
		PROBLEM "polynomial degree should be 0, 1, 2 or 3" ERROR;
	}
	if (d[id]->polynomial_degree > 0) { 
		/* we're doing polynomials through degree: */
		if (id > 0) {
			PROBLEM "polynomial degree will only work for a single variable" ERROR;
		} if (n_X > 1) {
			PROBLEM "polynomial degree only works when no other predictors are given" ERROR;
		}
		setup_polynomial_X(d[id]); /* standardized coordinate polynomials */
	}
	d[id]->vdist = INTEGER(vdist)[0];
	assert(n_X <= d[id]->n_X);
	current.X = (double *) emalloc(d[id]->n_X * sizeof(double));

	SET_POINT(&current);
	current.u.stratum = 0;
	current.attr = current.x = current.y = current.z = 0.0;
	for (i = 0; i < n; i++) { /* loop over points */
		current.attr = y[i];
		current.x = locs[i];
		if (dim >= 2)
			current.y = locs[n + i];
		if (dim >= 3)
			current.z = locs[2 * n + i];
		/* track min/max coordinates, also for z, for the qtree bbox */
		if (i == 0) {
			d[id]->maxX = d[id]->minX = current.x;
			d[id]->maxY = d[id]->minY = current.y;
			d[id]->maxZ = d[id]->minZ = current.z;
		} else {
			d[id]->minX = MIN(d[id]->minX, current.x);
			d[id]->maxX = MAX(d[id]->maxX, current.x);
			d[id]->minY = MIN(d[id]->minY, current.y);
			d[id]->maxY = MAX(d[id]->maxY, current.y);
			d[id]->minZ = MIN(d[id]->minZ, current.z);
			d[id]->minZ = MIN(d[id]->minZ, current.z);
		}
		for (j = 0; j < n_X; j++)
			current.X[j] = X[j * n + i];
		if (w != NULL)
			current.variance = 1.0/(w[i]);
		push_point(d[id], &current);
	}
	check_global_variables();
	d[id]->n_original = d[id]->n_list;
	efree(current.X);
	return(sy);
}

SEXP gstat_new_dummy_data(SEXP loc_dim, SEXP has_intercept, SEXP beta, 
		SEXP nmax, SEXP nmin, SEXP maxdist, SEXP vfn, SEXP is_projected,
		SEXP vdist) {
	int i, id, dim, has_int;
	char name[20];
	DATA **d = NULL;

	S_EVALUATOR
	dim = INTEGER(loc_dim)[0];
	if (dim <= 0)
		PROBLEM "dimension value impossible: %d", dim ERROR;
	if (dim > 3)
		PROBLEM "too many dimensions: %d", dim ERROR;
	assert(LENGTH(beta) > 0);

	id = get_n_vars();
	sprintf(name, "var%d", id);
	which_identifier(name);
	d = get_gstat_data();
	d[id]->id = id;

	d[id]->n_list = d[id]->n_max = 0;
	d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
	d[id]->x_coord = "x";
	d[id]->y_coord = "y";
	d[id]->z_coord = "z";
	d[id]->variable = "R data";
	d[id]->fname = "R data";
	has_int = INTEGER(has_intercept)[0];
	for (i = d[id]->n_X = 0; i < LENGTH(beta); i++)
		data_add_X(d[id], i + (has_int ? 0 : 1));
	assert(d[id]->n_X == LENGTH(beta));
	d[id]->dummy = 1;
	for (i = 0; i < LENGTH(beta); i++)
		d[id]->beta = push_d_vector(REAL(beta)[i], d[id]->beta);
	if (INTEGER(nmax)[0] > 0) /* leave default (large) if < 0 */
		d[id]->sel_max = INTEGER(nmax)[0];
/* I doubt whether using nmin for dummy data _ever_ can have a
 * meaning, but hey, let's add it anyway. */
	if (INTEGER(nmin)[0] > 0) /* leave default (0) if <= 0 */
		d[id]->sel_min = INTEGER(nmin)[0];
	if (REAL(maxdist)[0] > 0.0)
		d[id]->sel_rad = REAL(maxdist)[0];
	switch(INTEGER(vfn)[0]) {
		case 1: /* d[id]->variance_fn = v_identity; -> leave NULL */ break;
		case 2: d[id]->variance_fn = v_mu; break;
		case 3: d[id]->variance_fn = v_bin; break;
		case 4: d[id]->variance_fn = v_mu2; break;
		case 5: d[id]->variance_fn = v_mu3; break;
		default: PROBLEM "unknown variance function %d", 
				 	INTEGER(vfn)[0] ERROR;
	}
	gl_longlat = (INTEGER(is_projected)[0] == 0);
	d[id]->vdist = INTEGER(vdist)[0];
	d[id]->mode = X_BIT_SET | V_BIT_SET;
	if (dim > 1)
		d[id]->mode = d[id]->mode | Y_BIT_SET;
	if (dim > 2)
		d[id]->mode = d[id]->mode | Z_BIT_SET;
	set_norm_fns(d[id]); /* so gstat can calculate distances */
	check_global_variables();
	d[id]->n_original = d[id]->n_list;
	return(loc_dim);
}

SEXP gstat_predict(SEXP sn, SEXP slocs, SEXP sX, SEXP block_cols, SEXP block, 
			SEXP weights, SEXP nsim, SEXP blue) {
	double *locs, **est_all, *X;
	long i, j, k, n, nvars, nest, dim, n_X, ncols_block, 
		nrows_block, pos;
	DPOINT current, *bp = NULL;
	DATA **d = NULL, *vd = NULL, *area = NULL;
	SEXP ret;
	SEXP retvector;
	SEXP retvector_dim;
	extern unsigned int n_pred_locs; /* predict.c, used in msim.c */
	float ***msim = NULL;

	S_EVALUATOR

	nvars = get_n_vars();
	nest = nvars + (nvars * (nvars + 1))/2;
	n = INTEGER(sn)[0];
	if (n <= 0 || LENGTH(slocs) == 0 || LENGTH(sX) == 0)
		ErrMsg(ER_IMPOSVAL, "newdata empty or only NA's");
	if (LENGTH(slocs) % n != 0)
		PROBLEM "dimensions do not match: locations %d, nrows in X %ld",
			LENGTH(slocs), n ERROR;
	dim = LENGTH(slocs) / n;
	if (dim > 3)
		PROBLEM "too many spatial dimensions: %ld", dim ERROR;
	if (dim <= 0)
		PROBLEM "too few spatial dimensions: %ld", dim ERROR;
	locs = REAL(slocs);
	if (LENGTH(sX) % n != 0)
		PROBLEM "dimensions do not match: X %d and data %ld",
			LENGTH(sX), n ERROR;
	n_X = LENGTH(sX) / n;

	current.attr = current.x = current.y = current.z = 0.0;
	current.bitfield = 0;
	/* assuming points ... */
	SET_POINT(&current);
	/* and then do the block thing: */
	if (LENGTH(block_cols) == 0) {
		bp = get_block_p();
		bp->x = bp->y = bp->z = 0.0; /* obsolete, I'd guess */
		if (LENGTH(block) >= 1) {
			bp->x = REAL(block)[0];
			SET_BLOCK(&current);
		}
		if (LENGTH(block) >= 2)
			bp->y = REAL(block)[1];
		if (LENGTH(block) >= 3)
			bp->z = REAL(block)[2];
		if (LENGTH(block) > 3)
			pr_warning("block dimension can only be 3; using the first 3");
	} else if (LENGTH(block_cols) == 1) { /* if > 1, block contains multiple 2D blocks */
		ncols_block = INTEGER(block_cols)[0];
		if (ncols_block < 1 || ncols_block > 3)
			ErrMsg(ER_IMPOSVAL, "block dimensions should be in [1..3]");
		nrows_block = LENGTH(block) / ncols_block; /* nr of rows */
		if (nrows_block > 0) {
			area = create_data_area();
			area->colnvariance = 0;
			area->n_list = area->n_max = 0;
			area->id = ID_OF_AREA;
			area->mode = X_BIT_SET;
			if (ncols_block > 1)
				area->mode = area->mode & Y_BIT_SET;
			if (ncols_block > 2)
				area->mode = area->mode & Z_BIT_SET;
			for (i = 0; i < nrows_block; i++) {
				current.x = REAL(block)[i];
				if (ncols_block > 1)
					current.y = REAL(block)[nrows_block + i];
				if (ncols_block > 2)
					current.z = REAL(block)[2 * nrows_block + i];
				if (LENGTH(weights) > 0) {
					area->colnvariance = 1;
					current.variance = REAL(weights)[i];
				}
				push_point(area, &current);
			}
			SET_BLOCK(&current);
		}
		if (DEBUG_FORCE)
			print_data_list(area);
	}

	X = REAL(sX);
	assert(n_X > 0);
	current.X = (double *) emalloc(n_X * sizeof(double));
	current.u.stratum = 0;
	d = get_gstat_data();
	est_all = (double **) emalloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		est_all[i] = (double *) emalloc(nest * sizeof(double));
	/* 
	 * the following is to fake gstat's default method handling: 
	 * we got to suggest that we'll go through a list of prediction
	 * locations, a la the gstat ``data(): ... ;'' command.
	 */
	vd = get_dataval(); 
	vd->id = ID_OF_VALDATA; 
	vd->mode = d[0]->mode;
	/* set min/max[XYZ] */
	vd->minY = vd->maxY = vd->minZ = vd->maxZ = 0.0;
	vd->minX = vd->maxX = locs[0];
	for (i = 1; i < n; i++) {
		vd->minX = MIN(vd->minX, locs[i]);
		vd->maxX = MAX(vd->maxX, locs[i]);
	}
	if (dim >= 2) {
		vd->minY = vd->maxY = locs[n];
		for (i = 1; i < n; i++) {
			vd->minY = MIN(vd->minY, locs[n + i]);
			vd->maxY = MAX(vd->maxY, locs[n + i]);
		}
	}
	if (dim >= 3) {
		vd->minZ = vd->maxZ = locs[2 * n];
		for (i = 1; i < n; i++) {
			vd->minZ = MIN(vd->minZ, locs[2 * n + i]);
			vd->maxZ = MAX(vd->maxZ, locs[2 * n + i]);
		}
	}

	/* fill, and standardize coordinate predictors from degree = x */
	for (i = 0; i < nvars; i++) 
		setup_data_minmax(d[i]);
	setup_data_minmax(vd);
	for (i = 0; i < nvars; i++) 
		calc_polynomials(d[i]);
	/* calc_polynomials(vd); */ /* still no data in fake vd */

	vd->polynomial_degree = d[0]->polynomial_degree;
	if (vd->polynomial_degree > 0) {
		setup_polynomial_X(vd); /* standardized coordinate polynomials */
		current.X = (double *) erealloc(current.X, vd->n_X * sizeof(double));
	}

	/* so far for the faking; now let's see what gstat makes out of this: */
	if (INTEGER(nsim)[0] == 0) {
		if (INTEGER(blue)[0] == 0) { /* FALSE */
			if (get_method() == NSP) /* choose default */
				set_method(get_default_method());
		} else 
			set_method(LSLM);
	} else {
		if (INTEGER(nsim)[0] < 0) {
			gl_nsim = -(INTEGER(nsim)[0]);
			set_method(ISI);
		} else {
			gl_nsim = INTEGER(nsim)[0];
			set_method(GSI);
		}
		n_pred_locs = n;
		if (gl_nsim > 1)
			init_simulations(d);
		if (get_n_beta_set() != get_n_vars())
			setup_beta(d, get_n_vars(), gl_nsim);
	}
	set_mode();  /* simple, stratified, multivariable? */
	check_global_variables(); /* it's there, better do it now */
	if (debug_level)
		Rprintf("[%s]\n", method_string(get_method()));
#ifdef USING_R
# ifdef WIN32
	R_FlushConsole();
	R_ProcessEvents();
# endif
#endif
	for (i = 0; i < n; i++) {
		print_progress(i, n);
		if (LENGTH(block_cols) > 1)
			gstat_set_block(i, block, block_cols, &current);
		current.x = locs[i];
		if (dim >= 2)
			current.y = locs[n + i];
		if (dim >= 3)
			current.z = locs[2 * n + i];
		for (j = 0; j < n_X; j++)
			current.X[j] = X[j * n + i];
		/* transform standardized coordinate polynomial here */
		if (vd->polynomial_degree)
			calc_polynomial_point(vd, &current);
		for (j = 0; j < get_n_vars(); j++)
			select_at(d[j], &current);
		get_est(d, get_method(), &current, est_all[i]);
#ifdef USING_R
# ifdef WIN32
		R_ProcessEvents(); /* avoid terminal freeze in R/Win */
# endif
		R_CheckUserInterrupt();
#endif
	}
	print_progress(100, 100);
	PROTECT(ret = allocVector(VECSXP, 1));
	PROTECT(retvector_dim = allocVector(REALSXP, 2));
	REAL(retvector_dim)[0] = n; /* nrows */
	if (gl_nsim > 1) {
		PROTECT(retvector = allocVector(REALSXP, gl_nsim * nvars * n));
		msim = get_msim();
		for (i = pos = 0; i < nvars; i++)
			for (j = 0; j < gl_nsim; j++) 
				for (k = 0; k < n; k++)
					REAL(retvector)[pos++] = msim[i][k][j];
		REAL(retvector_dim)[1] = nvars * gl_nsim; /* ncols */
	} else {
		PROTECT(retvector = allocVector(REALSXP, n * nest));
		for (j = pos = 0; j < nest; j++) {
			for (i = 0; i < n; i++) {
				if (is_mv_double(&(est_all[i][j])))
#ifdef USING_R /* avoid NaN's to be returned */
					REAL(retvector)[pos] = NA_REAL;
#else
					na_set(REAL(retvector) + pos, S_MODE_DOUBLE);
					/* the documentation says it should be DOUBLE */
#endif
				else
					REAL(retvector)[pos] = est_all[i][j];
				pos++;
			}
		}
		REAL(retvector_dim)[1] = nest; /* ncols */
	}
	if (gl_nsim > 0)
		free_simulations();
	/* SET_DIM(retvector, retvector_dim); */
	setAttrib(retvector, R_DimSymbol, retvector_dim);
	SET_VECTOR_ELT(ret, 0, retvector);
	for (i = 0; i < n; i++)
		efree(est_all[i]);
	efree(est_all);
	efree(current.X);
	UNPROTECT(3);
	return(ret);
}

static void gstat_set_block(long i, SEXP block, SEXP block_cols, DPOINT *current) {
	DATA *area;
	VARIOGRAM *v;
	long nrows_block, start, end, j;

	if (i >= LENGTH(block_cols) || i < 0)
		ErrMsg(ER_IMPOSVAL, "block_cols length less than nr of prediction locations");
	nrows_block = LENGTH(block) / 2; /* nr of rows */
	start = INTEGER(block_cols)[i];
	if (i == LENGTH(block_cols) - 1)
		end = nrows_block;
	else
		end = INTEGER(block_cols)[i+1] - 1;
	area = get_data_area();
	if (area != NULL)
		free_data(area);
	area = create_data_area();
	area->n_list = area->n_max = 0;
	area->id = ID_OF_AREA;
	area->mode = X_BIT_SET & Y_BIT_SET;
	for (j = start - 1; j < end; j++) {
		current->x = REAL(block)[j];
		current->y = REAL(block)[nrows_block + j];
		push_point(area, current);
	}
	SET_BLOCK(current);
	if (DEBUG_FORCE)
		print_data_list(area); 
	for (j = 0; j < get_n_vgms(); j++) {
		v = get_vgm(j);
		if (v != NULL)
			v->block_semivariance_set = v->block_covariance_set = 0; /* don't store under these circumstances! */
	}
	return;
}

SEXP gstat_variogram(SEXP s_ids, SEXP cutoff, SEXP width, SEXP direction, 
		SEXP cressie, SEXP dX, SEXP boundaries, SEXP grid, SEXP cov,
		SEXP pseudo) {
	SEXP ret;
	SEXP np; 
	SEXP dist;
	SEXP gamma;
	SEXP sx;
	SEXP sy;
	SEXP ev_parameters;
	/* SEXP y; */
	long i, id1, id2, nest;
	VARIOGRAM *vgm;
	DATA **d;

	GRIDMAP *m;
	unsigned int row, col, n;

	S_EVALUATOR

	id1 = INTEGER(s_ids)[0];
	if (LENGTH(s_ids) > 1)
		id2 = INTEGER(s_ids)[1];
	else
		id2 = id1;
	vgm = get_vgm(LTI(id1,id2));
	vgm->id = LTI(id1,id2);
	vgm->id1 = id1;
	vgm->id2 = id2;
	if (INTEGER(cov)[0] == 0)
		vgm->ev->evt = (id1 == id2 ? SEMIVARIOGRAM : CROSSVARIOGRAM);
	else if (INTEGER(cov)[0] == 1)
		vgm->ev->evt = (id1 == id2 ? COVARIOGRAM : CROSSCOVARIOGRAM);
	else {
		if (id1 != id2)
			ErrMsg(ER_IMPOSVAL,
			"cannot compute pairwise relative cross semivariogram");
		if (INTEGER(cov)[0] == 2)
			vgm->ev->evt = PRSEMIVARIOGRAM;
	}
	/* vgm->ev->is_asym = INTEGER(asym)[0]; */
	vgm->ev->pseudo = INTEGER(pseudo)[0];
	vgm->ev->recalc = 1;
	vgm->fname = NULL;
	if (LENGTH(cutoff) > 0)
		gl_cutoff = REAL(cutoff)[0];
	if (LENGTH(width) > 0)
		gl_iwidth = REAL(width)[0];
	gl_alpha = REAL(direction)[0];
	gl_beta = REAL(direction)[1];
	gl_tol_hor = REAL(direction)[2];
	gl_tol_ver = REAL(direction)[3];
	gl_cressie = INTEGER(cressie)[0];
	if (LENGTH(dX) > 0) {
		d = get_gstat_data();
		d[id1]->dX = REAL(dX)[0];
		d[id2]->dX = REAL(dX)[0];
		/* printf("dX1: %g ", d[id1]->dX);
		printf("dX2: %g\n", d[id2]->dX); */
	} 
	for (i = 0; i < LENGTH(boundaries); i++) /* does nothing if LENGTH is 0 */
		push_bound(REAL(boundaries)[i]);
	switch (LENGTH(grid)) {
		case 0: case 1: break;
		case 6: vgm->ev->S_grid = gstat_S_fillgrid(grid); break;
		default: PROBLEM "unrecognized grid length in gstat_variogram" ERROR;
			break;
	}

	calc_variogram(vgm, NULL);

	if (vgm->ev->S_grid != NULL) {
		PROTECT(ret = allocVector(VECSXP, 4));
		m = vgm->ev->map;
		n = m->rows * m->cols;
		PROTECT(np = allocVector(REALSXP, n));
		PROTECT(gamma = allocVector(REALSXP, n));
		PROTECT(sx = allocVector(REALSXP, n));
		PROTECT(sy = allocVector(REALSXP, n));

		for (row = i = 0; row < m->rows; row++) {
			for (col = 0; col < m->cols; col++) {
				map_rowcol2xy(m, row, col, &(REAL(sx)[i]), 
								&(REAL(sy)[i]));
				REAL(np)[i] = vgm->ev->nh[i];
				if (vgm->ev->nh[i] > 0)
					REAL(gamma)[i] = vgm->ev->gamma[i];
				else 
#ifdef USING_R /* avoid NaN's to be returned */
					REAL(gamma)[i] = NA_REAL;
#else
					na_set(REAL(gamma) + i, S_MODE_DOUBLE);
					/* the documentation says it should be DOUBLE */
#endif
				i++;
			}
		}
		SET_VECTOR_ELT(ret, 0, sx);
		SET_VECTOR_ELT(ret, 1, sy);
		SET_VECTOR_ELT(ret, 2, np);
		SET_VECTOR_ELT(ret, 3, gamma);
		UNPROTECT(5);
	} else {
		if (vgm->ev->cloud)
			nest = vgm->ev->n_est;
		else {
			/* Rprintf("[zero: %d]\n", vgm->ev->zero); */
			if (vgm->ev->zero == ZERO_SPECIAL)
				nest = vgm->ev->n_est;
			else 
				nest = vgm->ev->n_est - 1;
		}
		PROTECT(ret = allocVector(VECSXP, 4));
		if (nest <= 0)
			return(ret);
		PROTECT(np = allocVector(REALSXP, nest));
		PROTECT(dist = allocVector(REALSXP, nest));
		PROTECT(gamma = allocVector(REALSXP, nest));
		PROTECT(ev_parameters = allocVector(REALSXP, 4));
		REAL(ev_parameters)[0] = vgm->ev->cutoff;
		REAL(ev_parameters)[1] = vgm->ev->iwidth;
		REAL(ev_parameters)[2] = vgm->ev->pseudo;
		REAL(ev_parameters)[3] = vgm->ev->is_asym;
		for (i = 0; i < nest; i++) {
			REAL(np)[i] = vgm->ev->nh[i];
			REAL(dist)[i] = vgm->ev->dist[i];
			REAL(gamma)[i] = vgm->ev->gamma[i];
		}
		SET_VECTOR_ELT(ret, 0, np);
		SET_VECTOR_ELT(ret, 1, dist);
		SET_VECTOR_ELT(ret, 2, gamma);
		SET_VECTOR_ELT(ret, 3, ev_parameters);
		UNPROTECT(5);
	}
	return(ret);
}

SEXP gstat_load_variogram(SEXP s_ids, SEXP s_model, SEXP s_sills, SEXP s_ranges, 
		SEXP s_kappas, SEXP s_anis_all, SEXP s_table) 
{
	VARIOGRAM *vgm;
	long i, n, id1, id2, max_id;
	double anis[5] = {0.0, 0.0, 0.0, 1.0, 1.0}, rpars[2], *sills, *ranges, 
		*kappas, *anis_all;
	const char *model;

	sills = REAL(s_sills);
	ranges = REAL(s_ranges);
	kappas = REAL(s_kappas);
	anis_all = REAL(s_anis_all);

	id1 = INTEGER(s_ids)[0];
	id2 = INTEGER(s_ids)[1];
	max_id = MAX(id1, id2);

	if (get_n_vars() == 0)
		which_identifier("xx"); /* at least "load" one dummy var */
	if (max_id >= get_n_vars())
		ErrMsg(ER_IMPOSVAL,
			"gstat_load_variogram has been called with max_id >= n_vars");

	vgm = get_vgm(LTI(id1,id2));
	assert(vgm != NULL);

	vgm->id = LTI(id1,id2);
	vgm->id1 = id1;
	vgm->id2 = id2;
	vgm->n_models = vgm->n_fit = 0;

	n = LENGTH(s_sills);
	for (i = 0; i < n; i++) {
#ifdef USING_R
		model = CHAR(STRING_ELT(s_model, i));
#else
		model = STRING_POINTER(s_model)[i];
#endif
		anis[0] = anis_all[0 * n + i];
		anis[1] = anis_all[1 * n + i];
		anis[2] = anis_all[2 * n + i];
		anis[3] = anis_all[3 * n + i];
		anis[4] = anis_all[4 * n + i];
		rpars[0] = ranges[i];
		rpars[1] = kappas[i];
		if (LENGTH(s_table) > 0)
			push_to_v_table(vgm, rpars[0], 
					LENGTH(s_table), REAL(s_table),
					(anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis);
		else
			push_to_v(vgm, model, sills[i], rpars, 2,
				(anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis, 1, 1);
	}
	update_variogram(vgm);
	if (DEBUG_DUMP)
		logprint_variogram(vgm, 1); 
	return(s_model);
}

SEXP gstat_variogram_values(SEXP ids, SEXP pars, SEXP covariance, SEXP dist_values) {
	double from, to, n, d, x = 1.0, y = 0.0, z = 0.0;
	int i, id1, id2, cov = 0, ndist = 0;
	VARIOGRAM *vgm;
	SEXP dist;
	SEXP gamma;
	SEXP ret;

	S_EVALUATOR

	if (LENGTH(pars) != 3 && LENGTH(pars) != 6)
		PROBLEM "supply three or six distance parameters" ERROR;
	from = REAL(pars)[0];
	to = REAL(pars)[1];
	n = REAL(pars)[2];
	ndist = LENGTH(dist_values);
	cov = INTEGER(covariance)[0];
	if (LENGTH(pars) == 6) {
		x = REAL(pars)[3];
		y = REAL(pars)[4];
		z = REAL(pars)[5];
	}

	id1 = INTEGER(ids)[0];
	id2 = INTEGER(ids)[1];
	vgm = get_vgm(LTI(id1,id2));

	if (ndist > 0) {
		PROTECT(dist = allocVector(REALSXP, ndist));
		PROTECT(gamma = allocVector(REALSXP, ndist));
		for (i = 0; i < ndist; i++) {
			d = REAL(dist_values)[i];
			REAL(dist)[i] = d;
			REAL(gamma)[i] = (cov ? 
				get_covariance(vgm, d * x, d * y, d * z) : 
				get_semivariance(vgm, d * x, d * y, d * z));
		}
	} else {
		PROTECT(dist = allocVector(REALSXP, n));
		PROTECT(gamma = allocVector(REALSXP, n));
		for (i = 0; i < n; i++) {
			d = from;
			if (i > 0) /* implies n > 1 */
				d += (i/(n-1))*(to-from);
			REAL(dist)[i] = d;
			REAL(gamma)[i] = (cov ? 
				get_covariance(vgm, d * x, d * y, d * z) : 
				get_semivariance(vgm, d * x, d * y, d * z));
		}
	}
	PROTECT(ret = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0, dist);
	SET_VECTOR_ELT(ret, 1, gamma);
	UNPROTECT(3);
	return(ret);
}

// Added by Paul Hiemstra, 30-06-2008
SEXP get_covariance_list(SEXP ids, SEXP covariance, SEXP dist_list) {
	double d, x = 1.0, y = 0.0, z = 0.0;
	int i, id1, id2, cov = 0;
	VARIOGRAM *vgm;
	SEXP dist;
	SEXP gamma;
	SEXP ret;
	int length_list = LENGTH(dist_list);

	S_EVALUATOR

	cov = INTEGER(covariance)[0];

	id1 = INTEGER(ids)[0];
	id2 = INTEGER(ids)[1];
	vgm = get_vgm(LTI(id1,id2));

	PROTECT(dist = allocVector(REALSXP, length_list));
	PROTECT(gamma = allocVector(REALSXP, length_list));
	for (i = 0; i < length_list; i++) {
		d = REAL(dist_list)[i];
		REAL(dist)[i] = d;
		REAL(gamma)[i] = (cov ? 
			get_covariance(vgm, d * x, d * y, d * z) : 
			get_semivariance(vgm, d * x, d * y, d * z));
	}
	PROTECT(ret = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0, dist);
	SET_VECTOR_ELT(ret, 1, gamma);
	UNPROTECT(3);
	return(ret);
}

SEXP gstat_get_variogram_models(SEXP dolong) {
	SEXP ret;
	int i, n = 0, do_long;
	
	for (i = 1; v_models[i].model != NOT_SP; i++)
		n++;

	do_long = INTEGER(dolong)[0];
	PROTECT(ret = allocVector(STRSXP, n));
	for (i = 1; v_models[i].model != NOT_SP; i++)
#ifdef USING_R
		SET_STRING_ELT(ret, i-1, 
				mkChar(do_long ? v_models[i].name_long : v_models[i].name));
#else
		STRING_POINTER(ret)[i-1] = 
					string_dup(do_long ? v_models[i].name_long : v_models[i].name);
#endif
	UNPROTECT(1);
	return(ret);
}

SEXP gstat_load_command(SEXP commands) {
	int i;
	const char *cmd;
	SEXP error;

	PROTECT(error = allocVector(INTSXP, 1));
	INTEGER(error)[0] = 0;
	for (i = 0; i < LENGTH(commands); i++) {
#ifdef USING_R
		cmd = CHAR(STRING_ELT(commands, i));
#else
		cmd = CHARACTER_POINTER(commands)[i];
#endif
		if (parse_cmd(cmd, NULL)) {
			Rprintf("internal gstat string parse error on [%s]", cmd);
			INTEGER(error)[0] = i+1;
			UNPROTECT(1);
			return(error); 
		}
	}
	UNPROTECT(1);
	return(error);
}

void s_gstat_progress(unsigned int current, unsigned int total) {
	static int perc_last = -1, sec_last = -1;
	int perc, sec;
	static time_t start;

#ifdef USING_R
	R_CheckUserInterrupt(); /* allow for user interrupt */
#endif

	if (total <= 0 || DEBUG_SILENT)
		return;

	if (sec_last == -1) {
		start = time(NULL);
		sec_last = 0;
	}
	perc = floor(100.0 * current / total);
	if (perc != perc_last) { /* another percentage -> calculate time: */
		if (current == total) { /* 100% done, reset: */
			Rprintf("\r%3d%% done\n", 100);
			perc_last = sec_last = -1;
		} else {
			sec = difftime(time(NULL), start);
			if (sec != sec_last) { /* another second -- don't print too often */
				Rprintf("\r%3d%% done", perc);
				perc_last = perc;
				sec_last = sec;
			}
		}
	}
}

void s_gstat_error(const char *mess, int level) {
	/*	PROBLEM error_messages[level], mess ERROR; */
	if (mess == NULL)
		PROBLEM "<NULL> message" ERROR
	else
		PROBLEM "%s", mess ERROR
}

void s_gstat_warning(const char *mess) {

	print_to_logfile_if_open(mess);

	if (DEBUG_SILENT)
		return;
	Rprintf("%s\n", mess);
	return;
}

void s_gstat_printlog(const char *mess) {

	if (DEBUG_SILENT)
		return;

	print_to_logfile_if_open(mess);
	Rprintf("%s", mess);
}


SEXP gstat_load_ev(SEXP np, SEXP dist, SEXP gamma) {

	int i, cloud = 1;
	VARIOGRAM *vgm;

	S_EVALUATOR

	which_identifier("xx");
	/*
	 * vgm = get_vgm(LTI(INTEGER(id)[0], INTEGER(id)[1]));
	 * */
	vgm = get_vgm(LTI(0, 0));
	if (vgm->ev == NULL)
		vgm->ev = init_ev();
	vgm->ev->evt = SEMIVARIOGRAM;
	vgm->ev->n_est = LENGTH(np);
	vgm->ev->n_max = LENGTH(np);
	vgm->ev->gamma = (double *) emalloc (sizeof(double) * vgm->ev->n_max);
	vgm->ev->dist = (double *) emalloc (sizeof(double) * vgm->ev->n_max);
	vgm->ev->nh = (unsigned long *) emalloc (sizeof(long) * vgm->ev->n_max);
	for (i = 0; i < vgm->ev->n_est; i++) {
		vgm->ev->nh[i] = REAL(np)[i];
		vgm->ev->dist[i] = REAL(dist)[i];
		vgm->ev->gamma[i] = REAL(gamma)[i];
		if (cloud && vgm->ev->nh[i] > 1)
			cloud = 0;
	}
	vgm->ev->cloud = cloud;
	if (DEBUG_VGMFIT)
		fprint_sample_vgm(NULL, vgm->ev);
	return(np);
}

SEXP gstat_fit_variogram(SEXP fit, SEXP fit_sill, SEXP fit_range) {
	int i;
	VARIOGRAM *vgm;
	SEXP ret;
	SEXP sills;
	SEXP ranges;
	SEXP SSErr;
	SEXP fit_is_singular;

	vgm = get_vgm(LTI(0, 0));
	vgm->ev->fit = INTEGER(fit)[0];
	for (i = 0; i < vgm->n_models; i++) {
		vgm->part[i].fit_sill = INTEGER(fit_sill)[i];
		vgm->part[i].fit_range = INTEGER(fit_range)[i];
	}
	update_variogram(vgm);
	if (DEBUG_VGMFIT)
		logprint_variogram(vgm, 1);
	fit_variogram(vgm);
	if (DEBUG_VGMFIT)
		logprint_variogram(vgm, 1);

	PROTECT(sills = allocVector(REALSXP, vgm->n_models));
	PROTECT(ranges = allocVector(REALSXP, vgm->n_models));
	for (i = 0; i < vgm->n_models; i++) {
		REAL(sills)[i] = vgm->part[i].sill;
		REAL(ranges)[i] = vgm->part[i].range[0];
	}

	PROTECT(ret = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ret, 0, sills);
	SET_VECTOR_ELT(ret, 1, ranges);

	PROTECT(fit_is_singular = allocVector(REALSXP, 1));
	REAL(fit_is_singular)[0] = vgm->fit_is_singular;
	SET_VECTOR_ELT(ret, 2, fit_is_singular);

	PROTECT(SSErr = allocVector(REALSXP, 1));
	REAL(SSErr)[0] = vgm->SSErr;
	SET_VECTOR_ELT(ret, 3, SSErr);

	UNPROTECT(5);
	return(ret);
}

SEXP gstat_debug_level(SEXP level) {
	debug_level = INTEGER(level)[0];
	return(level);
}

double s_r_uniform(void) {
	double u;

	S_EVALUATOR
#ifdef USING_R
	if (!seed_is_in) PROBLEM "s_r_uniform(): seed is not read" ERROR; 
#else /* S-Plus needs RANDIN/RANDOUT to be inside the function calling a RNG */
	RANDIN;
#endif
	u = R_UNIFORM;
#ifndef USING_R
	RANDOUT;
#endif
	return(u);
}

double s_r_normal(void) {
	double r;

	S_EVALUATOR
#ifdef USING_R
	if (!seed_is_in) PROBLEM "s_r_normal(): seed is not read" ERROR; 
#else
	RANDIN;
#endif
	r = R_NORMAL;
#ifndef USING_R
	RANDOUT;
#endif
	return(r);
}

static DATA_GRIDMAP *gstat_S_fillgrid(SEXP gridparams) {
	double x_ul, y_ul, cellsizex, cellsizey;
	unsigned int rows, cols;

	cellsizex = REAL(gridparams)[2];
	cellsizey = REAL(gridparams)[3];
	rows = (unsigned int) REAL(gridparams)[5];
	cols = (unsigned int) REAL(gridparams)[4];
	x_ul = REAL(gridparams)[0] - 0.5 * cellsizex;
	y_ul = REAL(gridparams)[1] + (rows - 0.5) * cellsizey;
	/*
	printf("%g %g %g %g %u %u\n", x_ul, y_ul, cellsizex, cellsizey, rows, cols);
	fflush(stdout);
	*/
	return gsetup_gridmap(x_ul, y_ul, cellsizex, cellsizey, rows, cols);
}

static void S_no_progress(unsigned int current, unsigned int total) {
#ifdef USING_R
	R_CheckUserInterrupt();
#endif
}
