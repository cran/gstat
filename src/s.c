/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2003 (C) Edzer J. Pebesma

    Edzer J. Pebesma, e.pebesma@geog.uu.nl
    Department of physical geography, Utrecht University
    P.O. Box 80.115, 3508 TC Utrecht, The Netherlands

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

#include "config.h" /* may define USING_R */

#ifdef USING_R
# include <R.h>
# include <Rdefines.h>
/* # include <Rinternals.h> */
# define R_UNIFORM unif_rand()
# define R_NORMAL  norm_rand()
# define RANDIN seed_in((long *) NULL)
# define RANDOUT seed_out((long *) NULL)
# define S_EVALUATOR
#else /* some S-Plus version; assuming >= 6 for now: */
# include "S.h"
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

#include "data.h"
#include "utils.h"
#include "userio.h"
#include "vario.h"
#include "lex.h"
#include "sem.h"
#include "glvars.h"
#include "debug.h"
#include "mapio.h"
#include "msim.h"
#include "random.h"

void no_progress(unsigned int current, unsigned int total);
void s_gstat_error(const char *mess, int level);
void s_gstat_printlog(const char *mess);
double s_r_uniform(void);
double s_r_normal(void);
static int seed_is_in = 0;

SEXP gstat_init(SEXP s_debug_level) {

	S_EVALUATOR

	remove_all();
	init_userio(1);  
	/* 1: set up for stdio */
	set_gstat_progress_handler(no_progress);
	set_gstat_error_handler(s_gstat_error);
	set_gstat_log_handler(s_gstat_printlog);
	setup_meschach_error_handler();
	init_global_variables();
	RANDIN; /* load R/S seed into rng */
	seed_is_in = 1;
	set_rng_functions(s_r_uniform, s_r_normal, "S/R random number generator");
	debug_level = INTEGER_POINTER(s_debug_level)[0];
	return(s_debug_level);
}

SEXP gstat_exit(SEXP x) {

	S_EVALUATOR

	RANDOUT; /* write seed back to R/S engine */
	seed_is_in = 0;
	remove_all();
	init_global_variables();
	return(x);
}

SEXP gstat_new_data(SEXP sy, SEXP slocs, SEXP sX, SEXP has_intercept, 
			SEXP beta, SEXP nmax, SEXP vfn) {
	double *y, *locs, *X;
	long i, j, id, n, dim, n_X, has_int;
	DPOINT current;
	DATA **d;
	char name[20];

	S_EVALUATOR

	sy = AS_NUMERIC(sy);
	n = LENGTH(sy);
	y = NUMERIC_POINTER(sy);
	if (n == 0)
		ErrMsg(ER_IMPOSVAL, "no data read");

	if (LENGTH(slocs) % n == 0)
		dim = LENGTH(slocs) / n;
	else
		PROBLEM
			"dimensions do not match: locations %d and data %d",
			LENGTH(slocs), n
		ERROR;
	if (dim <= 0 || dim > 3)
		PROBLEM
			"too many dimensions: %d", dim
		ERROR;
	locs = NUMERIC_POINTER(slocs);

	if (LENGTH(sX) % n == 0)
		n_X = LENGTH(sX) / n;
	else
		PROBLEM
			"dimensions do not match: X %d and data %d",
			LENGTH(sX), n
		ERROR;
	X = NUMERIC_POINTER(sX);

	assert(n_X > 0);
	current.X = (double *) emalloc(n_X * sizeof(double));
	current.z = 0.0;
	current.bitfield = 0;

	id = get_n_vars();
	sprintf(name, "var%d", id);
	which_identifier(name);
	d = get_gstat_data();
	d[id]->id = id;

	d[id]->n_list = d[id]->n_max = 0;
	d[id]->colnx = d[id]->colny = d[id]->colnvalue = d[id]->colnz = 0;
	d[id]->x_coord = string_dup("x (S-plus)");
	d[id]->y_coord = string_dup("y (S-plus)");
	d[id]->z_coord = string_dup("z (S-plus)");
	d[id]->variable = string_dup("S-plus data");
	d[id]->fname = string_dup("S-plus data");
	has_int = INTEGER_POINTER(has_intercept)[0];
	/* increase d[id]->n_X and set d[id]->colX[i]: */
	for (i = d[id]->n_X = 0; i < n_X; i++) 
		data_add_X(d[id], i + (has_int ? 0 : 1)); 
	assert(d[id]->n_X == n_X);
	for (i = 0; i < LENGTH(beta); i++) /* do nothing if beta is numeric(0) */
		d[id]->beta = push_to_vector(NUMERIC_POINTER(beta)[i], d[id]->beta);
	if (INTEGER_POINTER(nmax)[0] > 0) /* leave default (large) if < 0 */
		d[id]->sel_max = INTEGER_POINTER(nmax)[0];
	switch(INTEGER_POINTER(vfn)[0]) {
		case 1: /* d[id]->variance_fn = v_identity; */ break;
		case 2: d[id]->variance_fn = v_mu; break;
		case 3: d[id]->variance_fn = v_bin; break;
		default: PROBLEM "unknown variance function %d", 
				 	INTEGER_POINTER(vfn)[0] ERROR;
	}
	d[id]->mode = X_BIT_SET | V_BIT_SET;
	if (dim > 1)
		d[id]->mode = d[id]->mode | Y_BIT_SET;
	if (dim > 2)
		d[id]->mode = d[id]->mode | Z_BIT_SET;
	set_norm_fns(d[id]); /* so gstat can calculate distances */
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
		push_point(d[id], &current);
	}
	check_global_variables();
	d[id]->n_original = d[id]->n_list;
	return(sy);
}

SEXP gstat_new_dummy_data(SEXP loc_dim, SEXP has_intercept, SEXP beta, 
		SEXP nmax, SEXP vfn) {
	int i, id, dim, has_int;
	char name[20];
	DATA **d = NULL;

	S_EVALUATOR
	dim = INTEGER_POINTER(loc_dim)[0];
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
	d[id]->x_coord = string_dup("x (S-plus)");
	d[id]->y_coord = string_dup("y (S-plus)");
	d[id]->z_coord = string_dup("z (S-plus)");
	d[id]->variable = string_dup("S-plus data");
	d[id]->fname = string_dup("S-plus data");
	has_int = INTEGER_POINTER(has_intercept)[0];
	for (i = d[id]->n_X = 0; i < LENGTH(beta); i++)
		data_add_X(d[id], i + (has_int ? 0 : 1));
	assert(d[id]->n_X == LENGTH(beta));
	d[id]->dummy = 1;
	for (i = 0; i < LENGTH(beta); i++)
		d[id]->beta = push_to_vector(NUMERIC_POINTER(beta)[i], d[id]->beta);
	if (INTEGER_POINTER(nmax)[0] > 0) /* leave default (large) if < 0 */
		d[id]->sel_max = INTEGER_POINTER(nmax)[0];
	switch(INTEGER_POINTER(vfn)[0]) {
		case 1: /* d[id]->variance_fn = v_identity; */ break;
		case 2: d[id]->variance_fn = v_mu; break;
		case 3: d[id]->variance_fn = v_bin; break;
		default: PROBLEM "unknown variance function %d", 
				 	INTEGER_POINTER(vfn)[0] ERROR;
	}
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
			SEXP nsim) {
	double *locs, **est_all, *X;
	long i, j, k, n, nvars, nest, dim, n_X, ncols_block, nrows_block, pos;
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
	n = INTEGER_POINTER(sn)[0];
	if (n <= 0 || LENGTH(slocs) == 0 || LENGTH(sX) == 0)
		ErrMsg(ER_IMPOSVAL, "empty newdata");
	if (LENGTH(slocs) % n == 0)
		dim = LENGTH(slocs) / n;
	else
		PROBLEM
			"dimensions do not match: locations %d and data %d",
			LENGTH(slocs), n
		ERROR;
	if (dim <= 0 || dim > 3)
		PROBLEM
			"too many dimensions: %d", dim
		ERROR;
	locs = NUMERIC_POINTER(slocs);
	if (LENGTH(sX) % n == 0)
		n_X = LENGTH(sX) / n;
	else
		PROBLEM
			"dimensions do not match: X %d and data %d",
			LENGTH(sX), n
		ERROR;

	current.attr = current.x = current.y = current.z = 0.0;
	current.bitfield = 0;
	/* assuming points ... */
	SET_POINT(&current);
	/* and then do the block thing: */
	if (LENGTH(block_cols) == 0) {
		bp = get_block_p();
		bp->x = bp->y = bp->z = 0.0; /* obsolete, I'd guess */
		if (LENGTH(block) >= 1) {
			bp->x = NUMERIC_POINTER(block)[0];
			SET_BLOCK(&current);
		}
		if (LENGTH(block) >= 2)
			bp->y = NUMERIC_POINTER(block)[1];
		if (LENGTH(block) >= 3)
			bp->z = NUMERIC_POINTER(block)[2];
		if (LENGTH(block) > 3)
			pr_warning("block dimension can only be 3; using the first 3");
	} else {
		ncols_block = INTEGER_POINTER(block_cols)[0];
		if (ncols_block < 1 || ncols_block > 3)
			ErrMsg(ER_IMPOSVAL, "block dimensions should be in [1..3]");
		nrows_block = LENGTH(block) / ncols_block; /* nr of rows */
		if (nrows_block > 0) {
			area = create_data_area();
			area->n_list = area->n_max = 0;
			area->id = ID_OF_AREA;
			area->mode = X_BIT_SET;
			if (ncols_block > 1)
				area->mode = area->mode & Y_BIT_SET;
			if (ncols_block > 2)
				area->mode = area->mode & Z_BIT_SET;
			for (i = 0; i < nrows_block; i++) {
				current.x = NUMERIC_POINTER(block)[i];
				if (ncols_block > 1)
					current.y = NUMERIC_POINTER(block)[nrows_block + i];
				if (ncols_block > 2)
					current.z = NUMERIC_POINTER(block)[2 * nrows_block + i];
				push_point(area, &current);
			}
			SET_BLOCK(&current);
		}
		if (DEBUG_FORCE)
			print_data_list(area);
	}

	X = NUMERIC_POINTER(sX);
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
	/* so far for the faking; now let's see what gstat makes out of this: */
	if (INTEGER_POINTER(nsim)[0] == 0)
		set_method(get_default_method());
	else {
		if (INTEGER_POINTER(nsim)[0] < 0) {
			gl_nsim = -(INTEGER_POINTER(nsim)[0]);
			set_method(ISI);
		} else {
			gl_nsim = INTEGER_POINTER(nsim)[0];
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
		current.x = locs[i];
		if (dim >= 2)
			current.y = locs[n + i];
		if (dim >= 3)
			current.z = locs[2 * n + i];
		for (j = 0; j < n_X; j++)
			current.X[j] = X[j * n + i];
		for (j = 0; j < get_n_vars(); j++)
			select_at(d[j], &current);
		get_est(d, get_method(), &current, est_all[i]);
#ifdef USING_R
# ifdef WIN32
		R_ProcessEvents(); /* avoid terminal freeze */
# endif
#endif
	}
	PROTECT(ret = NEW_LIST(1));
	PROTECT(retvector_dim = NEW_NUMERIC(2));
	NUMERIC_POINTER(retvector_dim)[0] = n; /* nrows */
	if (gl_nsim > 1) {
		PROTECT(retvector = NEW_NUMERIC(gl_nsim * nvars * n));
		msim = get_msim();
		for (i = pos = 0; i < nvars; i++)
			for (j = 0; j < gl_nsim; j++) 
				for (k = 0; k < n; k++)
					NUMERIC_POINTER(retvector)[pos++] = msim[i][k][j];
		free_simulations();
		NUMERIC_POINTER(retvector_dim)[1] = nvars * gl_nsim; /* ncols */
	} else {
		PROTECT(retvector = NEW_NUMERIC(n * nest));
		for (j = pos = 0; j < nest; j++) {
			for (i = 0; i < n; i++) {
				if (is_mv_double(&(est_all[i][j])))
#ifdef USING_R /* avoid NaN's to be returned */
					NUMERIC_POINTER(retvector)[pos] = NA_REAL;
#else
					na_set(NUMERIC_POINTER(retvector) + pos, S_MODE_DOUBLE);
					/* the documentation says it should be DOUBLE */
#endif
				else
					NUMERIC_POINTER(retvector)[pos] = est_all[i][j];
				pos++;
			}
		}
		NUMERIC_POINTER(retvector_dim)[1] = nest; /* ncols */
	}
	SET_DIM(retvector, retvector_dim);
	SET_ELEMENT(ret, 0, retvector);
	for (i = 0; i < n; i++)
		efree(est_all[i]);
	efree(est_all);
	efree(current.X);
	UNPROTECT(3);
	return(ret);
}

SEXP gstat_variogram(SEXP s_ids, SEXP cutoff, SEXP width, SEXP direction, 
		SEXP cressie, SEXP dX, SEXP boundaries) {
	SEXP ret;
	SEXP np; 
	SEXP dist;
	SEXP gamma;
	long i, id1, id2;
	VARIOGRAM *vgm;
	DATA **d;

	S_EVALUATOR

	id1 = INTEGER_POINTER(s_ids)[0];
	if (LENGTH(s_ids) > 1)
		id2 = INTEGER_POINTER(s_ids)[1];
	else
		id2 = id1;
	vgm = get_vgm(LTI(id1,id2));
	vgm->id = LTI(id1,id2);
	vgm->id1 = id1;
	vgm->id2 = id2;
	vgm->ev->evt = (id1 == id2 ? SEMIVARIOGRAM : CROSSVARIOGRAM);
	vgm->ev->recalc = 1;
	vgm->fname = NULL;
	if (LENGTH(cutoff) > 0)
		gl_cutoff = NUMERIC_POINTER(cutoff)[0];
	if (LENGTH(width) > 0)
		gl_iwidth = NUMERIC_POINTER(width)[0];
	gl_alpha = NUMERIC_POINTER(direction)[0];
	gl_beta = NUMERIC_POINTER(direction)[1];
	gl_tol_hor = NUMERIC_POINTER(direction)[2];
	gl_tol_ver = NUMERIC_POINTER(direction)[3];
	gl_cressie = INTEGER_POINTER(cressie)[0];
	if (LENGTH(dX) > 0) {
		d = get_gstat_data();
		d[id1]->dX = NUMERIC_POINTER(dX)[0];
	} 
	for (i = 0; i < LENGTH(boundaries); i++) /* do nothing if LENGTH is 0 */
		push_bound(NUMERIC_POINTER(boundaries)[i]);

	calc_variogram(vgm, NULL);

	ret = NEW_LIST(3);
	if (vgm->ev->n_est <= 1)
		return(ret);
	np = NEW_NUMERIC(vgm->ev->n_est - 1);
	dist = NEW_NUMERIC(vgm->ev->n_est - 1);
	gamma = NEW_NUMERIC(vgm->ev->n_est - 1);
	for (i = 0; i < vgm->ev->n_est - 1; i++) {
		NUMERIC_POINTER(np)[i] = vgm->ev->nh[i];
		NUMERIC_POINTER(dist)[i] = vgm->ev->dist[i];
		NUMERIC_POINTER(gamma)[i] = vgm->ev->gamma[i];
	}
	SET_ELEMENT(ret, 0, np);
	SET_ELEMENT(ret, 1, dist);
	SET_ELEMENT(ret, 2, gamma);
	return(ret);
}

void Cgstat_load_variogram(int *ids, int *n_models, 
		char **model, double *sills, double *ranges, double *kappas,
		double *anis_all) 
{
	char *vgm_model;
	VARIOGRAM *vgm;
	int i, n, id1, id2, max_id;
	double anis[5] = {0.0, 0.0, 0.0, 1.0, 1.0}, rpars[2];

	id1 = ids[0];
	id2 = ids[1];
	max_id = MAX(id1, id2);

	if (get_n_vars() == 0)
		which_identifier("xx"); /* at least "load" one dummy var */
	if (max_id >= get_n_vars())
		ErrMsg(ER_IMPOSVAL,
			"Cgstat_load_variogram has been called with max_id >= n_vars");

	vgm = get_vgm(LTI(id1,id2));
	assert(vgm != NULL);

	vgm->id = LTI(id1,id2);
	vgm->id1 = id1;
	vgm->id2 = id2;
	vgm->n_models = vgm->n_fit = 0;

	n = *n_models;
	for (i = 0; i < n; i++) {
		anis[0] = anis_all[0 * n + i];
		anis[1] = anis_all[1 * n + i];
		anis[2] = anis_all[2 * n + i];
		anis[3] = anis_all[3 * n + i];
		anis[4] = anis_all[4 * n + i];
		rpars[0] = ranges[i];
		rpars[1] = kappas[i];
		push_to_v(vgm, model[i], sills[i], rpars, 2,
			(anis[3] == 1.0 && anis[4] == 1.0) ? NULL : anis, 1, 1);
	}
	update_variogram(vgm);
	if (DEBUG_DUMP)
		fprint_variogram(stdout, vgm, 1); 
}

SEXP gstat_variogram_values(SEXP ids, SEXP pars)
{
	double from, to, n, d, x = 1.0, y = 0.0, z = 0.0;
	int i, id1, id2;
	VARIOGRAM *vgm;
	SEXP dist;
	SEXP gamma;
	SEXP ret;

	S_EVALUATOR

	if (LENGTH(pars) != 3 && LENGTH(pars) != 6)
		PROBLEM "supply three or six distance parameters" ERROR;
	from = NUMERIC_POINTER(pars)[0];
	to = NUMERIC_POINTER(pars)[1];
	n = NUMERIC_POINTER(pars)[2];
	if (LENGTH(pars) == 6) {
		x = NUMERIC_POINTER(pars)[3];
		y = NUMERIC_POINTER(pars)[4];
		z = NUMERIC_POINTER(pars)[5];
	}

	id1 = INTEGER_POINTER(ids)[0];
	id2 = INTEGER_POINTER(ids)[1];
	vgm = get_vgm(LTI(id1,id2));

	dist = NEW_NUMERIC(n);
	gamma = NEW_NUMERIC(n);
	for (i = 0; i < n; i++) {
		d = from + (i/(n-1))*(to-from);
		NUMERIC_POINTER(dist)[i] = d;
		NUMERIC_POINTER(gamma)[i] = get_semivariance(vgm, d * x, d * y, d * z);
	}
	ret = NEW_LIST(2);
	SET_ELEMENT(ret, 0, dist);
	SET_ELEMENT(ret, 1, gamma);
	return(ret);
}

SEXP gstat_get_n_variogram_models(SEXP x) {
	SEXP n;

	S_EVALUATOR

	n = NEW_INTEGER(1);
	INTEGER_POINTER(n)[0] = get_n_variogram_models();
	return(n);
}

void Cgstat_get_variogram_models(char **names) {
	int i, n;

	for (i = 1; v_models[i].model != NOT_SP; i++)
		names[i-1] = string_dup(v_models[i].name);
}

void Cload_gstat_command(char **commands, int *n, int *error) {
	int i;

	*error = 0;
	for (i = 0; i < *n; i++) {
		if (parse_cmd(commands[i], NULL)) {
			*error = i+1;
			return;
		}
	}
	return;
}

void no_progress(unsigned int current, unsigned int total) {
}

void s_gstat_error(const char *mess, int level) {
		PROBLEM error_messages[level], mess ERROR;
}

void s_gstat_warning(const char *mess) {

	print_to_logfile_if_open(mess);

	Rprintf("gstat warning: %s\n", mess);
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
	 * vgm = get_vgm(LTI(INTEGER_POINTER(id)[0], INTEGER_POINTER(id)[1]));
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
		vgm->ev->nh[i] = NUMERIC_POINTER(np)[i];
		vgm->ev->dist[i] = NUMERIC_POINTER(dist)[i];
		vgm->ev->gamma[i] = NUMERIC_POINTER(gamma)[i];
		if (cloud && vgm->ev->nh[i] > 1)
			cloud = 0;
	}
	vgm->ev->cloud = cloud;
	if (DEBUG_VGMFIT)
		fprint_sample_vgm(stdout, vgm->ev);
	return(np);
}

SEXP gstat_fit_variogram(SEXP fit, SEXP fit_sill, SEXP fit_range) {
	int i;
	VARIOGRAM *vgm;
	SEXP ret;
	SEXP sills;
	SEXP ranges;
	SEXP SSErr;

	vgm = get_vgm(LTI(0, 0));
	vgm->ev->fit = INTEGER_POINTER(fit)[0];
	for (i = 0; i < vgm->n_models; i++) {
		vgm->part[i].fit_sill = INTEGER_POINTER(fit_sill)[i];
		vgm->part[i].fit_range = INTEGER_POINTER(fit_range)[i];
	}
	update_variogram(vgm);
	if (DEBUG_VGMFIT)
		fprint_variogram(stdout, vgm, 1);
	fit_variogram(vgm);
	if (DEBUG_VGMFIT)
		fprint_variogram(stdout, vgm, 1);
	ret = NEW_LIST(3);
	sills = NEW_NUMERIC(vgm->n_models);
	ranges = NEW_NUMERIC(vgm->n_models);
	SSErr = NEW_NUMERIC(1);
	for (i = 0; i < vgm->n_models; i++) {
		NUMERIC_POINTER(sills)[i] = vgm->part[i].sill;
		NUMERIC_POINTER(ranges)[i] = vgm->part[i].range[0];
	}
	NUMERIC_POINTER(SSErr)[0] = vgm->SSErr;
	SET_ELEMENT(ret, 0, sills);
	SET_ELEMENT(ret, 1, ranges);
	SET_ELEMENT(ret, 2, SSErr);
	return(ret);
}

SEXP gstat_debug_level(SEXP level) {
	debug_level = INTEGER_POINTER(level)[0];
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

