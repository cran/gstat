#include "R_ext/Lapack.h"
#include "debug.h"
#include "mtrx.h"

static void mtrx_chol0(int *N, double *A, int *info);
static void mtrx_backsolve0(int trans, int *k, int *ncb, double *rr, int *nrr, double *b);

MAT *CHfactor(MAT *m, int *info) {
    if (m->m != m->n) 
		error("CHfactor: 'm' must be a square matrix");
	m = m_transp(m, m); /* in situ */
	for (int i = 1; i < m->m; i++)
		for (int j = 0; j < i; j++)
			m->me[j][i] = 0.0; /* zero lower triangle of transposed Fortran order, = upper tri */
	mtrx_chol0((int *)&(m->n), m->v, info);
	return(m);
}

static void mtrx_chol0(int *N, double *A, int *info) {

	F77_CALL(dpotrf)("Upper", N, A, N, info);
	if (*info != 0) {
	    if (*info > 0 && DEBUG_COV)
			warning("the leading minor of order %d is not positive definite", *info);
	    if (*info < 0)
			error("argument %d of Lapack routine %s had invalid value", -(*info), "dpotrf");
	}
}

static void mtrx_backsolve0(int trans, /* transpose? */
		int *k,  /* size of mat */
		int *ncb,  /* nr of columns in b (1) */
		double *rr,  /* pointer to mat */
		int *nrr, /* nr of rows rr */
		double *b /* rhs, answer */ ) {

    size_t incr = *nrr + 1;
	for (int i = 0; i < *k; i++) /* check for zeros on diagonal */
		if (rr[i * incr] == 0.0)
			error("singular matrix in 'backsolve'. First zero in diagonal [%d]", i + 1);
	double one = 1.0;
	F77_CALL(dtrsm)("L", "U", trans ? "T" : "N", "N", k, ncb, &one, rr, nrr, b, k);
}

VEC *CHsolve(MAT *m, VEC *b, VEC *out) {
	if (m->m != m->n) 
		error("CHsolve: 'm' must be a square matrix");
	if (m->m != b->dim) 
		error("CHsolve: b dim does not match m");
    /* copy b to out */
	out = v_copy(b, v_resize(out, b->dim));
	/* mtrx_backsolve(A, mtrx_backsolve(A, b, ncols(A), 1), ncols(A), 0); */
	/* first */
	int one = 1;
	mtrx_backsolve0(1, (int *)&(m->m), &one, m->v, (int *)&(m->m), out->ve);
	/* second */
	mtrx_backsolve0(0, (int *)&(m->m), &one, m->v, (int *)&(m->m), out->ve);
	return(out);
}
