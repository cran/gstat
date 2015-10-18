/* interface roughly follows meschach; implementation rewritten from scratch */
#include <string.h> /* memset */
#include "defs.h" /* CDECL */
#include "utils.h" /* efree, emalloc */
#include "userio.h" /* ErrMsg */
#include "mtrx.h"

void m_free(MAT *m) {
	efree(m->me);
	efree(m->v);
	efree(m);
}

void v_free(VEC *v) {
	efree(v->ve);
	efree(v);
}

void iv_free(IVEC *iv) {
	efree(iv->ive);
	efree(iv);
}

MAT *m_init(void) {
	MAT *m = emalloc(sizeof(MAT));
	m->n = m->m = m->max_n = m->max_m = 0;
	m->me = (double **) NULL;
	m->v = (double *) NULL;
	return(m);
}

MAT *m_resize(MAT *m, size_t new_r, size_t new_c) {
	if (m == MNULL)
		m = m_init();
	if (new_r * new_c > m->max_n * m->max_m)
		m->v = (double *) erealloc(m->v, new_r * new_c * sizeof(double)); /* takes care of NULL m */
	if (new_r > m->max_m) {
		m->me = (double **) erealloc(m->me, new_r * sizeof(double *)); /* takes care of NULL m */
		m->max_m = new_r;
	}
	if (new_c > m->max_n)
		m->max_n = new_c;
	for (size_t i = 0; i < new_r; i++)
		m->me[i] = m->v + (i * new_c);
	m->m = new_r;
	m->n = new_c;
	return(m);
}

VEC *v_init(void) {
	VEC *v = emalloc(sizeof(VEC));
	v->dim = v->max_n = 0;
	v->ve = NULL;
	return(v);
}

VEC *v_resize(VEC *v, size_t new_n) {
	if (v == NULL)
		v = v_init();
	if (new_n > v->max_n) {
		v->ve = erealloc(v->ve, new_n * sizeof(double));
		v->max_n = new_n;
	}
	v->dim = new_n;
	return(v);
}

PERM *p_init(void) {
	PERM *p = emalloc(sizeof(PERM));
	p->size = p->max_size = 0;
	p->pe = (size_t *) NULL;
	return(p);
}

PERM *px_resize(PERM *p, size_t new_n) {
	if (p == PNULL)
		p = p_init();
	if (new_n > p->max_size) {
		p->pe = erealloc(p->pe, new_n * sizeof(size_t));
		p->max_size = new_n;
	}
	p->size = new_n;
	return(p);
}

IVEC *iv_init(void) {
	IVEC *iv = emalloc(sizeof(IVEC));
	iv->size = iv->max_size = 0;
	iv->ive = (int *) NULL;
	return(iv);
}

IVEC *iv_resize(IVEC *iv, size_t new_n) {
	if (iv == IVNULL)
		iv = iv_init();
	if (new_n > iv->max_size) {
		iv->ive = erealloc(iv->ive, new_n * sizeof(int));
		iv->max_size = new_n;
	}
	iv->size = new_n;
	return(iv);
}

MAT *m_zero(MAT *m) {
	if (m != MNULL)
		for (int i = 0; i < m->m; i++)
			for (int j = 0; j < m->n; j++)
				m->me[i][j] = 0.0;
	return(m);
}

VEC *v_zero(VEC *v) {
	if (v != VNULL)
		for (int i = 0; i < v->dim; i++)
			v->ve[i] = 0.0;
	return(v);
}

void set_col(MAT *M, size_t i, VEC *col) {
	if (i > M->n || col->dim != M->m)
		ErrMsg(ER_IMPOSVAL, "dimension mismatch in set_col");
	for (size_t j = 0; j < col->dim; j++)
		M->me[j][i] = col->ve[j];
}

VEC *get_col(MAT *M, size_t i, VEC *to) {
	if (i > M->n || to->dim != M->m)
		ErrMsg(ER_IMPOSVAL, "dimension mismatch in get_col");
	for (size_t j = 0; j < to->dim; j++)
		to->ve[j] = M->me[j][i];
	return(to);
}

MAT *m_inverse(MAT *in, MAT *out, int *info) { /* out <- in^{-1} */
	out = m_resize(out, in->m, in->n);
	MAT *m_cp = m_copy(in, MNULL);
	m_cp = CHfactor(m_cp, info);
	if (*info != 0) /* singular */
		return(in);
	VEC *v1 = v_zero(v_resize(VNULL, in->m));
	VEC *v2 = v_zero(v_resize(VNULL, in->m));
	for (int i = 0; i < out->m; i++) {
		v_zero(v1);
		v1->ve[i] = 1.0; /* diag */
		if (v1->dim != m_cp->m)
			ErrMsg(ER_IMPOSVAL, "dimension mismatch in m_inverse");
		CHsolve(m_cp, v1, v2);
		set_col(out, i, v2);
	}
	v_free(v1);
	v_free(v2);
	m_free(m_cp);
	return(out);
}

MAT *LDLfactor(MAT *M, int *info) { /* do Choleski, as the Lapack R/C interface doesn't LDL' */
	return(CHfactor(M, info));
}

VEC *LDLsolve(MAT *M, VEC *b, VEC *x) { /* same thing here */
	return(CHsolve(M, b, x));
}

VEC *vm_mlt(MAT *m, VEC *v, VEC *out) { /* out <- v m */
	if (m->m != v->dim)
		ErrMsg(ER_IMPOSVAL, "vm_mlt: dimensions");
	out = v_zero(v_resize(out, m->n));
    for (size_t i = 0; i < m->n; i++)
        for (size_t j = 0; j < v->dim; j++)
			out->ve[i] += v->ve[j] * m->me[j][i];
	return(out);
}

VEC *mv_mlt(MAT *m, VEC *v, VEC *out) { /* out <- m v */
	if (v == out)
		ErrMsg(ER_IMPOSVAL, "mv_mlt in situ");
	if (m->n != v->dim)
		ErrMsg(ER_IMPOSVAL, "mv_mlt non-matching sizes");
	out = v_zero(v_resize(out, m->m));
	for (int j = 0; j < m->m; j++)
		for (int i = 0; i < m->n; i++)
			out->ve[j] += m->me[j][i] * v->ve[i];
	return(out);
}

MAT *m_mlt(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 m2 */
	if (m1->n != m2->m)
		ErrMsg(ER_IMPOSVAL, "mv_mlt non-matching sizes");
	out = m_zero(m_resize(out, m1->m, m2->n));
	for (int i = 0; i < m1->m; i++)
		for (int j = 0; j < m2->n; j++)
			for (int k = 0; k < m1->n; k++)
				out->me[i][j] += m1->me[i][k] * m2->me[k][j];
	return(out);
}

VEC *sv_mlt(double s, VEC *v, VEC *out) { /* out <- s * v */
	out = v_resize(out, v->dim);
	for (int i = 0; i < v->dim; i++)
		out->ve[i] = s * v->ve[i];
	return(out);
}

double v_norm2(VEC *v) { /* 2-norm  */
	double sum = 0.0;
	for (int i = 0; i < v->dim; i++)
		sum += v->ve[i] * v->ve[i];
	return(sum);
}

MAT *ms_mltadd(MAT *m1, MAT *m2, double s, MAT *out) { /* out <- m1 + s * m2 */
	/* return m1 + s * m2 */
	if (m1->m != m2->m || m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "ms_mltadd: dimension mismatch");
	out = m_resize(out, m1->m, m1->n);
	for (int i = 0; i < m1->m; i++)
		for (int j = 0; j < m1->n; j++)
			out->me[i][j] = m1->me[i][j] + s * m2->me[i][j];
	return(out);
}

MAT *mtrm_mlt(MAT *m1, MAT *m2, MAT *out) { /* out <- m1'm2 */
	if (m1->m != m2->m)
		ErrMsg(ER_IMPOSVAL, "mtrm_mlt non-matching m arrays");
	out = m_zero(m_resize(out, m1->n, m2->n));
	for (int i = 0; i < m1->n; i++)
		for (int j = 0; j < m2->n; j++)
			for (int k = 0; k < m1->m; k++)
				out->me[i][j] += m1->me[k][i] * m2->me[k][j];
	return(out);
}

MAT *mmtr_mlt(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 m2' */
	if (m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "mmtr_mlt non-matching m arrays");
	out = m_zero(m_resize(out, m1->m, m2->m));
	for (int i = 0; i < m1->m; i++)
		for (int j = 0; j < m2->m; j++)
			out->me[i][j] = __ip__(m1->me[i], m2->me[j], m1->n);
			/* for (int k = 0; k < m1->n; k++)
				out->me[i][j] += m1->me[i][k] * m2->me[j][k];
			*/
	return(out);
}

MAT *m_copy(MAT *in, MAT *out){
	out = m_resize(out, in->m, in->n);
	memcpy(out->v, in->v, in->m * in->n * sizeof(double));
	return(out);
}

VEC *v_copy(VEC *in, VEC *out){
	out = v_resize(out, in->dim);
	memcpy(out->ve, in->ve, in->dim * sizeof(double));
	return(out);
}

VEC *get_row(MAT *m, size_t i, VEC *out) { /* out <- i-th row of m */
	if (i > m->m || out->dim != m->n)
		ErrMsg(ER_IMPOSVAL, "dimension mismatch in set_col");
	for (int j = 0; j < m->n; j++)
		out->ve[j] = m->me[i][j];
	return(out);
}

void QRfactor(MAT *a, VEC *b) {
	ErrMsg(ER_IMPOSVAL, "QRfactor not yet implemented");
}

double QRcondest(MAT *a) {
	ErrMsg(ER_IMPOSVAL, "QRcondest not yet implemented");
	return(0.0);
}

void QRsolve(MAT *m, VEC *a, VEC *b, VEC *c) {
	ErrMsg(ER_IMPOSVAL, "QRsolve not yet implemented");
}

double in_prod(VEC *a, VEC *b) { /* a'b */
	if (a->dim != b->dim)
		ErrMsg(ER_IMPOSVAL, "in_prod: dimensions don't match");
	return(__ip__(a->ve, b->ve, a->dim));
}

double __ip__(double *p1, double *p2, int n) { /* inner product of p1 and p2 */
	double d = 0.0;
	for (int i = 0; i < n; i++)
		d += p1[i] * p2[i];
	return(d);
}

MAT *sm_mlt(double s, MAT *m1, MAT *out) { /* out <- s * m1 */
	out = m_resize(out, m1->m, m1->n);
	for (int i = 0; i < m1->m; i++)
		for (int j = 0; j < m1->n; j++)
			out->me[i][j] = s * m1->me[i][j];
	return(out);
}

VEC *v_add(VEC *v1, VEC *v2, VEC *out) { /* out = v1 + v2 */
	if (v1->dim != v2->dim)
		ErrMsg(ER_IMPOSVAL, "v_sub size mismatch");
	out = v_resize(out, v1->dim);
	for (int i = 0; i < out->dim; i++)
		out->ve[i] = v1->ve[i] + v2->ve[i];
	return(out);
}

VEC *v_sub(VEC *v1, VEC *v2, VEC *out) { /* out = v1 - v2 */
	if (v1->dim != v2->dim)
		ErrMsg(ER_IMPOSVAL, "v_sub size mismatch");
	out = v_resize(out, v1->dim);
	for (int i = 0; i < out->dim; i++)
		out->ve[i] = v1->ve[i] - v2->ve[i];
	return(out);
}

MAT *m_add(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 + m2 */
	if (m1->m != m2->m || m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "m_add size mismatch");
	out = m_resize(out, m1->m, m1->n);
	for (int i = 0; i < m1->m; i++)
		for (int j = 0; j < m1->n; j++)
			out->me[i][j] = m1->me[i][j] + m2->me[i][j];
	return(out);
}

MAT *m_sub(MAT *m1, MAT *m2, MAT *out) { /* out <- m1 - m2 */
	if (m1->m != m2->m || m1->n != m2->n)
		ErrMsg(ER_IMPOSVAL, "m_sub size mismatch");
	out = m_resize(out, m1->m, m1->n);
	for (int i = 0; i < m1->m; i++)
		for (int j = 0; j < m1->n; j++)
			out->me[i][j] = m1->me[i][j] - m2->me[i][j];
	return(out);
}

MAT *m_transp(MAT *in, MAT *out) { /* out <- in' ; may be in situ */
	if (in == out) {
		for (int i = 0; i < in->m; i++)
			for (int j = 0; j < in->n; j++)
				out->me[j][i] = in->me[i][j];
	} else {
		double tmp;
		out = m_resize(out, in->n, in->m);
		for (int i = 1; i < in->n; i++)
			for (int j = 0; j < i; j++) { /* swap: */
				tmp = in->me[i][j];
				in->me[i][j] = in->me[j][i];
				in->me[j][i] = tmp;
			}
	}
	return(out);
}
