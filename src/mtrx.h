#define MATRIXH
/* interface copied from meschach; implementation rewritten from scratch */
typedef struct {
	size_t n, m, max_n, max_m;
	double **me, *v;
} MAT;

typedef	struct	{
	size_t dim, max_n;
	double *ve;
} VEC;

typedef	struct	{
	size_t size, max_size;
	size_t *pe;
} PERM;

typedef	struct	{
	size_t size, max_size;
	int *ive;
} IVEC;

#define PNULL (PERM *) NULL
#define MNULL (MAT *) NULL
#define VNULL (VEC *) NULL
#define IVNULL (IVEC *) NULL

#define M_FREE(x) { m_free(x); x = NULL; }
#define V_FREE(x) { v_free(x); x = NULL; }
void m_free(MAT *m);
void v_free(VEC *v);
void iv_free(IVEC *v);
void px_free(PERM *v);
#define m_get(i,j) m_resize(MNULL, i, j)
#define v_get(i) v_resize(VNULL, i)

MAT *m_resize(MAT *m, size_t new_r, size_t new_c);
VEC *v_resize(VEC *v, size_t new_n);
PERM *px_resize(PERM *p, size_t new_n);
IVEC *iv_resize(IVEC *v, size_t new_n);
MAT *m_zero(MAT *m);
VEC *v_zero(VEC *v);
void set_col(MAT *M, size_t i, VEC *col);
VEC *get_col(MAT *M, size_t i, VEC *to);
MAT *m_inverse(MAT *in, MAT *out, int *info);
MAT *LDLfactor(MAT *M, int *info);
VEC *LDLsolve(MAT *M, VEC *b, VEC *x);
VEC *vm_mlt(MAT *m, VEC *v, VEC *out);
VEC *mv_mlt(MAT *m, VEC *v, VEC *out);
MAT *m_mlt(MAT *m1, MAT *m2, MAT *out);
MAT *mtrm_mlt(MAT *m1, MAT *m2, MAT *out);
VEC *v_sub(VEC *v1, VEC *v2, VEC *out);
MAT *m_sub(MAT *m1, MAT *m2, MAT *out);
VEC *v_add(VEC *v1, VEC *v2, VEC *out);
VEC *sv_mlt(double s, VEC *v1, VEC *v2);
MAT *m_add(MAT *m1, MAT *m2, MAT *out);
MAT *m_copy(MAT *in, MAT *out);
VEC *v_copy(VEC *in, VEC *out);
double v_norm2(VEC *v);
VEC *CHsolve(MAT *A, VEC *b, VEC*out);
MAT *CHfactor(MAT *A, int *info);
double in_prod(VEC *a, VEC *b);
void QRfactor(MAT *a, VEC *b);
double QRcondest(MAT *a);
void QRsolve(MAT *m, VEC *a, VEC *b, VEC *c);
MAT *sm_mlt(double s, MAT *m1, MAT *out);
VEC *get_row(MAT *m, size_t i, VEC *out);
MAT *ms_mltadd(MAT *m1, MAT *m2, double s, MAT *out);
MAT *mmtr_mlt(MAT *m1, MAT *m2, MAT *out);
double __ip__(double *p1, double *p2, int n);
MAT *m_transp(MAT *in, MAT *out);

/* TODO: */
#define MACHEPS 1.19209e-07
