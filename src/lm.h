#ifndef LM_H
# define LM_H
void pred_lm(DATA **data, int n_vars, DPOINT *where, double *est);
void make_residuals_lm(DATA *d);
# ifdef MATRIXH /* MAT, VEC definitions */
MAT *get_X(DATA **d, MAT *X, int nvars);
MAT *get_X0(DATA **d, MAT *X0, DPOINT *where, int nvars);
double calc_mu(const DATA *d, const DPOINT *pt);
VEC *get_y(DATA **d, VEC *y, int nvars);
int is_singular(MAT *X, double epsilon);
void m_logoutput(MAT *a);
void v_logoutput(VEC *x);
# endif
# ifdef LM_TYPE_H
LM *calc_lm(LM *lm);
LM *init_lm(LM *lm);
void logprint_lm(DATA *d, LM *lm);
# endif
void free_lm(void *lm);
#endif /* LM_H */
