#ifndef LM_TYPE_H /* avoid multiple inclusion */
# define LM_TYPE_H
typedef struct {
	VEC *beta, /* parameter vector */
		*y, /* data vector */
		*Xty, /* X'y */
		*weights; /* weights in a WLS model: V-1, 1/sigma^2_i */
	MAT *X, /* design matrix */
		*Chol; /* Choleski decomposition of X'X */
	double MSErr, /* Mean Square Error */
		MSReg, /* Mean Square due to regression */
		SSErr, /* Sum of Squares error */
		SSReg, /* Sum of Squares regression */
		cn_max; /* max. allowed condition number; < 0 => don't check */
	int dfE, /* degrees of freedom error */
		dfReg, /* degrees of freedom regression */
		is_singular, /* flag if X'X is singular */
		no_variances, /* for first order estimates only */
		has_intercept; /* model has intercept, J is part of X */
} LM ;
#endif
