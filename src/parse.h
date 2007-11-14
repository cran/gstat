#ifndef BISON_Y_TAB_H
# define BISON_Y_TAB_H

#ifndef YYSTYPE
typedef union {
	int ival;
	unsigned int uval;
	double dval;
	char *sval;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	INT	257
# define	UINT	258
# define	REAL	259
# define	QSTR	260
# define	IDENT	261
# define	ID_DATA	262
# define	ID_X	263
# define	ID_VARIOGRAM	264
# define	ID_PREDICTIONS	265
# define	ID_VARIANCES	266
# define	ID_COVARIANCES	267
# define	ID_OUTPUT	268
# define	ID_MASKS	269
# define	ID_EDGES	270
# define	ID_SET	271
# define	ID_MERGE	272
# define	ID_AREA	273
# define	ID_BLOCK	274
# define	ID_METHOD	275
# define	ID_BOUNDS	276
# define	ID_MARGINALS	277


extern YYSTYPE gstat_yylval;

#endif /* not BISON_Y_TAB_H */
