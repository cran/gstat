/* A Bison parser, made from parse.y
   by GNU bison 1.35.  */

#define YYBISON 1  /* Identify Bison output.  */

#define yyparse gstat_yyparse
#define yylex gstat_yylex
#define yyerror gstat_yyerror
#define yylval gstat_yylval
#define yychar gstat_yychar
#define yydebug gstat_yydebug
#define yynerrs gstat_yynerrs
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

#line 1 "parse.y"

/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2003 (C) Edzer J. Pebesma

    Edzer J. Pebesma, e.pebesma@geog.uu.nl
    Department of physical geography, Utrecht University
    P.O. Box 80.115, 3508 TC Utrecht, The Netherlands

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

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
 * parse.y: LALR(1) grammar for the gstat command syntax.
 * to make parse.c, type ``make parse.c'', it will use bison or yacc.
 *
 * If you fail (or don't have bison or yacc), then copy the file parse.c_
 * to parse.c. All this can be prevented by running configure while NO_YACC
 * is defined.
 * 
 * The parser assumes that in the function yylex() each identifier is
 * duplicated to ylval.sval, not just a pointer-copy. (some memory loss
 * will occur as a result)
 *
 * hints to extend the parser: 
 * o add a command: copy all from the most similar available command
 *   (add a %token and %type declaration, add a rule, add a return value
 *   from yylex() -> see the IDENT actions in lex.l)
 * o add a ``set'' variable: modify is_set_expr(), glvars.[ch] and defaults.h
 * o add a data() command: modify data.[ch] and is_data_expr()
 * o add a variogram model: vario*.[ch]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"

#ifdef HAVE_UNISTD_H
# include <unistd.h> /* isatty() */
#endif

#include "data.h"
#include "vario.h"
#include "debug.h"
#include "glvars.h"
#include "userio.h"
#include "utils.h"
#include "lex.h"

static DATA *d = NULL, **dpp = NULL;
static DPOINT *bp = NULL;
static VARIOGRAM *v = NULL;
static int id = -1, id1 = -1, id2 = -1, col1 = -1, col2 = -1,
	fit_sill = 0, fit_range = 0, nrangepars = 1,
	vector_only = 0, allow_vector_only = 0;
static double range[NRANGEPARS], anis[5];
static char **ofn = NULL, *boundary_file = NULL;
static VARIOGRAM *parse_variogram = NULL;
static D_VECTOR *sd_vector = NULL;

#ifdef YYBISON
# ifndef __STDC__
#  define __STDC__
/* or else all const's will be defined empty */
# endif
#endif

typedef struct {
	const char *name;
	void *ptr;
	enum { 
		UNKNOWN, 
		IS_INT, 
		IS_UINT, 
		IS_REAL, 
		IS_STRING, 
		IS_D_VECTOR, 
		NO_ARG 
	} what;
	enum { 
		NOLIMIT, 
		GEZERO, 
		GTZERO 
	} limit;
} GSTAT_EXPR;

GSTAT_EXPR expr = { NULL, NULL, UNKNOWN, NOLIMIT };

static void push_data_X(DATA *d, int id);
static int is_data_expr(DATA *d, GSTAT_EXPR *expr, const char *fld);
static int is_set_expr(GSTAT_EXPR *expr, const char *fld);
static int is_block_expr(GSTAT_EXPR *expr, const char *s);
static void push_marginal(char *name, double val);
static void check_assign_expr(GSTAT_EXPR *expr);
static void reset_parser(void);
static void verify_data(DATA *d);

#define gstat_yyerror(s) lex_error()


#line 118 "parse.y"
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
#ifndef YYDEBUG
# define YYDEBUG 0
#endif



#define	YYFINAL		244
#define	YYFLAG		-32768
#define	YYNTBASE	36

/* YYTRANSLATE(YYLEX) -- Bison token number corresponding to YYLEX. */
#define YYTRANSLATE(x) ((unsigned)(x) <= 277 ? yytranslate[x] : 80)

/* YYTRANSLATE[YYLEX] -- Bison token number corresponding to YYLEX. */
static const char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    32,     2,
      30,    31,     2,    33,    27,    34,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    26,    24,
       2,    25,     2,     2,    35,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    28,     2,    29,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23
};

#if YYDEBUG
static const short yyprhs[] =
{
       0,     0,     1,     4,     6,     8,    11,    14,    17,    20,
      23,    26,    29,    32,    35,    38,    41,    44,    46,    48,
      50,    52,    54,    56,    60,    62,    66,    69,    71,    73,
      75,    77,    79,    81,    83,    85,    87,    89,    91,    93,
      95,    97,    99,   101,   103,   107,   109,   114,   118,   122,
     124,   128,   129,   133,   135,   139,   143,   147,   151,   153,
     155,   159,   161,   163,   165,   167,   171,   173,   177,   181,
     183,   187,   189,   193,   199,   207,   211,   217,   221,   226,
     230,   235,   242,   244,   247,   252,   258,   265,   272,   274,
     276,   280,   286,   294,   306,   320,   322,   325,   328,   330,
     333,   337,   344,   348,   355,   359,   368,   371,   373,   377,
     381,   385,   389,   391,   395,   399,   401,   405,   409,   411,
     415,   420,   431,   435,   439,   441,   444,   448,   452,   454,
     458,   460
};
static const short yyrhs[] =
{
      -1,    36,    37,     0,    42,     0,    24,     0,    45,    24,
       0,    58,    24,     0,    75,    24,     0,    71,    24,     0,
      73,    24,     0,    52,    24,     0,    56,    24,     0,    66,
      24,     0,    70,    24,     0,    76,    24,     0,    78,    24,
       0,    67,    24,     0,     3,     0,     5,     0,    25,     0,
      26,     0,    26,     0,    27,     0,    28,    42,    29,     0,
      43,     0,    42,    27,    43,     0,    42,    43,     0,    38,
       0,     7,     0,    19,     0,    20,     0,    22,     0,    13,
       0,     8,     0,    23,     0,    15,     0,    21,     0,    14,
       0,    11,     0,    17,     0,    12,     0,    10,     0,     9,
       0,    16,     0,    46,    26,    47,     0,    46,     0,     8,
      30,    44,    31,     0,     8,    30,    31,     0,     8,    30,
       1,     0,    48,     0,    47,    27,    48,     0,     0,     9,
      25,    49,     0,     6,     0,    51,    25,     3,     0,    51,
      25,     5,     0,    51,    25,     6,     0,    51,    25,    41,
       0,    51,     0,    50,     0,    49,    32,    50,     0,     7,
       0,     3,     0,     7,     0,    20,     0,    20,    26,    53,
       0,    54,     0,    53,    27,    54,     0,    55,    25,    38,
       0,     7,     0,    57,    26,    47,     0,    19,     0,    19,
      30,    31,     0,    59,    26,     6,    40,    60,     0,    59,
      26,     6,    40,     6,    40,    60,     0,    59,    26,     6,
       0,    59,    26,     6,    40,     6,     0,    59,    26,    60,
       0,    59,    26,     1,    24,     0,    10,    30,    31,     0,
      10,    30,    44,    31,     0,    10,    30,    44,    27,    44,
      31,     0,    61,     0,    60,    61,     0,    64,    62,    30,
      31,     0,    64,    62,    30,    63,    31,     0,    33,    64,
      62,    30,    63,    31,     0,    34,    64,    62,    30,    63,
      31,     0,     7,     0,    65,     0,    65,    27,    38,     0,
      65,    27,    38,    27,    38,     0,    65,    27,    38,    27,
      38,    27,    38,     0,    65,    27,    38,    27,    38,    27,
      38,    27,    38,    27,    38,     0,    65,    27,    38,    27,
      38,    27,    38,    27,    38,    27,    38,    27,    38,     0,
      38,     0,    35,    38,     0,    38,    35,     0,    38,     0,
      35,    38,     0,    14,    25,     6,     0,    11,    30,    44,
      31,    26,     6,     0,    11,    26,     6,     0,    12,    30,
      44,    31,    26,     6,     0,    12,    26,     6,     0,    13,
      30,    44,    27,    44,    31,    26,     6,     0,    17,    68,
       0,    68,     0,    69,    39,     3,     0,    69,    39,     4,
       0,    69,    39,     5,     0,    69,    39,     6,     0,     7,
       0,    21,    26,     7,     0,    15,    26,    72,     0,     6,
       0,    72,    27,     6,     0,    16,    26,    74,     0,     6,
       0,    74,    27,     6,     0,    18,    44,     7,    44,     0,
      18,    44,    30,     3,    31,     7,    44,    30,     3,    31,
       0,    22,    26,     6,     0,    22,    26,    77,     0,    38,
       0,    77,    38,     0,    77,    27,    38,     0,    23,    26,
      79,     0,    38,     0,    79,    27,    38,     0,     6,     0,
      79,    27,     6,     0
};

#endif

#if YYDEBUG
/* YYRLINE[YYN] -- source line where rule number YYN was defined. */
static const short yyrline[] =
{
       0,   142,   143,   144,   147,   148,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   162,   163,   165,
     166,   169,   170,   173,   176,   177,   178,   181,   189,   189,
     189,   189,   189,   189,   190,   190,   190,   190,   190,   191,
     191,   191,   191,   191,   195,   196,   199,   205,   209,   212,
     213,   216,   217,   218,   219,   227,   235,   242,   253,   261,
     262,   265,   279,   282,   290,   294,   297,   298,   301,   307,
     310,   313,   317,   323,   324,   325,   326,   327,   328,   334,
     340,   346,   358,   359,   362,   366,   369,   372,   378,   385,
     390,   395,   403,   412,   421,   433,   434,   435,   438,   439,
     442,   443,   448,   456,   461,   469,   481,   482,   485,   493,
     500,   507,   516,   519,   533,   536,   537,   540,   543,   544,
     547,   562,   581,   582,   585,   586,   587,   590,   593,   594,
     595,   596
};
#endif


#if (YYDEBUG) || defined YYERROR_VERBOSE

/* YYTNAME[TOKEN_NUM] -- String name of the token TOKEN_NUM. */
static const char *const yytname[] =
{
  "$", "error", "$undefined.", "INT", "UINT", "REAL", "QSTR", "IDENT", 
  "ID_DATA", "ID_X", "ID_VARIOGRAM", "ID_PREDICTIONS", "ID_VARIANCES", 
  "ID_COVARIANCES", "ID_OUTPUT", "ID_MASKS", "ID_EDGES", "ID_SET", 
  "ID_MERGE", "ID_AREA", "ID_BLOCK", "ID_METHOD", "ID_BOUNDS", 
  "ID_MARGINALS", "';'", "'='", "':'", "','", "'['", "']'", "'('", "')'", 
  "'&'", "'+'", "'-'", "'@'", "input", "command", "val", "assign", 
  "comcol", "d_vector", "d_list", "d_val", "any_id", "data_cmd", 
  "data_decl", "data_cont", "data_exp", "data_X", "data_X_what", 
  "data_what", "block_cmd", "block_cont", "block_exp", "block_lhs", 
  "area_cmd", "area_decl", "vgm_cmd", "vgm_decl", "vgm_cont", "vgm_model", 
  "vgm_model_type", "vgm_range", "sill_val", "range_val", "output_cmd", 
  "set_cmd", "set_exp", "set_lhs", "method_cmd", "mask_cmd", "mask_cont", 
  "edges_cmd", "edges_cont", "merge_cmd", "bounds_cmd", "bounds_exp", 
  "marginals_cmd", "marginals_cont", 0
};
#endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives. */
static const short yyr1[] =
{
       0,    36,    36,    36,    37,    37,    37,    37,    37,    37,
      37,    37,    37,    37,    37,    37,    37,    38,    38,    39,
      39,    40,    40,    41,    42,    42,    42,    43,    44,    44,
      44,    44,    44,    44,    44,    44,    44,    44,    44,    44,
      44,    44,    44,    44,    45,    45,    46,    46,    46,    47,
      47,    48,    48,    48,    48,    48,    48,    48,    48,    49,
      49,    50,    50,    51,    52,    52,    53,    53,    54,    55,
      56,    57,    57,    58,    58,    58,    58,    58,    58,    59,
      59,    59,    60,    60,    61,    61,    61,    61,    62,    63,
      63,    63,    63,    63,    63,    64,    64,    64,    65,    65,
      66,    66,    66,    66,    66,    66,    67,    67,    68,    68,
      68,    68,    69,    70,    71,    72,    72,    73,    74,    74,
      75,    75,    76,    76,    77,    77,    77,    78,    79,    79,
      79,    79
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN. */
static const short yyr2[] =
{
       0,     0,     2,     1,     1,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     1,     1,     1,
       1,     1,     1,     3,     1,     3,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     1,     4,     3,     3,     1,
       3,     0,     3,     1,     3,     3,     3,     3,     1,     1,
       3,     1,     1,     1,     1,     3,     1,     3,     3,     1,
       3,     1,     3,     5,     7,     3,     5,     3,     4,     3,
       4,     6,     1,     2,     4,     5,     6,     6,     1,     1,
       3,     5,     7,    11,    13,     1,     2,     2,     1,     2,
       3,     6,     3,     6,     3,     8,     2,     1,     3,     3,
       3,     3,     1,     3,     3,     1,     3,     3,     1,     3,
       4,    10,     3,     3,     1,     2,     3,     3,     1,     3,
       1,     3
};

/* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
   doesn't specify something else to do.  Zero means the default is an
   error. */
static const short yydefact[] =
{
       1,    17,    18,     0,    27,     3,    24,   112,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    71,    64,
       0,     0,     0,     4,     2,     0,    45,     0,     0,     0,
       0,     0,     0,     0,   107,     0,     0,     0,     0,     0,
       0,     0,     0,    26,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   106,    28,    33,    42,    41,    38,
      40,    32,    37,    35,    43,    39,    29,    30,    36,    31,
      34,     0,     0,     0,     0,     0,     0,     5,    51,    10,
      11,    51,     6,     0,    12,    16,    19,    20,     0,    13,
       8,     9,     7,    14,    15,    25,    48,    47,     0,    79,
       0,   102,     0,   104,     0,     0,   100,   115,   114,   118,
     117,     0,     0,    72,    69,    65,    66,     0,   113,   122,
     124,   123,   130,   128,   127,    53,    63,     0,    44,    49,
      58,    70,     0,    75,     0,     0,     0,    95,    77,    82,
       0,   108,   109,   110,   111,    46,     0,    80,     0,     0,
       0,     0,     0,   120,     0,     0,     0,     0,   125,     0,
       0,    51,     0,    78,    21,    22,     0,     0,     0,    96,
      97,    83,    88,     0,     0,     0,     0,     0,   116,   119,
       0,    67,    68,   126,   131,   129,    62,    61,    52,    59,
      50,    54,    55,    56,     0,    57,    76,    73,     0,     0,
       0,    81,   101,   103,     0,     0,     0,     0,     0,     0,
       0,    84,     0,    98,     0,    89,     0,     0,    60,    23,
      74,     0,     0,    99,    85,     0,   105,     0,    86,    87,
      90,     0,     0,   121,    91,     0,    92,     0,     0,     0,
      93,     0,    94,     0,     0
};

static const short yydefgoto[] =
{
       3,    24,   137,    88,   166,   195,     5,     6,    71,    25,
      26,   128,   129,   188,   189,   130,    27,   115,   116,   117,
      28,    29,    30,    31,   138,   139,   173,   214,   140,   215,
      32,    33,    34,    35,    36,    37,   108,    38,   110,    39,
      40,   121,    41,   124
};

static const short yypact[] =
{
      69,-32768,-32768,   127,-32768,    14,-32768,-32768,    -9,    -1,
     -12,    38,     2,    53,     9,    59,    90,   207,    78,    98,
     104,   105,   106,-32768,-32768,    92,   107,   128,   129,   132,
     130,   134,   131,   137,-32768,    33,   138,   139,   152,   159,
     160,   161,    69,-32768,    79,   158,   176,   207,   180,   207,
     207,   181,   182,   184,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,    39,   162,   185,   188,   101,   112,-32768,   113,-32768,
  -32768,   113,-32768,    10,-32768,-32768,-32768,-32768,   108,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,   165,-32768,
      40,-32768,   166,-32768,   167,   164,-32768,-32768,   172,-32768,
     174,   207,   199,-32768,-32768,   177,-32768,   183,-32768,-32768,
  -32768,    76,-32768,-32768,   178,-32768,-32768,   186,   179,-32768,
     206,   179,   189,    56,    25,    25,    69,   198,    21,-32768,
     227,-32768,-32768,-32768,-32768,-32768,   207,-32768,   210,   212,
     207,   234,   236,-32768,   213,   185,    69,    69,-32768,   120,
      70,   113,    34,-32768,-32768,-32768,    17,   227,   227,-32768,
  -32768,-32768,-32768,   215,   216,   237,   240,   217,-32768,-32768,
     242,-32768,-32768,-32768,-32768,-32768,-32768,-32768,   218,-32768,
  -32768,-32768,-32768,-32768,    69,-32768,    56,    21,   221,   222,
      22,-32768,-32768,-32768,   228,   207,    70,     7,    21,    28,
      28,-32768,    69,-32768,   224,   226,   250,   229,-32768,-32768,
      21,   230,   231,-32768,-32768,    69,-32768,   254,-32768,-32768,
     233,   232,    69,-32768,   238,    69,   239,    69,   241,    69,
     243,    69,-32768,   258,-32768
};

static const short yypgoto[] =
{
  -32768,-32768,     0,-32768,    68,-32768,    73,    -4,   -41,-32768,
  -32768,   190,   111,-32768,    63,-32768,-32768,-32768,   118,-32768,
  -32768,-32768,-32768,-32768,  -159,  -136,  -120,  -144,    -6,-32768,
  -32768,-32768,   259,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768
};


#define	YYLAST		275


static const short yytable[] =
{
       4,    43,   171,    98,   100,     4,   102,   197,   104,   105,
       1,   132,     2,     1,    46,     2,   133,     1,    47,     2,
       1,    44,     2,   196,     1,     1,     2,     2,     1,    45,
       2,     1,    50,     2,    42,    52,   219,   191,    95,   192,
     193,    42,     4,   134,   135,   136,   111,   198,   199,   220,
     134,   135,   136,   211,   134,   135,   136,   212,    86,    87,
     136,   171,   194,   212,    48,   221,   222,   146,    49,   112,
     153,   147,     1,   186,     2,   120,   123,   187,    51,     1,
      96,     2,   164,   165,   171,    53,    55,    56,    57,    58,
      59,    60,    61,    62,    63,    64,    65,     7,    66,    67,
      68,    69,    70,   157,     1,   174,     2,   119,    72,   177,
      97,   141,   142,   143,   144,     1,    77,     2,   122,   125,
     126,   158,   127,     1,    73,     2,   184,   243,   167,   168,
      74,    75,    76,    78,     7,     8,   169,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    79,    80,    82,    84,   182,   183,    81,   185,
      83,    85,    89,    90,   217,    55,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    91,    66,    67,    68,
      69,    70,   101,    92,    93,    94,   103,   106,   107,    99,
     109,   150,   114,   113,     4,   118,   145,   148,   149,   151,
     213,   152,   154,    43,   155,   159,   161,     4,   156,   213,
     213,   160,   223,   163,    55,    56,    57,    58,    59,    60,
      61,    62,    63,    64,    65,   230,    66,    67,    68,    69,
      70,   162,   234,   170,   172,   236,   175,   238,   176,   240,
     178,   242,   179,   202,   180,   200,   203,   201,   204,   205,
     206,   209,   210,   225,   216,   224,   226,   231,   244,   227,
     232,   228,   229,   233,   208,   235,   237,   207,   239,   218,
     241,   131,   190,   181,     0,    54
};

static const short yycheck[] =
{
       0,     5,   138,    44,    45,     5,    47,   166,    49,    50,
       3,     1,     5,     3,    26,     5,     6,     3,    30,     5,
       3,    30,     5,     6,     3,     3,     5,     5,     3,    30,
       5,     3,    30,     5,    27,    26,    29,     3,    42,     5,
       6,    27,    42,    33,    34,    35,     7,   167,   168,   208,
      33,    34,    35,    31,    33,    34,    35,    35,    25,    26,
      35,   197,    28,    35,    26,   209,   210,    27,    30,    30,
     111,    31,     3,     3,     5,    75,    76,     7,    25,     3,
       1,     5,    26,    27,   220,    26,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,     7,    19,    20,
      21,    22,    23,    27,     3,   146,     5,     6,    30,   150,
      31,     3,     4,     5,     6,     3,    24,     5,     6,     6,
       7,   121,     9,     3,    26,     5,     6,     0,   134,   135,
      26,    26,    26,    26,     7,     8,   136,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    24,    24,    24,    24,   156,   157,    26,   159,
      26,    24,    24,    24,   205,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    24,    19,    20,    21,
      22,    23,     6,    24,    24,    24,     6,     6,     6,    31,
       6,    27,     7,    31,   194,     7,    31,    31,    31,    27,
     200,    27,     3,   207,    27,    27,    27,   207,    25,   209,
     210,    25,   212,    24,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,   225,    19,    20,    21,    22,
      23,    25,   232,    35,     7,   235,    26,   237,    26,   239,
       6,   241,     6,     6,    31,    30,     6,    31,    31,     7,
      32,    30,    30,    27,    26,    31,     6,     3,     0,    30,
      27,    31,    31,    31,   196,    27,    27,   194,    27,   206,
      27,    81,   161,   155,    -1,    16
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison/bison.simple"

/* Skeleton output parser for bison,

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software
   Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser when
   the %semantic_parser declaration is not specified in the grammar.
   It was written by Richard Stallman by simplifying the hairy parser
   used when %semantic_parser is specified.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

#if ! defined (yyoverflow) || defined (YYERROR_VERBOSE)

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || defined (YYERROR_VERBOSE) */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYLTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
# if YYLSP_NEEDED
  YYLTYPE yyls;
# endif
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAX (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# if YYLSP_NEEDED
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE) + sizeof (YYLTYPE))	\
      + 2 * YYSTACK_GAP_MAX)
# else
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAX)
# endif

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAX;	\
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif


#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");			\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).

   When YYLLOC_DEFAULT is run, CURRENT is set the location of the
   first token.  By default, to implement support for ranges, extend
   its range to the last symbol.  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)       	\
   Current.last_line   = Rhs[N].last_line;	\
   Current.last_column = Rhs[N].last_column;
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#if YYPURE
# if YYLSP_NEEDED
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, &yylloc, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval, &yylloc)
#  endif
# else /* !YYLSP_NEEDED */
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval)
#  endif
# endif /* !YYLSP_NEEDED */
#else /* !YYPURE */
# define YYLEX			yylex ()
#endif /* !YYPURE */


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)
/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif

#ifdef YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif
#endif

#line 315 "/usr/share/bison/bison.simple"


/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
#  define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL
# else
#  define YYPARSE_PARAM_ARG YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
# endif
#else /* !YYPARSE_PARAM */
# define YYPARSE_PARAM_ARG
# define YYPARSE_PARAM_DECL
#endif /* !YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
# ifdef YYPARSE_PARAM
int yyparse (void *);
# else
int yyparse (void);
# endif
#endif

/* YY_DECL_VARIABLES -- depending whether we use a pure parser,
   variables are global, or local to YYPARSE.  */

#define YY_DECL_NON_LSP_VARIABLES			\
/* The lookahead symbol.  */				\
int yychar;						\
							\
/* The semantic value of the lookahead symbol. */	\
YYSTYPE yylval;						\
							\
/* Number of parse errors so far.  */			\
int yynerrs;

#if YYLSP_NEEDED
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES			\
						\
/* Location data for the lookahead symbol.  */	\
YYLTYPE yylloc;
#else
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES
#endif


/* If nonreentrant, generate the variables here. */

#if !YYPURE
YY_DECL_VARIABLES
#endif  /* !YYPURE */

int
yyparse (YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  /* If reentrant, generate the variables here. */
#if YYPURE
  YY_DECL_VARIABLES
#endif  /* !YYPURE */

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yychar1 = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack. */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;

#if YYLSP_NEEDED
  /* The location stack.  */
  YYLTYPE yylsa[YYINITDEPTH];
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;
#endif

#if YYLSP_NEEDED
# define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
# define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  YYSIZE_T yystacksize = YYINITDEPTH;


  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
#if YYLSP_NEEDED
  YYLTYPE yyloc;
#endif

  /* When reducing, the number of symbols on the RHS of the reduced
     rule. */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;
#if YYLSP_NEEDED
  yylsp = yyls;
#endif
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  */
# if YYLSP_NEEDED
	YYLTYPE *yyls1 = yyls;
	/* This used to be a conditional around just the two extra args,
	   but that might be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yyls1, yysize * sizeof (*yylsp),
		    &yystacksize);
	yyls = yyls1;
# else
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);
# endif
	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);
# if YYLSP_NEEDED
	YYSTACK_RELOCATE (yyls);
# endif
# undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
#if YYLSP_NEEDED
      yylsp = yyls + yysize - 1;
#endif

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yychar1 = YYTRANSLATE (yychar);

#if YYDEBUG
     /* We have to keep this `#if YYDEBUG', since we use variables
	which are defined only if `YYDEBUG' is set.  */
      if (yydebug)
	{
	  YYFPRINTF (stderr, "Next token is %d (%s",
		     yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise
	     meaning of a token, for further debugging info.  */
# ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
# endif
	  YYFPRINTF (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %d (%s), ",
	      yychar, yytname[yychar1]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to the semantic value of
     the lookahead token.  This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

#if YYLSP_NEEDED
  /* Similarly for the default location.  Let the user run additional
     commands if for instance locations are ranges.  */
  yyloc = yylsp[1-yylen];
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
#endif

#if YYDEBUG
  /* We have to keep this `#if YYDEBUG', since we use variables which
     are defined only if `YYDEBUG' is set.  */
  if (yydebug)
    {
      int yyi;

      YYFPRINTF (stderr, "Reducing via rule %d (line %d), ",
		 yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (yyi = yyprhs[yyn]; yyrhs[yyi] > 0; yyi++)
	YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
      YYFPRINTF (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif

  switch (yyn) {

case 1:
#line 142 "parse.y"
{ ; }
    break;
case 2:
#line 143 "parse.y"
{ reset_parser(); }
    break;
case 3:
#line 144 "parse.y"
{ vector_only = 1; }
    break;
case 4:
#line 147 "parse.y"
{ ; }
    break;
case 5:
#line 148 "parse.y"
{ ; }
    break;
case 6:
#line 149 "parse.y"
{ update_variogram(v); }
    break;
case 7:
#line 150 "parse.y"
{ ; }
    break;
case 8:
#line 151 "parse.y"
{ ; }
    break;
case 9:
#line 152 "parse.y"
{ ; }
    break;
case 10:
#line 153 "parse.y"
{ ; }
    break;
case 11:
#line 154 "parse.y"
{ ; }
    break;
case 12:
#line 155 "parse.y"
{ ; }
    break;
case 13:
#line 156 "parse.y"
{ ; }
    break;
case 14:
#line 157 "parse.y"
{ ; }
    break;
case 15:
#line 158 "parse.y"
{ ; }
    break;
case 16:
#line 159 "parse.y"
{ ; }
    break;
case 17:
#line 162 "parse.y"
{ yyval.dval = (double) yyvsp[0].ival; }
    break;
case 19:
#line 165 "parse.y"
{ ; }
    break;
case 20:
#line 166 "parse.y"
{ ; }
    break;
case 21:
#line 169 "parse.y"
{ ; }
    break;
case 22:
#line 170 "parse.y"
{ ; }
    break;
case 23:
#line 173 "parse.y"
{ ; }
    break;
case 27:
#line 181 "parse.y"
{ 
			if (d == NULL)
				sd_vector = push_to_vector(yyvsp[0].dval, sd_vector);
			else
				d->beta = push_to_vector(yyvsp[0].dval, d->beta);
		}
    break;
case 43:
#line 192 "parse.y"
{ ; }
    break;
case 44:
#line 195 "parse.y"
{ verify_data(d); }
    break;
case 45:
#line 196 "parse.y"
{ d->dummy = 1; }
    break;
case 46:
#line 199 "parse.y"
{
			id = which_identifier(yyvsp[-1].sval);
			dpp = get_gstat_data();
			d = dpp[id];
			d->id = id;
		}
    break;
case 47:
#line 205 "parse.y"
{
			d = get_dataval();
			d->id = ID_OF_VALDATA;
		}
    break;
case 48:
#line 209 "parse.y"
{ ErrMsg(ER_SYNTAX, "invalid identifier"); }
    break;
case 51:
#line 216 "parse.y"
{ ; }
    break;
case 53:
#line 218 "parse.y"
{ d->fname = yyvsp[0].sval; }
    break;
case 54:
#line 219 "parse.y"
{
			switch (expr.what) { 
				case IS_INT: *((int *)expr.ptr) = yyvsp[0].ival; break;
				case IS_REAL: *((double *)expr.ptr) = (double) yyvsp[0].ival; break;
				default: lex_error(); YYERROR; break;
			}
			check_assign_expr(&expr);
		}
    break;
case 55:
#line 227 "parse.y"
{
			if (expr.what != IS_REAL) {
				lex_error();
				YYERROR;
			}
			*((double *)expr.ptr) = yyvsp[0].dval;
			check_assign_expr(&expr);
		}
    break;
case 56:
#line 235 "parse.y"
{
			if (expr.what != IS_STRING) {
				lex_error();
				YYERROR;
			}
			*((char **)expr.ptr) = yyvsp[0].sval;
		}
    break;
case 57:
#line 242 "parse.y"
{
			if (expr.what != IS_D_VECTOR) {
				lex_error();
				YYERROR;
			}
			/*
			*((D_VECTOR **)expr.ptr) = sd_vector; 
			printf("[[ %d ]]\n", sd_vector->size); 
			sd_vector = NULL;
			*/
		}
    break;
case 58:
#line 253 "parse.y"
{
			if (expr.what != NO_ARG) {
				lex_error();
				YYERROR;
			}
		}
    break;
case 61:
#line 265 "parse.y"
{
			for (id = 0; id < N_POLY; id++) {
				if (almost_equals(yyvsp[0].sval, polynomial[id].name)) {
					id += POLY_MIN;
					break; /* i-loop */
				}
			}
			if (id < 0)
				data_add_X(d, id);
			else {
				lex_error();
				YYERROR;
			}
		}
    break;
case 62:
#line 279 "parse.y"
{ push_data_X(d, yyvsp[0].ival); }
    break;
case 63:
#line 282 "parse.y"
{
			if (! is_data_expr(d, &expr, yyvsp[0].sval)) {
				lex_error();
				YYERROR;
			}
		}
    break;
case 64:
#line 290 "parse.y"
{
			bp = get_block_p();
			bp->x = -1.0; /* will be set to grid cell size in predict.c */
		}
    break;
case 68:
#line 301 "parse.y"
{
			*((double *)expr.ptr) = yyvsp[0].dval;
			check_assign_expr(&expr);
		}
    break;
case 69:
#line 307 "parse.y"
{ if (! is_block_expr(&expr, yyvsp[0].sval)) { lex_error(); YYERROR; }}
    break;
case 70:
#line 310 "parse.y"
{ ; }
    break;
case 71:
#line 313 "parse.y"
{
			d = create_data_area();
			d->id = ID_OF_AREA;
		}
    break;
case 72:
#line 317 "parse.y"
{
			d = create_data_area();
			d->id = ID_OF_AREA;
		}
    break;
case 73:
#line 323 "parse.y"
{ v->fname = yyvsp[-2].sval; }
    break;
case 74:
#line 324 "parse.y"
{ v->fname = yyvsp[-4].sval; v->fname2 = yyvsp[-2].sval; }
    break;
case 75:
#line 325 "parse.y"
{v->fname = yyvsp[0].sval; }
    break;
case 76:
#line 326 "parse.y"
{v->fname = yyvsp[-2].sval; v->fname2 = yyvsp[0].sval; }
    break;
case 78:
#line 328 "parse.y"
{ 
			/* this will eat the ';' as well, but we're bailing out anyway: */
			YYERROR; 
		}
    break;
case 79:
#line 334 "parse.y"
{ 
			/* only allow this when called through read_variogram(): */
			assert(parse_variogram != NULL);
			v = parse_variogram;
			v->n_models = v->n_fit = 0;
		}
    break;
case 80:
#line 340 "parse.y"
{
			id = which_identifier(yyvsp[-1].sval);
			v = get_vgm(LTI(id,id));
			v->id = v->id1 = v->id2 = id;
			v->n_models = v->n_fit = 0;
		}
    break;
case 81:
#line 346 "parse.y"
{
			id1 = which_identifier(yyvsp[-3].sval);
			id2 = which_identifier(yyvsp[-1].sval);
			id = LTI(id1,id2);
			v = get_vgm(id);
			v->id = id;
			v->id1 = id1;
			v->id2 = id2;
			v->n_models = v->n_fit = 0;
		}
    break;
case 84:
#line 362 "parse.y"
{
			range[0] = 0.0;
			push_to_v(v, yyvsp[-2].sval, yyvsp[-3].dval, range, 1, NULL, fit_sill, fit_range);
		}
    break;
case 85:
#line 366 "parse.y"
{
			push_to_v(v, yyvsp[-3].sval, yyvsp[-4].dval, range, nrangepars, anis, fit_sill, fit_range);
		}
    break;
case 86:
#line 369 "parse.y"
{
			push_to_v(v, yyvsp[-3].sval, yyvsp[-4].dval, range, nrangepars, anis, fit_sill, fit_range);
		}
    break;
case 87:
#line 372 "parse.y"
{
			push_to_v(v, yyvsp[-3].sval, -1.0 * yyvsp[-4].dval, range, nrangepars, anis, 
				fit_sill, fit_range);
		}
    break;
case 88:
#line 378 "parse.y"
{ 
			if (which_variogram_model(yyvsp[0].sval) == NOT_SP) {
				lex_error(); YYERROR;
			}
	}
    break;
case 89:
#line 385 "parse.y"
{ 
			range[0] = yyvsp[0].dval; 
			nrangepars = 1;
			anis[0] = -9999.0; 
		}
    break;
case 90:
#line 390 "parse.y"
{
			range[0] = yyvsp[-2].dval;
			range[1] = yyvsp[0].dval;
			nrangepars = 2;
		}
    break;
case 91:
#line 395 "parse.y"
{
			range[0] = yyvsp[-4].dval;
			nrangepars = 1;
			anis[0] = yyvsp[-2].dval;
			anis[3] = yyvsp[0].dval;
			anis[1] = anis[2] = 0.0;
			anis[4] = 1.0;
		}
    break;
case 92:
#line 403 "parse.y"
{
			range[0] = yyvsp[-6].dval;
			range[1] = yyvsp[-4].dval;
			nrangepars = 2;
			anis[0] = yyvsp[-2].dval;
			anis[3] = yyvsp[0].dval;
			anis[1] = anis[2] = 0.0;
			anis[4] = 1.0;
		}
    break;
case 93:
#line 412 "parse.y"
{
			range[0] = yyvsp[-10].dval;
			nrangepars = 1;
			anis[0] = yyvsp[-8].dval;
			anis[1] = yyvsp[-6].dval;
			anis[2] = yyvsp[-4].dval;
			anis[3] = yyvsp[-2].dval;
			anis[4] = yyvsp[0].dval;
		}
    break;
case 94:
#line 421 "parse.y"
{
			range[0] = yyvsp[-12].dval;
			range[1] = yyvsp[-10].dval;
			nrangepars = 2;
			anis[0] = yyvsp[-8].dval;
			anis[1] = yyvsp[-6].dval;
			anis[2] = yyvsp[-4].dval;
			anis[3] = yyvsp[-2].dval;
			anis[4] = yyvsp[0].dval;
		}
    break;
case 95:
#line 433 "parse.y"
{ fit_sill = 1; }
    break;
case 96:
#line 434 "parse.y"
{ fit_sill = 0; yyval.dval = yyvsp[0].dval; }
    break;
case 97:
#line 435 "parse.y"
{ fit_sill = 0; yyval.dval = yyvsp[-1].dval; }
    break;
case 98:
#line 438 "parse.y"
{ fit_range = 1; }
    break;
case 99:
#line 439 "parse.y"
{ fit_range = 0; yyval.dval = yyvsp[0].dval; }
    break;
case 100:
#line 442 "parse.y"
{ o_filename = yyvsp[0].sval; }
    break;
case 101:
#line 443 "parse.y"
{ 
			id = which_identifier(yyvsp[-3].sval);
			ofn = (char **) get_outfile_name();
			ofn[2 * id] = yyvsp[0].sval;
		}
    break;
case 102:
#line 448 "parse.y"
{ 
			if (get_n_vars() == 0) {
				lex_error();
				ErrMsg(ER_SYNTAX, "define data first");
			}
			ofn = (char **) get_outfile_name();
			ofn[0] = yyvsp[0].sval;
		}
    break;
case 103:
#line 456 "parse.y"
{ 
			id = which_identifier(yyvsp[-3].sval);
			ofn = (char **) get_outfile_name();
			ofn[2 * id + 1] = yyvsp[0].sval;
		}
    break;
case 104:
#line 461 "parse.y"
{ 
			if (get_n_vars() == 0) {
				lex_error();
				ErrMsg(ER_SYNTAX, "define data first");
			}
			ofn = (char **) get_outfile_name();
			ofn[1] = yyvsp[0].sval;
		}
    break;
case 105:
#line 469 "parse.y"
{ 
			id = get_n_vars();
			id1 = which_identifier(yyvsp[-5].sval);
			id2 = which_identifier(yyvsp[-3].sval);
			if (id != get_n_vars())	
				ErrMsg(ER_SYNTAX, "define all data(..) before covariances(..,..)");
			ofn = (char **) get_outfile_name();
			id = 2 * id + LTI2(id1, id2);
			ofn[id] = yyvsp[0].sval;
		}
    break;
case 108:
#line 485 "parse.y"
{
			switch (expr.what) { 
				case IS_INT: *((int *)expr.ptr) = yyvsp[0].ival; break;
				case IS_REAL: *((double *)expr.ptr) = (double) yyvsp[0].ival; break;
				default: lex_error(); YYERROR;
			}
			check_assign_expr(&expr);
		}
    break;
case 109:
#line 493 "parse.y"
{
			switch (expr.what) { 
				case IS_UINT: *((unsigned int *)expr.ptr) = yyvsp[0].uval; break;
				default: lex_error(); YYERROR;
			}
			check_assign_expr(&expr);
		}
    break;
case 110:
#line 500 "parse.y"
{
			switch (expr.what) {
				case IS_REAL: *((double *)expr.ptr) = yyvsp[0].dval; break;
				default: lex_error(); YYERROR; break;
			}
			check_assign_expr(&expr);
		}
    break;
case 111:
#line 507 "parse.y"
{
			if (expr.what != IS_STRING) {
				lex_error();
				YYERROR;
			}
			*((char **) expr.ptr) = yyvsp[0].sval;
		}
    break;
case 112:
#line 516 "parse.y"
{ if (! is_set_expr(&expr, yyvsp[0].sval)) { lex_error(); YYERROR; }}
    break;
case 113:
#line 519 "parse.y"
{
			for (id = 1; methods[id].name != NULL; id++) {
				if (almost_equals(yyvsp[0].sval, methods[id].name)) {
					set_method(methods[id].m);
					break; /* id-loop */
				}
			}
			if (methods[id].m == NSP) {
				lex_error();
				YYERROR;
			}
		}
    break;
case 115:
#line 536 "parse.y"
{ push_mask_name(yyvsp[0].sval); }
    break;
case 116:
#line 537 "parse.y"
{ push_mask_name(yyvsp[0].sval); }
    break;
case 118:
#line 543 "parse.y"
{ push_edges_name(yyvsp[0].sval); }
    break;
case 119:
#line 544 "parse.y"
{ push_edges_name(yyvsp[0].sval); }
    break;
case 120:
#line 547 "parse.y"
{
			if (!almost_equals(yyvsp[-1].sval, "w$ith"))
				lex_error();
			id1 = which_identifier(yyvsp[-2].sval);
			id2 = which_identifier(yyvsp[0].sval);
			col1 = col2 = 0;
			dpp = get_gstat_data();
			if (id1 < id2) { /* swap id's */
				id = id1; id1 = id2; id2 = id;
			}
			if (push_to_merge_table(dpp[id1], id2, col1, col2)) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
			}
		}
    break;
case 121:
#line 562 "parse.y"
{
			if (!almost_equals(yyvsp[-4].sval, "w$ith"))
				lex_error();
			id1 = which_identifier(yyvsp[-8].sval);
			id2 = which_identifier(yyvsp[-3].sval);
			col1 = yyvsp[-6].ival;
			col2 = yyvsp[-1].ival;
			dpp = get_gstat_data();
			if (id1 < id2) { /* swap id and col */
				id = id1; id1 = id2; id2 = id;
				id = col1; col1 = col2; col2 = id;
			}
			if (push_to_merge_table(dpp[id1], id2, col1, col2)) {
				lex_error();
				ErrMsg(ER_IMPOSVAL, "attempt to merge failed");
			}
		}
    break;
case 122:
#line 581 "parse.y"
{ boundary_file = yyvsp[0].sval; }
    break;
case 124:
#line 585 "parse.y"
{ push_bound(yyvsp[0].dval); }
    break;
case 125:
#line 586 "parse.y"
{ push_bound(yyvsp[0].dval); }
    break;
case 126:
#line 587 "parse.y"
{ push_bound(yyvsp[0].dval); }
    break;
case 128:
#line 593 "parse.y"
{ push_marginal(NULL, yyvsp[0].dval); }
    break;
case 129:
#line 594 "parse.y"
{ push_marginal(NULL, yyvsp[0].dval); }
    break;
case 130:
#line 595 "parse.y"
{ push_marginal(yyvsp[0].sval, -1.0); }
    break;
case 131:
#line 596 "parse.y"
{ push_marginal(yyvsp[0].sval, -1.0); }
    break;
}

#line 705 "/usr/share/bison/bison.simple"


  yyvsp -= yylen;
  yyssp -= yylen;
#if YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;
#if YYLSP_NEEDED
  *++yylsp = yyloc;
#endif

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("parse error, unexpected ") + 1;
	  yysize += yystrlen (yytname[YYTRANSLATE (yychar)]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "parse error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[YYTRANSLATE (yychar)]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exhausted");
	}
      else
#endif /* defined (YYERROR_VERBOSE) */
	yyerror ("parse error");
    }
  goto yyerrlab1;


/*--------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action |
`--------------------------------------------------*/
yyerrlab1:
  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;
      YYDPRINTF ((stderr, "Discarding token %d (%s).\n",
		  yychar, yytname[yychar1]));
      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;


/*-------------------------------------------------------------------.
| yyerrdefault -- current state does not do anything special for the |
| error token.                                                       |
`-------------------------------------------------------------------*/
yyerrdefault:
#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */

  /* If its default is to accept any token, ok.  Otherwise pop it.  */
  yyn = yydefact[yystate];
  if (yyn)
    goto yydefault;
#endif


/*---------------------------------------------------------------.
| yyerrpop -- pop the current state because it cannot handle the |
| error token                                                    |
`---------------------------------------------------------------*/
yyerrpop:
  if (yyssp == yyss)
    YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#if YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "Error: state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

/*--------------.
| yyerrhandle.  |
`--------------*/
yyerrhandle:
  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

/*---------------------------------------------.
| yyoverflowab -- parser overflow comes here.  |
`---------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}
#line 599 "parse.y"


static int is_data_expr(DATA *d, GSTAT_EXPR *expr, const char *fld) {
#define TABLE_SIZE 32
	GSTAT_EXPR data_options[TABLE_SIZE];
	int i = 0;
#define fill_table(n, p, w, l) \
 data_options[i].name = n; data_options[i].ptr = p; \
 data_options[i].what = w; data_options[i].limit = l; i++;

/* set up table: */
	fill_table("x",          &(d->colnx),        IS_INT, GTZERO  )
	fill_table("y",          &(d->colny),        IS_INT, GTZERO  )
	fill_table("z",          &(d->colnz),        IS_INT, GTZERO  )
	fill_table("v",          &(d->colnvalue),    IS_INT, GTZERO  )
	fill_table("V",          &(d->colnvariance), IS_INT, GTZERO  )
	fill_table("d",          &(d->polynomial_degree), IS_INT, GEZERO  )
	fill_table("max",        &(d->sel_max),      IS_INT, GEZERO  )
	fill_table("omax",       &(d->oct_max),      IS_INT, GEZERO  )
	fill_table("min",        &(d->sel_min),      IS_INT, GEZERO  )
	fill_table("n$max",      &(d->init_max),     IS_INT, GEZERO  )
	fill_table("togrid",     &(d->togrid),       IS_INT,  GEZERO )
	fill_table("I",          &(d->Icutoff),      IS_REAL, NOLIMIT )
	fill_table("mv",         &(d->mv),           IS_REAL, NOLIMIT )
	fill_table("rad$ius",    &(d->sel_rad),      IS_REAL, GTZERO  )
	fill_table("dX",         &(d->dX),           IS_REAL, GEZERO  )
	fill_table("b$eta",      &(d->beta),         IS_D_VECTOR, NOLIMIT )
	fill_table("stan$dard",  &(d->standard),     NO_ARG, NOLIMIT )
	fill_table("log",        &(d->log),          NO_ARG, NOLIMIT )
	fill_table("av$erage",   &(d->average),      NO_ARG, NOLIMIT )
	fill_table("re$gion",    &(d->region),       NO_ARG, NOLIMIT )
	fill_table("du$mmy",     &(d->dummy),        NO_ARG, NOLIMIT )
	fill_table("res$idual",  &(d->calc_residuals), NO_ARG, NOLIMIT )
	fill_table("vdist",      &(d->vdist),        NO_ARG, NOLIMIT )
	fill_table("force",      &(d->force),        NO_ARG, NOLIMIT )
	fill_table("Cat$egory",   &(d->Category),    IS_STRING, NOLIMIT )
	fill_table("ID",         &(d->coln_id),      IS_INT, GTZERO  )
	fill_table("VarF$unction", &(d->var_fn_str), IS_STRING, NOLIMIT  )
	fill_table("nscore",     &(d->nscore_table), IS_STRING, NOLIMIT )
	fill_table("every",      &(d->every),        IS_INT, GTZERO )
	fill_table("offset",     &(d->offset),       IS_INT, NOLIMIT )
	fill_table("prob",       &(d->prob),         IS_REAL, GTZERO )
	fill_table(NULL, NULL, IS_INT, NOLIMIT )

/* check TABLE_SIZE was set correctly... */
	assert(i == TABLE_SIZE); 

	expr->ptr = NULL;
	expr->what = UNKNOWN;
	expr->limit = NOLIMIT;

	for (i = 0; data_options[i].name != NULL; i++) {
		if (almost_equals(fld, data_options[i].name)) {
			expr->name = fld;
			expr->ptr = data_options[i].ptr;
			expr->what = data_options[i].what;
			expr->limit = data_options[i].limit;
			if (expr->what == NO_ARG)
				*((int *) expr->ptr) = 1;
			return 1;
		}
	}

	/* non-standard cases not in data_options[] table: */
	if (almost_equals(fld, "s$tratum")) {
		if (d->id != ID_OF_VALDATA)
			return 0;
		expr->ptr = &(d->colns); 
		expr->what = IS_INT; 
		expr->limit = GTZERO;
		d->what_is_u = U_ISSTRATUM;
	} else if (almost_equals(fld, "av$erage")) {
		d->average = 1; expr->what = NO_ARG;
	} else if (almost_equals(fld, "noav$erage")) {
		d->average = 0; expr->what = NO_ARG;
	} else if (almost_equals(fld, "nores$idual")) {
		d->calc_residuals = 0; expr->what = NO_ARG;
	} else if (almost_equals(fld, "square")) {
		d->square = 1; expr->what = NO_ARG;
	} else if (almost_equals(fld, "c")) {
		pr_warning("use `v' instead of `c' in data definition");
	} else if (almost_equals(fld, "sk_mean")) { /* move it to beta: */
		d->beta = NULL;
		d->beta = push_to_vector(-9999.0, d->beta);
		expr->ptr = &(d->beta->val[0]);
		expr->what = IS_REAL; 
		expr->limit = NOLIMIT;
	} 

	return (expr->what != UNKNOWN);
}

static int is_set_expr(GSTAT_EXPR *expr, const char *name) {
/*
 * parse sequences like `set zmap = 50.0;' or `set zmap = 50, idp = 2.5;'
 * (int, float or string)
 */
	int i;

	const GSTAT_EXPR set_options[] = {
	{ "cn$_max",        &gl_cn_max,       IS_REAL, GTZERO  },
	{ "co$incide",      &gl_coincide,     IS_INT,  GEZERO  },
	{ "Cr$essie",       &gl_cressie,      IS_INT,  GEZERO  },
	{ "a$lpha",         &gl_alpha,        IS_REAL, GEZERO  },
	{ "b$eta",          &gl_beta,         IS_REAL, GEZERO  },
	{ "c$utoff",        &gl_cutoff,       IS_REAL, GTZERO  },
	{ "de$bug",         &debug_level,     IS_INT,  GEZERO  },
	{ "display",        &gl_display,      IS_STRING, NOLIMIT },
	{ "do$ts",          &gl_dots,         IS_INT,  GEZERO  },
	{ "fit",            &gl_fit,          IS_INT,  GEZERO  },
	{ "fit_l$imit",     &gl_fit_limit,    IS_REAL, GTZERO  },
	{ "fo$rmat",        &gl_format,       IS_STRING, NOLIMIT },
	{ "fr$action",      &gl_fraction,     IS_REAL, GTZERO  },
	{ "gcv",            &gl_gcv,          IS_REAL, GTZERO  },
	{ "gls$_residuals", &gl_gls_residuals, IS_INT, GEZERO  },
	{ "gnuplot",        &gl_gnuplot,      IS_STRING, NOLIMIT },
	{ "gnuplot35",      &gl_gnuplot35,    IS_STRING, NOLIMIT },
	{ "gpt$erm",        &gl_gpterm,       IS_STRING, NOLIMIT  },
	{ "id$p",           &gl_idp,          IS_REAL, GEZERO  },
	{ "in$tervals",     &gl_n_intervals,  IS_INT,  GTZERO  },
	{ "it$er",          &gl_iter,         IS_INT,  GEZERO  },
	{ "j$graph",        &gl_jgraph,       IS_INT,  GEZERO  },
	{ "lhs",            &gl_lhs,          IS_INT,  GEZERO  },
	{ "log$file",       &logfile_name,    IS_STRING, NOLIMIT },
	{ "mvbeta",   		&gl_mvbeta,       IS_INT, GEZERO },
	{ "mv$string",		&gl_mv_string,    IS_STRING, NOLIMIT },
	{ "n_uk",           &gl_n_uk,         IS_INT,  GEZERO  },
	{ "numbers",        &gl_numbers,      IS_INT,  GEZERO  },
	{ "nb$lockdiscr",   &gl_nblockdiscr,  IS_INT,  GTZERO  },
	{ "no$check",       &gl_nocheck,      IS_INT,  GEZERO  },
	{ "ns$im",          &gl_nsim,         IS_INT,  GTZERO  },
	{ "o$utputfile",    &o_filename,      IS_STRING, NOLIMIT },
	{ "or$der",         &gl_order,        IS_INT,  GEZERO },
	{ "pag$er",         &gl_pager,        IS_STRING, NOLIMIT },
	{ "pl$otfile",      &gl_plotfile,     IS_STRING, NOLIMIT },
	{ "q$uantile",      &gl_quantile,     IS_REAL, GTZERO  },
	{ "rp",             &gl_rp,           IS_INT,  GEZERO  },
	{ "sec$ure",        &gl_secure,       IS_INT,  GTZERO  },
	{ "see$d",          &gl_seed,         IS_INT,  GTZERO  },
	{ "useed",          &gl_seed,         IS_UINT,  GEZERO  },
	{ "spa$rse",        &gl_sparse,       IS_INT,  GEZERO  },
	{ "spi$ral",        &gl_spiral,       IS_INT,  GEZERO  },
	{ "spl$it",         &gl_split,        IS_INT,  GTZERO  },
	{ "sy$mmetric",     &gl_sym_ev,       IS_INT,  GEZERO  },
	{ "tol_h$or",       &gl_tol_hor,      IS_REAL, GEZERO  },
	{ "tol_v$er",       &gl_tol_ver,      IS_REAL, GEZERO  },
	{ "v$erbose",       &debug_level,     IS_INT,  GEZERO  },
	{ "w$idth",         &gl_iwidth,       IS_REAL, GEZERO  },
	{ "x$valid",        &gl_xvalid,       IS_INT,  GEZERO  },
	{ "zero_di$st",     &gl_zero_est,     IS_INT,  GEZERO  },
	{ "zero",           &gl_zero,         IS_REAL, GEZERO  },
	{ "zm$ap",          &gl_zmap,         IS_REAL, NOLIMIT },
	{ "plotw$eights",   &gl_plotweights,  IS_INT, GEZERO   },
	{ NULL, NULL, 0, 0 }
	};

	for (i = 0; set_options[i].name; i++)
		if (almost_equals(name, set_options[i].name))
			break; /* break out i-loop */
	if (set_options[i].name == NULL)
		return 0;

	if (almost_equals((const char *)name,"nb$lockdiscr"))
		gl_gauss = 0; /* side effect */

	expr->name = name;
	expr->ptr = set_options[i].ptr;
	expr->what = set_options[i].what;
	expr->limit = set_options[i].limit;

	return 1;
}

static void check_assign_expr(GSTAT_EXPR *expr) {
/* for INT and REAL expressions, check range */
	double val;

	switch(expr->what) {
		case IS_INT: 
			val = (double) (*((int *)(expr->ptr)));
			break;
		case IS_REAL: 
			val = (*((double *)(expr->ptr)));
			break;
		default:
			return;
	}
	if (expr->limit == GEZERO && val < 0.0) {
		lex_error();
		pr_warning("value should be non-negative");
		ErrMsg(ER_IMPOSVAL, expr->name);
	}
	if (expr->limit == GTZERO && val <= 0.0) {
		lex_error();
		pr_warning("value should be positive");
		ErrMsg(ER_IMPOSVAL, expr->name);
	}
}

static void push_data_X(DATA *d, int id) {
	if (id == -1) { /* remove default intercept */
		if (d->n_X > 1) {
			lex_error();
			ErrMsg(ER_SYNTAX, "-1 only as first argument following X="); 
		}
		d->n_X = 0;
	} else if (id == 0) {
		lex_error();
		ErrMsg(ER_SYNTAX, "intercept is default"); 
	} else /* id > 0 */
		data_add_X(d, id);
}

static int is_block_expr(GSTAT_EXPR *expr, const char *s) {
	DPOINT *bp;

	bp = get_block_p();
	expr->name = s;
	expr->limit = GEZERO;
	expr->what = IS_REAL;
	if (almost_equals(s, "dx"))
		expr->ptr = &(bp->x);
	else if (almost_equals(s, "dy"))
		expr->ptr = &(bp->y);
	else if (almost_equals(s, "dz"))
		expr->ptr = &(bp->z);
	else
		return 0;
	return 1;
}

static void push_marginal(char *name, double value) {
	static int names = -1;

	if (names == -1)
		names = (name != NULL);

	if (name) {
		if (!names) {
			lex_error();
			ErrMsg(ER_SYNTAX, "only real values allowed"); 
		}
		gl_marginal_names = (char **) erealloc(gl_marginal_names,
			++gl_n_marginals * sizeof(char *));
		gl_marginal_names[gl_n_marginals - 1] = name;
	} else {
		if (names) {
			lex_error();
			ErrMsg(ER_SYNTAX, "only quoted strings allowed"); 
		}
		gl_marginal_values = (double *) erealloc (gl_marginal_values,
			++gl_n_marginals * sizeof(double));
		gl_marginal_values[gl_n_marginals - 1] = value;
	}
	return;
}

static void reset_parser(void) {
/* savety first: reset all static globals (should be unnessesary) */
	v = NULL;
	d = NULL;
	bp = NULL;
	ofn = NULL;
	expr.ptr = NULL;
	expr.what =  UNKNOWN;
	expr.limit =  NOLIMIT;
	id = id1 = id2 = col1 = col2 = -1;
}

int parse_cmd(const char *cmd, const char *fname) {
	set_lex_source(cmd, fname);
	reset_parser();
	return yyparse();
}

int parse_file(const char *fname) {
/* 
 * parse commands in file fname
 */
	int stdin_isatty = 1;
	char *cp;

	if (fname == NULL || strcmp(fname, "-") == 0) {
#ifdef HAVE_UNISTD_H
		stdin_isatty = isatty(fileno(stdin));
#endif
		if (stdin_isatty)
			cp = string_prompt("gstat> ");
		else
			cp = string_file(NULL);
	} else /* read from file */
		cp = string_file(fname);

	if (parse_cmd(cp, fname))
		ErrMsg(ER_SYNTAX, fname);
	efree(cp);

	if (boundary_file != NULL) {
		cp = string_file(boundary_file);
		if (parse_cmd(cp, boundary_file))
			ErrMsg(ER_SYNTAX, boundary_file);
		efree(cp);
	}

	if (vector_only && !allow_vector_only)
		ErrMsg(ER_SYNTAX, fname);

	return 0;
}

void parse_gstatrc(void) {
	char *fname = NULL, *cp;

	if ((fname = getenv(GSTATRC)) != NULL) {
		if (! file_exists(fname)) {
			message("environment variable %s:\n", GSTATRC);
			ErrMsg(ER_READ, fname);
		}
		parse_file(fname);
	} else if ((cp = getenv("HOME")) != NULL) {
		fname = (char *) emalloc(strlen(cp) + strlen(HOMERCFILE) + 2);
		sprintf(fname, "%s/%s", cp, HOMERCFILE);
		if (file_exists(fname))
			parse_file(fname);
		efree(fname);
	}
	return;
}

int read_variogram(VARIOGRAM *v, const char *source) {
	char *cp;
	int rval;

	parse_variogram = v;
	cp = (char *) emalloc((strlen(source) + 20) * sizeof(char));
	sprintf(cp, "variogram(): %s;", source);
	rval = parse_cmd(cp, NULL);
	parse_variogram = NULL; /* for savety */
	efree(cp);
	return rval;
}

int read_vector(D_VECTOR *d, char *fname) {
	int rval;

	assert(d != NULL);
	sd_vector = d;

	allow_vector_only = 1;

	rval = parse_file(fname);

	if (! vector_only)  {
		message("stat: only numeric input allowed -- \n");
		ErrMsg(ER_IMPOSVAL, fname);
	}

	return rval;
}

static void verify_data(DATA *d) { /* declaration : contents */

	if (d->var_fn_str != NULL) {
		if (almost_equals(d->var_fn_str, "mu"))
			d->variance_fn = v_mu;
		else if (almost_equals(d->var_fn_str, "mu(1-mu)"))
			d->variance_fn = v_bin;
		else if (almost_equals(d->var_fn_str, "identity"))
			d->variance_fn = v_identity;
		else {
			lex_error();
			message("variance function %s not supported:\n", d->var_fn_str);
			ErrMsg(ER_SYNTAX, d->var_fn_str);
		}
	}
}
