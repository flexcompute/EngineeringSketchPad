/*
 * $Id: c3d_global.h,v 1.10 2012/04/30 20:16:36 maftosmi Exp $
 *
 */
#ifndef __C3D_GLOBAL_H_
#define __C3D_GLOBAL_H_

/* ------------ system include files --------- */
#include <stdio.h>
#include <math.h>
/* ------------------------------------------- */

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MAX3(a,b,c) ((((a)>(b))&&((a)>(c))) ? (a) : (((b)>(c)) ? (b) : (c)))
#define MIN3(a,b,c) ((((a)<(b))&&((a)<(c))) ? (a) : (((b)<(c)) ? (b) : (c)))

#define TWO_TO_THE(EXP)  (1L << (EXP))                    /* type is long */
#define ONE_MB 1048576
#define ONE_KB 1024

#define IS_NAN(A)    ((A) != (A)) 

#define DIM 3

#define X  0                   /*      ...dimension names */
#define Y  1
#define Z  2


#define CUBE_OF(X)   ((X)*(X)*(X))
#define ABS(X)       (((X) < 0. ) ? (-(X)) : (X))
#define SIGN(A)      (((A) < 0. ) ?    -1. : 1.)
#define DOT(A,B) ((A)[X]*(B)[X] + (A)[Y]*(B)[Y] + (A)[Z]*(B)[Z])
/* recursive allowed? MAG = sqrt(DOT) */
#define MAGNITUDE(A) (sqrt((A)[X]*(A)[X] + (A)[Y]*(A)[Y] + (A)[Z]*(A)[Z]))
#define SQUARE(A) ((A)*(A))
#define SUB( A, B, C ) { \
  (C)[X] =  (A)[X] - (B)[X]; \
  (C)[Y] =  (A)[Y] - (B)[Y]; \
  (C)[Z] =  (A)[Z] - (B)[Z]; \
   }
#define ADDVEC( A, B, C ) { \
  (C)[X] =  (A)[X] + (B)[X]; \
  (C)[Y] =  (A)[Y] + (B)[Y]; \
  (C)[Z] =  (A)[Z] + (B)[Z]; \
   }

#define STRING_LEN 511
#define FILENAME_LEN 256


#define MACHINE_EPSILON     1.E-14
#define SINGLE_EPS           3.e-7       /*  (~ 1./(TWO_TO_THE(22))).  */
#define REAL_INFINITY         1.E12     

/* unstructured mesh constants */
#define NO_CELL_FLAG_INDX     -1


/*                    ----- other flags            -----  */
#define UNSET                 -888888
#define BAD_SHORT              65535
#define BAD_INDEX      -17

/*                    ----- define Error Codes      ----  */
#define FILE_ERROR      -1
#define PARSE_ERROR     -3
#define ASSERT_ERROR    -5

#define ERR  printf(" ===> ERROR:  ");printf /* standardize IO msgs */
#define WARN printf(" ===> WARNING:");printf
#define NOTE printf("\r    o  "     );printf
#define CONT printf("\r     . "     );printf
#define ATTN printf(" ===> ATTENTION: ");printf
#define INFO printf("# INFO: "  );    printf


/*                    ----- ASSERT from S.Maguire ISBN:1-55615-551-4 ----- */
/*                          the wrapper  prevents the compiler from
			    complaining if you include it in several places
			                                         -M. Aftosmis*/
#include <assert.h>

/* _Assert is defined in memory_util.c */
extern void _Assert(const char *, const unsigned);

#ifdef DEBUG
#
#  define ASSERT(f)                \
      if (!(f)) _Assert(__FILE__, __LINE__)
#else
#     define ASSERT(f)   
#endif

#endif /* __GLOBAL_H_ */
