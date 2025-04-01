
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: c3d_global.h,v 1.14 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __C3D_GLOBAL_H_
#define __C3D_GLOBAL_H_

#include <stdio.h>
#include <math.h>

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

#define PI 3.14159265358979

#define CUBE_OF(X)   ((X)*(X)*(X))
#define ABS(X)       (((X) < 0. ) ? (-(X)) : (X))
#define SIGN(A)      (((A) < 0. ) ?    -1. : 1.)
#define DOT(A,B) ((A)[X]*(B)[X] + (A)[Y]*(B)[Y] + (A)[Z]*(B)[Z])
/* recursive allowed? MAG = sqrt(DOT) */
#define MAGNITUDE(A) (sqrt((A)[X]*(A)[X] + (A)[Y]*(A)[Y] + (A)[Z]*(A)[Z]))
#define LIN_MAGNITUDE(A,B,C) (((B)[X]*(C)[X] + (B)[Y]*(C)[Y] + (B)[Z]*(C)[Z])/(A))
#define SQUARE(A) ((A)*(A))
#define SUB( A, B, C )     { \
  (C)[X] =  (A)[X] - (B)[X]; \
  (C)[Y] =  (A)[Y] - (B)[Y]; \
  (C)[Z] =  (A)[Z] - (B)[Z]; \
   }
#define ADDVEC( A, B, C ) {  \
  (C)[X] =  (A)[X] + (B)[X]; \
  (C)[Y] =  (A)[Y] + (B)[Y]; \
  (C)[Z] =  (A)[Z] + (B)[Z]; \
   }
#define SCALEVEC( A, B, C ){ \
  (C)[X] =  (A) * (B)[X];    \
  (C)[Y] =  (A) * (B)[Y];    \
  (C)[Z] =  (A) * (B)[Z];    \
   }
#define LIN_SCALEVEC( A, DA, B, DB, C ){ \
  (C)[X] =  (A) * (DB)[X] + (DA) * (B)[X]; \
  (C)[Y] =  (A) * (DB)[Y] + (DA) * (B)[Y]; \
  (C)[Z] =  (A) * (DB)[Z] + (DA) * (B)[Z]; \
   }
#define CROSSVEC( A, B, C ) {                 \
  (C)[X] = (A)[Y] * (B)[Z] - (A)[Z] * (B)[Y]; \
  (C)[Y] = (A)[Z] * (B)[X] - (A)[X] * (B)[Z]; \
  (C)[Z] = (A)[X] * (B)[Y] - (A)[Y] * (B)[X]; \
   }
#define LIN_CROSSVEC( A, DA, B, DB, C ) {     \
  (C)[X] = ((A)[Y] * (DB)[Z] + (DA)[Y] * (B)[Z]) - ((A)[Z] * (DB)[Y] + (DA)[Z] * (B)[Y]); \
  (C)[Y] = ((A)[Z] * (DB)[X] + (DA)[Z] * (B)[X]) - ((A)[X] * (DB)[Z] + (DA)[X] * (B)[Z]); \
  (C)[Z] = ((A)[X] * (DB)[Y] + (DA)[X] * (B)[Y]) - ((A)[Y] * (DB)[X] + (DA)[Y] * (B)[X]); \
   }

#define STRING_LEN   511
#define FILENAME_LEN 256

#define MACHINE_EPSILON     1.E-14
#define SINGLE_EPS           3.e-7       /*  (~ 1./(TWO_TO_THE(22))).  */
#define REAL_INFINITY        1.E12

/* unstructured mesh constants */
#define NO_CELL_FLAG_INDX     -1

/*                    ----- other flags            -----  */
#define UNSET      -888888
#define BAD_SHORT    65535
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


/* ----- ASSERT from S.Maguire ISBN:1-55615-551-4 -----
 *       the wrapper  prevents the compiler from
 *       complaining if you include it in several places
 */
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/* _Assert is defined in memory_util.c */
extern void _Assert(const char *, const unsigned);

#ifdef DEBUG
# define ASSERT(f) if (!(f)) _Assert(__FILE__, __LINE__)
#else
# define ASSERT(f)
#endif

#ifdef __cplusplus
}
#endif

#endif /* __GLOBAL_H_ */
