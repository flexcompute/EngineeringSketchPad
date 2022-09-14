/*
 ************************************************************************
 *                                                                      *
 * common.h -- common header file                                       *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2022  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

#ifndef _COMMON_H_
#define _COMMON_H_

/* macros for error checking */
#define CHECK_STATUS(X)                                                 \
    if (status < SUCCESS) {                                             \
        printf( "ERROR:: BAD STATUS = %d from %s (called from %s:%d)\n", status, #X, routine, __LINE__); \
        goto cleanup;                                                   \
    }
#define SET_STATUS(STAT,X)                                              \
    status = STAT;                                                      \
    printf( "ERROR:: BAD STATUS = %d from %s (called from %s:%d)\n", status, #X, routine, __LINE__); \
    goto cleanup;
#define MALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR != NULL) {                                                  \
        printf("ERROR:: MALLOC overwrites for %s=%llx (called from %s:%d)\n", #PTR, (long long)PTR, routine, __LINE__); \
        status = BAD_MALLOC;                                            \
        goto cleanup;                                                   \
    }                                                                   \
    PTR = (TYPE *) EG_alloc((SIZE) * sizeof(TYPE));                     \
    if (PTR == NULL) {                                                  \
        printf("ERROR:: MALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
        status = BAD_MALLOC;                                            \
        goto cleanup;                                                   \
    }
#define RALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR == NULL) {                                                  \
        MALLOC(PTR,TYPE,SIZE);                                          \
    } else {                                                            \
       realloc_temp = EG_reall(PTR, (SIZE) * sizeof(TYPE));            \
       if (PTR == NULL) {                                               \
           printf("ERROR:: RALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
           status = BAD_MALLOC;                                         \
           goto cleanup;                                                \
       } else {                                                         \
           PTR = (TYPE *)realloc_temp;                                  \
       }                                                                \
    }
#define FREE(PTR)                                               \
    if (PTR != NULL) {                                          \
        EG_free(PTR);                                           \
    }                                                           \
    PTR = NULL;
#define STRLEN(A)   (int)strlen(A)

#define ROUTINE(NAME) char routine[] = #NAME ;\
    if (routine[0] == '\0') printf("bad routine(%s)\n", routine);

/* macros for status printing */
#define SPRINT0(OL,FORMAT)                                   \
    if (outLevel >= OL) {                                    \
        printf(FORMAT); printf("\n");                        \
    }
#define SPRINT0x(OL,FORMAT)                                  \
    if (outLevel >= OL) {                                    \
        printf(FORMAT); fflush(stdout);                      \
    }
#define SPRINT1(OL,FORMAT,A)                                 \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A); printf("\n");                      \
    }
#define SPRINT1x(OL,FORMAT,A)                                \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A); fflush(stdout);                    \
    }
#define SPRINT2(OL,FORMAT,A,B)                               \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B); printf("\n");                    \
    }
#define SPRINT2x(OL,FORMAT,A,B)                              \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B); fflush(stdout);                  \
    }
#define SPRINT3(OL,FORMAT,A,B,C)                             \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C); printf("\n");                  \
    }
#define SPRINT3x(OL,FORMAT,A,B,C)                            \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C); fflush(stdout);                \
    }
#define SPRINT4(OL,FORMAT,A,B,C,D)                           \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D); printf("\n");                \
    }
#define SPRINT4x(OL,FORMAT,A,B,C,D)                          \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D); fflush(stdout);              \
    }
#define SPRINT5(OL,FORMAT,A,B,C,D,E)                         \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E); printf("\n");              \
    }
#define SPRINT6(OL,FORMAT,A,B,C,D,E,F)                       \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F); printf("\n");            \
    }
#define SPRINT7(OL,FORMAT,A,B,C,D,E,F,G)                     \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G); printf("\n");          \
    }
#define SPRINT8(OL,FORMAT,A,B,C,D,E,F,G,H)                   \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H); printf("\n");        \
    }
#define SPRINT9(OL,FORMAT,A,B,C,D,E,F,G,H,I)                 \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I); printf("\n");      \
    }
#define SPRINT10(OL,FORMAT,A,B,C,D,E,F,G,H,I,J)              \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I,J); printf("\n");    \
    }
#define SPRINT11(OL,FORMAT,A,B,C,D,E,F,G,H,I,J,K)            \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I,J,K); printf("\n");  \
    }
#define SPRINT12(OL,FORMAT,A,B,C,D,E,F,G,H,I,J,K,L)          \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I,J,K,L); printf("\n"); \
    }
#define SPRINT13(OL,FORMAT,A,B,C,D,E,F,G,H,I,J,K,L,M)        \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I,J,K,L,M); printf("\n");  \
    }
#define SPRINT14(OL,FORMAT,A,B,C,D,E,F,G,H,I,J,K,L,M,N)      \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I,J,K,L,M,N); printf("\n");  \
    }
#define SPRINT15(OL,FORMAT,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O)    \
    if (outLevel >= OL) {                                    \
        printf(FORMAT,A,B,C,D,E,F,G,H,I,J,K,L,M,N,O); printf("\n");  \
    }

#define SPLINT_CHECK_FOR_NULL(X)                                        \
    if ((X) == NULL) {                                                  \
        printf("ERROR:: SPLINT found %s is NULL (called from %s:%d)\n", #X, routine, __LINE__); \
        status = -9999;                                                 \
        goto cleanup;                                                   \
    }

/* error codes */
#define           BAD_MALLOC      -900

/*  miscellaneous macros */
#define           PI              3.1415926535897931159979635
#define           TWOPI           6.2831853071795862319959269
#define           PIo2            1.5707963267948965579989817
#define           PIo4            0.7853981633974482789994909
#define           PIo180          0.0174532925199432954743717

#define           HUGEQ           99999999.0
#define           HUGEI           9999999
#define           EPS03           1.0e-03
#define           EPS06           1.0e-06
#define           EPS09           1.0e-09
#define           EPS12           1.0e-12
#define           EPS20           1.0e-20
#define           SQR(A)          ((A) * (A))
#define           NINT(A)         (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define           MINMAX(A,B,C)   MIN(MAX((A), (B)), (C))
#define           SIGN(A)         (((A) < 0) ? -1 : (((A) > 0) ? +1 : 0))
#define           FSIGN(A,B)      ((B) >= 0 ? fabs(A) : -fabs(A))


#endif   /* _COMMON_H_ */
