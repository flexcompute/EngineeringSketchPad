/*
 ************************************************************************
 *                                                                      *
 * udpUtilities.h -- header for udp*.c                                  *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
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

#ifndef _UDPUTILITIES_H_
#define _UDPUTILITIES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "egads.h"
#include "udp.h"

//$$$#define DEBUG_UDP

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define PI       3.1415926535897931159979635
#define TWOPI    6.2831853071795862319959269
#define NINT(A)  (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))
#define SQR(A)   ((A)*(A))

#define ROUTINE(NAME) char routine[] = #NAME ;\
    if (routine[0] == '\0') printf("bad routine(%s)\n", routine);
#define CHECK_STATUS(X)                                                 \
    if (status < EGADS_SUCCESS) {                                       \
        printf( "ERROR:: BAD STATUS = %d from %s (called from %s:%d)\n", status, #X, routine, __LINE__); \
        goto cleanup;                                                   \
    }
#define SET_STATUS(STAT,X)                                              \
    status = STAT;                                                      \
    printf( "ERROR:: BAD STATUS = %d from %s (called from %s:%d)\n", status, #X, routine, __LINE__); \
    goto cleanup;
#define MALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR != NULL) {                                                  \
        printf("ERROR:: MALLOC overwrites for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
        status = EGADS_MALLOC;                                          \
        goto cleanup;                                                   \
    }                                                                   \
    PTR = (TYPE *) EG_alloc((SIZE) * sizeof(TYPE));                     \
    if (PTR == NULL) {                                                  \
        printf("ERROR:: MALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
        status = EGADS_MALLOC;                                          \
        goto cleanup;                                                   \
    }
#define RALLOC(PTR,TYPE,SIZE)                                           \
    realloc_temp = EG_reall(PTR, (SIZE) * sizeof(TYPE));                \
    if (PTR == NULL) {                                                  \
        printf("ERROR:: RALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
        status = EGADS_MALLOC;                                          \
        goto cleanup;                                                   \
    } else {                                                            \
        PTR = (TYPE *)realloc_temp;                                     \
    }
#define FREE(PTR)                                               \
    if (PTR != NULL) {                                          \
        EG_free(PTR);                                           \
    }                                                           \
    PTR = NULL;


#define SPLINT_CHECK_FOR_NULL(X)                                        \
    if ((X) == NULL) {                                                  \
        printf("ERROR:: SPLINT found %s is NULL (called from %s:%d)\n", #X, routine, __LINE__); \
        status = OCSM_UNKNOWN;                                          \
        goto cleanup;                                                   \
    }

/* definition of the udp structures */
typedef struct {
    void      *val;
    void      *dot;
    int       size;
} udparg_T;

typedef struct {
    ego        ebody;
    udparg_T   arg[NUMUDPARGS];
    int        *bodyList;
    void       *data;         /* private data */
} udp_T;

/* storage for UDPs */
static int    numUdp = 0;     /* number of UDPs */
static udp_T  *udps  = NULL;  /* array  of UDPs */

#endif  /* _UDPUTILITIES_H_ */

