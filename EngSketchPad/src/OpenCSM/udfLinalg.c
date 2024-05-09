/*
 ************************************************************************
 *                                                                      *
 * udpLinalg -- perform linear algeraic calculations                    *
 *                                                                      *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                      John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2024  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS       4
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"
#include "common.h"

/* shorthands for accessing argument values and velocities */
#define OPER(IUDP    ) ((char   *) (udps[IUDP].arg[0].val))
#define M1(  IUDP,I,J) ((double *) (udps[IUDP].arg[1].val))[(I)*udps[IUDP].arg[1].ncol+(J)]
#define M2(  IUDP,I,J) ((double *) (udps[IUDP].arg[2].val))[(I)*udps[IUDP].arg[2].ncol+(J)]
#define ANS( IUDP,I,J) ((double *) (udps[IUDP].arg[3].val))[(I)*udps[IUDP].arg[3].ncol+(J)]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"oper",     "m1",     "m2",     "ans",     };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRREAL, ATTRREAL, -ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,        0,        0,         };
static double argDdefs[NUMUDPARGS] = {0.,         0.,       0.,       0.,        };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
#include <assert.h>

static int matsol(double A[], double b[], int n, int m, double x[]);


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  emodel,                 /* (in)  Model containing Body */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     i, j, k;
    int     oclass, mtype, nchild, *senses;
    double  data[18];
    char    *message=NULL;
    ego     context, eref, *ebodys;
    void    *realloc_temp = NULL;            /* used by RALLOC macro */
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 100, "expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("M1=\n");
    for (i = 0; i < udps[0].arg[1].nrow; i++) {
        for (j = 0; j < udps[0].arg[1].ncol; j++) {
            printf(" %10.5f", M1(0,i,j));
        }
        printf("\n");
    }
    printf("M2=\n");
    for (i = 0; i < udps[0].arg[2].nrow; i++) {
        for (j = 0; j < udps[0].arg[2].ncol; j++) {
            printf(" %10.5f", M2(0,i,j));
        }
        printf("\n");
    }
#endif

    /* check arguments (and perform matrix operations) */
    if        (strcmp(OPER(0), "add") == 0 || strcmp(OPER(0), "ADD") == 0) {
        /* scalar addition:  M1(scalar) + M2(matrix) */
        if (udps[0].arg[1].nrow == 1 && udps[0].arg[1].ncol == 1) {
            udps[0].arg[3].size = udps[0].arg[2].size;
            udps[0].arg[3].nrow = udps[0].arg[2].nrow;
            udps[0].arg[3].ncol = udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            for (i = 0; i < udps[0].arg[3].nrow; i++) {
                for (j = 0; j < udps[0].arg[3].ncol; j++) {
                    ANS(0,i,j) = M1(0,0,0) + M2(0,i,j);
                }
            }

        /* matrix addition:  M1(matrix) + M2(matrix) */
        } else if (udps[0].arg[1].nrow == udps[0].arg[2].nrow &&
                   udps[0].arg[1].ncol == udps[0].arg[2].ncol   ) {
            udps[0].arg[3].size = udps[0].arg[2].size;
            udps[0].arg[3].nrow = udps[0].arg[2].nrow;
            udps[0].arg[3].ncol = udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            for (i = 0; i < udps[0].arg[3].nrow; i++) {
                for (j = 0; j < udps[0].arg[3].ncol; j++) {
                    ANS(0,i,j) = M1(0,i,j) + M2(0,i,j);
                }
            }

        } else {
            snprintf(message, 100, "M1 is not scalar nor are M1 and M2 the same shape");
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }

    } else if (strcmp(OPER(0), "sub") == 0 || strcmp(OPER(0), "SUB") == 0) {
        /* scalar subtraction:  M1(scalar) - M2(matrix)  */
        if (udps[0].arg[1].nrow == 1 && udps[0].arg[1].ncol == 1) {
            udps[0].arg[3].size = udps[0].arg[2].size;
            udps[0].arg[3].nrow = udps[0].arg[2].nrow;
            udps[0].arg[3].ncol = udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            for (i = 0; i < udps[0].arg[3].nrow; i++) {
                for (j = 0; j < udps[0].arg[3].ncol; j++) {
                    ANS(0,i,j) = M1(0,0,0) - M2(0,i,j);
                }
            }

        /* matrix subrtaction:  M1(matrix) - M2(matrix) */
        } else if (udps[0].arg[1].nrow == udps[0].arg[2].nrow &&
                   udps[0].arg[1].ncol == udps[0].arg[2].ncol   ) {
            udps[0].arg[3].size = udps[0].arg[2].size;
            udps[0].arg[3].nrow = udps[0].arg[2].nrow;
            udps[0].arg[3].ncol = udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            for (i = 0; i < udps[0].arg[3].nrow; i++) {
                for (j = 0; j < udps[0].arg[3].ncol; j++) {
                    ANS(0,i,j) = M1(0,i,j) - M2(0,i,j);
                }
            }

        } else {
            snprintf(message, 100, "M1 is not scalar nor are M1 and M2 the same shape");
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }

    } else if (strcmp(OPER(0), "mult") == 0 || strcmp(OPER(0), "MULT") == 0) {
        /* scalar multiplication:  M1(scalar) * M2(matrix)  */
        if (udps[0].arg[1].nrow == 1 && udps[0].arg[1].ncol == 1) {
            udps[0].arg[3].size = udps[0].arg[2].size;
            udps[0].arg[3].nrow = udps[0].arg[2].nrow;
            udps[0].arg[3].ncol = udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            for (i = 0; i < udps[0].arg[3].nrow; i++) {
                for (j = 0; j < udps[0].arg[3].ncol; j++) {
                    ANS(0,i,j) = M1(0,0,0) * M2(0,i,j);
                }
            }

        /* matrix multiplication:  M1(matrix) * M2(matrix) */
        } else if (udps[0].arg[1].ncol == udps[0].arg[2].nrow) {
            udps[0].arg[3].size = udps[0].arg[1].nrow * udps[0].arg[2].ncol;
            udps[0].arg[3].nrow = udps[0].arg[1].nrow;
            udps[0].arg[3].ncol =                       udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            for (i = 0; i < udps[0].arg[1].nrow; i++) {
                for (j = 0; j < udps[0].arg[2].ncol; j++) {
                    ANS(0,i,j) = 0;
                    for (k = 0; k < udps[0].arg[1].ncol; k++) {
                        ANS(0,i,j) += M1(0,i,k) * M2(0,k,j);
                    }
                }
            }

        } else {
            snprintf(message, 100, "M1 is not scalar nor are M1 and M2 compatable shapes");
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }


    } else if (strcmp(OPER(0), "div"  ) == 0 || strcmp(OPER(0), "DIV"  ) == 0 ||
               strcmp(OPER(0), "solve") == 0 || strcmp(OPER(0), "SOLVE") == 0) {
        /* matrix division (matrix solve):  M1inv(matrix) * M2(matrix) */
        if (udps[0].arg[1].nrow == udps[0].arg[1].ncol &&
            udps[0].arg[1].nrow == udps[0].arg[2].nrow   ) {
            udps[0].arg[3].size = udps[0].arg[1].ncol * udps[0].arg[2].ncol;
            udps[0].arg[3].nrow = udps[0].arg[1].ncol;
            udps[0].arg[3].ncol =                       udps[0].arg[2].ncol;

            RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

            status = matsol(&(M1(0,0,0)), &(M2(0,0,0)), udps[0].arg[1].nrow, udps[0].arg[2].ncol, &(ANS(0,0,0)));
            CHECK_STATUS(matsol);

        } else {
            snprintf(message, 100, "M1 and M2 compatable shapes");
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }

    } else if (strcmp(OPER(0), "trans") == 0 || strcmp(OPER(0), "TRANS") == 0) {
        /* matrix transpose:  M1trans(any) */
        udps[0].arg[3].size = udps[0].arg[1].size;
        udps[0].arg[3].nrow = udps[0].arg[1].ncol;
        udps[0].arg[3].ncol = udps[0].arg[1].nrow;

        RALLOC(udps[0].arg[3].val, double, udps[0].arg[3].size);

        for (i = 0; i < udps[0].arg[3].nrow; i++) {
            for (j = 0; j < udps[0].arg[3].ncol; j++) {
                ANS(0,i,j) = M1(0,j,i);
            }
        }

    } else {
        snprintf(message, 100, "OPER must be ADD, SUB, MULT, DIV, SOLVE, or TRANS");
        status = OCSM_UDP_ERROR2;
        goto cleanup;
    }

#ifdef DEBUG
    printf("ANS=\n");
    for (i = 0; i < udps[0].arg[3].nrow; i++) {
        for (j = 0; j < udps[0].arg[3].ncol; j++) {
            printf(" %10.5f", ANS(0,i,j));
        }
        printf("\n");
    }
#endif

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

    /* make a copy of the Body (so that it does not get removed
     when OpenCSM deletes emodel) */
    status = EG_copyObject(ebodys[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* add a special Attribute to the Body to tell OpenCSM that there
       is no topological change and hence it should not adjust the
       Attributes on the Body in finishBody() */
    status = EG_attributeAdd(*ebody, "__noTopoChange__", ATTRSTRING,
                             0, NULL, NULL, "udfEditAttr");
    CHECK_STATUS(EG_attributeAdd);

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

    if (strlen(message) > 0) {
        *string = message;
        printf("%s\n", message);
    } else if (status != EGADS_SUCCESS) {
        FREE(message);
        *string = udpErrorStr(status);
    } else {
        FREE(message);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                      *
 ************************************************************************
 */

int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
   /*@unused@*/int    npnt,             /* (in)  number of points */
   /*@unused@*/int    entType,          /* (in)  OCSM entity type */
   /*@unused@*/int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
   /*@unused@*/double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d)\n",
           (long long)ebody, npnt, entType, entIndex);
#endif

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        status = EGADS_NOTMODEL;
        goto cleanup;
    }


    /* this routine is not written yet */
    status = EGADS_NOLOAD;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   matsol - Gaussian elimination with partial pivoting                *
 *                                                                      *
 ************************************************************************
 */

static int
matsol(double    A[],                   /* (in)  matrix to be solved (stored rowwise) */
                                        /* (out) upper-triangular form of matrix */
       double    b[],                   /* (in)  right hand side */
                                        /* (out) right-hand side after swapping */
       int       n,                     /* (in)  size of matrix */
       int       m,                     /* (in)  number of right-hand sides */
       double    x[])                   /* (out) solution of A*x=b */
{
    int       status = SUCCESS;         /* (out) return status */

    int       ir, jc, kc, imax;
    double    amax, swap, fact;

    ROUTINE(matsol);

    /* --------------------------------------------------------------- */

    /* reduce each column of A */
    for (kc = 0; kc < n; kc++) {

        /* find pivot element */
        imax = kc;
        amax = fabs(A[kc*n+kc]);

        for (ir = kc+1; ir < n; ir++) {
            if (fabs(A[ir*n+kc]) > amax) {
                imax = ir;
                amax = fabs(A[ir*n+kc]);
            }
        }

        /* check for possibly-singular matrix (ie, near-zero pivot) */
        if (amax < EPS12) {
            status = OCSM_UDP_ERROR3;
            goto cleanup;
        }

        /* if diagonal is not pivot, swap rows in A and b */
        if (imax != kc) {
            for (jc = 0; jc < n; jc++) {
                swap         = A[kc  *n+jc];
                A[kc  *n+jc] = A[imax*n+jc];
                A[imax*n+jc] = swap;
            }

            for (jc = 0; jc < m; jc++) {
                swap         = b[kc  *m+jc];
                b[kc  *m+jc] = b[imax*m+jc];
                b[imax*m+jc] = swap;
            }
        }

        /* row-reduce part of matrix to the bottom of and right of [kc,kc] */
        for (ir = kc+1; ir < n; ir++) {
            fact = A[ir*n+kc] / A[kc*n+kc];

            for (jc = kc+1; jc < n; jc++) {
                A[ir*n+jc] -= fact * A[kc*n+jc];
            }

            for (jc = 0; jc < m; jc++) {
                b[ir*m+jc] -= fact * b[kc*m+jc];
            }

            A[ir*n+kc] = 0;
        }
    }

    /* back-substitution pass */
    for (jc = 0; jc < m; jc++) {
        x[(n-1)*m+jc] = b[(n-1)*m+jc] / A[(n-1)*n+(n-1)];

        for (ir = n-2; ir >= 0; ir--) {
            x[ir*m+jc] = b[ir*m+jc];
            for (kc = ir+1; kc < n; kc++) {
                x[ir*m+jc] -= A[ir*n+kc] * x[kc*m+jc];
            }
            x[ir*m+jc] /= A[ir*n+ir];
        }
    }

cleanup:
    return status;
}
