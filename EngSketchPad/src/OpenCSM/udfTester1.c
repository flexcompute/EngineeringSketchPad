/*
 ************************************************************************
 *                                                                      *
 * udpTester1 -- simple test udp                                        *
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

#define NUMUDPARGS       3
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"
#include "common.h"

/* shorthands for accessing argument values and velocities */
#define AMAT(    IUDP,I,J) ((double *) (udps[IUDP].arg[0].val))[(I)*udps[IUDP].arg[0].ncol+(J)]
#define BMAT(    IUDP,I,J) ((double *) (udps[IUDP].arg[1].val))[(I)*udps[IUDP].arg[1].ncol+(J)]
#define BMAT_DOT(IUDP,I,J) ((double *) (udps[IUDP].arg[1].dot))[(I)*udps[IUDP].arg[1].ncol+(J)]
#define CMAT(    IUDP,I,J) ((double *) (udps[IUDP].arg[2].val))[(I)*udps[IUDP].arg[2].ncol+(J)]
#define CMAT_DOT(IUDP,I,J) ((double *) (udps[IUDP].arg[2].dot))[(I)*udps[IUDP].arg[2].ncol+(J)]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"amat",   "bmat",      "cmat",       };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREALSEN, -ATTRREALSEN, };
static int    argIdefs[NUMUDPARGS] = {0,        0,           0,            };
static double argDdefs[NUMUDPARGS] = {0.,       0.,          0.,           };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
#include <assert.h>


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

    int     irow, icol, k;
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

    /* check arguments */
    if (udps[0].arg[0].ncol != udps[0].arg[1].nrow) {
        snprintf(message, 100, "amat.ncol (%d) != mat.nrow (%d)", udps[0].arg[0].ncol, udps[0].arg[1].nrow);
        status = OCSM_UDP_ERROR1;
        goto cleanup;
    }

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    for (irow = 0; irow < udps[0].arg[0].nrow; irow++) {
        for (icol = 0; icol < udps[0].arg[0].ncol; icol++) {
            printf("amat(     0,%d,%d) = %f\n", irow, icol, AMAT(    0,irow,icol));
        }
    }
    for (irow = 0; irow < udps[0].arg[1].nrow; irow++) {
        for (icol = 0; icol < udps[0].arg[1].ncol; icol++) {
            printf("bmat(     0,%d,%d) = %f\n", irow, icol, BMAT(    0,irow,icol));
            printf("bmat_dot( 0,%d,%d) = %f\n", irow, icol, BMAT_DOT(0,irow,icol));
        }
    }
#endif

    /* make enough room for the matrix product */
    udps[0].arg[2].size = udps[0].arg[0].nrow * udps[0].arg[1].ncol;
    udps[0].arg[2].nrow = udps[0].arg[0].nrow;
    udps[0].arg[2].ncol =                       udps[0].arg[1].ncol;

    RALLOC(udps[0].arg[2].val, double, udps[0].arg[2].size);
    RALLOC(udps[0].arg[2].dot, double, udps[0].arg[2].size);

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    for (irow = 0; irow < udps[numUdp].arg[0].nrow; irow++) {
        for (icol = 0; icol < udps[numUdp].arg[0].ncol; icol++) {
            printf("amat(     %d,%d,%d) = %f\n", numUdp, irow, icol, AMAT(    numUdp,irow,icol));
        }
    }
    for (irow = 0; irow < udps[numUdp].arg[1].nrow; irow++) {
        for (icol = 0; icol < udps[numUdp].arg[1].ncol; icol++) {
            printf("bmat(     %d,%d,%d) = %f\n", numUdp, irow, icol, BMAT(    numUdp,irow,icol));
            printf("bmat_dot( %d,%d,%d) = %f\n", numUdp, irow, icol, BMAT_DOT(numUdp,irow,icol));
        }
    }
#endif

    /* perform the matrix multiplication and set the output values */
#ifndef __clang_analyzer__
    for (irow = 0; irow < udps[numUdp].arg[0].nrow; irow++) {
        for (icol = 0; icol < udps[numUdp].arg[1].ncol; icol++) {
            CMAT(    numUdp,irow,icol) = 0;
            CMAT_DOT(numUdp,irow,icol) = 0;
            for (k = 0; k < udps[0].arg[0].ncol; k++) {
                CMAT(numUdp,irow,icol) += AMAT(numUdp,irow,k) * BMAT(numUdp,k,icol);
            }
        }
    }
#endif

#ifdef DEBUG
    for (irow = 0; irow < udps[numUdp].arg[2].nrow; irow++) {
        for (icol = 0; icol < udps[numUdp].arg[2].ncol; icol++) {
            printf("C[%2d,%2d]=%12.6f\n", irow, icol, CMAT(numUdp,irow, icol));
        }
    }
#endif

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
    int    irow, icol, k;

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

    /* compute the sensitivity of the matrix product */
    for (irow = 0; irow < udps[iudp].arg[0].nrow; irow++) {
        for (icol = 0; icol < udps[iudp].arg[1].ncol; icol++) {
            CMAT_DOT(iudp,irow,icol) = 0;
            for (k = 0; k < udps[iudp].arg[0].ncol; k++) {
                CMAT_DOT(iudp,irow,icol) += AMAT(iudp,irow,k) * BMAT_DOT(iudp,k,icol);
            }
        }
    }

    status = EGADS_SUCCESS;

cleanup:
    return status;
}
