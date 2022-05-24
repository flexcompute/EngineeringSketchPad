/*
 ************************************************************************
 *                                                                      *
 * udfNuscale -- performs non-uniform scaling on a B-spline Body        *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
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

#define NUMUDPARGS       6
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define XSCALE(IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
#define YSCALE(IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define ZSCALE(IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define XCENT( IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define YCENT( IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define ZCENT( IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"xscale", "yscale", "zscale", "xcent",  "ycent",  "zcent",  };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {1,        1,        1,        0,        0,        0,        };
static double argDdefs[NUMUDPARGS] = {1.,       1.,       1.,       0.,       0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define  EPS06   1.0e-6

/* unpublished routine in OpenCSM.c */
extern int convertToBSplines(ego inbody, double mat[], ego *ebody);


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

    int     oclass, mtype, nchild, *senses, ndum;
    double  xyz[18], mat[12];
    ego     *ebodys, eref, *edum;

    char    *message=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("xscale(0) = %f\n", XSCALE(0));
    printf("yscale(0) = %f\n", YSCALE(0));
    printf("zscale(0) = %f\n", ZSCALE(0));
    printf("xcent( 0) = %f\n", XCENT( 0));
    printf("ycent( 0) = %f\n", YCENT( 0));
    printf("zcent( 0) = %f\n", ZCENT( 0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            xyz, &nchild, &ebodys, &senses);
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

    status = EG_getTopology(ebodys[0], &eref, &oclass, &mtype,
                            xyz, &ndum, &edum, &senses);
    CHECK_STATUS(EG_getTopology);

    /* check arguments */
    if        (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "xscale should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "yscale should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "zscale should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        snprintf(message, 100, "xcent should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size > 1) {
        snprintf(message, 100, "ycent should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        snprintf(message, 100, "zcent should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (XSCALE(0) < EPS06 && mtype == SOLIDBODY) {
        snprintf(message, 100, "xscale must be positive for a SolidBody");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (YSCALE(0) < EPS06 && mtype == SOLIDBODY) {
        snprintf(message, 100, "yscale must be positive for a SolidBody");;
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (ZSCALE(0) < EPS06 && mtype == SOLIDBODY) {
        snprintf(message, 100, "zscale must be positive for a SOlidBody");;
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (fabs(XSCALE(0)) < EPS06 && fabs(YSCALE(0)) < EPS06) {
        snprintf(message, 100, "xscale and yscale cannot both be 0");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (fabs(YSCALE(0)) < EPS06 && fabs(ZSCALE(0)) < EPS06) {
        snprintf(message, 100, "yscale and zscale cannot both be 0");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (fabs(ZSCALE(0)) < EPS06 && fabs(XSCALE(0)) < EPS06) {
        snprintf(message, 100, "zscale and xscale cannot both be 0");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("xscale(%d) = %f\n", numUdp, XSCALE(numUdp));
    printf("yscale(%d) = %f\n", numUdp, YSCALE(numUdp));
    printf("zscale(%d) = %f\n", numUdp, ZSCALE(numUdp));
    printf("xcent( %d) = %f\n", numUdp, XCENT( numUdp));
    printf("ycent( %d) = %f\n", numUdp, YCENT( numUdp));
    printf("zcent( %d) = %f\n", numUdp, ZCENT( numUdp));
#endif

    /* set up transformation matrix */
    mat[ 0] = XSCALE(0);   mat[ 1] = 0;           mat[ 2] = 0;           mat[ 3] = XCENT(0) * (1-XSCALE(0));
    mat[ 4] = 0;           mat[ 5] = YSCALE(0);   mat[ 6] = 0;           mat[ 7] = YCENT(0) * (1-YSCALE(0));
    mat[ 8] = 0;           mat[ 9] = 0;           mat[10] = ZSCALE(0);   mat[11] = ZCENT(0) * (1-ZSCALE(0));

    status = convertToBSplines(ebodys[0], mat, ebody);
    CHECK_STATUS(convertToBSplines);

    /* set the output value */

    /* the Body is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
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
    int iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    /* this routine is not written yet */
    return EGADS_NOLOAD;
}

