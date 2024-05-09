/*
 ************************************************************************
 *                                                                      *
 * udpTester2 -- simple test udp                                        *
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

#define NUMUDPARGS 5
#include "udpUtilities.h"
#include "common.h"

/* shorthands for accessing argument values and velocities */
#define CG(      IUDP,I)   ((double *) (udps[IUDP].arg[0].val))[I]
#define CG_DOT(  IUDP,I)   ((double *) (udps[IUDP].arg[0].dot))[I]
#define RAD(     IUDP)     ((double *) (udps[IUDP].arg[1].val))[0]
#define RAD_DOT( IUDP)     ((double *) (udps[IUDP].arg[1].dot))[0]
#define SCALE(   IUDP)     ((double *) (udps[IUDP].arg[2].val))[0]
#define VOL(     IUDP)     ((double *) (udps[IUDP].arg[3].val))[0]
#define VOL_DOT( IUDP)     ((double *) (udps[IUDP].arg[3].dot))[0]
#define BBOX(    IUDP,I,J) ((double *) (udps[IUDP].arg[4].val))[3*(I)+(J)]
#define BBOX_DOT(IUDP,I,J) ((double *) (udps[IUDP].arg[4].dot))[3*(I)+(J)]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"cg",        "rad",       "scale",  "vol",        "bbox",       };
static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRREAL, -ATTRREALSEN, -ATTRREALSEN, };
static int    argIdefs[NUMUDPARGS] = {0,           0,           0,        0,            0,            };
static double argDdefs[NUMUDPARGS] = {0.,          1.,          1.,       0.,           0.,           };

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
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    double  data[18];
    char    *message=NULL;
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

    /* check arguments */
    if (udps[0].arg[0].size == 1 && CG(0,0) == 0) {
        udps[0].arg[0].size = 3;
        udps[0].arg[0].nrow = 3;
        udps[0].arg[0].ncol = 1;

        RALLOC(udps[0].arg[0].val, double, udps[0].arg[0].size);
        RALLOC(udps[0].arg[0].dot, double, udps[0].arg[0].size);

        CG(0,0) = 0;
        CG(0,1) = 0;
        CG(0,2) = 0;

        CG_DOT(0,0) = 0;
        CG_DOT(0,1) = 0;
        CG_DOT(0,2) = 0;
    }

    if (udps[0].arg[0].size != 3) {
        snprintf(message, 100, "if \"cg\" is specified, it should have 3 values");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "\"rad\" should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (RAD(0) <= 0) {
        snprintf(message, 100, "\"rad\" (=%f) should be positive", RAD(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "\"scale\" should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (SCALE(0) <= 0) {
        snprintf(message, 100, "\"scale\" (=%f) should be positive", SCALE(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    }

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("cg(     0) = %f %f %f\n", CG(0,0), CG(0,1), CG(0,2));
    printf("rad(    0) = %f\n", RAD(    0));
    printf("rad_dot(0) = %f\n", RAD_DOT(0));
    printf("scale(  0) = %f\n", SCALE(  0));
#endif

    /* make enough room for the bounding box */
    udps[0].arg[4].size = 6;
    udps[0].arg[4].nrow = 2;
    udps[0].arg[4].ncol = 3;

    RALLOC(udps[0].arg[4].val, double, udps[0].arg[4].size);
    RALLOC(udps[0].arg[4].dot, double, udps[0].arg[4].size);

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("cg(     %d) = %f %f %f\n", numUdp, CG(numUdp,0), CG(numUdp,1), CG(numUdp,2));
    printf("rad(    %d) = %f\n", numUdp, RAD(    numUdp));
    printf("rad_dot(%d) = %f\n", numUdp, RAD_DOT(numUdp));
    printf("scale(  %d) = %f\n", numUdp, SCALE(  numUdp));
#endif

    /* generate the centered box */
    data[0] = CG(0,0) - RAD(0) * SCALE(0);
    data[1] = CG(0,1) - RAD(0) * SCALE(0);
    data[2] = CG(0,2) - RAD(0) * SCALE(0);
    data[3] =       2 * RAD(0) * SCALE(0);
    data[4] =       2 * RAD(0) * SCALE(0);
    data[5] =       2 * RAD(0) * SCALE(0);

    status = EG_makeSolidBody(context, BOX, data, ebody);
    CHECK_STATUS(EG_makeSolidBody);

    assert(*ebody != NULL);

    status = EG_getMassProperties(*ebody, data);
    CHECK_STATUS(EG_getMassProperties);

    /* set the output value(s) associated with the box */
    VOL(     numUdp    ) = pow(2*RAD(numUdp)*SCALE(numUdp), 3);
    BBOX(    numUdp,0,0) = CG(numUdp,0) - RAD(numUdp) * SCALE(numUdp);
    BBOX(    numUdp,0,1) = CG(numUdp,1) - RAD(numUdp) * SCALE(numUdp);
    BBOX(    numUdp,0,2) = CG(numUdp,2) - RAD(numUdp) * SCALE(numUdp);
    BBOX(    numUdp,1,0) = CG(numUdp,0) + RAD(numUdp) * SCALE(numUdp);
    BBOX(    numUdp,1,1) = CG(numUdp,1) + RAD(numUdp) * SCALE(numUdp);
    BBOX(    numUdp,1,2) = CG(numUdp,2) + RAD(numUdp) * SCALE(numUdp);

    VOL_DOT( numUdp    ) = 0;
    BBOX_DOT(numUdp,0,0) = 0;
    BBOX_DOT(numUdp,0,1) = 0;
    BBOX_DOT(numUdp,0,2) = 0;
    BBOX_DOT(numUdp,1,0) = 0;
    BBOX_DOT(numUdp,1,1) = 0;
    BBOX_DOT(numUdp,1,2) = 0;

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
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp, ipnt, nnode, nedge, nface;
    int    oclass, mtype, nchild, *senses;
    double data[18];
    ego    *enodes, *eedges, *efaces, eent, eref, *echilds;

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

    /* find the ego entity */
    if (entType == OCSM_NODE) {
        status = EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        eent = enodes[entIndex-1];

        EG_free(enodes);
    } else if (entType == OCSM_EDGE) {
        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
        CHECK_STATUS(EG_getBodyTopos);

        eent = eedges[entIndex-1];

        EG_free(eedges);
    } else if (entType == OCSM_FACE) {
        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        eent = efaces[entIndex-1];

        EG_free(efaces);
    } else {
        printf("udpSensitivity: bad entType=%d\n", entType);
        status = EGADS_ATTRERR;
        goto cleanup;
    }

    /* loop through the points */
    for (ipnt = 0; ipnt < npnt; ipnt++) {

        /* find the physical coordinates */
        if        (entType == OCSM_NODE) {
            status = EG_getTopology(eent, &eref, &oclass, &mtype,
                                    data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);
        } else if (entType == OCSM_EDGE) {
            status = EG_evaluate(eent, &(uvs[ipnt]), data);
            CHECK_STATUS(EG_evaluate);
        } else if (entType == OCSM_FACE) {
            status = EG_evaluate(eent, &(uvs[2*ipnt]), data);
            CHECK_STATUS(EG_evaluate);
        }

        /* compute the sensitivity */
        vels[3*ipnt  ] = CG_DOT(iudp,0) + RAD_DOT(iudp) * (data[0]-CG(iudp,0)) / RAD(iudp);
        vels[3*ipnt+1] = CG_DOT(iudp,1) + RAD_DOT(iudp) * (data[1]-CG(iudp,1)) / RAD(iudp);
        vels[3*ipnt+2] = CG_DOT(iudp,2) + RAD_DOT(iudp) * (data[2]-CG(iudp,2)) / RAD(iudp);
    }

    /* compute the sensitivity of the volume */
    /*@-evalorder@*/
    VOL_DOT(iudp) = 24 * pow(SCALE(iudp),3) * pow(RAD(iudp),2) * RAD_DOT(iudp);
    /*@+evalorder@*/

    /* compute the sensitivity of the bounding box */
    BBOX_DOT(iudp,0,0) = CG_DOT(iudp,0) - RAD_DOT(iudp) * SCALE(iudp);
    BBOX_DOT(iudp,0,1) = CG_DOT(iudp,1) - RAD_DOT(iudp) * SCALE(iudp);
    BBOX_DOT(iudp,0,2) = CG_DOT(iudp,2) - RAD_DOT(iudp) * SCALE(iudp);
    BBOX_DOT(iudp,1,0) = CG_DOT(iudp,0) + RAD_DOT(iudp) * SCALE(iudp);
    BBOX_DOT(iudp,1,1) = CG_DOT(iudp,1) + RAD_DOT(iudp) * SCALE(iudp);
    BBOX_DOT(iudp,1,2) = CG_DOT(iudp,2) + RAD_DOT(iudp) * SCALE(iudp);

    status = EGADS_SUCCESS;

cleanup:
    return status;
}
