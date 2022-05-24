/*
 ************************************************************************
 *                                                                      *
 * udpStag -- udp file for a simple turbomachinery airfoil generator    *
 *                                                                      *
 *            Written by Bridget Dixon @ Syracuse University            *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                      John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2022  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 9
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define RAD1(     IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
#define RAD1_DOT( IUDP)  ((double *) (udps[IUDP].arg[0].dot))[0]
#define BETA1(    IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define BETA1_DOT(IUDP)  ((double *) (udps[IUDP].arg[1].dot))[0]
#define GAMA1(    IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define GAMA1_DOT(IUDP)  ((double *) (udps[IUDP].arg[2].dot))[0]
#define RAD2(     IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define RAD2_DOT( IUDP)  ((double *) (udps[IUDP].arg[3].dot))[0]
#define BETA2(    IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define BETA2_DOT(IUDP)  ((double *) (udps[IUDP].arg[4].dot))[0]
#define GAMA2(    IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]
#define GAMA2_DOT(IUDP)  ((double *) (udps[IUDP].arg[5].dot))[0]
#define ALFA(     IUDP)  ((double *) (udps[IUDP].arg[6].val))[0]
#define ALFA_DOT( IUDP)  ((double *) (udps[IUDP].arg[6].val))[0]
#define XFRNT(    IUDP)  ((double *) (udps[IUDP].arg[7].val))[0]
#define XFRNT_DOT(IUDP)  ((double *) (udps[IUDP].arg[7].val))[0]
#define XREAR(    IUDP)  ((double *) (udps[IUDP].arg[8].val))[0]
#define XREAR_DOT(IUDP)  ((double *) (udps[IUDP].arg[8].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"rad1",   "beta1",  "gama1",
                                      "rad2",   "beta2",  "gama2",
                                      "alfa",   "xfrnt",  "xrear",  };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRREAL,
                                      ATTRREAL, ATTRREAL, ATTRREAL,
                                      ATTRREAL, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,        0,        0,
                                      0,        0,        0,
                                      0,        0,        0,        };
static double argDdefs[NUMUDPARGS] = {0.10,     30.0,     10.0,
                                      0.05,    -40.0,      5.0,
                                      -30.0,    0.333,    0.667     };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define   PI        3.1415926535897931159979635
#define   TWOPI     6.2831853071795862319959269
#define   PIo2      1.5707963267948965579989817
#define   PIo4      0.7853981633974482789994909
#define   PIo180    0.0174532925199432954743717


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

    int     header[3], senses[4];
    double  x[7], y[7], theta3, theta4, theta5, theta6;
    double  xa, ya, xb, yb, xc, yc, xd, yd, data[18];
    ego     enodes[5], ecurve, eedges[8], eloop, eface;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */
    
#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("rad1(0)       = %f\n", RAD1(     0));
    printf("rad1_dot(0)   = %f\n", RAD1_DOT( 0));
    printf("beta1(0)      = %f\n", BETA1(    0));
    printf("beta1_dot(0)  = %f\n", BETA1_DOT(0));
    printf("gama1(0)      = %f\n", GAMA1(    0));
    printf("gama1_dot(0)  = %f\n", GAMA1_DOT(0));
    printf("rad2(0)       = %f\n", RAD2(     0));
    printf("rad2_dot(0)   = %f\n", RAD2_DOT( 0));
    printf("beta2(0)      = %f\n", BETA2(    0));
    printf("beta2_dot(0)  = %f\n", BETA2_DOT(0));
    printf("gama2(0)      = %f\n", GAMA2(    0));
    printf("gama2_dot(0)  = %f\n", GAMA2_DOT(0));
    printf("alfa(0)       = %f\n", ALFA(     0));
    printf("alfa_dot(0)   = %f\n", ALFA_DOT( 0));
    printf("xfrnt(0)      = %f\n", XFRNT(    0));
    printf("xfrnt_dot(0)  = %f\n", XFRNT_DOT(0));
    printf("xrear(0)      = %f\n", XREAR(    0));
    printf("xrear_dot(0)  = %f\n", XREAR_DOT(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        printf(" udpExecute: rad1 should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (RAD1(0) <= 0) {
        printf(" udpExecute: rad1 should be positive\n");
        status  = EGADS_RANGERR;
        goto cleanup;
        

    } else if (udps[0].arg[1].size > 1) {
        printf(" udpExecute: beta1 should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        printf(" udpExecute: gama1 should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        printf(" udpExecute: rad2 should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (RAD2(0) <= 0) {
        printf(" udpExecute: rad2 should be positive\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size > 1) {
        printf(" udpExecute: beta2 should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        printf(" udpExecute: gama2 should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[6].size > 1) {
        printf(" udpExecute: alfa should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[7].size > 1) {
        printf(" udpExecute: xfrnt should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (XFRNT(0) <= 0 || XFRNT(0) >= XREAR(0)) {
        printf(" udpExecute: xfrnt should be between 0 and xrear\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[8].size > 1) {
        printf(" udpExecute: xrear should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (XREAR(0) <= XFRNT(0) || XREAR(0) >= 1) {
        printf(" udpExecute: xrear should be between xfrnt and 1\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("rad1(%d)       = %f\n", numUdp, RAD1(     numUdp));
    printf("rad1_dot(%d)   = %f\n", numUdp, RAD1_DOT( numUdp));
    printf("beta1(%d)      = %f\n", numUdp, BETA1(    numUdp));
    printf("beta1_dot(%d)  = %f\n", numUdp, BETA1_DOT(numUdp));
    printf("gama1(%d)      = %f\n", numUdp, GAMA1(    numUdp));
    printf("gama1_dot(%d)  = %f\n", numUdp, GAMA1_DOT(numUdp));
    printf("rad2(%d)       = %f\n", numUdp, RAD2(     numUdp));
    printf("rad2_dot(%d)   = %f\n", numUdp, RAD2_DOT( numUdp));
    printf("beta2(%d)      = %f\n", numUdp, BETA2(    numUdp));
    printf("beta2_dot(%d)  = %f\n", numUdp, BETA2_DOT(numUdp));
    printf("gama2(%d)      = %f\n", numUdp, GAMA2(    numUdp));
    printf("gama2_dot(%d)  = %f\n", numUdp, GAMA2_DOT(numUdp));
    printf("alfa(%d)       = %f\n", numUdp, ALFA(     numUdp));
    printf("alfa_dot(%d)   = %f\n", numUdp, ALFA_DOT( numUdp));
    printf("xfrnt(%d)      = %f\n", numUdp, XFRNT(    numUdp));
    printf("xfrnt_dot(%d)  = %f\n", numUdp, XFRNT_DOT(numUdp));
    printf("xrear(%d)      = %f\n", numUdp, XREAR(    numUdp));
    printf("xrear_dot(%d)  = %f\n", numUdp, XREAR_DOT(numUdp));
#endif

    /* centers of leading and trailing edge circles */
    x[1] = RAD1(numUdp);
    y[1] = 0.0;

    x[2] = 1.0 - RAD2(numUdp);
    y[2] = y[1] + (x[2] - x[1]) * tan(ALFA(numUdp)*PIo180);

    /* leading edge tangency points */
    theta3 = (BETA1(numUdp) + 90 + GAMA1(numUdp)) * PIo180;
    x[3]   = x[1] + RAD1(numUdp) * cos(theta3);
    y[3]   = y[1] + RAD1(numUdp) * sin(theta3);
    theta6 = (BETA1(numUdp) - 90 - GAMA1(numUdp)) * PIo180;
    x[6]   = x[1] + RAD1(numUdp) * cos(theta6);
    y[6]   = y[1] + RAD1(numUdp) * sin(theta6);

    /* trailing edge tangency points */
    theta4 = (BETA2(numUdp) + 90 - GAMA2(numUdp)) * PIo180;
    x[4]   = x[2] + RAD2(numUdp) * cos(theta4);
    y[4]   = y[2] + RAD2(numUdp) * sin(theta4);
    theta5 = (BETA2(numUdp) - 90 + GAMA2(numUdp)) * PIo180;
    x[5]   = x[2] + RAD2(numUdp) * cos(theta5);
    y[5]   = y[2] + RAD2(numUdp) * sin(theta5);

    /* upper bezier */
    xa = XFRNT(numUdp);
    ya = y[3] + (xa - x[3]) * tan(theta3+PIo2);

    xb = XREAR(numUdp);
    yb = y[4] + (xb - x[4]) * tan(theta4+PIo2);

    /* lower bezier */
    xc = XFRNT(numUdp);
    yc = y[6] + (xc - x[6]) * tan(theta6+PIo2);

    xd = XREAR(numUdp);
    yd = y[5] + (xd - x[5]) * tan(theta5+PIo2);

    /* make Nodes (starting at upper trailing edge) */
    data[0] = x[4];
    data[1] = y[4];
    data[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 9, data, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    data[0] = x[3];
    data[1] = y[3];
    data[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 9, data, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    data[0] = x[6];
    data[1] = y[6];
    data[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 9, data, 0, NULL, NULL, &(enodes[2]));
    CHECK_STATUS(EG_makeTopology);

    data[0] = x[5];
    data[1] = y[5];
    data[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 9, data, 0, NULL, NULL, &(enodes[3]));
    CHECK_STATUS(EG_makeTopology);

    enodes[4] = enodes[0];

    /* make the Edges (starting at upper surface) */
    header[0] = 0;
    header[1] = 3;
    header[2] = 4;
    data[  0] = x[4];    data[  1] = y[4];    data[  2] = 0;
    data[  3] = xb;      data[  4] = yb;      data[  5] = 0;
    data[  6] = xa;      data[  7] = ya;      data[  8] = 0;
    data[  9] = x[3];    data[ 10] = y[3];    data[ 11] = 0;

    status = EG_makeGeometry(context, CURVE, BEZIER, NULL, header, data, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    data[0] = 0;
    data[1] = 1;
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, data, 2, &(enodes[0]), NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

    data[ 0] = x[1];           data[ 1] = y[1];           data[ 2] = 0;
    data[ 3] = x[3] - x[1];    data[ 4] = y[3] - y[1];    data[ 5] = 0;
    data[ 6] = -data[4];       data[ 7] = +data[3];       data[ 8] = 0;
    data[ 9] = RAD1(numUdp);

    status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    data[0] = 0;
    data[1] = PI - 2 * GAMA1(numUdp) * PIo180;
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, data, 2, &(enodes[1]), NULL, &(eedges[1]));
    CHECK_STATUS(EG_makeTopology);

    header[0] = 0;
    header[1] = 3;
    header[2] = 4;
    data[  0] = x[6];    data[  1] = y[6];    data[  2] = 0;
    data[  3] = xc;      data[  4] = yc;      data[  5] = 0;
    data[  6] = xd;      data[  7] = yd;      data[  8] = 0;
    data[  9] = x[5];    data[ 10] = y[5];    data[ 11] = 0;

    status = EG_makeGeometry(context, CURVE, BEZIER, NULL, header, data, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    data[0] = 0;
    data[1] = 1;
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, data, 2, &(enodes[2]), NULL, &(eedges[2]));
    CHECK_STATUS(EG_makeTopology);

    data[ 0] = x[2];           data[ 1] = y[2];           data[ 2] = 0;
    data[ 3] = x[5] - x[2];    data[ 4] = y[5] - y[2];    data[ 5] = 0;
    data[ 6] = -data[4];       data[ 7] = +data[3];       data[ 8] = 0;
    data[ 9] = RAD2(numUdp);

    status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    data[0] = 0;
    data[1] = PI - 2 * GAMA2(numUdp) * PIo180;
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, data, 2, &(enodes[3]), NULL, &(eedges[3]));
    CHECK_STATUS(EG_makeTopology);


    /* make Loop from the Edges */
    senses[0] = SFORWARD;
    senses[1] = SFORWARD;
    senses[2] = SFORWARD;
    senses[3] = SFORWARD;
    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, eedges, senses, &eloop);
    CHECK_STATUS(EG_makeTopology);

    /* make Face from the Loop */
    status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
    CHECK_STATUS(EG_makeFace);

    /* create the FaceBody (which will be returned) */
    status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* set the output value(s) */
    status = EG_getMassProperties(*ebody, data);
    CHECK_STATUS(EG_getMassProperties);

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
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
    int    iudp, judp;

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
