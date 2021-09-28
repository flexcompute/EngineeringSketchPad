/*
 ************************************************************************
 *                                                                      *
 * udpBox -- udp file to generate a box                                 *
 *                                                                      *
 *            Written by Bridget Dixon @ Syracuse University            *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                      John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2021  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 6
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define DX(    IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
#define DX_DOT(IUDP)  ((double *) (udps[IUDP].arg[0].dot))[0]
#define DY(    IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define DY_DOT(IUDP)  ((double *) (udps[IUDP].arg[1].dot))[0]
#define DZ(    IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define DZ_DOT(IUDP)  ((double *) (udps[IUDP].arg[2].dot))[0]
#define RAD(   IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define AREA(  IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define VOLUME(IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"dx",        "dy",        "dz",        "rad",    "area",    "volume", };
static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRREALSEN, ATTRREAL, -ATTRREAL, -ATTRREAL,};
static int    argIdefs[NUMUDPARGS] = {0,           0,           0,           0,        0,         0,        };
static double argDdefs[NUMUDPARGS] = {0.,          0.,          0.,          0.,       0.,        0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define  EPS06   1.0e-6


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
    int     wire, sense[8], nedges, add=1;
    double  node1[3], node2[3], node3[3], node4[3], node5[3], node6[3], node7[3], node8[3];
    double  cent1[3], cent2[3], cent3[3], cent4[3], axis1[3], axis2[3];
    double  data[18], trange[2];
    char    *message=NULL;
    ego     enodes[9], ecurve, eedges[8], eloop, eface, etemp, *eedges2;

    ROUTINE(udpExecute);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("dx(0)      = %f\n", DX(    0));
    printf("dx_dot(0)  = %f\n", DX_DOT(0));
    printf("dy(0)      = %f\n", DY(    0));
    printf("dy_dot(0)  = %f\n", DY_DOT(0));
    printf("dz(0)      = %f\n", DZ(    0));
    printf("dz_dot(0)  = %f\n", DZ_DOT(0));
    printf("rad(0)     = %f\n", RAD(   0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "dx should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (DX(0) < 0) {
        snprintf(message, 100, "dx = %f < 0", DX(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "dy should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (DY(0) < 0) {
        snprintf(message, 100, "dy = %f < 0", DY(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "dz should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (DZ(0) < 0) {
        snprintf(message, 100, "dz = %f < 0", DZ(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        snprintf(message, 100, "rad should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (RAD(0) < 0) {
        snprintf(message, 100, "rad = %f < 0", RAD(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (DX(0) <= 0 && DY(0) <= 0 && DZ(0) <= 0) {
        snprintf(message, 100, "dx=dy=dz=0");
        status  = EGADS_GEOMERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("dx[%d]     = %f\n", numUdp, DX(    numUdp));
    printf("dx_dot[%d] = %f\n", numUdp, DX_DOT(numUdp));
    printf("dy[%d]     = %f\n", numUdp, DY(    numUdp));
    printf("dy_dot[%d] = %f\n", numUdp, DY_DOT(numUdp));
    printf("dz[%d]     = %f\n", numUdp, DZ(    numUdp));
    printf("dz_dot[%d] = %f\n", numUdp, DZ_DOT(numUdp));
    printf("rad[%d]    = %f\n", numUdp, RAD(   numUdp));
#endif

    /* check for 3D solid body (and make if requested) */
    if (DX(0) > 0 && DY(0) > 0 && DZ(0) > 0)  {

        /* make sure that radius is not too big */
        if (2*RAD(0) >= DX(0) || 2*RAD(0) >= DY(0) || 2*RAD(0) >= DZ(0)) {
            snprintf(message, 100, "radius must be less half of all side lengths");
            status = EGADS_GEOMERR;
            goto cleanup;
        }

        data[0] = -DX(0) / 2;
        data[1] = -DY(0) / 2;
        data[2] = -DZ(0) / 2;
        data[3] =  DX(0);
        data[4] =  DY(0);
        data[5] =  DZ(0);

        /* Make Solid Body */
        status = EG_makeSolidBody(context, BOX, data, ebody);
        CHECK_STATUS(EG_makeSolidBody);
        if (*ebody == NULL) goto cleanup;    // needed for splint

        /* apply rounded edges (if required) */
        if (RAD(0) > 0) {
            etemp = *ebody;

            status = EG_getBodyTopos(etemp, NULL, EDGE, &nedges, &eedges2);
            CHECK_STATUS(EG_getBodyTopos);

            status = EG_filletBody(etemp, nedges, eedges2, RAD(0), ebody, NULL);
            CHECK_STATUS(EG_filletBody);

            EG_free(eedges2);
        }

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;
    }

    /* otherwise this is either a 1D wirebody or a 2D facebody */
    wire = 0;            /* flag to tell if wire or face */

    /* geometry start */
    node1[0] = -DX(0) / 2;
    node1[1] = -DY(0) / 2;
    node1[2] = -DZ(0) / 2;

    /* wire body dz > 0 */
    if        (DX(0) == 0 && DY(0) == 0) {
        wire = 1;

        node2[0] =    0;
        node2[1] =    0;
        node2[2] = DZ(0) / 2;

    /* wire body dy > 0 */
    } else if (DX(0) == 0 && DZ(0) == 0) {
        wire = 1;

        node2[0] =    0;
        node2[1] = DY(0) / 2;
        node2[2] =    0;

    /* wire body dx > 0 */
    } else if (DY(0) == 0 && DZ(0) == 0) {
        wire = 1;

        node2[0] = DX(0) / 2;
        node2[1] =    0;
        node2[2] =    0;

    /* face in the x-y plane */
    } else if (DZ(0) == 0) {
        node2[0] = +DX(0) / 2;
        node2[1] = -DY(0) / 2;
        node2[2] =     0;

        node3[0] = +DX(0) / 2;
        node3[1] = +DY(0) / 2;
        node3[2] =     0;

        node4[0] = -DX(0) / 2;
        node4[1] = +DY(0) / 2;
        node4[2] =     0;

    /* face in the y-z plane */
    } else if (DX(0) == 0) {
        node2[0] =     0;
        node2[1] = +DY(0) / 2;
        node2[2] = -DZ(0) / 2;

        node3[0] =     0;
        node3[1] = +DY(0) / 2;
        node3[2] = +DZ(0) / 2;

        node4[0] =     0;
        node4[1] = -DY(0) / 2;
        node4[2] = +DZ(0) / 2;

    /* face in the x-z plane */
    } else {
        node2[0] = -DX(0) / 2;
        node2[1] =     0;
        node2[2] = +DZ(0) / 2;

        node3[0] = +DX(0) / 2;
        node3[1] =     0;
        node3[2] = +DZ(0) / 2;

        node4[0] = +DX(0) / 2;
        node4[1] =     0;
        node4[2] = -DZ(0) / 2;
    }

    /* make 1D wire body */
    if (wire == 1) {

        /* check to ensure that rad was not set */
        if (RAD(0) > 0) {
            snprintf(message, 100, "rad cannot be set for wirebody");
            status = EGADS_GEOMERR;
            goto cleanup;
        }

        /* make Nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 1 */
        data[0] = node1[0];
        data[1] = node1[1];
        data[2] = node1[2];
        data[3] = node2[0] - node1[0];
        data[4] = node2[1] - node1[1];
        data[5] = node2[2] - node1[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node1, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node2, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* make Loop from this Edge */
        sense[0] = SFORWARD;
        sense[1] = SFORWARD;
        status = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, 1, eedges, sense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        /* create the WireBody (which will be returned) */
        status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);
        if (*ebody == NULL) goto cleanup;

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;

    /* make 2D face body (without rounded corners) */
    } else if (RAD(0) == 0) {

        /*        y,z,x
                    ^
                    |
              4     |    3
                    |
                    +----------> x,y,z

              1----->    2
        */

        /* make Nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node3, 0, NULL, NULL, &(enodes[2]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node4, 0, NULL, NULL, &(enodes[3]));
        CHECK_STATUS(EG_makeTopology);

        /* first and last Nodes are the same */
        enodes[4] = enodes[0];

        /* make the Line 1 */
        data[0] = node1[0];
        data[1] = node1[1];
        data[2] = node1[2];
        data[3] = node2[0] - node1[0];
        data[4] = node2[1] - node1[1];
        data[5] = node2[2] - node1[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node1, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node2, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 2 */
        data[0] = node2[0];
        data[1] = node2[1];
        data[2] = node2[2];
        data[3] = node3[0] - node2[0];
        data[4] = node3[1] - node2[1];
        data[5] = node3[2] - node2[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node2, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node3, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 3 */
        data[0] = node3[0];
        data[1] = node3[1];
        data[2] = node3[2];
        data[3] = node4[0] - node3[0];
        data[4] = node4[1] - node3[1];
        data[5] = node4[2] - node3[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node3, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node4, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[2]), NULL, &(eedges[2]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 4 */
        data[0] = node4[0];
        data[1] = node4[1];
        data[2] = node4[2];
        data[3] = node1[0] - node4[0];
        data[4] = node1[1] - node4[1];
        data[5] = node1[2] - node4[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node4, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node1, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[3]), NULL, &(eedges[3]));
        CHECK_STATUS(EG_makeTopology);

        /* make Loop from this Edge */
        sense[0] = SFORWARD;
        sense[1] = SFORWARD;
        sense[2] = SFORWARD;
        sense[3] = SFORWARD;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, eedges, sense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        /* make Face from the loop */
        status = EG_makeFace(eloop, SREVERSE, NULL, &eface);
        CHECK_STATUS(EG_makeFace);

        /* since this will make a PLANE, we need to add an Attribute
           to tell OpenCSM to scale the UVs when computing sensitivities */
        status = EG_attributeAdd(eface, "_scaleuv", ATTRINT, 1,
                                 &add, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        /* create the FaceBody (which will be returned) */
        status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);
        if (*ebody == NULL) goto cleanup;    // needed for splint

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;

    /* make 2D face body (with rounded corners) */
    } else {

        /*        y,z,x
                    ^
                6   |    5
              7 c4  |   c3 4
                    |
                    +----------> x,y,z
              8 c1      c2 3
                1----->  2
        */

        /* make Nodes */
        if (DX(0) > 0 && DY(0) > 0 && DZ(0) == 0){

            /* make sure that radius is not too big */
            if (2*RAD(0) >= DX(0) || 2*RAD(0) >= DY(0)) {
                snprintf(message, 100, "radius cannot be greater than half of any side length");
                status = EGADS_GEOMERR;
                goto cleanup;
            }

            /* in x-y plane */
            node1[0] = -DX(0) / 2 + RAD(0);
            node1[1] = -DY(0) / 2;
            node1[2] =  0;

            node2[0] = +DX(0) / 2 - RAD(0);
            node2[1] = -DY(0) / 2;
            node2[2] =  0;

            node3[0] = +DX(0) / 2;
            node3[1] = -DY(0) / 2 + RAD(0);
            node3[2] =  0;

            node4[0] = +DX(0) / 2;
            node4[1] = +DY(0) / 2 - RAD(0);
            node4[2] =  0;

            node5[0] = +DX(0) / 2 - RAD(0);
            node5[1] = +DY(0) / 2;
            node5[2] =  0;

            node6[0] = -DX(0) / 2 + RAD(0);
            node6[1] = +DY(0) / 2;
            node6[2] =  0;

            node7[0] = -DX(0) / 2;
            node7[1] = +DY(0) / 2 - RAD(0);
            node7[2] =  0;

            node8[0] = -DX(0) / 2;
            node8[1] = -DY(0) / 2 + RAD(0);
            node8[2] =  0;

            cent1[0] = -DX(0) / 2 + RAD(0);
            cent1[1] = -DY(0) / 2 + RAD(0);
            cent1[2] =  0;

            cent2[0] = +DX(0) / 2 - RAD(0);
            cent2[1] = -DY(0) / 2 + RAD(0);
            cent2[2] =  0;

            cent3[0] = +DX(0) / 2 - RAD(0);
            cent3[1] = +DY(0) / 2 - RAD(0);
            cent3[2] =  0;

            cent4[0] = -DX(0) / 2 + RAD(0);
            cent4[1] = +DY(0) / 2 - RAD(0);
            cent4[2] =  0;

            axis1[0] = 1;
            axis1[1] = 0;
            axis1[2] = 0;

            axis2[0] = 0;
            axis2[1] = 1;
            axis2[2] = 0;

        } else if (DX(0) > 0 && DY(0) == 0 && DZ(0) > 0) {

            /* make sure that radius is not too big */
            if (2*RAD(0) >= DX(0) || 2*RAD(0) >= DZ(0)) {
                snprintf(message, 100, "radius cannot be greater than half of any side length");
                status = EGADS_GEOMERR;
                goto cleanup;
            }

            /* in z-x plane */
            node1[0] = -DX(0) / 2;
            node1[1] =  0;
            node1[2] = -DZ(0) / 2 + RAD(0);

            node2[0] = -DX(0) / 2;
            node2[1] =  0;
            node2[2] = +DZ(0) / 2 - RAD(0);

            node3[0] = -DX(0) / 2 + RAD(0);
            node3[1] =  0;
            node3[2] = +DZ(0) / 2;

            node4[0] = +DX(0) / 2 - RAD(0);
            node4[1] =  0;
            node4[2] = +DZ(0) / 2;

            node5[0] = +DX(0) / 2;
            node5[1] =  0;
            node5[2] = +DZ(0) / 2 - RAD(0);

            node6[0] = +DX(0) / 2;
            node6[1] =  0;
            node6[2] = -DZ(0) / 2 + RAD(0);

            node7[0] = +DX(0) / 2 - RAD(0);
            node7[1] =  0;
            node7[2] = -DZ(0) / 2;

            node8[0] = -DX(0) / 2 + RAD(0);
            node8[1] =  0;
            node8[2] = -DZ(0) / 2;

            cent1[0] = -DX(0) / 2 + RAD(0);
            cent1[1] =  0;
            cent1[2] = -DZ(0) / 2 + RAD(0);

            cent2[0] = -DX(0) / 2 + RAD(0);
            cent2[1] =  0;
            cent2[2] = +DZ(0) / 2 - RAD(0);

            cent3[0] = +DX(0) / 2 - RAD(0);
            cent3[1] =  0;
            cent3[2] = +DZ(0) / 2 - RAD(0);

            cent4[0] = +DX(0) / 2 - RAD(0);
            cent4[1] =  0;
            cent4[2] = -DZ(0) / 2 + RAD(0);

            axis1[0] = 0;
            axis1[1] = 0;
            axis1[2] = 1;

            axis2[0] = 1;
            axis2[1] = 0;
            axis2[2] = 0;

        } else {

            /* make sure that radius is not too big */
            if (2*RAD(0) >= DY(0) || 2*RAD(0) >= DZ(0)) {
                snprintf(message, 100, "radius cannot be greater than half of any side length");
                status = EGADS_GEOMERR;
                goto cleanup;
            }

            /* in y-z plane */
            node1[0] =  0;
            node1[1] = -DY(0) / 2 + RAD(0);
            node1[2] = -DZ(0) / 2;

            node2[0] =  0;
            node2[1] = +DY(0) / 2 - RAD(0);
            node2[2] = -DZ(0) / 2;

            node3[0] =  0;
            node3[1] = +DY(0) / 2;
            node3[2] = -DZ(0) / 2 + RAD(0);

            node4[0] =  0;
            node4[1] = +DY(0) / 2;
            node4[2] = +DZ(0) / 2 - RAD(0);

            node5[0] =  0;
            node5[1] = +DY(0) / 2 - RAD(0);
            node5[2] = +DZ(0) / 2;

            node6[0] =  0;
            node6[1] = -DY(0) / 2 + RAD(0);
            node6[2] = +DZ(0) / 2;

            node7[0] =  0;
            node7[1] = -DY(0) / 2;
            node7[2] = +DZ(0) / 2 - RAD(0);

            node8[0] =  0;
            node8[1] = -DY(0) / 2;
            node8[2] = -DZ(0) / 2 + RAD(0);

            cent1[0] =  0;
            cent1[1] = -DY(0) / 2 + RAD(0);
            cent1[2] = -DZ(0) / 2 + RAD(0);

            cent2[0] =  0;
            cent2[1] = +DY(0) / 2 - RAD(0);
            cent2[2] = -DZ(0) / 2 + RAD(0);

            cent3[0] =  0;
            cent3[1] = +DY(0) / 2 - RAD(0);
            cent3[2] = +DZ(0) / 2 - RAD(0);

            cent4[0] =  0;
            cent4[1] = -DY(0) / 2 + RAD(0);
            cent4[2] = +DZ(0) / 2 - RAD(0);

            axis1[0] = 0;
            axis1[1] = 1;
            axis1[2] = 0;

            axis2[0] = 0;
            axis2[1] = 0;
            axis2[2] = 1;
        }

        /* make the nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node3, 0, NULL, NULL, &(enodes[2]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node4, 0, NULL, NULL, &(enodes[3]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node5, 0, NULL, NULL, &(enodes[4]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node6, 0, NULL, NULL, &(enodes[5]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node7, 0, NULL, NULL, &(enodes[6]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node8, 0, NULL, NULL, &(enodes[7]));
        CHECK_STATUS(EG_makeTopology);

        /* first and last Nodes are the same */
        enodes[8] = enodes[0];

        /* make the Line 1 */
        data[0] = node1[0];
        data[1] = node1[1];
        data[2] = node1[2];
        data[3] = node2[0] - node1[0];
        data[4] = node2[1] - node1[1];
        data[5] = node2[2] - node1[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node1, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node2, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* make rounded corner 1 */
        data[0] = cent2[0];
        data[1] = cent2[1];
        data[2] = cent2[2];
        data[3] = axis1[0];
        data[4] = axis1[1];
        data[5] = axis1[2];
        data[6] = axis2[0];
        data[7] = axis2[1];
        data[8] = axis2[2];
        data[9] = RAD(0);

        status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node2, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node3, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[0] > trange[1]){
            trange[1] += TWOPI;
        }

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 2 */
        data[0] = node3[0];
        data[1] = node3[1];
        data[2] = node3[2];
        data[3] = node4[0] - node3[0];
        data[4] = node4[1] - node3[1];
        data[5] = node4[2] - node3[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node3, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node4, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[2]), NULL, &(eedges[2]));
        CHECK_STATUS(EG_makeTopology);

        /* make rounded corner 2 */
        data[0] = cent3[0];
        data[1] = cent3[1];
        data[2] = cent3[2];
        data[3] = axis1[0];
        data[4] = axis1[1];
        data[5] = axis1[2];
        data[6] = axis2[0];
        data[7] = axis2[1];
        data[8] = axis2[2];
        data[9] = RAD(0);

        status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node4, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node5, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[0] > trange[1]){
            trange[1] += TWOPI;
        }

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[3]), NULL, &(eedges[3]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 3 */
        data[0] = node5[0];
        data[1] = node5[1];
        data[2] = node5[2];
        data[3] = node6[0] - node5[0];
        data[4] = node6[1] - node5[1];
        data[5] = node6[2] - node5[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node5, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node6, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[4]), NULL, &(eedges[4]));
        CHECK_STATUS(EG_makeTopology);

        /* make rounded corner 3 */
        data[0] = cent4[0];
        data[1] = cent4[1];
        data[2] = cent4[2];
        data[3] = axis1[0];
        data[4] = axis1[1];
        data[5] = axis1[2];
        data[6] = axis2[0];
        data[7] = axis2[1];
        data[8] = axis2[2];
        data[9] = RAD(0);

        status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node6, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node7, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[0] > trange[1]){
            trange[1] += TWOPI;
        }

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[5]), NULL, &(eedges[5]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line 4 */
        data[0] = node7[0];
        data[1] = node7[1];
        data[2] = node7[2];
        data[3] = node8[0] - node7[0];
        data[4] = node8[1] - node7[1];
        data[5] = node8[2] - node7[2];
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node7, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node8, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[6]), NULL, &(eedges[6]));
        CHECK_STATUS(EG_makeTopology);

        /* make rounded corner 4 */
        data[0] = cent1[0];
        data[1] = cent1[1];
        data[2] = cent1[2];
        data[3] = axis1[0];
        data[4] = axis1[1];
        data[5] = axis1[2];
        data[6] = axis2[0];
        data[7] = axis2[1];
        data[8] = axis2[2];
        data[9] = RAD(0);

        status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* get the parameter range */
        status = EG_invEvaluate(ecurve, node8, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, node1, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[0] > trange[1]){
            trange[1] += TWOPI;
        }

        /* make edge */
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, &(enodes[7]), NULL, &(eedges[7]));
        CHECK_STATUS(EG_makeTopology);

        /* make Loop from this Edge */
        sense[0] = SFORWARD;
        sense[1] = SFORWARD;
        sense[2] = SFORWARD;
        sense[3] = SFORWARD;
        sense[4] = SFORWARD;
        sense[5] = SFORWARD;
        sense[6] = SFORWARD;
        sense[7] = SFORWARD;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 8, eedges, sense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        /* make Face from the loop */
        status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
        CHECK_STATUS(EG_makeFace);

        /* since this will make a PLANE, we need to add an Attribute
           to tell OpenCSM to scale the UVs when computing sensitivities */
        status = EG_attributeAdd(eface, "_scaleuv", ATTRINT, 1,
                                 &add, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        /* create the FaceBody (which will be returned) */
        status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);
        if (*ebody == NULL) goto cleanup;    // needed for splint

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;
    }

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

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
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp, ipnt, nnode, nedge, nface, nchild, oclass, mtype, *senses;
    double data[18];
    ego    eref, *echilds, *enodes, *eedges, *efaces, eent;

    ROUTINE(udpSensitivity);
    
#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f)\n",
           (long long)ebody, npnt, entType, entIndex, uvs[0]);
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

    /* cannot compute sensitivity for RAD */
    if (RAD(iudp) > 0) {
        printf("udpSensitivity: dx=%f, dy=%f, dz=%f, rad=%f\n",
               DX(iudp), DY(iudp), DZ(iudp), RAD(iudp));
        return EGADS_ATTRERR;
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
        if (fabs(DX(iudp)) > EPS06) {
            vels[3*ipnt  ] = data[0] / DX(iudp) * DX_DOT(iudp);
        } else {
            vels[3*ipnt  ] = 0;
        }
        if (fabs(DY(iudp)) > EPS06) {
            vels[3*ipnt+1] = data[1] / DY(iudp) * DY_DOT(iudp);
        } else {
            vels[3*ipnt+1] = 0;
        }
        if (fabs(DZ(iudp)) > EPS06) {
            vels[3*ipnt+2] = data[2] / DZ(iudp) * DZ_DOT(iudp);
        } else {
            vels[3*ipnt+2] = 0;
        }
    }

    status = EGADS_SUCCESS;

cleanup:
    return status;
}
