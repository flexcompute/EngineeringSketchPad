/*
 ************************************************************************
 *                                                                      *
 * udpSample -- sample UDP                                              *
 *                                                                      *
 *              this makes a box, plate, or wire (centered at origin)   *
 *              and returns its area and volume (for demo purposes)     *
 *                                                                      *
 *              sensitivities can be computed for DX, DY, and/or DZ     *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
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

/* the number of "input" Bodys

   this only needs to be specified if this is a UDF (user-defined
   function) that consumes Bodys from OpenCSM's stack.  (the default
   value is 0).  the Bodys are passed into udpExecute via its
   emodel orgument */
#define NUMUDPINPUTBODYS 0

/*  the name of the routine to clean up private data

    this only needs to be specified if the UDP hangs private
    data onto udps[iudp].data (which is a void*), and that data
    could not be freed by a simple call to EG_free().  cases where
    this happens is when udps[iudp].data points to a structure that
    contains components that had be allocated via separate calls
    to EG_alloc (or malloc) or which used an allocator other
    than EG_alloc */
#define FREEUDPDATA(A) freePrivateData(A)
static int freePrivateData(void *data);

/* the number of arguments (specified below) */
#define NUMUDPARGS 6

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define DX(        IUDP  )  ((double *) (udps[IUDP].arg[0].val))[0]
#define DX_DOT(    IUDP  )  ((double *) (udps[IUDP].arg[0].dot))[0]
#define DY(        IUDP  )  ((double *) (udps[IUDP].arg[1].val))[0]
#define DY_DOT(    IUDP  )  ((double *) (udps[IUDP].arg[1].dot))[0]
#define DZ(        IUDP  )  ((double *) (udps[IUDP].arg[2].val))[0]
#define DZ_DOT(    IUDP  )  ((double *) (udps[IUDP].arg[2].dot))[0]
#define CENTER(    IUDP,I)  ((double *) (udps[IUDP].arg[3].val))[I]
#define CENTER_DOT(IUDP,I)  ((double *) (udps[IUDP].arg[3].dot))[I]
#define CENTER_SIZ(IUDP  )               udps[IUDP].arg[3].size
#define AREA(      IUDP  )  ((double *) (udps[IUDP].arg[4].val))[0]
#define VOLUME(    IUDP  )  ((double *) (udps[IUDP].arg[5].val))[0]

/* data about possible arguments
      argNames: argument name
      argTypes: argument type: ATTRSTRING   string
                               ATTRINT      integer input
                              -ATTRINT      integer output
                               ATTRREAL     double  input
                              -ATTRREAL     double  output
                               ATTRREALSENS double input (for which a sensitivity can be calculated)
                               ATTRFILE     file input
      argIdefs: default value for ATTRINT
      argDdefs: default value for ATTRREAL or ATTRREALSENS */
static char  *argNames[NUMUDPARGS] = {"dx",        "dy",        "dz",        "center",    "area",    "volume", };
static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, ATTRREALSEN, ATTRREALSEN, -ATTRREAL, -ATTRREAL,};
static int    argIdefs[NUMUDPARGS] = {0,           0,           0,           0,           0,         0,        };
static double argDdefs[NUMUDPARGS] = {0.,          0.,          0.,          0.,          0.,        0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"


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
    int     wire, sense[4];
    double  node1[3], node2[3], node3[3], node4[3];
    double  data[18], trange[2];
    char    *message=NULL;
    ego     enodes[9], ecurve, eedges[8], eloop, eface;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("dx(0)         = %f\n", DX(    0));
    printf("dx_dot(0)     = %f\n", DX_DOT(0));
    printf("dy(0)         = %f\n", DY(    0));
    printf("dy_dot(0)     = %f\n", DY_DOT(0));
    printf("dz(0)         = %f\n", DZ(    0));
    printf("dz_dot(0)     = %f\n", DZ_DOT(0));
    if (CENTER_SIZ(0) == 3) {
        printf("center(0)     = %f %f %f\n", CENTER(    0,0), CENTER(    0,1), CENTER(    0,2));
        printf("center_dot(0) = %f %f %f\n", CENTER_DOT(0,0), CENTER_DOT(0,1), CENTER_DOT(0,2));
    }
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

    } else if (DX(0) <= 0 && DY(0) <= 0 && DZ(0) <= 0) {
        snprintf(message, 100, "dx=dy=dz=0");
        status  = EGADS_GEOMERR;
        goto cleanup;

    } else if (CENTER_SIZ(0) == 1) {
        // CENTER is not used

    } else if (CENTER_SIZ(0) != 3) {
        snprintf(message, 100, "center should contain 3 entries");
        status  = EGADS_GEOMERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("dx[%d]         = %f\n", numUdp, DX(    numUdp));
    printf("dx_dot[%d]     = %f\n", numUdp, DX_DOT(numUdp));
    printf("dy[%d]         = %f\n", numUdp, DY(    numUdp));
    printf("dy_dot[%d]     = %f\n", numUdp, DY_DOT(numUdp));
    printf("dz[%d]         = %f\n", numUdp, DZ(    numUdp));
    printf("dz_dot[%d]     = %f\n", numUdp, DZ_DOT(numUdp));
    if (CENTER_SIZ(numUdp) == 3) {
        printf("center[%d]     = %f %f %f\n", numUdp, CENTER(    numUdp,0), CENTER(    numUdp,1), CENTER(    numUdp,2));
        printf("center_dot[%d] = %f %f %f\n", numUdp, CENTER_DOT(numUdp,0), CENTER_DOT(numUdp,1), CENTER_DOT(numUdp,2));
    }
#endif

    /* make private data (not needed here, but included to
       show how one would do this */
    MALLOC(udps[numUdp].data, char, 30);

    strcpy((char*)(udps[numUdp].data), "this is test private data");

    /* check for 3D SolidBody (and make if requested) */
    if (DX(0) > 0 && DY(0) > 0 && DZ(0) > 0)  {

        /*
                ^ Y
                |
                4----11----8
               /:         /|
              3 :        7 |
             /  4       /  8
            3----12----7   |
            |   2-----9|---6  --> X
            |  '       |  /
            2 1        6 5
            |'         |/
            1----10----5
           /
          Z

        */

        data[0] = -DX(0) / 2;
        data[1] = -DY(0) / 2;
        data[2] = -DZ(0) / 2;
        data[3] =  DX(0);
        data[4] =  DY(0);
        data[5] =  DZ(0);

        /* move the Body if CENTER has 3 values */
        if (CENTER_SIZ(0) == 3) {
            data[0] += CENTER(0,0);
            data[1] += CENTER(0,1);
            data[2] += CENTER(0,2);
        }

        /* Make SolidBody */
        status = EG_makeSolidBody(context, BOX, data, ebody);
        CHECK_STATUS(EG_makeSolidBody);
        if (*ebody == NULL) goto cleanup;    // needed for splint

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (Body) */
        udps[numUdp].ebody = *ebody;

        goto cleanup;
    }

    /* otherwise this is either a 1D WireBody or a 2D FaceBody */

    /* geometry start */
    node1[0] = - DX(0) / 2;
    node1[1] = - DY(0) / 2;
    node1[2] = - DZ(0) / 2;

    /* WireBody dz > 0 */
    if        (DX(0) == 0 && DY(0) == 0) {
        wire = 1;
        node2[0] = 0;            node2[1] = 0;            node2[2] = DZ(0) / 2;

    /* WireBody dy > 0 */
    } else if (DX(0) == 0 && DZ(0) == 0) {
        wire = 1;
        node2[0] = 0;            node2[1] = DY(0) / 2;    node2[2] = 0;

    /* WireBody dx > 0 */
    } else if (DY(0) == 0 && DZ(0) == 0) {
        wire = 1;
        node2[0] = DX(0) / 2;    node2[1] = 0;            node2[2] = 0;

    /* Face in the x-y plane */
    } else if (DZ(0) == 0) {
        wire = 0;
        node2[0] = +DX(0) / 2;   node2[1] = -DY(0) / 2;   node2[2] = 0;
        node3[0] = +DX(0) / 2;   node3[1] = +DY(0) / 2;   node3[2] = 0;
        node4[0] = -DX(0) / 2;   node4[1] = +DY(0) / 2;   node4[2] = 0;

    /* Face in the y-z plane */
    } else if (DX(0) == 0) {
        wire = 0;
        node2[0] = 0;            node2[1] = +DY(0) / 2;   node2[2] = -DZ(0) / 2;
        node3[0] = 0;            node3[1] = +DY(0) / 2;   node3[2] = +DZ(0) / 2;
        node4[0] = 0;            node4[1] = -DY(0) / 2;   node4[2] = +DZ(0) / 2;

    /* Face in the x-z plane */
    } else {
        wire = 0;
        node2[0] = -DX(0) / 2;   node2[1] = 0;            node2[2] = +DZ(0) / 2;
        node3[0] = +DX(0) / 2;   node3[1] = 0;            node3[2] = +DZ(0) / 2;
        node4[0] = +DX(0) / 2;   node4[1] = 0;            node4[2] = -DZ(0) / 2;
    }

    /* make 1D WireBody */
    if (wire == 1) {

        /* move the Nodes if CENTER has 3 values */
        if (CENTER_SIZ(0) == 3) {
            node1[0] += CENTER(0,0);
            node1[1] += CENTER(0,1);
            node1[2] += CENTER(0,2);

            node2[0] += CENTER(0,0);
            node2[1] += CENTER(0,1);
            node2[2] += CENTER(0,2);
        }

        /* make Nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* make the Line */
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

        /* make Edge */
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
        if (*ebody == NULL) goto cleanup;    // needed for splint

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (Body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;

    /* make 2D Face Body */
    } else {

        /*
                  y,z,x
                    ^
                    :
              4----<3----3
              |     :    |
              4v    +- - 2^ - -> x,y,z
              |          |
              1----1>----2
        */

        /* move the Nodes if CENTER has 3 values */
        if (CENTER_SIZ(0) == 3) {
            node1[0] += CENTER(0,0);
            node1[1] += CENTER(0,1);
            node1[2] += CENTER(0,2);

            node2[0] += CENTER(0,0);
            node2[1] += CENTER(0,1);
            node2[2] += CENTER(0,2);

            node3[0] += CENTER(0,0);
            node3[1] += CENTER(0,1);
            node3[2] += CENTER(0,2);

            node4[0] += CENTER(0,0);
            node4[1] += CENTER(0,1);
            node4[2] += CENTER(0,2);
        }

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

        /* make Edge */
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

        /* make Edge */
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

        /* make Edge */
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

        /* make Edge */
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

        /* create the FaceBody (which will be returned) */
        status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);
        if (*ebody == NULL) goto cleanup;    // needed for splint

        /* set the output value(s) */
        status = EG_getMassProperties(*ebody, data);
        CHECK_STATUS(EG_getMassProperties);

        AREA(0)   = data[1];
        VOLUME(0) = data[0];

        /* remember this model (Body) */
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

    int    iudp, judp, i, inode, iedge, iface;
    double dx_dot, dy_dot, dz_dot;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
           (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
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

    /* WireBody */
    if        (DY(iudp) <= 0 && DZ(iudp) <= 0) {
        if (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = 0;
            } else if (inode == 2) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if (iedge == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        for (i = 0; i < npnt; i++) {
            vels[3*i  ] = dx_dot;
            vels[3*i+1] = dy_dot;
            vels[3*i+2] = dz_dot;
        }
    } else if (DZ(iudp) <= 0 && DX(iudp) <= 0) {
        if (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = 0;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (inode == 2) {
                dx_dot = 0;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if (iedge == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        for (i = 0; i < npnt; i++) {
            vels[3*i  ] = dx_dot;
            vels[3*i+1] = dy_dot;
            vels[3*i+2] = dz_dot;
        }
    } else if (DX(iudp) <= 0 && DY(iudp) <= 0) {
        if (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 2) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = +DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if (iedge == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        for (i = 0; i < npnt; i++) {
            vels[3*i  ] = dx_dot;
            vels[3*i+1] = dy_dot;
            vels[3*i+2] = dz_dot;
        }

    /* SheetBody in XY plane (since the velocity on each Node and Edge is
                              a constant, we need to just compute one value) */
    } else if (DZ(iudp) <= 0) {
        if        (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (inode == 2) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (inode == 3) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (inode == 4) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if        (iedge == 1) {
                dx_dot = 0;                 dy_dot = -DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (iedge == 2) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = 0;                 dz_dot = 0;
            } else if (iedge == 3) {
                dx_dot = 0;                 dy_dot = +DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (iedge == 4) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = 0;                 dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_FACE) {
            iface = entType;

            if (iface == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iface=%d\n", iface);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        for (i = 0; i < npnt; i++) {
            vels[3*i  ] = dx_dot;
            vels[3*i+1] = dy_dot;
            vels[3*i+2] = dz_dot;
        }

    /* SheetBody in YZ plane (since the velocity on each Node and Edge is
                              a constant, we need to just compute one value) */
    } else if (DX(iudp) <= 0) {
        if        (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = 0;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 2) {
                dx_dot = 0;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 3) {
                dx_dot = 0;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 4) {
                dx_dot = 0;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if        (iedge == 1) {
                dx_dot = 0;   dy_dot = 0;                 dz_dot = -DZ_DOT(iudp)/2;
            } else if (iedge == 2) {
                dx_dot = 0;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (iedge == 3) {
                dx_dot = 0;   dy_dot = 0;                 dz_dot = +DZ_DOT(iudp)/2;
            } else if (iedge == 4) {
                dx_dot = 0;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_FACE) {
            iface = entType;

            if (iface == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iface=%d\n", iface);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        for (i = 0; i < npnt; i++) {
            vels[3*i  ] = dx_dot;
            vels[3*i+1] = dy_dot;
            vels[3*i+2] = dz_dot;
        }

    /* SheetBody in ZX plane (since the velocity on each Node and Edge is
                              a constant, we need to just compute one value) */
    } else if (DY(iudp) <= 0) {
        if        (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 2) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 3) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 4) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = -DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if        (iedge == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = 0;
            } else if (iedge == 2) {
                dx_dot = 0;                 dy_dot = 0;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (iedge == 3) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = 0;   dz_dot = 0;
        } else if (iedge == 4) {
                dx_dot = 0;                 dy_dot = 0;   dz_dot = -DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_FACE) {
            iface = entType;

            if (iface == 1) {
                dx_dot = 0;   dy_dot = 0;   dz_dot = 0;
            } else {
                printf("udpSensitivity: bad iface=%d\n", iface);
                status = EGADS_INDEXERR;
                goto cleanup;
            }

        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        for (i = 0; i < npnt; i++) {
            vels[3*i  ] = dx_dot;
            vels[3*i+1] = dy_dot;
            vels[3*i+2] = dz_dot;
        }

    /* SolidBody (since the velocity on each Node, Edge, and Face is a
                  constant, we need just to compute one vaue) */
    } else {
        if        (entType == OCSM_NODE) {
            inode = entIndex;

            if        (inode == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 2) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 3) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 4) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 5) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 6) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (inode == 7) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (inode == 8) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad inode=%d\n", inode);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_EDGE) {
            iedge = entIndex;

            if        (iedge == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot =  0;
            } else if (iedge == 2) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot =  0;                dz_dot = +DZ_DOT(iudp)/2;
            } else if (iedge == 3) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot =  0;
            } else if (iedge == 4) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot =  0;                dz_dot = -DZ_DOT(iudp)/2;
            } else if (iedge == 5) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = -DY_DOT(iudp)/2;   dz_dot =  0;
            } else if (iedge == 6) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot =  0;                dz_dot = +DZ_DOT(iudp)/2;
            } else if (iedge == 7) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = +DY_DOT(iudp)/2;   dz_dot =  0;
            } else if (iedge == 8) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot =  0;                dz_dot = -DZ_DOT(iudp)/2;
            } else if (iedge == 9) {
                dx_dot =  0;                dy_dot = -DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (iedge == 10) {
                dx_dot =  0;                dy_dot = -DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else if (iedge == 11) {
                dx_dot =  0;                dy_dot = +DY_DOT(iudp)/2;   dz_dot = -DZ_DOT(iudp)/2;
            } else if (iedge == 12) {
                dx_dot =  0;                dy_dot = +DY_DOT(iudp)/2;   dz_dot = +DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad iedge=%d\n", iedge);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else if (entType == OCSM_FACE) {
            iface = entIndex;

            if        (iface == 1) {
                dx_dot = -DX_DOT(iudp)/2;   dy_dot = 0;                 dz_dot = 0;
            } else if (iface == 2) {
                dx_dot = +DX_DOT(iudp)/2;   dy_dot = 0;                 dz_dot = 0;
            } else if (iface == 3) {
                dx_dot = 0;                 dy_dot = -DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (iface == 4) {
                dx_dot = 0;                 dy_dot = +DY_DOT(iudp)/2;   dz_dot = 0;
            } else if (iface == 5) {
                dx_dot = 0;                 dy_dot = 0;                 dz_dot = -DZ_DOT(iudp)/2;
            } else if (iface == 6) {
                dx_dot = 0;                 dy_dot = 0;                 dz_dot = +DZ_DOT(iudp)/2;
            } else {
                printf("udpSensitivity: bad iface=%d\n", iface);
                status = EGADS_INDEXERR;
                goto cleanup;
            }
        } else {
            printf("udpSensitivity: bad entType=%d\n", entType);
            status = EGADS_ATTRERR;
            goto cleanup;
        }

        /* account for the movement of the CENTER if it has 3 values */
        if (CENTER_SIZ(iudp) == 3) {
            for (i = 0; i < npnt; i++) {
                vels[3*i  ] = CENTER_DOT(iudp,0) + dx_dot;
                vels[3*i+1] = CENTER_DOT(iudp,1) + dy_dot;
                vels[3*i+2] = CENTER_DOT(iudp,2) + dz_dot;
            }
        }
    }

    status = EGADS_SUCCESS;

cleanup:
    return status;
}



/*
 ************************************************************************
 *                                                                      *
 *   freePrivateData - free private data (just an example)              *
 *                                                                      *
 ************************************************************************
 */

static int
freePrivateData(void  *data)            /* (in)  pointer to private data */
{
    int    status = EGADS_SUCCESS;

    printf("freePrivateData(%s)\n", (char*)(data));

    /* note: this function would not be necessary if we are only calling
             EG_free (and then FREEUDPDATA would not be defined above).
             It is simply included here to show how one would write
             such a function if the allocation was more complicated
             than a simple EG_alloc() */

    EG_free(data);

//cleanup:
    return status;
}
