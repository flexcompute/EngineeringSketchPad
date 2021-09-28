/*
 *************************************************************************
 *                                                                       *
 * udpParaboloid --- create a paraboloid or parabola                     *
 *                                                                       *
 * written  by William Saueressig                                        *
 * modified by John Dannenhoffer                                         *
 *                                                                       *
 *************************************************************************
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

/* the number of "input" Bodys

   this only needs to be specified if this is a UDF (user-defined
   function) that consumes Bodys from OpenCSM's stack.  (the default
   value is 0).  the Bodys are passed into udpExecute via its
   emodel orgument */

#define NUMUDPINPUTBODYS 0

/* the number of arguments (specified below) */
#define NUMUDPARGS 3

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define XLENGTH(IUDP) ((double *) (udps[IUDP].arg[0].val))[0]
#define YRADIUS(IUDP) ((double *) (udps[IUDP].arg[1].val))[0]
#define ZRADIUS(IUDP) ((double *) (udps[IUDP].arg[2].val))[0]

/* data about possible arguments
      argNames: argument name
      argTypes: argument type: ATTRSTRING   string
                               ATTRINT      integer input
                              -ATTRINT      integer output
                               ATTRREAL     double  input
                              -ATTRREAL     double  output
                               ATTRREALSENS double input (for which a sensitivity can be calculated)
      argIdefs: default value for ATTRINT
      argDdefs: default value for ATTRREAL or ATTRREALSENS */
static char  *argNames[NUMUDPARGS] = {"xlength", "yradius", "zradius",};
static int    argTypes[NUMUDPARGS] = {ATTRREAL,  ATTRREAL,  ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,         0,         0,        };
static double argDdefs[NUMUDPARGS] = {0.,        0.,        0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototype for function defined below */
static int convertToBSpline(ego inbody, double matrix[], ego *outbody);


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

    /* declare status to EGADS_SUCCESS (0) */
    int     status = EGADS_SUCCESS;

    /* declare variables */
    int     psens[3];
    int     oclass, mtype, nchild, *senses;

    double  node0[3], node1[3], node2[3], ymirror[12], zmirror[12];
    double  lineGeom[6], paraGeom[10], planGeom[9], axisGeom[6], trange[2], tempf[3];
    double  xyz[18], mat[12];
    char    *message=NULL;

    ego     enodes[3], ecurve[3], eedges[3], eloop, esurf, eface;
    ego     equarter[2], ehalf[2], efull, *ebodys, eref, echild[3], etran;

    ROUTINE(udpExecute);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("xlength(0)    = %f\n", XLENGTH(0));
    printf("yradius(0)    = %f\n", YRADIUS(0));
    printf("zradius(0)    = %f\n", ZRADIUS(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "udpParaboloid: xlength should be a scalar.\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (XLENGTH(0) <= 0) {
        snprintf(message, 100, "udpParaboloid: xlength has to be greater than 0\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "udpParaboloid: yradius should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (YRADIUS(0) < 0) {
        snprintf(message, 100, "udpParaboloid: yradius cannot be negative\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "udpParaboloid: xlength should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (ZRADIUS(0) < 0) {
        snprintf(message, 100, "udpParaboloid: zradius cannot be negative\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (YRADIUS(0) <= 0 && ZRADIUS(0) <= 0) {
        snprintf(message, 100, "udpParaboloid: yradius and zradius cannot both be zero");
        status  = EGADS_GEOMERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("xlength[%d]    = %f\n", numUdp, XLENGTH(numUdp));
    printf("yradius[%d]    = %f\n", numUdp, YRADIUS(numUdp));
    printf("zradius[%d]    = %f\n", numUdp, ZRADIUS(numUdp));
#endif

    /*
     ******************************************************************************
     *    Make 3D Solid Body if xlength, yradius, and zradius are all positive    *
     ******************************************************************************
     */

    if (XLENGTH(0) > 0 && YRADIUS(0) > 0 && ZRADIUS(0) > 0)  {

        /* set values of the Nodes */
        // Node at origin
        node0[0] = 0;
        node0[1] = 0;
        node0[2] = 0;

        // Node along x-axis
        node1[0] = XLENGTH(0);
        node1[1] = 0;
        node1[2] = 0;

        // Node on parabola
        node2[0] = XLENGTH(0);
        node2[1] = 1;
        node2[2] = 0;

        /* make the Nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node0, 0, NULL, NULL, &enodes[0]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &enodes[1]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &enodes[2]);
        CHECK_STATUS(EG_makeTopology);

        /* geometry for Line along global x-axis */
        lineGeom[0] = node0[0];   lineGeom[1] = node0[1];   lineGeom[2] = node0[2];
        lineGeom[3] = node1[0];   lineGeom[4] = node1[1];   lineGeom[5] = node1[2];

        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, lineGeom, &ecurve[0]);
        CHECK_STATUS(EG_makeGeometry);

        /* geometry for Line parallel to global y-axis */
        lineGeom[0] = node1[0];
        lineGeom[1] = node1[1];
        lineGeom[2] = node1[2];

        lineGeom[3] = node2[0] - node1[0];
        lineGeom[4] = node2[1] - node1[1];
        lineGeom[5] = node2[2] - node1[2];

        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, lineGeom, &ecurve[1]);
        CHECK_STATUS(EG_makeGeometry);

        /* geometry for the PARABOLA */
        paraGeom[0] = 0;   paraGeom[1] = 0;   paraGeom[2] = 0;
        paraGeom[3] = 1;   paraGeom[4] = 0;   paraGeom[5] = 0;
        paraGeom[6] = 0;   paraGeom[7] = 1;   paraGeom[8] = 0;

        /* focus = 1 / (4 * a) */
        paraGeom[9] = 0.25 / XLENGTH(0);

        /* make the Curve */
        status = EG_makeGeometry(context, CURVE, PARABOLA, NULL, NULL, paraGeom, &ecurve[2]);
        CHECK_STATUS(EG_makeGeometry);

        /* get the trange for CURVE 0 */
        status = EG_invEvaluate(ecurve[0], node0, &(trange[0]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve[0], node1, &(trange[1]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        /* set the children */
        echild[0] = enodes[0];
        echild[1] = enodes[1];

#ifdef DEBUG
        printf("Make Edge 0, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* make Edge along global x-axis */
        status = EG_makeTopology(context, ecurve[0], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[0]);
        CHECK_STATUS(EG_makeTopology);

        /* get the trange for Curve 1 */
        status = EG_invEvaluate(ecurve[1], node1, &trange[0], tempf);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve[1], node2, &trange[1], tempf);
        CHECK_STATUS(EG_invEvaluate);

        /* set the children */
        echild[0] = enodes[1];
        echild[1] = enodes[2];

#ifdef DEBUG
        printf("Make Edge 1, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* make Edge parallel to global y-axis */
        status = EG_makeTopology(context, ecurve[1], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[1]);
        CHECK_STATUS(EG_makeTopology);

        /* get the trange for Curve 3 */
        trange[0] = 0;
        trange[1] = 1;

        /* set the children */
        echild[0] = enodes[0];
        echild[1] = enodes[2];

#ifdef DEBUG
        printf("Make Edge 2, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* Make the edge for parabola section */
        status = EG_makeTopology(context, ecurve[2], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[2]);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Done making EdgeS.\n");
#endif

#ifdef DEBUG
        printf("Make Loops.\n");
#endif

        /* Make a Loop */
        echild[0] = eedges[0];
        echild[1] = eedges[1];
        echild[2] = eedges[2];

        psens[0] = SFORWARD;
        psens[1] = SFORWARD;
        psens[2] = SREVERSE;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 3, echild, psens, &eloop);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Make Surface.\n");
#endif

        /* Make the Surface */
        planGeom[0] = 0;   planGeom[1] = 0;   planGeom[2] = 0;
        planGeom[3] = 1;   planGeom[4] = 0;   planGeom[5] = 0;
        planGeom[6] = 0;   planGeom[7] = 1;   planGeom[8] = 0;

        status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, planGeom, &esurf);
        CHECK_STATUS(EG_makeGeometry);

#ifdef DEBUG
        printf("Make Face.\n");
#endif

        /* Create the Face */
        status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL, 1, &eloop, psens, &eface);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Rotate to create the body.\n");
#endif

        /* Rotate to create a Body */
        axisGeom[0] = 0;   axisGeom[1] = 0;   axisGeom[2] = 0;
        axisGeom[3] = 1;   axisGeom[4] = 0;   axisGeom[5] = 0;

        status = EG_rotate(eface, 90, axisGeom, &equarter[0]);
        CHECK_STATUS(EG_rotate);

#ifdef DEBUG
        printf("Converting the quarter paraboloid into the full paraboloid.\n");
#endif

        /* mirror the quarter parabola over the y axis */
        ymirror[ 0] =  1;   ymirror[ 1] =  0;   ymirror[ 2] =  0;   ymirror[ 3] =  0;
        ymirror[ 4] =  0;   ymirror[ 5] = -1;   ymirror[ 6] =  0;   ymirror[ 7] =  0;
        ymirror[ 8] =  0;   ymirror[ 9] =  0;   ymirror[10] =  1;   ymirror[11] =  0;

        status = EG_makeTransform(context, ymirror, &etran);
        CHECK_STATUS(EG_makeTransform);

        /* copy the Object */
        status = EG_copyObject(equarter[0], etran, &equarter[1]);
        CHECK_STATUS(EG_copyObject);

        /* fuse the two quarters together */
        status = EG_generalBoolean(equarter[0], equarter[1], FUSION, 0.0, &ehalf[0]);
        CHECK_STATUS(EG_generalBoolean);

        status = EG_deleteObject(equarter[0]);
        CHECK_STATUS(EG_deleteObject);

        status = EG_deleteObject(equarter[1]);
        CHECK_STATUS(EG_deleteObject);

        /* mirror the half parabola across the z axis */
        zmirror[ 0] =  1;   zmirror[ 1] =  0;   zmirror[ 2] =  0;   zmirror[ 3] =  0;
        zmirror[ 4] =  0;   zmirror[ 5] =  1;   zmirror[ 6] =  0;   zmirror[ 7] =  0;
        zmirror[ 8] =  0;   zmirror[ 9] =  0;   zmirror[10] = -1;   zmirror[11] =  0;

        status = EG_makeTransform(context, zmirror, &etran);
        CHECK_STATUS(EG_makeTransform);

        /* copy the Object */
        status = EG_copyObject(ehalf[0], etran, &ehalf[1]);
        CHECK_STATUS(EG_copyObject);

        /* fuse the two halfs together */
        status = EG_generalBoolean(ehalf[0], ehalf[1], FUSION, 0.0, &efull);
        CHECK_STATUS(EG_generalBoolean);

        status = EG_deleteObject(ehalf[0]);
        CHECK_STATUS(EG_deleteObject);

        status = EG_deleteObject(ehalf[1]);
        CHECK_STATUS(EG_deleteObject);

        /* scale Y and Z radii */
        status = EG_getTopology(efull, &eref, &oclass, &mtype,
                                xyz, &nchild, &ebodys, &senses);
        CHECK_STATUS(EG_getTopology);

        /* set up transformation matrix */
        mat[ 0] = 1;           mat[ 1] = 0;            mat[ 2] = 0;            mat[ 3] = 0;
        mat[ 4] = 0;           mat[ 5] = YRADIUS(0);   mat[ 6] = 0;            mat[ 7] = 0;
        mat[ 8] = 0;           mat[ 9] = 0;            mat[10] = ZRADIUS(0);   mat[11] = 0;

        status = convertToBSpline(ebodys[0], mat, ebody);
        CHECK_STATUS(EG_getTopology);

     /*
      **********************************************
      *    Create a 2D parabola in the X-Z plane   *
      **********************************************
     */
    } else if (YRADIUS(0) == 0 ) {

        /* set values of the Nodes */

        // Node on top end of parabolas
        node0[0] =  XLENGTH(0);
        node0[1] =  0;
        node0[2] =  ZRADIUS(0);

        // Node on bottom end of parabolas
        node1[0] =  XLENGTH(0);
        node1[1] =  0;
        node1[2] = -ZRADIUS(0);

        // Node at the origin
        node2[0] =  XLENGTH(0);
        node2[1] =  0;
        node2[2] =  0;

        /* make the Nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node0, 0, NULL, NULL, &enodes[0]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &enodes[1]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &enodes[2]);
        CHECK_STATUS(EG_makeTopology);

        /* geometry for Line along global y-axis */
        lineGeom[0] = node1[0];
        lineGeom[1] = node1[1];
        lineGeom[2] = node1[2];

        lineGeom[3] = node0[0] - node1[0];
        lineGeom[4] = node0[1] - node1[1];
        lineGeom[5] = node0[2] - node1[2];

        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, lineGeom, &ecurve[0]);
        CHECK_STATUS(EG_makeGeometry);

        /* geometry for the PARABOLA */
        paraGeom[0] = 0;   paraGeom[1] = 0;   paraGeom[2] = 0;
        paraGeom[3] = 1;   paraGeom[4] = 0;   paraGeom[5] = 0;
        paraGeom[6] = 0;   paraGeom[7] = 0;   paraGeom[8] = 1;

        // focus = r^2 / (4 * a)
        paraGeom[9] = ZRADIUS(0)*ZRADIUS(0) / (4*XLENGTH(0));

        /* make the Curve */
        status = EG_makeGeometry(context, CURVE, PARABOLA, NULL, NULL, paraGeom, &ecurve[1]);
        CHECK_STATUS(EG_makeGeometry);

        /* get the trange for Curve 0 */
        status = EG_invEvaluate(ecurve[0], node1, &(trange[0]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve[0], node0, &(trange[1]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        /* set the children */
        echild[0] = enodes[1];
        echild[1] = enodes[0];

#ifdef DEBUG
        printf("Make Edge 0, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* make Edge along global Z-axis */
        status = EG_makeTopology(context, ecurve[0], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[0]);
        CHECK_STATUS(EG_makeTopology);

        /* get the trange for Curve 1 */
        status = EG_invEvaluate(ecurve[1], node1, &(trange[0]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve[1], node0, &(trange[1]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        /* set the children */
        echild[0] = enodes[1];
        echild[1] = enodes[0];

#ifdef DEBUG
        printf("Make Edge 1, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* Make the edge for parabola section */
        status = EG_makeTopology(context, ecurve[1], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[1]);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Done making EdgeS.\n");
#endif

#ifdef DEBUG
        printf("Make Loops.\n");
#endif

        /* Make a Loop */
        echild[0] = eedges[0];
        echild[1] = eedges[1];

        psens[0] = SFORWARD;
        psens[1] = SREVERSE;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 2, echild, psens, &eloop);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Make Surface.\n");
#endif

        /* Make the Surface */
        planGeom[0] = 0;   planGeom[1] = 0;   planGeom[2] = 0;
        planGeom[3] = 1;   planGeom[4] = 0;   planGeom[5] = 0;
        planGeom[6] = 0;   planGeom[7] = 0;   planGeom[8] = 1;

        status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, planGeom, &esurf);
        CHECK_STATUS(EG_makeGeometry);

#ifdef DEBUG
        printf("Make Face.\n");
#endif

        /* Create the face */
        status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL, 1, &eloop, psens, &eface);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

     /*
      **********************************************
      *    Create a 2D parabola in the X-Y plane   *
      **********************************************
     */
    } else if (ZRADIUS(0) == 0 ) {

        /* set values of the Nodes */

        // Node on top end of parabolas
        node0[0] =  XLENGTH(0);
        node0[1] =  YRADIUS(0);
        node0[2] =  0;

        // Node on bottom end of parabolas
        node1[0] =  XLENGTH(0);
        node1[1] = -YRADIUS(0);
        node1[2] =  0;

        // Node at the origin
        node2[0] =  XLENGTH(0);
        node2[1] =  0;
        node2[2] =  0;

        /* make the Nodes */
        status = EG_makeTopology(context, NULL, NODE, 0, node0, 0, NULL, NULL, &enodes[0]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &enodes[1]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &enodes[2]);
        CHECK_STATUS(EG_makeTopology);


        /* geometry for Line along global y-axis */
        lineGeom[0] = node1[0];
        lineGeom[1] = node1[1];
        lineGeom[2] = node1[2];

        lineGeom[3] = node0[0] - node1[0];
        lineGeom[4] = node0[1] - node1[1];
        lineGeom[5] = node0[2] - node1[2];

        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, lineGeom, &ecurve[0]);
        CHECK_STATUS(EG_makeGeometry);

        /* geometry for the PARABOLA */
        paraGeom[0] = 0;   paraGeom[1] = 0;   paraGeom[2] = 0;
        paraGeom[3] = 1;   paraGeom[4] = 0;   paraGeom[5] = 0;
        paraGeom[6] = 0;   paraGeom[7] = 1;   paraGeom[8] = 0;

        // focus = r^2 / (4 * a)
        paraGeom[9] = YRADIUS(0)*YRADIUS(0) / (4*XLENGTH(0));

        /* make the Curve */
        status = EG_makeGeometry(context, CURVE, PARABOLA, NULL, NULL, paraGeom, &ecurve[1]);
        CHECK_STATUS(EG_makeGeometry);

        /* get the trange for Curve 0 */
        status = EG_invEvaluate(ecurve[0], node1, &(trange[0]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve[0], node0, &(trange[1]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        /* set the children */
        echild[0] = enodes[1];
        echild[1] = enodes[0];

#ifdef DEBUG
        printf("Make Edge 0, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* make Edge along global Y-axis */
        status = EG_makeTopology(context, ecurve[0], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[0]);
        CHECK_STATUS(EG_makeTopology);

        /* get the trange for Curve 1 */
        status = EG_invEvaluate(ecurve[1], node1, &(trange[0]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve[1], node0, &(trange[1]), tempf);
        CHECK_STATUS(EG_invEvaluate);

        /* set the children */
        echild[0] = enodes[1];
        echild[1] = enodes[0];

#ifdef DEBUG
        printf("Make Edge 1, trange[0]=%f, trange[1]=%f\n", trange[0], trange[1]);
#endif

        /* Make the edge for parabola section */
        status = EG_makeTopology(context, ecurve[1], EDGE, TWONODE, trange, 2, echild, NULL, &eedges[1]);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Done making EdgeS.\n");
#endif

#ifdef DEBUG
        printf("Make Loops.\n");
#endif

        /* Make a Loop */
        echild[0] = eedges[0];
        echild[1] = eedges[1];

        psens[0] = SFORWARD;
        psens[1] = SREVERSE;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 2, echild, psens, &eloop);
        CHECK_STATUS(EG_makeTopology);

#ifdef DEBUG
        printf("Make Surface.\n");
#endif

        /* Make the Surface */
        planGeom[0] = 0;   planGeom[1] = 0;   planGeom[2] = 0;
        planGeom[3] = 1;   planGeom[4] = 0;   planGeom[5] = 0;
        planGeom[6] = 0;   planGeom[7] = 1;   planGeom[8] = 0;

        status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, planGeom, &esurf);
        CHECK_STATUS(EG_makeGeometry);

#ifdef DEBUG
        printf("Make Face.\n");
#endif

        /* Create the face */
        status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL, 1, &eloop, psens, &eface);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

    }

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;

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
   /*@unused@*/int    npnt,             /* (in)  number of points */
   /*@unused@*/int    entType,          /* (in)  OCSM entity type */
   /*@unused@*/int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
   /*@unused@*/double vels[])           /* (out) velocities */
{
    int iudp, judp;

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


/*
 ************************************************************************
 *                                                                      *
 *   convertToBSpline - conert Body to BSpline with transformation      *
 *                                                                      *
 ************************************************************************
 */

static int
convertToBSpline(ego    inbody,         /* (in)  input Body */
                 double mat[],          /* (in)  transformation matrix */
                 ego    *ebody)         /* (out) pointer to output Body */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, ichild, nchild, *senses;
    int     inode, nnode, iedge, nedge, iface, nface, iloop, nloop, ishell, nshell, *ivals, i, sense;
    int     oclass1, mtype1, nchild1, *senses1, jedge, periodic;
    double  xyz[18], xyz2[18], *rvals, trange[2], xyzClose[3];
    double  uv[18], range[4], tt[2], swap;
    ego     context, eref, *echilds, *enodes_orig=NULL, *eedges_orig=NULL, *efaces_orig=NULL, *eloops_orig=NULL, *eshells_orig=NULL;
    ego     *enodes=NULL, *ecurves=NULL, *eedges=NULL, *esurfaces=NULL, *efaces=NULL, *eshells=NULL;
    ego     *eloops2=NULL, *eloops=NULL, *efoo=NULL;
    ego     ecurve, esurface, etemp[2];
            ego *echilds1;

    ROUTINE(convertToBspline);

    /* --------------------------------------------------------------- */

    status = EG_getContext(inbody, &context);
    CHECK_STATUS(EG_getContext);

#ifdef DEBUG
    printf("convertToBSpline(mat=%10.5f %10.5f %10.5f %10.5f\n", mat[ 0], mat[ 1], mat[ 2], mat[ 3]);
    printf("                     %10.5f %10.5f %10.5f %10.5f\n", mat[ 4], mat[ 5], mat[ 6], mat[ 7]);
    printf("                     %10.5f %10.5f %10.5f %10.5f\n", mat[ 8], mat[ 9], mat[10], mat[11]);

    (void) ocsmPrintEgo(inbody);
#endif

    /* get the Nodes, Edges, and Faces associated with the input Body */
    status = EG_getBodyTopos(inbody, NULL, NODE, &nnode, &enodes_orig);
    CHECK_STATUS(EG_getBodyTopos);

#ifdef DEBUG
    printf("error isnt the node");
#endif

    status = EG_getBodyTopos(inbody, NULL, EDGE, &nedge, &eedges_orig);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(inbody, NULL, FACE, &nface, &efaces_orig);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(inbody, NULL, SHELL, &nshell, &eshells_orig);
    CHECK_STATUS(EG_getBodyTopos);

#ifdef DEBUG
    printf("nnode=%d, nedge=%d, nface=%d, nshell=%d\n", nnode, nedge, nface, nshell);
#endif

    /* get storage for the new Nodes, Curves, Edges, Surfaces, Faces, and Shells */
    MALLOC(enodes,    ego, nnode );
    MALLOC(ecurves,   ego, nedge );
    MALLOC(eedges,    ego, nedge );
    MALLOC(esurfaces, ego, nface );
    MALLOC(efaces,    ego, nface );
    MALLOC(eshells,   ego, nshell);

    /* make new Nodes by transforming the original Nodes */
    for (inode = 0; inode < nnode; inode++) {
        SPLINT_CHECK_FOR_NULL(enodes_orig );

        status = EG_getTopology(enodes_orig[inode], &eref, &oclass, &mtype,
                                xyz2, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        xyz[0] = mat[ 0] * xyz2[0] + mat[ 1] * xyz2[1] + mat[ 2] * xyz2[2] + mat[ 3];
        xyz[1] = mat[ 4] * xyz2[0] + mat[ 5] * xyz2[1] + mat[ 6] * xyz2[2] + mat[ 7];
        xyz[2] = mat[ 8] * xyz2[0] + mat[ 9] * xyz2[1] + mat[10] * xyz2[2] + mat[11];

        status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, &enodes[inode]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeDup(enodes_orig[inode], enodes[inode]);
        CHECK_STATUS(EG_attributeDup);
    }

    /* special processing for NodeBody */
    if (nnode == 1) {
        MALLOC(eloops, ego, 1);

        trange[0] = 0;
        trange[1] = 1;
        sense     = SFORWARD;

        status = EG_makeTopology(context, NULL, EDGE, DEGENERATE, trange, 1, enodes, &sense, &eedges[0]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 1, eedges, &sense, eloops);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, eloops, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        FREE(eloops);

        goto cleanup;
    }

    /* make new Edges by transforming the original Edges */
    for (iedge = 0; iedge < nedge; iedge++) {
        SPLINT_CHECK_FOR_NULL(eedges_orig );

        status = EG_getTopology(eedges_orig[iedge], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if (mtype == DEGENERATE) {
            eedges[iedge] = NULL;
            continue;
        }

        status = EG_convertToBSpline(eedges_orig[iedge], &ecurve);
        CHECK_STATUS(EG_convertToBSpline);

        status = EG_getGeometry(ecurve, &oclass, &mtype, &eref, &ivals, &rvals);
        CHECK_STATUS(EG_getGeometry);

        for (i = ivals[3]; i < ivals[3]+3*ivals[2]; i+=3) {
            xyz2[0] = rvals[i  ];
            xyz2[1] = rvals[i+1];
            xyz2[2] = rvals[i+2];

            rvals[i  ] = mat[ 0] * xyz2[0] + mat[ 1] * xyz2[1] + mat[ 2] * xyz2[2] + mat[ 3];
            rvals[i+1] = mat[ 4] * xyz2[0] + mat[ 5] * xyz2[1] + mat[ 6] * xyz2[2] + mat[ 7];
            rvals[i+2] = mat[ 8] * xyz2[0] + mat[ 9] * xyz2[1] + mat[10] * xyz2[2] + mat[11];
        }

        status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, ivals, rvals, &ecurves[iedge]);
        CHECK_STATUS(EG_makeGeometry);

        EG_free(ivals);
        EG_free(rvals);

        status = EG_getTopology(eedges_orig[iedge], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        etemp[0] = enodes[EG_indexBodyTopo(inbody, echilds[       0])-1];
        etemp[1] = enodes[EG_indexBodyTopo(inbody, echilds[nchild-1])-1];

        status = EG_getTopology(etemp[0], &eref, &oclass1, &mtype1,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_invEvaluate(ecurves[iedge], xyz, &trange[0], xyzClose);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_getTopology(etemp[1], &eref, &oclass1, &mtype1,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_invEvaluate(ecurves[iedge], xyz, &trange[1], xyzClose);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[0] > trange[1]) {
            trange[0] -= TWOPI;
        }

        status = EG_makeTopology(context, ecurves[iedge], oclass, mtype, trange,
                                 2, etemp, NULL, &eedges[iedge]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeDup(eedges_orig[iedge], eedges[iedge]);
        CHECK_STATUS(EG_attributeDup);
    }

    /* special processing for WireBody */
    if (nface == 0) {
        status = EG_getBodyTopos(inbody, NULL, LOOP, &nloop, &eloops_orig);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(eloops_orig);

        MALLOC(eloops, ego, nloop);

        status = EG_getTopology(eloops_orig[0], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        MALLOC(efoo, ego, nchild);

        for (iedge = 0; iedge < nchild; iedge++) {
            efoo[iedge] = eedges[EG_indexBodyTopo(inbody, echilds[iedge])-1];
        }

        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, efoo, senses, &eloops[0]);
        CHECK_STATUS(EG_makeTopology);

        FREE(efoo  );

        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, eloops, senses, ebody);
        CHECK_STATUS(EG_makeTopology);

        FREE(eloops);

        goto cleanup;
    }

    /* make new Faces by transforming the original Faces */
    for (iface = 0; iface < nface; iface++) {
        SPLINT_CHECK_FOR_NULL(efaces_orig );

        status = EG_getTopology(efaces_orig[iface], &eref, &oclass, &mtype,
                                xyz, &nloop, &eloops2, &senses);
        CHECK_STATUS(EG_getTopology);

        SPLINT_CHECK_FOR_NULL(eloops2);

        status = EG_convertToBSpline(efaces_orig[iface], &esurface);
        CHECK_STATUS(EG_convertToBSpline);

        status = EG_getGeometry(esurface, &oclass, &mtype, &eref, &ivals, &rvals);
        CHECK_STATUS(EG_getGeometry);

        for (i = ivals[3]+ivals[6]; i < ivals[3]+ivals[6]+3*ivals[2]*ivals[5]; i+=3) {
            xyz2[0] = rvals[i  ];
            xyz2[1] = rvals[i+1];
            xyz2[2] = rvals[i+2];

            rvals[i  ] = mat[ 0] * xyz2[0] + mat[ 1] * xyz2[1] + mat[ 2] * xyz2[2] + mat[ 3];
            rvals[i+1] = mat[ 4] * xyz2[0] + mat[ 5] * xyz2[1] + mat[ 6] * xyz2[2] + mat[ 7];
            rvals[i+2] = mat[ 8] * xyz2[0] + mat[ 9] * xyz2[1] + mat[10] * xyz2[2] + mat[11];
        }

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, ivals, rvals, &esurfaces[iface]);
        CHECK_STATUS(EG_makeGeometry);

        EG_free(ivals);
        EG_free(rvals);

        /* process the Loops */
        MALLOC(eloops, ego, nloop);

        for (iloop = 0; iloop < nloop; iloop++) {
            status = EG_getTopology(eloops2[iloop], &eref, &oclass, &mtype,
                                    xyz, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            MALLOC(efoo, ego, 2*nchild);

            /* process each Edge in each Loop (but store NULLs for DEGENERATE Edges */
            for (iedge = 0; iedge < nchild; iedge++) {
                status = EG_getTopology(echilds[iedge], &eref, &oclass1, &mtype1,
                                        xyz, &nchild1, &echilds1, &senses1);
                CHECK_STATUS(EG_getTopology);

                if (mtype1 != DEGENERATE) {
                    efoo[iedge] = eedges[EG_indexBodyTopo(inbody, echilds[iedge])-1];

                    status = EG_otherCurve(esurfaces[iface], efoo[iedge], 0, &efoo[iedge+nchild]);
                    CHECK_STATUS(EG_otherCurve);
                } else {
                    efoo[iedge       ] = NULL;
                    efoo[iedge+nchild] = NULL;
                }
            }

            /* if any of the Edges were NULL (DEGENERATE), proces them now */
            for (iedge = 0; iedge < nchild; iedge++) {
                if (efoo[iedge] != NULL) continue;

                /* Edge before iedge */
                if (iedge > 0) {
                    jedge = iedge - 1;
                } else {
                    jedge = nchild - 1;
                }

                status = EG_getTopology(efoo[jedge], &eref, &oclass1, &mtype1,
                                        xyz, &nchild1, &echilds1, &senses1);
                CHECK_STATUS(EG_getTopology);

                status = EG_getRange(efoo[jedge], tt, &periodic);
                CHECK_STATUS(EG_getRange);

                if (senses[jedge] == SFORWARD) {
                    etemp[0] = echilds1[1];
                    status = EG_evaluate(efoo[jedge+nchild], &tt[1], &uv[0]);
                    CHECK_STATUS(EG_evaluate);
                } else {
                    etemp[0] = echilds1[0];
                    status = EG_evaluate(efoo[jedge+nchild], &tt[0], &uv[0]);
                    CHECK_STATUS(EG_evaluate);
                }

                /* Edge after iedge */
                if (iedge < nchild-1) {
                    jedge = iedge + 1;
                } else {
                    jedge = 0;
                }

                status = EG_getRange(efoo[jedge], tt, &periodic);
                CHECK_STATUS(EG_getRange);

                if (senses[jedge] == SFORWARD) {
                    status = EG_evaluate(efoo[jedge+nchild], &tt[0], &uv[2]);
                    CHECK_STATUS(EG_evaluate);
                } else {
                    status = EG_evaluate(efoo[jedge+nchild], &tt[1], &uv[2]);
                    CHECK_STATUS(EG_evaluate);
                }

                /* make the DEGENERATE Edge and associated PCurve */
                if (senses[iedge] != SFORWARD) {
                    swap  = uv[2];
                    uv[2] = uv[0];
                    uv[0] = swap;

                    swap  = uv[3];
                    uv[3] = uv[1];
                    uv[1] = swap;
                }

                uv[2] -= uv[0];
                uv[3] -= uv[1];

                range[0] = 0;
                range[1] = sqrt(uv[2]*uv[2] + uv[3]*uv[3]);
                sense    = SFORWARD;

                status = EG_makeTopology(context, NULL, EDGE, DEGENERATE, range, 1, etemp, &sense, &efoo[iedge]);
                CHECK_STATUS(EG_makeTopology);

                status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, uv, &efoo[iedge+nchild]);
                CHECK_STATUS(EG_makeGeometry);
            }

            /* make the Loop */
            status = EG_makeTopology(context, esurfaces[iface], oclass, mtype,
                                     NULL, nchild, efoo, senses, &eloops[iloop]);
            CHECK_STATUS(EG_makeTopology);

            FREE(efoo);
        }

        /* make the Face */
        status = EG_getTopology(efaces_orig[iface], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_makeTopology(context, esurfaces[iface], oclass, mtype,
                                 xyz, nloop, eloops, senses, &efaces[iface]);
        CHECK_STATUS(EG_makeTopology);

        FREE(eloops);

        status = EG_attributeDup(efaces_orig[iface], efaces[iface]);
        CHECK_STATUS(EG_attributeDup);
    }

    /* make the Shells */
    for (ishell = 0; ishell < nshell; ishell++) {
        SPLINT_CHECK_FOR_NULL(eshells_orig);

        status = EG_getTopology(eshells_orig[ishell], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        MALLOC(efoo, ego, nchild);

        for (ichild = 0; ichild < nchild; ichild++) {
            iface = EG_indexBodyTopo(inbody, echilds[ichild]) - 1;
            efoo[ichild] = efaces[iface];
        }

        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, efoo, senses, &eshells[ishell]);
        CHECK_STATUS(EG_makeTopology);
    }

    /* finally make the Body */
    status = EG_getTopology(inbody, &eref, &oclass, &mtype,
                            xyz, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (mtype == SHEETBODY || mtype == SOLIDBODY) {
        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, eshells, senses, ebody);
        CHECK_STATUS(EG_makeTopology);
    } else if (mtype == FACEBODY && nchild == 1) {
        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, efaces, senses, ebody);
        CHECK_STATUS(EG_makeTopology);
    } else {
        printf(" convertToBSpline: oclass=%d, mtype=%d, nchild=%d\n", oclass, mtype, nchild);
        status = EGADS_TOPOERR;
        CHECK_STATUS(EG_makeTopology);
    }

    status = EG_attributeDup(inbody, *ebody);
    CHECK_STATUS(EG_attributeDup);

cleanup:
    if (enodes_orig  != NULL) EG_free(enodes_orig );
    if (eedges_orig  != NULL) EG_free(eedges_orig );
    if (eloops_orig  != NULL) EG_free(eloops_orig );
    if (efaces_orig  != NULL) EG_free(efaces_orig );
    if (eshells_orig != NULL) EG_free(eshells_orig);

    FREE(enodes   );
    FREE(ecurves  );
    FREE(eedges   );
    FREE(esurfaces);
    FREE(efaces   );
    FREE(eshells  );
    FREE(eloops   );
    FREE(efoo     );

    return status;
}
