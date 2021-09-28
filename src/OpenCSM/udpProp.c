/*
 ****************************************************************************
 *                                                                          *
 * udpProp -- udp file to generate a propeller                              *
 *                                                                          *
 *            Written by Daniel Oluwalana @ Syracuse University             *
 *                                                                          *
 * Calculations based on                                                    *
 * Adkins, C. N., & Liebeck, R. H. (1994). Design of optimum propellers.    *
 *               Journal of Propulsion and Power, 10(5), 676-682.           *
 *                                                                          *
 ****************************************************************************
 */

/*
 * Copyright (C) 2021  John F. Dannenhoffer, III (Syracuse University)
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

/* the number of arguments */
#define NUMUDPARGS 16

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define NBLADE( IUDP)   ((int    *) (udps[IUDP].arg[ 0].val))[0]
#define CPOWER( IUDP)   ((double *) (udps[IUDP].arg[ 1].val))[0]
#define LAMBDA( IUDP)   ((double *) (udps[IUDP].arg[ 2].val))[0]
#define REYR(   IUDP)   ((double *) (udps[IUDP].arg[ 3].val))[0]
#define RTIP(   IUDP)   ((double *) (udps[IUDP].arg[ 4].val))[0]
#define RHUB(   IUDP)   ((double *) (udps[IUDP].arg[ 5].val))[0]
#define CLIFT(  IUDP)   ((double *) (udps[IUDP].arg[ 6].val))[0]
#define CDRAG(  IUDP)   ((double *) (udps[IUDP].arg[ 7].val))[0]
#define ALFA(   IUDP)   ((double *) (udps[IUDP].arg[ 8].val))[0]

#define SHDIAM( IUDP)   ((double *) (udps[IUDP].arg[ 9].val))[0]
#define SHXMIN( IUDP)   ((double *) (udps[IUDP].arg[10].val))[0]
#define SHXMAX( IUDP)   ((double *) (udps[IUDP].arg[11].val))[0]
#define SPDIAM( IUDP)   ((double *) (udps[IUDP].arg[12].val))[0]
#define SPXMIN( IUDP)   ((double *) (udps[IUDP].arg[13].val))[0]

#define CTHRUST(IUDP)   ((double *) (udps[IUDP].arg[14].val))[0]
#define EFF(    IUDP)   ((double *) (udps[IUDP].arg[15].val))[0]

/* data about possible arguments */
static char*    argNames[NUMUDPARGS] = {"nblade",  "cpower",  "lambda",  "reyr",    "rtip",    "rhub",
                                        "clift",   "cdrag",   "alfa",
                                        "shdiam",  "shxmin",  "shxmax",  "spdiam",  "spxmin",
                                        "cthrust", "eff",                                               };
static int      argTypes[NUMUDPARGS] = {ATTRINT,   ATTRREAL,   ATTRREAL, ATTRREAL,  ATTRREAL,  ATTRREAL,
                                        ATTRREAL,  ATTRREAL,   ATTRREAL,
                                        ATTRREAL,  ATTRREAL,   ATTRREAL, ATTRREAL,  ATTRREAL,
                                        -ATTRREAL,-ATTRREAL                                             };
static int      argIdefs[NUMUDPARGS] = {2,         0,          0,        0,         0,         0,
                                        0,         0,          0,
                                        0,         0,          0,        0,         0,
                                        0,         0,                                                   };
static double   argDdefs[NUMUDPARGS] = {0.,        0.,         0.,       0.,        0.,        0.,
                                        0.,        0.,         0.,
                                        0.,        0.,         0.,       0.,        0.,
                                        0.,        0.,                                                  };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototypes for routine defined below */
static int    adkins(int iudp, int nsect, double r[], double chord[], double beta[], double *Tc);
static int    naca(double m, double p, double t, int npnt, double pnt[]);

// number of airfoil cross-sections
#define NSECT 13
// number of points around airfoil
#define NPNT  101

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

    int     ipnt, sizes[2], periodic, sense[3], isect, iblade;
    int     oclass, mtype, nchild, *senses;
    double  Tc, *radius=NULL, *chord=NULL, *beta=NULL;
    double  camber, locmax, thick, sparrad, ang, tiptreat[2];
    double  tle, tdata[4], dxytol=1e06, *pnt=NULL, result[3], data[18], xform[12], cosang, sinang;
    char    *message=NULL;
    ego     ecurve, eline, enodes[5], eedges[3], eloop, eairfoil, eshaft, exform, esects[NSECT+2], eblade;
    ego     etemp1, etemp2, emodel, eref, *echilds, eface;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("nblade(0) = %d\n", NBLADE(0));
    printf("cpower(0) = %f\n", CPOWER(0));
    printf("lambda(0) = %f\n", LAMBDA(0));
    printf("reyr(  0) = %f\n", REYR(  0));
    printf("rtip(  0) = %f\n", RTIP(  0));
    printf("rhub(  0) = %f\n", RHUB(  0));
    printf("clift( 0) = %f\n", CLIFT( 0));
    printf("cdrag( 0) = %f\n", CDRAG( 0));
    printf("alfa(  0) = %f\n", ALFA(  0));

    printf("shdiam(0) = %f\n", SHDIAM(0));
    printf("shxmin(0) = %f\n", SHXMIN(0));
    printf("shxmax(0) = %f\n", SHXMAX(0));
    printf("spdiam(0) = %f\n", SPDIAM(0));
    printf("spxmin(0) = %f\n", SPXMIN(0));
#endif

    /* default return values */
    *ebody = NULL;
    *nMesh = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (NBLADE(0) < 2) {
        snprintf(message, 100, "nblade=%d should be > 1", NBLADE(0));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("nblade(%d) = %d\n", numUdp, NBLADE(numUdp));
    printf("cpower(%d) = %f\n", numUdp, CPOWER(numUdp));
    printf("lambda(%d) = %f\n", numUdp, LAMBDA(numUdp));
    printf("reyr(  %d) = %f\n", numUdp, REYR(  numUdp));
    printf("rtip(  %d) = %f\n", numUdp, RTIP(  numUdp));
    printf("rhub(  %d) = %f\n", numUdp, RHUB(  numUdp));
    printf("clift( %d) = %f\n", numUdp, CLIFT( numUdp));
    printf("cdrag( %d) = %f\n", numUdp, CDRAG( numUdp));
    printf("alfa(  %d) = %f\n", numUdp, ALFA(  numUdp));

    printf("shdiam(%d) = %f\n", numUdp, SHDIAM(numUdp));
    printf("shxmin(%d) = %f\n", numUdp, SHXMIN(numUdp));
    printf("shxmax(%d) = %f\n", numUdp, SHXMAX(numUdp));
    printf("spdiam(%d) = %f\n", numUdp, SPDIAM(numUdp));
    printf("spxmin(%d) = %f\n", numUdp, SPXMIN(numUdp));
#endif

    /* perform design */
    MALLOC(radius, double, NSECT);
    MALLOC(chord,  double, NSECT);
    MALLOC(beta,   double, NSECT);

    status = adkins(numUdp, NSECT, radius, chord, beta, &Tc);
    CHECK_STATUS(adkins);

    /* since last chord is 0, make it half of the second-last chord */
    chord[NSECT-1] = 0.5 * chord[NSECT-2];

    printf("    radius      chord  beta(deg)\n");
    for (isect = 0; isect < NSECT; isect++) {
        printf("%10.5f %10.5f %10.5f\n", radius[isect], chord[isect], beta[isect]*180/PI);
    }

    /* make the airfoil cross-section  and an associated Loop */
    MALLOC(pnt, double, 3*NPNT);

    camber = 0.04;
    locmax = 0.40;
    thick  = 0.15;
    status = naca(camber, locmax, thick, NPNT, pnt);
    CHECK_STATUS(naca);

    /* create Node at upper trailing Edge */
    ipnt = 0;
    status = EG_makeTopology(context, NULL, NODE, 0, &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    /* create Node at leading Edge */
    ipnt = (NPNT - 1) / 2;

    status = EG_makeTopology(context, NULL, NODE, 0, &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create Node at lower trailing Edge */
    ipnt = NPNT - 1;
    status = EG_makeTopology(context, NULL, NODE, 0, &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[2]));
    CHECK_STATUS(EG_makeTopology);

    enodes[3] = enodes[0];

    /* create spline curve from upper TE, to LE, to lower TE */
    sizes[0] = NPNT;
    sizes[1] = 0;
    status = EG_approximate(context, 0, dxytol, sizes, pnt, &ecurve);
    CHECK_STATUS(EG_approximate);

    /* find t at leading Edge point */
    ipnt = (NPNT - 1) / 2;
    status = EG_invEvaluate(ecurve, &pnt[3*ipnt], &(tle), result);
    CHECK_STATUS(EG_invEvaluate);

    /* make Edge for upper surface */
    status = EG_getRange(ecurve, tdata, &periodic);
    CHECK_STATUS(EG_getRange);
    tdata[1] = tle;

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

    /* make Edge for lower surface */
    status = EG_getRange(ecurve, tdata, &periodic);
    CHECK_STATUS(EG_getRange);
    tdata[0] = tle;

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create line segment at trailing Edge */
    ipnt = NPNT - 1;
    data[0] = pnt[3*ipnt  ];
    data[1] = pnt[3*ipnt+1];
    data[2] = pnt[3*ipnt+2];
    data[3] = pnt[0] - data[0];
    data[4] = pnt[1] - data[1];
    data[5] = pnt[2] - data[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
    CHECK_STATUS(EG_makeGeometry);

    /* make Edge for this line */
    status = EG_invEvaluate(eline, data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    data[0] = pnt[0];
    data[1] = pnt[1];
    data[2] = pnt[2];
    status = EG_invEvaluate(eline, data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

    status = EG_makeTopology(context, eline, EDGE, TWONODE, tdata, 2, &(enodes[2]), NULL, &(eedges[2]));
    CHECK_STATUS(EG_makeTopology);

    /* create loop of the two or three Edges */
    sense[0] = SFORWARD;
    sense[1] = SFORWARD;
    sense[2] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 3, eedges, sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeFace(eloop, SFORWARD, NULL, &eairfoil);
    CHECK_STATUS(EG_makeFace);

    /* make the spar */
    sparrad = sqrt(thick / 2 / PI) * chord[0];     // roughly same area as airfoil
    
    for (ipnt = 0; ipnt < NPNT; ipnt++) {
        ang = 0.05 + (double)(ipnt)/(double)(NPNT-1) * (2*PI - 0.10);
        pnt[3*ipnt  ] = sparrad * cos(ang);
        pnt[3*ipnt+1] = sparrad * sin(ang);
        pnt[3*ipnt+2] = 0;
    }

    /* create Node at upper trailing Edge */
    ipnt = 0;
    status = EG_makeTopology(context, NULL, NODE, 0, &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    /* create Node at leading Edge */
    ipnt = (NPNT - 1) / 2;

    status = EG_makeTopology(context, NULL, NODE, 0, &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create Node at lower trailing Edge */
    ipnt = NPNT - 1;
    status = EG_makeTopology(context, NULL, NODE, 0, &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[2]));
    CHECK_STATUS(EG_makeTopology);

    enodes[3] = enodes[0];

    /* create spline curve from upper TE, to LE, to lower TE */
    sizes[0] = NPNT;
    sizes[1] = 0;
    status = EG_approximate(context, 0, dxytol, sizes, pnt, &ecurve);
    CHECK_STATUS(EG_approximate);

    /* find t at leading Edge point */
    ipnt = (NPNT - 1) / 2;
    status = EG_invEvaluate(ecurve, &pnt[3*ipnt], &(tle), result);
    CHECK_STATUS(EG_invEvaluate);

    /* make Edge for upper surface */
    status = EG_getRange(ecurve, tdata, &periodic);
    CHECK_STATUS(EG_getRange);
    tdata[1] = tle;

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

    /* make Edge for lower surface */
    status = EG_getRange(ecurve, tdata, &periodic);
    CHECK_STATUS(EG_getRange);
    tdata[0] = tle;

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create line segment at trailing Edge */
    ipnt = NPNT - 1;
    data[0] = pnt[3*ipnt  ];
    data[1] = pnt[3*ipnt+1];
    data[2] = pnt[3*ipnt+2];
    data[3] = pnt[0] - data[0];
    data[4] = pnt[1] - data[1];
    data[5] = pnt[2] - data[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
    CHECK_STATUS(EG_makeGeometry);

    /* make Edge for this line */
    status = EG_invEvaluate(eline, data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    data[0] = pnt[0];
    data[1] = pnt[1];
    data[2] = pnt[2];
    status = EG_invEvaluate(eline, data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

    status = EG_makeTopology(context, eline, EDGE, TWONODE, tdata, 2, &(enodes[2]), NULL, &(eedges[2]));
    CHECK_STATUS(EG_makeTopology);

    /* create loop of the two or three Edges */
    sense[0] = SFORWARD;
    sense[1] = SFORWARD;
    sense[2] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 3, eedges, sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeFace(eloop, SFORWARD, NULL, &eshaft);
    CHECK_STATUS(EG_makeFace);

    /* stack two circles at the axis and half-way to the hub */
    for (isect = 0; isect < 2; isect++) {

        /* make a copy of eshaft that has been transformed so that it is
           on the stacking line */
        cosang = cos(PI/2 - beta[0]);
        sinang = sin(PI/2 - beta[0]);

        xform[ 0] = +cosang;
        xform[ 1] = -sinang;
        xform[ 2] = 0;
        xform[ 3] = 0;

        xform[ 4] = +sinang;
        xform[ 5] = +cosang;
        xform[ 6] = 0;
        xform[ 7] = 0;

        xform[ 8] = 0;
        xform[ 9] = 0;
        xform[10] = 1;
        xform[11] = (double)(isect)/(double)(2) * radius[0];

        status = EG_makeTransform(context, xform, &exform);
        CHECK_STATUS(EG_makeTransform);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_copyObject(eshaft, exform, &(esects[isect]));
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(exform);
        CHECK_STATUS(EG_deleteObject);
    }

    /* stack airfoil, from hub to tip */
    for (isect = 0; isect < NSECT; isect++) {

        /* make a copy of eairfoil that has been transformed so that its quarter-chord
           is on the stacking line */
        cosang = cos(PI/2 - beta[isect]);
        sinang = sin(PI/2 - beta[isect]);

        xform[ 0] = +chord[isect] * cosang;
        xform[ 1] = -chord[isect] * sinang;
        xform[ 2] = 0;
        xform[ 3] = -chord[isect] * cosang / 4;

        xform[ 4] = +chord[isect] * sinang;
        xform[ 5] = +chord[isect] * cosang;
        xform[ 6] = 0;
        xform[ 7] = -chord[isect] * sinang / 4;

        xform[ 8] = 0;
        xform[ 9] = 0;
        xform[10] =  chord[isect];
        xform[11] =  radius[isect];

        status = EG_makeTransform(context, xform, &exform);
        CHECK_STATUS(EG_makeTransform);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_copyObject(eairfoil, exform, &(esects[isect+2]));
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(exform);
        CHECK_STATUS(EG_deleteObject);
    }

    /* now make the blend */
    tiptreat[0] = 0;
    tiptreat[1] = 4;

    status = EG_blend(NSECT+2, esects, NULL, tiptreat, &eblade);
    CHECK_STATUS(EG_blend);

    /* make the shaft */
    if (SHDIAM(numUdp) > 0) {
        data[0] = SHXMIN(numUdp);
        data[1] = 0;
        data[2] = 0;
        data[3] = SHXMAX(numUdp);
        data[4] = 0;
        data[5] = 0;
        data[6] = SHDIAM(numUdp) / 2;

        status = EG_makeSolidBody(context, CYLINDER, data, &etemp1);
        CHECK_STATUS(EG_makeSolidBody);
    } else {
        etemp1 = NULL;
    }

    /* now assemble NBLADEs together */
    for (iblade = 0; iblade < NBLADE(numUdp); iblade++) {
        ang = 2 * PI * (double)(iblade) / (double)(NBLADE(numUdp));

        xform[ 0] = 1;   xform[ 1] = 0;           xform[ 2] = 0;           xform[ 3] = 0;
        xform[ 4] = 0;   xform[ 5] = +cos(ang);   xform[ 6] = -sin(ang);   xform[ 7] = 0;
        xform [8] = 0;   xform[ 9] = +sin(ang);   xform[10] = +cos(ang);   xform[11] = 0;

        status = EG_makeTransform(context, xform, &exform);
        CHECK_STATUS(EG_makeTransform);

        status = EG_copyObject(eblade, exform, &etemp2);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(exform);
        CHECK_STATUS(EG_deleteObject);

        if (etemp1 != NULL) {
            status = EG_generalBoolean(etemp1, etemp2, FUSION, 0, &emodel);
            CHECK_STATUS(EG_generalBoolean);

            status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                                    &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_copyObject(echilds[0], NULL, &etemp1);
            CHECK_STATUS(EG_copyObject);

            status = EG_deleteObject(emodel);
            CHECK_STATUS(EG_deleteObject);
        } else {
            status = EG_copyObject(eblade, NULL, &etemp1);
            CHECK_STATUS(EG_deleteObject);
        }
    }

    if (etemp1 == NULL) {
        status = -999;
        goto cleanup;
    }

    /* possibly add in a spinner */
    if (SHDIAM(numUdp) > 0 && SPDIAM(numUdp) > SHDIAM(numUdp) && SPXMIN(numUdp) < SHXMIN(numUdp)) {

        /* centerline */
        data[0] = SPXMIN(numUdp);
        data[1] = 0;
        data[2] = 0;
        status = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        data[0] = SHXMIN(numUdp);
        data[1] = 0;
        data[2] = 0;
        status = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        data[0] = SPXMIN(numUdp);
        data[1] = 0;
        data[2] = 0;
        data[3] = SHXMIN(numUdp) - SPXMIN(numUdp);
        data[4] = 0;
        data[5] = 0;
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        data[0]  = 0;
        data[1]  = SHXMIN(numUdp) - SPXMIN(numUdp);
        sense[0] = SFORWARD;
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, data, 2, &(enodes[0]), sense, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* downstream */
        data[0] = SHXMIN(numUdp);
        data[1] = SPDIAM(numUdp) / 2;
        data[2] = 0;
        status = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL, &(enodes[2]));
        CHECK_STATUS(EG_makeTopology);

        data[0] = SHXMIN(numUdp);
        data[1] = 0;
        data[2] = 0;
        data[3] = 0;
        data[4] = SPDIAM(numUdp) / 2;
        data[5] = 0;
        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        data[0]  = 0;
        data[1]  = SPDIAM(numUdp) / 2;
        sense[0] = SFORWARD;
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, data, 2, &(enodes[1]), sense, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        /* parabola */
        enodes[3] = enodes[0];
        enodes[4] = enodes[2];

        data[0] = SPXMIN(numUdp);   data[1] = 0;   data[2] = 0;
        data[3] = 1;                data[4] = 0;   data[5] = 0;
        data[6] = 0;                data[7] = 1;   data[8] = 0;
        data[9] = SPDIAM(numUdp) * SPDIAM(numUdp) / (SHXMIN(numUdp) - SPXMIN(numUdp)) / 16;
        status = EG_makeGeometry(context, CURVE, PARABOLA, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        tdata[0] = 0;
        tdata[1] = SPDIAM(numUdp) / 2;

        sense[0] = SFORWARD;
        status = EG_makeTopology(context, ecurve, EDGE, TWONODE, tdata, 2, &(enodes[3]), sense, &(eedges[2]));
        CHECK_STATUS(EG_makeTopology);

        /* make Loop, Face, and Body */
        status = EG_makeLoop(3, eedges, NULL, 0, &eloop);
        CHECK_STATUS(EG_makeLoop);

        status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
        CHECK_STATUS(EG_makeFace);

        /* add in first half body of revolution */
        data[0] = 0;
        data[1] = 0;
        data[2] = 0;
        data[3] = 1;
        data[4] = 0;
        data[5] = 0;
        status = EG_rotate(eface, +180, data, &etemp2);
        CHECK_STATUS(EG_rotate);

        status = EG_generalBoolean(etemp1, etemp2, FUSION, 0, &emodel);
        CHECK_STATUS(EG_generalBoolean);

        status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_copyObject(echilds[0], NULL, &etemp1);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(emodel);
        CHECK_STATUS(EG_deleteObject);

        /* add in second half body of revolution */
        status = EG_rotate(eface, -180, data, &etemp2);
        CHECK_STATUS(EG_rotate);

        status = EG_generalBoolean(etemp1, etemp2, FUSION, 0, &emodel);
        CHECK_STATUS(EG_generalBoolean);

        status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_copyObject(echilds[0], NULL, &etemp1);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(emodel);
        CHECK_STATUS(EG_deleteObject);
    }

    *ebody = etemp1;

    /* set the output value(s) */
    CTHRUST(0) = Tc;
    EFF(    0) = Tc / CPOWER(0);

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(radius);
    FREE(chord );
    FREE(beta  );

    FREE(pnt);

    if (pnt != NULL) EG_free(pnt);

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
        status = EGADS_NOTMODEL;
        goto cleanup;
    }

    status = EGADS_SUCCESS;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   adkins - calculation procedure described in:                       *
 *            Adkins, C.N., and Liebeck, R.H.,                          *
 *            "Design of Optimum Propellers"                            *
 *            JPP, Vol 10, No 5, 1994, pp 676--682                      *
 *                                                                      *
 ************************************************************************
 */

static int
adkins(int     iudp,                    /* (in)  iudp index */
       int     nsect,                   /* (in)  number of sections */
       double  radius[],                /* (out) radial position */
       double  chord[],                 /* (out) chord distribution */
       double  beta[],                  /* (out) twist angle distribution */
       double  *Tc)                     /* (out) thrust coefficient */
{
    int    status = EGADS_SUCCESS;

    int    iter, isect;
    double zeta, phi_t, I1, I2, J1, J2, I1prime, I2prime, J1prime, J2prime, zeta_new;
    double xi, f, F, phi, G, Wc_V, eps, a, W_V;

    ROUTINE(adkins);

    /* --------------------------------------------------------------- */

    /* initial guess for displacement velocity ratio */
    zeta = 0;

    /* iteration loop for zeta */
    for (iter = 0; iter < 10; iter++) {

#ifdef DEBUG
        printf("zeta=%10.5f\n", zeta);
#endif

        /* initialize integrals */
        I1 = 0;
        I2 = 0;
        J1 = 0;
        J2 = 0;

        /* equation 20 */
        phi_t = atan(LAMBDA(iudp) * (1 + zeta/2));

        /* loop over the sections */
        for (isect = 0; isect < nsect; isect++) {

            /* evenly space radial positions */
            radius[isect] = RHUB(iudp) + (double)(isect)/(double)(nsect-1) * (RTIP(iudp) - RHUB(iudp));

            /* definition */
            xi = radius[isect] / RTIP(iudp);

            /* equation 19 */
            f = NBLADE(iudp) / 2 * (1 - xi) / sin(phi_t);
            if (f < 0) f = 0;

            /* equation 18 */
            F = 2 / PI * acos(exp(-f));

            /* equation 21 */
            phi = atan(tan(phi_t) / xi);

            /* equation 5 */
            G = F * xi / LAMBDA(iudp) * cos(phi) * sin(phi);

            /* equation 16 */
            Wc_V = 4 * PI * LAMBDA(iudp) * G * RTIP(iudp) * zeta / (CLIFT(iudp) * NBLADE(iudp));

            /* definition */
//$$$            Re = Wc_V * REYR(iudp) / RTIP(iudp);

            /* airfoil data */
            eps = CDRAG(iudp) / CLIFT(iudp);

            /* equation 7a */
            a = zeta / 2 * cos(phi)*cos(phi) * (1 - eps * tan(phi));

            /* equation 7b */
//$$$            aprime = zeta / (2*xi/LAMBDA(iudp)) * cos(phi) * sin(phi) * (1 + eps/tan(phi));

            /* equation 17 */
            W_V = (1 + a) / sin(phi);

            /* chord */
            chord[isect] = Wc_V / W_V;

            /* definition */
            beta[isect] = ALFA(iudp)*PI/180 + phi;

            /* equation 11a */
            I1prime = 4 * xi * G * (1 - eps * tan(phi));

            /* equation 11b */
            I2prime = LAMBDA(iudp) * (I1prime/2/xi) * (1 + eps/tan(phi)) * sin(phi) * cos(phi);

            /* equation 11c */
            J1prime = 4 * xi * G * (1 + eps/tan(phi));

            /* equation 11d */
            J2prime = J1prime/2 * (1 - eps*tan(phi)) * cos(phi) * cos(phi);

            /* trapezoidal integration */
            if (isect > 0 && isect < nsect-1) {
                I1 += I1prime;
                I2 += I2prime;
                J1 += J1prime;
                J2 += J2prime;
            } else {
                I1 += I1prime / 2;
                I2 += I2prime / 2;
                J1 += J1prime / 2;
                J2 += J2prime / 2;
            }
        }

        /* normalize the integrals */
        I1 *= (1 - RHUB(iudp)/RTIP(iudp)) / (nsect - 1);
        I2 *= (1 - RHUB(iudp)/RTIP(iudp)) / (nsect - 1);
        J1 *= (1 - RHUB(iudp)/RTIP(iudp)) / (nsect - 1);
        J2 *= (1 - RHUB(iudp)/RTIP(iudp)) / (nsect - 1);

        /* equation 14 */
        zeta_new = - (J1 / 2 / J2) + sqrt((J1/2/J2)*(J1/2/J2) + CPOWER(iudp)/J2);

        /* equation 15 */
        *Tc = (I1 - I2 * zeta_new) * zeta_new;

        /* convergence check */
        if (fabs(zeta_new-zeta) < 1e-5) {
            goto cleanup;
        }

        /* update zeta for next iteration */
        zeta = zeta_new;
    }

    /* getting here means that we did not converge */
    status = EGADS_RANGERR;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   naca - generate points arount NACA 4-series airfoil                *
 *                                                                      *
 ************************************************************************
 */

static int
naca(double m,                          /* (in)  maximum camber/chord */
     double p,                          /* (in)  location of max camber/chord */
     double t,                          /* (in)  maximum thickness to chord */
     int    npnt,                       /* (in)  number of points (TE-upper-LE-lower-TE) */
     double pnt[])                      /* (out) x-, y-, z-coordinates */
{
    int     status = EGADS_SUCCESS;

    int     ipnt;
    double  zeta, s, yt, yc, theta;

    ROUTINE(naca);

    /* --------------------------------------------------------------- */

    /* points around airfoil (upper and then lower) */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        zeta = TWOPI * ipnt / (npnt-1);
        s    = (1 + cos(zeta)) / 2;
        yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));

        if (s < p) {
            yc    =      m / (  p)/(  p) * (s * (2*p -   s));
            theta = atan(m / (  p)/(  p) * (    (2*p - 2*s)));
        } else {
            yc    =      m / (1-p)/(1-p) * ((1-2*p) + s * (2*p -   s));
            theta = atan(m / (1-p)/(1-p) * (              (2*p - 2*s)));
        }

        if (ipnt < npnt/2) {
            pnt[3*ipnt  ] = s  - yt * sin(theta);
            pnt[3*ipnt+1] = yc + yt * cos(theta);
            pnt[3*ipnt+2] = 0;
        } else if (ipnt == npnt/2) {
            pnt[3*ipnt  ] = 0;
            pnt[3*ipnt+1] = 0;
            pnt[3*ipnt+2] = 0;
        } else {
            pnt[3*ipnt  ] = s  + yt * sin(theta);
            pnt[3*ipnt+1] = yc - yt * cos(theta);
            pnt[3*ipnt+2] = 0;
        }
    }

//cleanup:
    return status;
}
