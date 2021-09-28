/*
 ************************************************************************
 *                                                                      *
 * udpTire -- udp file to generate a tire                               *
 *                                                                      *
 *            Written by Jake Miller @ Syracuse University              *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                      John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2020  John F. Dannenhoffer, III (Syracuse University)
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
#ifdef UDP
   #define WIDTH(IUDP)      ((double *) (udps[IUDP].arg[0].val))[0]
   #define MINRAD(IUDP)     ((double *) (udps[IUDP].arg[1].val))[0]
   #define MAXRAD(IUDP)     ((double *) (udps[IUDP].arg[2].val))[0]
   #define FILLETRAD(IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
   #define PLATETHICK(IUDP) ((double *) (udps[IUDP].arg[4].val))[0]
   #define PATTERN(IUDP)    ((double *) (udps[IUDP].arg[5].val))[0]
   #define BOLTS(IUDP)      ((double *) (udps[IUDP].arg[6].val))[0]
   #define BOLTRAD(IUDP)    ((double *) (udps[IUDP].arg[7].val))[0]
   #define VOLUME(IUDP)     ((double *) (udps[IUDP].arg[8].val))[0]

   /* data about possible arguments */
   static char  *argNames[NUMUDPARGS] = {"width",  "minrad", "maxrad", "fillrad", "platethick", "patternrad", "bolts",   "boltrad", "volume", };
   static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL,  ATTRREAL,     ATTRREAL,      ATTRREAL, ATTRREAL,  -ATTRREAL,};
   static int    argIdefs[NUMUDPARGS] = {0,        0,        0,        0,         0,            0,             0,        0,         0,        };
   static double argDdefs[NUMUDPARGS] = {0.,       0.,       0.,       0.,        0.,           0.,            0.,       0.,        0.,       };

   /* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                            udpGet, udpVel, udpClean, udpMesh */
   #include "udpUtilities.c"
#else
   #define WIDTH(IUDP)      5.0
   #define MINRAD(IUDP)     8.0
   #define MAXRAD(IUDP)    12.0
   #define FILLETRAD(IUDP)  2.0
   #define PLATETHICK(IUDP) 0.5
   #define PATTERN(IUDP)    4.0
   #define BOLTS(IUDP)      5.0
   #define BOLTRAD(IUDP)    1.0
   #define VOLUME(IUDP)     myVolume

/*@null@*/ static char *
udpErrorStr(int stat)                   /* (in)  status number */
{
    char *string;                       /* (out) error message */

    string = EG_alloc(25*sizeof(char));
    if (string == NULL) {
        return string;
    }
    snprintf(string, 25, "EGADS status = %d", stat);

    return string;
}

void printEgo(/*@unused@*/ /*@null@*/ego obj) {};

int
main(/*@unused@*/ int  argc,                    /* (in)  number of arguments */
     /*@unused@*/ char *argv[])                 /* (in)  array of arguments */
{
    int  status, nMesh;
    char *string;
    ego  context, ebody, emodel;

    /* dummy call to prevent compiler warnings */
    printEgo(NULL);

    /* define a context */
    status = EG_open(&context);
    printf("EG_open -> status=%d\n", status);
    if (status < 0) exit(EXIT_FAILURE);

    /* call the execute routine */
    status = udpExecute(context, &ebody, &nMesh, &string);
    printf("udpExecute -> status=%d\n", status);
    if (status < 0) exit(EXIT_FAILURE);

    EG_free(string);

    /* make and dump the model */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL,
                             1, &ebody, NULL, &emodel);
    printf("EG_makeTopology -> status=%d\n", status);
    if (status < 0) exit(EXIT_FAILURE);

    status = EG_saveModel(emodel, "tire.egads");
    printf("EG_saveModel -> status=%d\n", status);
    if (status < 0) exit(EXIT_FAILURE);

    /* cleanup */
    status = EG_deleteObject(emodel);
    printf("EG_close -> status=%d\n", status);

    status = EG_close(context);
    printf("EG_close -> status=%d\n", status);

    return EXIT_SUCCESS;
}
#endif


/*
 ************************************************************************
 *                                                                      *
 *   periodicSeam - build a PCurve along a periodic seam                *
 *                                                                      *
 ************************************************************************
 */

static int
periodicSeam(ego eedge,                 /* (in)  Edge associated with seam */
             int sense,                 /* (in)  sense of PCurve */
             ego *pcurve)               /* (out) PCurve created */
{
    int    status = EGADS_SUCCESS;
    int    per;
    double range[2], data[4];
    ego    context;

    status = EG_getContext(eedge, &context);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_getRange(eedge, range, &per);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set up u and v at ends */
    data[0] = TWOPI;
    data[1] = range[0];
    data[2] = 0;
    data[3] = sense;
    if (sense == -1) data[1] = range[1];

    /* make (linear) PCurve */
    status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, pcurve);

cleanup:
    return status;
}


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
    int     sense[20], oclass, mtype, nchild, *senses, i;
    double  node1[3], node2[3], node3[3], node4[3], node5[3], node6[3], node7[3], node8[3];
    double  cent1[3], cent2[3], axis1[3], axis2[3], axis3[3];
    double  data[18], trange[2];
    ego     enodes[8], ecurve[16], eedges[16], eloop, efaces[8], eshell;
    ego     esurface[4], epcurve[4], ebody1, ebody2, ebody3, ebody4;
    ego     elist[20], emodel, *echilds2, source, *echilds, eref;

#ifndef UDP
    double  myVolume;
#endif

#ifdef DEBUG
    printf("udpExecute(context=%lx)\n", (long)context);
    printf("width(0)       = %f\n", WIDTH(    0));
    printf("minrad(0)      = %f\n", MINRAD(   0));
    printf("maxrad(0)      = %f\n", MAXRAD(   0));
    printf("fillrad(0)     = %f\n", FILLETRAD(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

#ifdef UDP
    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        printf(" udpExecute: width should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (WIDTH(0) <= 0) {
        printf(" udpExecute: width = %f < 0\n", WIDTH(0));
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[1].size > 1) {
        printf(" udpExecute: minrad should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (MINRAD(0) <= 0) {
        printf(" udpExecute: minrad = %f < 0\n", MINRAD(0));
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[2].size > 1) {
        printf(" udpExecute: maxrad should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (MAXRAD(0) <= 0) {
        printf(" udpExecute: maxrad = %f < 0\n", MAXRAD(0));
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (FILLETRAD(0) < 0) {
        printf(" udpExecute: width = %f < 0\n", WIDTH(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[3].size > 1) {
        printf(" udpExecute: minrad should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (PLATETHICK(0) < 0) {
        printf(" udpExecute: width = %f < 0\n", WIDTH(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[4].size > 1) {
        printf(" udpExecute: minrad should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (PATTERN(0) < 0) {
        printf(" udpExecute: width = %f < 0\n", WIDTH(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[5].size > 1) {
        printf(" udpExecute: minrad should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (BOLTS(0) < 0) {
        printf(" udpExecute: width = %f < 0\n", WIDTH(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[6].size > 1) {
        printf(" udpExecute: minrad should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (BOLTRAD(0) < 0) {
        printf(" udpExecute: width = %f < 0\n", WIDTH(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[7].size > 1) {
        printf(" udpExecute: minrad should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (PATTERN(0) > MINRAD(0)) {
        printf(" udpExecute: patternrad must be less than minrad");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (WIDTH(0) <= 0 && MINRAD(0) <= 0 && MAXRAD(0) <= 0) {
        printf(" udpExecute: width=minrad=maxrad=0\n");
        status  = EGADS_GEOMERR;
        goto cleanup;
    } else if (MINRAD(0) > MAXRAD(0)) {
        printf(" udpExecute: minrad cannot be bigger than maxrad\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    if (status < 0) {
        printf(" udpExecute: problem caching arguments\n");
        goto cleanup;
    }
#endif

#ifdef DEBUG
    printf("width[     %d] = %f\n", numUdp, WIDTH(     numUdp));
    printf("minrad[    %d] = %f\n", numUdp, MINRAD(    numUdp));
    printf("maxrad[    %d] = %f\n", numUdp, MAXRAD(    numUdp));
    printf("fillrad[   %d] = %f\n", numUdp, FILLETRAD( numUdp));
    printf("platethick[%d] = %f\n", numUdp, PLATETHICK(numUdp));
    printf("bolts[     %d] = %f\n", numUdp, BOLTS(     numUdp));
    printf("patternrad[%d] = %f\n", numUdp, PATTERN(   numUdp));
    printf("boltrad[   %d] = %f\n", numUdp, BOLTRAD(   numUdp));
#endif

    /* Node locations */
    node1[0] = -MINRAD(0);  node1[1] = 0;  node1[2] = -WIDTH(0) / 2;
    node2[0] = -MINRAD(0);  node2[1] = 0;  node2[2] =  WIDTH(0) / 2;
    node3[0] = -MAXRAD(0);  node3[1] = 0;  node3[2] =  WIDTH(0) / 2;
    node4[0] = -MAXRAD(0);  node4[1] = 0;  node4[2] = -WIDTH(0) / 2;
    node5[0] =  MINRAD(0);  node5[1] = 0;  node5[2] = -WIDTH(0) / 2;
    node6[0] =  MAXRAD(0);  node6[1] = 0;  node6[2] = -WIDTH(0) / 2;
    node7[0] =  MAXRAD(0);  node7[1] = 0;  node7[2] =  WIDTH(0) / 2;
    node8[0] =  MINRAD(0);  node8[1] = 0;  node8[2] =  WIDTH(0) / 2;

    /* make Nodes */
    status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &(enodes[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node2, 0, NULL, NULL, &(enodes[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node3, 0, NULL, NULL, &(enodes[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node4, 0, NULL, NULL, &(enodes[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node5, 0, NULL, NULL, &(enodes[4]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node6, 0, NULL, NULL, &(enodes[5]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node7, 0, NULL, NULL, &(enodes[6]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, NODE, 0, node8, 0, NULL, NULL, &(enodes[7]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 1 */
    data[0] = node1[0];
    data[1] = node1[1];
    data[2] = node1[2];
    data[3] = node2[0] - node1[0];
    data[4] = node2[1] - node1[1];
    data[5] = node2[2] - node1[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[0], node1, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[0], node2, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[0];
    elist[1] = enodes[1];
    status = EG_makeTopology(context, ecurve[0], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 2 */
    data[0] = node2[0];
    data[1] = node2[1];
    data[2] = node2[2];
    data[3] = node3[0] - node2[0];
    data[4] = node3[1] - node2[1];
    data[5] = node3[2] - node2[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[1], node2, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[1], node3, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[1];
    elist[1] = enodes[2];
    status = EG_makeTopology(context, ecurve[1], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 3 */
    data[0] = node3[0];
    data[1] = node3[1];
    data[2] = node3[2];
    data[3] = node4[0] - node3[0];
    data[4] = node4[1] - node3[1];
    data[5] = node4[2] - node3[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[2], node3, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[2], node4, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[2];
    elist[1] = enodes[3];
    status = EG_makeTopology(context, ecurve[2], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 4 */
    data[0] = node4[0];
    data[1] = node4[1];
    data[2] = node4[2];
    data[3] = node1[0] - node4[0];
    data[4] = node1[1] - node4[1];
    data[5] = node1[2] - node4[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[3], node4, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[3], node1, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[3];
    elist[1] = enodes[0];
    status = EG_makeTopology(context, ecurve[3], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 5 */
    data[0] = node5[0];
    data[1] = node5[1];
    data[2] = node5[2];
    data[3] = node6[0] - node5[0];
    data[4] = node6[1] - node5[1];
    data[5] = node6[2] - node5[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[4]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[4], node5, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[4], node6, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[4];
    elist[1] = enodes[5];
    status = EG_makeTopology(context, ecurve[4], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[4]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 6 */
    data[0] = node6[0];
    data[1] = node6[1];
    data[2] = node6[2];
    data[3] = node7[0] - node6[0];
    data[4] = node7[1] - node6[1];
    data[5] = node7[2] - node6[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[5]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[5], node6, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[5], node7, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[5];
    elist[1] = enodes[6];
    status = EG_makeTopology(context, ecurve[5], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[5]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 7 */
    data[0] = node7[0];
    data[1] = node7[1];
    data[2] = node7[2];
    data[3] = node8[0] - node7[0];
    data[4] = node8[1] - node7[1];
    data[5] = node8[2] - node7[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[6]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[6], node7, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[6], node8, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[6];
    elist[1] = enodes[7];
    status = EG_makeTopology(context, ecurve[6], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[6]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (linear) Edge 8 */
    data[0] = node8[0];
    data[1] = node8[1];
    data[2] = node8[2];
    data[3] = node5[0] - node8[0];
    data[4] = node5[1] - node8[1];
    data[5] = node5[2] - node8[2];
    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurve[7]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[7], node8, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[7], node5, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    elist[0] = enodes[7];
    elist[1] = enodes[4];
    status = EG_makeTopology(context, ecurve[7], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[7]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* data used in creating the arcs */
    axis1[0] = 1;   axis1[1] = 0;   axis1[2] = 0;
    axis2[0] = 0;   axis2[1] = 1;   axis2[2] = 0;
    axis3[0] = 0;   axis3[1] = 0;   axis3[2] = 1;

    cent1[0] = 0;   cent1[1] = 0;   cent1[2] = -WIDTH(0) / 2;
    cent2[0] = 0;   cent2[1] = 0;   cent2[2] =  WIDTH(0) / 2;

    /* make (circular) Edge 9 */
    data[0] = cent1[0];   data[1] = cent1[1];   data[2] = cent1[2];
    data[3] = axis1[0];   data[4] = axis1[1];   data[5] = axis1[2];
    data[6] = axis2[0];   data[7] = axis2[1];   data[8] = axis2[2];   data[9] = MINRAD(0);
    status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &(ecurve[8]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[8], node5, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[8], node1, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[4];
    elist[1] = enodes[0];
    status = EG_makeTopology(context, ecurve[8], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[8]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 10 */
    data[0] = cent2[0];   data[1] = cent2[1];   data[2] = cent2[2];
    data[3] = axis1[0];   data[4] = axis1[1];   data[5] = axis1[2];
    data[6] = axis2[0];   data[7] = axis2[1];   data[8] = axis2[2];   data[9] = MINRAD(0);
    status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &(ecurve[9]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[9], node8, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[9], node2, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[7];
    elist[1] = enodes[1];
    status = EG_makeTopology(context, ecurve[9], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[9]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 11 */
    data[0] = cent1[0];   data[1] = cent1[1];   data[2] = cent1[2];
    data[3] = axis1[0];   data[4] = axis1[1];   data[5] = axis1[2];
    data[6] = axis2[0];   data[7] = axis2[1];   data[8] = axis2[2];   data[9] = MAXRAD(0);
    status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &(ecurve[10]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[10], node6, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[10], node4, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[5];
    elist[1] = enodes[3];
    status = EG_makeTopology(context, ecurve[10], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[10]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 12 */
    data[0] = cent2[0];   data[1] = cent2[1];   data[2] = cent2[2];
    data[3] = axis1[0];   data[4] = axis1[1];   data[5] = axis1[2];
    data[6] = axis2[0];   data[7] = axis2[1];   data[8] = axis2[2];   data[9] = MAXRAD(0);
    status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, data, &(ecurve[11]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[11], node7, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[11], node3, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[6];
    elist[1] = enodes[2];
    status = EG_makeTopology(context, ecurve[11], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[11]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the outer cylindrical surface */
    data[0]  = cent1[0];   data[1]  = cent1[1];   data[2]  = cent1[2];
    data[3]  = axis1[0];   data[4]  = axis1[1];   data[5]  = axis1[2];
    data[6]  = axis2[0];   data[7]  = axis2[1];   data[8]  = axis2[2];
    data[9]  = axis3[0];   data[10] = axis3[1];   data[11] = axis3[2];   data[12] = MAXRAD(0);
    status = EG_makeGeometry(context, SURFACE, CYLINDRICAL, NULL, NULL, data, &(esurface[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeGeometry(context, SURFACE, CYLINDRICAL, NULL, NULL, data, &(esurface[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the inner cylindrical surface */
    data[0]  = cent1[0];   data[1]  = cent1[1];   data[2]  = cent1[2];
    data[3]  = axis1[0];   data[4]  = axis1[1];   data[5]  = axis1[2];
    data[6]  = axis2[0];   data[7]  = axis2[1];   data[8]  = axis2[2];
    data[9]  = axis3[0];   data[10] = axis3[1];   data[11] = axis3[2];   data[12] = MINRAD(0);
    status = EG_makeGeometry(context, SURFACE, CYLINDRICAL, NULL, NULL, data, &(esurface[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeGeometry(context, SURFACE, CYLINDRICAL, NULL, NULL, data, &(esurface[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (planar) Face 1 */
    sense[0] = SFORWARD;    sense[1] = SREVERSE;    sense[2] = SFORWARD;    sense[3] = SFORWARD;
    elist[0] = eedges[3];   elist[1] = eedges[8];   elist[2] = eedges[4];   elist[3] = eedges[10];
    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeFace(eloop, SFORWARD, NULL, &(efaces[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (planar) Face 2 */
    sense[0] = SFORWARD;    sense[1] = SREVERSE;    sense[2] = SFORWARD;    sense[3] = SFORWARD;
    elist[0] = eedges[1];   elist[1] = eedges[11];  elist[2] = eedges[6];   elist[3] = eedges[9];
    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeFace(eloop, SFORWARD, NULL, &(efaces[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (cylindrical) Face 3 */
    status = EG_otherCurve(esurface[0], ecurve[2], 0, &(epcurve[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[0], ecurve[10], 0, &(epcurve[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[0], ecurve[5], 0, &(epcurve[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[0], ecurve[11], 0, &(epcurve[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    sense[0] = SFORWARD;    sense[1] = SREVERSE;    sense[2] = SFORWARD;    sense[3] = SFORWARD;
    elist[0] = eedges[2];   elist[1] = eedges[10];  elist[2] = eedges[5];   elist[3] = eedges[11];
    elist[4] = epcurve[0];  elist[5] = epcurve[1];  elist[6] = epcurve[2];  elist[7] = epcurve[3];
    status = EG_makeTopology(context, esurface[0], LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, esurface[0], FACE, SREVERSE, NULL, 1, &eloop, sense, &(efaces[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (cylindrical) Face 4 */
    status = EG_otherCurve(esurface[1], ecurve[0], 0, &(epcurve[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[1], ecurve[9], 0, &(epcurve[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[1], ecurve[7], 0, &(epcurve[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[1], ecurve[8], 0, &(epcurve[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    sense[0] = SFORWARD;    sense[1] = SREVERSE;    sense[2] = SFORWARD;    sense[3] = SFORWARD;
    elist[0] = eedges[0];   elist[1] = eedges[9];   elist[2] = eedges[7];   elist[3] = eedges[8];
    elist[4] = epcurve[0];  elist[5] = epcurve[1];  elist[6] = epcurve[2];  elist[7] = epcurve[3];
    status = EG_makeTopology(context, esurface[1], LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, esurface[1], FACE, SFORWARD, NULL, 1, &eloop, sense, &(efaces[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 13 */
    ecurve[12] = ecurve[8];   /* reuse */

    status = EG_invEvaluate(ecurve[12], node1, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[12], node5, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[0];
    elist[1] = enodes[4];
    status = EG_makeTopology(context, ecurve[12], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[12]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 14 */
    ecurve[13] = ecurve[9];   /* reuse */

    status = EG_invEvaluate(ecurve[13], node2, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[13], node8, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[1];
    elist[1] = enodes[7];
    status = EG_makeTopology(context, ecurve[13], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[13]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 15 */
    ecurve[14] = ecurve[10];

    status = EG_invEvaluate(ecurve[14], node4, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[14], node6, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[3];
    elist[1] = enodes[5];
    status = EG_makeTopology(context, ecurve[14], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[14]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (circular) Edge 16 */
    ecurve[15] = ecurve[11];

    status = EG_invEvaluate(ecurve[15], node3, &(trange[0]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_invEvaluate(ecurve[15], node7, &(trange[1]), data);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (trange[0] > trange[1]) trange[1] += TWOPI;   /* ensure trange[1] > trange[0] */

    elist[0] = enodes[2];
    elist[1] = enodes[6];
    status = EG_makeTopology(context, ecurve[15], EDGE, TWONODE, trange, 2, elist, NULL, &(eedges[15]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (planar) Face 5 */
    sense[0] = SFORWARD;    sense[1] = SREVERSE;    sense[2] = SREVERSE;    sense[3] = SREVERSE;
    elist[0] = eedges[14];  elist[1] = eedges[4];   elist[2] = eedges[12];  elist[3] = eedges[3];
    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeFace(eloop, SFORWARD, NULL, &(efaces[4]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (planar) Face 6 */
    sense[0] = SFORWARD;    sense[1] = SREVERSE;    sense[2] = SFORWARD;    sense[3] = SFORWARD;
    elist[0] = eedges[6];   elist[1] = eedges[13];  elist[2] = eedges[1];   elist[3] = eedges[15];
    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeFace(eloop, SFORWARD, NULL, &(efaces[5]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (cylindrical) Face 7 */
    status = EG_otherCurve(esurface[2], ecurve[2], 0, &(epcurve[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[2], ecurve[14], 0, &(epcurve[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* EG_otherCurve cannot be used along periodic seam */
//  status = EG_otherCurve(esurface[2], ecurve[5], 0, &(epcurve[2]));
    status = periodicSeam(eedges[5], SFORWARD,        &(epcurve[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[2], ecurve[15], 0, &(epcurve[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    sense[0] = SFORWARD;    sense[1] = SFORWARD;    sense[2] = SFORWARD;    sense[3] = SREVERSE;
    elist[0] = eedges[2];   elist[1] = eedges[14];  elist[2] = eedges[5];   elist[3] = eedges[15];
    elist[4] = epcurve[0];  elist[5] = epcurve[1];  elist[6] = epcurve[2];  elist[7] = epcurve[3];
    status = EG_makeTopology(context, esurface[2], LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, esurface[2], FACE, SFORWARD, NULL, 1, &eloop, sense, &(efaces[6]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make (cylindrical) Face 8 */
    status = EG_otherCurve(esurface[3], ecurve[0], 0, &(epcurve[0]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[3], ecurve[13], 0, &(epcurve[1]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* EG_otherCurve cannot be used along periodic seam */
//  status = EG_otherCurve(esurface[3], ecurve[7], 0, &(epcurve[2]));
    status = periodicSeam(eedges[7], SREVERSE,        &(epcurve[2]));
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_otherCurve(esurface[3], ecurve[12], 0, &(epcurve[3]));
    if (status != EGADS_SUCCESS) goto cleanup;

    sense[0] = SFORWARD;    sense[1] = SFORWARD;    sense[2] = SFORWARD;    sense[3] = SREVERSE;
    elist[0] = eedges[0];   elist[1] = eedges[13];  elist[2] = eedges[7];   elist[3] = eedges[12];
    elist[4] = epcurve[0];  elist[5] = epcurve[1];  elist[6] = epcurve[2];  elist[7] = epcurve[3];
    status = EG_makeTopology(context, esurface[3], LOOP, CLOSED, NULL, 4, elist, sense, &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, esurface[3], FACE, SREVERSE, NULL, 1, &eloop, sense, &(efaces[7]));
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the shell and initial Body */
    status = EG_makeTopology(context, NULL, SHELL, CLOSED, NULL, 8, efaces, NULL, &eshell);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1, &eshell, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* add fillets if desired (result is ebody2) */
    if (FILLETRAD(0) > 0) {
        elist[0] = eedges[10];   elist[1] = eedges[11];   elist[2] = eedges[14];   elist[3] = eedges[15];
        status = EG_filletBody(ebody1, 4, elist, FILLETRAD(0), &ebody2, NULL);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_deleteObject(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;
    } else {
        ebody2 = ebody1;
    }

    /* add wheel if desired (result is source) */
    if (PLATETHICK(0) > 0) {
        data[0] = 0;   data[1] = 0;   data[2] =  PLATETHICK(0) / 2;
        data[3] = 0;   data[4] = 0;   data[5] = -PLATETHICK(0) / 2;
        data[6] = (MINRAD(0) + MAXRAD(0)) / 2;
        status = EG_makeSolidBody(context, CYLINDER, data, &ebody3);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_solidBoolean(ebody2, ebody3, FUSION, &emodel);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_deleteObject(ebody2);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_deleteObject(ebody3);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_getTopology(emodel, &eref, &oclass, &mtype, data, &nchild, &echilds, &senses);
        if (status != EGADS_SUCCESS) goto cleanup;

        if (oclass != MODEL || nchild != 1) {
            printf("You didn't input a model or are returning more than one body ochild = %d, nchild = %d/n", oclass, nchild);
            status = -999;
            goto cleanup;
        }

        status = EG_copyObject(echilds[0], NULL, &source);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_deleteObject(emodel);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* add bolt holes */
        for (i = 0; i < NINT(BOLTS(0)); i++) {
            data[0] = PATTERN(0) * cos(i * (2 * PI / BOLTS(0)));
            data[1] = PATTERN(0) * sin(i * (2 * PI / BOLTS(0)));
            data[2] = PLATETHICK(0) / 2;
            data[3] = PATTERN(0) * cos(i * (2 * PI / BOLTS(0)));
            data[4] = PATTERN(0) * sin(i * (2 * PI / BOLTS(0)));
            data[5] = -PLATETHICK(0) / 2;
            data[6] = BOLTRAD(0);

            status = EG_makeSolidBody(context, CYLINDER, data, &ebody4);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_solidBoolean(source, ebody4, SUBTRACTION, &emodel);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_deleteObject(source);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_deleteObject(ebody4);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_getTopology(emodel, &eref, &oclass, &mtype, data, &nchild, &echilds2, &senses);
            if (status != EGADS_SUCCESS) goto cleanup;

            if (oclass != MODEL || nchild != 1) {
                printf("You didn't input a model or are returning more than one body ochild = %d, nchild = %d/n", oclass, nchild);
                status = -999;
                goto cleanup;
            }

            status = EG_copyObject(echilds2[0], NULL, &source);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_deleteObject(emodel);
            if (status != EGADS_SUCCESS) goto cleanup;
        }
        *ebody = source;
    } else {
        *ebody = ebody2;
    }

    /* set the output value(s) */
    status = EG_getMassProperties(*ebody, data);
    if (status != EGADS_SUCCESS) {
        goto cleanup;
    }

    VOLUME(0) = data[0];

    /* remember this model (body) */
#ifdef UDP
    udps[numUdp].ebody = *ebody;
#else
    printf("myVolume = %f\n", myVolume);
#endif

#ifdef DEBUG
    printf("udpExecute -> *ebody=%lx\n", (long)(*ebody));
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

