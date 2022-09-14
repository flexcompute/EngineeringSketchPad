/*
 ************************************************************************
 *                                                                      *
 * udpRadwaf -- udp file to generate radial waffle                      *
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

#define NUMUDPARGS 4
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define YSIZE( IUDP)      ((double *) (udps[IUDP].arg[0].val))[0]
#define ZSIZE( IUDP)      ((double *) (udps[IUDP].arg[1].val))[0]
#define NSPOKE(IUDP)      ((int    *) (udps[IUDP].arg[2].val))[0]
#define XFRAME(IUDP,I)    ((double *) (udps[IUDP].arg[3].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"ysize",  "zsize",  "nspoke", "xframe", };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRINT,  ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,        0,        0,        -1,       };
static double argDdefs[NUMUDPARGS] = {0.,       0.,       0.,       -1.,      };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define EPS03  0.001

static int makeNode(ego context, double xyz[], ego *enode);
static int makeEdge(ego enode1, ego enode2, ego *eedge);
static int makeFrame(ego context, int iframe, double xframe, double ysize, double zsize,
                     int nspoke, ego enodes[], ego eedges[], ego efaces[]);


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

    int     iframe, nframe, ispoke, nface, iface, senList[4], faceattr[2];
    double  area;
    ego     *enodes0=NULL, *eedges0=NULL, *efaces0=NULL;
    ego     *enodes1=NULL, *eedges1=NULL, *efaces1=NULL;
    ego     edgList[4], *facList=NULL, eloop, eshell;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("ysize( 0)   = %f\n", YSIZE( 0));
    printf("zsize( 0)   = %f\n", ZSIZE( 0));
    printf("nspoke(0)   = %d\n", NSPOKE(0));
    for (iframe = 1; iframe < udps[0].arg[3].size; iframe++) {
        printf("xframe(0,%d) = %f\n", iframe, XFRAME(0,iframe));
    }
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    nframe = udps[0].arg[3].size;

    if (udps[0].arg[0].size > 1) {
        printf(" udpExecute: ysize should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (YSIZE(0) <= 0) {
        printf(" udpExecute: ysize = %f <= 0\n", YSIZE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        printf(" udpExecute: zsize should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (ZSIZE(0) <= 0) {
        printf(" udpExecute: zsize = %f <= 0\n", ZSIZE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        printf(" udpExecute: nspoke should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (NSPOKE(0) < 2) {
        printf(" udpExecute: nspoke = %d < 2\n", NSPOKE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size < 2) {
        printf(" udpExecute: xframe should contain at least 2 values\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else {
        for (iframe = 1; iframe < nframe; iframe++) {
            if (XFRAME(0,iframe) <= XFRAME(0,iframe-1)) {
                printf(" udpExecute: xframe should be in ascending order\n");
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("ysize( %d)   = %f\n", numUdp, YSIZE( numUdp));
    printf("zsize( %d)   = %f\n", numUdp, ZSIZE( numUdp));
    printf("nspoke(%d)   = %d\n", numUdp, NSPOKE(numUdp));
    for (iframe = 0; iframe < nframe; iframe++) {
        printf("xframe(%d,%d) = %f\n", numUdp, iframe, XFRAME(numUdp,iframe));
    }
#endif

    MALLOC(enodes0, ego, (NSPOKE(numUdp)+1));
    MALLOC(eedges0, ego, (NSPOKE(numUdp)+1));
    MALLOC(efaces0, ego, (NSPOKE(numUdp)+1));
    MALLOC(enodes1, ego, (NSPOKE(numUdp)+1));
    MALLOC(eedges1, ego, (NSPOKE(numUdp)+1));
    MALLOC(efaces1, ego, (NSPOKE(numUdp)+1));

    MALLOC(facList, ego, (2*NSPOKE(numUdp)*nframe));

    eedges0[0] = NULL;
    efaces0[0] = NULL;
    eedges1[0] = NULL;
    efaces1[0] = NULL;

    /* make star-shaped pattern at x=xframe(0) */
    status = makeFrame(context, 1, XFRAME(numUdp,0), YSIZE(numUdp), ZSIZE(numUdp), NSPOKE(numUdp),
                       enodes0, eedges0, efaces0);
    CHECK_STATUS(makeFrame);

    nface = 0;
    for (iface = 1; iface <= NSPOKE(numUdp); iface++) {
        facList[nface++] = efaces0[iface];
    }

    /* make star-shaped pattern at rest of xframes */
    for (iframe = 1; iframe < nframe; iframe++) {
        status = makeFrame(context, iframe+2, XFRAME(numUdp,iframe), YSIZE(numUdp), ZSIZE(numUdp), NSPOKE(numUdp),
                           enodes1, eedges1, efaces1);
        CHECK_STATUS(makeFrame);

        for (iface = 1; iface <= NSPOKE(numUdp); iface++) {
            facList[nface++] = efaces1[iface];
        }

        /* make the Faces connecting the two star-shaped patterns */
        status = makeEdge(enodes0[0], enodes1[0], &edgList[0]);
        senList[0] = SFORWARD;
        CHECK_STATUS(makeEdge);

        for (ispoke = 1; ispoke <= NSPOKE(numUdp); ispoke++) {
            edgList[1] = eedges1[ispoke];
            senList[1] = SFORWARD;

            status = makeEdge(enodes1[ispoke], enodes0[ispoke], &edgList[2]);
            senList[2] = SFORWARD;
            CHECK_STATUS(makeEdge);

            edgList[3] = eedges0[ispoke];
            senList[3] = SREVERSE;

            status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL,
                                     4, edgList, senList, &eloop);
            CHECK_STATUS(EG_makeTopology);

            status = EG_getArea(eloop, NULL, &area);
            CHECK_STATUS(EG_getArea);

            if (area < 0) {
                status = EG_makeFace(eloop, SFORWARD, NULL, &facList[nface++]);
            } else {
                status = EG_makeFace(eloop, SREVERSE, NULL, &facList[nface++]);
            }
            CHECK_STATUS(EG_makeFace);

            faceattr[0] = ispoke;
            faceattr[1] = iframe + 1;

            status = EG_attributeAdd(facList[nface-1], "spoke", ATTRINT,
                                     2, faceattr, NULL, NULL);
            CHECK_STATUS(EG_attributeAdd);
        }

        /* set up for next pass through this loop */
        for (ispoke = 0; ispoke <= NSPOKE(numUdp); ispoke++) {
            enodes0[ispoke] = enodes1[ispoke];
            eedges0[ispoke] = eedges1[ispoke];
            efaces0[ispoke] = efaces1[ispoke];
        }
    }

    /* make the Shell and then the SheetBody */
    status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL,
                             nface, facList, NULL, &eshell);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL,
                             1, &eshell, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (enodes0 != NULL) free(enodes0);
    if (eedges0 != NULL) free(eedges0);
    if (efaces0 != NULL) free(efaces0);
    if (enodes1 != NULL) free(enodes1);
    if (eedges1 != NULL) free(eedges1);
    if (efaces1 != NULL) free(efaces1);
    if (facList != NULL) free(facList);

    if (status < 0) {
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


/*
 ************************************************************************
 *                                                                      *
 *   makeNode - make a Node                                             *
 *                                                                      *
 ************************************************************************
 */

static int
makeNode(ego    context,
         double xyz[],
         ego    *enode)
{
    int    status = EGADS_SUCCESS;

#ifdef DEBUG
    printf("makeMode(%10.5f %10.5f %10.5f)\n", xyz[0], xyz[1], xyz[2]);
#endif

    status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, enode);
    
//cleanup:
#ifdef DEBUG
    ocsmPrintEgo(*enode);
#endif

    if (status != EGADS_SUCCESS) {
        printf(" udpExecute: problem in makeNode(%f, %f, %f)\n", xyz[0], xyz[1], xyz[2]);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   makeEdge - make Edge between Nodes                                 *
 *                                                                      *
 ************************************************************************
 */

static int
makeEdge(ego     enode1,
         ego     enode2,
         ego     *eedge)
{
    int    status = EGADS_SUCCESS;

    int    oclass, mtype, nchild, *senses;
    double xyz1[3], xyz2[3], data[6], trange[2];
    ego    context, eref, *echilds, ecurve, enodes[2];

    ROUTINE(makeEdge);
    
#ifdef DEBUG
    printf("makeEdge(%lx %lx)\n", (long)enode1, (long)enode2);
#endif

    status = EG_getContext(enode1, &context);
    CHECK_STATUS(EG_getContext);

    status = EG_getTopology(enode1, &eref, &oclass, &mtype,
                            xyz1, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(enode2, &eref, &oclass, &mtype,
                            xyz2, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    data[0] = xyz1[0];
    data[1] = xyz1[1];
    data[2] = xyz1[2];
    data[3] = xyz2[0] - xyz1[0];
    data[4] = xyz2[1] - xyz1[1];
    data[5] = xyz2[2] - xyz1[2];

    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    enodes[0] = enode1;
    enodes[1] = enode2;

    trange[0] = 0;
    trange[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, enodes, NULL, eedge);
    CHECK_STATUS(EG_makeTopology);

cleanup:
#ifdef DEBUG
    ocsmPrintEgo(*eedge);
#endif

    if (status != EGADS_SUCCESS) {
        printf(" uedExecute: problem in makeEdge\n");
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   makeFrame - make a frame (consisting of Nodes, Edges, and Faces)   *
 *                                                                      *
 ************************************************************************
 */

static int
makeFrame(ego    context,
          int    iframe,
          double xframe,
          double ysize,
          double zsize,
          int    nspoke,
          ego    enodes[],
          ego    eedges[],
          ego    efaces[])
{
    int    status = EGADS_SUCCESS;

    int    ispoke, jspoke, senList[6], iedge, faceattr[2];
    double xyz[3], theta0, theta1, theta2, theta3, theta, dtheta, area;
    ego    edgList[6], ebeg, eend, eloop;

    ROUTINE(makeFrame);

    /* make the Node at the center */
    xyz[0] = xframe;
    xyz[1] = 0;
    xyz[2] = 0;
    status = makeNode(context, xyz, &enodes[0]);
    CHECK_STATUS(makeNode);

    /* make the spokes (Nodes and Edges) */
    theta0 = atan(ysize / zsize);
    theta1 = PI - theta0;
    theta2 = PI + theta0;
    theta3 = TWOPI - theta0;

    dtheta = TWOPI / nspoke;

    theta = 0;
    for (ispoke = 1; ispoke <= nspoke; ispoke++) {
        if (theta <= 0) {
            xyz[0] = xframe;
            xyz[1] = 0;
            xyz[2] = zsize/2;
        } else if (theta <= theta0) {
            xyz[0] =  xframe;
            xyz[1] =  zsize/2 * tan(theta);
            xyz[2] =  zsize/2;
        } else if (theta <= theta1) {
            xyz[0] =  xframe;
            xyz[1] =  ysize/2;
            xyz[2] =  ysize/2 / tan(theta);
        } else if (theta <= theta2) {
            xyz[0] =  xframe;
            xyz[1] = -zsize/2 * tan(theta);
            xyz[2] = -zsize/2;
        } else if (theta <= theta3) {
            xyz[0] =  xframe;
            xyz[1] = -ysize/2;
            xyz[2] = -ysize/2 / tan(theta);
        } else {
            xyz[0] =  xframe;
            xyz[1] =  zsize/2 * tan(theta);
            xyz[2] =  zsize/2;
        }
        status = makeNode(context, xyz, &enodes[ispoke]);
        CHECK_STATUS(makeNode);

        status = makeEdge(enodes[0], enodes[ispoke], &eedges[ispoke]);
        CHECK_STATUS(makeEdge);

        theta += dtheta;
    }

    /* fill in the Faces between the spokes */
    theta = 0;
    for (ispoke = 1; ispoke <= nspoke; ispoke++) {
        if (ispoke < nspoke) {
            jspoke = ispoke + 1;
        } else {
            jspoke = 1;
        }

        iedge = 0;

        edgList[iedge] = eedges[ispoke];
        senList[iedge] = SFORWARD;
        iedge++;

        ebeg = enodes[ispoke];

        if (theta < theta0-EPS03 && (theta+dtheta) > theta0+EPS03) {
            xyz[0] =  xframe;
            xyz[1] = +ysize/2;
            xyz[2] = +zsize/2;
            status = makeNode(context, xyz, &eend);
            CHECK_STATUS(makeNode);

            status = makeEdge(ebeg, eend, &edgList[iedge]);
            senList[iedge] = SFORWARD;
            CHECK_STATUS(makeEdge);
            iedge++;

            ebeg = eend;
        }

        if (theta < theta1-EPS03 && (theta+dtheta) > theta1+EPS03) {
            xyz[0] =  xframe;
            xyz[1] = +ysize/2;
            xyz[2] = -zsize/2;
            status = makeNode(context, xyz, &eend);
            CHECK_STATUS(makeNode);

            status = makeEdge(ebeg, eend, &edgList[iedge]);
            senList[iedge] = SFORWARD;
            CHECK_STATUS(makeEdge);
            iedge++;

            ebeg = eend;
        }

        if (theta < theta2-EPS03 && (theta+dtheta) > theta2+EPS03) {
            xyz[0] =  xframe;
            xyz[1] = -ysize/2;
            xyz[2] = -zsize/2;
            status = makeNode(context, xyz, &eend);
            CHECK_STATUS(makeNode);

            status = makeEdge(ebeg, eend, &edgList[iedge]);
            senList[iedge] = SFORWARD;
            CHECK_STATUS(makeEdge);
            iedge++;

            ebeg = eend;
        }

        if (theta < theta3-EPS03 && (theta+dtheta) > theta3+EPS03) {
            xyz[0] =  xframe;
            xyz[1] = -ysize/2;
            xyz[2] = +zsize/2;
            status = makeNode(context, xyz, &eend);
            CHECK_STATUS(makeNode);

            status = makeEdge(ebeg, eend, &edgList[iedge]);
            senList[iedge] = SFORWARD;
            CHECK_STATUS(makeEdge);
            iedge++;

            ebeg = eend;
        }

        status = makeEdge(ebeg, enodes[jspoke], &edgList[iedge]);
        senList[iedge] = SFORWARD;
        CHECK_STATUS(makeEdge);
        iedge++;

        edgList[iedge] = eedges[jspoke];
        senList[iedge] = SREVERSE;
        iedge++;

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL,
                                 iedge, edgList, senList, &eloop);
        CHECK_STATUS(EG_makeTopology);

        status = EG_getArea(eloop, NULL, &area);
        CHECK_STATUS(EG_getArea);

        if (area < 0) {
            status = EG_makeFace(eloop, SFORWARD, NULL, &efaces[ispoke]);
        } else {
            status = EG_makeFace(eloop, SREVERSE, NULL, &efaces[ispoke]);
        }
        CHECK_STATUS(EG_makeFace);

        faceattr[0] = iframe;
        faceattr[1] = ispoke;
        status = EG_attributeAdd(efaces[ispoke], "frame", ATTRINT,
                                 2, faceattr, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        theta +=  dtheta;
    }

cleanup:
    if (status != EGADS_SUCCESS) {
        printf(" udpExecute: problem in makeFrame(%f)\n", xframe);
    }

    return status;
}
