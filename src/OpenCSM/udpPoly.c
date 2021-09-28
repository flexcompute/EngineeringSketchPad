/*
 ************************************************************************
 *                                                                      *
 * udpPoly -- udp file to generate general polyherda                    *
 *                                                                      *
 *            Written by Paul Mokotoff @ Syracuse University            *
 *            Modified by John Dannenhoffer @ Syracuse University       *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2019/2021  John F. Dannenhoffer, III (Syracuse University)
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

/* define necessary UDP arguments */
#define MAXUDPINPUTBODYS 0
#define NUMUDPARGS       1

/* define tolerance for Node comparison */
#define EPS06  1.0e-6

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define POINTS(IUDP,I,J)  ((double *) (udps[IUDP].arg[0].val))[(I)+3*(J)]

static char  *argNames[NUMUDPARGS] = {"points",    };
//$$$static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, };
static int    argTypes[NUMUDPARGS] = {ATTRREAL,    };
static int    argIdefs[NUMUDPARGS] = {0,           };
static double argDdefs[NUMUDPARGS] = {0.,          };

/* get utility routines */
#include "udpUtilities.c"

/* prototype functions defined below */
static int PM_makeNode(ego context, int inode, double nodeval[], ego enode[], int *nnode);
static int PM_makeEdge(ego context, int iedge, ego enodeA, ego enodeB, ego eedges[], int *nedge);
static int PM_makeFace(ego context, ego eedgeA, ego eedgeB, ego eedgeC, ego eedgeD, int senses[], ego efaces[], int *nface);

/******************************************************************************/
/*********************** Node, Edge, and Face Locations ***********************/
/******************************************************************************/

/*
      Node locations:            Edge locations:            Face locations:

              ^ J                       ^ J                        ^ J           0: imin
              |                         |                          |             1: imax
              3----------2              |----1------               |-----------  2: jmin
             /|         /|             /|         /|              /|         /|  3: jmax
            / |        / |            / 4        / |             / |    3   / |  4: kmin
           /  |       /  |          10  |       11 5            /  |       /  |  5: kmax
          7----------6   |          /-----3----/   |           /----------/ 4 |
          |   0------|---1  --> I   |   ----0--|----  --> I    | 0 -------|----  --> I
          |  /       |  /           6  /       7  /            |  /       | 1/
          | /        | /            | 8        | 9             | /    2   | /
          |/         |/             |/         |/              |/  5      |/
          4----------5              -----2------               ------------
         /                         /                          /
        K                         K                          K
*/


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
    int status = EGADS_SUCCESS;

    int    i, senses[4], nnode, nedge, nface, oclass, mtype, nchild, *senses2;
    double nodeval[24], trange[2], data[4];
    ego    enodes[8], eedges[12], eloop, efaces[6], eshell, emodel, eref, *echilds;

    ROUTINE(udpExecute);
    
#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    for (i = 0; i < udps[0].arg[0].size/3; i++) {
        printf("points[%2d]= %10.5f %10.5f %10.5f\n", i, POINTS(0,0,i), POINTS(0,1,i), POINTS(0,2,i));
    }
#endif
        
    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if (udps[0].arg[0].size < 3) {
        printf(" udpExecute: \"points\" should contain at least 3 entries\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[0].size%3 != 0) {
        printf(" udpExecute: \"points\" should contain a multiple of 3 entries\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

    /* put POINTS into full array */
    if       (udps[0].arg[0].size == 3) {   // point
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,0);   nodeval[ 4] = POINTS(0,1,0);   nodeval[ 5] = POINTS(0,2,0);
        nodeval[ 6] = POINTS(0,0,0);   nodeval[ 7] = POINTS(0,1,0);   nodeval[ 8] = POINTS(0,2,0);
        nodeval[ 9] = POINTS(0,0,0);   nodeval[10] = POINTS(0,1,0);   nodeval[11] = POINTS(0,2,0);
        nodeval[12] = POINTS(0,0,0);   nodeval[13] = POINTS(0,1,0);   nodeval[14] = POINTS(0,2,0);
        nodeval[15] = POINTS(0,0,0);   nodeval[16] = POINTS(0,1,0);   nodeval[17] = POINTS(0,2,0);
        nodeval[18] = POINTS(0,0,0);   nodeval[19] = POINTS(0,1,0);   nodeval[20] = POINTS(0,2,0);
        nodeval[21] = POINTS(0,0,0);   nodeval[22] = POINTS(0,1,0);   nodeval[23] = POINTS(0,2,0);
    } else if (udps[0].arg[0].size == 6) {   // line
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,0);   nodeval[ 4] = POINTS(0,1,0);   nodeval[ 5] = POINTS(0,2,0);
        nodeval[ 6] = POINTS(0,0,0);   nodeval[ 7] = POINTS(0,1,0);   nodeval[ 8] = POINTS(0,2,0);
        nodeval[ 9] = POINTS(0,0,0);   nodeval[10] = POINTS(0,1,0);   nodeval[11] = POINTS(0,2,0);
        nodeval[12] = POINTS(0,0,1);   nodeval[13] = POINTS(0,1,1);   nodeval[14] = POINTS(0,2,1);
        nodeval[15] = POINTS(0,0,1);   nodeval[16] = POINTS(0,1,1);   nodeval[17] = POINTS(0,2,1);
        nodeval[18] = POINTS(0,0,1);   nodeval[19] = POINTS(0,1,1);   nodeval[20] = POINTS(0,2,1);
        nodeval[21] = POINTS(0,0,1);   nodeval[22] = POINTS(0,1,1);   nodeval[23] = POINTS(0,2,1);
    } else if (udps[0].arg[0].size == 9) {   // triangle
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,1);   nodeval[ 4] = POINTS(0,1,1);   nodeval[ 5] = POINTS(0,2,1);
        nodeval[ 6] = POINTS(0,0,2);   nodeval[ 7] = POINTS(0,1,2);   nodeval[ 8] = POINTS(0,2,2);
        nodeval[ 9] = POINTS(0,0,2);   nodeval[10] = POINTS(0,1,2);   nodeval[11] = POINTS(0,2,2);
        nodeval[12] = POINTS(0,0,0);   nodeval[13] = POINTS(0,1,0);   nodeval[14] = POINTS(0,2,0);
        nodeval[15] = POINTS(0,0,1);   nodeval[16] = POINTS(0,1,1);   nodeval[17] = POINTS(0,2,1);
        nodeval[18] = POINTS(0,0,2);   nodeval[19] = POINTS(0,1,2);   nodeval[20] = POINTS(0,2,2);
        nodeval[21] = POINTS(0,0,2);   nodeval[22] = POINTS(0,1,2);   nodeval[23] = POINTS(0,2,2);
    } else if (udps[0].arg[0].size == 12) {  // quadrilateral
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,1);   nodeval[ 4] = POINTS(0,1,1);   nodeval[ 5] = POINTS(0,2,1);
        nodeval[ 6] = POINTS(0,0,2);   nodeval[ 7] = POINTS(0,1,2);   nodeval[ 8] = POINTS(0,2,2);
        nodeval[ 9] = POINTS(0,0,3);   nodeval[10] = POINTS(0,1,3);   nodeval[11] = POINTS(0,2,3);
        nodeval[12] = POINTS(0,0,0);   nodeval[13] = POINTS(0,1,0);   nodeval[14] = POINTS(0,2,0);
        nodeval[15] = POINTS(0,0,1);   nodeval[16] = POINTS(0,1,1);   nodeval[17] = POINTS(0,2,1);
        nodeval[18] = POINTS(0,0,2);   nodeval[19] = POINTS(0,1,2);   nodeval[20] = POINTS(0,2,2);
        nodeval[21] = POINTS(0,0,3);   nodeval[22] = POINTS(0,1,3);   nodeval[23] = POINTS(0,2,3);
    } else if (udps[0].arg[0].size == 15) {  // pyramid
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,1);   nodeval[ 4] = POINTS(0,1,1);   nodeval[ 5] = POINTS(0,2,1);
        nodeval[ 6] = POINTS(0,0,2);   nodeval[ 7] = POINTS(0,1,2);   nodeval[ 8] = POINTS(0,2,2);
        nodeval[ 9] = POINTS(0,0,3);   nodeval[10] = POINTS(0,1,3);   nodeval[11] = POINTS(0,2,3);
        nodeval[12] = POINTS(0,0,4);   nodeval[13] = POINTS(0,1,4);   nodeval[14] = POINTS(0,2,4);
        nodeval[15] = POINTS(0,0,4);   nodeval[16] = POINTS(0,1,4);   nodeval[17] = POINTS(0,2,4);
        nodeval[18] = POINTS(0,0,4);   nodeval[19] = POINTS(0,1,4);   nodeval[20] = POINTS(0,2,4);
        nodeval[21] = POINTS(0,0,4);   nodeval[22] = POINTS(0,1,4);   nodeval[23] = POINTS(0,2,4);
    } else if (udps[0].arg[0].size == 18) {  // wedge
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,1);   nodeval[ 4] = POINTS(0,1,1);   nodeval[ 5] = POINTS(0,2,1);
        nodeval[ 6] = POINTS(0,0,2);   nodeval[ 7] = POINTS(0,1,2);   nodeval[ 8] = POINTS(0,2,2);
        nodeval[ 9] = POINTS(0,0,2);   nodeval[10] = POINTS(0,1,2);   nodeval[11] = POINTS(0,2,2);
        nodeval[12] = POINTS(0,0,3);   nodeval[13] = POINTS(0,1,3);   nodeval[14] = POINTS(0,2,3);
        nodeval[15] = POINTS(0,0,4);   nodeval[16] = POINTS(0,1,4);   nodeval[17] = POINTS(0,2,4);
        nodeval[18] = POINTS(0,0,5);   nodeval[19] = POINTS(0,1,5);   nodeval[20] = POINTS(0,2,5);
        nodeval[21] = POINTS(0,0,5);   nodeval[22] = POINTS(0,1,5);   nodeval[23] = POINTS(0,2,5);
    } else if (udps[0].arg[0].size == 24) {  // hexahedron
        nodeval[ 0] = POINTS(0,0,0);   nodeval[ 1] = POINTS(0,1,0);   nodeval[ 2] = POINTS(0,2,0);
        nodeval[ 3] = POINTS(0,0,1);   nodeval[ 4] = POINTS(0,1,1);   nodeval[ 5] = POINTS(0,2,1);
        nodeval[ 6] = POINTS(0,0,2);   nodeval[ 7] = POINTS(0,1,2);   nodeval[ 8] = POINTS(0,2,2);
        nodeval[ 9] = POINTS(0,0,3);   nodeval[10] = POINTS(0,1,3);   nodeval[11] = POINTS(0,2,3);
        nodeval[12] = POINTS(0,0,4);   nodeval[13] = POINTS(0,1,4);   nodeval[14] = POINTS(0,2,4);
        nodeval[15] = POINTS(0,0,5);   nodeval[16] = POINTS(0,1,5);   nodeval[17] = POINTS(0,2,5);
        nodeval[18] = POINTS(0,0,6);   nodeval[19] = POINTS(0,1,6);   nodeval[20] = POINTS(0,2,6);
        nodeval[21] = POINTS(0,0,7);   nodeval[22] = POINTS(0,1,7);   nodeval[23] = POINTS(0,2,7);
    } else {
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* make the Nodes */
    nnode = 0;

    for (i = 0; i < 8; i++) {
        status = PM_makeNode(context, i, nodeval, enodes, &nnode);
        CHECK_STATUS(PM_makeNode);
    }

#ifdef DEBUG
    for (i = 0; i < 8; i++) {
        printf("enodes[%d]=%llx\n", i, (long long)enodes[i]);
//        ocsmPrintEgo(enodes[i]);
    }
#endif

    /* if a NodeBody is being made, make it now */
    if (nnode == 1) {
        trange[0] = 0;   trange[1] = 1;
        status = EG_makeTopology(context, NULL, EDGE, DEGENERATE, trange, 1, &enodes[0], NULL, &eedges[0]);
        CHECK_STATUS(EG_makeTopology);

        senses[0] = SFORWARD;
        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 1, &eedges[0], &senses[0], &eloop);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

        goto cleanup;
    }

    /* make the Curves and Edges */
    nedge = 0;

    /* imin -> imax Edges */
    status = PM_makeEdge(context,  0, enodes[0], enodes[1], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  1, enodes[3], enodes[2], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  2, enodes[4], enodes[5], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  3, enodes[7], enodes[6], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    /* jmin -> jmax Edges */
    status = PM_makeEdge(context,  4, enodes[0], enodes[3], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  5, enodes[1], enodes[2], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  6, enodes[4], enodes[7], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  7, enodes[5], enodes[6], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    /* kmin -> kmax Edges */
    status = PM_makeEdge(context,  8, enodes[0], enodes[4], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context,  9, enodes[3], enodes[7], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context, 10, enodes[1], enodes[5], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

    status = PM_makeEdge(context, 11, enodes[2], enodes[6], eedges, &nedge);
    CHECK_STATUS(PM_makeEdge);

#ifdef DEBUG
    for (i = 0; i < 12; i++) {
        printf("eedges[%2d]=%llx\n", i, (long long)eedges[i]);
//        ocsmPrintEgo(eedges[i]);
    }
#endif

    /* if a WireBody is being made, the body can now be built. */
    if (nnode == 2 && nedge == 1) {
        for (i = 0; i < 12; i++) {
            if (eedges[i] != NULL) {
                senses[0] = SFORWARD;
                status = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, 1, &eedges[i], &senses[0], &eloop);
                CHECK_STATUS(EG_makeTopology);

                status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, ebody);
                CHECK_STATUS(EG_makeTopology);

                /* remember this model (body) */
                udps[numUdp].ebody = *ebody;
                
                goto cleanup;
            }
        }
    }

    /* make the Faces */
    nface = 0;

    senses[0] = SFORWARD;
    senses[1] = SFORWARD;
    senses[2] = SREVERSE;
    senses[3] = SREVERSE;

    /* make imin Face */
    status = PM_makeFace(context, eedges[ 8], eedges[ 6], eedges[ 9], eedges[ 4], senses, efaces, &nface);
    CHECK_STATUS(PM_makeFace);

    /* make imax Face (if distinct from imin Face) */
    if (enodes[0] != enodes[1] || enodes[3] != enodes[2] || enodes[4] != enodes[5] || enodes[7] != enodes[6]) {
        status = PM_makeFace(context, eedges[ 5], eedges[11], eedges[ 7], eedges[10], senses, efaces, &nface);
        CHECK_STATUS(PM_makeFace);
    }

    /* make jmin Face */
    status = PM_makeFace(context, eedges[ 0], eedges[10], eedges[ 2], eedges[ 8], senses, efaces, &nface);
    CHECK_STATUS(PM_makeFace);

    /* make jmax Face (if distinct from jmin Face) */
    if (enodes[0] != enodes[3] || enodes[1] != enodes[2] || enodes[4] != enodes[7] || enodes[5] != enodes[6]) {
        status = PM_makeFace(context, eedges[ 9], eedges[ 3], eedges[11], eedges[ 1], senses, efaces, &nface);
        CHECK_STATUS(PM_makeFace);
    }

    /* make kmin Face */
    status = PM_makeFace(context, eedges[ 4], eedges[ 1], eedges[ 5], eedges[ 0], senses, efaces, &nface);
    CHECK_STATUS(PM_makeFace);

    /* make kmax Face (if distinct from kmin Face) */
    if (enodes[0] != enodes[4] || enodes[1] != enodes[5] || enodes[2] != enodes[6] || enodes[3] != enodes[7]) {
        status = PM_makeFace(context, eedges[ 2], eedges[ 7], eedges[ 3], eedges[ 6], senses, efaces, &nface);
        CHECK_STATUS(PM_makeFace);
    }

#ifdef DEBUG
    for (i = 0; i < 6; i++) {
        printf("efaces[%d]=%llx\n", i, (long long)efaces[i]);
    }
#endif

    /* SheetBody (triangle or quadrilateral) */
    if (nface == 1) {
        eshell = NULL;
        for (i = 0; i < 6; i++) {
            if (efaces[i] != NULL) {
                status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL, 1, &efaces[i], NULL, &eshell);
                CHECK_STATUS(EG_makeTopology);

                break;
            }
        }

        if (eshell != NULL) {
            status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1, &eshell, NULL, ebody);
            CHECK_STATUS(EG_makeTopology);
        } else {
            printf(" udpExecute: we could not make a SHELL\n");
            status = -999;
            goto cleanup;
        }

    /* SolidBody */
    } else {
        status = EG_sewFaces(nface, efaces, 0, 0, &emodel);
        CHECK_STATUS(EG_sewFaces);

        status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses2);
        CHECK_STATUS(EG_getTopology);

        status = EG_copyObject(echilds[0], NULL, ebody);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(emodel);
        CHECK_STATUS(EG_deleteObject);
    }

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
    }
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   PM_makeNode - make a unique Node                                   *
 *                                                                      *
 ************************************************************************
 */

static int
PM_makeNode(ego    context,             /* (in)  EGADS context */
            int    inode,               /* (in)  Node being generated */
            double nodeval[],           /* (in)  array of 8 corner coordinates */
            ego    enodes[],            /* (in)  enodes associated with 0<=i<inode */
                                        /* (out) enode  associated with inode */
            int    *nnode)              /* (both)number of unique Nodes */
{
    int status = EGADS_SUCCESS;

    int i;

    /* reuse a previous Node if coordinates are the same */
    for (i = 0; i < inode; i++) {
        if (fabs(nodeval[3*inode  ]-nodeval[3*i  ]) < EPS06 &&
            fabs(nodeval[3*inode+1]-nodeval[3*i+1]) < EPS06 &&
            fabs(nodeval[3*inode+2]-nodeval[3*i+2]) < EPS06   ) {
            enodes[inode] = enodes[i];
            goto cleanup;
        }
    }

    /* otherwise create a new Node */
    (*nnode)++;

    status = EG_makeTopology(context, NULL, NODE, 0, &(nodeval[3*inode]), 0, NULL, NULL, &enodes[inode]);

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   PM_makeEdge - make a unique Edge                                   *
 *                                                                      *
 ************************************************************************
 */

static int
PM_makeEdge(ego    context,             /* (in)  EGADS context */
            int    iedge,               /* (in)  Edge beging generated */
            ego    enodeA,              /* (in)  Node at beginning */
            ego    enodeB,              /* (in)  Node at end */
            ego    eedges[],            /* (in)  eedges associated with 0<=i<iedge */
                                        /* (out) eedge  associated with iedge (or NULL) */
            int    *nedge)              /* (both)number of unique Edges */
{
    int status = EGADS_SUCCESS;

    int    i, oclass, mtype, nchild, *senses;
    double nodeA[3], nodeB[3], data[18], trange[2], fill[18];
    ego    eref, ecurve, *echilds, echild[2];

    ROUTINE(PM_makeEdge);

    /* if Nodes are the same, return a NULL */
    if (enodeA == enodeB) {
        eedges[iedge] = NULL;
        goto cleanup;
    }

    /* check if edge has already been made. */
    for (i = 0; i < iedge; i++) {
        if (eedges[i] == NULL) continue;

        status = EG_getTopology(eedges[i], &eref, &oclass, &mtype, data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if        (enodeA == echilds[0] && enodeB == echilds[1]) {
            eedges[iedge] = eedges[i];
            goto cleanup;
        } else if (enodeA == echilds[1] && enodeB == echilds[0]) {
            eedges[iedge] = eedges[i];
            goto cleanup;
        }
    }

    /* getting here means that we need to make the Edge */
    (*nedge)++;

    /* make the Curve */
    status = EG_getTopology(enodeA, &eref, &oclass, &mtype, nodeA, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(enodeB, &eref, &oclass, &mtype, nodeB, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    data[0] = nodeA[0];
    data[1] = nodeA[1];
    data[2] = nodeA[2];
    data[3] = nodeB[0] - nodeA[0];
    data[4] = nodeB[1] - nodeA[1];
    data[5] = nodeB[2] - nodeA[2];

    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    /* find trange values */
    status = EG_invEvaluate(ecurve, nodeA, &trange[0], fill);
    CHECK_STATUS(EG_invEvaluate);

    status = EG_invEvaluate(ecurve, nodeB, &trange[1], fill);
    CHECK_STATUS(EG_invEvaluate);

    /* make the Edge */
    echild[0] = enodeA;
    echild[1] = enodeB;
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, echild, NULL, &(eedges[iedge]));
    CHECK_STATUS(EG_makeTopology);

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   PM_makeFace - make a Face                                          *
 *                                                                      *
 ************************************************************************
 */

static int
PM_makeFace(ego    context,             /* (in)  EGADS context */
            ego    eedgeA,              /* (in)  first  eedge */
            ego    eedgeB,              /* (in)  second eedge */
            ego    eedgeC,              /* (in)  third  eedge */
            ego    eedgeD,              /* (in)  fourth eedge */
            int    senses[],            /* (in)  array of senses (for Loop) */
            ego    efaces[],            /* (in)  efaces associated with 0<=i<*nface */
                                        /* (out) eface  associated with *nface */
            int    *nface)              /* (both)number of unique Faces */
{
    int status = EGADS_SUCCESS;

    int    nedge, periodic, header[7], ndata=0, nnode, oclass, mtype, nchild;
    double trange[4], data[20], xyzsw[3], xyzse[3], xyznw[3], xyzne[3];
    ego    myEedges[8], eloop, esurface, eref, *enodes, *echilds;

    ROUTINE(PM_makeFace);

    /* count the number of non-NULL Edges */
    nedge = 0;
    if (eedgeA != NULL) {
        myEedges[nedge] = eedgeA;
        nedge++;
    }
    if (eedgeB != NULL) {
        myEedges[nedge] = eedgeB;
        nedge++;
    }
    if (eedgeC != NULL) {
        myEedges[nedge] = eedgeC;
        nedge++;
    }
    if (eedgeD != NULL) {
        myEedges[nedge] = eedgeD;
        nedge++;
    }

    /* if there are 3 Edges, make a planar Face */
    if (nedge == 3) {

        /* make a Loop of the 3 Edges */
        status = EG_makeLoop(3, myEedges, NULL, 0, &eloop);
        CHECK_STATUS(EG_makeLoop);

        /* make a Face from the Loop */
        status = EG_makeFace(eloop, SFORWARD, NULL, &(efaces[*nface]));
        CHECK_STATUS(EG_makeFace);

        (*nface)++;

    /* if there are 4 Edges, make a warped surfaces and its Face */
    } else if (nedge == 4) {
        status = EG_getTopology(myEedges[0], &eref, &oclass, &mtype,
                                data, &nnode, &enodes, &senses);
        CHECK_STATUS(EG_getTopology);


        status = EG_getTopology(enodes[0], &eref, &oclass, &mtype,
                                xyzsw, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getTopology(enodes[1], &eref, &oclass, &mtype,
                                xyzse, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);
            
        status = EG_getTopology(myEedges[2], &eref, &oclass, &mtype,
                                data, &nnode, &enodes, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getTopology(enodes[0], &eref, &oclass, &mtype,
                                xyznw, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzne,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);
            
        header[0] = 0;              // bitflag
        header[1] = 1;              // u degree
        header[2] = 2;              // u CP
        header[3] = 4;              // u knot
        header[4] = 1;              // v degree
        header[5] = 2;              // v CP
        header[6] = 4;              // v knot

        data[ndata++] = 0;          // u knots
        data[ndata++] = 0;
        data[ndata++] = 1;
        data[ndata++] = 1;

        data[ndata++] = 0;          // v knots
        data[ndata++] = 0;
        data[ndata++] = 1;
        data[ndata++] = 1;

        data[ndata++] = xyzsw[0];   // sw CP
        data[ndata++] = xyzsw[1];
        data[ndata++] = xyzsw[2];

        data[ndata++] = xyzse[0];   // se CP
        data[ndata++] = xyzse[1];
        data[ndata++] = xyzse[2];

        data[ndata++] = xyznw[0];   // nw CP
        data[ndata++] = xyznw[1];
        data[ndata++] = xyznw[2];

        data[ndata++] = xyzne[0];   // ne CP
        data[ndata++] = xyzne[1];
        data[ndata++] = xyzne[2];

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, data, &esurface);
        CHECK_STATUS(EG_makeGeometry);

        status = EG_getRange(esurface, trange, &periodic);
        CHECK_STATUS(EG_getRange);

        /* make the Face from the Loop */
        status = EG_makeFace(esurface, SFORWARD, trange, &(efaces[*nface]));
        CHECK_STATUS(EG_makeFace);

        (*nface)++;
    }

cleanup:
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
    int iudp, judp, ipnt;

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

    /* these are dummy returns */
//$$$    if        (entType == OCSM_NODE) {
//$$$        for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$            vels[3*ipnt  ] = entIndex;
//$$$            vels[3*ipnt+1] = 0;
//$$$            vels[3*ipnt+2] = 0;
//$$$        }
//$$$    } else if (entType == OCSM_EDGE) {
//$$$        for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$            vels[3*ipnt  ] = 0;
//$$$            vels[3*ipnt+1] = entIndex;
//$$$            vels[3*ipnt+2] = 0;
//$$$        }
//$$$    } else if (entType == OCSM_FACE) {
//$$$        for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$            vels[3*ipnt  ] = 0;
//$$$            vels[3*ipnt+1] = 0;
//$$$            vels[3*ipnt+2] = entIndex;
//$$$        }
//$$$    } else {
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            vels[3*ipnt  ] = 0;
            vels[3*ipnt+1] = 0;
            vels[3*ipnt+2] = 0;
        }
//$$$    }

    return EGADS_SUCCESS;
}
