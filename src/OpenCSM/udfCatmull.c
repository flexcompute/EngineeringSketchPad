/*
 ************************************************************************
 *                                                                      *
 * udfCatmull -- generate Catmull-Clark subdivision surfaces            *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
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

#define NUMUDPINPUTBODYS 1
#define NUMUDPARGS       4

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define NSUBDIV( IUDP)  ((int    *) (udps[IUDP].arg[0].val))[0]
#define PROGRESS(IUDP)  ((int    *) (udps[IUDP].arg[1].val))[0]
#define AREA(    IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define VOLUME(  IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]

static char  *argNames[NUMUDPARGS] = {"nsubdiv", "progress", "area",    "volume", };
static int    argTypes[NUMUDPARGS] = {ATTRINT,   ATTRINT,    -ATTRREAL, -ATTRREAL,};
static int    argIdefs[NUMUDPARGS] = {1,         0,          0,         0,        };
static double argDdefs[NUMUDPARGS] = {1.,        0.,         0.,        0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include <assert.h>

/*
 ************************************************************************
 *                                                                      *
 * structures                                                           *
 *                                                                      *
 ************************************************************************
 */

typedef struct {
    double  xyz[3];           /* coordinates of original Node */
    double  dxyz[3];          /* change in coordinates */

    int     limit;            /* limit on Node movement */
    int     nedge;            /* number of incident Edges */
    int     nface;            /* number of incident Faces */

    ego     enode;            /* ego */
} node_TT;

typedef struct {
    int     ibeg;             /* Node at beginning (bias-0) */
    int     iend;             /* Node at end       (bias-0) */

    int     ileft;            /* Face on left      (bias-0) */
    int     irite;            /* Face of rite      (bias-0) */

    int     limit;            /* limit of Edge movement */
    int     inext;            /* next subEdge or -1 */

    ego     eedge;            /* ego */
} edge_TT;

typedef struct {
    int     isouth;           /* Edge on south     (bias-0) */
    int     ieast;            /* Edge on east      (bias-0) */
    int     inorth;           /* Edge on north     (bias-0) */
    int     iwest;            /* Edge on west      (bias-0) */

    int     isw;              /* Node at southwest (bias-0) */
    int     ise;              /* Node at southeast (bias-0) */
    int     ine;              /* Node at northeast (bias-0) */
    int     inw;              /* Node at northwest (bias-0) */

    int     js;               /* sense of south Edge */
    int     je;               /* sense of east  Edge */
    int     jn;               /* sense of north Edge */
    int     jw;               /* sense of west  Edge */

    int     ic;               /* new Face Node */
    int     limit;            /* limit on Face movement */

    ego     eface;            /* ego */
} face_TT;

typedef struct {
    int     mnode;            /* maximum number of Nodes */
    int     nnode;            /* current number of Nodes */
    node_TT *nodes;           /* array          of nodes */

    int     medge;            /* maximum number of Edges */
    int     nedge;            /* current number of Edges */
    edge_TT *edges;           /* array          of Edges */

    int     mface;            /* maximum number of Faces */
    int     nface;            /* current number of Faces */
    face_TT *faces;           /* array          of Faces */

    ego     context;          /* EGADS context */
} poly_TT;

static void *realloc_temp=NULL;              /* used by RALLOC macro */

/*
 ************************************************************************
 *                                                                      *
 * prototypes for functions defined below                               *
 *                                                                      *
 ************************************************************************
 */

static int addNode(poly_TT *poly, double x, double y, double z);
static int addEdge(poly_TT *poly, int ibeg, int iend);
static int addFace(poly_TT *poly, int isouth, int ieast, int inorth, int west, int limit);
static int subdivide(poly_TT *poly);
static int makeBrep(poly_TT *poly, ego *ebody);
static int printPoly(poly_TT *poly);

#define MAX(A,B) (((A) < (B)) ? (B) : (A))


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

    int     oclass, mtype, nchild, nloop, *senses;
    int     nnode, inode, nedge, iedge, nface, iface;
    int     ibeg, iend, is, ie, in, iw;
    int     isubdiv, one=1, limit, atype, alen;
    CINT    *tempIlist;
    double  data[18];
    CDOUBLE *tempRlist;
    CCHAR   *tempClist;
    ego     context, *ebodys, eref, *echilds, *eloops, rgeom;
    ego     *enodes=NULL, *eedges=NULL, *efaces=NULL;

    poly_TT poly;

    ROUTINE(udpExecute);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("nsubdiv  = %d\n", NSUBDIV( 0));
    printf("progress = %d\n", PROGRESS(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* initialize the poly_TT structure */
    poly.mnode = 0;
    poly.nnode = 0;
    poly.nodes = NULL;

    poly.medge = 0;
    poly.nedge = 0;
    poly.edges = NULL;

    poly.mface = 0;
    poly.nface = 0;
    poly.faces = NULL;

    poly.context = NULL;

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        printf(" udpExecute: expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        printf(" udpExecute: expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    poly.context = context;

    /* check arguments */

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("nsubdiv( %d) = %d\n", numUdp, NSUBDIV( numUdp));
    printf("progress(%d) = %d\n", numUdp, PROGRESS(numUdp));
#endif

    /* get pointers to the input Body's Nodes, Edges, and Faces */
    status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, &enodes);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebodys[0], NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    if (nnode <= 0) goto cleanup;       // needed for scan-build
    if (nedge <= 0) goto cleanup;       // needed for scan-build
    if (nface <= 0) goto cleanup;       // needed for scan-build

    if (enodes == NULL) goto cleanup;   // needed for splint
    if (eedges == NULL) goto cleanup;   // needed for splint
    if (efaces == NULL) goto cleanup;   // needed for splint

    /* check that the input Body his a solid and that all Edges
       are lines */
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getTopology(eedges[iedge], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getGeometry(eref, &oclass, &mtype, &rgeom, NULL, NULL);
        CHECK_STATUS(EG_getGeometry);

        if (oclass != CURVE || mtype != LINE) {
            printf(" udpExecute: expecting solid Body.  oclass=%d, mtype=%d\n", oclass, mtype);
            status = EGADS_GEOMERR;
            goto cleanup;
        }
    }

    /* create the initial poly by extracting from the input Body */
    for (inode = 0; inode < nnode; inode++) {
        status = EG_getTopology(enodes[inode], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = addNode(&poly, data[0], data[1], data[2]);
        CHECK_STATUS(addNode);
    }

    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getTopology(eedges[iedge], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        ibeg = EG_indexBodyTopo(ebodys[0], echilds[0]) - 1;
        iend = EG_indexBodyTopo(ebodys[0], echilds[1]) - 1;

        status = addEdge(&poly, ibeg, iend);
        CHECK_STATUS(addEdge);
    }

    for (iface = 0; iface < nface; iface++) {
        status = EG_getTopology(efaces[iface], &eref, &oclass, &mtype,
                                data, &nloop, &eloops, &senses);
        CHECK_STATUS(EG_getTopology);

        if (nloop != 1) {
            printf(" udpExecute: expecting Face %d to have one Loop %d\n", iface, nloop);
            status = EGADS_GEOMERR;
            goto cleanup;
        }

        status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if (nchild != 4) {
            printf(" udpExecute: expecting Face %d to have four Edges %d\n", iface, nchild);
            status = -999;
            goto cleanup;
        }

        if        (senses[0] == SFORWARD && senses[1] == SFORWARD &&
                   senses[2] == SREVERSE && senses[3] == SREVERSE   ) {
            is = EG_indexBodyTopo(ebodys[0], echilds[0]) - 1;
            ie = EG_indexBodyTopo(ebodys[0], echilds[1]) - 1;
            in = EG_indexBodyTopo(ebodys[0], echilds[2]) - 1;
            iw = EG_indexBodyTopo(ebodys[0], echilds[3]) - 1;
        } else if (senses[0] == SREVERSE && senses[1] == SFORWARD &&
                   senses[2] == SFORWARD && senses[3] == SREVERSE   ) {
            is = EG_indexBodyTopo(ebodys[0], echilds[1]) - 1;
            ie = EG_indexBodyTopo(ebodys[0], echilds[2]) - 1;
            in = EG_indexBodyTopo(ebodys[0], echilds[3]) - 1;
            iw = EG_indexBodyTopo(ebodys[0], echilds[0]) - 1;
        } else if (senses[0] == SREVERSE && senses[1] == SREVERSE &&
                   senses[2] == SFORWARD && senses[3] == SFORWARD   ) {
            is = EG_indexBodyTopo(ebodys[0], echilds[2]) - 1;
            ie = EG_indexBodyTopo(ebodys[0], echilds[3]) - 1;
            in = EG_indexBodyTopo(ebodys[0], echilds[0]) - 1;
            iw = EG_indexBodyTopo(ebodys[0], echilds[1]) - 1;
        } else if (senses[0] == SFORWARD && senses[1] == SREVERSE &&
                   senses[2] == SREVERSE && senses[3] == SFORWARD   ) {
            is = EG_indexBodyTopo(ebodys[0], echilds[3]) - 1;
            ie = EG_indexBodyTopo(ebodys[0], echilds[0]) - 1;
            in = EG_indexBodyTopo(ebodys[0], echilds[1]) - 1;
            iw = EG_indexBodyTopo(ebodys[0], echilds[2]) - 1;
        } else {
            printf(" udpExecute: expecting Face %d to have senses %d %d %d %d\n",
                   iface, senses[0], senses[1], senses[2], senses[3]);
            status = EGADS_GEOMERR;
            goto cleanup;
        }

        limit = 0;
        status = EG_attributeRet(efaces[iface], "limitFace", &atype, &alen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status == EGADS_SUCCESS) {
            if (atype == ATTRREAL && alen == 1) {
                limit |= NINT(tempRlist[0]);
            }
        }

        status = addFace(&poly, is, ie, in, iw, limit);
        CHECK_STATUS(addFace);
    }

    if        (PROGRESS(numUdp) == 1) {
        printf("      initial   nnode=%5d, nedge=%5d, nface=%5d\n",
               poly.nnode, poly.nedge, poly.nface);
    } else if (PROGRESS(numUdp) == 2) {
        status = printPoly(&poly);
        CHECK_STATUS(printPoly);
    }

    if (poly.edges == NULL) goto cleanup;    // needed for splint

    /* make sure that we have a good polyhedron */
    for (iedge = 0; iedge < poly.nedge; iedge++) {
        if (poly.edges[iedge].ileft < 0 ||
            poly.edges[iedge].irite < 0   ) {
            printf(" udpExecute: initial poly is not watertight. iedge=%d, ileft=%d, irite=%d\n",
                   iedge, poly.edges[iedge].ileft, poly.edges[iedge].irite);
            status = EGADS_GEOMERR;
            goto cleanup;
        }
    }

    /* perform the subdivisions */
    for (isubdiv = 0; isubdiv < NSUBDIV(numUdp); isubdiv++) {
        status = subdivide(&poly);
        CHECK_STATUS(subdivide);

        if        (PROGRESS(numUdp) == 1) {
            printf("      sdiv %3d  nnode=%5d, nedge=%5d, nface=%5d\n",
                   isubdiv+1, poly.nnode, poly.nedge, poly.nface);
        } else if (PROGRESS(numUdp) == 2) {
            status = printPoly(&poly);
            CHECK_STATUS(printPoly);
        }

    }

    /* make the Brep */
    status = makeBrep(&poly, ebody);
    CHECK_STATUS(makeBrep);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* set the output value(s) */
    status = EG_getMassProperties(*ebody, data);
    CHECK_STATUS(EG_getMassProperties);

    AREA(0)   = data[1];
    VOLUME(0) = data[0];

    /* tell OpenCSM that the Faces do not have a _body attribute */
    status = EG_attributeAdd(*ebody, "__markFaces__", ATTRINT, 1,
                             &one, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);

    /* remove data from poly structure */
    if (poly.nodes != NULL) EG_free(poly.nodes);
    if (poly.edges != NULL) EG_free(poly.edges);
    if (poly.faces != NULL) EG_free(poly.faces);

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
 *   addNode - add a Node to poly                                       *
 *                                                                      *
 ************************************************************************
 */

static int
addNode(poly_TT *poly,                  /* (in)  pointer to poly_TT */
        double  x,                      /* (in)  x-coordinate */
        double  y,                      /* (in)  y-coordinate */
        double  z)                      /* (in)  z-coordinate */
{
    int    status = EGADS_SUCCESS;

    ROUTINE(addNode);
    
#ifdef DEBUG
    printf("call addNode(%f, %f, %f) -> %d\n", x, y, z, poly->nnode);
#endif

    /* make room for new Node */
    if (poly->nnode >= poly->mnode-1) {
        poly->mnode += 50;

        if (poly->nnode == 0) {
            MALLOC(poly->nodes, node_TT, poly->mnode);
        } else {
            RALLOC(poly->nodes, node_TT, poly->mnode);
        }
    }

    /* store the coordinates */
    poly->nodes[poly->nnode].xyz[0] = x;
    poly->nodes[poly->nnode].xyz[1] = y;
    poly->nodes[poly->nnode].xyz[2] = z;

    /* initialize the change in coordinates */
    poly->nodes[poly->nnode].dxyz[0] = 0;
    poly->nodes[poly->nnode].dxyz[1] = 0;
    poly->nodes[poly->nnode].dxyz[2] = 0;

    /* initialize other information */
    poly->nodes[poly->nnode].nedge  = 0;
    poly->nodes[poly->nnode].nface  = 0;
    poly->nodes[poly->nnode].limit  = 0;
    poly->nodes[poly->nnode].enode  = NULL;

    /* increment number of Nodes */
    poly->nnode++;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   addEdge - add a Edge to poly                                       *
 *                                                                      *
 ************************************************************************
 */

static int
addEdge(poly_TT *poly,                  /* (in)  pointer to poly_TT */
        int     ibeg,                   /* (in)  beg Node (bias-0) */
        int     iend)                   /* (in)  end Node (bias-0) */
{
    int    status = EGADS_SUCCESS;

    ROUTINE(addEdge);

#ifdef DEBUG
    printf("call addEdge(ibeg=%d, iend=%d) -> %d\n", ibeg, iend, poly->nedge);
#endif

    /* make room for new Edge */
    if (poly->nedge >= poly->medge-1) {
        poly->medge += 50;

        if (poly->nedge == 0) {
            MALLOC(poly->edges, edge_TT, poly->medge);
        } else {
            RALLOC(poly->edges, edge_TT, poly->medge);
        }
    }

    /* store the end points */
    poly->edges[poly->nedge].ibeg = ibeg;
    poly->edges[poly->nedge].iend = iend;

    /* initialize other information */
    poly->edges[poly->nedge].ileft = -1;
    poly->edges[poly->nedge].irite = -1;
    poly->edges[poly->nedge].inext = -1;
    poly->edges[poly->nedge].limit =  0;
    poly->edges[poly->nedge].eedge = NULL;

    /* increment the valence of the two Nodes */
    poly->nodes[ibeg].nedge ++;
    poly->nodes[iend].nedge ++;

    /* increment number of Edges */
    poly->nedge++;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   addFace - add a Face to poly                                       *
 *                                                                      *
 ************************************************************************
 */

static int
addFace(poly_TT *poly,                  /* (in)  pointer to poly_TT */
        int     isouth,                 /* (in)  south Edge (bias-0) */
        int     ieast,                  /* (in)  east  Edge (bias-0) */
        int     inorth,                 /* (in)  north Edge (bias-0) */
        int     iwest,                  /* (in)  west  Edge (bias-0) */
        int     limit)                  /* (in)  limits on point movement */
                                        /*       =1 to limit X */
                                        /*       =2 to limit Y */
                                        /*       =4 to limit Z */
{
    int    status = EGADS_SUCCESS;

    ROUTINE(addFace);

#ifdef DEBUG
    printf("call addFace(isouth=%d, ieast=%d, inorth=%d, iwest=%d) -> %d\n", isouth, ieast, inorth, iwest, poly->nface);
#endif

    /* make room for new Face */
    if (poly->nface >= poly->mface-1) {
        poly->mface += 50;

        if (poly->nface == 0) {
            MALLOC(poly->faces, face_TT, poly->mface);
        } else {
            RALLOC(poly->faces, face_TT, poly->mface);
        }
    }

    /* store Edges */
    poly->faces[poly->nface].isouth = isouth;
    poly->faces[poly->nface].ieast  = ieast;
    poly->faces[poly->nface].inorth = inorth;
    poly->faces[poly->nface].iwest  = iwest;

    /* store the corner points */
    poly->faces[poly->nface].isw   = poly->edges[isouth].ibeg;
    poly->faces[poly->nface].ise   = poly->edges[ieast ].ibeg;
    poly->faces[poly->nface].ine   = poly->edges[inorth].iend;
    poly->faces[poly->nface].inw   = poly->edges[iwest ].iend;

    /* initialize other information */
    poly->faces[poly->nface].ic    = -1;
    poly->faces[poly->nface].limit = limit;
    poly->faces[poly->nface].eface = NULL;

    /* increment valence of 4 Nodes */
    poly->nodes[poly->faces[poly->nface].isw].nface ++;
    poly->nodes[poly->faces[poly->nface].ise].nface ++;
    poly->nodes[poly->faces[poly->nface].ine].nface ++;
    poly->nodes[poly->faces[poly->nface].inw].nface ++;

    /* set up pointers from Edges to this Face */
    poly->edges[isouth].ileft = poly->nface;
    poly->edges[ieast ].ileft = poly->nface;
    poly->edges[inorth].irite = poly->nface;
    poly->edges[iwest ].irite = poly->nface;

    /* propagate limits to its Nodes and Edges */
    poly->nodes[poly->faces[poly->nface].isw].limit |= limit;
    poly->nodes[poly->faces[poly->nface].ise].limit |= limit;
    poly->nodes[poly->faces[poly->nface].ine].limit |= limit;
    poly->nodes[poly->faces[poly->nface].inw].limit |= limit;

    poly->edges[isouth].limit |= limit;
    poly->edges[ieast ].limit |= limit;
    poly->edges[inorth].limit |= limit;
    poly->edges[iwest ].limit |= limit;

    /* increment number of Faces */
    poly->nface++;

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   subdivide - perform one Catmull-Clark subdivision                  *
 *                                                                      *
 ************************************************************************
 */

static int
subdivide(poly_TT *poly)                /* (in)  pointer to poly_TT */
{
    int    status = EGADS_SUCCESS;

    int    inode, nnode, iedge, nedge, iface, nface;
    int    ibeg, iend, ileft, irite, inew, isw, ise, ine, inw;
    double xedge, yedge, zedge, xface, yface, zface;

    ROUTINE(subdivide);

#ifdef DEBUG
    printf("calling subdivide: nnode=%d, nedge=%d, nface=%d\n", poly->nnode, poly->nedge, poly->nface);
#endif

    /* remember the numer of Nodes, Eges, and Faces */
    nnode = poly->nnode;
    nedge = poly->nedge;
    nface = poly->nface;

    /* initialize changes at the (original) Nodes */
    for (inode = 0; inode < nnode; inode++) {
        poly->nodes[inode].dxyz[0] = 0;
        poly->nodes[inode].dxyz[1] = 0;
        poly->nodes[inode].dxyz[2] = 0;
    }

    /* visit each Face */
    for (iface = 0; iface < nface; iface++) {
        isw = poly->faces[iface].isw;
        ise = poly->faces[iface].ise;
        ine = poly->faces[iface].ine;
        inw = poly->faces[iface].inw;

        /* new Face points are the average of its original Nodes */
        xface = ( poly->nodes[isw].xyz[0] + poly->nodes[ise].xyz[0]
                + poly->nodes[ine].xyz[0] + poly->nodes[inw].xyz[0]) / 4;
        yface = ( poly->nodes[isw].xyz[1] + poly->nodes[ise].xyz[1]
                + poly->nodes[ine].xyz[1] + poly->nodes[inw].xyz[1]) / 4;
        zface = ( poly->nodes[isw].xyz[2] + poly->nodes[ise].xyz[2]
                + poly->nodes[ine].xyz[2] + poly->nodes[inw].xyz[2]) / 4;

        /* add a new Node */
        status = addNode(poly, xface, yface, zface);
        CHECK_STATUS(addNode);

        poly->faces[iface].ic = poly->nnode - 1;

        /* add Face point coordinates into change at Nodes */
        poly->nodes[isw].dxyz[0] += xface;
        poly->nodes[isw].dxyz[1] += yface;
        poly->nodes[isw].dxyz[2] += zface;

        poly->nodes[ise].dxyz[0] += xface;
        poly->nodes[ise].dxyz[1] += yface;
        poly->nodes[ise].dxyz[2] += zface;

        poly->nodes[ine].dxyz[0] += xface;
        poly->nodes[ine].dxyz[1] += yface;
        poly->nodes[ine].dxyz[2] += zface;

        poly->nodes[inw].dxyz[0] += xface;
        poly->nodes[inw].dxyz[1] += yface;
        poly->nodes[inw].dxyz[2] += zface;
    }

    /* visit each Edge */
    for (iedge = 0; iedge < nedge; iedge++) {
        ibeg  = poly->edges[iedge].ibeg;
        iend  = poly->edges[iedge].iend;
        ileft = poly->faces[poly->edges[iedge].ileft].ic;
        irite = poly->faces[poly->edges[iedge].irite].ic;

        /* new Edge points are the average of its original Nodes and adjacent Face points */
        if ((poly->edges[iedge].limit & 1) == 0) {
            xedge = (  poly->nodes[ibeg ].xyz[0] + poly->nodes[iend ].xyz[0]
                     + poly->nodes[ileft].xyz[0] + poly->nodes[irite].xyz[0]) / 4;
        } else {
            xedge = (  poly->nodes[ibeg ].xyz[0] + poly->nodes[iend ].xyz[0]) / 2;
        }
        if ((poly->edges[iedge].limit & 2) == 0) {
            yedge = (  poly->nodes[ibeg ].xyz[1] + poly->nodes[iend ].xyz[1]
                     + poly->nodes[ileft].xyz[1] + poly->nodes[irite].xyz[1]) / 4;
        } else {
            yedge = (  poly->nodes[ibeg ].xyz[1] + poly->nodes[iend ].xyz[1]) / 2;
        }
        if ((poly->edges[iedge].limit & 4) == 0) {
            zedge = (  poly->nodes[ibeg ].xyz[2] + poly->nodes[iend ].xyz[2]
                     + poly->nodes[ileft].xyz[2] + poly->nodes[irite].xyz[2]) / 4;
        } else {
            zedge = (  poly->nodes[ibeg ].xyz[2] + poly->nodes[iend ].xyz[2]) / 2;
        }

        /* add a new Node */
        status = addNode(poly, xedge, yedge, zedge);
        CHECK_STATUS(addNode);

        inew = poly->nnode - 1;

        /* add twice Edge midpoint into change at Nodes */
        poly->nodes[ibeg].dxyz[0] += poly->nodes[ibeg].xyz[0] + poly->nodes[iend].xyz[0];
        poly->nodes[ibeg].dxyz[1] += poly->nodes[ibeg].xyz[1] + poly->nodes[iend].xyz[1];
        poly->nodes[ibeg].dxyz[2] += poly->nodes[ibeg].xyz[2] + poly->nodes[iend].xyz[2];

        poly->nodes[iend].dxyz[0] += poly->nodes[ibeg].xyz[0] + poly->nodes[iend].xyz[0];
        poly->nodes[iend].dxyz[1] += poly->nodes[ibeg].xyz[1] + poly->nodes[iend].xyz[1];
        poly->nodes[iend].dxyz[2] += poly->nodes[ibeg].xyz[2] + poly->nodes[iend].xyz[2];

        /* add a new Edge (for second half) */
        status = addEdge(poly, inew, iend);
        CHECK_STATUS(addEdge);

        /* modify old Edge */
        poly->edges[iedge].iend  = inew;
        poly->edges[iedge].inext = poly->nedge - 1;

        /* modify valence for effected Nodes */
        poly->nodes[inew].nedge ++;
        poly->nodes[iend].nedge --;
    }

    /* visit each Face */
    for (iface = 0; iface < nface; iface++) {

        /* create 4 Edges */
        status = addEdge(poly, poly->edges[poly->faces[iface].isouth].iend, poly->faces[iface].ic);
        CHECK_STATUS(addEdge);

        status = addEdge(poly, poly->faces[iface].ic, poly->edges[poly->faces[iface].ieast ].iend);
        CHECK_STATUS(addEdge);

        status = addEdge(poly, poly->faces[iface].ic, poly->edges[poly->faces[iface].inorth].iend);
        CHECK_STATUS(addEdge);

        status = addEdge(poly, poly->edges[poly->faces[iface].iwest ].iend, poly->faces[iface].ic);
        CHECK_STATUS(addEdge);

        /* add three new Faces */
        status = addFace(poly,
                         poly->edges[poly->faces[iface].isouth].inext,
                         poly->faces[iface].ieast,
                         poly->nedge - 3,
                         poly->nedge - 4,
                         poly->faces[iface].limit);
        CHECK_STATUS(addFace);

        status = addFace(poly,
                         poly->nedge - 3,
                         poly->edges[poly->faces[iface].ieast ].inext,
                         poly->edges[poly->faces[iface].inorth].inext,
                         poly->nedge - 2,
                         poly->faces[iface].limit);
        CHECK_STATUS(addFace);

        status = addFace(poly,
                         poly->nedge - 1,
                         poly->nedge - 2,
                         poly->faces[iface].inorth,
                         poly->edges[poly->faces[iface].iwest ].inext,
                         poly->faces[iface].limit);
        CHECK_STATUS(addFace);

        /* modify old Face */
        poly->faces[iface].ieast  = poly->nedge - 4;
        poly->faces[iface].inorth = poly->nedge - 1;

        poly->faces[iface].ise = poly->faces[poly->nface-3].isw;
        poly->faces[iface].ine = poly->faces[poly->nface-2].isw;
        poly->faces[iface].inw = poly->faces[poly->nface-1].isw;

        poly->edges[poly->nedge-4].ileft = iface;
        poly->edges[poly->nedge-1].irite = iface;
    }

    /* update original Node points */
    for (inode = 0;  inode < nnode; inode++) {
        if ((poly->nodes[inode].limit & 1) == 0) {
            poly->nodes[inode].xyz[0] = ((poly->nodes[inode].nedge - 3) * poly->nodes[inode].xyz[0] + poly->nodes[inode].dxyz[0] / poly->nodes[inode].nedge) / poly->nodes[inode].nedge;
        }
        if ((poly->nodes[inode].limit & 2) == 0) {
            poly->nodes[inode].xyz[1] = ((poly->nodes[inode].nedge - 3) * poly->nodes[inode].xyz[1] + poly->nodes[inode].dxyz[1] / poly->nodes[inode].nedge) / poly->nodes[inode].nedge;
        }
        if ((poly->nodes[inode].limit & 4) == 0) {
            poly->nodes[inode].xyz[2] = ((poly->nodes[inode].nedge - 3) * poly->nodes[inode].xyz[2] + poly->nodes[inode].dxyz[2] / poly->nodes[inode].nedge) / poly->nodes[inode].nedge;
        }
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   makeBrep - make an EGADS Brep                                      *
 *                                                                      *
 ************************************************************************
 */

static int
makeBrep(poly_TT *poly,                 /* (in)  pointer to poly_TT */
         ego     *ebody)                /* (out) pointer to EGADS Body */
{
    int        status = EGADS_SUCCESS;

    int    inode, iedge, iface, *senses=NULL;
    int    oclass, mtype, nchild;
    double begXYZ[3], endXYZ[3], data[18], trange[2];
    ego    begNode, endNode, sEdge, eEdge, nEdge, wEdge;
    ego    eref, *echilds, ecurve, enodes[2], esurf, *eedges=NULL, eloop, *efaces=NULL, eshell;

    ROUTINE(makeBrep);
    
    /* make the Nodes */
    for (inode = 0; inode < poly->nnode; inode++) {
        status = EG_makeTopology(poly->context, NULL, NODE, 0,
                                 poly->nodes[inode].xyz, 0, NULL, NULL, &(poly->nodes[inode].enode));
        CHECK_STATUS(EG_makeTopology);
    }

    /* make the Edges */
    for (iedge = 0; iedge < poly->nedge; iedge++) {

        /* get two Nodes and context */
        begNode = poly->nodes[poly->edges[iedge].ibeg].enode;
        endNode = poly->nodes[poly->edges[iedge].iend].enode;

        /* get the coordinates of begNode and endNode */
        status = EG_getTopology(begNode, &eref, &oclass, &mtype,
                                begXYZ, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getTopology(endNode, &eref, &oclass, &mtype,
                                endXYZ, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        /* create a lime between begNode and endNode */
        data[0] = begXYZ[0];
        data[1] = begXYZ[1];
        data[2] = begXYZ[2];

        data[3] = endXYZ[0] - begXYZ[0];
        data[4] = endXYZ[1] - begXYZ[1];
        data[5] = endXYZ[2] - begXYZ[2];

        status = EG_makeGeometry(poly->context, CURVE, LINE, NULL, NULL, data, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        /* find the parametric coordinates at endpoints */
        status = EG_invEvaluate(ecurve, begXYZ, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurve, endXYZ, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        /* make the Edge */
        enodes[0] = begNode;
        enodes[1] = endNode;

        status = EG_makeTopology(poly->context, ecurve, EDGE, TWONODE,
                                 trange, 2, enodes, NULL, &(poly->edges[iedge].eedge));
        CHECK_STATUS(EG_makeTopology);
    }

    /* make the Faces */
    senses = NULL;

    MALLOC(eedges, ego, 8                 );
    MALLOC(efaces, ego, poly->nface       );
    MALLOC(senses, int, MAX(8,poly->nface));

    for (iface = 0; iface < poly->nface; iface++) {

        /* get four Faces and context */
        sEdge = poly->edges[poly->faces[iface].isouth].eedge;
        eEdge = poly->edges[poly->faces[iface].ieast ].eedge;
        nEdge = poly->edges[poly->faces[iface].inorth].eedge;
        wEdge = poly->edges[poly->faces[iface].iwest ].eedge;

        /* make a temporary Loop from the Edges */
        eedges[0] = sEdge;
        eedges[1] = eEdge;
        eedges[2] = nEdge;
        eedges[3] = wEdge;

        status = EG_makeLoop(4, eedges, NULL, 0, &eloop);
        CHECK_STATUS(EG_makeLoop);

        /* make a surface from the loop */
        status = EG_isoCline(eloop, 0, 0, &esurf);
        CHECK_STATUS(EG_isoCline);

        /* remove the loop */
        status = EG_deleteObject(eloop);
        CHECK_STATUS(EG_deleteObject);

        /* make a new loop (that references esurf) */
        eedges[0] = sEdge;
        senses[0] = SFORWARD;
        status = EG_otherCurve(esurf, eedges[0], 0, &(eedges[4]));
        CHECK_STATUS(EG_otherCurve);

        eedges[1] = eEdge;
        senses[1] = SFORWARD;
        status = EG_otherCurve(esurf, eedges[1], 0, &(eedges[5]));
        CHECK_STATUS(EG_otherCurve);

        eedges[2] = nEdge;
        senses[2] = SREVERSE;
        status = EG_otherCurve(esurf, eedges[2], 0, &(eedges[6]));
        CHECK_STATUS(EG_otherCurve);

        eedges[3] = wEdge;
        senses[3] = SREVERSE;
        status = EG_otherCurve(esurf, eedges[3], 0, &(eedges[7]));
        CHECK_STATUS(EG_otherCurve);

        status = EG_makeTopology(poly->context, esurf, LOOP, CLOSED,
                                 NULL, 4, eedges, senses, &eloop);
        CHECK_STATUS(EG_makeTopology);

        /* make the Face */
        status = EG_makeTopology(poly->context, esurf, FACE, SFORWARD,
                                 NULL, 1, &eloop, NULL, &(poly->faces[iface].eface));
        CHECK_STATUS(EG_makeTopology);
    }

    /* make a Shell and solid Body */
    for (iface = 0; iface < poly->nface; iface++) {
        efaces[iface] = poly->faces[iface].eface;
        senses[iface] = SFORWARD;
    }

    status = EG_makeTopology(poly->context, NULL, SHELL, CLOSED,
                             NULL, poly->nface, efaces, senses, &eshell);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(poly->context, NULL, BODY, SOLIDBODY,
                             NULL, 1, &eshell, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

cleanup:
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);
    if (senses != NULL) EG_free(senses);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   printPoly - print poly data structure                              *
 *                                                                      *
 ************************************************************************
 */

static int
printPoly(poly_TT *poly)                /* (in)  pointer to poly_TT */
{
    int    status = EGADS_SUCCESS;

    int    inode, iedge, iface;

    printf("inode     x          y          z      nedge nface limit\n");
    for (inode = 0; inode < poly->nnode; inode++) {
        printf("%5d %10.5f %10.5f %10.5f %5d %5d %5d\n", inode,
               poly->nodes[inode].xyz[0],
               poly->nodes[inode].xyz[1],
               poly->nodes[inode].xyz[2],
               poly->nodes[inode].nedge,
               poly->nodes[inode].nface,
               poly->nodes[inode].limit);
    }

    printf("iedge  ibeg  iend ileft irite inext limit\n");
    for (iedge = 0; iedge < poly->nedge; iedge++) {
        printf("%5d %5d %5d %5d %5d %5d %5d\n", iedge,
               poly->edges[iedge].ibeg,
               poly->edges[iedge].iend,
               poly->edges[iedge].ileft,
               poly->edges[iedge].irite,
               poly->edges[iedge].inext,
               poly->edges[iedge].limit);
    }

    printf("iface    is    ie    in    iw   isw   ise   ine   inw    ic limit\n");
    for (iface = 0; iface < poly->nface; iface++) {
        printf("%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n", iface,
               poly->faces[iface].isouth,
               poly->faces[iface].ieast,
               poly->faces[iface].inorth,
               poly->faces[iface].iwest,
               poly->faces[iface].isw,
               poly->faces[iface].ise,
               poly->faces[iface].ine,
               poly->faces[iface].inw,
               poly->faces[iface].ic,
               poly->faces[iface].limit);
    }

//cleanup:
    return status;
}
