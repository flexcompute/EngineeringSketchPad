/*
 ************************************************************************
 *                                                                      *
 * udpCreatePoly -- create a tetgen .poly input file for two Bodys      *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2021  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 2
#define NUMUDPINPUTBODYS 2
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)    ((char   *) (udps[IUDP].arg[0].val))
#define HOLE(    IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename", "hole",   };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,        };
static double argDdefs[NUMUDPARGS] = {0.,         0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define           CINT            const int
#define           CDOUBLE         const double
#define           CCHAR           const char

/* prototype for function defined below */


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

    int     oclass, mtype, nchild, *senses, stat, npnt;
    int     i, ifirst, ilast, inode, iedge, iface, ipnt, ntri, itri;
    int     npid, npid0, npid1, npid2;
    int     nnode_inner, nedge_inner, nface_inner, npid_inner, neid_inner, made_inner=0;
    int     nnode_outer, nedge_outer, nface_outer, npid_outer, neid_outer, made_outer=0;
    int     *nodePID_inner=NULL, **edgePID_inner=NULL, **facePID_inner=NULL;
    int     *nodePID_outer=NULL, **edgePID_outer=NULL, **facePID_outer=NULL;
    CINT    *ptype, *pindx, *tris, *tric;
    double  data[4], data0[14], data1[14], bbox[6], hole[3], size, params[3];
    CDOUBLE *xyz, *uv;
    ego     context, eref, *ebodys, *echilds;
    ego     ebody_inner, etess_inner=NULL, *enodes_inner=NULL, *eedges_inner=NULL, *efaces_inner=NULL;
    ego     ebody_outer, etess_outer=NULL, *enodes_outer=NULL, *eedges_outer=NULL, *efaces_outer=NULL;
    ego     etemp, topref, next, prev, etess;

    FILE    *fpPoly=NULL;

    ROUTINE(udpExecute);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("filename(0) = %s\n", FILENAME(0));
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf("hole(0)[%d]   = %f\n", i, HOLE(0,i));
    }
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check that Model was input that contains on Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        printf(" udpExecute: expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 2) {
        printf(" udpExecute: expecting Model to contain two Bodys (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* check arguments */

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf("hole(%d)[%d]   = %f\n", numUdp, i, HOLE(0,i));
    }
#endif

    /* the outer Body is the one with the larger volume */
    status = EG_getMassProperties(ebodys[0], data0);
    CHECK_STATUS(EG_getMassProperties);

    status = EG_getMassProperties(ebodys[0], data1);
    CHECK_STATUS(EG_getMassProperties);

    if (data0[0] > data1[0]) {
        ebody_outer = ebodys[0];
        ebody_inner = ebodys[1];
    } else {
        ebody_outer = ebodys[1];
        ebody_inner = ebodys[0];
    }

    /* find a point that is inside the inner Body (needed as hole
       point for tetgen)*/
    if (udps[0].arg[1].size == 3) {
        hole[0] = HOLE(0,0);
        hole[1] = HOLE(0,1);
        hole[2] = HOLE(0,2);
    } else {
        status = EG_getBoundingBox(ebody_inner, bbox);
        CHECK_STATUS(EG_getBoundingBox);

        ifirst = 52;
        ilast  = -1;
        for (i = 0; i < 51; i++) {
            hole[0] = (bbox[0] + bbox[3]) / 2;
            hole[1] = (bbox[1] + bbox[4]) / 2;
            hole[2] =  bbox[2] + (bbox[5] - bbox[2]) * i / 50.;

            status = EG_inTopology(ebody_inner, hole);
            if (status == 0) {
                ifirst = MIN(ifirst, i);
                ilast  = MAX(ilast,  i);
            }
        }

        if (ifirst > ilast) {
            printf("unable to find hole point\n");
            status = EGADS_NOTFOUND;
            goto cleanup;
        } else {
            hole[2] = bbox[2] + (bbox[5] - bbox[2]) * (ifirst + ilast) / 100.;
        }
    }

    /* look for tessellation associated with the inner Body */
    status = EG_getContext(ebody_inner, &etemp);
    CHECK_STATUS(EG_getContext);

    status = EG_getInfo(etemp, &oclass, &mtype, &topref, &prev, &next);
    CHECK_STATUS(EG_getInfo);

    while (next != NULL) {
        etemp  = next;
        status = EG_getInfo(etemp, &oclass, &mtype, &topref, &prev, &next);
        CHECK_STATUS(EG_getInfo);

        if (oclass == TESSELLATION) {
            etess  = etemp;
            status = EG_statusTessBody(etess, &etemp, &stat, &npnt);
            CHECK_STATUS(EG_statusTessBody);

            /* if we found the correct tessellation, use it */
            if (etemp == ebody_inner) {
                etess_inner = etess;
#ifdef DEBUG
                printf("skipping tessellation for inner Body\n");
#endif
                break;
            }
        }
    }

    /* tessellaate the inner Body (if it does not already exist) */
    if (etess_inner == NULL) {
#ifdef DEBUG
        printf("tessellating inner Body\n");
#endif
        status = EG_getBoundingBox(ebody_inner, bbox);
        CHECK_STATUS(EG_getBoundingBox);

        size = sqrt( (bbox[3] - bbox[0]) * (bbox[3] - bbox[0])
                    +(bbox[4] - bbox[1]) * (bbox[4] - bbox[1])
                    +(bbox[5] - bbox[2]) * (bbox[5] - bbox[2]));

        params[0] = 0.0250 * size;
        params[1] = 0.0010 * size;
        params[2] = 15.0;

        status = EG_makeTessBody(ebody_inner, params, &etess_inner);
        CHECK_STATUS(EG_makeTessBody);

        made_inner = 1;

        status = EG_attributeAdd(ebody_inner, "_tParams",
                                 ATTRREAL, 3, NULL, params, NULL);
        CHECK_STATUS(EG_attributeAdd);
    }

    /* look for tessellation associated with the outer Body */
    status = EG_getContext(ebody_outer, &etemp);
    CHECK_STATUS(EG_getContext);

    status = EG_getInfo(etemp, &oclass, &mtype, &topref, &prev, &next);
    CHECK_STATUS(EG_getInfo);

    while (next != NULL) {
        etemp  = next;
        status = EG_getInfo(etemp, &oclass, &mtype, &topref, &prev, &next);
        CHECK_STATUS(EG_getInfo);

        if (oclass == TESSELLATION) {
            etess  = etemp;
            status = EG_statusTessBody(etess, &etemp, &stat, &npnt);
            CHECK_STATUS(EG_statusTessBody);

            /* if we found the correct tessellation, use it */
            if (etemp == ebody_outer) {
                etess_outer = etess;
#ifdef DEBUG
                printf("skipping tessellation for outer Body\n");
#endif
                break;
            }
        }
    }

    /* tessellate the outer Body */
    if (etess_outer == NULL) {
#ifdef DEBUG
        printf("tessellating outer Body\n");
#endif
        status = EG_getBoundingBox(ebody_outer, bbox);
        CHECK_STATUS(EG_getBoundingBox);

        size = sqrt( (bbox[3] - bbox[0]) * (bbox[3] - bbox[0])
                    +(bbox[4] - bbox[1]) * (bbox[4] - bbox[1])
                    +(bbox[5] - bbox[2]) * (bbox[5] - bbox[2]));

        params[0] = 0.0250 * size;
        params[1] = 0.0010 * size;
        params[2] = 15.0;

        status = EG_makeTessBody(ebody_outer, params, &etess_outer);
        CHECK_STATUS(EG_makeTessBody);

        made_outer = 1;

        status = EG_attributeAdd(ebody_outer, "_tParams",
                                 ATTRREAL, 3, NULL, params, NULL);
        CHECK_STATUS(EG_attributeAdd);
    }

    /* get the statistics for the inner Body */
    status = EG_getBodyTopos(ebody_inner, NULL, NODE, &nnode_inner, &enodes_inner);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody_inner, NULL, EDGE, &nedge_inner, &eedges_inner);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody_inner, NULL, FACE, &nface_inner, &efaces_inner);
    CHECK_STATUS(EG_getBodyTopos);

    if ( etess_inner == NULL) goto cleanup;  // needed for splint
    if (enodes_inner == NULL) goto cleanup;  // needed for splint
    if (eedges_inner == NULL) goto cleanup;  // needed for splint
    if (efaces_inner == NULL) goto cleanup;  // needed for splint

    npid_inner = nnode_inner;
    neid_inner = 0;

    for (iedge = 1; iedge <= nedge_inner; iedge++) {
        status = EG_getTessEdge(etess_inner, iedge,
                                &npnt, &xyz, &uv);
        CHECK_STATUS(EG_getTessEdge);

        npid_inner += (npnt - 2);
    }

    for (iface = 1; iface <= nface_inner; iface++) {
        status = EG_getTessFace(etess_inner, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (ptype[ipnt] == -1) npid_inner++;
        }

        neid_inner += ntri;
    }
#ifdef DEBUG
    printf("inner Body: nnode=%5d, nedge=%5d, nface=%5d, npid=%8d, neid=%8d\n",
           nnode_inner, nedge_inner, nface_inner, npid_inner, neid_inner);
#endif

    /* get the statistics for the outer Body */
    status = EG_getBodyTopos(ebody_outer, NULL, NODE, &nnode_outer, &enodes_outer);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody_outer, NULL, EDGE, &nedge_outer, &eedges_outer);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody_outer, NULL, FACE, &nface_outer, &efaces_outer);
    CHECK_STATUS(EG_getBodyTopos);

    if ( etess_outer == NULL) goto cleanup;  // needed for splint
    if (enodes_outer == NULL) goto cleanup;  // needed for splint
    if (eedges_outer == NULL) goto cleanup;  // needed for splint
    if (efaces_outer == NULL) goto cleanup;  // needed for splint

    npid_outer = nnode_outer;
    neid_outer = 0;

    for (iedge = 1; iedge <= nedge_outer; iedge++) {
        status = EG_getTessEdge(etess_outer, iedge,
                                &npnt, &xyz, &uv);
        CHECK_STATUS(EG_getTessEdge);

        npid_outer += (npnt - 2);
    }

    for (iface = 1; iface <= nface_outer; iface++) {
        status = EG_getTessFace(etess_outer, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (ptype[ipnt] == -1) npid_outer++;
        }

        neid_outer += ntri;
    }
#ifdef DEBUG
    printf("outer Body: nnode=%5d, nedge=%5d, nface=%5d, npid=%8d, neid=%8d\n",
           nnode_outer, nedge_outer, nface_outer, npid_outer, neid_outer);
#endif

    /* open the output file */
    fpPoly = fopen(FILENAME(0), "w");
    if (fpPoly == NULL) {
        status = EGADS_WRITERR;
        goto cleanup;
    }

    npid = 0;

    /* part 1: write the points to the .poly file */
    fprintf(fpPoly, "# Part 1 - node list\n");
    fprintf(fpPoly, "# Node count,  3 dim,  no attributes, no boundary markers\n");
    fprintf(fpPoly, "%10d  3  0  0\n", npid_inner+npid_outer);
    fprintf(fpPoly, "# Node index, node coordinates\n");

    /* write out the points associated with the inner Body */
    MALLOC(nodePID_inner, int,  (nnode_inner+1));
    MALLOC(edgePID_inner, int*, (nedge_inner+1));
    MALLOC(facePID_inner, int*, (nface_inner+1));

    for (iedge = 0; iedge <= nedge_inner; iedge++) {
        edgePID_inner[iedge] = NULL;
    }

    for (iface = 0; iface <= nface_inner; iface++) {
        facePID_inner[iface] = NULL;
    }

    for (inode = 1; inode <= nnode_inner; inode++) {
        fprintf(fpPoly, "# inner body, node %d\n", inode);

        npid++;
        nodePID_inner[inode] = npid;

        status = EG_getTopology(enodes_inner[inode-1], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
                npid, data[0], data[1], data[2]);
    }

    for (iedge = 1; iedge <= nedge_inner; iedge++) {
        fprintf(fpPoly, "# inner body, edge %d\n", iedge);

        status = EG_getTessEdge(etess_inner, iedge,
                                &npnt, &xyz, &uv);
        CHECK_STATUS(EG_getTessEdge);

        MALLOC(edgePID_inner[iedge], int, npnt);

        for (ipnt = 1; ipnt < npnt-1; ipnt++) {
            npid++;
            edgePID_inner[iedge][ipnt] = npid;

            fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
                    npid, xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2]);
        }
    }

    for (iface = 1; iface <= nface_inner; iface++) {
        fprintf(fpPoly, "# inner body, face %d\n", iface);

        status = EG_getTessFace(etess_inner, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        MALLOC(facePID_inner[iface], int, npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (ptype[ipnt] == 0) {
                facePID_inner[iface][ipnt] = nodePID_inner[pindx[ipnt]];
            } else if (ptype[ipnt] > 0) {
                facePID_inner[iface][ipnt] = edgePID_inner[pindx[ipnt]][ptype[ipnt]-1];
            } else {
                npid++;
                facePID_inner[iface][ipnt] = npid;

                fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
                        npid, xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2]);
            }
        }
    }

    /* write out the points associated with the outer Body */
    MALLOC(nodePID_outer, int,  (nnode_outer+1));
    MALLOC(edgePID_outer, int*, (nedge_outer+1));
    MALLOC(facePID_outer, int*, (nface_outer+1));

    for (iedge = 0; iedge <= nedge_outer; iedge++) {
        edgePID_outer[iedge] = NULL;
    }

    for (iface = 0; iface <= nface_outer; iface++) {
        facePID_outer[iface] = NULL;
    }

    for (inode = 1; inode <= nnode_outer; inode++) {
        fprintf(fpPoly, "# outer body, node %d\n", inode);

        npid++;
        nodePID_outer[inode] = npid;

        status = EG_getTopology(enodes_outer[inode-1], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
                npid, data[0], data[1], data[2]);
    }

    for (iedge = 1; iedge <= nedge_outer; iedge++) {
        fprintf(fpPoly, "# outer body, edge %d\n", iedge);

        status = EG_getTessEdge(etess_outer, iedge,
                                &npnt, &xyz, &uv);
        CHECK_STATUS(EG_getTessEdge);

        MALLOC(edgePID_outer[iedge], int, npnt);

        for (ipnt = 1; ipnt < npnt-1; ipnt++) {
            npid++;
            edgePID_outer[iedge][ipnt] = npid;

            fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
                    npid, xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2]);
        }
    }

    for (iface = 1; iface <= nface_outer; iface++) {
        fprintf(fpPoly, "# outer body, face %d\n", iface);

        status = EG_getTessFace(etess_outer, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        MALLOC(facePID_outer[iface], int, npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (ptype[ipnt] == 0) {
                facePID_outer[iface][ipnt] = nodePID_outer[pindx[ipnt]];
            } else if (ptype[ipnt] > 0) {
                facePID_outer[iface][ipnt] = edgePID_outer[pindx[ipnt]][ptype[ipnt]-1];
            } else {
                npid++;
                facePID_outer[iface][ipnt] = npid;

                fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
                        npid, xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2]);
            }
        }
    }

    /* part 2: write the facets to the .poly file */
    fprintf(fpPoly, "# Part 2 - facet list\n");
    fprintf(fpPoly, "# Facet count,  1 boundary marker\n");
    fprintf(fpPoly, "%8d  1\n", neid_inner+neid_outer);
    fprintf(fpPoly, "# facets\n");

    /* write out the triangles associated with the inner Body */
    for (iface = 1; iface <= nface_inner; iface++) {
        fprintf(fpPoly, "# inner body, face %d\n", iface);

        status = EG_getTessFace(etess_inner, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        for (itri = 0; itri < ntri; itri++) {
            npid0 = facePID_inner[iface][tris[3*itri  ]-1];
            npid1 = facePID_inner[iface][tris[3*itri+1]-1];
            npid2 = facePID_inner[iface][tris[3*itri+2]-1];

            fprintf(fpPoly, "1  0  1\n");
            fprintf(fpPoly, "%8d  %8d  %8d  %8d\n",
                    3, npid0, npid1, npid2);
        }
    }

    /* write out the triangles associated with the outer Body */
    for (iface = 1; iface <= nface_outer; iface++) {
        fprintf(fpPoly, "# outer body, face %d\n", iface);

        status = EG_getTessFace(etess_outer, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        for (itri = 0; itri < ntri; itri++) {
            npid0 = facePID_outer[iface][tris[3*itri  ]-1];
            npid1 = facePID_outer[iface][tris[3*itri+1]-1];
            npid2 = facePID_outer[iface][tris[3*itri+2]-1];

            fprintf(fpPoly, "1  0  2\n");
            fprintf(fpPoly, "%8d  %8d  %8d  %8d\n",
                    3, npid0, npid1, npid2);
        }
    }

    /* part 3: write the hole (in the middle of the inner Body) to the .poly file */
    fprintf(fpPoly, "# Part 3 - hole list\n");
    fprintf(fpPoly, "# Hole count\n");
    fprintf(fpPoly, "1\n");
    fprintf(fpPoly, "%8d  %20.10f %20.10f %20.10f\n",
            1, hole[0], hole[1], hole[2]);

    /* part 4: write the region info to the .poly file */
    fprintf(fpPoly, "# Part 4 - region list\n");
    fprintf(fpPoly, "# Region count\n");
    fprintf(fpPoly, "0\n");

    /* finialize and close the file */
    fclose(fpPoly);
    fpPoly = NULL;

    /* return a copy of the inner Body */
    status = EG_copyObject(ebody_inner, NULL, ebody);
    CHECK_STATUS(EG_copyObject);

cleanup:
    if (fpPoly != NULL) {
        fprintf(fpPoly, "$$$ error encountered during writePoly\n");
        fclose(fpPoly);
    }

    /* inner Body */
    if (facePID_inner != NULL) {
        for (iface = 1; iface <= nface_inner; iface++) {
            if (facePID_inner[iface] != NULL) free(facePID_inner[iface]);
        }
        free(facePID_inner);
    }

    if (edgePID_inner != NULL) {
        for (iedge = 1; iedge <= nedge_inner; iedge++) {
            if (edgePID_inner[iedge] != NULL) free(edgePID_inner[iedge]);
        }
        free(edgePID_inner);
    }

    if (nodePID_inner != NULL) free(nodePID_inner);

    EG_free(enodes_inner);
    EG_free(eedges_inner);
    EG_free(efaces_inner);

    if (made_inner == 1 && etess_inner != NULL) EG_deleteObject(etess_inner);

    /* outer Body */
    if (facePID_outer != NULL) {
        for (iface = 1; iface <= nface_outer; iface++) {
            if (facePID_outer[iface] != NULL) free(facePID_outer[iface]);
        }
        free(facePID_outer);
    }

    if (edgePID_outer != NULL) {
        for (iedge = 1; iedge <= nedge_outer; iedge++) {
            if (edgePID_outer[iedge] != NULL) free(edgePID_outer[iedge]);
        }
        free(edgePID_outer);
    }

    if (nodePID_outer != NULL) free(nodePID_outer);

    EG_free(enodes_outer);
    EG_free(eedges_outer);
    EG_free(efaces_outer);

    if (made_outer == 1 && etess_outer != NULL) EG_deleteObject(etess_outer);

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
