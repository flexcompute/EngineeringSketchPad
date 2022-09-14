/*
 ************************************************************************
 *                                                                      *
 * udpCreateBEM -- create a BEM and write a file for a Body             *
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

#define NUMUDPARGS       5
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define SPACE(   IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define IMIN(    IUDP)  ((int    *) (udps[IUDP].arg[2].val))[0]
#define IMAX(    IUDP)  ((int    *) (udps[IUDP].arg[3].val))[0]
#define NOCROD(  IUDP)  ((int    *) (udps[IUDP].arg[4].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename",  "space",  "imin",  "imax",  "nocrod",};
static int    argTypes[NUMUDPARGS] = {ATTRSTRING,  ATTRREAL, ATTRINT, ATTRINT, ATTRINT, };
static int    argIdefs[NUMUDPARGS] = {0,           0,        3,       5,       0,       };
static double argDdefs[NUMUDPARGS] = {0.,          0.,       0.,      0.,      0.,      };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

/* prototype for function defined below */
static int createBemFile(ego ebody, char filename[], double space, int imin, int imax);


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

    int     oclass, mtype, nchild, *senses;
    double  data[4];
    ego     context, eref, *ebodys;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */
    
#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("filename(0) = %s\n", FILENAME(0));
    printf("imin(    0) = %d\n", IMIN(    0));
    printf("imax(    0) = %d\n", IMAX(    0));
    printf("nocrod(  0) = %d\n", NOCROD(  0));
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
    } else if (nchild != 1) {
        printf(" udpExecute: expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* check arguments */
    if        (udps[0].arg[1].size >= 2) {
        printf(" udpExecute: space should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (SPACE(0) <= 0) {
        printf(" udpExecute: space = %f <= 0\n", SPACE(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size >= 2) {
        printf(" udpExecute: imin should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (IMIN(0) <= 0) {
        printf(" udpExecute: imin = %d <= 0\n", IMIN(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size >= 2) {
        printf(" udpExecute: imax should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (IMAX(0) < IMIN(0)) {
        printf(" udpExecute: imax = %d < imin = %d\n", IMAX(0), IMIN(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
    printf("space(   %d) = %f\n", numUdp, SPACE(   numUdp));
    printf("imin(    %d) = %d\n", numUdp, IMIN(    numUdp));
    printf("imax(    %d) = %d\n", numUdp, IMAX(    numUdp));
    printf("nocrod(  %d) = %d\n", numUdp, NOCROD(  numUdp));
#endif

    /* make a copy of the Body (so that it does not get removed
     when OpenCSM deletes emodel) */
    status = EG_copyObject(ebodys[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* annotate the Body and create the BEM file */
    status = createBemFile(*ebody, FILENAME(numUdp), SPACE(numUdp),
                          IMIN(numUdp), IMAX(numUdp));
    CHECK_STATUS(createBemFile);

    /* the copy of the Body that was annotated is returned */
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
 *   createBemFile - create a BEM and write to file                     *
 *                                                                      *
 ************************************************************************
 */

static int
createBemFile(ego    ebody,             /* (in)  EGADS Body */
              char   filename[],        /* (in)  file to be written */
              double space,             /* (in)  nominal spacing along Edges */
              int    imin,              /* (in)  minimum points along any Edge */
              int    imax)              /* (in)  maximum points along any Edge */
{
    int       status = 0;               /* (out) return status */

    int       periodic, nchange, oclass, mtype, nchild, *senses;
    int       ncid, ngid, npid, neid, ipid, atype, alen, nattr, iattr;
    int       inode, iedge, iface, npnt, ntri, ipnt, itri, i, one=1;
    int       *npnts=NULL, *isouth=NULL, *ieast=NULL, *inorth=NULL, *iwest=NULL;
    int       *ignoreNode=NULL, *ignoreEdge=NULL, *ignoreFace=NULL;
    int       *nodeGID=NULL, **edgeGID=NULL, **faceGID=NULL;
    int       itemp, ntemp;
    const int *ptype, *pindx, *tris, *tric, *tempIlist;
    double    params[3], range[2], arclen, data[4];
    double    *rpos=NULL;
    const double *xyz, *uv, *tempRlist;
    const char *tempClist, *aname;
    FILE      *fpBEM=NULL;
    int       nnode, nedge, nface, ibeg, iend;
    ego       *enodes=NULL, *eedges=NULL, *efaces=NULL, eref, *echilds, *etemp, eloop;
    ego       etess, newTess;

    ROUTINE(createBemFile);
    
    /* --------------------------------------------------------------- */

    /* get number of Nodes, Edges, and Faces in ebody */
    status = EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    if (enodes == NULL) goto cleanup;   // needed for splint
    if (eedges == NULL) goto cleanup;   // needed for splint
    if (efaces == NULL) goto cleanup;   // needed for splint

    /* determine the nominal number of points along each Edge */
    MALLOC(npnts, int,    (nedge+1));
    MALLOC(rpos,  double, (imax   ));

    /* add .tParams to the Body */
    params[0] = 2 * space;
    params[1] =     space;
    params[2] = 30.0;

    status = EG_attributeAdd(ebody, ".tParams", ATTRREAL, 3,
                             NULL, params, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* determine the number of points (including endpoints) for each Edge */
    for (iedge = 1; iedge <= nedge; iedge++) {
        status = EG_getRange(eedges[iedge-1], range, &periodic);
        CHECK_STATUS(EG_getRange);

        status = EG_arcLength(eedges[iedge-1], range[0], range[1], &arclen);
        CHECK_STATUS(EG_arcLength);

        npnts[iedge] = MIN(MAX(MAX(imin,2), (int)(1+arclen/space)), imax);

        /* half as many points (with possible round-up) because of quadder */
        npnts[iedge] = 1 + npnts[iedge] / 2;
    }

    /* make arrays to remember if Node, Edge, or Face is to be ignored */
    MALLOC(ignoreNode, int, nnode);
    MALLOC(ignoreEdge, int, nedge);
    MALLOC(ignoreFace, int, nface);

    /* keep track of whether or not each Node, Edge, and Face is
       ignored or not */
    for (inode = 0; inode < nnode; inode++) {
        ignoreNode[inode] = 0;

        status = EG_attributeRet(enodes[inode], "ignoreNode",
                                 &atype, &alen, &tempIlist, &tempRlist, &tempClist);
        if (status == EGADS_SUCCESS && atype == ATTRSTRING) {
            if (strcmp(tempClist, "true") == 0) {
                ignoreNode[inode] = 1;
            }
        }
    }

    for (iedge = 0; iedge < nedge; iedge++) {
        ignoreEdge[iedge] = 0;

        status = EG_attributeRet(eedges[iedge], "ignoreEdge",
                                 &atype, &alen, &tempIlist, &tempRlist, &tempClist);
        if (status == EGADS_SUCCESS && atype == ATTRSTRING) {
            if (strcmp(tempClist, "true") == 0) {
                ignoreEdge[iedge] = 1;
            }
        }

        /* if the Edge is not ignored, we cannot ignore its Nodes either */
        if (ignoreEdge[iedge] == 0) {
            status = EG_getBodyTopos(ebody, eedges[iedge], NODE, &ntemp, &etemp);
            CHECK_STATUS(EG_getBodyTopos);

            for (itemp = 0; itemp < ntemp; itemp++) {
                inode = EG_indexBodyTopo(ebody, etemp[itemp]);
                if (inode > 0) {
                    ignoreNode[inode-1] = 0;
                }
            }

            EG_free(etemp);
        }
    }

    for (iface = 0; iface < nface; iface++) {
        ignoreFace[iface] = 0;

        status = EG_attributeRet(efaces[iface], "ignoreFace",
                                 &atype, &alen, &tempIlist, &tempRlist, &tempClist);
        if (status == EGADS_SUCCESS && atype == ATTRSTRING) {
            if (strcmp(tempClist, "true") == 0) {
                ignoreFace[iface] = 1;
            }
        }

        /* if the Face is not ignored, we cannot ignore its Nodes or Edges either */
        if (ignoreFace[iface] == 0) {
            status = EG_getBodyTopos(ebody, efaces[iface], NODE, &ntemp, &etemp);
            CHECK_STATUS(EG_getBodyTopos);

            for (itemp = 0; itemp < ntemp; itemp++) {
                inode = EG_indexBodyTopo(ebody, etemp[itemp]);
                if (inode > 0) {
                    ignoreNode[inode-1] = 0;
                }
            }

            EG_free(etemp);

            status = EG_getBodyTopos(ebody, efaces[iface], EDGE, &ntemp, &etemp);
            CHECK_STATUS(EG_getBodyTopos);

            for (itemp = 0; itemp < ntemp; itemp++) {
                iedge = EG_indexBodyTopo(ebody, etemp[itemp]);
                if (iedge > 0) {
                    ignoreEdge[iedge-1] = 0;
                }
            }

            EG_free(etemp);
        }
    }

    /* make arrays for "opposite" sides of four-sided Faces (with only one loop) */
    MALLOC(isouth, int, (nface+1));
    MALLOC(ieast,  int, (nface+1));
    MALLOC(inorth, int, (nface+1));
    MALLOC(iwest,  int, (nface+1));

    for (iface = 1; iface <= nface; iface++) {
        isouth[iface] = 0;
        ieast[ iface] = 0;
        inorth[iface] = 0;
        iwest[ iface] = 0;

        status = EG_getTopology(efaces[iface-1],
                                &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if (nchild != 1) continue;

        eloop = echilds[0];
        status = EG_getTopology(eloop, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if (nchild != 4) continue;

        isouth[iface] = status = EG_indexBodyTopo(ebody, echilds[0]);
        CHECK_STATUS(EG_indexBodyTopo);

        ieast[iface]  = status = EG_indexBodyTopo(ebody, echilds[1]);
        CHECK_STATUS(EG_indexBodyTopo);

        inorth[iface] = status = EG_indexBodyTopo(ebody, echilds[2]);
        CHECK_STATUS(EG_indexBodyTopo);

        iwest[iface]  = status = EG_indexBodyTopo(ebody, echilds[3]);
        CHECK_STATUS(EG_indexBodyTopo);
    }

    /* make "opposite" sides of four-sided Faces (with only one loop) match */
    nchange = 1;
    for (i = 0; i < 20; i++) {
        nchange = 0;

        for (iface = 1; iface <= nface; iface++) {
            if (ignoreFace[iface-1] == 1) continue;

            if (isouth[iface] <= 0 || ieast[iface] <= 0 ||
                inorth[iface] <= 0 || iwest[iface] <= 0   ) continue;

            if        (npnts[iwest[iface]] < npnts[ieast[iface]]) {
                npnts[iwest[iface]] = npnts[ieast[iface]];
                nchange++;
            } else if (npnts[ieast[iface]] < npnts[iwest[iface]]) {
                npnts[ieast[iface]] = npnts[iwest[iface]];
                nchange++;
            }

            if        (npnts[isouth[iface]] < npnts[inorth[iface]]) {
                npnts[isouth[iface]] = npnts[inorth[iface]];
                nchange++;
            } else if (npnts[inorth[iface]] < npnts[isouth[iface]]) {
                npnts[inorth[iface]] = npnts[isouth[iface]];
                nchange++;
            }
        }
        if (nchange == 0) break;
    }
    if (nchange > 0) {
        status = -999;
        goto cleanup;
    }

    /* mark the Edges with npnts[iedge] evenly-spaced points */
    for (iedge = 1; iedge <= nedge; iedge++) {
        for (i = 1; i < npnts[iedge]-1; i++) {
            rpos[i-1] = (double)(i) / (double)(npnts[iedge]-1);
        }

        if (npnts[iedge] == 2) {
            i = 0;
            status = EG_attributeAdd(eedges[iedge-1],
                                     ".rPos", ATTRINT, 1, &i, NULL, NULL);
            CHECK_STATUS(EG_attributeAdd);
        } else {
            status = EG_attributeAdd(eedges[iedge-1],
                                     ".rPos", ATTRREAL, npnts[iedge]-2, NULL, rpos, NULL);
            CHECK_STATUS(EG_attributeAdd);
        }
    }

    /* make the new tessellation */
    status = EG_makeTessBody(ebody, params, &etess);
    CHECK_STATUS(EG_makeTessBody);

    status = EG_attributeAdd(ebody, "_tParams",
                             ATTRREAL, 3, NULL, params, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* convert the triangles to quads */
    status = EG_quadTess(etess, &newTess);
    if (status == EGADS_SUCCESS) {
        EG_deleteObject(etess);

        etess = newTess;
    }

    /* put attribute on Body so that OpenCSM makes quads too */
    status = EG_attributeAdd(ebody, "_makeQuads", ATTRINT, 1, &one, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* initialize the CIDs, GIDs, PIDs, and EIDs for the output file */
    ncid = 0;   // coordinate  IDs
    ngid = 0;   // global grid IDs
    npid = 0;   // property    IDs
    neid = 0;   // element     IDs

    /* open the output file */
    fpBEM = fopen(filename, "w");
    if (fpBEM == NULL) {
        status = EGADS_WRITERR;
        goto cleanup;
    }

    MALLOC(nodeGID, int,  (nnode+1));
    MALLOC(edgeGID, int*, (nedge+1));
    MALLOC(faceGID, int*, (nface+1));

    /* create arrays to hold global IDs for all points in tessellation */
    for (iedge = 0; iedge <= nedge; iedge++) {
        edgeGID[iedge] = NULL;
    }

    for (iface = 0; iface <= nface; iface++) {
        faceGID[iface] = NULL;
    }

    /* write out the GRID cards and keep track of the associated GIDs */
    for (inode = 1; inode <= nnode; inode++) {
        if (ignoreNode[inode-1] == 1) continue;

        fprintf(fpBEM, "$ node %d", inode);
        status = EG_attributeNum(enodes[inode-1], &nattr);
        CHECK_STATUS(EG_attributeNum);

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(enodes[inode-1], iattr,
                                     &aname, &atype, &alen,
                                     &tempIlist, &tempRlist, &tempClist);
            CHECK_STATUS(EG_attributeGet);

            if        (atype == ATTRINT) {
                fprintf(fpBEM, " %s=", aname);
                for (i = 0; i < alen; i++) {
                    fprintf(fpBEM, "%d;", tempIlist[i]);
                }
            } else if (atype == ATTRREAL) {
                fprintf(fpBEM, " %s=", aname);
                for (i = 0; i < alen; i++) {
                    fprintf(fpBEM, "%f;", tempRlist[i]);
                }
            } else if (atype == ATTRSTRING) {
                fprintf(fpBEM, " %s=%s", aname, tempClist);
            }
        }
        fprintf(fpBEM, "\n");

        ngid++;
        nodeGID[inode] = ngid;

        status = EG_getTopology(enodes[inode-1], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        fprintf(fpBEM, "GRID    %8d %7d %7.4f %7.4f %7.4f\n",
                ngid, ncid, data[0],
                            data[1],
                            data[2]);
    }

    for (iedge = 1; iedge <= nedge; iedge++) {
        if (ignoreEdge[iedge-1] == 1) continue;

        fprintf(fpBEM, "$ edge %d", iedge);
        status = EG_attributeNum(eedges[iedge-1], &nattr);
        CHECK_STATUS(EG_attributeNum);

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(eedges[iedge-1], iattr,
                                     &aname, &atype, &alen,
                                     &tempIlist, &tempRlist, &tempClist);
            CHECK_STATUS(EG_attributeGet);

            if        (atype == ATTRINT) {
                fprintf(fpBEM, " %s=", aname);
                for (i = 0; i < alen; i++) {
                    fprintf(fpBEM, "%d;", tempIlist[i]);
                }
            } else if (atype == ATTRREAL) {
                fprintf(fpBEM, " %s=", aname);
                for (i = 0; i < alen; i++) {
                    fprintf(fpBEM, "%f;", tempRlist[i]);
                }
            } else if (atype == ATTRSTRING) {
                fprintf(fpBEM, " %s=%s", aname, tempClist);
            }
        }
        fprintf(fpBEM, "\n");

        status = EG_getTessEdge(etess, iedge,
                                &npnt, &xyz, &uv);
        CHECK_STATUS(EG_getTessEdge);

        if (npnt <= 0) continue;

        MALLOC(edgeGID[iedge], int, npnt);

        status = EG_getTopology(eedges[iedge-1], &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        ibeg = status = EG_indexBodyTopo(ebody, echilds[0]);
        CHECK_STATUS(EG_indexBodyTopo);

        iend = status = EG_indexBodyTopo(ebody, echilds[1]);
        CHECK_STATUS(EG_indexBodyTopo);

        edgeGID[iedge][0] = nodeGID[ibeg];

        for (ipnt = 1; ipnt < npnt-1; ipnt++) {
            ngid++;
            edgeGID[iedge][ipnt] = ngid;

            fprintf(fpBEM, "GRID    %8d %7d %7.4f %7.4f %7.4f\n",
                    ngid, ncid, xyz[3*ipnt  ],
                                xyz[3*ipnt+1],
                                xyz[3*ipnt+2]);
        }


        edgeGID[iedge][npnt-1] = nodeGID[iend];
    }

    for (iface = 1; iface <= nface; iface++) {
        if (ignoreFace[iface-1] == 1) continue;

        fprintf(fpBEM, "$ face %d", iface);
        status = EG_attributeNum(efaces[iface-1], &nattr);
        CHECK_STATUS(EG_attributeNum);

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(efaces[iface-1], iattr,
                                     &aname, &atype, &alen,
                                     &tempIlist, &tempRlist, &tempClist);
            CHECK_STATUS(EG_attributeGet);

            if        (atype == ATTRINT) {
                fprintf(fpBEM, " %s=", aname);
                for (i = 0; i < alen; i++) {
                    fprintf(fpBEM, "%d;", tempIlist[i]);
                }
            } else if (atype == ATTRREAL) {
                fprintf(fpBEM, " %s=", aname);
                for (i = 0; i < alen; i++) {
                    fprintf(fpBEM, "%f;", tempRlist[i]);
                }
            } else if (atype == ATTRSTRING) {
                fprintf(fpBEM, " %s=%s", aname, tempClist);
            }
        }
        fprintf(fpBEM, "\n");

//$$$        status = EG_getQuads(etess, iface,
//$$$                             &npnt, &xyz, &uv, &ptype, &pindx,
//$$$                             &npatch);
//$$$        if (status < EGADS_SUCCESS) goto cleanup;
//$$$
//$$$        if (npatch <= 0) {
        status = EG_getTessFace(etess, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);
//$$$        }

        MALLOC(faceGID[iface], int, npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (ptype[ipnt] == 0) {
                faceGID[iface][ipnt] = nodeGID[pindx[ipnt]];
            } else if (ptype[ipnt] > 0) {
                faceGID[iface][ipnt] = edgeGID[pindx[ipnt]][ptype[ipnt]-1];
            } else {
                ngid++;
                faceGID[iface][ipnt] = ngid;

                fprintf(fpBEM, "GRID    %8d %7d %7.4f %7.4f %7.4f\n",
                        ngid, ncid, xyz[3*ipnt  ],
                                    xyz[3*ipnt+1],
                                    xyz[3*ipnt+2]);
            }
        }
    }

    /* write out the CROD cards associated with all Edges */
    if (NOCROD(numUdp) == 0) {
        for (iedge = 1; iedge <= nedge; iedge++) {
            if (ignoreEdge[iedge-1] == 1) continue;

            npid++;
            fprintf(fpBEM, "$ edge %d\n", iedge);

            status = EG_getTessEdge(etess, iedge,
                                    &npnt, &xyz, &uv);
            CHECK_STATUS(EG_getTessEdge);

            for (ipnt = 1; ipnt < npnt; ipnt++) {
                neid++;
#ifndef __clang_analyzer__
                fprintf(fpBEM, "CROD    %8d %7d %7d %7d\n",
                        neid, npid, edgeGID[iedge][ipnt-1],
                                    edgeGID[iedge][ipnt  ]);
#endif
            }
        }
    }

    /* write out the CQUAD4 cards associated with all Faces */
    for (iface = 1; iface <= nface; iface++) {
        if (ignoreFace[iface-1] == 1) continue;

        npid++;
        fprintf(fpBEM, "$ face %d\n", iface);

        status = EG_getTessFace(etess, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        for (itri = 0; itri < ntri; itri+=2) {
            neid++;
            fprintf(fpBEM, "CQUAD4  %8d %7d %7d %7d %7d %7d\n",
                        neid, npid, faceGID[iface][tris[3*itri  ]-1],
                                    faceGID[iface][tris[3*itri+1]-1],
                                    faceGID[iface][tris[3*itri+2]-1],
                                    faceGID[iface][tris[3*itri+5]-1]);
        }
    }

    /* write out the CQUAD4 or CTRIA3 cards associated with all Faces */
//$$$    for (iface = 1; iface <= nface; iface++) {
//$$$        npid++;
//$$$        fprintf(fpBEM, "$ face %d\n", iface);
//$$$
//$$$        status = EG_getQuads(etess, iface,
//$$$                             &npnt, &xyz, &uv, &ptype, &pindx,
//$$$                             &npatch);
//$$$        if (status < EGADS_SUCCESS) goto cleanup;
//$$$
//$$$        if (npatch > 0) {
//$$$            for (ipatch = 1; ipatch <= npatch; ipatch++) {
//$$$                status = EG_getPatch(etess, iface, ipatch,
//$$$                                     &n1, &n2, &pvindex, &pbounds);
//$$$                if (status < EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                for (i2 = 1; i2 < n2; i2++) {
//$$$                    for (i1 = 1; i1 < n1; i1++) {
//$$$                        neid++;
//$$$                        fprintf(fpBEM, "CQUAD4  %8d %7d %7d %7d %7d %7d\n",
//$$$                                neid, npid, faceGID[iface][pvindex[(i1-1)+n1*(i2-1)]-1],
//$$$                                            faceGID[iface][pvindex[(i1  )+n1*(i2-1)]-1],
//$$$                                            faceGID[iface][pvindex[(i1  )+n1*(i2  )]-1],
//$$$                                            faceGID[iface][pvindex[(i1-1)+n1*(i2  )]-1]);
//$$$                    }
//$$$                }
//$$$
//$$$            }
//$$$        } else {
//$$$            status = EG_getTessFace(etess, iface,
//$$$                                    &npnt, &xyz, &uv, &ptype, &pindx,
//$$$                                    &ntri, &tris, &tric);
//$$$            if (status < EGADS_SUCCESS) goto cleanup;
//$$$
//$$$            for (itri = 0; itri < ntri; itri++) {
//$$$                neid++;
//$$$                fprintf(fpBEM, "CTRIA3  %8d %7d %7d %7d %7d\n",
//$$$                        neid, npid, faceGID[iface][tris[3*itri  ]-1],
//$$$                                    faceGID[iface][tris[3*itri+1]-1],
//$$$                                    faceGID[iface][tris[3*itri+2]-1]);
//$$$            }
//$$$        }
//$$$    }

    /* write out the PSHELL cards associated with all Faces */
    fprintf(fpBEM, "$ properties and materials\n");

    if (NOCROD(numUdp) == 0) {
        for (ipid = nedge+1; ipid <= npid; ipid++) {
            fprintf(fpBEM, "PSHELL  %8d       1     1.0\n", ipid);
        }
    } else {
        for (ipid = 1; ipid <= npid; ipid++) {
            fprintf(fpBEM, "PSHELL  %8d       1     1.0\n", ipid);
        }
    }

    /* write out the MAT1 card */
    fprintf(fpBEM, "$\n");
    fprintf(fpBEM, "MAT1           1     30.      9.     1.0     1.0\n");

    /* finialize and close the file */
    fclose(fpBEM);
    fpBEM = NULL;

cleanup:
    FREE(enodes);
    FREE(eedges);
    FREE(efaces);

    if (fpBEM != NULL) {
        fprintf(fpBEM, "$$$ error encountered during writeBEM\n");
        fclose(fpBEM);
    }

    if (faceGID != NULL) {
        for (iface = 1; iface <= nface; iface++) {
            FREE(faceGID[iface]);
        }
        FREE(faceGID);
    }

    if (edgeGID != NULL) {
        for (iedge = 1; iedge <= nedge; iedge++) {
            FREE(edgeGID[iedge]);
        }
        FREE(edgeGID);
    }

    FREE(nodeGID);

    FREE(iwest );
    FREE(inorth);
    FREE(ieast );
    FREE(isouth);
    FREE(rpos  );
    FREE(npnts );

    FREE(ignoreNode);
    FREE(ignoreEdge);
    FREE(ignoreFace);

    return status;
}
