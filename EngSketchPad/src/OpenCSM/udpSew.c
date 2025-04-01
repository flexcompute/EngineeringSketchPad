/*
 ************************************************************************
 *                                                                      *
 * udpSew -- udp file to sew faces in a STEP file into a SOLIDBODY      *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2025  John F. Dannenhoffer, III (Syracuse University)
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
#define FILENAME(  IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define TOLER(     IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define BODYNUM(   IUDP)  ((int    *) (udps[IUDP].arg[2].val))[0]
#define REMOVEDUPS(IUDP)  ((int    *) (udps[IUDP].arg[3].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename",  "toler",  "bodynum", "removedups", };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING,  ATTRREAL, ATTRINT,   ATTRINT,      };
static int    argIdefs[NUMUDPARGS] = {0,           0,        0,         0,            };
static double argDdefs[NUMUDPARGS] = {0.,          0.,       0.,        0.,           };

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

    int     oclass, mtype, nchild1, ichild1, *senses, nface, mface;
    int     nchild2, ichild2, nchild3, ichild3, nchild4, ichild4, nchild5;
    int     nman, nnon, iedge, nedge, ibody, iface, jface;
    int     *dupface=NULL;
    double  data[4], bboxi[6], bboxj[6], massi[14], massj[14];
    ego     emodel1, emodel2, eref, *eedges;
    ego     *ebodys1, *ebodys2, *ebodys3, *ebodys4, *ebodys5;
    ego     *faceList=NULL;
    ego     topRef, prev, next;
    udp_T   *udps = *Udps;
    void    *realloc_temp=NULL;              /* used by RALLOC macro */

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename      = %s\n", FILENAME(  0));
    printf("toler(0)      = %f\n", TOLER(     0));
    printf("bodynum(0)    = %d\n", BODYNUM(   0));
    printf("removedups(0) = %d\n", REMOVEDUPS(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if (udps[0].arg[1].size > 1) {
        printf(" udpExecute: toler should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (TOLER(0) < 0) {
        printf(" udpExecute: toler = %f < 0\n", TOLER(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (BODYNUM(0) < -1) {
        printf(" udpExecute: bodynum = %d < -1\n", BODYNUM(0));
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d)    = %s\n", numUdp, FILENAME(  numUdp));
    printf("toler(%d)       = %f\n", numUdp, TOLER(     numUdp));
    printf("bodynum(%d)     = %d\n", numUdp, BODYNUM(   numUdp));
    printf("removedups(%d)  = %d\n", numUdp, REMOVEDUPS(numUdp));
#endif

    /* load the model */
    status = EG_loadModel(context, 0, FILENAME(0), &emodel1);
    if ((status != EGADS_SUCCESS) || (emodel1 == NULL)) {
        goto cleanup;
    }

    /* make a list of the Faces in emodel1 */
    mface = 50;
    nface = 0;

    MALLOC(faceList, ego, mface);
    MALLOC(dupface,  int, mface);

    status = EG_getTopology(emodel1, &eref, &oclass, &mtype,
                            data, &nchild1, &ebodys1, &senses);
    if (status != EGADS_SUCCESS) printf("EG_getTopology -> status=%d\n", status);

    for (ichild1 = 0; ichild1 < nchild1; ichild1++) {
        status = EG_getTopology(ebodys1[ichild1], &eref, &oclass, &mtype,
                                data, &nchild2, &ebodys2, &senses);
        if (status != EGADS_SUCCESS) printf("EG_getTopology -> status=%d\n", status);

        if (oclass == FACE) {
            if (nface >= mface) {
                mface += 50;
                RALLOC(faceList, ego, mface);
                RALLOC(dupface,  int, mface);
            }
            faceList[nface++] = ebodys1[ichild1];
            continue;
        }

        for (ichild2 = 0; ichild2 < nchild2; ichild2++) {
            status = EG_getTopology(ebodys2[ichild2], &eref, &oclass, &mtype,
                                    data, &nchild3, &ebodys3, &senses);
            if (status != EGADS_SUCCESS) printf("EG_getTopology -> status=%d\n", status);

            if (oclass == FACE) {
                if (nface >= mface) {
                    mface += 50;
                    RALLOC(faceList, ego, mface);
                    RALLOC(dupface,  int, mface);
                }
                faceList[nface++] = ebodys2[ichild2];
                continue;
            }

            for (ichild3 = 0; ichild3 < nchild3; ichild3++) {
                status = EG_getTopology(ebodys3[ichild3], &eref, &oclass, &mtype,
                                        data, &nchild4, &ebodys4, &senses);
                if (status != EGADS_SUCCESS) printf("EG_getTopology -> status=%d\n", status);

                if (oclass == FACE) {
                    if (nface >= mface) {
                        mface += 50;
                        RALLOC(faceList, ego, mface);
                        RALLOC(dupface,  int, mface);
                    }
                    faceList[nface++] = ebodys3[ichild3];
                    continue;
                }

                for (ichild4 = 0; ichild4 < nchild4; ichild4++) {
                    status = EG_getTopology(ebodys4[ichild4], &eref, &oclass, &mtype,
                                            data, &nchild5, &ebodys5, &senses);
                    if (status != EGADS_SUCCESS) printf("EG_getTopology -> status=%d\n", status);

                    if (oclass == FACE) {
                        if (nface >= mface) {
                            mface += 50;
                            RALLOC(faceList, ego, mface);
                            RALLOC(dupface,  int, mface);
                        }
                        faceList[nface++] = ebodys4[ichild4];
                        continue;
                    }
                }
            }
        }
    }
    printf(" udpExecute: there are %4d Faces to sew with toler=%f\n", nface, TOLER(0));

    /* check if any Faces are duplicated */
    if (REMOVEDUPS(0) > 0) {

        /* initialize the list of duplicate Faces */
        for (iface = 0; iface < nface; iface++) {
            dupface[iface] = 0;
        }

        /* duplicate Face have same bounding box and area */
        for (iface = 0; iface < nface; iface++) {
            status = EG_getBoundingBox(faceList[iface], bboxi);
            CHECK_STATUS(EG_getBoundingBox);

            for (jface = iface+1; jface < nface; jface++) {
                if (dupface[iface] > 0) continue;

                status = EG_getBoundingBox(faceList[jface], bboxj);
                CHECK_STATUS(EG_getBoundingBox);

                if (fabs(bboxi[0]-bboxj[0]) < EPS12 &&
                    fabs(bboxi[1]-bboxj[1]) < EPS12 &&
                    fabs(bboxi[2]-bboxj[2]) < EPS12 &&
                    fabs(bboxi[3]-bboxj[3]) < EPS12 &&
                    fabs(bboxi[4]-bboxj[4]) < EPS12 &&
                    fabs(bboxi[5]-bboxj[5]) < EPS12   ) {

                    status = EG_getMassProperties(faceList[iface], massi);
                    CHECK_STATUS(EG_getMassProperties);

                    status = EG_getMassProperties(faceList[jface], massj);
                    CHECK_STATUS(EG_getMassProperties);

                    if (fabs(massi[1]-massj[1]) < EPS09) {
                        printf("   Faces %3d and %3d are equivalent and will be removed\n", iface, jface);

                        dupface[iface] = 1;
                        dupface[jface] = 1;

                        break;
                    }
                }
            }
        }

        /* remove the duplicates from the list */
        for (iface = nface-1; iface >= 0; iface--) {
            if (dupface[iface] > 0) {
                faceList[iface] = faceList[nface-1];
                nface--;
            }
        }
    }

    /* sew the Faces into a new Model */
    status = EG_sewFaces(nface, faceList, TOLER(0), 1, &emodel2);
    if (status != EGADS_SUCCESS || emodel2 == NULL) {
        printf(" udpExecute: error while sewing Faces\n");
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;
    }

    status = EG_getTopology(emodel2, &eref, &oclass, &mtype,
                            data, &nchild1, &ebodys1, &senses);
    if (status != EGADS_SUCCESS || oclass != MODEL) {
        printf(" udpExecute: sewing failed\n");
        status = EGADS_NODATA;
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;

    } else {
        for (ibody = 0; ibody < nchild1; ibody++) {
            status = EG_getBodyTopos(ebodys1[ibody], NULL, FACE, &nchild2, NULL);
            if (status != EGADS_SUCCESS) {
                status = EGADS_NODATA;
                goto cleanup;
            }

            printf("             body %3d contains %5d Faces\n", ibody+1, nchild2);
        }
    }

    if (BODYNUM(0) == -1) {
        *ebody = emodel2;

        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }

        udps[numUdp].ebody = *ebody;

        goto cleanup;

    } else if (BODYNUM(0) == 0 && nchild1 > 1) {
        printf(" udpExecute: expecting emodel2 to have one child  (nchild1=%d)\n", nchild1);
        printf("             try re-running with increased toler\n");
        status = EGADS_NODATA;
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;
    } else if (BODYNUM(0) > nchild1) {
        printf(" udpExecute: (bodynum=%d) should not exceed (nchild1=%d)\n", BODYNUM(0), nchild1);
        printf("             try re-running with increased toler\n");
        status = EGADS_NODATA;
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;
    } else if (nchild1 != 1) {
        ibody = BODYNUM(0) - 1;
    } else {
        ibody = 0;
    }

    printf("             body %3d selected for processing\n", ibody+1);
    status = EG_getTopology(ebodys1[ibody], &eref, &oclass, &mtype,
                            data, &nchild2, &ebodys2, &senses);
    if (status != EGADS_SUCCESS || oclass != BODY || nchild2 != 1) {
        printf(" udpExecute: expecting ebodys1[%d] to have one child  (nchild2=%d)\n", ibody, nchild2);
        status = EGADS_NODATA;
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;
    }

    /* copy the single Body from the new Model */
    status = EG_copyObject(ebodys1[ibody], NULL, ebody);
    if (status != EGADS_SUCCESS) {
        printf(" udpExecute: problem copying BODY\n");
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;
    }
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* determine the number of maniford and non-manifold Edges */
    nman = 0;
    nnon = 0;

    status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
    if (status != EGADS_SUCCESS) {
        printf(" udpExecute: problem getting Edge information\n");
        if (emodel1 != NULL) {
            EG_deleteObject(emodel1);
            emodel1 = NULL;
        }
        if (emodel2 != NULL) {
            EG_deleteObject(emodel2);
            emodel2 = NULL;
        }
        goto cleanup;
    }

    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getInfo(eedges[iedge], &oclass, &mtype, &topRef, &prev, &next);
        if (status != EGADS_SUCCESS) {
            printf("EG_getInfo -> status=%d\n", status);
            continue;
        }

        if (mtype == DEGENERATE) continue;

        status = EG_getBodyTopos(*ebody, eedges[iedge], FACE, &nface, NULL);
        if (status == EGADS_SUCCESS) {
            if (nface == 2) {
                nman++;
            } else {
                nnon++;
            }
        } else {
            nnon++;      // problem Edges are counted as non-manifold
        }
    }

    if (eedges != NULL) EG_free(eedges);

    printf("             there are %4d manifold     Edges in Body %3d after sewing\n", nman, ibody+1);
    printf("             there are %4d non-manifold Edges in Body %3d after sewing\n", nnon, ibody+1);

    /* clean up models */
    if (emodel1 != NULL) {
        EG_deleteObject(emodel1);
        emodel1 = NULL;
    }
    if (emodel2 != NULL) {
        EG_deleteObject(emodel2);
        emodel2 = NULL;
    }

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(faceList);
    FREE(dupface );

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
