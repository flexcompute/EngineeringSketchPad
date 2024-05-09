/*
 ************************************************************************
 *                                                                      *
 * udfDeform -- deform Bsplines on Body                                 *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2024  John F. Dannenhoffer, III (Syracuse University)
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
#define IFACE(IUDP,I)  ((int    *) (udps[IUDP].arg[0].val))[I]
#define IU(   IUDP,I)  ((int    *) (udps[IUDP].arg[1].val))[I]
#define IV(   IUDP,I)  ((int    *) (udps[IUDP].arg[2].val))[I]
#define DIST( IUDP,I)  ((double *) (udps[IUDP].arg[3].val))[I]
#define TOLER(IUDP  )  ((double *) (udps[IUDP].arg[4].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"iface", "iu",    "iv",    "dist",   "toler",  };
static int    argTypes[NUMUDPARGS] = {ATTRINT, ATTRINT, ATTRINT, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,        -1,     -1,      0,        0,        };
static double argDdefs[NUMUDPARGS] = {0.,       0.,      0.,     0.,       0.,       };

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
udpExecute(ego  emodel,                 /* (in)  Model containing Body */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, nchild, *senses, *header=NULL;
    int     i, j, ij, ij2, nu, nv, nnode, nedge, nface, iface;
    int     *ireplace=NULL;
    double  data[18], *rdata=NULL, *norm=NULL, uvec[3], vvec[3], size;
    ego     *ebodys, *enodes=NULL, *eedges=NULL, *efaces=NULL, esurface, eref, *echilds, context, enew;
    ego     *ereplace=NULL;
    udp_T   *udps = *Udps;

    char    *message=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("iface(0) =");
    for (i = 0; i < udps[0].arg[0].size; i++) {
        printf(" %10d", IFACE(0,i));
    }
    printf("\n");
    printf("iu(   0) =");
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf(" %10d", IU(0,i));
    }
    printf("\n");
    printf("iv(   0) =");
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf(" %10d", IV(0,i));
    }
    printf("\n");
    printf("idist(0) =");
    for (i = 0; i < udps[0].arg[0].size; i++) {
        printf(" %10.6f", DIST(0,i));
    }
    printf("\n");
    printf("toler(0) = %e\n", TOLER(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 100, "expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, &enodes);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebodys[0], NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    /* check arguments */
    if        (udps[0].arg[0].size != udps[0].arg[1].size) {
        snprintf(message, 100, "\"iface\" and \"iu\" should be same length");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[0].size != udps[0].arg[2].size) {
        snprintf(message, 100, "\"iface\" and \"iv\" should be same length");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[0].size != udps[0].arg[3].size) {
        snprintf(message, 100, "\"iface\" and \"dist\" should be same length");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (TOLER(0) < 0) {
        snprintf(message, 100, "\"toler\" must be non-negative");
        status = EGADS_RANGERR;
        goto cleanup;

    }

    SPLINT_CHECK_FOR_NULL(enodes);
    SPLINT_CHECK_FOR_NULL(eedges);
    SPLINT_CHECK_FOR_NULL(efaces);

    for (i = 0; i < udps[0].arg[0].size; i++) {
        iface = IFACE(0,i);

        if (iface < 1 || iface > nface) {
            snprintf(message, 100, "\"iface[%d]\" (%d) is not between 1 and %d", i, iface, nface);
            status = EGADS_RANGERR;
            goto cleanup;
        }

        status = EG_getTopology(efaces[iface-1], &esurface, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getGeometry(esurface, &oclass, &mtype, &eref, &header, &rdata);
        CHECK_STATUS(EG_getGeometry);

        if (mtype != BSPLINE) {
            snprintf(message, 100, "\"iface[%d]\" (%d) is not a BSPLINE", i+1, iface);
            status = EGADS_RANGERR;
            goto cleanup;
        }

        SPLINT_CHECK_FOR_NULL(header);
        SPLINT_CHECK_FOR_NULL(rdata );

        if        (IU(0,i) < 1 || IU(0,i) > header[2]) {
            snprintf(message, 100, "\"iu[%d]\" (%d) is not between 1 and %d", i+1, IU(0,i), header[2]);
            status = EGADS_RANGERR;
            goto cleanup;
        } else if (IV(0,i) < 1 || IV(0,i) > header[5]) {
            snprintf(message, 100, "\"iv[%d]\" (%d) is not between 1 and %d", i+1, IV(0,i), header[5]);
            status = EGADS_RANGERR;
            goto cleanup;
        }

        if (header != NULL) {EG_free(header); header = NULL;}
        if (rdata  != NULL) {EG_free(rdata ); rdata  = NULL;}
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("iface(%d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[0].size; i++) {
        printf(" %10d", IFACE(numUdp,i));
    }
    printf("\n");
    printf("iu(   %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf(" %10d", IU(numUdp,i));
    }
    printf("\n");
    printf("iv(   %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[2].size; i++) {
        printf(" %10d", IV(numUdp,i));
    }
    printf("\n");
    printf("idist(%d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[3].size; i++) {
        printf(" %10.6f", DIST(numUdp,i));
    }
    printf("\n");
    printf("toler(%d) = %e\n", numUdp, TOLER(numUdp));
#endif

    /* get the context associated with the model */
    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* make a list to store the original and replacement Faces */
    MALLOC(ireplace, int, nface);
    MALLOC(ereplace, ego, nface);

    /* loop through IFACE to find the Faces that should be
       replaced */
    for (iface = 1; iface <= nface; iface++) {
        ireplace[iface-1] = 0;
    }

    for (i = 0; i < udps[numUdp].arg[0].size; i++) {
        ireplace[IFACE(numUdp,i)-1]++;
    }

    /* create the replacement Faces */
    for (iface = 1; iface <= nface; iface++) {

        /* use the current Face if not modified */
        if (ireplace[iface-1] == 0) {
            status = EG_copyObject(efaces[iface-1], NULL, &ereplace[iface-1]);
            CHECK_STATUS(EG_copyObject);

            continue;
        }

        /* get the surface associated with iface */
        status = EG_getTopology(efaces[iface-1], &esurface, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getGeometry(esurface, &oclass, &mtype, &eref, &header, &rdata);
        CHECK_STATUS(EG_getGeometry);

        SPLINT_CHECK_FOR_NULL(header);
        SPLINT_CHECK_FOR_NULL(rdata );

        nu = header[2];
        nv = header[5];

#ifdef DEBUG
        int k, m=0;
        for (k = 0; k < header[3]; k++) {
            printf("uknot[%2d] = %12.5f\n", k, rdata[m++]);
        }
        for (k = 0; k < header[6]; k++) {
            printf("vknot[%2d] = %12.5f\n", k, rdata[m++]);
        }
        for (i = 0; i < nu; i++) {
            for (j = 0; j < nv; j++) {
                printf("cp[%2d,%2d] = %12.5f %12.5f %12.5f\n", i, j, rdata[m], rdata[m+1], rdata[m+2]);
                m += 3;
            }
        }
#endif

        /* get the normals to the control net */
        MALLOC(norm, double, 3*nu*nv);

        for (i = 0; i < nu; i++) {
            for (j = 0; j < nv; j++) {
                ij  = 3 * (i + j * nu);
                ij2 = ij + header[3] + header[6];

                if        (i == 0) {
                    uvec[0] = rdata[ij2+3] - rdata[ij2  ];
                    uvec[1] = rdata[ij2+4] - rdata[ij2+1];
                    uvec[2] = rdata[ij2+5] - rdata[ij2+2];
                } else if (i == nu-1) {
                    uvec[0] = rdata[ij2  ] - rdata[ij2-3];
                    uvec[1] = rdata[ij2+1] - rdata[ij2-2];
                    uvec[2] = rdata[ij2+2] - rdata[ij2-1];
                } else {
                    uvec[0] = rdata[ij2+3] - rdata[ij2-3];
                    uvec[1] = rdata[ij2+4] - rdata[ij2-2];
                    uvec[2] = rdata[ij2+5] - rdata[ij2-1];
                }

                if        (j == 0) {
                    vvec[0] = rdata[ij2+3*nu  ] - rdata[ij2       ];
                    vvec[1] = rdata[ij2+3*nu+1] - rdata[ij2     +1];
                    vvec[2] = rdata[ij2+3*nu+2] - rdata[ij2     +2];
                } else if (j == nv-1) {
                    vvec[0] = rdata[ij2       ] - rdata[ij2-3*nu  ];
                    vvec[1] = rdata[ij2     +1] - rdata[ij2-3*nu+1];
                    vvec[2] = rdata[ij2     +2] - rdata[ij2-3*nu+2];
                } else {
                    vvec[0] = rdata[ij2+3*nu  ] - rdata[ij2-3*nu  ];
                    vvec[1] = rdata[ij2+3*nu+1] - rdata[ij2-3*nu+1];
                    vvec[2] = rdata[ij2+3*nu+2] - rdata[ij2-3*nu+2];
                }

                norm[ij  ] = uvec[1] * vvec[2] - uvec[2] * vvec[1];
                norm[ij+1] = uvec[2] * vvec[0] - uvec[0] * vvec[2];
                norm[ij+2] = uvec[0] * vvec[1] - uvec[1] * vvec[0];

                size = sqrt(norm[ij]*norm[ij] + norm[ij+1]*norm[ij+1] + norm[ij+2]*norm[ij+2]);

                norm[ij  ] /= size;
                norm[ij+1] /= size;
                norm[ij+2] /= size;
#ifdef DEBUG
                printf("norm[%2d,%2d] = %12.5f %12.5f %12.5f\n", i, j, norm[ij], norm[ij+1], norm[ij+2]);
#endif
            }
        }

        /* loop through the inputs and modify the appropriate control points */
        for (i = 0; i < udps[numUdp].arg[0].size; i++) {
            if (IFACE(numUdp,i) != iface) continue;

            ij  = 3 * ((IU(numUdp,i)-1) + nu * (IV(numUdp,i)-1));
            ij2 = ij + header[3] + header[6];

            rdata[ij2  ] += norm[ij  ] * DIST(numUdp,i);
            rdata[ij2+1] += norm[ij+1] * DIST(numUdp,i);
            rdata[ij2+2] += norm[ij+2] * DIST(numUdp,i);
        }

        FREE(norm);

        /* make the new surface and Face */
        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rdata, &esurface);
        CHECK_STATUS(EG_makeGeometry);

        status = EG_makeFace(esurface, SFORWARD, data, &ereplace[iface-1]);
        CHECK_STATUS(EG_makeFace);

        status = EG_attributeDup(efaces[iface-1], ereplace[iface-1]);
        CHECK_STATUS(EG_attributeDup);

        if (header != NULL) {EG_free(header); header = NULL;}
        if (rdata  != NULL) {EG_free(rdata ); rdata  = NULL;}
    }

    /* make the new Body with the replacement Faces */
    status = EG_sewFaces(nface, ereplace, TOLER(numUdp), 0, &enew);
    CHECK_STATUS(EG_sewFaces);

    status = EG_getTopology(enew, &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (nchild != 1) {
        snprintf(message, 100, "sewing Faces yielded %d Bodys (not 1)", nchild);
        status = OCSM_UDP_ERROR1;
        goto cleanup;
    }

    /* remember the Body and delete the Model created by sewing */
    status = EG_copyObject(echilds[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);

    status = EG_deleteObject(enew);
    CHECK_STATUS(EG_deleteObject);

    SPLINT_CHECK_FOR_NULL(*ebody);

    /* tell OpenCSM to put _body, _brch, and Branch Attributes on the Faces */
    status = EG_attributeAdd(*ebody, "__markFaces__", ATTRSTRING, 1, NULL, NULL, "true");
    CHECK_STATUS(EG_attributeAdd);

    /* set the output value */

    /* the Body is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(ireplace);
    FREE(ereplace);

    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (efaces != NULL) EG_free(efaces);

    if (header != NULL) EG_free(header);
    if (rdata != NULL)  EG_free(rdata );

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

