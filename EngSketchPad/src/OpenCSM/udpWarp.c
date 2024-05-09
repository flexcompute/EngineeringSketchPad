/*
 ************************************************************************
 *                                                                      *
 * udpWarp -- warp one Face of a Body via control point movement        *
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

#define NUMUDPARGS       4
#define NUMUDPINPUTBODYS 0
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define EGADSFILE(IUDP)   ((char   *) (udps[IUDP].arg[0].val))
#define IFACE(    IUDP)   ((int    *) (udps[IUDP].arg[1].val))[0]
#define DIST(     IUDP,I) ((double *) (udps[IUDP].arg[2].val))[I]
#define DIST_DOT( IUDP,I) ((double *) (udps[IUDP].arg[2].dot))[I]
#define TOLER(    IUDP)   ((double *) (udps[IUDP].arg[3].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"egadsfile", "iface", "dist",      "toler",  };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,    ATTRINT, ATTRREALSEN, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,           0,       0,           0,        };
static double argDdefs[NUMUDPARGS] = {0.,          0.,      0.,          0.,       };

/* routine to clean up private data */
#define FREEUDPDATA(A) freePrivateData(A)
static int freePrivateData(void *data);

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "egads_dot.h"


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

    int     oclass, mtype, nchild, *senses, *header=NULL;
    int     i, j, k, ij, ij2, nu, nv, nnode, inode, nedge, iedge, nface, iface;
    double  data[18], *rdata=NULL, *norm=NULL, uvec[3], vvec[3], size;
    ego     *enodes=NULL, *eedges=NULL, *efaces=NULL, esurface, eref, *echilds, emodel, enew;
    ego     *ereplace=NULL;
    char    *message=NULL;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("egadsfile(0) = %s\n", EGADSFILE(0));
    printf("iface(    0) = %d\n", IFACE(    0));
    printf("dist(     0) =");
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf(" %10.6f", DIST(0,i));
    }
    printf("\n");
    printf("toler(    0) = %f\n", TOLER(    0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check some of the arguments */
    if        (udps[0].arg[0].size == 0) {
        snprintf(message, 100, "\"egadefile\" must be specified");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (strstr(EGADSFILE(0), ".egads") == NULL) {
        snprintf(message, 100, "\"%s\" does not have \".egads\" suffix", EGADSFILE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size != 1) {
        snprintf(message, 100, "\"toler\" must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;
        
    }

    /* read the .egads file and get its Body and its Nodes, Edges, and Faces (echilds[0]) */
    status = EG_loadModel(context, 0, EGADSFILE(0), &emodel);
    if (status != EGADS_SUCCESS) {
        snprintf(message, 100, "\"%s\" could not be opened", EGADSFILE(0));
        goto cleanup;
    }

    status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (nchild != 1) {
        (void) EG_deleteObject(emodel);

        snprintf(message, 100, "\"%s\" contains %d Bodys (not 1)", EGADSFILE(0), nchild);
        status = EGADS_RANGERR;
        goto cleanup;
    }

    status = EG_getBodyTopos(echilds[0], NULL, NODE, &nnode, &enodes);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(echilds[0], NULL, EDGE, &nedge, &eedges);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(echilds[0], NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(enodes);
    SPLINT_CHECK_FOR_NULL(eedges);
    SPLINT_CHECK_FOR_NULL(efaces);

    /* remove the __trace__ attributes from the Nodes, Edges, and Faces */
    for (inode = 0; inode < nnode; inode++) {
        status = EG_attributeDel(enodes[inode], "__trace__");
        CHECK_STATUS(EG_attributeDel);
    }
    
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_attributeDel(eedges[iedge], "__trace__");
        CHECK_STATUS(EG_attributeDel);
    }
    
    for (iface = 0; iface < nface; iface++) {
        status = EG_attributeDel(efaces[iface], "__trace__");
        CHECK_STATUS(EG_attributeDel);

        status = EG_attributeDel(efaces[iface], "_hist");
        CHECK_STATUS(EG_attributeDel);
    }
    
    /* check more of the arguments */
    if (udps[0].arg[1].size != 1) {
        snprintf(message, 100, "\"iface\" must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (IFACE(0) < 1 || IFACE(0) > nface) {
        snprintf(message, 100, "\"iface\" = %d (should be between 1 and %d)", IFACE(0), nface);
        status = EGADS_RANGERR;
        goto cleanup;

    }

    /* get the Surface associated with IFACE(0) */
    status = EG_getTopology(efaces[IFACE(0)-1], &esurface, &oclass, &mtype, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getGeometry(esurface, &oclass, &mtype, &eref, &header, &rdata);
    CHECK_STATUS(EG_getGeometry);

    if (mtype != BSPLINE) {
        snprintf(message, 100, "Face[%d] is not a BSPLINE", IFACE(0));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    SPLINT_CHECK_FOR_NULL(header);
    SPLINT_CHECK_FOR_NULL(rdata );

    /* check the length of the dist argument */
    if (udps[0].arg[2].size != (header[2]-2) * (header[5]-2)) {
        snprintf(message, 100, "\"dist\" has %d entries (should be %d*%d=%d)", udps[0].arg[2].size,
                 header[2]-2, header[5]-2, (header[2]-2)*(header[5]-2));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    EG_free(header);
    EG_free(rdata );

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("egadsfile(%d) = %s\n", numUdp, EGADSFILE(numUdp));
    printf("iface(    %d) = %d\n", numUdp, IFACE(    numUdp));
    printf("dist(     %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[2].size; i++) {
        printf(" %10.6f", DIST(numUdp,i));
    }
    printf("\n");
    printf("toler(    %d) = %f\n", numUdp, TOLER(    numUdp));
#endif

    /* make a list to store the original and replacement Faces */
    MALLOC(ereplace, ego, nface);

    /* create the replacement Faces */
    for (iface = 1; iface <= nface; iface++) {

        /* use the current Face if not modified */
        if (iface != IFACE(0)) {
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
        int m=0;
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
        MALLOC(udps[numUdp].data, double, 3*nu*nv);
        norm = (double *) udps[numUdp].data;

        for (i = 0; i < nu; i++) {
            for (j = 0; j < nv; j++) {
                ij  = 3 * (i * nu + j);
                ij2 = ij + header[3] + header[6];

                if        (i == 0) {
                    vvec[0] = rdata[ij2+3*nu  ] - rdata[ij2       ];
                    vvec[1] = rdata[ij2+3*nu+1] - rdata[ij2     +1];
                    vvec[2] = rdata[ij2+3*nu+2] - rdata[ij2     +2];
                } else if (i == nu-1) {
                    vvec[0] = rdata[ij2       ] - rdata[ij2-3*nu  ];
                    vvec[1] = rdata[ij2     +1] - rdata[ij2-3*nu+1];
                    vvec[2] = rdata[ij2     +2] - rdata[ij2-3*nu+2];
                } else {
                    vvec[0] = rdata[ij2+3*nu  ] - rdata[ij2-3*nu  ];
                    vvec[1] = rdata[ij2+3*nu+1] - rdata[ij2-3*nu+1];
                    vvec[2] = rdata[ij2+3*nu+2] - rdata[ij2-3*nu+2];
                }

                if        (j == 0) {
                    uvec[0] = rdata[ij2+3] - rdata[ij2  ];
                    uvec[1] = rdata[ij2+4] - rdata[ij2+1];
                    uvec[2] = rdata[ij2+5] - rdata[ij2+2];
                } else if (j == nv-1) {
                    uvec[0] = rdata[ij2  ] - rdata[ij2-3];
                    uvec[1] = rdata[ij2+1] - rdata[ij2-2];
                    uvec[2] = rdata[ij2+2] - rdata[ij2-1];
                } else {
                    uvec[0] = rdata[ij2+3] - rdata[ij2-3];
                    uvec[1] = rdata[ij2+4] - rdata[ij2-2];
                    uvec[2] = rdata[ij2+5] - rdata[ij2-1];
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
        k = 0;
        for (i = 1; i < nu-1; i++) {
            for (j = 1; j < nv-1; j++) {
                ij  = 3 * (i *nu + j);
                ij2 = ij + header[3] + header[6];

                rdata[ij2  ] += norm[ij  ] * DIST(numUdp,k);
                rdata[ij2+1] += norm[ij+1] * DIST(numUdp,k);
                rdata[ij2+2] += norm[ij+2] * DIST(numUdp,k);

                k++;
            }
        }

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
    int status = EGADS_SUCCESS;
    
    int     iudp, judp, ipnt, nface, oclass, mtype, nchild, *senses, *header=NULL, ndata;
    int     i, j, k, m, ij;
    double  *rdata=NULL, *rdata_dot=NULL, data[18], data_dot[18], uv[2], uv_dot[2], *norm;
    ego     *efaces=NULL, esurface, *echilds, eref;

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
        status = EGADS_NOTMODEL;
        goto cleanup;
    }

    norm = (double *)(udps[iudp].data);

    /* sensitivities on IFACE */
    if (entType == OCSM_FACE && entIndex == IFACE(iudp)) {

        /* set up dots on the Surface */
        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(efaces);

        status = EG_getTopology(efaces[entIndex-1], &esurface, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_getGeometry(esurface, &oclass, &mtype, &eref, &header, &rdata);
        CHECK_STATUS(EG_getGeometry);

        SPLINT_CHECK_FOR_NULL(header);
        SPLINT_CHECK_FOR_NULL(rdata );

        ndata = header[3] + header[6] + 3 * header[2] * header[5];

        MALLOC(rdata_dot, double, ndata);

        k = 0;
        for (i = 0; i < header[3]; i++) {
            rdata_dot[k++] = 0;
        }
        for (j = 0; j < header[6]; j++) {
            rdata_dot[k++] = 0;
        }
        m = 0;
        for (i = 0; i < header[5]; i++) {
            for (j = 0; j < header[2]; j++) {
                if (i > 0 && i < header[2]-1 && j > 0 && j < header[5]-1) {
                    ij  = 3 * (i * header[2] + j);
                    rdata_dot[k++] = norm[ij  ] * DIST_DOT(iudp,m);
                    rdata_dot[k++] = norm[ij+1] * DIST_DOT(iudp,m);
                    rdata_dot[k++] = norm[ij+2] * DIST_DOT(iudp,m);
                    m++;
                } else {
                    rdata_dot[k++] = 0;
                    rdata_dot[k++] = 0;
                    rdata_dot[k++] = 0;
                }
            }
        }

        status = EG_setGeometry_dot(esurface, oclass, mtype, header, rdata, rdata_dot);
        CHECK_STATUS(EG_setGeometry_dot);

        /* evaluate the velocities */
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            uv[    0] = uvs[2*ipnt  ];
            uv[    1] = uvs[2*ipnt+1];
            uv_dot[0] = 0;
            uv_dot[1] = 0;
            
            status = EG_evaluate_dot(esurface, uv, uv_dot, data, data_dot);
            CHECK_STATUS(EG_evaluate_dot);

            vels[3*ipnt  ] = data_dot[0];
            vels[3*ipnt+1] = data_dot[1];
            vels[3*ipnt+2] = data_dot[2];
        }

    /* all other sensitivities */
    } else {
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            vels[3*ipnt  ] = 0;
            vels[3*ipnt+1] = 0;
            vels[3*ipnt+2] = 0;
        }
    }

cleanup:
    FREE(rdata_dot);

    if (efaces != NULL) EG_free(efaces);
    if (header != NULL) EG_free(header);
    if (rdata  != NULL) EG_free(rdata );
    
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   freePrivateData - free private data                                *
 *                                                                      *
 ************************************************************************
 */

static int
freePrivateData(void  *data)            /* (in)  pointer to private data */
{
    int    status = EGADS_SUCCESS;

    if (data != NULL) {
        EG_free(data);
    }

//cleanup:
    return status;
}
