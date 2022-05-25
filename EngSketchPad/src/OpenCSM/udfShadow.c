/*
 ************************************************************************
 *                                                                      *
 * udfShadow -- find mass properties for shadow of a Body               *
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

#define NUMUDPARGS 7
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define NUMPTS(IUDP)  ((int    *) (udps[IUDP].arg[0].val))[0]
#define AREA(  IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define XCENT( IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define YCENT( IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define IXX(   IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define IXY(   IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]
#define IYY(   IUDP)  ((double *) (udps[IUDP].arg[6].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"numpts", "area",   "xcent",    "ycent",   "ixx",     "ixy",     "iyy",     };
static int    argTypes[NUMUDPARGS] = {ATTRINT,  -ATTRREAL, -ATTRREAL, -ATTRREAL, -ATTRREAL, -ATTRREAL, -ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {1001,      0,         0,         0,         0,         0,        0,         };
static double argDdefs[NUMUDPARGS] = {0.,        0.,        0.,        0.,        0.,        0.,       0.,        };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"

#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define           EPS06           1.0e-6


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

    int     oclass, mtype, nchild, *senses, i, j, k;
    int     ip0, ix0, iy0, ip1, ix1, iy1, ip2, ix2, iy2, ix3, iy3, iswap;
    int     nface, iface, itri, ntri, npnt;
    CINT    *ptype, *pindx, *tris, *tric;
    double  data[18], bbox[6], size, params[3];
    double  slp_left, slp_rite, x_left, x_rite;
    double  dx, dy, area, xcent, ycent, Ixx, Ixy, Iyy;
    CDOUBLE *xyz, *uv;
    ego     context, eref, *ebodys, etess=NULL;
    char    *gp=NULL, *message=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("numpts(0) = %d\n", NUMPTS(0));
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
    if (udps[0].arg[0].size != 1) {
        snprintf(message, 100, "numpts should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (NUMPTS(0) < 10) {
        snprintf(message, 100, "numpts should be greater than 10");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("numpts(%d) = %d\n", numUdp, NUMPTS(numUdp));
#endif

    /* make a copy of the Body (so that it does not get removed
     when OpenCSM deletes emodel) */
    status = EG_copyObject(ebodys[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* find the bounding box forebodys[0] */
    status = EG_getBoundingBox(ebodys[0], bbox);
    CHECK_STATUS(EG_getBoundingBox);

    /* get a piece of GraphPaper on which the Triangles will be drawn */
    MALLOC(gp, char, NUMPTS(0)*NUMPTS(0));

    /* initialize the GraphPaper */
    k = 0;
    for (j = 0; j < NUMPTS(0); j++) {
        for (i = 0; i < NUMPTS(0); i++) {
            gp[k++] = 0;
        }
    }

    /* tessellate ebody */
    size = sqrt(SQR(bbox[3]-bbox[0]) + SQR(bbox[4]-bbox[1]) + SQR(bbox[5]-bbox[2]));

    params[0] = 0.0250  * size;
    params[1] = 0.0075  * size;
    params[2] = 20;
    status = EG_makeTessBody(ebodys[0], params, &etess);
    CHECK_STATUS(EG_makeTessBody);

    SPLINT_CHECK_FOR_NULL(etess);

    /* loop through all the Triangles of all the Face and draw the projection
       of the Triangle on the Graphpaper */
    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, NULL);
    CHECK_STATUS(EG_getBodyTopos);

    for (iface = 1; iface <= nface; iface++) {
        status = EG_getTessFace(etess, iface,
                                &npnt, &xyz, &uv, &ptype, &pindx,
                                &ntri, &tris, &tric);
        CHECK_STATUS(EG_getTessFace);

        for (itri = 0; itri < ntri; itri++) {
            ip0 = tris[3*itri  ] - 1;
            ip1 = tris[3*itri+1] - 1;
            ip2 = tris[3*itri+2] - 1;

            ix0 = (NUMPTS(0) - 1) * (xyz[3*ip0  ] - bbox[0]) / (bbox[3] - bbox[0]);
            ix1 = (NUMPTS(0) - 1) * (xyz[3*ip1  ] - bbox[0]) / (bbox[3] - bbox[0]);
            ix2 = (NUMPTS(0) - 1) * (xyz[3*ip2  ] - bbox[0]) / (bbox[3] - bbox[0]);

            iy0 = (NUMPTS(0) - 1) * (xyz[3*ip0+1] - bbox[1]) / (bbox[4] - bbox[1]);
            iy1 = (NUMPTS(0) - 1) * (xyz[3*ip1+1] - bbox[1]) / (bbox[4] - bbox[1]);
            iy2 = (NUMPTS(0) - 1) * (xyz[3*ip2+1] - bbox[1]) / (bbox[4] - bbox[1]);

            /* sort points in increasing iy */
            if (iy1 < iy0) {
                iswap = ix0;
                ix0   = ix1;
                ix1   = iswap;

                iswap = iy0;
                iy0   = iy1;
                iy1   = iswap;
            }
            if (iy2 < iy0) {
                iswap = ix0;
                ix0   = ix2;
                ix2   = iswap;

                iswap = iy0;
                iy0   = iy2;
                iy2   = iswap;
            }
            if (iy2 < iy1) {
                iswap = ix1;
                ix1   = ix2;
                ix2   = iswap;

                iswap = iy1;
                iy1   = iy2;
                iy2   = iswap;
            }

            /* collapsed triangle */
            if (iy0 == iy2) {

            /* triangle with flat bottom (0,1=2) */
            } else if (iy1 == iy2) {
                if (ix1 > ix2) {
                    iswap = ix1;
                    ix1   = ix2;
                    ix2   = iswap;
                }

                slp_left = (double)(ix1 - ix0) / (double)(iy1 - iy0);
                slp_rite = (double)(ix2 - ix0) / (double)(iy2 - iy0);

                x_left = ix0;
                x_rite = ix0;

                for (j = MAX(iy0,0); j <= MIN(iy1,NUMPTS(0)-1); j++) {
                    for (i = MAX((int)x_left,0); i <= MIN((int)x_rite,NUMPTS(0)-1); i++) {
                        gp[i+j*NUMPTS(0)] = 1;
                    }

                    x_left += slp_left;
                    x_rite += slp_rite;
                }

            /* triangle with flat top (0=1,2) */
            } else if (iy0 == iy1) {
                if (ix0 < ix1) {
                    iswap = ix0;
                    ix0   = ix1;
                    ix1   = iswap;
                }

                slp_left = (double)(ix2 - ix1) / (double)(iy2 - iy1);
                slp_rite = (double)(ix2 - ix0) / (double)(iy2 - iy0);

                x_left = ix2;
                x_rite = ix2;

                for (j = MAX(iy2,0); j >= MIN(iy0,NUMPTS(0)-1); j--) {
                    for (i = MAX((int)x_left,0); i <= MIN((int)x_rite,NUMPTS(0)-1); i++) {
                        gp[i+j*NUMPTS(0)] = 1;
                    }

                    x_left -= slp_left;
                    x_rite -= slp_rite;
                }

            /* general triangle */
            } else {
                /* point on line between 0 and 2 */
                ix3 = ix0 + (double)(iy1 - iy0) / (double)(iy2 - iy0) * (double)(ix2 - ix0);
                iy3 = iy1;

                /* top half-triangle with flat bottom (0,1=3) */
                if (ix1 > ix3) {
                    iswap = ix1;
                    ix1   = ix3;
                    ix3   = iswap;
                }

                slp_left = (double)(ix1 - ix0) / (double)(iy1 - iy0);
                slp_rite = (double)(ix3 - ix0) / (double)(iy3 - iy0);

                x_left = ix0;
                x_rite = ix0;

                for (j = MAX(iy0,0); j <= MIN(iy1, NUMPTS(0)-1); j++) {
                    for (i = MAX((int)x_left,0); i <= MIN((int)x_rite,NUMPTS(0)-1); i++) {
                        gp[i+j*NUMPTS(0)] = 1;
                    }

                    x_left += slp_left;
                    x_rite += slp_rite;
                }

                /* bottom half-triangle with flat top (3=1,2) */
                slp_left = (double)(ix2 - ix1) / (double)(iy2 - iy1);
                slp_rite = (double)(ix2 - ix3) / (double)(iy2 - iy3);

                x_left = ix2;
                x_rite = ix2;

                for (j = MAX(iy2,0); j >= MIN(iy3,NUMPTS(0)-1); j--) {
                    for (i = MAX((int)x_left,0); i <= MIN((int)x_rite,NUMPTS(0)-1); i++) {
                        gp[i+j*NUMPTS(0)] = 1;
                    }

                    x_left -= slp_left;
                    x_rite -= slp_rite;
                }
            }
        }
    }

#ifdef DEBUG
    /* print the GraphPaper */
    k = 0;
    for (j = 0; j < NUMPTS(0); j++) {
        for (i = 0; i < NUMPTS(0); i++) {
            if (gp[k++] == 0) {
                printf(".");
            } else {
                printf("X");
            }
        }
        printf("\n");
    }
#endif

    /* accumulate the mass property sums */
    area  = 0;
    xcent = 0;
    ycent = 0;
    Ixx   = 0;
    Ixy   = 0;
    Iyy   = 0;

    k = 0;
    for (j = 0; j < NUMPTS(0); j++) {
        for (i = 0; i < NUMPTS(0); i++) {
            if (gp[k++] != 0) {
                area  += 1;
                xcent += i;
                ycent += j;
                Ixx   += j * j;
                Ixy   -= i * j;
                Iyy   += i * i;
            }
        }
    }

    /* compute the mass properties */
    dx = (bbox[3] - bbox[0]) / (NUMPTS(0) - 1);
    dy = (bbox[4] - bbox[1]) / (NUMPTS(0) - 1);

    area  *= dx * dy;
    xcent *= dx * dy * dx;
    ycent *= dx * dy * dy;
    Ixx   *= dx * dy * dy * dy;
    Ixy   *= dx * dy * dx * dy;
    Iyy   *= dx * dy * dx * dx;

    xcent = xcent / area + bbox[0];
    ycent = ycent / area + bbox[1];
    Ixx   = Ixx   - area * (ycent-bbox[1]) * (ycent-bbox[1]);
    Ixy   = Ixy   + area * (xcent-bbox[0]) * (ycent-bbox[1]);
    Iyy   = Iyy   - area * (xcent-bbox[0]) * (xcent-bbox[0]);

    /* set the output values */
    AREA( 0) = area;
    XCENT(0) = xcent;
    YCENT(0) = ycent;
    IXX(  0) = Ixx;
    IXY(  0) = Ixy;
    IYY(  0) = Iyy;

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(gp);

    if (etess != NULL) EG_deleteObject(etess);

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
