/*
 ************************************************************************
 *                                                                      *
 * udpNaca456 -- udp file to generate a 4-, 5-, or 6-series airfoil     *
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

#define NUMUDPARGS 9
#include "udpUtilities.h"

/* codes for various airfoils

   mptt    -> thkcode=4,   toc=11/100,                     camcode=2,  cmax=m/100, xmaxc=p/10
   mptt-lx -> thkcode=4M,  toc=tt/100, leindex=l, xmaxt=x, camcode=2,  cmax=m/100, xmaxc=p/10
   mp0tt   -> thkcode=4,   toc=tt/100,                     camcode=3,  cl=m*.15,   xmaxc=p/20
   mp1tt   -> thkcode=4,   toc=tt/100,                     camcode=3r, cl=m*.15,   xmaxc=p/20
   63-mtt  -> thkcode=63,  toc=tt/100,                     camcode=6,  cl=m/10,    a=??
   63Amtt  -> thkcode=63a, toc=tt/100,                     camcode=6m, cl=m/10,    a=0.8
   64-mtt  -> thkcode=64,  toc=tt/100,                     camcode=6,  cl=m/10,    a=??
   64Amtt  -> thkcode=64a, toc/tt/100,                     camcode=6m, cl=m/10,    a=0.8
   65-mtt  -> thkcode=65,  toc=tt/100,                     camcode=6,  cl=m/10,    a=??
   65Amtt  -> thkcode=65a, toc/tt/100,                     camcode=6m, cl=m/10,    a=0.8
   66-mtt  -> thkcode=66,  toc=tt/100,                     camcode=6,  cl=m/10,    a=??
   67-mtt  -> thkcode=67,  toc=tt/100,                     camcode=6,  cl=m/10,    a=??
*/

/* declaration for fortran90 routine */
#ifdef WIN32
void NACA456 (int*   ithkcode, double* toc,  double* xmaxt, double* leindex,
              int*   icamcode, double* cmax, double* xmaxc, double* cl,      double* a,
              int*   nairfoil, double xairfoil[], double yairfoil[]);
#else
void naca456_(int*   ithkcode, double* toc,  double* xmaxt, double* leindex,
              int*   icamcode, double* cmax, double* xmaxc, double* cl,      double* a,
              int*   nairfoil, double xairfoil[], double yairfoil[]);
#endif

/* shorthands for accessing argument values and velocities */
#define THKCODE(      IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define TOC(          IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define XMAXT(        IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define LEINDEX(      IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define CAMCODE(      IUDP)  ((char   *) (udps[IUDP].arg[4].val))
#define CMAX(         IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]
#define XMAXC(        IUDP)  ((double *) (udps[IUDP].arg[6].val))[0]
#define CL(           IUDP)  ((double *) (udps[IUDP].arg[7].val))[0]
#define A(            IUDP)  ((double *) (udps[IUDP].arg[8].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"thkcode", "toc",     "xmaxt",  "leindex", "camcode",  "cmax",   "xmaxc",  "cl",     "a",      };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRREAL, ATTRREAL, ATTRREAL,  ATTRSTRING, ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,          0,        0,        0,         0,          0,        0,        0,        0,        };
static double argDdefs[NUMUDPARGS] = {0.,         0.,       0.,       0.,        0.,         0.,       0.,       0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define           HUGEQ           99999999.0
#define           PIo2            1.5707963267948965579989817
#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))


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

    int     ithkcode=0, icamcode=0, bluntte;
    int     nairfoil, nle;
    int     npnt, ipnt, sense[3], sizes[2], periodic;
    double  *xairfoil=NULL, *yairfoil=NULL, *pts=NULL;
    double  data[18], tdata[2], result[3], range[4], eval[18], norm[3];
    double  dxytol = 1.0e-6;
    ego     enodes[4], eedges[3], ecurve, eline, eloop, eface, enew;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */
    
#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("thkcode(0) = %s\n", THKCODE(0));
    printf("toc(0)     = %f\n", TOC(    0));
    printf("xmaxt(0)   = %f\n", XMAXT(  0));
    printf("leindex(0) = %f\n", LEINDEX(0));
    printf("camcode(0) = %s\n", CAMCODE(0));
    printf("cmax(0)    = %f\n", CMAX(   0));
    printf("xmaxc(0)   = %f\n", XMAXC(  0));
    printf("cl(0)      = %f\n", CL(     0));
    printf("a(0)       = %f\n", A(      0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* find the thickness code */
    if        (strcmp(THKCODE(0), "4"  ) == 0) {
        ithkcode = 4;
    } else if (strcmp(THKCODE(0), "4M" ) == 0) {
        ithkcode = 41;
    } else if (strcmp(THKCODE(0), "63" ) == 0) {
        ithkcode = 63;
    } else if (strcmp(THKCODE(0), "63A") == 0) {
        ithkcode = 631;
    } else if (strcmp(THKCODE(0), "64" ) == 0) {
        ithkcode = 64;
    } else if (strcmp(THKCODE(0), "64A") == 0) {
        ithkcode = 641;
    } else if (strcmp(THKCODE(0), "65" ) == 0) {
        ithkcode = 65;
    } else if (strcmp(THKCODE(0), "65A") == 0) {
        ithkcode = 651;
    } else if (strcmp(THKCODE(0), "66" ) == 0) {
        ithkcode = 66;
    } else if (strcmp(THKCODE(0), "67" ) == 0) {
        ithkcode = 67;
    } else {
        printf(" udpExecute: thkcode should be 4, 4M, 63, 63A, 64, 64A, 65, 65A, 66, or 67\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* find the camber code */
    if        (strcmp(CAMCODE(0), "0" ) == 0) {
        icamcode = 0;
    } else if (strcmp(CAMCODE(0), "2" ) == 0) {
        icamcode = 2;
    } else if (strcmp(CAMCODE(0), "3" ) == 0) {
        icamcode = 3;
    } else if (strcmp(CAMCODE(0), "3R") == 0) {
        icamcode = 31;
    } else if (strcmp(CAMCODE(0), "6" ) == 0) {
        icamcode = 6;
    } else if (strcmp(CAMCODE(0), "6M") == 0) {
        icamcode = 61;
    } else {
        printf(" udpExecute: camcode should be 0, 2, 3, 3R, 6, 6M\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* check other arguments */
    if        (udps[0].arg[1].size > 1) {
        printf(" udpExecute: toc should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (TOC(0) < 0) {
        printf(" udpExecute: toc = %f < 0\n", TOC(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        printf(" udpExecute: xmaxt should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (XMAXT(0) < 0) {
        printf(" udpExecute: xmaxt = %f < 0\n", XMAXT(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        printf(" udpExecute: leindex should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (LEINDEX(0) < 0) {
        printf(" udpExecute: leindex = %f < 0\n", LEINDEX(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        printf(" udpExecute: cmax should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[6].size > 1) {
        printf(" udpExecute: xmaxc should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (XMAXC(0) < 0) {
        printf(" udpExecute: xmaxc = %f < 0\n", XMAXC(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[7].size > 1) {
        printf(" udpExecute: cl should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (CL(0) < 0) {
        printf(" udpExecute: cl = %f < 0\n", CL(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[8].size > 1) {
        printf(" udpExecute: a should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (A(0) < 0) {
        printf(" udpExecute: a = %f < 0\n", A(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("thkcode(%d) = %s\n", numUdp, THKCODE(numUdp));
    printf("toc(%d)     = %f\n", numUdp, TOC(    numUdp));
    printf("xmaxt(%d)   = %f\n", numUdp, XMAXT(  numUdp));
    printf("leindex(%d) = %f\n", numUdp, LEINDEX(numUdp));
    printf("camcode(%d) = %s\n", numUdp, CAMCODE(numUdp));
    printf("cmax(%d)    = %f\n", numUdp, CMAX(   numUdp));
    printf("xmaxc(%d)   = %f\n", numUdp, XMAXC(  numUdp));
    printf("cl(%d)      = %f\n", numUdp, CL(     numUdp));
    printf("a(%d)       = %f\n", numUdp, A(      numUdp));
#endif

    /* mallocs required by Windows compiler */
    npnt = 128;
    MALLOC(xairfoil, double, npnt);
    MALLOC(yairfoil, double, npnt);

    /* call naca456 */
#ifdef WIN32
    NACA456 (&ithkcode, udps[numUdp].arg[1].val, udps[numUdp].arg[2].val,
                        udps[numUdp].arg[3].val,
             &icamcode, udps[numUdp].arg[5].val, udps[numUdp].arg[6].val,
                        udps[numUdp].arg[7].val, udps[numUdp].arg[8].val,
             &nairfoil, xairfoil, yairfoil);
#else
    naca456_(&ithkcode, udps[numUdp].arg[1].val, udps[numUdp].arg[2].val,
                        udps[numUdp].arg[3].val,
             &icamcode, udps[numUdp].arg[5].val, udps[numUdp].arg[6].val,
                        udps[numUdp].arg[7].val, udps[numUdp].arg[8].val,
             &nairfoil, xairfoil, yairfoil);
#endif

    /* always expecting at least 5 airfoil points */
    if (nairfoil < 5) {
        printf("naca456 returned nairfoil=%d\n", nairfoil);
        status = EGADS_NODATA;
        goto cleanup;
    } else if (nairfoil > npnt) {
        printf("naca456 returned nairfoil=%d but there is only room for %d\n", nairfoil, npnt);
        status = EGADS_NODATA;
        goto cleanup;
    }


#ifdef DEBUG
    printf("in udpExecute\n");
    for (ipnt = 0; ipnt < nairfoil; ipnt++) {
        printf("%4d  %12.6f %12.6f\n", ipnt, xairfoil[ipnt], yairfoil[ipnt]);
    }
#endif

    nle = (nairfoil - 1) / 2;

    /* mallocs required by Windows compiler */
    npnt = nairfoil;
    MALLOC(pts, double, 3*npnt);

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        pts[3*ipnt  ] = xairfoil[ipnt];
        pts[3*ipnt+1] = yairfoil[ipnt];
        pts[3*ipnt+2] = 0.0;
    }

    /* if just a camberline, create a WireBody of the meanline */
    if (fabs(TOC(numUdp)) < EPS06) {

        /* create Node at leading edge */
        ipnt = nle;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 data, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at trailing edge */
        ipnt = npnt - 1;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 data, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* create spline curve from LE to TE */
        sizes[0] = nle + 1;
        sizes[1] = 0;
        status = EG_approximate(context, 0, dxytol, sizes, &(pts[3*nle]), &ecurve);
        CHECK_STATUS(EG_approximate);

        /* make Edge for camberline */
        ipnt = nle;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_invEvaluate(ecurve, data, &(tdata[0]), result);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_getRange(ecurve, range, &periodic);
        CHECK_STATUS(EG_getRange);
        tdata[1] = range[1];

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* create loop of the Edges */
        sense[0] = SFORWARD;

        status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                 NULL, 1, eedges, sense, &eloop);
        CHECK_STATUS(EG_makeTopology);

        /* create the WireBody (which will be returned) */
        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

    /* create FaceBody of the airfoil */
    } else {

        /* determine if the airfoil has a sharp or blunt trailing edge */
        if (fabs(xairfoil[0]-xairfoil[nairfoil-1]) < EPS06 &&
            fabs(yairfoil[0]-yairfoil[nairfoil-1]) < EPS06   ) {
            bluntte = 0;
        } else {
            bluntte = 1;
        }

        /* create Node at upper trailing edge */
        ipnt = 0;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 data, 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at leading edge */
        ipnt = nle;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 data, 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* if a blunt trailing edge, create Node at lower trailing edge */
        if (bluntte) {
            ipnt = npnt - 1;
            data[0] = pts[3*ipnt  ];
            data[1] = pts[3*ipnt+1];
            data[2] = pts[3*ipnt+2];
            status = EG_makeTopology(context, NULL, NODE, 0,
                                     data, 0, NULL, NULL, &(enodes[2]));
            CHECK_STATUS(EG_makeTopology);

            enodes[3] = enodes[0];
        } else {
            enodes[2] = enodes[0];
        }

        /* create spline curve from upper TE, to LE, to lower TE */
        sizes[0] = npnt;
        sizes[1] = 0;
        status = EG_approximate(context, 0, dxytol, sizes, pts, &ecurve);
        CHECK_STATUS(EG_approximate);

        /* make Edge for upper surface */
        status = EG_getRange(ecurve, range, &periodic);
        CHECK_STATUS(EG_getRange);
        tdata[0] = range[0];

        ipnt = nle;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_invEvaluate(ecurve, data, &(tdata[1]), result);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* make Edge for lower surface */
        ipnt = nle;
        data[0] = pts[3*ipnt  ];
        data[1] = pts[3*ipnt+1];
        data[2] = pts[3*ipnt+2];
        status = EG_invEvaluate(ecurve, data, &(tdata[0]), result);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_getRange(ecurve, range, &periodic);
        CHECK_STATUS(EG_getRange);
        tdata[1] = range[1];

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        /* if a blunt trailing edge, create line segment at trailing edge */
        if (bluntte) {
            ipnt = npnt - 1;
            data[0] = pts[3*ipnt  ];
            data[1] = pts[3*ipnt+1];
            data[2] = pts[3*ipnt+2];
            data[3] = pts[0] - data[0];
            data[4] = pts[1] - data[1];
            data[5] = pts[2] - data[2];
            status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
            CHECK_STATUS(EG_makeGeometry);

            /* make Edge for this line */
            status = EG_invEvaluate(eline, data, &(tdata[0]), result);
            CHECK_STATUS(EG_invEvaluate);

            data[0] = pts[0];
            data[1] = pts[1];
            data[2] = pts[2];
            status = EG_invEvaluate(eline, data, &(tdata[1]), result);
            CHECK_STATUS(EG_invEvaluate);

            status = EG_makeTopology(context, eline, EDGE, TWONODE,
                                     tdata, 2, &(enodes[2]), NULL, &(eedges[2]));
            CHECK_STATUS(EG_makeTopology);
        }

        /* create loop of the edges Edges */
        sense[0] = SFORWARD;
        sense[1] = SFORWARD;
        sense[2] = SFORWARD;

        if (bluntte) {
            status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                     NULL, 3, eedges, sense, &eloop);
        } else {
            status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                     NULL, 2, eedges, sense, &eloop);
        }
        CHECK_STATUS(EG_makeTopology);

        /* make Face from the loop */
        status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
        CHECK_STATUS(EG_makeFace);

        /* find the direction of the Face normal */
        status = EG_getRange(eface, range, &periodic);
        CHECK_STATUS(EG_getRange);

        range[0] = (range[0] + range[1]) / 2;
        range[1] = (range[2] + range[3]) / 2;

        status = EG_evaluate(eface, range, eval);
        CHECK_STATUS(EG_evaluate);

        norm[0] = eval[4] * eval[8] - eval[5] * eval[7];
        norm[1] = eval[5] * eval[6] - eval[3] * eval[8];
        norm[2] = eval[3] * eval[7] - eval[4] * eval[6];

        /* if the normal is not positive, flip the Face */
        if (norm[2] < 0) {
            status = EG_flipObject(eface, &enew);
            CHECK_STATUS(EG_flipObject);
            eface = enew;
        }

        /* create the FaceBody (which will be returned) */
        status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                                 NULL, 1, &eface, sense, ebody);
        CHECK_STATUS(EG_makeTopology);
    }

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (pts      != NULL) free(pts     );
    if (xairfoil != NULL) free(xairfoil);
    if (yairfoil != NULL) free(yairfoil);

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
    int    iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
           (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
#endif

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
