/*
 ************************************************************************
 *                                                                      *
 * udpKulfan -- udp file to generate Kulfan airfoils                    *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 *            Based upon prescription by B.M. Kulfan and J.E. Bussoletti*
 *              "Fundamental Parametric Geometry Representations for    *
 *               Aircraft Component Shapes", AIAA-2006-6948             *
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

#define NUMUDPARGS 5
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define CLASS( IUDP,I)  ((double *) (udps[IUDP].arg[0].val))[I]
#define ZTAIL( IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]
#define AUPPER(IUDP,I)  ((double *) (udps[IUDP].arg[2].val))[I]
#define ALOWER(IUDP,I)  ((double *) (udps[IUDP].arg[3].val))[I]
#define NUMPTS(IUDP  )  ((int    *) (udps[IUDP].arg[4].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"class",  "ztail",  "aupper", "alower", "numpts", };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, ATTRINT,  };
static int    argIdefs[NUMUDPARGS] = {0,        0,        0,        0,        101,      };
static double argDdefs[NUMUDPARGS] = {0.,       0.,       0.,       0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
 udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define           HUGEQ           99999999.0
#define           PIo2            1.5707963267948965579989817
#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

#ifdef GRAFIC
    #include "grafic.h"
#endif

/*
 ************************************************************************
 *                                                                      *
 *   factorial - factorial of a number (with cache)                     *
 *                                                                      *
 ************************************************************************
 */

static long
factorial(int n)
{
    long result;
    static long table[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                           362880, 3628800, 39916800, 479001600};

    /* bad input */
    if        (n < 0) {
        result = 0;

    /* use (cached) result from table */
    } else if (n < 13) {
        result = table[n];

    /* compute the result if not in the table */
    } else {
        result = n;
        while (--n > 1) {
            result *= n;
        }
    }

    return result;
}


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

    int     r, n, ipnt, sense[3], sizes[2], periodic, nedge;
    double  zeta, shape, s, K, pow1, pow2, *pnts=NULL;
    double  data[18], tdata[2], result[3], range[4], eval[18], norm[3];
    double  dxytol = 1.0e-6;
    ego     enodes[4], eedges[3], ecurve, eline, eloop, eface, enew;
#ifdef GRAFIC
    float   *xplot=NULL, *yplot=NULL;
#endif

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if        (udps[0].arg[0].size < 2) {
        printf(" udpExecute: class should contain 2 values (nose,tail)\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size < 2) {
        printf(" udpExecute: ztail should contain 2 values (upper,lower)\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size < 2) {
        printf(" udpExecute: aupper should contain at least 1 value\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size < 2) {
        printf(" udpExecute: aupper should contain at least 1 value\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size != 1 || NUMPTS(0) < 11) {
        printf(" udpExecute: numpts should contain one positive number greater than 11\n");
        status = EGADS_RANGERR;
        goto cleanup;

    }

#ifdef DEBUG
    {
        int i;

        printf("class( 0)[0]   = %f  (nose)\n", CLASS(0,0));
        printf("class( 0)[1]   = %f  (tail)\n", CLASS(0,1));
        printf("ztail( 0)[0]   = %f  (upper)\n", ZTAIL(0,0));
        printf("ztail( 0)[1]   = %f  (lower)\n", ZTAIL(0,1));
        for (i = 0; i < udps[0].arg[2].size; i++) {
            printf("aupper(0)[%d]  = %f\n", i, AUPPER(0,i));
        }
        for (i = 0; i < udps[0].arg[3].size; i++) {
            printf("alower(0)[%d]  = %f\n", i, ALOWER(0,i));
        }
        printf("numpts(0)      = %d\n", NUMPTS( 0));
    }
#endif

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

    /* mallocs required by Windows compiler */
    MALLOC(pnts, double, (3*NUMPTS(numUdp)));

    /* points around airfoil ( upper and lower) */
    for (ipnt = 0; ipnt < NUMPTS(numUdp); ipnt++) {
        zeta = TWOPI * ipnt / (NUMPTS(numUdp)-1);
        s    = (1 + cos(zeta)) / 2;

        /* upper surface */
        if ( ipnt < (NUMPTS(numUdp)-1)/2){
            n = udps[numUdp].arg[2].size - 1;

            shape = 0;
            for (r = 0; r <= n; r++) {
                pow1 = pow(1-s, n-r);
                pow2 = pow(  s,   r);
                K    = factorial(n) / factorial(r) / factorial(n-r);
                shape += AUPPER(numUdp,r) * K * pow1 * pow2;
            }

            pow1 = pow(  s, CLASS(numUdp,0));
            pow2 = pow(1-s, CLASS(numUdp,1));
            pnts[3*ipnt  ] = s;
            pnts[3*ipnt+1] = pow1 * pow2 * shape + ZTAIL(numUdp,0) * s;
            pnts[3*ipnt+2] = 0;

        /* leading edge */
        } else if (ipnt == (NUMPTS(numUdp)-1)/2) {
            pnts[3*ipnt  ] = 0;
            pnts[3*ipnt+1] = 0;
            pnts[3*ipnt+2] = 0;

        /* lower surface */
        } else if (ipnt > (NUMPTS(numUdp)-1)/2) {
            n = udps[numUdp].arg[3].size - 1;

            shape = 0;
            for (r = 0; r <= n; r++) {
                pow1 = pow(1-s, n-r);
                pow2 = pow(  s,   r);
                K = factorial(n) / factorial(r) / factorial(n-r);
                shape += ALOWER(numUdp,r) * K * pow1 * pow2;
            }

            pow1 = pow(  s, CLASS(numUdp,0));
            pow2 = pow(1-s, CLASS(numUdp,1));
            pnts[3*ipnt  ] = s;
            pnts[3*ipnt+1] = pow1 * pow2 * shape + ZTAIL(numUdp,1) * s;
            pnts[3*ipnt+2] = 0;
        }
    }

    /* create Node at upper trailing edge */
    ipnt = 0;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    /* create Node at leading edge */
    ipnt = (NUMPTS(numUdp) - 1) / 2;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create Node at lower trailing edge */
    ipnt = NUMPTS(numUdp) - 1;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &(enodes[2]));
    CHECK_STATUS(EG_makeTopology);

    enodes[3] = enodes[0];

#ifdef DEBUG
    {
        int i;

        for (i = 0; i < NUMPTS(numUdp); i++) {
            printf("%3d  %10.4f %10.4f %10.4f\n", i, pnts[3*i], pnts[3*i+1], pnts[3*i+2]);
        }
    }
#endif

    /* create spline curve from upper TE, to LE, to lower TE */
    sizes[0] = NUMPTS(numUdp);
    sizes[1] = 0;
    status = EG_approximate(context, 0, dxytol, sizes, pnts, &ecurve);
    CHECK_STATUS(EG_approximate);

#ifdef GRAFIC
    #define NUMEVAL 5000

    if (1) {
        int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        int    nline=0, npnt=0, ilin[NUMEVAL+3], isym[NUMEVAL+3], nper[NUMEVAL+3], i;
        double tt, evaldata[18];

        MALLOC(xplot, float, NUMPTS(numUdp)+3*NUMEVAL);
        MALLOC(yplot, float, NUMPTS(numUdp)+3*NUMEVAL);

        for (i = 0; i < NUMPTS(numUdp); i++) {
            xplot[npnt] = pnts[3*i  ];
            yplot[npnt] = pnts[3*i+1];
            npnt++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] =  NUMPTS(numUdp);
        nline++;

        for (i = 0; i < NUMEVAL; i++) {
            tt = (double)(i) / (double)(NUMEVAL-1);
            (void) EG_evaluate(ecurve, &tt, evaldata);

            xplot[npnt] = evaldata[0];
            yplot[npnt] = evaldata[1];
            npnt++;
        }
        ilin[nline] = +GR_SOLID;
        isym[nline] = -GR_PLUS;
        nper[nline] =  NUMEVAL;
        nline++;

        for (i = 0; i < NUMEVAL; i++) {
            tt = (double)(i) / (double)(NUMEVAL+1);
            (void) EG_evaluate(ecurve, &tt, evaldata);

            xplot[npnt] = evaldata[0];
            yplot[npnt] = evaldata[1];
            npnt++;

            xplot[npnt] = evaldata[0] + evaldata[4] / 100;
            yplot[npnt] = evaldata[1] - evaldata[3] / 100;
            npnt++;

            ilin[nline] = GR_SOLID;
            isym[nline] = 0;
            nper[nline] = 2;
            nline++;
        }

        grinit_(&io_kbd, &io_scr, "udpKulfan", STRLEN("udpKulfan"));
        grline_(ilin, isym, &nline,                 "~x~y~O=points, line=fit",
                &indgr, xplot, yplot, nper, STRLEN( "~x~y~O=points, line=fit"));
    }
#endif

    /* make Edge for upper surface */
    tdata[0] = 0;   // because EG_approximate was used

    ipnt = (NUMPTS(numUdp) - 1) / 2;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_invEvaluate(ecurve, data, &(tdata[1]), result);
    CHECK_STATUS(EG_invEvaluate);

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             tdata, 2, &(enodes[0]), NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

    /* make Edge for lower surface */
    ipnt = (NUMPTS(numUdp) - 1) / 2;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_invEvaluate(ecurve, data, &(tdata[0]), result);
    CHECK_STATUS(EG_invEvaluate);

    tdata[1] = 1;   // because EG_approximate was used

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create line segment at trailing edge */
    ipnt = NUMPTS(numUdp) - 1;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    data[3] = pnts[0] - data[0];
    data[4] = pnts[1] - data[1];
    data[5] = pnts[2] - data[2];

    if (fabs(data[3]) > EPS06 || fabs(data[4]) > EPS06 || fabs(data[5]) > EPS06) {
        nedge = 3;

        status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
        CHECK_STATUS(EG_makeGeometry);

        /* make Edge for this line */
        status = EG_invEvaluate(eline, data, &(tdata[0]), result);
        CHECK_STATUS(EG_invEvaluate);

        data[0] = pnts[0];
        data[1] = pnts[1];
        data[2] = pnts[2];
        status = EG_invEvaluate(eline, data, &(tdata[1]), result);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeTopology(context, eline, EDGE, TWONODE,
                                 tdata, 2, &(enodes[2]), NULL, &(eedges[2]));
        CHECK_STATUS(EG_makeTopology);
    } else {
        nedge = 2;

        status = EG_deleteObject(eedges[1]);
        CHECK_STATUS(EG_deleteObject);

        status = EG_deleteObject(enodes[2]);
        CHECK_STATUS(EG_deleteObject);

        enodes[2] = enodes[0];

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 tdata, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);
    }

    /* create loop of either two or three Edges */
    sense[0] = SFORWARD;
    sense[1] = SFORWARD;
    sense[2] = SFORWARD;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                             NULL, nedge, eedges, sense, &eloop);
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

    /* set the output value(s) */

    udps[numUdp].ebody = *ebody;
#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (pnts != NULL) free(pnts);

#ifdef GRAFIC
    FREE(xplot);
    FREE(yplot);
#endif

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
//$$$    int    status = EGADS_SUCCESS;

    int    iudp, judp;

//$$$    int    ipnt, nedge, r, n;
//$$$    double xyz[18], shape, shape_dot, s, K, temp;
//$$$    ego    *eedges=NULL;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

//$$$#ifdef DEBUG
//$$$    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
//$$$           (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
//$$$#endif

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

    /* velocities for Nodes */
//$$$    if (entType == OCSM_NODE) {
//$$$        if        (entIndex == 1) {
//$$$            vels[0] = 0;
//$$$            vels[1] = ZTAIL_DOT(iudp,0);
//$$$            vels[2] = 0;
//$$$        } else if (entIndex == 2) {
//$$$            vels[0] = 0;
//$$$            vels[1] = 0;
//$$$            vels[2] = 0;
//$$$        } else if (entIndex == 3) {
//$$$            vels[0] = 0;
//$$$            vels[1] = ZTAIL_DOT(iudp,1);
//$$$            vels[2] = 0;
//$$$        } else {
//$$$            printf(" udpSensitivity: bad node index\n");
//$$$            return EGADS_GEOMERR;
//$$$        }
//$$$
    /* velocities for Edges */
//$$$    } else if (entType == OCSM_EDGE) {
//$$$        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
//$$$        if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$        if        (entIndex == 1) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedges[0], &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                n = udps[iudp].arg[2].size - 1;
//$$$                s = xyz[0];
//$$$
//$$$                shape     = 0;
//$$$                shape_dot = 0;
//$$$                for (r = 0; r <= n; r++) {
//$$$                    K = factorial(n) / factorial(r) / factorial(n-r);
//$$$                    shape     += AUPPER(    iudp,r) * K * pow(1-s,n-r) * pow(s,r);
//$$$                    shape_dot += AUPPER_DOT(iudp,r) * K * pow(1-s,n-r) * pow(s,r);
//$$$                }
//$$$                temp = pow(s,CLASS(iudp,0)) * pow(1-s,CLASS(iudp,1));
//$$$
//$$$                vels[3*ipnt  ] = 0;
//$$$                vels[3*ipnt+1] = temp * shape_dot + ZTAIL_DOT(iudp,0) * s;
//$$$                if (s   > 0) vels[3*ipnt+1] += temp * shape * log(  s) * CLASS_DOT(iudp,0);
//$$$                if (1-s > 0) vels[3*ipnt+1] += temp * shape * log(1-s) * CLASS_DOT(iudp,1);
//$$$                vels[3*ipnt+2] = 0;
//$$$            }
//$$$        } else if (entIndex == 2) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedges[1], &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                n = udps[iudp].arg[3].size - 1;
//$$$                s = xyz[0];
//$$$
//$$$                shape     = 0;
//$$$                shape_dot = 0;
//$$$                for (r = 0; r <= n; r++) {
//$$$                    K = factorial(n) / factorial(r) / factorial(n-r);
//$$$                    shape     += ALOWER(    iudp,r) * K * pow(1-s,n-r) * pow(s,r);
//$$$                    shape_dot += ALOWER_DOT(iudp,r) * K * pow(1-s,n-r) * pow(s,r);
//$$$                }
//$$$                temp = pow(s,CLASS(iudp,0)) * pow(1-s,CLASS(iudp,1));
//$$$
//$$$                vels[3*ipnt  ] = 0;
//$$$                vels[3*ipnt+1] = temp * shape_dot + ZTAIL_DOT(iudp,1) * s;
//$$$                if (s   > 0) vels[3*ipnt+1] += temp * shape * log(  s) * CLASS_DOT(iudp,0);
//$$$                if (1-s > 0) vels[3*ipnt+1] += temp * shape * log(1-s) * CLASS_DOT(iudp,1);
//$$$                vels[3*ipnt+2] = 0;
//$$$            }
//$$$        } else if (entIndex == 3) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                status = EG_evaluate(eedges[2], &(uvs[ipnt]), xyz);
//$$$                if (status != EGADS_SUCCESS) goto cleanup;
//$$$
//$$$                s = (xyz[1] - ZTAIL(iudp,1)) / (ZTAIL(iudp,0) - ZTAIL(iudp,1));
//$$$
//$$$                vels[3*ipnt  ] = 0;
//$$$                vels[3*ipnt+1] = (1-s) * ZTAIL_DOT(iudp,1) + s * ZTAIL_DOT(iudp,0);
//$$$                vels[3*ipnt+2] = 0;
//$$$            }
//$$$        } else {
//$$$            printf(" udpSensitivity: bad index\n");
//$$$            return EGADS_GEOMERR;
//$$$        }
//$$$
    /* velocities for Face */
//$$$    } else if (entType == OCSM_FACE) {
//$$$        if (entIndex == 1) {
//$$$            for (ipnt = 0; ipnt < npnt; ipnt++) {
//$$$                vels[3*ipnt  ] = 0;
//$$$                vels[3*ipnt+1] = 0;
//$$$                vels[3*ipnt+2] = 0;
//$$$            }
//$$$        } else {
//$$$            printf(" udpSensitivity: bad Face index\n");
//$$$            return EGADS_GEOMERR;
//$$$        }
//$$$
//$$$    } else {
//$$$        printf(" udpSensitivity: bad entity type\n");
//$$$        return EGADS_GEOMERR;
//$$$    }
//$$$
//$$$    status = EGADS_SUCCESS;
//$$$
//$$$cleanup:
//$$$    if (eedges != NULL) EG_free(eedges);
//$$$
//$$$    return status;
}
