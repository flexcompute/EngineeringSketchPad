/*
 ************************************************************************
 *                                                                      *
 * udpNaca -- udp file to generate a NACAmptt airfoil                   *
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

#define NUMUDPARGS 6
#include "udpUtilities.h"

#undef  OFFSET

/* shorthands for accessing argument values and velocities */
#define SERIES(       IUDP)  ((int    *) (udps[IUDP].arg[0].val))[0]
#define THICKNESS(    IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define THICKNESS_DOT(IUDP)  ((double *) (udps[IUDP].arg[1].dot))[0]
#define CAMBER(       IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define CAMBER_DOT(   IUDP)  ((double *) (udps[IUDP].arg[2].dot))[0]
#define MAXLOC(       IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
//$$$#define MAXLOC_DOT(   IUDP)  ((double *) (udps[IUDP].arg[3].dot))[0]
#define OFFSET(       IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define SHARPTE(      IUDP)  ((int    *) (udps[IUDP].arg[5].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"series", "thickness", "camber",    "maxloc", "offset", "sharpte", };
static int    argTypes[NUMUDPARGS] = {ATTRINT,  ATTRREALSEN, ATTRREALSEN, ATTRREAL, ATTRREAL, ATTRINT,   };
static int    argIdefs[NUMUDPARGS] = {12,       0,           0,           0,        0,        0,         };
static double argDdefs[NUMUDPARGS] = {0.,       0.,          0.,          0.40,     0.,       0.,        };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#ifdef GRAFIC
   #include "grafic.h"
   #define DEBUG
#endif

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

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

    int     npnt, ipnt, jpnt, kpnt, ile, sense[3], sizes[2], periodic, sharpte, shorten, i;
    int     add=1;
    double  m, p, t, zeta, s, yt, yc, theta, *pnt=NULL, *pnt_save=NULL;
    double  data[18], trange[2], tle, result[3], range[4], eval[18], norm[3];
    double  dx, dy, ds, x1, y1, x2, y2, x3, y3, x4, y4, dd, ss, tt, xx, yy, frac, dold;
    double  dxytol = 1.0e-6;
    char    *message=NULL;
    ego     enodes[4], eedges[3], ecurve, eline, eloop, eface, enew;

    ROUTINE(udpExecute);

#ifdef GRAFIC
    float   xplot[1203], yplot[1203];
    int     io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
    int     nline=0, nplot=0, ilin[3], isym[3], nper[3];
    char    pltitl[80];
#endif

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("series(        0) = %d\n", SERIES(       0));
    printf("thickness(     0) = %f\n", THICKNESS(    0));
    printf("thickness_dotT(0) = %f\n", THICKNESS_DOT(0));
    printf("camber(        0) = %f\n", CAMBER(       0));
    printf("camber_dot(    0) = %f\n", CAMBER_DOT(   0));
    printf("maxloc(        0) = %f\n", MAXLOC(       0));
//$$$    printf("maxloc_dot(    0) = %f\n", MAXLOC_DOT(   0));
    printf("offset(        0) = %f\n", OFFSET(       0));
    printf("sharpte(       0) = %d\n", SHARPTE(      0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "series should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (SERIES(0) <= 0) {
        snprintf(message, 100, "series = %d <= 0", SERIES(0));
        status  =  EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "thickness should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)THICKNESS(0) < 0) {
        snprintf(message, 100, "thickness = %f < 0", THICKNESS(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "camber should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        snprintf(message, 100, "maxloc should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)MAXLOC(0) <=0) {
        snprintf(message, 100, "maxloc = %f <= 0", MAXLOC(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if ((double)MAXLOC(0) >= 1) {
        snprintf(message, 100, "maxloc = %f >= 1", MAXLOC(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size > 1) {
        snprintf(message, 100, "offset should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        snprintf(message, 100, "sharpte should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (SHARPTE(0) != 0 && SHARPTE(0) != 1   ) {
        snprintf(message, 100, "sharpte should be 0 or 1");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("series(       %d) = %d\n", numUdp, SERIES(       numUdp));
    printf("thickness(    %d) = %f\n", numUdp, THICKNESS(    numUdp));
    printf("thickness_dot(%d) = %f\n", numUdp, THICKNESS_DOT(numUdp));
    printf("camber(       %d) = %f\n", numUdp, CAMBER(       numUdp));
    printf("camber_dot(   %d) = %f\n", numUdp, CAMBER_DOT(   numUdp));
    printf("maxloc(       %d) = %f\n", numUdp, MAXLOC(       numUdp));
//$$$    printf("maxloc_dot(   %d) = %f\n", numUdp, MAXLOC_DOT(   numUdp));
    printf("offset(       %d) = %f\n", numUdp, OFFSET(       numUdp));
    printf("sharpte(      %d) = %d\n", numUdp, SHARPTE(      numUdp));
#endif

    /* if thickness == camber == 0 and maxloc == 0.40, then use series */
    if (THICKNESS(0) == 0  && CAMBER(0) == 0 && fabs(MAXLOC(0)-0.40) < EPS06) {
        /* break series into camber (m), locn of max camber (p), and thickness (t) */
        m = (int)( SERIES(0)           / 1000);
        p = (int)((SERIES(0) - 1000*m) /  100);
        t =        SERIES(0)           %  100;

        if (p == 0) {
            p = 4;
        }

        m /= 100;
        p /=  10;
        t /= 100;
    } else {
        m = CAMBER(0);
        p = MAXLOC(0);
        t = THICKNESS(0);
    }

    /* if t > 0, generate FaceBody for profile */
    if (t > 0) {

        /* mallocs required by Windows compiler */
        npnt = 101;
        MALLOC(pnt,      double, 3*npnt);
        MALLOC(pnt_save, double, 3*npnt);

        sharpte = SHARPTE(0);
        ile     = (npnt - 1) / 2;

        /* points around airfoil (upper and then lower) */
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            zeta = TWOPI * ipnt / (npnt-1);
            s    = (1 + cos(zeta)) / 2;
            if (sharpte == 0) {
                yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
            } else {
                yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
            }

            if (s < p) {
                yc    =      m / (  p)/(  p) * (s * (2*p -   s));
                theta = atan(m / (  p)/(  p) * (    (2*p - 2*s)));
            } else {
                yc    =      m / (1-p)/(1-p) * ((1-2*p) + s * (2*p -   s));
                theta = atan(m / (1-p)/(1-p) * (              (2*p - 2*s)));
            }

            if (ipnt < npnt/2) {
                pnt[3*ipnt  ] = s  - yt * sin(theta);
                pnt[3*ipnt+1] = yc + yt * cos(theta);
                pnt[3*ipnt+2] = 0;
            } else if (ipnt == npnt/2) {
                pnt[3*ipnt  ] = 0;
                pnt[3*ipnt+1] = 0;
                pnt[3*ipnt+2] = 0;
            } else {
                pnt[3*ipnt  ] = s  + yt * sin(theta);
                pnt[3*ipnt+1] = yc - yt * cos(theta);
                pnt[3*ipnt+2] = 0;
            }
        }

#ifdef GRAFIC
        /* original points */
        if (npnt > 101) exit(0);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            xplot[nplot] = pnt[3*ipnt  ];
            yplot[nplot] = pnt[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = GR_PLUS;
        nper[nline] = npnt;
        nline++;
#endif

        /* create offset if required */
        if (OFFSET(0) != 0) {
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                pnt_save[3*ipnt  ] = pnt[3*ipnt  ];
                pnt_save[3*ipnt+1] = pnt[3*ipnt+1];
                pnt_save[3*ipnt+2] = pnt[3*ipnt+2];
            }

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                if (ipnt == 0) {
                    dx = pnt_save[3] - pnt_save[0];
                    dy = pnt_save[4] - pnt_save[1];
                } else if (ipnt == npnt-1) {
                    dx = pnt_save[3*ipnt  ] - pnt_save[3*ipnt-3];
                    dy = pnt_save[3*ipnt+1] - pnt_save[3*ipnt-2];
                } else {
                    dx = pnt_save[3*ipnt+3] - pnt_save[3*ipnt-3];
                    dy = pnt_save[3*ipnt+4] - pnt_save[3*ipnt-2];
                }

                ds = sqrt(dx * dx + dy * dy);
                pnt[3*ipnt  ] += OFFSET(0) * dy / ds;
                pnt[3*ipnt+1] -= OFFSET(0) * dx / ds;
            }
        }

        /* if offset is positive, extrapolate by offset distance by moving first and last points */
        if (OFFSET(0) > 0) {
            ipnt = 0;
            dold = sqrt(SQR(pnt[3*ipnt+3]-pnt[3*ipnt  ])
                       +SQR(pnt[3*ipnt+4]-pnt[3*ipnt+1])
                       +SQR(pnt[3*ipnt+5]-pnt[3*ipnt+2]));
            frac = (dold + OFFSET(0)) / dold;

            pnt[3*ipnt  ] = pnt[3*ipnt+3] + frac * (pnt[3*ipnt  ] - pnt[3*ipnt+3]);
            pnt[3*ipnt+1] = pnt[3*ipnt+4] + frac * (pnt[3*ipnt+1] - pnt[3*ipnt+4]);
            pnt[3*ipnt+2] = pnt[3*ipnt+5] + frac * (pnt[3*ipnt+2] - pnt[3*ipnt+5]);

            ipnt = npnt - 1;
            dold = sqrt(SQR(pnt[3*ipnt-3]-pnt[3*ipnt  ])
                       +SQR(pnt[3*ipnt-2]-pnt[3*ipnt+1])
                       +SQR(pnt[3*ipnt-1]-pnt[3*ipnt+2]));
            frac = (dold + OFFSET(0)) / dold;

            pnt[3*ipnt  ] = pnt[3*ipnt-3] - frac * (pnt[3*ipnt-3] - pnt[3*ipnt  ]);
            pnt[3*ipnt+1] = pnt[3*ipnt-2] - frac * (pnt[3*ipnt-2] - pnt[3*ipnt+1]);
            pnt[3*ipnt+2] = pnt[3*ipnt-1] - frac * (pnt[3*ipnt-1] - pnt[3*ipnt+2]);
        }

        /* if offset is negative, look for self-intersections near leading Edge */
        if (OFFSET(0) < 0) {
            shorten = 0;
            for (ipnt = npnt/4; ipnt < 3*npnt/4; ipnt++) {
                x1 = pnt[3*ipnt  ];
                y1 = pnt[3*ipnt+1];
                x2 = pnt[3*ipnt+3];
                y2 = pnt[3*ipnt+4];

                for (jpnt = ipnt+2; jpnt < 3*npnt/4; jpnt++) {
                    x3 = pnt[3*jpnt  ];
                    y3 = pnt[3*jpnt+1];
                    x4 = pnt[3*jpnt+3];
                    y4 = pnt[3*jpnt+4];

                    dd = (x2 - x1) * (y3 - y4) - (x3 - x4) * (y2 - y1);
                    if (fabs(dd) < EPS12) continue;

                    ss = ((x3 - x1) * (y3 - y4) - (x3 - x4) * (y3 - y1)) / dd;
                    tt = ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / dd;

                    if (ss >= 0 && ss <= 1 && tt >= 0 && tt <= 1) {
                        pnt[3*ipnt+3] = x1 * (1 - ss) + x2 * ss;
                        pnt[3*ipnt+4] = y1 * (1 - ss) + y2 * ss;
                        ipnt++;
                        ile = ipnt;

                        for (kpnt = jpnt+1; kpnt < npnt; kpnt++) {
                            pnt[3*ipnt+3] = pnt[3*kpnt  ];
                            pnt[3*ipnt+4] = pnt[3*kpnt+1];
                            pnt[3*ipnt+5] = pnt[3*kpnt+2];
                            ipnt++;
                        }
                        npnt = ipnt - 1;
                        shorten = 1;
                        break;
                    }
                }
                if (shorten == 1) break;
            }
        }

        /* if offset is negative, shorten to eliminate any self-intersection
           near trailing Edge */
        if (OFFSET(0) < 0) {
            shorten = 0;

            for (ipnt = 0; ipnt < npnt/4; ipnt++) {
                x1 = pnt[3*ipnt  ];
                y1 = pnt[3*ipnt+1];
                x2 = pnt[3*ipnt+3];
                y2 = pnt[3*ipnt+4];

                for (jpnt = npnt*3/4; jpnt < npnt-1; jpnt++) {
                    x3 = pnt[3*jpnt  ];
                    y3 = pnt[3*jpnt+1];
                    x4 = pnt[3*jpnt+3];
                    y4 = pnt[3*jpnt+4];

                    dd = (x2 - x1) * (y3 - y4) - (x3 - x4) * (y2 - y1);
                    if (fabs(dd) < EPS12) continue;

                    ss = ((x3 - x1) * (y3 - y4) - (x3 - x4) * (y3 - y1)) / dd;
                    tt = ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / dd;

                    if (ss >= 0 && ss <= 1 && tt >= 0 && tt <= 1) {
                        xx = (1 - ss) * x1 + ss * x2;
                        yy = (1 - ss) * y1 + ss * y2;

                        pnt[       0] = xx;
                        pnt[       1] = yy;

                        for (i = 1; i <= ipnt; i++) {
                            frac  = (double)(i) / (double)(ipnt+1);
                            pnt[3*i  ] = (1-frac) * xx + frac * pnt[3*ipnt+3];
                            pnt[3*i+1] = (1-frac) * yy + frac * pnt[3*ipnt+4];
                        }

                        pnt[3*npnt-3] = xx;
                        pnt[3*npnt-2] = yy;

                        for (i = 1; i < npnt-jpnt-1; i++) {
                            frac = (double)(i) / (double)(npnt-jpnt-1);
                            pnt[3*(npnt-i-1)  ] = (1-frac) * xx + frac * pnt[3*jpnt  ];
                            pnt[3*(npnt-i-1)+1] = (1-frac) * yy + frac * pnt[3*jpnt+1];
                        }

                        shorten = 1;
                        sharpte = 1;     /* offset curve will have sharp te */
                        break;
                    }
                }
                if (shorten) break;
            }

        /* if offset is positive, then we will not have a sharp te */
        } else if (OFFSET(0) > 0) {
            sharpte = 0;
        }

#ifdef GRAFIC
        /* offset points */
        if (OFFSET(0) != 0) {
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                xplot[nplot] = pnt[3*ipnt  ];
                yplot[nplot] = pnt[3*ipnt+1];
                nplot++;
            }
            ilin[nline] = -GR_DOTTED;
            isym[nline] = GR_CIRCLE;
            nper[nline] = npnt;
            nline++;
        }
#endif

        /* create Node at upper trailing edge */
        ipnt = 0;
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at leading edge */
        ipnt = ile;
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at lower trailing edge */
        if (sharpte == 0) {
            ipnt = npnt - 1;
            status = EG_makeTopology(context, NULL, NODE, 0,
                                     &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[2]));
            CHECK_STATUS(EG_makeTopology);

            enodes[3] = enodes[0];
        } else {
            enodes[2] = enodes[0];
        }

        /* create spline curve from upper TE, to LE, to lower TE */
        sizes[0] = npnt;
        sizes[1] = 0;
        status = EG_approximate(context, 0, dxytol, sizes, pnt, &ecurve);
        CHECK_STATUS(EG_approximate);

#ifdef GRAFIC
        /* curve fit */
        status = EG_getRange(ecurve, trange, &periodic);
        CHECK_STATUS(EG_getRange);

        for (ipnt = 0; ipnt < 1001; ipnt++) {
            double tt = trange[0] + (double)(ipnt) / (double)(1001-1) * (trange[1] - trange[0]);
            status = EG_evaluate(ecurve, &tt, data);
            CHECK_STATUS(EG_evaluate);

            xplot[nplot] = data[0];
            yplot[nplot] = data[1];
            nplot++;
        }
        ilin[nline] = GR_SOLID;
        isym[nline] = -GR_SQUARE;
        nper[nline] = 1001;
        nline++;

        sprintf(pltitl, "~x~y~airfoil: offset=%6.3f +=orig, O=offset", OFFSET(0));
        grinit_(&io_kbd, &io_scr, "udpNaca", strlen("udpNaca"));
        grline_(ilin, isym, &nline, pltitl, &indgr, xplot, yplot, nper, strlen(pltitl));
#endif

        /* find t at leading edge point */
        ipnt = (npnt - 1) / 2;
        status = EG_invEvaluate(ecurve, &(pnt[3*ipnt]), &(tle), result);
        CHECK_STATUS(EG_invEvaluate);

        /* make Edge for upper surface */
        status = EG_getRange(ecurve, trange, &periodic);
        CHECK_STATUS(EG_getRange);
        trange[1] = tle;

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 trange, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        /* make Edge for lower surface */
        status = EG_getRange(ecurve, trange, &periodic);
        CHECK_STATUS(EG_getRange);
        trange[0] = tle;

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 trange, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        /* create line segment at trailing edge */
        if (sharpte == 0) {
            ipnt = npnt - 1;
            data[0] = pnt[3*ipnt  ];
            data[1] = pnt[3*ipnt+1];
            data[2] = pnt[3*ipnt+2];
            data[3] = pnt[0] - data[0];
            data[4] = pnt[1] - data[1];
            data[5] = pnt[2] - data[2];
            status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
            CHECK_STATUS(EG_makeGeometry);

            /* make Edge for this line */
            status = EG_invEvaluate(eline, data, &(trange[0]), result);
            CHECK_STATUS(EG_invEvaluate);

            status = EG_invEvaluate(eline, &(pnt[0]), &(trange[1]), result);
            CHECK_STATUS(EG_invEvaluate);

            status = EG_makeTopology(context, eline, EDGE, TWONODE,
                                     trange, 2, &(enodes[2]), NULL, &(eedges[2]));
            CHECK_STATUS(EG_makeTopology);
        }

        /* create loop of the two or three Edges */
        sense[0] = SFORWARD;
        sense[1] = SFORWARD;
        sense[2] = SFORWARD;

        if (sharpte == 0) {
            status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                     NULL, 3, eedges, sense, &eloop);
            CHECK_STATUS(EG_makeTopology);
        } else {
            status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                     NULL, 2, eedges, sense, &eloop);
            CHECK_STATUS(EG_makeTopology);
        }

        /* make Face from the loop */
        status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
        CHECK_STATUS(EG_makeFace);

        /* since this will make a PLANE, we need to add an Attribute
           to tell OpenCSM to scale the UVs when computing sensitivities */
        status = EG_attributeAdd(eface, "_scaleuv", ATTRINT, 1,
                                 &add, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

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

    /* t == 0, so create WireBody of the meanline */
    } else {

        /* mallocs required by Windows compiler */
        npnt = 51;
        MALLOC(pnt, double, (3*npnt    ));

        /* points along meanline (leading edge to trailing edge) */
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            zeta = PI * ipnt / (npnt-1);
            s    = (1 - cos(zeta)) / 2;

            if (s < p) {
                yc =  m / (  p)/(  p) * (          s * (2*p - s));
            } else {
                yc =  m / (1-p)/(1-p) * ((1-2*p) + s * (2*p - s));
            }

            pnt[3*ipnt  ] = s;
            pnt[3*ipnt+1] = yc;
            pnt[3*ipnt+2] = 0;
        }

        /* create Node at leading edge */
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(pnt[0]), 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at trailing edge */
        ipnt = npnt - 1;
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* create spline curve from LE to TE */
        sizes[0] = npnt;
        sizes[1] = 0;
        status = EG_approximate(context, 0, dxytol, sizes, pnt, &ecurve);
        CHECK_STATUS(EG_approximate);

        /* make Edge for camberline */
        status = EG_getRange(ecurve, trange, &periodic);
        CHECK_STATUS(EG_getRange);

        status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                 trange, 2, &(enodes[0]), NULL, &(eedges[0]));
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
    }

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    FREE(pnt     );
    FREE(pnt_save);

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
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp, ipnt, nnode, nedge, nface, nchild, oclass, mtype, *senses;
    int    sharpte, iter;
    double x, dx, dx_ds, dx_dt, y, dy, dy_ds, dy_dt, dyt_ds, dyt_dt, dyc_ds, th, dth_ds, D;
    double s, t, t_dot, m, m_dot, p;
    double yt, yt_dot, yc, yc_dot, theta, theta_dot;
    double data[18], temp1, temp2;
    ego    eref, *echilds, *enodes, *eedges, *efaces, eent;

    ROUTINE(udpSensitivity);

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

    /* get parameters and their derivatives */
    if (THICKNESS(iudp) == 0  && CAMBER(iudp) == 0 && fabs(MAXLOC(iudp)-0.40) < EPS06) {
        /* break series into camber (m), locn of max camber (p), and thickness (t) */
        m = (int)( SERIES(iudp)           / 1000);
        p = (int)((SERIES(iudp) - 1000*m) /  100);
        t =        SERIES(iudp)           %  100;

        if (p == 0) {
            p = 4;
        }

        t      /= 100;
        t_dot   = 0;
        m      /= 100;
        m_dot   = 0;
        p      /=  10;
        sharpte = SHARPTE(      iudp);
    } else {
        t       = THICKNESS(    iudp);
        t_dot   = THICKNESS_DOT(iudp);
        m       = CAMBER(       iudp);
        m_dot   = CAMBER_DOT(   iudp);
        p       = MAXLOC(       iudp);
        sharpte = SHARPTE(      iudp);
    }

    /* find the ego entity */
    if (entType == OCSM_NODE) {
        status = EG_getBodyTopos(ebody, NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        eent = enodes[entIndex-1];

        EG_free(enodes);
    } else if (entType == OCSM_EDGE) {
        status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, &eedges);
        CHECK_STATUS(EG_getBodyTopos);

        eent = eedges[entIndex-1];

        EG_free(eedges);
    } else if (entType == OCSM_FACE) {
        status = EG_getBodyTopos(ebody, NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        eent = efaces[entIndex-1];

        EG_free(efaces);
    } else {
        printf("udpSensitivity: bad entType=%d\n", entType);
        status = EGADS_ATTRERR;
        goto cleanup;
    }

    /* loop through the points */
    for (ipnt = 0; ipnt < npnt; ipnt++) {

        /* find the physical coordinates */
        if        (entType == OCSM_NODE) {
            status = EG_getTopology(eent, &eref, &oclass, &mtype,
                                    data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);
        } else if (entType == OCSM_EDGE) {
            status = EG_evaluate(eent, &(uvs[ipnt]), data);
            CHECK_STATUS(EG_evaluate);
        } else if (entType == OCSM_FACE) {
            status = EG_evaluate(eent, &(uvs[2*ipnt]), data);
            CHECK_STATUS(EG_evaluate);
        }

        /* special case for leading edge */
        if (fabs(data[0]) < EPS06 && fabs(data[1]) < EPS06) {
            vels[3*ipnt  ] = 0;
            vels[3*ipnt+1] = 0;
            vels[3*ipnt+2] = 0;
            continue;
        }

        /* find the s and thick associated with data[] via 
           a Newton search */
        s = fabs(data[0]);
        if        (s < EPS06) {
            s = EPS06;
        } else if (s > 1) {
            s = 1;
        }

        for (iter = 0; iter < 30; iter++) {

            if (sharpte == 0) {
                yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
            } else {
                yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
            }

            if (s < p) {
                yc =      m / (  p)/(  p) * (s * (2*p -   s));
                th = atan(m / (  p)/(  p) * (    (2*p - 2*s)));
            } else {
                yc =      m / (1-p)/(1-p) * ((1-2*p) + s * (2*p -   s));
                th = atan(m / (1-p)/(1-p) * (              (2*p - 2*s)));
            }

            x = s  - yt * sin(th);
            y = yc + yt * cos(th);

            dx = x - data[0];
            dy = y - data[1];

            if (fabs(dx) < EPS06 && fabs(dy) < EPS12) break;

            if (sharpte == 0) {
                dyt_ds =   t/0.20 * (0.2969/2 / sqrt(s) +     (-0.1260 + s * (-0.3516*2 + s * ( 0.2843*3 + s * (-0.1015*4)))));
                dyt_dt = 1.0/0.20 * (0.2969   * sqrt(s) + s * (-0.1260 + s * (-0.3516   + s * ( 0.2843   + s * (-0.1015  )))));
            } else {
                dyt_ds =   t/0.20 * (0.2969/2 / sqrt(s) +     (-0.1260 + s * (-0.3516*2 + s * ( 0.2843*3 + s * (-0.1036*4)))));
                dyt_dt = 1.0/0.20 * (0.2969   * sqrt(s) + s * (-0.1260 + s * (-0.3516   + s * ( 0.2843   + s * (-0.1036  )))));
            }

            if (s < p) {
                dyc_ds = m / (  p)/(  p) * (2*p - 2*s);
                dth_ds = 1 / (1 + dyc_ds * dyc_ds);
            } else {
                dyc_ds = m / (1-p)/(1-p) * (2*p - 2*s);
                dth_ds = 1 / (1 + dyc_ds * dyc_ds);
            }

            dx_ds = 1      - sin(th) * dyt_ds - yt * cos(th) * dth_ds;
            dx_dt =        - sin(th) * dyt_dt;
            dy_ds = dyc_ds + cos(th) * dyt_ds - yt * sin(th) * dth_ds;
            dy_dt =        + cos(th) * dyt_dt;
            
            D   =  dx_ds * dy_dt - dy_ds * dx_dt;
            s  -= (dx    * dy_dt - dy    * dx_dt) / D;
            t  -= (dx_ds * dy    - dy_ds * dx   ) / D;

            if        (s < EPS06) {
                s = EPS06;
            } else if (s > 1) {
                s = 1;
            }
        }

        /* evaluate at this s and t */
        if (sharpte == 0) {
            yt     = t     / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
            yt_dot = t_dot / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
        } else {
            yt     = t     / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
            yt_dot = t_dot / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
        }

        if (t < 0) yt_dot = -yt_dot;

        if (s < p) {
            temp1     = (s * (2*p -   s)) / (p) / (p);
            temp2     = (    (2*p - 2*s)) / (p) / (p);
            theta     = atan(m * temp2);
            yc_dot    = m_dot * temp1;
            theta_dot = m_dot * temp2 / (1 + m*m*temp2*temp2);
        } else {
            temp1     = ((1-2*p) + s * (2*p -   s)) / (1-p) / (1-p);
            temp2     = (              (2*p - 2*s)) / (1-p) / (1-p);
            theta     = atan(m * temp2);
            yc_dot    = m_dot * temp1;
            theta_dot = m_dot * temp2 / (1 + m*m*temp2*temp2);
        }

        vels[3*ipnt  ] =        - yt_dot * sin(theta) - theta_dot * yt * cos(theta);
        vels[3*ipnt+1] = yc_dot + yt_dot * cos(theta) - theta_dot * yt * sin(theta);
        vels[3*ipnt+2] = 0;
    }

cleanup:
    return status;
}
