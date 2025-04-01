/*
 ************************************************************************
 *                                                                      *
 * udfNaca6mc -- udp file to generate a NACA 6 series with multi-comber *
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

#define NUMUDPINPUTBODYS  1
#define NUMUDPARGS        2
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define CLT(IUDP,I)  ((double *) (udps[IUDP].arg[0].val))[I]
#define A(  IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"clt",    "a",      };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {12,       0,        };
static double argDdefs[NUMUDPARGS] = {0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#ifdef GRAFIC
   #include "grafic.h"
   #define DEBUG 1
#endif

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

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
udpExecute(ego  emodel,                 /* (in)  input model */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     npnt, ipnt, sense[3], sizes[2], periodic, i, j;
    int     oclass, mtype, nchild, *senses, nnode, sharpte;
    double  *pnt=NULL;
    double  clt, a, g, h, term1, term1p, term2, term2p, term3, term3p, term4, term4p, xc, yc;
    double  tlo, flo, tmid, fmid, thi, fhi;
    double  zeta, s, yt, ycp, theta, tle, result[18], data[18], eval[18];
    double  trange[4], range[4], norm[3];
    double  dxytol = EPS06;
    char    *message=NULL;
    ego     context, *ebodys, eref, enodes[4], eedges[3], ecurve, eline, eloop, eface, enew;
    ego     *echilds=NULL;
    udp_T   *udps = *Udps;

#ifdef GRAFIC
    float   xplot[5000], yplot[5000];
    int     io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
    int     nline=0, nplot=0, ilin[200], isym[200], nper[200];
    char    pltitl[80];
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("clt(  0) =");
    for (i = 0; i < udps[0].arg[0].size; i++) {
        printf(" %f", CLT(0,i));
    }
    printf("\n");
    printf("a(    0) =");
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf(" %f", A(  0,i));
    }
    printf("\n");
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size != udps[0].arg[1].size) {
        snprintf(message, 100, "clt and a should be the same length");
        status  = EGADS_RANGERR;
        goto cleanup;

    }

    /* check that Model that was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 100, "expecting Model to contain one Body (not %d)", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("clt(  %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[0].size; i++) {
        printf(" %f", CLT(numUdp,i));
    }
    printf("\n");
    printf("a(    %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf(" %f", A(  numUdp,i));
    }
    printf("\n");
#endif

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* get the number of Nodes in the input Body */
    status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, NULL);
    CHECK_STATUS(EG_getBodyTopos);

#ifdef DEBUG
    ocsmPrintEgo(ebodys[0]);
    printf("nnode=%d\n", nnode);
#endif

    /* if input Body does has only one Node, just generate a WireBody */
    if (nnode == 1) {

        /* mallocs required by Windows compiler */
        npnt = 101;
        MALLOC(pnt, double, (3*npnt    ));

        /* create points to represent the composite camber line */
        j = 0;
        pnt[3*j  ] = 0;
        pnt[3*j+1] = 0;
        pnt[3*j+2] = 0;

        for (j = 1; j < npnt-1; j++) {
            xc = (double)(j) / (double)(npnt-1);
            yc = 0;

#ifndef __clang_analyzer__
            for (i = 0; i < udps[0].arg[0].size; i++) {
                a   = A(  numUdp,i);
                clt = CLT(numUdp,i);

                g = -1/(1-a) * (a*a * (0.5 * log(a)-0.25) + 0.25);
                h =  1/(1-a) * (0.5 * (1-a)*(1-a) * log(1-a) - 0.25*(1-a)*(1-a)) + g;

                if (fabs(a-xc) < EPS20) {
                    term1 = 0;
                    term2 = 0;
                } else {
                    term1 = 0.50 * (a-xc) * (a-xc) * log(fabs(a-xc));
                    term2 = 0.25 * (a-xc) * (a-xc);
                }
                if (fabs(1-xc) < EPS20) {
                    term3 = 0;
                    term4 = 0;
                } else {
                    term3 = 0.50 * (1-xc) * (1-xc) * log(     1-xc);
                    term4 = 0.25 * (1-xc) * (1-xc);
                }
                yc += clt/(TWOPI*(a+1)) * (1/(1-a) * (term1 - term3 + term4 - term2) - xc*log(xc) + g - h*xc);
            }
#endif

            pnt[3*j  ] = xc;
            pnt[3*j+1] = yc;
            pnt[3*j+2] = 0;
        }

        j = npnt - 1;
        pnt[3*j  ] = 1;
        pnt[3*j+1] = 0;
        pnt[3*j+2] = 0;

#ifdef DEBUG
        for (j = 0; j < npnt; j++) {
            printf("%3d %10.5f %10.5f %10.5f\n", j, pnt[3*j], pnt[3*j+1], pnt[3*j+2]);
        }
#endif

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

     /* if input Body has multiple Nodes, use the first Edge to get the thickness */
    } else {

        /* mallocs required by Windows compiler */
        npnt = 201;
        MALLOC(pnt, double, 3*npnt);

        status = EG_getBodyTopos(ebodys[0], NULL, EDGE, &nchild, &echilds);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(echilds);

        status = EG_getRange(echilds[0], trange, &periodic);
        CHECK_STATUS(EG_getRange);

        /* points around airfoil (upper and then lower) */
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            zeta = TWOPI * ipnt / (npnt-1);
            s    = (1 + cos(zeta)) / 2;

            /* get the thickness (via a binary search) */
            if (s <= EPS06) {
                yt = 0;
            } else {
                tlo = trange[0];
                status = EG_evaluate(echilds[0], &tlo, data);
                CHECK_STATUS(EG_evaluate);
                flo = data[0] - s;

                if (fabs(flo) < EPS12) {
                    yt = data[1];

                } else {
                    thi = trange[1];
                    status = EG_evaluate(echilds[0], &thi, data);
                    CHECK_STATUS(EG_evaluate);
                    fhi = data[0] - s;

                    if (fabs(fhi) < EPS12) {
                        yt = data[1];

                    } else {
                        if (flo*fhi > 0) {
                            snprintf(message, 100, "%e and %e do not bracket root", flo, fhi);
                            status = EGADS_NOTFOUND;
                            goto cleanup;
                        }

                        while (thi-tlo >= EPS12) {
                            tmid = (tlo + thi) / 2;
                            status = EG_evaluate(echilds[0], &tmid, data);
                            CHECK_STATUS(EG_evaluate);
                            fmid = data[0] - s;

                            if (fabs(fmid) < EPS12) {
                                break;
                            } else if (flo*fmid <= 0) {
                                thi = tmid;
                            } else {
                                tlo = tmid;
                                flo = fmid;
                            }
                        }

                        status = EG_evaluate(echilds[0], &tmid, data);
                        CHECK_STATUS(EG_evaluate);

                        yt = data[1];
                    }
                }
            }

            /* create the camberline */
            yc  = 0;
            ycp = 0;
#ifndef __clang_analyzer__
            for (i = 0; i < udps[0].arg[0].size; i++) {
                a   = A(  numUdp,i);
                clt = CLT(numUdp,i);

                g = -1/(1-a) * (a*a * (0.5 * log(a)-0.25) + 0.25);
                h =  1/(1-a) * (0.5 * (1-a)*(1-a) * log(1-a) - 0.25*(1-a)*(1-a)) + g;

                if (fabs(a-s ) < EPS20) {
                    term1  = 0;
                    term1p = 0;
                    term2  = 0;
                    term2p = 0;
                } else if (a-s < 0) {
                    term1  = 0.50 * (s-a) * (s-a) * log(s-a);
                    term1p =        (s-a) * (log(s-a) + 0.5);
                    term2  = 0.25 * (s-a) * (s-a);
                    term2p = 0.50 * (a-s);
                } else {
                    term1  = 0.50 * (a-s) * (a-s) * log(a-s);
                    term1p =      - (a-s) * (log(a-s) + 0.5);
                    term2  = 0.25 * (a-s) * (a-s);
                    term2p = -.50 * (a-s);
                }
                if (fabs(1-s) < EPS20) {
                    term3  = 0;
                    term3p = 0;
                    term4  = 0;
                    term4p = 0;
                } else {
                    term3  = 0.50 * (1-s) * (1-s) * log(1-s );
                    term3p =      - (1-s) * (log(1-s) + 0.5);
                    term4  = 0.25 * (1-s) * (1-s);
                    term4p = -.50 * (1-s);
                }
                yc  += clt/(TWOPI*(a+1)) * ((term1  - term3  + term4  - term2 ) / (1-a) - s *log(MAX(s, EPS12))     + g - h*s);
                ycp += clt/(TWOPI*(a+1)) * ((term1p - term3p + term4p - term2p) / (1-a) -    log(MAX(s, EPS12)) - 1     - h  );
            }
#endif
            theta = atan(ycp);

#ifdef DEBUG
            printf("ipnt=%3d, s=%10.6f, yt=%10.6f, yc=%10.6f (%10.6f), theta=%10.6f\n", ipnt, s, yt, yc, ycp, theta);
#endif

            if (ipnt < npnt/2) {
                pnt[3*ipnt  ] = s  - yt * sin(theta);
                pnt[3*ipnt+1] = yc + yt * cos(theta);
                pnt[3*ipnt+2] = 0;

#ifdef GRAFIC
                xplot[nplot] = s;
                yplot[nplot] = yc;
                nplot++;

                xplot[nplot] = pnt[3*ipnt  ];
                yplot[nplot] = pnt[3*ipnt+1];
                nplot++;

                ilin[nline] = +GR_DASHED;
                isym[nline] = +GR_CIRCLE;
                nper[nline] = 2;
                nline++;
#endif
            } else if (ipnt == npnt/2) {
                pnt[3*ipnt  ] = 0;
                pnt[3*ipnt+1] = 0;
                pnt[3*ipnt+2] = 0;
            } else {
                pnt[3*ipnt  ] = s  + yt * sin(theta);
                pnt[3*ipnt+1] = yc - yt * cos(theta);
                pnt[3*ipnt+2] = 0;

#ifdef GRAFIC
                xplot[nplot] = s;
                yplot[nplot] = yc;
                nplot++;

                xplot[nplot] = pnt[3*ipnt  ];
                yplot[nplot] = pnt[3*ipnt+1];
                nplot++;

                ilin[nline] = +GR_DASHED;
                isym[nline] = +GR_CIRCLE;
                nper[nline] = 2;
                nline++;
#endif
            }
        }

#ifdef DEBUG
        for (j = 0; j < npnt; j++) {
            printf("%3d %10.5f %10.5f %10.5f\n", j, pnt[3*j], pnt[3*j+1], pnt[3*j+2]);
        }
#endif
#ifdef GRAFIC
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            xplot[nplot] = pnt[3*ipnt  ];
            yplot[nplot] = pnt[3*ipnt+1];
            nplot++;
        }
        ilin[nline] = GR_SOLID;
        isym[nline] = GR_PLUS;
        nper[nline] = npnt;
        nline++;

        snprintf(pltitl, 79, "~x~y~airfoil");
        grinit_(&io_kbd, &io_scr, "udpNaca6mc", strlen("udpNaca6mc"));
        grline_(ilin, isym, &nline, pltitl, &indgr, xplot, yplot, nper, strlen(pltitl));
#endif

        /* create Node at upper trailing edge */
        ipnt = 0;
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[0]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at leading edge */
        ipnt = (npnt - 1) / 2;
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 &(pnt[3*ipnt]), 0, NULL, NULL, &(enodes[1]));
        CHECK_STATUS(EG_makeTopology);

        /* create Node at lower trailing edge */
        if (fabs(pnt[0]-pnt[3*npnt-3]) < EPS06 &&
            fabs(pnt[1]-pnt[3*npnt-2]) < EPS06 &&
            fabs(pnt[2]-pnt[3*npnt-1]) < EPS06   ) {
            sharpte = 1;
        } else {
            sharpte = 0;
        }

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

        SPLINT_CHECK_FOR_NULL(*ebody);

        /* tell OpenCSM to put _body, _brch, and Branch Attributes on the Faces */
        status = EG_attributeAdd(*ebody, "__markFaces__", ATTRSTRING, 1, NULL, NULL, "true");
        CHECK_STATUS(EG_attributeAdd);
    }

#ifdef DEBUG
    printf("*ebody\n");
    ocsmPrintEgo(*ebody);
#endif

    /* set the output value(s) */

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    FREE(pnt);

    if (echilds != NULL) EG_free(echilds);

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
