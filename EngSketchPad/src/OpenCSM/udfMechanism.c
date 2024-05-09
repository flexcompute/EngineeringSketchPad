/*
 ************************************************************************
 *                                                                      *
 * udfMchanism -- solve mechanism equations and move Bodys              *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Wil Martin                *
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

#define NUMUDPARGS 1
#define NUMUDPINPUTBODYS -99
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FIXED(IUDP) ((char   *) udps[IUDP].arg[0].val)

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"fixed",     };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING,  };
static int    argIdefs[NUMUDPARGS] = {0,           };
static double argDdefs[NUMUDPARGS] = {0.,          };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"

typedef struct {
    char   name[10];          /* point name */
    int    stat;              /* =0 not fixed, =1 fixed */
    double x;                 /* x-location (when fixed) */
    double y;                 /* y-location (when fixed) */
} pnt_T;

typedef struct {
    int    npnt;              /* number of points */
    int    pnt[10];           /* point indicies */
    double x[10];             /* x-locations */
    double y[10];             /* y-locations */
    int    nfix;              /* number of fixed points */
    double dx;                /* x-translation */
    double dy;                /* y-translation */
    double xrot;              /* x-rotation center */
    double yrot;              /* y-rotation center */
    double ang;               /* rotation angle */
} bar_T;

#ifdef DEBUG
static void printall(int npnt, pnt_T pnt[], int nbar, bar_T bar[]);
#endif


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

    int     oclass, mtype, ichild, nchild, *senses, nchange;
    int     iattr, nattr, attrType, attrLen;
    int     iter, ipnt, jpnt, kpnt, ibar, ibar1, ipnt1, jpnt1, ibar2, ipnt2, jpnt2;
    CINT    *tempIlist;
    double  data[4], L0, L1, L2, phi0, phi1, phi2, alfa, beta, theta1, theta2;
    double  xavg, yavg, theta1a, xa, ya, dista, theta1b, xb, yb, distb;
    double  cosang, sinang, xform11, xform12, xform13, xform21, xform22, xform23, xold, yold, M[12];
    CDOUBLE *tempRlist;
    char    *message=NULL, *fixed=NULL, *name=NULL;
    CCHAR   *tempClist, *attrName;
    ego     context, eref, *ebodys, exform, *newBodys=NULL;
#ifdef DEBUG
    udp_T   *udps = *Udps;
#else
    udp_T   *udps;
#endif

    int     npnt;
    pnt_T   *pnt = NULL;
    int     nbar;
    bar_T   *bar = NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("fixed(0) = %s\n", FIXED(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains more than one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, " udpExecute: expecting a Model");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild <= 1) {
        snprintf(message, 100, " udpExecute: expecting Model to contain multiple Bodys (not %d)", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("fixed(%d) = %s\n", numUdp, FIXED(numUdp));
#endif

    MALLOC(pnt, pnt_T, 3*nchild);
    MALLOC(bar, bar_T,   nchild);

    for (ipnt = 0; ipnt < 3*nchild; ipnt++) {
        pnt[ipnt].stat = 0;
        pnt[ipnt].x    = 0;
        pnt[ipnt].y    = 0;
    }
    for (ibar = 0; ibar < nchild; ibar++) {
        bar[ibar].npnt = 0;
        bar[ibar].nfix = 0;
        bar[ibar].dx   = 0;
        bar[ibar].dy   = 0;
        bar[ibar].xrot = 0;
        bar[ibar].yrot = 0;
        bar[ibar].ang  = 0;
    }

    npnt = 0;
    nbar = 0;

    /* loop through input Bodys and set up pnt and bar entries */
    for (ichild = 0; ichild < nchild; ichild++) {
        status = EG_attributeNum(ebodys[ichild], &nattr);
        CHECK_STATUS(EG_attributeNum);

        for (iattr = 1; iattr <= nattr; iattr++) {
            status = EG_attributeGet(ebodys[ichild], iattr, &attrName, &attrType, &attrLen,
                                     &tempIlist, &tempRlist, &tempClist);
            CHECK_STATUS(EG_attributeGet);

            if (attrType != ATTRCSYS) continue;

            jpnt = -1;
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                if (strcmp(pnt[ipnt].name, attrName) == 0) {
                    jpnt = ipnt;
                    break;
                }
            }

            if (jpnt < 0) {
                strcpy(pnt[npnt++].name, attrName);
                jpnt = npnt-1;
            }

            bar[nbar].pnt[bar[nbar].npnt] = jpnt;
            bar[nbar].x[  bar[nbar].npnt] = tempRlist[0];
            bar[nbar].y[  bar[nbar].npnt] = tempRlist[1];
            bar[nbar].npnt++;
        }
        nbar++;
    }

#ifdef DEBUG
    printf("\ninitial setup\n");
    printall(npnt, pnt, nbar, bar);
#endif

    /* set the "fixed" flag for all known points */
    MALLOC(fixed, char, strlen(FIXED(0))+4);
    MALLOC(name,  char, strlen(FIXED(0))+4);

    snprintf(fixed, strlen(FIXED(0))+3, ";%s;", FIXED(0));

    for (ipnt = 0; ipnt < npnt; ipnt++) {
        snprintf(name, strlen(FIXED(0) )+3, ";%s;", pnt[ipnt].name);

        if (strstr(fixed, name) == NULL) continue;

        for (ibar = 0; ibar < nbar; ibar++) {
            for (jpnt = 0; jpnt < bar[ibar].npnt; jpnt++) {
                if (bar[ibar].pnt[jpnt] == ipnt) {
                    pnt[ipnt].stat = 1;
                    pnt[ipnt].x    = bar[ibar].x[jpnt];
                    pnt[ipnt].y    = bar[ibar].y[jpnt];

                    bar[ibar].nfix++;
                }
            }
        }
    }

#ifdef DEBUG
    printf("\nafter fixing known points\n");
    printall(npnt, pnt, nbar, bar);
#endif

    /* any bar with 2 or more fixed points should mark all
       its points as fixed */
    for (ibar = 0; ibar < nbar; ibar++) {
        if (bar[ibar].nfix <  2             ) continue;
        if (bar[ibar].nfix == bar[ibar].npnt) continue;

        for (ipnt = 0; ipnt < bar[ibar].npnt; ipnt++) {
            if (pnt[bar[ibar].pnt[ipnt]].stat == 0) {
                pnt[bar[ibar].pnt[ipnt]].stat = 1;
                pnt[bar[ibar].pnt[ipnt]].x    = bar[ibar].x[ipnt];
                pnt[bar[ibar].pnt[ipnt]].y    = bar[ibar].y[ipnt];

                bar[ibar].nfix++;
            }
        }
    }

#ifdef DEBUG
    printf("\nafter fixing points in bars with 2 fixed points\n");
    printall(npnt, pnt, nbar, bar);
#endif

    /* iteratively move the bars until nothing else can be done */
    for (iter = 0; iter < 2*nbar; iter++) {
        nchange = 0;

        /* for each point that is currently fixed, look for a bar that
           is currently free (nfix==0) but which uses the point */
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (pnt[ipnt].stat != 1) continue;

            for (ibar = 0; ibar < nbar; ibar++) {
                if (bar[ibar].nfix != 0) continue;

                for (jpnt = 0; jpnt < bar[ibar].npnt; jpnt++) {
                    if (bar[ibar].pnt[jpnt] == ipnt) {
#ifdef DEBUG
                        printf("translating bar %d because of point %d (%s)\n", ibar, ipnt, pnt[ipnt].name);
#endif
                        bar[ibar].dx = pnt[ipnt].x - bar[ibar].x[jpnt];
                        bar[ibar].dy = pnt[ipnt].y - bar[ibar].y[jpnt];

                        /* translate the points in the bar */
                        for (kpnt = 0; kpnt < bar[ibar].npnt; kpnt++) {
                            bar[ibar].x[kpnt] += bar[ibar].dx;
                            bar[ibar].y[kpnt] += bar[ibar].dy;
                        }

                        bar[ibar].nfix = 1;
                        nchange++;

                        break;
                    }
                }
            }
        }

#ifdef DEBUG
        printf("\nafter loop 2 for iter=%d, nchange=%d\n", iter, nchange);
        printall(npnt, pnt, nbar, bar);
#endif

        /* for each point that is currently not fixed, look for two
           bars that use the point and which are also only constrained
           at one point each (nfix==1) */
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            if (pnt[ipnt].stat != 0) continue;

            ibar1 = -1;
            ibar2 = -1;
            jpnt1 = -1;
            jpnt2 = -1;
            for (ibar = 0; ibar < nbar; ibar++) {
                if (bar[ibar].nfix != 1) continue;

                for (jpnt = 0; jpnt < bar[ibar].npnt; jpnt++) {
                    if (bar[ibar].pnt[jpnt] == ipnt) {
                        ibar1 = ibar2;
                        jpnt1 = jpnt2;
                        ibar2 = ibar;
                        jpnt2 = jpnt;
                        break;
                    }
                }
            }
            if (ibar1 < 0) continue;

            /* remember the index in ibar1 and ibar2 for ipnt */
            ipnt1 = -1;
            for (jpnt = 0; jpnt < bar[ibar1].npnt; jpnt++) {
                if (pnt[bar[ibar1].pnt[jpnt]].stat == 1) {
                    ipnt1 = jpnt;
                    break;
                }
            }

            ipnt2 = -1;
            for (jpnt = 0; jpnt < bar[ibar2].npnt; jpnt++) {
                if (pnt[bar[ibar2].pnt[jpnt]].stat == 1) {
                    ipnt2 = jpnt;
                    break;
                }
            }

            /* we have found two bars (ibar1 and ibar2), so rotate ibar about ipnt1
               and rotate ibar2 about ipnt2 such that jpnt1 in ibar1 and
               jpnt2 in ibar2 coincide.  this is done by a little
               newton solve */
#ifdef DEBUG
            printf("rotating bar %d about %d (%s) to match %d (%s)\n",
                   ibar1, ipnt1, pnt[bar[ibar1].pnt[ipnt1]].name, jpnt1, pnt[bar[ibar1].pnt[jpnt1]].name);
            printf("     and bar %d about %d (%s) to match %d (%s)\n",
                   ibar2, ipnt2, pnt[bar[ibar2].pnt[ipnt2]].name, jpnt2, pnt[bar[ibar2].pnt[jpnt2]].name);
#endif

            L0 = sqrt((bar[ibar1].y[ipnt1]-bar[ibar2].y[ipnt2]) * (bar[ibar1].y[ipnt1]-bar[ibar2].y[ipnt2])
                     +(bar[ibar1].x[ipnt1]-bar[ibar2].x[ipnt2]) * (bar[ibar1].x[ipnt1]-bar[ibar2].x[ipnt2]));
            L1 = sqrt((bar[ibar1].y[ipnt1]-bar[ibar1].y[jpnt1]) * (bar[ibar1].y[ipnt1]-bar[ibar1].y[jpnt1])
                     +(bar[ibar1].x[ipnt1]-bar[ibar1].x[jpnt1]) * (bar[ibar1].x[ipnt1]-bar[ibar1].x[jpnt1]));
            L2 = sqrt((bar[ibar2].y[ipnt2]-bar[ibar2].y[jpnt2]) * (bar[ibar2].y[ipnt2]-bar[ibar2].y[jpnt2])
                     +(bar[ibar2].x[ipnt2]-bar[ibar2].x[jpnt2]) * (bar[ibar2].x[ipnt2]-bar[ibar2].x[jpnt2]));

            if (L0 >= L1+L2 || L1 >= L2+L0 || L2 >= L0+L1) {
                printf("bar[%2d].pnt[%d] (%s) .x=%10.5f .y=%10.5f\n", ibar1, ipnt1, pnt[bar[ibar1].pnt[ipnt1]].name, bar[ibar1].x[ipnt1], bar[ibar1].y[ipnt1]);
                printf("bar[%2d].pnt[%d] (%s) .x=%10.5f .y=%10.5f  L1=%10.5f\n", ibar1, jpnt1, pnt[bar[ibar1].pnt[jpnt1]].name, bar[ibar1].x[jpnt1], bar[ibar1].y[jpnt1], L1);
                printf("bar[%2d].pnt[%d] (%s) .x=%10.5f .y=%10.5f\n", ibar2, ipnt2, pnt[bar[ibar2].pnt[ipnt2]].name, bar[ibar2].x[ipnt2], bar[ibar2].y[ipnt2]);
                printf("bar[%2d].pnt[%d] (%s) .x=%10.5f .y=%10.5f  L2=%10.5f\n", ibar2, jpnt2, pnt[bar[ibar2].pnt[jpnt2]].name, bar[ibar2].x[jpnt2], bar[ibar2].y[jpnt2], L2);
                printf("L0 (%s) to (%s)=%10.5f\n", pnt[bar[ibar1].pnt[ipnt1]].name, pnt[bar[ibar2].pnt[ipnt2]].name, L0);

                snprintf(message, 100, "incompatable distances for \"%s\", \"%s\", and \"%s\"",
                         pnt[bar[ibar1].pnt[ipnt1]].name,
                         pnt[bar[ibar1].pnt[jpnt1]].name,
                         pnt[bar[ibar2].pnt[ipnt2]].name);
                status = OCSM_UDP_ERROR2;
                goto cleanup;
            }

            phi0 = atan2(bar[ibar2].y[ipnt2]-bar[ibar1].y[ipnt1], bar[ibar2].x[ipnt2]-bar[ibar1].x[ipnt1]);
            phi1 = atan2(bar[ibar1].y[jpnt1]-bar[ibar1].y[ipnt1], bar[ibar1].x[jpnt1]-bar[ibar1].x[ipnt1]);
            phi2 = atan2(bar[ibar2].y[jpnt2]-bar[ibar2].y[ipnt2], bar[ibar2].x[jpnt2]-bar[ibar2].x[ipnt2]);

            alfa = acos((L2*L2 + L0*L0 - L1*L1) / (2 * L2 * L0));
            beta = acos((L0*L0 + L1*L1 - L2*L2) / (2 * L0 * L1));

            /* try both left and right of line and choose the one which
               most closely matches the average of jpnt1 and jpnt2 */
            xavg = (bar[ibar1].x[jpnt1] + bar[ibar2].x[jpnt2]) / 2;
            yavg = (bar[ibar1].y[jpnt1] + bar[ibar2].y[jpnt2]) / 2;

            theta1a = phi0 - phi1 - beta;
            xa      = bar[ibar1].x[ipnt1] + L1 * cos(phi1+theta1a);
            ya      = bar[ibar1].y[ipnt1] + L1 * sin(phi1+theta1a);
            dista   = sqrt((xa-xavg)*(xa-xavg) + (ya-yavg)*(ya-yavg));

            theta1b = phi0 - phi1 + beta;
            xb      = bar[ibar1].x[ipnt1] + L1 * cos(phi1+theta1b);
            yb      = bar[ibar1].y[ipnt1] + L1 * sin(phi1+theta1b);
            distb   = sqrt((xb-xavg)*(xb-xavg) + (yb-yavg)*(yb-yavg));

            if (dista < distb) {
                theta1 = phi0 - phi1 - beta;
                theta2 = phi0 - phi2 - (PI - alfa);
            } else {
                theta1 = phi0 - phi1 + beta;
                theta2 = phi0 - phi2 + (PI - alfa);
            }

            bar[ibar1].xrot = bar[ibar1].x[ipnt1];
            bar[ibar1].yrot = bar[ibar1].y[ipnt1];
            bar[ibar1].ang  = theta1;

            bar[ibar2].xrot = bar[ibar2].x[ipnt2];
            bar[ibar2].yrot = bar[ibar2].y[ipnt2];
            bar[ibar2].ang  = theta2;

            /* rotate and fix all free points in bar1 */
            cosang = cos(theta1);
            sinang = sin(theta1);

            xform11 = + cosang;
            xform12 = - sinang;
            xform13 = + (1 - cosang) * bar[ibar1].xrot + (  + sinang) * bar[ibar1].yrot;
            xform21 = + sinang;
            xform22 = + cosang;
            xform23 = + (  - sinang) * bar[ibar1].xrot + (1 - cosang) * bar[ibar1].yrot;

            for (jpnt = 0; jpnt < bar[ibar1].npnt; jpnt++) {
                xold = bar[ibar1].x[jpnt];
                yold = bar[ibar1].y[jpnt];

                bar[ibar1].x[jpnt] = xform11 * xold + xform12 * yold + xform13;
                bar[ibar1].y[jpnt] = xform21 * xold + xform22 * yold + xform23;

                if (pnt[bar[ibar1].pnt[jpnt]].stat == 0) {
                    pnt[bar[ibar1].pnt[jpnt]].x    = bar[ibar1].x[jpnt];
                    pnt[bar[ibar1].pnt[jpnt]].y    = bar[ibar1].y[jpnt];
                    pnt[bar[ibar1].pnt[jpnt]].stat = 1;
                }
            }
            bar[ibar1].nfix = bar[ibar1].npnt;

            /* rotate and fix all free points in bar2 */
            cosang = cos(theta2);
            sinang = sin(theta2);

            xform11 = + cosang;
            xform12 = - sinang;
            xform13 = + (1 - cosang) * bar[ibar2].xrot + (  + sinang) * bar[ibar2].yrot;
            xform21 = + sinang;
            xform22 = + cosang;
            xform23 = + (  - sinang) * bar[ibar2].xrot + (1 - cosang) * bar[ibar2].yrot;

            for (jpnt = 0; jpnt < bar[ibar2].npnt; jpnt++) {
                xold = bar[ibar2].x[jpnt];
                yold = bar[ibar2].y[jpnt];

                bar[ibar2].x[jpnt] = xform11 * xold + xform12 * yold + xform13;
                bar[ibar2].y[jpnt] = xform21 * xold + xform22 * yold + xform23;

                if (pnt[bar[ibar2].pnt[jpnt]].stat == 0) {
                    pnt[bar[ibar2].pnt[jpnt]].x    = bar[ibar2].x[jpnt];
                    pnt[bar[ibar2].pnt[jpnt]].y    = bar[ibar2].y[jpnt];
                    pnt[bar[ibar2].pnt[jpnt]].stat = 1;
                }
            }
            bar[ibar2].nfix = bar[ibar2].npnt;

            /* check that jpnt1 and jpnt2 actually coincide */
#ifdef DEBUG
            printf("error = %10.5f %10.5f\n", (bar[ibar1].x[jpnt1])-(bar[ibar2].x[jpnt2]),
                                              (bar[ibar1].y[jpnt1])-(bar[ibar2].y[jpnt2]));
#endif
            nchange++;

            break;
        }

#ifdef DEBUG
        printf("\nafter loop 3 for iter=%d, nchange=%d\n", iter, nchange);
        printall(npnt, pnt, nbar, bar);
#endif

        /* if we have not made any changes, break out now */
        if (nchange == 0) break;
    }

    /* make sure that whole structure is solved */
    for (ibar = 0; ibar < nbar; ibar++) {
        if (bar[ibar].nfix != bar[ibar].npnt) {
            snprintf(message, 100, "bar %d could not be placed", ibar);
            status = OCSM_UDP_ERROR1;
            goto cleanup;
        }
    }

    /* make copies of the input Bars with the appropriate transformations */
    MALLOC(newBodys, ego, nbar);

    for (ibar = 0; ibar < nbar; ibar++) {
        cosang = cos(bar[ibar].ang);
        sinang = sin(bar[ibar].ang);

        M[ 0] = + cosang;
        M[ 1] = - sinang;
        M[ 2] = 0;
        M[ 3] = + cosang * bar[ibar].dx - sinang * bar[ibar].dy + (1 - cosang) * bar[ibar].xrot + (  + sinang) * bar[ibar].yrot;

        M[ 4] = + sinang;
        M[ 5] = + cosang;
        M[ 6] = 0;
        M[ 7] = + sinang * bar[ibar].dx + cosang * bar[ibar].dy + (  - sinang) * bar[ibar].xrot + (1 - cosang) * bar[ibar].yrot;

        M[ 8] = 0;
        M[ 9] = 0;
        M[10] = 1;
        M[11] = 0;

        status = EG_makeTransform(context, M, &exform);
        CHECK_STATUS(EG_makeTransform);

        status = EG_copyObject(ebodys[ibar], exform, &newBodys[ibar]);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(exform);
        CHECK_STATUS(EG_deleteObject);
    }

    /* make a Model with the moved Bars */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, nbar, newBodys, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* set the output value */

    /* return the transformed Model */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(fixed);
    FREE(name );
    FREE(newBodys);
    FREE(pnt);
    FREE(bar);

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


/*
 ************************************************************************
 *                                                                      *
 *   printAll - print point and bar information                         *
 *                                                                      *
 ************************************************************************
 */

#ifdef DEBUG
static void
printall(int    npnt,                   /* (in)  number of points */
         pnt_T  pnt[],                  /* (in)  array  of points */
         int    nbar,                   /* (in)  number of bars */
         bar_T  bar[])                  /* (in)  array  of bars */
{
    int     ipnt, ibar;

    /* --------------------------------------------------------------- */

    printf(" ipnt name     stat      x          y\n");
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        printf("%5d %-10s %2d %10.5f %10.5f\n", ipnt,
               pnt[ipnt].name, pnt[ipnt].stat, pnt[ipnt].x, pnt[ipnt].y);
    }

    printf(" ibar  npnt  nfix        dx         dy        xrot       yrot       ang\n");
    for (ibar = 0; ibar < nbar; ibar++) {
        printf("%5d %5d %5d   %10.5f %10.5f %10.5f %10.5f %10.5f\n", ibar,
               bar[ibar].npnt, bar[ibar].nfix, bar[ibar].dx, bar[ibar].dy, bar[ibar].xrot, bar[ibar].yrot, bar[ibar].ang);
        for (ipnt = 0; ipnt < bar[ibar].npnt; ipnt++) {
            printf("                                                                          %2d (%s) %10.5f %10.5f\n",
                   bar[ibar].pnt[ipnt], pnt[bar[ibar].pnt[ipnt]].name, bar[ibar].x[ipnt], bar[ibar].y[ipnt]);
        }
    }
}
#endif

