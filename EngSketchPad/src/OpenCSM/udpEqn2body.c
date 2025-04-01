/*
 ************************************************************************
 *                                                                      *
 * udpEqn2body -- udp file to generate WireBody or SheetBody from eqn   *
 *                                                                      *
 *            Written by Sydney Jud @ Syracuse University               *
 *            Modified by John Dannenhoffer                             *
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

#define NUMUDPARGS 8
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define XEQN(  IUDP)    ((char   *) (udps[IUDP].arg[0].val))
#define YEQN(  IUDP)    ((char   *) (udps[IUDP].arg[1].val))
#define ZEQN(  IUDP)    ((char   *) (udps[IUDP].arg[2].val))
#define URANGE(IUDP,I)  ((double *) (udps[IUDP].arg[3].val))[I]
#define VRANGE(IUDP,I)  ((double *) (udps[IUDP].arg[4].val))[I]
#define TOLER( IUDP)    ((double *) (udps[IUDP].arg[5].val))[0]
#define NPNT(  IUDP)    ((int    *) (udps[IUDP].arg[6].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"xeqn",     "yeqn",     "zeqn",     "urange",  "vrange",  "toler",  "npnt",  "rebuild",   };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRSTRING, ATTRSTRING, ATTRREAL,  ATTRREAL,  ATTRREAL, ATTRINT, ATTRREBUILD, };
static int    argIdefs[NUMUDPARGS] = {0,          0,          0,          0,         0,         0,        101,     0,           };
static double argDdefs[NUMUDPARGS] = {1.,         0.,         0.,         0.,        0.,        1.e-5,    0.,      0.,          };


/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototype for function defined below */
static int  makeWireBody(ego context, int iudp, char message[], ego *ebody, int *NumUdp, udp_T *udps);
static int  makeSheetBody(ego context, int iudp, char message[], ego *ebody, int *NumUdp, udp_T *udps);


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

#ifdef DEBUG
    int     i;
#endif
    char    *message=NULL;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("XEQN(  0) = %s\n", XEQN(0));
    printf("YEQN(  0) = %s\n", YEQN(0));
    printf("ZEQN(  0) = %s\n", ZEQN(0));
    printf("URANGE(0) = ");
    for (i = 0; i < udps[0].arg[3].size; i++) {
        printf(" %f", URANGE(0,i));
    }
    printf("\n");
    printf("VRANGE(0) = ");
    for (i = 0; i < udps[0].arg[4].size; i++) {
        printf(" %f", VRANGE(0,i));
    }
    printf("\n");
    printf("TOLER( 0) = %f\n", TOLER(0));
    printf("NPNT(  0) = %d\n", NPNT( 0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 1024);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[3].size <= 1) {
        snprintf(message, 1024, "URANGE must have 2 values");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (URANGE(0,1) <= URANGE(0,0)) {
        snprintf(message, 1024, "URANGE must specify a positive interval");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[4].size == 1 && VRANGE(0,0) != 0) {
        snprintf(message, 1024, "VRANGE must specify a positive interval (if specified)");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[4].size >= 2 && VRANGE(0,1) <= VRANGE(0,0)) {
        snprintf(message, 1024, "VRANGE must specify a positive interval (if specified)");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[5].size != 1) {
        snprintf(message, 1024, "TOLER must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (TOLER(0) < 0) {
        snprintf(message, 1024, "TOLER (=%f) must be non-negative", TOLER(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[6].size != 1) {
        snprintf(message, 1024, "NPNT must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (NPNT(0) < 5) {
        snprintf(message, 1024, "NPNT (=%d) must be >= 5", NPNT(0));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("XEQN(  %d) = %s\n", numUdp, XEQN(numUdp));
    printf("YEQN(  %d) = %s\n", numUdp, YEQN(numUdp));
    printf("ZEQN(  %d) = %s\n", numUdp, ZEQN(numUdp));
    printf("URANGE(%d) = ",     numUdp);
    for (i = 0; i < udps[numUdp].arg[3].size; i++) {
        printf(" %f", URANGE(numUdp,i));
    }
    printf("\n");
    printf("VRANGE(%d) = ",     numUdp);
    for (i = 0; i < udps[numUdp].arg[4].size; i++) {
        printf(" %f", VRANGE(numUdp,i));
    }
    printf("\n");
    printf("TOLER( %d) = %f\n", numUdp, TOLER(numUdp));
    printf("NPNT(  %d) = %d\n", numUdp, NPNT( numUdp));
#endif

    /* determine if a WireBody or SheetBody is to be made */
    if (udps[0].arg[4].size == 1) {
        status = makeWireBody(context, 0, message, ebody, NumUdp, udps);
        CHECK_STATUS(makeWireBody);
    } else {
        status = makeSheetBody(context, 0, message, ebody, NumUdp, udps);
        CHECK_STATUS(makeSheetBody);
    }

    /* set the output value(s) */

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
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
 *   makeWireBody - make a WireBody by evaluating the eqns              *
 *                                                                      *
 ************************************************************************
 */

static int
makeWireBody(ego     context,           /* (in)  EGADS context */
             int     iudp,              /* (in)  udp index */
             char    message[],         /* (out) error message */
             ego     *ebody,            /* (out) created WireBody */
 /*@unused@*/int     *NumUdp,
             udp_T   *udps)
{
    int      status = EGADS_SUCCESS;    /* (out) return status */

    int      ipmtr, ipnt, sizes[2], periodic, sense, status2;
    double   *xyz=NULL, uu, vel, urange[2];
    char     tempStr[256];
    ego      ecurve, enodes[2], eedge, eloop;

    void     *modl;
    modl_T   *MODL=NULL;

    ROUTINE(makeWireBody);

    /* --------------------------------------------------------------- */

    /* get a pointer to the MODL */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)modl;

    /* check to see if "u" already exists as a paramater */
    status = ocsmFindPmtr(MODL, "u", 0, 0, 0, &ipmtr);
    CHECK_STATUS(ocsmFindPmtr);

    if (ipmtr != 0) {
        snprintf(message, 1024, "Parameter \"u\" already exists");
        status = OCSM_NAME_ALREADY_DEFINED;
        goto cleanup;
    } else {
        status = ocsmFindPmtr(MODL, "u", OCSM_LOCALVAR, 1, 1, &ipmtr);
        CHECK_STATUS(ocsmFindPmtr);
    }

    /* get storage for the points */
    MALLOC(xyz, double, 3*NPNT(iudp));

    /* create the array of points by evaluating the given expressions */
    for (ipnt = 0; ipnt < NPNT(iudp); ipnt++) {
        uu = URANGE(iudp,0) + ((double)ipnt)/((double)(NPNT(iudp)-1)) * (URANGE(iudp,1) - URANGE(iudp,0));

        status = ocsmSetValuD(MODL, ipmtr, 1, 1, uu);
        CHECK_STATUS(ocsmSetValuD);

        if (udps[iudp].arg[0].size == 0) {
            xyz[3*ipnt  ] = 0;
        } else {
            status = ocsmEvalExpr(MODL, XEQN(iudp), &(xyz[3*ipnt  ]), &vel, tempStr);
            CHECK_STATUS(ocsmEvalExpr);

            if (strlen(tempStr) > 0) {
                snprintf(message, 1024, "XEQN should evaluate to a number");
                status = OCSM_ILLEGAL_TYPE;
                goto cleanup;
            }
        }

        if (udps[iudp].arg[1].size == 0) {
            xyz[3*ipnt+1] = 0;
        } else {
            status = ocsmEvalExpr(MODL, YEQN(iudp), &(xyz[3*ipnt+1]), &vel, tempStr);
            CHECK_STATUS(ocsmEvalExpr);

            if (strlen(tempStr) > 0) {
                snprintf(message, 1024, "YEQN should evaluate to a number");
                status = OCSM_ILLEGAL_TYPE;
                goto cleanup;
            }
        }

        if (udps[iudp].arg[2].size == 0) {
            xyz[3*ipnt+2] = 0;
        } else {
            status = ocsmEvalExpr(MODL, ZEQN(iudp), &(xyz[3*ipnt+2]), &vel, tempStr);
            CHECK_STATUS(ocsmEvalExpr);

            if (strlen(tempStr) > 0) {
                snprintf(message, 1024, "ZEQN should evaluate to a number");
                status = OCSM_ILLEGAL_TYPE;
                goto cleanup;
            }
        }

#ifdef DEBUG
        printf("%10.5f   %10.5f %10.5f %10.5f\n", uu, xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2]);
#endif
    }

    /* generate the Curve */
    sizes[0] = NPNT(iudp);
    sizes[1] = 0;

    status = EG_approximate(context, 0, TOLER(iudp), sizes, xyz, &ecurve);
    CHECK_STATUS(EG_approximate);

    status = EG_getRange(ecurve, urange, &periodic);
    CHECK_STATUS(EG_getRange);

    /* generate the Nodes at the ends of the Curve */
    status = EG_makeTopology(context, NULL, NODE, 0,
                             &(xyz[0]), 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, NODE, 0,
                             &(xyz[3*NPNT(iudp)-3]), 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    /* create the Edge, Loop, and WireBody */
    sense = SFORWARD;
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             urange, 2, enodes, &sense, &eedge);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, LOOP, OPEN,
                             NULL, 1, &eedge, &sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                             NULL, 1, &eloop, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

cleanup:
    /* remove the temporary LOCALVAR */
    status2 = ocsmFindPmtr(MODL, "u", 0, 0, 0, &ipmtr);
    if (status2 == EGADS_SUCCESS && ipmtr > 0) {
        (void) ocsmDelPmtr(MODL, ipmtr);
    }

    FREE(xyz);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   makeSheetBody - make a SheetBody by evaluating the eqns            *
 *                                                                      *
 ************************************************************************
 */

static int
makeSheetBody(ego     context,          /* (in)  EGADS context */
              int     iudp,              /* (in)  udp index */
              char    message[],         /* (out) error message */
              ego     *ebody,            /* (out) created SheetBody */
  /*@unused@*/int     *NumUdp,
              udp_T   *udps)
{
    int      status = EGADS_SUCCESS;    /* (out) return status */

    int      ipmtru, ipmtrv, ipntu, ipntv, kpnt, sizes[2], periodic, status2;
    double   *xyz=NULL, uu, vv, vel, uvrange[4];
    char     tempStr[256];
    ego      esurface, eface, eshell;

    void     *modl;
    modl_T   *MODL=NULL;

    ROUTINE(makeSheetBody);

    /* --------------------------------------------------------------- */

    /* get a pointer to the MODL */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)modl;

    /* check to see if "u" or "v" already exists as a paramater */
    status = ocsmFindPmtr(MODL, "u", 0, 0, 0, &ipmtru);
    CHECK_STATUS(ocsmFindPmtr);

    if (ipmtru != 0) {
        snprintf(message, 1024, "Parameter \"u\" already exists");
        status = OCSM_NAME_ALREADY_DEFINED;
        goto cleanup;
    } else {
        status = ocsmFindPmtr(MODL, "u", OCSM_LOCALVAR, 1, 1, &ipmtru);
        CHECK_STATUS(ocsmFindPmtr);
    }

    status = ocsmFindPmtr(MODL, "v", 0, 0, 0, &ipmtrv);
    CHECK_STATUS(ocsmFindPmtr);

    if (ipmtrv != 0) {
        snprintf(message, 1024, "Parameter \"v\" already exists");
        status = OCSM_NAME_ALREADY_DEFINED;
        goto cleanup;
    } else {
        status = ocsmFindPmtr(MODL, "v", OCSM_LOCALVAR, 1, 1, &ipmtrv);
        CHECK_STATUS(ocsmFindPmtr);
    }

    /* get storage for the points */
    MALLOC(xyz, double, 3*NPNT(iudp)*NPNT(iudp));

    /* create the array of points by evaluating the given expressions */
    kpnt = 0;
    for (ipntu = 0; ipntu < NPNT(iudp); ipntu++) {
        for (ipntv = 0; ipntv < NPNT(iudp); ipntv++) {
        
            uu = URANGE(iudp,0) + ((double)ipntu)/((double)(NPNT(iudp)-1)) * (URANGE(iudp,1) - URANGE(iudp,0));
            vv = VRANGE(iudp,0) + ((double)ipntv)/((double)(NPNT(iudp)-1)) * (VRANGE(iudp,1) - VRANGE(iudp,0));

            status = ocsmSetValuD(MODL, ipmtru, 1, 1, uu);
            CHECK_STATUS(ocsmSetValuD);

            status = ocsmSetValuD(MODL, ipmtrv, 1, 1, vv);
            CHECK_STATUS(ocsmSetValuD);

            if (udps[iudp].arg[0].size == 0) {
                xyz[3*kpnt  ] = 0;
            } else {
                status = ocsmEvalExpr(MODL, XEQN(iudp), &(xyz[3*kpnt  ]), &vel, tempStr);
                CHECK_STATUS(ocsmEvalExpr);

                if (strlen(tempStr) > 0) {
                    snprintf(message, 1024, "XEQN should evaluate to a number");
                    status = OCSM_ILLEGAL_TYPE;
                    goto cleanup;
                }
            }

            if (udps[iudp].arg[1].size == 0) {
                xyz[3*kpnt+1] = 0;
            } else {
                status = ocsmEvalExpr(MODL, YEQN(iudp), &(xyz[3*kpnt+1]), &vel, tempStr);
                CHECK_STATUS(ocsmEvalExpr);

                if (strlen(tempStr) > 0) {
                    snprintf(message, 1024, "YEQN should evaluate to a number");
                    status = OCSM_ILLEGAL_TYPE;
                    goto cleanup;
                }
            }

            if (udps[iudp].arg[2].size == 0) {
                xyz[3*kpnt+2] = 0;
            } else {
                status = ocsmEvalExpr(MODL, ZEQN(iudp), &(xyz[3*kpnt+2]), &vel, tempStr);
                CHECK_STATUS(ocsmEvalExpr);

                if (strlen(tempStr) > 0) {
                    snprintf(message, 1024, "ZEQN should evaluate to a number");
                    status = OCSM_ILLEGAL_TYPE;
                    goto cleanup;
                }
            }

#ifdef DEBUG
            printf("%10.5f %10.5f   %10.5f %10.5f %10.5f\n", uu, vv, xyz[3*kpnt], xyz[3*kpnt+1], xyz[3*kpnt+2]);
#endif
            kpnt++;
        }
    }

    /* generate the Curve */
    sizes[0] = NPNT(iudp);
    sizes[1] = NPNT(iudp);

    status = EG_approximate(context, 0, 0, sizes, xyz, &esurface);
    CHECK_STATUS(EG_approximate);

    status = EG_getRange(esurface, uvrange, &periodic);
    CHECK_STATUS(EG_getRange);

    /* generate a Face */
    status = EG_makeFace(esurface, SFORWARD, uvrange, &eface);
    CHECK_STATUS(EG_makeFace);

    /* make the Shell and SheetBody */
    status = EG_makeTopology(context, NULL, SHELL, OPEN,
                             NULL, 1, &eface, NULL, &eshell);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                             NULL, 1, &eshell, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

cleanup:
    /* remove the temporary LOCALVAR */
    status2 = ocsmFindPmtr(MODL, "", 0, 0, 0, &ipmtru);
    if (status2 == EGADS_SUCCESS && ipmtru > 0) {
        (void) ocsmDelPmtr(MODL, ipmtru);
    }

    status2 = ocsmFindPmtr(MODL, "", 0, 0, 0, &ipmtrv);
    if (status2 == EGADS_SUCCESS && ipmtrv > 0) {
        (void) ocsmDelPmtr(MODL, ipmtrv);
    }

    FREE(xyz);

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
