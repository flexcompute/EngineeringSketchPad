/*
 ************************************************************************
 *                                                                      *
 * udpBiconvex -- biconvex airfoil between (0,0) and (1,0)              *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

#define UDP    1
#define EPS06  1.0e-6

/*
 * Copyright (C) 2013/2021  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 2
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#ifdef UDP
    #define THICK(     IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
    #define THICK_DOT( IUDP)  ((double *) (udps[IUDP].arg[0].dot))[0]
    #define CAMBER(    IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
    #define CAMBER_DOT(IUDP)  ((double *) (udps[IUDP].arg[1].dot))[0]

    /* data about possible arguments */
    static char  *argNames[NUMUDPARGS] = {"thick",     "camber",    };
    static int    argTypes[NUMUDPARGS] = {ATTRREALSEN, ATTRREALSEN, };
    static int    argIdefs[NUMUDPARGS] = {0,           0,           };
    static double argDdefs[NUMUDPARGS] = {0.,          0.,          };

    /* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                             udpGet, udpVel, udpClean, udpMesh */
    #include "udpUtilities.c"
#else
    #define THICK(     IUDP)  0.12
    #define THICK_DOT( IUDP)  0.00
    #define CAMBER(    IUDP)  0.04
    #define CAMBER_DOT(IUDP)  0.00

/*@null@*/ char *
udpErrorStr(int stat)                   /* (in)  status number */
{
    char *string;                       /* (out) error message */

    MALLOC(string, char, 25);
    snprintf(string, 25, "EGADS status = %d", stat);

    return string;
}

int
main(int       argc,                    /* (in)  number of arguments */
     char      *argv[])                 /* (in)  array of arguments */
{
    int  status, nMesh;
    char *string;
    ego  context, ebody, emodel;

    /* dummy call to prevent compiler warnings */
    ocsmPrintEgo(NULL);

    /* define a context */
    status = EG_open(&context);
    printf("EG_open -> status=%d\n", status);
    if (status < 0) exit(1);

    /* call the execute routine */
    status = udpExecute(context, &ebody, &nMesh, &string);
    printf("udpExecute -> status=%d\n", status);
    if (status < 0) exit(1);

    EG_free(string);

    /* make and dump the model */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL,
                             1, &ebody, NULL, &emodel);
    printf("EG_makeTopology -> status=%d\n", status);
    if (status < 0) exit(1);

    status = EG_saveModel(emodel, "biconvex.egads");
    printf("EG_saveModel -> status=%d\n", status);
    if (status < 0) exit(1);

    /* cleanup */
    status = EG_deleteObject(emodel);
    printf("EG_close -> status=%d\n", status);

    status = EG_close(context);
    printf("EG_close -> status=%d\n", status);
}
#endif


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
    int     iedge, sense[2];
    double  L, H, R, nodeLE[3], nodeTE[3], circ[10], trange[2], data[18];
    ego     enodes[3], ecurves[2], eedges[2], eloop, eface;

    ROUTINE(udpExecute);
    
#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("thick(0)       = %f\n", THICK(     0));
    printf("thick_dot(0)   = %f\n", THICK_DOT (0));
    printf("camber(0)      = %f\n", CAMBER(    0));
    printf("camber_dot(0)  = %f\n", CAMBER_DOT(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

#ifdef UDP
    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        printf(" udpExecute: thick should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (THICK(0) <= 0) {
        printf(" udpExecute: thick=%f < 0\n", THICK(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        printf(" udpExecute: camber should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);
#endif

#ifdef DEBUG
    printf("thick[%d]      = %f\n", numUdp, THICK(     numUdp));
    printf("thick_dot[%d]  = %f\n", numUdp, THICK_DOT( numUdp));
    printf("camber[%d]     = %f\n", numUdp, CAMBER(    numUdp));
    printf("camber_dot[%d] = %f\n", numUdp, CAMBER_DOT(numUdp));
#endif

    /* compute the geometry */
    for (iedge = 0; iedge < 2; iedge++) {
        if (iedge == 0) {
            H = CAMBER(0) + THICK(0) / 2;
        } else {
            H = CAMBER(0) - THICK(0) / 2;
        }

        /* circle which is convex up */
        if (H > EPS06) {
            L = 0.50;
            R = (L*L + H*H) / (2*H);

            circ[0] = L;
            circ[1] = H - R;
            circ[2] = 0;
            circ[3] = 1;
            circ[4] = 0;
            circ[5] = 0;
            circ[6] = 0;
            circ[7] = 1;
            circ[8] = 0;
            circ[9] = R;

            status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, circ, &(ecurves[iedge]));
            CHECK_STATUS(EG_makeGeometry);

        } else if (H < -EPS06) {
            L = 0.50;
            R = -(L*L + H*H) / (2*H);

            circ[0] = L;
            circ[1] = H + R;
            circ[2] = 0;
            circ[3] = 1;
            circ[4] = 0;
            circ[5] = 0;
            circ[6] = 0;
            circ[7] = 1;
            circ[8] = 0;
            circ[9] = R;

            status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL, circ, &(ecurves[iedge]));
            CHECK_STATUS(EG_makeGeometry);

        } else if (iedge == 0) {
            circ[0] =  1;
            circ[1] =  0;
            circ[2] =  0;
            circ[3] = -1;
            circ[4] =  0;
            circ[5] =  0;

            status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, circ, &(ecurves[iedge]));
            CHECK_STATUS(EG_makeGeometry);

        } else {
            circ[0] =  0;
            circ[1] =  0;
            circ[2] =  0;
            circ[3] =  1;
            circ[4] =  0;
            circ[5] =  0;

            status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, circ, &(ecurves[iedge]));
            CHECK_STATUS(EG_makeGeometry);
        }
    }

    /* make Nodes */
    nodeLE[0] = 0;
    nodeLE[1] = 0;
    nodeLE[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 0, nodeLE, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    nodeTE[0] = 1;
    nodeTE[1] = 0;
    nodeTE[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 0, nodeTE, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    enodes[2] = enodes[0];

    /* make Edge for upper surface */
    if (CAMBER(0)+THICK(0)/2 >= 0) {
        status = EG_invEvaluate(ecurves[0], nodeTE, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurves[0], nodeLE, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[1] < trange[0]) trange[1] += TWOPI;

        status = EG_makeTopology(context, ecurves[0], EDGE, TWONODE, trange, 2, &(enodes[1]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        sense[0] = SFORWARD;
    } else {
        status = EG_invEvaluate(ecurves[0], nodeTE, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurves[0], nodeLE, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[1] < trange[0]) trange[1] += TWOPI;

        status = EG_makeTopology(context, ecurves[0], EDGE, TWONODE, trange, 2, &(enodes[0]), NULL, &(eedges[0]));
        CHECK_STATUS(EG_makeTopology);

        sense[0] = SREVERSE;
    }

    /* make Edge for lower surface */
    if (CAMBER(0)-THICK(0)/2 <= 0) {
        status = EG_invEvaluate(ecurves[1], nodeLE, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurves[1], nodeTE, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[1] < trange[0]) trange[1] += TWOPI;

        status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE, trange, 2, &(enodes[0]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        sense[1] = SFORWARD;
    } else {
        status = EG_invEvaluate(ecurves[1], nodeTE, &(trange[0]), data);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_invEvaluate(ecurves[1], nodeLE, &(trange[1]), data);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[1] < trange[0]) trange[1] += TWOPI;

        status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE, trange, 2, &(enodes[1]), NULL, &(eedges[1]));
        CHECK_STATUS(EG_makeTopology);

        sense[1] = SREVERSE;
    }

    /* make Loop from these Edges */
    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 2, eedges, sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    /* make Face from the loop */
    status = EG_makeFace(eloop, SFORWARD, NULL, &eface);
    CHECK_STATUS(EG_makeFace);

    /* create the FaceBody (which will be returned) */
    status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &eface, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* no output value(s) */

    /* remember this model (Body) */
#ifdef UDP
    udps[numUdp].ebody = *ebody;
#endif

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
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
    int    status = EGADS_SUCCESS;

    int    iudp, judp, ipnt;
    double H, H_dot, L, R, R_dot;
    double thbeg, thbeg_dot, thend, thend_dot, s, th, th_dot;

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

    /* velocity of a Node is zero */
    if (entType == OCSM_NODE) {
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            vels[3*ipnt  ] = 0;
            vels[3*ipnt+1] = 0;
            vels[3*ipnt+2] = 0;
        }

    /* velocity of a Face is zero */
    } else if (entType == OCSM_FACE) {
        for (ipnt = 0; ipnt < npnt; ipnt++) {
            vels[3*ipnt  ] = 0;
            vels[3*ipnt+1] = 0;
            vels[3*ipnt+2] = 0;
        }

        /* velocity of the Edges */
    } else {
        if (entIndex == 1) {
            H     = CAMBER(    iudp) + THICK(    iudp) / 2;
            H_dot = CAMBER_DOT(iudp) + THICK_DOT(iudp) / 2;
        } else {
            H     = CAMBER(    iudp) - THICK(    iudp) / 2;
            H_dot = CAMBER_DOT(iudp) - THICK_DOT(iudp) / 2;
        }

        /* circle which is convex up */
        if (H > EPS06) {
            L = 0.5;

            R     = (H*H + L*L) / (2*H);
            R_dot = H_dot * (1 - (H*H + L*L) / (2 * H*H));

            thbeg     = atan2(R-H, L);
            thbeg_dot = (R_dot - H_dot) * L / ((R-H) * (R-H) + L*L);

            thend     = PI - thbeg;
            thend_dot =    - thbeg_dot;

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                s = (uvs[ipnt] - thbeg) / (thend - thbeg);

                th     = thbeg     * (1-s) + thend     * s;
                th_dot = thbeg_dot * (1-s) + thend_dot * s;

                vels[3*ipnt  ] = R_dot * cos(th) - th_dot * R * sin(th);
                vels[3*ipnt+1] = R_dot * sin(th) + th_dot * R * cos(th) - R_dot + H_dot;
                vels[3*ipnt+2] = 0;
            }

        /* circle which is convex down */
        } else if (H < -EPS06) {
            L = 0.5;

            R     = -(H*H + L*L) / (2*H);
            R_dot = -H_dot * (1 - (H*H + L*L) / (2 * H*H));

            thbeg     = PI + atan2(R+H, L);
            thbeg_dot = (R_dot + H_dot) * L / ((R+H) * (R+H) + L*L);

            thend     = 3*PI - thbeg;
            thend_dot =      - thbeg_dot;

            for (ipnt = 0; ipnt < npnt; ipnt++) {
                s = (uvs[ipnt] - thbeg) / (thend - thbeg);

                th     = thbeg     * (1-s) + thend     * s;
                th_dot = thbeg_dot * (1-s) + thend_dot * s;

                vels[3*ipnt  ] = R_dot * cos(th) - th_dot * R * sin(th);
                vels[3*ipnt+1] = R_dot * sin(th) + th_dot * R * cos(th) + R_dot + H_dot;
                vels[3*ipnt+2] = 0;
            }

        /* straight line */
        } else {
            for (ipnt = 0; ipnt < npnt; ipnt++) {
                vels[3*ipnt  ] = 0;
                vels[3*ipnt+1] = 0;
                vels[3*ipnt+2] = 0;
            }
        }
    }

//cleanup:
    return status;
}
