/*
 ************************************************************************
 *                                                                      *
 * udpPod -- udp file to generate a OpenVSP pod                         *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2022  John F. Dannenhoffer, III (Syracuse University)
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

#define  UDP        1
#define  NUMUDPARGS 3
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#ifdef UDP
   #define LENGTH(  IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
   #define FINENESS(IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
   #define VOLUME(  IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]

   /* data about possible arguments */
   static char  *argNames[NUMUDPARGS] = {"length", "fineness", "volume" , };
   static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL,   -ATTRREAL, };
   static int    argIdefs[NUMUDPARGS] = {0,        0,          0,         };
   static double argDdefs[NUMUDPARGS] = {0.,       0.,         0.,        };

   /* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                            udpGet, udpVel, udpClean, udpMesh */
   #include "udpUtilities.c"
#else
   #define LENGTH(  IUDP)  5.0
   #define FINENESS(IUDP)  6.0
   #define VOLUME(  IUDP)  myVolume

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
    CHECK_STATUS(EG_open);

    /* call the execute routine */
    status = udpExecute(context, &ebody, &nMesh, &string);
    CHECK_STATUS(udpExecute);

    EG_free(string);

    /* make and dump the model */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL,
                             1, &ebody, NULL, &emodel);
    CHECK_STATUS(EG_makeTopology);

    status = EG_saveModel(emodel, "pod.egads");
    CHECK_STATUS(EG_saveModel);

    /* cleanup */
    status = EG_deleteObject(emodel);
    CHECK_STATUS(EG_deleteObject);

    status = EG_close(context);
    CHECK_STATUS(EG_close);

    return EXIT_SUCCESS;
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

    int     header[3], sense[1], oclass, mtype, nchild, *senses;
    double  node0[3], node1[3], depth, cp[21], data[18], trange[2];
    double  axis[6] = { 0, 0, 0, 1, 0, 0};
    ego     enodes[2], ecurve, eedge, eloop, eref, *echilds, elist[2];
    ego     ewire, esheet;

#ifndef UDP
    double  myVolume;
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("length(0)         = %f\n", LENGTH(0));
    printf("fineness(0)       = %f\n", FINENESS(0));
#endif

    /* default return values */
    *ebody = NULL;
    *nMesh = 0;
    *string = NULL;

#ifdef UDP
    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        printf(" udpExecute: length should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (LENGTH(0) <= 0) {
        printf(" udpExecute: length = %f <= 0\n", LENGTH(0));
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (udps[0].arg[1].size > 1) {
        printf(" udpExecute: fineness should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (FINENESS(0) <= 0) {
        printf(" udpExecute: fine_ratio = %f <= 0\n", FINENESS(0));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);
    
#endif

#ifdef DEBUG
    printf("length[  %d] = %f\n", numUdp, LENGTH(  numUdp));
    printf("fineness[%d] = %f\n", numUdp, FINENESS(numUdp));

#endif

    /* Node locations */
    depth = LENGTH(0) / FINENESS(0);

    /* make Nodes */
    node0[0] = 0;           node0[1] = 0;  node0[2] = 0;
    node1[0] = LENGTH(0);   node1[1] = 0;  node1[2] = 0;

    status = EG_makeTopology(context, NULL, NODE, 0, node0, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, NODE, 0, node1, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    /* make (bezier) Edge 0 */
    header[0] = 0;  // bitFlag
    header[1] = 3;  // degree
    header[2] = 7;  // nCP

    cp[ 0] = 0;                cp[ 1] = 0;            cp[ 2] = 0;
    cp[ 3] = 0.05*LENGTH(0);   cp[ 4] = 0.95*depth;   cp[ 5] = 0;
    cp[ 6] = 0.20*LENGTH(0);   cp[ 7] =      depth;   cp[ 8] = 0;
    cp[ 9] = 0.50*LENGTH(0);   cp[10] =      depth;   cp[11] = 0;
    cp[12] = 0.60*LENGTH(0);   cp[13] =      depth;   cp[14] = 0;
    cp[15] = 0.95*LENGTH(0);   cp[16] = 0.30*depth;   cp[17] = 0;
    cp[18] =      LENGTH(0);   cp[19] = 0;            cp[20] = 0;

    status = EG_makeGeometry(context, CURVE, BEZIER, NULL,
                             header, cp, &ecurve);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_invEvaluate(ecurve, node0, &(trange[0]), data);
    CHECK_STATUS(EG_invEvaluate);

    status = EG_invEvaluate(ecurve, node1, &(trange[1]), data);
    CHECK_STATUS(EG_invEvaluate);

    elist[0] = enodes[0];
    elist[1] = enodes[1];
    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, elist, NULL, &eedge);
    CHECK_STATUS(EG_makeTopology);

    /* make a WireBody */
    sense[0] = SFORWARD;
    status = EG_makeTopology(context, NULL, LOOP, OPEN,
                             NULL, 1, &eedge, sense, &eloop);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                             NULL, 1, &eloop, NULL, &ewire);
    CHECK_STATUS(EG_makeTopology);

    /* rotate WireBody into SheetBody */
    status = EG_rotate(ewire, +360., axis, &esheet);
    CHECK_STATUS(EG_rotate);

    /* convert SheetBody into SolidBody */
    status = EG_getTopology(esheet, &eref, &oclass, &mtype, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_makeTopology(context, NULL, BODY, SOLIDBODY,
                             NULL, 1, echilds, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    /* set the output value(s) */
    status = EG_getMassProperties(*ebody, data);
    CHECK_STATUS(EG_getMassProperties);

    VOLUME(0) = data[0];

    /* remember this model (body) */
#ifdef UDP
    udps[numUdp].ebody = *ebody;
#else
    printf("myVolume = %f\n", myVolume);
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

