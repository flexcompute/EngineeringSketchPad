/*
 ************************************************************************
 *                                                                      *
 * udfGuide -- move cross-section along a guide curve                   *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

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

#define NUMUDPINPUTBODYS 2
#define NUMUDPARGS       3

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define NXSECT(IUDP)     ((int    *) (udps[IUDP].arg[0].val))[0]
#define ORIGIN(IUDP,I)   ((double *) (udps[IUDP].arg[1].val))[I]
#define AXIS(  IUDP,I)   ((double *) (udps[IUDP].arg[2].val))[I]

static char  *argNames[NUMUDPARGS] = {"nxsect", "origin", "axis",   };
static int    argTypes[NUMUDPARGS] = {ATTRINT,  ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {5,        0,        0,        };
static double argDdefs[NUMUDPARGS] = {0.,       0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "OpenCSM.h"

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

    int     i, oclass, mtype, nchild, *senses, periodic, nloop, nedge;
    double  origin[3], data[18], trange[4], tt, theta, xform[16];
    ego     context, eref, *ebodys, *echilds, *eloops, *eedges;
    ego     eface, eloop, exform, *exsects=NULL;

    ROUTINE(udpExecute);
    
#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("nxsect    = %d\n", NXSECT(0));
    printf("origin(0) = ");
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf("%f ", ORIGIN(0,i));
    }
    printf("\n");
    printf("axis(  0) = ");
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf("%f ", AXIS(0,i));
    }
    printf("\n");
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check/process arguments */
    if (udps[0].arg[0].size > 1) {
        printf(" udpExecute: nxsect should be a scalar\n");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (NXSECT(0) <= 0) {
        printf(" udpExecute: nxsect = %d <= 0\n", NXSECT(0));
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size == 1 && AXIS(0,0) == 0) {
    } else if (udps[0].arg[2].size == 6) {
    } else {
        printf(" udpExecute: axis must be blank or have 6 elements\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* get cross-section origin */
    if (udps[0].arg[1].size == 3) {
        origin[0] = ORIGIN(0,0);
        origin[1] = ORIGIN(0,1);
        origin[2] = ORIGIN(0,2);
    } else {
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
    }

    /* check that Model was input that contains two Bodys */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        printf(" udpExecute: expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 2) {
        printf(" udpExecute: expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

#ifdef DEBUG
    printf("emodel\n");
    ocsmPrintEgo(emodel);
#endif

    /* make sure first Body (the cross-section) is a FaceBody, SheetBody, or WireBody
       and the second Body (the path) is a WireBody */
    status = EG_getTopology(ebodys[0], &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != BODY || (mtype != WIREBODY && mtype != FACEBODY && mtype != SHEETBODY)) {
            printf(" udpExecute: left Body must be WireBody or SheetBody\n");
            status = EGADS_NOTBODY;
            goto cleanup;
    }

    status = EG_getTopology(ebodys[1], &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != BODY || mtype != WIREBODY) {
            printf(" udpExecute: RITE Body must be WireBody\n");
            status = EGADS_NOTBODY;
            goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("nxsect(%d) = %d\n", numUdp, NXSECT(numUdp));
    printf("origin(%d) = ", numUdp);
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf("%f ", ORIGIN(numUdp,i));
    }
    printf("\n");
    printf("axis(  0) = ");
    for (i = 0; i < udps[numUdp].arg[2].size; i++) {
        printf("%f ", AXIS(numUdp,i));
    }
    printf("\n");
#endif

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* get an array to hold the cross-sections */
    MALLOC(exsects, ego, NXSECT(numUdp));

    /* for now, make sure guide curve is comprised of a single Edge */
    status = EG_getTopology(ebodys[1], &eref, &oclass, &mtype,
                            data, &nloop, &eloops, &senses);
    CHECK_STATUS(EG_getTopology);
    if (nloop != 1) {status = -990; goto cleanup;}

    status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                            data, &nedge, &eedges, &senses);
    CHECK_STATUS(EG_getTopology);
    if (nedge != 1) {status = -991; goto cleanup;}

    status = EG_getRange(eedges[0], trange, &periodic);
    CHECK_STATUS(EG_getRange);

    /* extract Face and Loop from ebodys[0] */
    status = EG_getTopology(ebodys[0], &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);
    eface = echilds[0];

    if        (mtype == WIREBODY) {
        eloop = eface;
    } else {
        status = EG_getTopology(eface, &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);
        eloop = echilds[0];
    }

    /* create nxsect Faces/Loops */
    for (i = 0; i < NXSECT(0); i++) {
        tt = trange[0] + (trange[1] - trange[0]) * (double)(i) / (double)(NXSECT(0)-1);
        status = EG_evaluate(eedges[0], &tt, data);
        CHECK_STATUS(EG_evaluate);

        /* this only works for translation or for axes aligned with x-, y-, or z- */
        if (udps[0].arg[2].size != 6) {
            xform[ 0] =  1;          xform[ 1] =  0;          xform[ 2] =  0;          xform[ 3] = data[0]-origin[0];
            xform[ 4] =  0;          xform[ 5] =  1;          xform[ 6] =  0;          xform[ 7] = data[1]-origin[1];
            xform[ 8] =  0;          xform[ 9] =  0;          xform[10] =  1;          xform[11] = data[2]-origin[2];
        } else if (AXIS(0,3) == 1 && AXIS(0,4) == 0 && AXIS(0,5) == 0) {
            theta     = atan2(data[2]-AXIS(0,2), data[1]-AXIS(0,1));
            xform[ 0] =  1;          xform[ 1] =  0;          xform[ 2] =  0;          xform[ 3] = data[0]-origin[0];
            xform[ 4] =  0;          xform[ 5] =  cos(theta); xform[ 6] = -sin(theta); xform[ 7] = data[1]-origin[1];
            xform[ 8] =  0;          xform[ 9] =  sin(theta); xform[10] =  cos(theta); xform[11] = data[2]-origin[2];
        } else if (AXIS(0,3) == 0 && AXIS(0,4) == 1 && AXIS(0,5) == 0) {
            theta     = atan2(data[0]-AXIS(0,0), data[2]-AXIS(0,2));
            xform[ 0] =  cos(theta); xform[ 1] =  0;          xform[ 2] =  sin(theta); xform[ 3] = data[0]-origin[0];
            xform[ 4] =  0;          xform[ 5] =  0;          xform[ 6] =  0;          xform[ 7] = data[1]-origin[1];
            xform[ 8] = -sin(theta); xform[ 9] =  0;          xform[10] =  cos(theta); xform[11] = data[2]-origin[2];
        } else if (AXIS(0,3) == 0 && AXIS(0,4) == 0 && AXIS(0,5) == 1) {
            theta     = atan2(data[1]-AXIS(0,1), data[0]-AXIS(0,0));
            xform[ 0] =  cos(theta); xform[ 1] = -sin(theta); xform[ 2] =  0;          xform[ 3] = data[0]-origin[0];
            xform[ 4] =  sin(theta); xform[ 5] =  cos(theta); xform[ 6] =  0;          xform[ 7] = data[1]-origin[1];
            xform[ 8] =  0;          xform[ 9] =  0;          xform[10] =  1;          xform[11] = data[2]-origin[2];
        } else {
            printf(" udpExecute: axis must be aligned with x-, y-, or z-axis\n");
            status = EGADS_RANGERR;
            goto cleanup;
        }

        status = EG_makeTransform(context, xform, &exform);
        CHECK_STATUS(EG_makeTransform);

        if (i == 0 || i == NXSECT(0)-1) {
            status = EG_copyObject(eface, exform, &exsects[i]);
        } else {
            status = EG_copyObject(eloop, exform, &exsects[i]);
        }
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(exform);
        CHECK_STATUS(EG_deleteObject);
    }

    /* create the blend */
    status = EG_blend(NXSECT(0), exsects, NULL, NULL, ebody);
    CHECK_STATUS(EG_blend);
    if (*ebody == NULL) goto cleanup;   // needed fro splint

#ifdef DEBUG
    printf("*ebody\n");
    (void) ocsmPrintEgo(*ebody);
#endif

    /* add __markFaces__ attribute so that sweep is treated as primitive */
    {
        int idum=0;

        status = EG_attributeAdd(*ebody, "__markFaces__", ATTRINT, 1,
                                 &idum, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);
    }

    /* no output value(s) */

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;
    goto cleanup;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    if (exsects != NULL) EG_free(exsects);

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
