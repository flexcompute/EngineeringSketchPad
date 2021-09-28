/*
 ************************************************************************
 *                                                                      *
 * udfStiffener -- Generate stiffeners from input Model                 *
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

#define NUMUDPARGS       4
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define BEG(    I,IUDP)  ((double *) (udps[IUDP].arg[0].val))[I]
#define END(    I,IUDP)  ((double *) (udps[IUDP].arg[1].val))[I]
#define DEPTH(    IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define ANGLE(    IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"beg",    "end",    "depth",  "angle",  };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,        0,        0,        0,        };
static double argDdefs[NUMUDPARGS] = {0.,       0.,       0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define           PIo180          0.0174532925199432954743717
#define           EPS06           1.0e-06


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

    int     oclass, mtype, nchild, *senses, senses2[4], periodic, nface, iedge;
    double  ubeg, vbeg, uend, vend, trange[2], xform[12], length, dt;
    double  data[18], uv[2], xyz0[18], xyz1[18], xyz2[18], xyz3[18], norm[4];
    ego     context, eref, *ebodys, enodes[4], eedges[8], eloop, eface, etemp[2], exform;
    ego     eshell, *echilds;
    ego     *efaces=NULL, epcurve, epcurve2, esurf, ecurves[4];

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("depth(0) = %f\n", DEPTH(0));
    printf("angle(0) = %f\n", ANGLE(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

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

    /* check that Model was input that contains one Face */
    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);
    if (efaces == NULL) goto cleanup;   // needed for splint

    if (nface != 1) {
        printf(" udpExecute: input Body should have one Face\n");
        status = EGADS_EMPTY;
        goto cleanup;
    }

    /* check arguments */
    if        (udps[0].arg[0].size == 2) {
        ubeg    = BEG(0,0);
        vbeg    = BEG(1,0);
    } else if (udps[0].arg[0].size == 3) {
        xyz0[0] = BEG(0,0);
        xyz0[1] = BEG(1,0);
        xyz0[2] = BEG(2,0);

        status = EG_invEvaluate(efaces[0], xyz0, uv, data);
        CHECK_STATUS(EG_invEvaluate);

        ubeg = uv[0];
        vbeg = uv[1];
    } else {
        printf(" udpExecute: \"beg\" should contain 2 or 3 values\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    if        (udps[0].arg[1].size == 2) {
        uend    = END(0,0);
        vend    = END(1,0);
    } else if (udps[0].arg[1].size == 3) {
        xyz0[0] = END(0,0);
        xyz0[1] = END(1,0);
        xyz0[2] = END(2,0);

        status = EG_invEvaluate(efaces[0], xyz0, uv, data);
        CHECK_STATUS(EG_invEvaluate);

        uend = uv[0];
        vend = uv[1];
    } else {
        printf(" udpExecute: \"end\" should contain 2 or 3 values\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    if (udps[0].arg[2].size != 1) {
        printf(" udpExecute: \"depth\" should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (fabs(DEPTH(0)) < EPS06) {
        printf(" udpExecute: \"depth\" should be non-zero\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    if (udps[0].arg[3].size != 1) {
        printf(" udpExecute: \"angle\" should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;
    } else if (fabs(DEPTH(0)) > 89) {
        printf(" udpExecute: \"angle\" should be less than 89\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("depth(%d) = %f\n", numUdp, DEPTH(numUdp));
    printf("angle(%d) = %f\n", numUdp, ANGLE(numUdp));
#endif

    /* make copy of Face so that Model can be removed */
    status = EG_copyObject(efaces[0], NULL, &eface);
    CHECK_STATUS(EG_copyObject);

    /* create a Bspline Pcurve on the Face (so that its extents are finite) */
    status = EG_getTopology(eface, &esurf, &oclass, &mtype, data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    data[0] = ubeg;
    data[1] = vbeg;
    data[2] = uend - ubeg;
    data[3] = vend - vbeg;

    status = EG_makeGeometry(context, PCURVE, LINE, esurf, NULL, data, &epcurve);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_invEvaluate(epcurve, data, uv, xyz0);
    CHECK_STATUS(EG_invEvaluate);
    trange[0] = uv[0];

    data[0] = uend;
    data[1] = vend;
    status = EG_invEvaluate(epcurve, data, uv, xyz1);
    CHECK_STATUS(EG_invEvaluate);
    trange[1] = uv[0];

    status = EG_convertToBSplineRange(epcurve, trange, &epcurve2);
    CHECK_STATUS(EG_convertToBSplineRange);

    /* use the Pcurve to create an Edge (eedges[0]) on the input Face */
    status = EG_otherCurve(esurf, epcurve2, 0, &(ecurves[0]));
    CHECK_STATUS(EG_otherCurve);

    status = EG_evaluate(ecurves[0], &(trange[0]), xyz0);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, xyz0, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_evaluate(ecurves[0], &(trange[1]), xyz1);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, xyz1, 0, NULL, NULL, &(enodes[1]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, ecurves[0], EDGE, TWONODE, trange, 2, enodes, NULL, &(eedges[0]));
    CHECK_STATUS(EG_makeTopology);

    /* find the "normal" direction at the center of the stiffener */
    uv[0] = (ubeg + uend) / 2;
    uv[1] = (vbeg + vend) / 2;

    status = EG_evaluate(esurf, uv, data);
    CHECK_STATUS(EG_evaluate);

    norm[0]  = data[4] * data[8] - data[5] * data[7];
    norm[1]  = data[5] * data[6] - data[3] * data[8];
    norm[2]  = data[3] * data[7] - data[4] * data[6];
    norm[3]  = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    norm[0] /= norm[3];
    norm[1] /= norm[3];
    norm[2] /= norm[3];

    /* make a copy of the Curve in eedges[0] */
    xform[ 0] = 1;   xform[ 1] = 0;   xform[ 2] = 0;  xform[ 3] = norm[0] * DEPTH(numUdp);
    xform[ 4] = 0;   xform[ 5] = 1;   xform[ 6] = 0;  xform[ 7] = norm[1] * DEPTH(numUdp);
    xform[ 8] = 0;   xform[ 9] = 0;   xform[10] = 1;  xform[11] = norm[2] * DEPTH(numUdp);

    status = EG_makeTransform(context, xform, &exform);
    CHECK_STATUS(EG_makeTransform);

    status = EG_copyObject(ecurves[0], exform, &(ecurves[2]));
    CHECK_STATUS(EG_copyObject);

    status = EG_deleteObject(exform);
    CHECK_STATUS(EG_deleteObject);

    /* make a new Edge that is "shorter" than copied Edge */
    status = EG_getRange(eedges[0], data, &periodic);
    CHECK_STATUS(EG_getRange);
    
    status = EG_arcLength(eedges[0], data[0], data[1], &length);
    CHECK_STATUS(EG_arcLength);

    dt = sqrt(xform[3]*xform[3] + xform[7]*xform[7] + xform[11]*xform[11]) / length
        * (data[1] - data[0]) * tan(ANGLE(numUdp) * PIo180);

    trange[0] = data[0] + dt;
    trange[1] = data[1] - dt;

    status = EG_evaluate(ecurves[2], &(trange[0]), xyz3);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, xyz3, 0, NULL, NULL, &(enodes[3]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_evaluate(ecurves[2], &(trange[1]), xyz2);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, xyz2, 0, NULL, NULL, &(enodes[2]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_invEvaluate(ecurves[2], xyz3, uv, data);
    CHECK_STATUS(EG_invEvaluate);
    trange[0] = uv[0];

    status = EG_invEvaluate(ecurves[2], xyz2, uv, data);
    CHECK_STATUS(EG_invEvaluate);
    trange[1] = uv[0];

    etemp[0] = enodes[3];
    etemp[1] = enodes[2];
    status = EG_makeTopology(context, ecurves[2], EDGE, TWONODE, trange, 2, etemp, NULL, &(eedges[2]));
    CHECK_STATUS(EG_makeTopology);

    /* make the "angled" Edges at the ends of the stiffener */
    data[0] = xyz0[0];
    data[1] = xyz0[1];
    data[2] = xyz0[2];
    data[3] = xyz3[0] - xyz0[0];
    data[4] = xyz3[1] - xyz0[1];
    data[5] = xyz3[2] - xyz0[2];

    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurves[3]));
    CHECK_STATUS(EG_makeGeometry);

    status = EG_invEvaluate(ecurves[3], xyz0, uv, data);
    CHECK_STATUS(EG_invEvaluate);
    trange[0] = uv[0];

    status = EG_invEvaluate(ecurves[3], xyz3, uv, data);
    CHECK_STATUS(EG_invEvaluate);
    trange[1] = uv[0];

    etemp[0] = enodes[0];
    etemp[1] = enodes[3];
    status = EG_makeTopology(context, ecurves[3], EDGE, TWONODE, trange, 2, etemp, NULL, &(eedges[3]));
    CHECK_STATUS(EG_makeTopology);

    data[0] = xyz1[0];
    data[1] = xyz1[1];
    data[2] = xyz1[2];
    data[3] = xyz2[0] - xyz1[0];
    data[4] = xyz2[1] - xyz1[1];
    data[5] = xyz2[2] - xyz1[2];

    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &(ecurves[1]));
    CHECK_STATUS(EG_makeGeometry);

    status = EG_invEvaluate(ecurves[1], xyz1, uv, data);
    CHECK_STATUS(EG_invEvaluate);
    trange[0] = uv[0];

    status = EG_invEvaluate(ecurves[1], xyz2, uv, data);
    CHECK_STATUS(EG_invEvaluate);
    trange[1] = uv[0];

    etemp[0] = enodes[1];
    etemp[1] = enodes[2];
    status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE, trange, 2, etemp, NULL, &(eedges[1]));
    CHECK_STATUS(EG_makeTopology);

    /* make a loop of the Edges */
    senses2[0] = SFORWARD;
    senses2[1] = SFORWARD;
    senses2[2] = SREVERSE;
    senses2[3] = SREVERSE;

    status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 4, eedges, senses2, &eloop);
    CHECK_STATUS(EG_makeTopology);

    /* create a Surface from the Loop (using the overloaded EG_isoCline) */
    status = EG_isoCline(eloop, 0, 0.0, &esurf);
    CHECK_STATUS(EG_isoCline);

    status = EG_deleteObject(eloop);
    CHECK_STATUS(EG_deleteObject);

    /* find the PCurves associated with the Edges */
    for (iedge = 0; iedge < 4; iedge++) {
        status = EG_otherCurve(esurf, eedges[iedge], 0.0, &(eedges[iedge+4]));
        CHECK_STATUS(EG_otherCurve);
    }

    /* make a new Loop and a Face using this Surface */
    status = EG_makeTopology(context, esurf, LOOP, CLOSED, NULL, 4, eedges, senses2, &eloop);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL, 1, &eloop, senses2, &eface);
    CHECK_STATUS(EG_makeTopology);

    /* make a Shell and SheetBody */
    status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL, 1, &eface, NULL, &eshell);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1, &eshell, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);
    if (*ebody == NULL) goto cleanup;   // needed for splint

    status = EG_attributeAdd(*ebody, "__markFaces__", ATTRSTRING, 1, NULL, NULL, "true");
    CHECK_STATUS(EG_attributeAdd);

    /* set the output value (none) */

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
    }

    EG_free(efaces);

    return status;
}


/*
 *****************************************************************************
 *                                                                           *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                           *
 *****************************************************************************
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
