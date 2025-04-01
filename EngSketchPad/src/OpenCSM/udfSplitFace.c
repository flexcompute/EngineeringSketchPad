/*
 ************************************************************************
 *                                                                      *
 * udfSplitFace -- splits a given Face                                  *
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

#define NUMUDPARGS       3
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define IFACE(IUDP  )  ((int    *) (udps[IUDP].arg[0].val))[0]
#define BEG(  IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]
#define END(  IUDP,I)  ((double *) (udps[IUDP].arg[2].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"iface", "beg",    "end",    };
static int    argTypes[NUMUDPARGS] = {ATTRINT, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,       0,        0,        };
static double argDdefs[NUMUDPARGS] = {0.,      0.,       0.,       };

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

    int     oclass, mtype, nchild, *senses, senses2[4], match;
    int     nnode, nedge0, iedge0, nedge1, iedge1, nface, iface;
    double  ubeg, vbeg, uend, vend, trange[2];
    double  data[18], uv[2], xyz0[18], xyz1[18];
    char    *message=NULL;
    ego     context, eref, *ebodys=NULL, enodes2[2], eedge, eloop, *echilds;
    ego     *enodes=NULL, *eedges0=NULL, *eedges1=NULL, *efaces=NULL;
    ego     epcurve, epcurve2, esurf, ecurve, enewmodel, ewirebody;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("iface(0) = %d\n", IFACE(0));
    if        (udps[0].arg[1].size == 1) {
        printf("beg(  0) = %f\n",       BEG(0,0));
    } else if (udps[0].arg[1].size == 2) {
        printf("beg(  0) = %f %f\n",    BEG(0,0), BEG(0,1));
    } else {
        printf("beg(  0) = %f %f %f\n", BEG(0,0), BEG(0,1), BEG(0,2));
    }
    if        (udps[0].arg[2].size == 1) {
        printf("end(  0) = %f\n",       END(0,0));
    } else if (udps[0].arg[2].size == 2) {
        printf("end(  0) = %f %f\n",    END(0,0), END(0,1));
    } else {
        printf("end(  0) = %f %f %f\n", END(0,0), END(0,1), END(0,2));
    }
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 100, "expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    SPLINT_CHECK_FOR_NULL(ebodys);

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* get the number of Faces in the Body */
    status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(efaces);

    /* check arguments */
    if (IFACE(0) < 1 || IFACE(0) > nface) {
        snprintf(message, 100, "\"iface\" should be between 1 and %d", nface);
        status = EGADS_RANGERR;
        goto cleanup;
    }

    status = EG_getTopology(efaces[IFACE(0)-1], &esurf, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if        (udps[0].arg[1].size == 1) {
        status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(enodes);

        if (NINT(BEG(0,0)) < 0 || NINT(BEG(0,0)) >= nnode) {
            snprintf(message, 100, "inode=%d is out of range", NINT(BEG(0,0)));
            status = EGADS_RANGERR;
            goto cleanup;
        }

        status = EG_evaluate(enodes[NINT(BEG(0,0))-1], 0, xyz0);
        CHECK_STATUS(EG_evaluate);

        status = EG_invEvaluate(esurf, xyz0, uv, data);
        CHECK_STATUS(EG_invEvaluate);

        ubeg = uv[0];
        vbeg = uv[1];

        EG_free(enodes);   enodes = NULL;
    } else if (udps[0].arg[1].size == 2) {
        ubeg    = BEG(0,0);
        vbeg    = BEG(0,1);
    } else if (udps[0].arg[1].size == 3) {
        xyz0[0] = BEG(0,0);
        xyz0[1] = BEG(0,1);
        xyz0[2] = BEG(0,2);

        status = EG_invEvaluate(esurf, xyz0, uv, data);
        CHECK_STATUS(EG_invEvaluate);

        ubeg = uv[0];
        vbeg = uv[1];
    } else {
        snprintf(message, 100, "\"beg\" should contain 1, 2, or 3 values\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    if        (udps[0].arg[2].size == 1) {
        status = EG_getBodyTopos(ebodys[0], NULL, NODE, &nnode, &enodes);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(enodes);

        if (NINT(END(0,0)) < 0 || NINT(END(0,0)) >= nnode) {
            snprintf(message, 100, "inode=%d is out of range", NINT(END(0,0)));
            status = EGADS_RANGERR;
            goto cleanup;
        }

        status = EG_evaluate(enodes[NINT(END(0,0))-1], 0, xyz0);
        CHECK_STATUS(EG_evaluate);

        status = EG_invEvaluate(esurf, xyz0, uv, data);
        CHECK_STATUS(EG_invEvaluate);

        uend = uv[0];
        vend = uv[1];

        EG_free(enodes);   enodes = NULL;
    } else if (udps[0].arg[2].size == 2) {
        uend    = END(0,0);
        vend    = END(0,1);
    } else if (udps[0].arg[2].size == 3) {
        xyz0[0] = END(0,0);
        xyz0[1] = END(0,1);
        xyz0[2] = END(0,2);

        status = EG_invEvaluate(esurf, xyz0, uv, data);
        CHECK_STATUS(EG_invEvaluate);

        uend = uv[0];
        vend = uv[1];
    } else {
        snprintf(message, 100, "\"end\" should contain 1, 2, or 3 values\n");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* make sure beg and end are within the Face */
    uv[0] = ubeg;
    uv[1] = vbeg;
    status = EG_inFace(efaces[IFACE(0)-1], uv);
    if (status != EGADS_SUCCESS) {
        snprintf(message, 100, "beg (%f,%f) is not in Face %d", uv[0], uv[1], IFACE(0));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    uv[0] = uend;
    uv[1] = vend;
    status = EG_inFace(efaces[IFACE(0)-1], uv);
    if (status != EGADS_SUCCESS) {
        snprintf(message, 100, "end (%f,%f) is not in Face %d", uv[0], uv[1], IFACE(0));
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("iface(%d) = %d\n", numUdp, IFACE(numUdp));
    if        (udps[0].arg[1].size == 1) {
        printf("beg(%d) = %f\n",       numUdp, BEG(numUdp,0));
    } else if (udps[0].arg[1].size == 2) {
        printf("beg(%d) = %f %f\n",    numUdp, BEG(numUdp,0), BEG(numUdp,1));
    } else {
        printf("beg(%d) = %f %f %f\n", numUdp, BEG(numUdp,0), BEG(numUdp,1), BEG(numUdp,2));
    }
    if        (udps[0].arg[2].size == 1) {
        printf("end(%d) = %f\n",       numUdp, END(numUdp,0));
    } else if (udps[0].arg[2].size == 2) {
        printf("end(%d) = %f %f\n",    numUdp, END(numUdp,0), END(numUdp,1));
    } else {
        printf("end(%d) = %f %f %f\n", numUdp, END(numUdp,0), END(numUdp,1), END(numUdp,2));
    }
#endif

    iface = IFACE(0);

    /* create a Bspline Pcurve on the Face (so that its extents are finite) */
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

    /* use the Pcurve to create an Edge (eedge) on the input Face */
    status = EG_otherCurve(esurf, epcurve2, 0, &ecurve);
    CHECK_STATUS(EG_otherCurve);

    status = EG_evaluate(ecurve, &(trange[0]), xyz0);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, xyz0, 0, NULL, NULL, &(enodes2[0]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_evaluate(ecurve, &(trange[1]), xyz1);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, xyz1, 0, NULL, NULL, &(enodes2[1]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE, trange, 2, enodes2, NULL, &eedge);
    CHECK_STATUS(EG_makeTopology);

    senses2[0] = SFORWARD;
    status = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, 1, &eedge, senses2, &eloop);
    CHECK_STATUS(EG_makeTopology);

    status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, &ewirebody);
    CHECK_STATUS(EG_makeTopology);

    /* scribe the input Body with the new WireBody (Edge) */
    status = EG_generalBoolean(ebodys[0], ewirebody, SPLITTER, 0, &enewmodel);
    CHECK_STATUS(EG_generalBoolean);

    status = EG_getTopology(enewmodel, &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (nchild != 1) {
        snprintf(message, 100, "splitter created %d Bodys (not 1)", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    *ebody = echilds[0];

    /* put an attribute on the Edges that came from the scribe */
    status = EG_getBodyTopos(ebodys[0], NULL, EDGE, &nedge0, &eedges0);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge1, &eedges1);
    CHECK_STATUS(EG_getBodyTopos);

    SPLINT_CHECK_FOR_NULL(eedges0);
    SPLINT_CHECK_FOR_NULL(eedges1);

    for (iedge1 = 0; iedge1 < nedge1; iedge1++) {
        match = 0;

        for (iedge0 = 0; iedge0 < nedge0; iedge0++) {
            if (EG_isSame(eedges0[iedge0], eedges1[iedge1]) == EGADS_SUCCESS) {
                match = 1;
                break;
            }
        }

        if (match == 0) {
            status = EG_attributeAdd(eedges1[iedge1], "__splitFace__",
                                     ATTRINT, 1, &iface, NULL, NULL);
            CHECK_STATUS(EG_attributeAdd);
        }
    }

    EG_free(eedges0);   eedges0 = NULL;
    EG_free(eedges1);   eedges1 = NULL;

    /* set the output value (none) */

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (*ebody == NULL && ebodys != NULL) {
        (void) EG_copyObject(ebodys[0], NULL, ebody);
    }

    if (strlen(message) > 0) {
        *string = message;
        printf("%s\n", message);
    } else if (status != EGADS_SUCCESS) {
        FREE(message);
        *string = udpErrorStr(status);
    } else {
        FREE(message);
    }

    if (enodes  != NULL) EG_free(enodes );
    if (eedges0 != NULL) EG_free(eedges0);
    if (eedges1 != NULL) EG_free(eedges1);
    if (efaces  != NULL) EG_free(efaces );

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
