/*
 ************************************************************************
 *                                                                      *
 * udfDroop -- droops leading and/or trailing edge of FaceBody          *
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

#define NUMUDPARGS       6
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define XLE(    IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
#define THETALE(IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define POWLE(  IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define XTE(    IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define THETATE(IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define POWTE(  IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"xle",    "thetale", "powle",  "xte",    "thetate", "powte",  };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL,  ATTRREAL, ATTRREAL, ATTRREAL,  ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {0,        0,         0,        0,        0,         0,        };
static double argDdefs[NUMUDPARGS] = {-100.,    0.,        1.,       100.,     0.,        1.,      };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#ifdef GRAFIC
   #include "grafic.h"
   #define DEBUG
#endif

#ifdef DEBUG
   #include "OpenCSM.h"
#endif

/* prototype for function defined below */

extern int EG_isPlanar(const ego object);

#define  NCP    11


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

    int     oclass, mtype, mtype2, nchild, *senses, *senses2, *idata=NULL;
    int     ncp, nknot, periodic, i, nnode, nedge, nloop, nface, iedge;
    double  bbox[6], scale, mat[12];
    double  data[18], trange[4], *tranges=NULL, tt, *rdata=NULL;
    double  xyz_beg[3], xyz_end[3], frac, xyz_out[3], dyle, dyte;
    char    *message=NULL;
    ego     context, eref, *ebodys, *echilds, eplane;
    ego     etemp1=NULL, etemp2=NULL, exform=NULL, *enodes, *eedges, *eloops=NULL, *efaces=NULL;
    ego     *ebsplines=NULL, *newnodes=NULL, *newedges=NULL, eloop, eface, topRef, prev, next;
    udp_T   *udps = *Udps;

#ifdef GRAFIC
    float   xplot[1000], yplot[1000];
    int     io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
    int     nline=0, nplot=0, ilin[10], isym[10], nper[10];
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("xle(    0) = %f\n", XLE(    0));
    printf("thetale(0) = %f\n", THETALE(0));
    printf("xte(    0) = %f\n", XTE(    0));
    printf("thetate(0) = %f\n", THETATE(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

#ifdef DEBUG
    printf("(input) emodel:\n");
    (void) ocsmPrintEgo(emodel);
#endif

    /* get the context associated with emodel (needed for subsequent
       constructions) */
    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* check that Model was input that contains one FaceBody or SheetBody */
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

    status = EG_getTopology(ebodys[0], &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != BODY || (mtype != FACEBODY && mtype != SHEETBODY)) {
        snprintf(message, 100, "expecting one SheetBody");
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    /* check arguments */
    if        (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "xle should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "thetale should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (THETALE(0) < -89 || THETALE(0) > 89) {
        snprintf(message, 100, "thetale = %f should be between -89 and +89", THETALE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "powle should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (POWLE(0) <= 0) {
        snprintf(message, 100, " powle = %f should be positive", POWLE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        snprintf(message, 100, "xte should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size > 1) {
        snprintf(message, 100, "thetate should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (THETATE(0) < -89 || THETATE(0) > 89) {
        snprintf(message, 100, "thetate = %f should be between -89 and +89", THETATE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        snprintf(message, 100, "powte should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (XTE(0) <= XLE(0)) {
        snprintf(message, 100, "powte = %f should be greater than powle = %f", POWTE(0), POWLE(0));
        status = EGADS_RANGERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("xle(    %d) = %f\n", numUdp, XLE(    numUdp));
    printf("thetale(%d) = %f\n", numUdp, THETALE(numUdp));
    printf("xte(    %d) = %f\n", numUdp, XTE(    numUdp));
    printf("thetate(%d) = %f\n", numUdp, THETATE(numUdp));
#endif

#ifdef GRAFIC
    /* plot pivot locations */
    if (XLE(0) > 0) {
        xplot[nplot] = XLE(0);
        yplot[nplot] = 0.0;
        nplot++;

        ilin[nline] = 0;
        isym[nline] = +GR_PLUS;
        nper[nline]  = 1;
        nline++;
    }

    if (XTE(0) < 1) {
        xplot[nplot] = XTE(0);
        yplot[nplot] = 0.0;
        nplot++;

        ilin[nline] = 0;
        isym[nline] = +GR_PLUS;
        nper[nline]  = 1;
        nline++;
    }
#endif

    /* get the bounding box associated with the input Body */
    status = EG_getBoundingBox(ebodys[0], bbox);
    CHECK_STATUS(EG_getBoundingBox);

    /* if it does not have roughly constant z, stop */
    if (fabs(bbox[5]-bbox[2]) > EPS06) {
        snprintf(message, 100, "Face does not have constant z");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* transform ebodys[0] so that it has unit chord and lleading edge is at x=z=0 */
    scale = 1.0 / (bbox[3] - bbox[0]);

    mat[ 0] = scale;  mat[ 1] = 0;      mat[ 2] = 0;      mat[ 3] = - bbox[0]            * scale;
    mat[ 4] = 0;      mat[ 5] = scale;  mat[ 6] = 0;      mat[ 7] = -(bbox[1]+bbox[4])/2 * scale;
    mat[ 8] = 0;      mat[ 9] = 0;      mat[10] = scale;  mat[11] = - bbox[2]            * scale;

    status = EG_makeTransform(context, mat, &exform);
    CHECK_STATUS(EG_makeTransform);

    SPLINT_CHECK_FOR_NULL(exform);

    status = EG_copyObject(ebodys[0], exform, &etemp1);
    CHECK_STATUS(EG_copyObject);

    SPLINT_CHECK_FOR_NULL(etemp1);

    status = EG_deleteObject(exform);
    CHECK_STATUS(EG_deleteObject);

    exform = NULL;

    /* get the Loop associated with the (modified) input Body */
    status = EG_getBodyTopos(etemp1, NULL, LOOP,
                             &nloop, &eloops);
    CHECK_STATUS(EG_getBodyTopos);
    if (eloops == NULL) goto cleanup;   // needed for splint

    if (nloop != 1) {
        snprintf(message, 100, "Body has %d Loops (expecting only 1)", nloop);
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* get the Face associated with the (modified) input Body */
    status = EG_getBodyTopos(etemp1, NULL, FACE,
                             &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);
    if (efaces == NULL) goto cleanup;   // needed for splint

    if (nface != 1) {
        snprintf(message, 100, "Body has %d Faces (expecting only 1)", nface);
        status = EGADS_RANGERR;
        goto cleanup;
    }

    status = EG_getTopology(efaces[0], &eref, &oclass, &mtype2,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    /* copy needs to be made so that emodel can be removed by OpenCSM */
    status = EG_copyObject(eref, NULL, &eplane);
    CHECK_STATUS(EG_copyObject);

    /* get the BSplines associated with the Edges in the Loop */
    status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                            data, &nedge, &eedges, &senses);
    CHECK_STATUS(EG_getTopology);

    MALLOC(ebsplines, ego,      nedge   );
    MALLOC(tranges,   double, 2*nedge   );
    MALLOC(newnodes,  ego,     (nedge+1));
    MALLOC(newedges,  ego,      nedge   );

    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getTopology(eedges[iedge], &ebsplines[iedge], &oclass, &mtype,
                                data, &nnode, &enodes, &senses2);
        CHECK_STATUS(EG_getTopology);

        status = EG_getInfo(ebsplines[iedge], &oclass, &mtype, &topRef, &prev, &next);
        CHECK_STATUS(EG_getInfo);

        tranges[2*iedge  ] = data[0];
        tranges[2*iedge+1] = data[1];

        status = EG_getGeometry(ebsplines[iedge], &oclass, &mtype,
                                &eref, &idata, &rdata);
        CHECK_STATUS(EG_getGeometry);

        if (mtype == LINE) {
            EG_free(idata);     idata = NULL;
            EG_free(rdata);     rdata = NULL;

            MALLOC(idata, int,     4       );
            MALLOC(rdata, double, (4*NCP+2));

            idata[0] = 0;
            idata[1] = 1;
            idata[2] = NCP;
            idata[3] = NCP + 2;

            rdata[0] = 0;
            for (i = 0; i < NCP; i++) {
                rdata[i+1] = (double)(i) / (double)(NCP-1);
            }
            rdata[NCP+1] = 1;

            status = EG_getTopology(enodes[0], &eref, &oclass, &mtype,
                                    xyz_beg, &nchild, &echilds, &senses2);
            CHECK_STATUS(EG_getTopology);

            status = EG_getTopology(enodes[1], &eref, &oclass, &mtype,
                                    xyz_end, &nchild, &echilds, &senses2);
            CHECK_STATUS(EG_getTopology);

            for (i = 0; i < NCP; i++) {
                frac = (double)(i) / (double)(NCP-1);
                rdata[NCP+3*i+2] = (1-frac) * xyz_beg[0] + frac * xyz_end[0];
                rdata[NCP+3*i+3] = (1-frac) * xyz_beg[1] + frac * xyz_end[1];
                rdata[NCP+3*i+4] = (1-frac) * xyz_beg[2] + frac * xyz_end[2];
            }

            status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                                     idata, rdata, &ebsplines[iedge]);
            CHECK_STATUS(EG_makeGeometry);

            tranges[2*iedge  ] = 0;
            tranges[2*iedge+1] = 1;
        } else if (mtype != BSPLINE) {
            status = EG_convertToBSpline(eedges[iedge], &ebsplines[iedge]);
            CHECK_STATUS(EG_convertToBSpline);

            status = EG_getRange(ebsplines[iedge], trange, &periodic);
            CHECK_STATUS(EG_getRange);

            tranges[2*iedge  ] = trange[0];
            tranges[2*iedge+1] = trange[1];
        }

        EG_free(idata);     idata = NULL;
        EG_free(rdata);     rdata = NULL;
    }

#ifdef GRAFIC
    /* plot the initial control points */
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getGeometry(ebsplines[iedge], &oclass, &mtype,
                                &eref, &idata, &rdata);
        CHECK_STATUS(EG_getGeometry);

        ncp     = idata[2];
        nknot   = idata[3];

        for (i = 0; i < ncp; i++) {
            xplot[nplot] = rdata[nknot+3*i  ];
            yplot[nplot] = rdata[nknot+3*i+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_CIRCLE;
        nper[nline] = ncp;
        nline++;

        EG_free(idata);     idata = NULL;
        EG_free(rdata);     rdata = NULL;
    }
#endif

    /* compute the modified y at the leading and trailing edge */
    if (XLE(0) > 0) {
        dyle = - (XLE(0)-0) * tan(THETALE(0) * PI/180);
    } else {
        dyle = 0;
    }

    if (XTE(0) < 1) {
        dyte = + (1-XTE(0)) * tan(THETATE(0) * PI/180);
    } else {
        dyte = 0;
    }

    /* modify the control points forward of XLE and aft of XTE and
       create  the new Bspline curve */
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getGeometry(ebsplines[iedge], &oclass, &mtype,
                                &eref, &idata, &rdata);
        CHECK_STATUS(EG_getGeometry);
        if (idata == NULL) goto cleanup;     // needed for splint
        if (rdata == NULL) goto cleanup;     // needed for splint

        ncp     = idata[2];
        nknot   = idata[3];

        for (i = 0; i < ncp; i++) {
            if        (rdata[nknot+3*i] < XLE(0)) {
                frac = (XLE(0) - rdata[nknot+3*i]) / (XLE(0) - 0);
                rdata[nknot+3*i+1] += dyle * pow(frac, POWLE(0));

            } else if (rdata[nknot+3*i] > XTE(0)) {
                frac = (rdata[nknot+3*i] - XTE(0)) / (1 - XTE(0));
                rdata[nknot+3*i+1] += dyte * pow(frac, POWTE(0));
            }
        }

        if (XLE(0) > 0 || XTE(0) < 1) {
            status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                                     idata, rdata, &(ebsplines[iedge]));
            CHECK_STATUS(EG_makeGeometry);
        }

        EG_free(idata);     idata = NULL;
        EG_free(rdata);     rdata = NULL;
    }

#ifdef GRAFIC
    /* plot the modified control points */
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_getGeometry(ebsplines[iedge], &oclass, &mtype,
                                &eref, &idata, &rdata);
        CHECK_STATUS(EG_getGeometry);

        ncp     = idata[2];
        nknot   = idata[3];

        for (i = 0; i < ncp; i++) {
            xplot[nplot] = rdata[nknot+3*i  ];
            yplot[nplot] = rdata[nknot+3*i+1];
            nplot++;
        }
        ilin[nline] = -GR_DASHED;
        isym[nline] = +GR_SQUARE;
        nper[nline] = ncp;
        nline++;

        EG_free(idata);     idata = NULL;
        EG_free(rdata);     rdata = NULL;
    }
#endif

#ifdef GRAFIC
    grinit_(&io_kbd, &io_scr, "udfDroop", strlen("udfDroop"));
    grline_(ilin, isym, &nline,                "~x~y~airfoil",
            &indgr, xplot, yplot, nper, strlen("~x~y~airfoil"));
#endif

    /* make the Nodes for the new Body */
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_evaluate(ebsplines[iedge], &(tranges[2*iedge]), data);
        CHECK_STATUS(EG_evaluate);

        status = EG_invEvaluate(ebsplines[iedge], data, &tt, xyz_out);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_makeTopology(context, NULL, NODE, 0,
                                 xyz_out, 0, NULL,NULL, &(newnodes[iedge]));
        CHECK_STATUS(EG_makeTopology);
    }
    newnodes[nedge] = newnodes[0];

    /* make the Edges for the new Body */
    for (iedge = 0; iedge < nedge; iedge++) {
        status = EG_makeTopology(context, ebsplines[iedge], EDGE, TWONODE,
                                 &(tranges[2*iedge]), 2, &(newnodes[iedge]),
                                 NULL, &(newedges[iedge]));
        CHECK_STATUS(EG_makeTopology);
    }

    /* make the Face and copy attributes from original Face */
    status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                             NULL, nedge, newedges, senses, &eloop);
    CHECK_STATUS(EG_makeTopology);

    senses[0] = SFORWARD;
    status = EG_makeTopology(context, eplane, FACE, mtype2,
                             NULL, 1, &eloop, senses, &eface);
    CHECK_STATUS(EG_makeTopology);

    status = EG_attributeDup(efaces[0], eface);
    CHECK_STATUS(EG_attributeDup);

    /* make the Body and copy attributes from original Body */
    status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                             NULL, 1, &eface, NULL, &etemp2);
    CHECK_STATUS(EG_makeTopology);
    if (etemp2 == NULL) goto cleanup;   // needed for splint

#ifdef DEBUG
    printf("(output) *ebody:\n");
    (void) ocsmPrintEgo(*ebody);
#endif

    /* rescale the Body */
    mat[ 0] = 1/scale;  mat[ 1] = 0;        mat[ 2] = 0;        mat[ 3] = + bbox[0];
    mat[ 4] = 0;        mat[ 5] = 1/scale;  mat[ 6] = 0;        mat[ 7] = +(bbox[1]+bbox[4])/2;;
    mat[ 8] = 0;        mat[ 9] = 0;        mat[10] = 1/scale;  mat[11] = + bbox[2];

    status = EG_makeTransform(context, mat, &exform);
    CHECK_STATUS(EG_makeTransform);

    SPLINT_CHECK_FOR_NULL(exform);

    status = EG_copyObject(etemp2, exform, ebody);
    CHECK_STATUS(EG_copyObject);

    SPLINT_CHECK_FOR_NULL(*ebody);

    status = EG_deleteObject(exform);
    CHECK_STATUS(EG_deleteObject);

    exform = NULL;

    status = EG_attributeDup(etemp1, *ebody);
    CHECK_STATUS(EG_attributeDup);

    /* there are no output values to set */

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (etemp1    != NULL) EG_deleteObject(etemp1);
    if (etemp2    != NULL) EG_deleteObject(etemp2);
    if (exform    != NULL) EG_deleteObject(exform);

    if (eloops    != NULL) EG_free(eloops   );
    if (efaces    != NULL) EG_free(efaces   );
    if (rdata     != NULL) EG_free(rdata    );
    if (idata     != NULL) EG_free(idata    );
    if (newnodes  != NULL) EG_free(newnodes );
    if (newedges  != NULL) EG_free(newedges );
    if (tranges   != NULL) EG_free(tranges  );
    if (ebsplines != NULL) EG_free(ebsplines);

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
