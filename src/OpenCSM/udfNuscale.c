/*
 ************************************************************************
 *                                                                      *
 * udfNuscale -- performs non-uniform scaling on a B-spline Body        *
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

#define NUMUDPARGS       6
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define XSCALE(IUDP)  ((double *) (udps[IUDP].arg[0].val))[0]
#define YSCALE(IUDP)  ((double *) (udps[IUDP].arg[1].val))[0]
#define ZSCALE(IUDP)  ((double *) (udps[IUDP].arg[2].val))[0]
#define XCENT( IUDP)  ((double *) (udps[IUDP].arg[3].val))[0]
#define YCENT( IUDP)  ((double *) (udps[IUDP].arg[4].val))[0]
#define ZCENT( IUDP)  ((double *) (udps[IUDP].arg[5].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"xscale", "yscale", "zscale", "xcent",  "ycent",  "zcent",  };
static int    argTypes[NUMUDPARGS] = {ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL, };
static int    argIdefs[NUMUDPARGS] = {1,        1,        1,        0,        0,        0,        };
static double argDdefs[NUMUDPARGS] = {1.,       1.,       1.,       0.,       0.,       0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define  EPS06   1.0e-6

/* prototype for function defined below */
static int convertToBSpline(ego inbody, double matrix[], ego *outbody);


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

    int     oclass, mtype, nchild, *senses;
    double  xyz[18], mat[12];
    ego     *ebodys, eref;

    char    *message=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("xscale(0) = %f\n", XSCALE(0));
    printf("yscale(0) = %f\n", YSCALE(0));
    printf("zscale(0) = %f\n", ZSCALE(0));
    printf("xcent( 0) = %f\n", XCENT( 0));
    printf("ycent( 0) = %f\n", YCENT( 0));
    printf("zcent( 0) = %f\n", ZCENT( 0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            xyz, &nchild, &ebodys, &senses);
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

    /* check arguments */
    if        (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "xscale should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "yscale should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "zscale should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        snprintf(message, 100, "xcent should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[4].size > 1) {
        snprintf(message, 100, "ycent should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size > 1) {
        snprintf(message, 100, "zcent should be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (fabs(XSCALE(0)) < EPS06 && fabs(YSCALE(0)) < EPS06) {
        snprintf(message, 100, "xscale and yscale cannot both be 0");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (fabs(YSCALE(0)) < EPS06 && fabs(ZSCALE(0)) < EPS06) {
        snprintf(message, 100, "yscale and zscale cannot both be 0");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (fabs(ZSCALE(0)) < EPS06 && fabs(ZSCALE(0)) < EPS06) {
        snprintf(message, 100, "zscale and xscale cannot both be 0");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("xscale(%d) = %f\n", numUdp, XSCALE(numUdp));
    printf("yscale(%d) = %f\n", numUdp, YSCALE(numUdp));
    printf("zscale(%d) = %f\n", numUdp, ZSCALE(numUdp));
    printf("xcent( %d) = %f\n", numUdp, XCENT( numUdp));
    printf("ycent( %d) = %f\n", numUdp, YCENT( numUdp));
    printf("zcent( %d) = %f\n", numUdp, ZCENT( numUdp));
#endif

    /* set up transformation matrix */
    mat[ 0] = XSCALE(0);   mat[ 1] = 0;           mat[ 2] = 0;           mat[ 3] = XCENT(0) * (1-XSCALE(0));
    mat[ 4] = 0;           mat[ 5] = YSCALE(0);   mat[ 6] = 0;           mat[ 7] = YCENT(0) * (1-YSCALE(0));
    mat[ 8] = 0;           mat[ 9] = 0;           mat[10] = ZSCALE(0);   mat[11] = ZCENT(0) * (1-ZSCALE(0));

    status = convertToBSpline(ebodys[0], mat, ebody);
    CHECK_STATUS(convertToBSpline);

    /* set the output value */

    /* the Body is returned */
    udps[numUdp].ebody = *ebody;

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
 *   convertToBSpline - conert Body to BSpline with transformation      *
 *                                                                      *
 ************************************************************************
 */

static int
convertToBSpline(ego    inbody,         /* (in)  input Body */
                 double mat[],          /* (in)  transformation matrix */
                 ego    *ebody)         /* (out) pointer to output Body */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, ichild, nchild, *senses;
    int     inode, nnode, iedge, nedge, iface, nface, iloop, nloop, ishell, nshell, *ivals, i, sense;
    int     oclass1, mtype1, nchild1, *senses1, jedge, periodic;
    double  xyz[18], xyz2[18], *rvals, trange[2], xyzClose[3];
    double  uv[18], range[4], tt[2], swap;
    ego     context, eref, *echilds, *enodes_orig=NULL, *eedges_orig=NULL, *efaces_orig=NULL, *eloops_orig=NULL, *eshells_orig=NULL;
    ego     *enodes=NULL, *ecurves=NULL, *eedges=NULL, *esurfaces=NULL, *efaces=NULL, *eshells=NULL;
    ego     *eloops2=NULL, *eloops=NULL, *efoo=NULL;
    ego     ecurve, esurface, etemp[2];
            ego *echilds1;

    ROUTINE(convertToBspline);

    /* --------------------------------------------------------------- */

    status = EG_getContext(inbody, &context);
    CHECK_STATUS(EG_getContext);

#ifdef DEBUG
    printf("convertToBSpline(mat=%10.5f %10.5f %10.5f %10.5f\n", mat[ 0], mat[ 1], mat[ 2], mat[ 3]);
    printf("                     %10.5f %10.5f %10.5f %10.5f\n", mat[ 4], mat[ 5], mat[ 6], mat[ 7]);
    printf("                     %10.5f %10.5f %10.5f %10.5f\n", mat[ 8], mat[ 9], mat[10], mat[11]);

    (void) ocsmPrintEgo(inbody);
#endif
    
    /* get the Nodes, Edges, and Faces associated with the input Body */
    status = EG_getBodyTopos(inbody, NULL, NODE, &nnode, &enodes_orig);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(inbody, NULL, EDGE, &nedge, &eedges_orig);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(inbody, NULL, FACE, &nface, &efaces_orig);
    CHECK_STATUS(EG_getBodyTopos);

    status = EG_getBodyTopos(inbody, NULL, SHELL, &nshell, &eshells_orig);
    CHECK_STATUS(EG_getBodyTopos);

#ifdef DEBUG
    printf("nnode=%d, nedge=%d, nface=%d, nshell=%d\n", nnode, nedge, nface, nshell);
#endif

    /* get storage for the new Nodes, Curves, Edges, Surfaces, Faces, and Shells */
    MALLOC(enodes,    ego, nnode );
    MALLOC(ecurves,   ego, nedge );
    MALLOC(eedges,    ego, nedge );
    MALLOC(esurfaces, ego, nface );
    MALLOC(efaces,    ego, nface );
    MALLOC(eshells,   ego, nshell);

    /* make new Nodes by transforming the original Nodes */
    for (inode = 0; inode < nnode; inode++) {
        SPLINT_CHECK_FOR_NULL(enodes_orig );

        status = EG_getTopology(enodes_orig[inode], &eref, &oclass, &mtype,
                                xyz2, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        xyz[0] = mat[ 0] * xyz2[0] + mat[ 1] * xyz2[1] + mat[ 2] * xyz2[2] + mat[ 3];
        xyz[1] = mat[ 4] * xyz2[0] + mat[ 5] * xyz2[1] + mat[ 6] * xyz2[2] + mat[ 7];
        xyz[2] = mat[ 8] * xyz2[0] + mat[ 9] * xyz2[1] + mat[10] * xyz2[2] + mat[11];

        status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, &enodes[inode]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeDup(enodes_orig[inode], enodes[inode]);
        CHECK_STATUS(EG_attributeDup);
    }

    /* special processing for NodeBody */
    if (nnode == 1) {
        MALLOC(eloops, ego, 1);

        trange[0] = 0;
        trange[1] = 1;
        sense     = SFORWARD;

        status = EG_makeTopology(context, NULL, EDGE, DEGENERATE, trange, 1, enodes, &sense, &eedges[0]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 1, eedges, &sense, eloops);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, eloops, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        FREE(eloops);

        goto cleanup;
    }

    /* make new Edges by transforming the original Edges */
    for (iedge = 0; iedge < nedge; iedge++) {
        SPLINT_CHECK_FOR_NULL(eedges_orig );

        status = EG_getTopology(eedges_orig[iedge], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        if (mtype == DEGENERATE) {
            eedges[iedge] = NULL;
            continue;
        }

        status = EG_convertToBSpline(eedges_orig[iedge], &ecurve);
        CHECK_STATUS(EG_convertToBSpline);

        status = EG_getGeometry(ecurve, &oclass, &mtype, &eref, &ivals, &rvals);
        CHECK_STATUS(EG_getGeometry);

        for (i = ivals[3]; i < ivals[3]+3*ivals[2]; i+=3) {
            xyz2[0] = rvals[i  ];
            xyz2[1] = rvals[i+1];
            xyz2[2] = rvals[i+2];

            rvals[i  ] = mat[ 0] * xyz2[0] + mat[ 1] * xyz2[1] + mat[ 2] * xyz2[2] + mat[ 3];
            rvals[i+1] = mat[ 4] * xyz2[0] + mat[ 5] * xyz2[1] + mat[ 6] * xyz2[2] + mat[ 7];
            rvals[i+2] = mat[ 8] * xyz2[0] + mat[ 9] * xyz2[1] + mat[10] * xyz2[2] + mat[11];
        }

        status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, ivals, rvals, &ecurves[iedge]);
        CHECK_STATUS(EG_makeGeometry);

        EG_free(ivals);
        EG_free(rvals);

        status = EG_getTopology(eedges_orig[iedge], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        etemp[0] = enodes[EG_indexBodyTopo(inbody, echilds[       0])-1];
        etemp[1] = enodes[EG_indexBodyTopo(inbody, echilds[nchild-1])-1];

        status = EG_getTopology(etemp[0], &eref, &oclass1, &mtype1,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_invEvaluate(ecurves[iedge], xyz, &trange[0], xyzClose);
        CHECK_STATUS(EG_invEvaluate);

        status = EG_getTopology(etemp[1], &eref, &oclass1, &mtype1,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_invEvaluate(ecurves[iedge], xyz, &trange[1], xyzClose);
        CHECK_STATUS(EG_invEvaluate);

        if (trange[0] > trange[1]) {
            trange[0] -= TWOPI;
        }

        status = EG_makeTopology(context, ecurves[iedge], oclass, mtype, trange,
                                 2, etemp, NULL, &eedges[iedge]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_attributeDup(eedges_orig[iedge], eedges[iedge]);
        CHECK_STATUS(EG_attributeDup);
    }

    /* special processing for WireBody */
    if (nface == 0) {
        status = EG_getBodyTopos(inbody, NULL, LOOP, &nloop, &eloops_orig);
        CHECK_STATUS(EG_getBodyTopos);

        SPLINT_CHECK_FOR_NULL(eloops_orig);

        MALLOC(eloops, ego, nloop);

        status = EG_getTopology(eloops_orig[0], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        MALLOC(efoo, ego, nchild);

        for (iedge = 0; iedge < nchild; iedge++) {
            efoo[iedge] = eedges[EG_indexBodyTopo(inbody, echilds[iedge])-1];
        }

        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, efoo, senses, &eloops[0]);
        CHECK_STATUS(EG_makeTopology);

        FREE(efoo  );

        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, eloops, senses, ebody);
        CHECK_STATUS(EG_makeTopology);

        FREE(eloops);

        goto cleanup;
    }

    /* make new Faces by transforming the original Faces */
    for (iface = 0; iface < nface; iface++) {
        SPLINT_CHECK_FOR_NULL(efaces_orig );

        status = EG_getTopology(efaces_orig[iface], &eref, &oclass, &mtype,
                                xyz, &nloop, &eloops2, &senses);
        CHECK_STATUS(EG_getTopology);

        SPLINT_CHECK_FOR_NULL(eloops2);

        status = EG_convertToBSpline(efaces_orig[iface], &esurface);
        CHECK_STATUS(EG_convertToBSpline);

        status = EG_getGeometry(esurface, &oclass, &mtype, &eref, &ivals, &rvals);
        CHECK_STATUS(EG_getGeometry);

        for (i = ivals[3]+ivals[6]; i < ivals[3]+ivals[6]+3*ivals[2]*ivals[5]; i+=3) {
            xyz2[0] = rvals[i  ];
            xyz2[1] = rvals[i+1];
            xyz2[2] = rvals[i+2];

            rvals[i  ] = mat[ 0] * xyz2[0] + mat[ 1] * xyz2[1] + mat[ 2] * xyz2[2] + mat[ 3];
            rvals[i+1] = mat[ 4] * xyz2[0] + mat[ 5] * xyz2[1] + mat[ 6] * xyz2[2] + mat[ 7];
            rvals[i+2] = mat[ 8] * xyz2[0] + mat[ 9] * xyz2[1] + mat[10] * xyz2[2] + mat[11];
        }

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, ivals, rvals, &esurfaces[iface]);
        CHECK_STATUS(EG_makeGeometry);

        EG_free(ivals);
        EG_free(rvals);

        /* process the Loops */
        MALLOC(eloops, ego, nloop);

        for (iloop = 0; iloop < nloop; iloop++) {
            status = EG_getTopology(eloops2[iloop], &eref, &oclass, &mtype,
                                    xyz, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            MALLOC(efoo, ego, 2*nchild);

            /* process each Edge in each Loop (but store NULLs for DEGENERATE Edges */
            for (iedge = 0; iedge < nchild; iedge++) {
                status = EG_getTopology(echilds[iedge], &eref, &oclass1, &mtype1,
                                        xyz, &nchild1, &echilds1, &senses1);
                CHECK_STATUS(EG_getTopology);

                if (mtype1 != DEGENERATE) {
                    efoo[iedge] = eedges[EG_indexBodyTopo(inbody, echilds[iedge])-1];

                    status = EG_otherCurve(esurfaces[iface], efoo[iedge], 0, &efoo[iedge+nchild]);
                    CHECK_STATUS(EG_otherCurve);
                } else {
                    efoo[iedge       ] = NULL;
                    efoo[iedge+nchild] = NULL;
                }
            }

            /* if any of the Edges were NULL (DEGENERATE), proces them now */
            for (iedge = 0; iedge < nchild; iedge++) {
                if (efoo[iedge] != NULL) continue;

                /* Edge before iedge */
                if (iedge > 0) {
                    jedge = iedge - 1;
                } else {
                    jedge = nchild - 1;
                }

                status = EG_getTopology(efoo[jedge], &eref, &oclass1, &mtype1,
                                        xyz, &nchild1, &echilds1, &senses1);
                CHECK_STATUS(EG_getTopology);

                status = EG_getRange(efoo[jedge], tt, &periodic);
                CHECK_STATUS(EG_getRange);

                if (senses[jedge] == SFORWARD) {
                    etemp[0] = echilds1[1];
                    status = EG_evaluate(efoo[jedge+nchild], &tt[1], &uv[0]);
                    CHECK_STATUS(EG_evaluate);
                } else {
                    etemp[0] = echilds1[0];
                    status = EG_evaluate(efoo[jedge+nchild], &tt[0], &uv[0]);
                    CHECK_STATUS(EG_evaluate);
                }

                /* Edge after iedge */
                if (iedge < nchild-1) {
                    jedge = iedge + 1;
                } else {
                    jedge = 0;
                }

                status = EG_getRange(efoo[jedge], tt, &periodic);
                CHECK_STATUS(EG_getRange);

                if (senses[jedge] == SFORWARD) {
                    status = EG_evaluate(efoo[jedge+nchild], &tt[0], &uv[2]);
                    CHECK_STATUS(EG_evaluate);
                } else {
                    status = EG_evaluate(efoo[jedge+nchild], &tt[1], &uv[2]);
                    CHECK_STATUS(EG_evaluate);
                }

                /* make the DEGENERATE Edge and associated PCurve */
                if (senses[iedge] != SFORWARD) {
                    swap  = uv[2];
                    uv[2] = uv[0];
                    uv[0] = swap;

                    swap  = uv[3];
                    uv[3] = uv[1];
                    uv[1] = swap;
                }

                uv[2] -= uv[0];
                uv[3] -= uv[1];

                range[0] = 0;
                range[1] = sqrt(uv[2]*uv[2] + uv[3]*uv[3]);
                sense    = SFORWARD;

                status = EG_makeTopology(context, NULL, EDGE, DEGENERATE, range, 1, etemp, &sense, &efoo[iedge]);
                CHECK_STATUS(EG_makeTopology);

                status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, uv, &efoo[iedge+nchild]);
                CHECK_STATUS(EG_makeGeometry);
            }

            /* make the Loop */
            status = EG_makeTopology(context, esurfaces[iface], oclass, mtype,
                                     NULL, nchild, efoo, senses, &eloops[iloop]);
            CHECK_STATUS(EG_makeTopology);

            FREE(efoo);
        }

        /* make the Face */
        status = EG_getTopology(efaces_orig[iface], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        status = EG_makeTopology(context, esurfaces[iface], oclass, mtype,
                                 xyz, nloop, eloops, senses, &efaces[iface]);
        CHECK_STATUS(EG_makeTopology);

        FREE(eloops);

        status = EG_attributeDup(efaces_orig[iface], efaces[iface]);
        CHECK_STATUS(EG_attributeDup);
    }

    /* make the SHELLs */
    for (ishell = 0; ishell < nshell; ishell++) {
        SPLINT_CHECK_FOR_NULL(eshells_orig);

        status = EG_getTopology(eshells_orig[ishell], &eref, &oclass, &mtype,
                                xyz, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        MALLOC(efoo, ego, nchild);

        for (ichild = 0; ichild < nchild; ichild++) {
            iface = EG_indexBodyTopo(inbody, echilds[ichild]) - 1;
            efoo[ichild] = efaces[iface];
        }

        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, efoo, senses, &eshells[ishell]);
        CHECK_STATUS(EG_makeTopology);
    }

    /* finally make the Body */
    status = EG_getTopology(inbody, &eref, &oclass, &mtype,
                            xyz, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (mtype == SHEETBODY || mtype == SOLIDBODY) {
        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, eshells, senses, ebody);
        CHECK_STATUS(EG_makeTopology);
    } else if (mtype == FACEBODY && nchild == 1) {
        status = EG_makeTopology(context, NULL, oclass, mtype,
                                 xyz, nchild, efaces, senses, ebody);
        CHECK_STATUS(EG_makeTopology);
    } else {
        printf(" convertToBSpline: oclass=%d, mtype=%d, nchild=%d\n", oclass, mtype, nchild);
        status = EGADS_TOPOERR;
        goto cleanup;
    }

    status = EG_attributeDup(inbody, *ebody);
    CHECK_STATUS(EG_attributeDup);

cleanup:
    if (enodes_orig  != NULL) EG_free(enodes_orig );
    if (eedges_orig  != NULL) EG_free(eedges_orig );
    if (eloops_orig  != NULL) EG_free(eloops_orig );
    if (efaces_orig  != NULL) EG_free(efaces_orig );
    if (eshells_orig != NULL) EG_free(eshells_orig);

    FREE(enodes   );
    FREE(ecurves  );
    FREE(eedges   );
    FREE(esurfaces);
    FREE(efaces   );
    FREE(eshells  );
    FREE(eloops   );
    FREE(efoo     );

    return status;
}
