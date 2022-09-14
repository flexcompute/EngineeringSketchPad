/*
 ************************************************************************
 *                                                                      *
 * udpFitcurve -- udp file to fit a Bspline curve to a set of points    *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
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

#define NUMUDPARGS 9
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)    ((char   *) (udps[IUDP].arg[0].val))
#define NCP(     IUDP)    ((int    *) (udps[IUDP].arg[1].val))[0]
#define ORDERED( IUDP)    ((int    *) (udps[IUDP].arg[2].val))[0]
#define PERIODIC(IUDP)    ((int    *) (udps[IUDP].arg[3].val))[0]
#define SPLIT(   IUDP,I)  ((int    *) (udps[IUDP].arg[4].val))[I]
#define XFORM(   IUDP,I)  ((double *) (udps[IUDP].arg[5].val))[I]
#define NPNT(    IUDP)    ((int    *) (udps[IUDP].arg[6].val))[0]
#define RMS(     IUDP)    ((double *) (udps[IUDP].arg[7].val))[0]
#define XYZ(     IUDP,I)  ((double *) (udps[IUDP].arg[8].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename", "ncp",   "ordered", "periodic", "split", "xform",  "npnt",   "rms",     "xyz", };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,   ATTRINT, ATTRINT,   ATTRINT,    ATTRINT, ATTRREAL, -ATTRINT, -ATTRREAL, 0,     };
static int    argIdefs[NUMUDPARGS] = {0,          0,       1,         0,          0,       0,        0,        0,         0,     };
static double argDdefs[NUMUDPARGS] = {0.,         0.,      1.,        0.,         0.,      0.,       0.,       0.,        0.,    };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
#include <assert.h>

#ifdef GRAFIC
   #include "grafic.h"
   #define DEBUG
#endif

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) > (B)) ? (A) : (B))
#define           EPS06           1.0e-6
#define           HUGEQ           1.0e+20

static int    EG_fitBspline(ego context,
                            int npnt, int bitflag, double xyz[],
                            int ncp, ego *ecurve, double *rms);
static int    fit1dCloud(int m, int ordered, double XYZcloud[],
                         int n, double cp[], double *normf);
static int    eval1dBspline(double T, int n, double cp[], double XYZ[],
                            /*@null@*/double dXYZdT[], /*@null@*/double dXYZdP[]);
static int    cubicBsplineBases(int ncp, double t, double N[], double dN[]);
static int    solveSparse(double SAv[], int SAi[], double b[], double x[],
                          int itol, double *errmax, int *iter);
static double L2norm(double f[], int n);

#ifdef GRAFIC
static int    plotCurve(int npnt, ego ecurve);
#endif

static void *realloc_temp=NULL;              /* used by RALLOC macro */


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

    int     *senses=NULL;
    int     ipnt, jpnt, mpnt, npnt, iedge, nedge, nument, idum, i, ieof;
    int     bitflag, wraparound, periodic;
    double  xin, yin, zin, rms, range[4], norm[3], range_save;
    double  data[18], uv_out[2], xyz_out[3], *cpdata=NULL;
    char    *message=NULL, *filename, *token;
    FILE    *fp=NULL;
    ego     ecurve, eloop, eface, enew;
    ego     *enodes=NULL, *eedges=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename(0) = %s\n", FILENAME(0));
    printf("ncp(     0) = %d\n", NCP(     0));
    printf("ordered( 0) = %d\n", ORDERED( 0));
    printf("periodic(0) = %d\n", PERIODIC(0));
    printf("split(   0) = %d",   SPLIT(   0,0));
    for (i = 1; i < udps[0].arg[4].size; i++) {
        printf(" %d", SPLIT(0,i));
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
    if (strlen(FILENAME(0)) == 0) {
        snprintf(message, 100, "filename must not be null\n");
        status  = EGADS_NODATA;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "ncp should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (NCP(0) < 3) {
        snprintf(message, 100, "ncp = %d < 3\n", NCP(0));
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[2].size > 1) {
        snprintf(message, 100, "ordered should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (ORDERED(0) != 0 && ORDERED(0) != 1) {
        snprintf(message, 100, "ordered should be 0 or 1\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[3].size > 1) {
        snprintf(message, 100, "periodic should be a scalar\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (PERIODIC(0) != 0 && PERIODIC(0) != 1) {
        snprintf(message, 100, "periodic should be 0 or 1\n");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[5].size != 1 && udps[0].arg[5].size != 12) {
        snprintf(message, 100, "xform should have 1 or 12 elements\n");
        status = EGADS_RANGERR;
        goto cleanup;

    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
    printf("ncp(     %d) = %d\n", numUdp, NCP(     numUdp));
    printf("ordered( %d) = %d\n", numUdp, ORDERED( numUdp));
    printf("periodic(%d) = %d\n", numUdp, PERIODIC(numUdp));
    printf("split(   %d) = %d",   numUdp, SPLIT(   numUdp,0));
    for (i = 1; i < udps[0].arg[4].size; i++) {
        printf(" %d", SPLIT(numUdp,i));
    }
    printf("\n");
#endif

    /* open the file or stream */
    filename = FILENAME(numUdp);

    if (strncmp(filename, "<<\n", 3) == 0) {
        token = strtok(filename, " \t\n");
        if (token == NULL) {
            ieof = 1;
        } else {
            ieof = 0;
        }
    } else {
        fp = fopen(filename, "r");
        if (fp == NULL) {
            snprintf(message, 100, "could not open file \"%s\"", filename);
            status = EGADS_NOTFOUND;
            goto cleanup;
        }

        ieof = feof(fp);
    }
    if (ieof == 1) {
        snprintf(message, 100, "premature enf-of-file found");
        status = EGADS_NODATA;
        goto cleanup;
    }

    mpnt = 100;
    RALLOC(udps[numUdp].arg[8].val, double, 3*mpnt);

    /* fill the table of points by re-reading the file */
    npnt  = 0;
    nedge = 1;
    while (1) {
        if (npnt >= mpnt-1) {
            mpnt += 100;
            RALLOC(udps[numUdp].arg[8].val, double, 3*mpnt);
        }

        if (fp != NULL) {
            nument = fscanf(fp, "%lf %lf %lf\n", &xin, &yin, &zin);
            if (nument != 3) break;
        } else {
            token = strtok(NULL, " \t\n");
            if (token == NULL) break;
            xin = strtod(token, NULL);

            token = strtok(NULL, " \t\n");
            if (token == NULL) break;
            yin = strtod(token, NULL);

            token = strtok(NULL, " \t\n");
            if (token == NULL) break;
            zin = strtod(token, NULL);
        }

        if (udps[numUdp].arg[5].size == 1) {
            XYZ(numUdp,3*npnt  ) = xin;
            XYZ(numUdp,3*npnt+1) = yin;
            XYZ(numUdp,3*npnt+2) = zin;
        } else {
            XYZ(numUdp,3*npnt  ) = XFORM(numUdp, 0) * xin + XFORM(numUdp, 1) * yin + XFORM(numUdp, 2) * zin + XFORM(numUdp, 3);
            XYZ(numUdp,3*npnt+1) = XFORM(numUdp, 4) * xin + XFORM(numUdp, 5) * yin + XFORM(numUdp, 6) * zin + XFORM(numUdp, 7);
            XYZ(numUdp,3*npnt+2) = XFORM(numUdp, 8) * xin + XFORM(numUdp, 9) * yin + XFORM(numUdp,10) * zin + XFORM(numUdp,11);
        }
        npnt++;

        if (npnt > 1 &&
            fabs(XYZ(numUdp,3*npnt-6)-XYZ(numUdp,3*npnt-3)) < EPS06 &&
            fabs(XYZ(numUdp,3*npnt-5)-XYZ(numUdp,3*npnt-2)) < EPS06 &&
            fabs(XYZ(numUdp,3*npnt-4)-XYZ(numUdp,3*npnt-1)) < EPS06   ) {
            nedge++;
        }

        for (i = 0; i < udps[0].arg[4].size; i++) {
            if (npnt == SPLIT(0,i)+i) {
                XYZ(numUdp,3*npnt  ) = XYZ(numUdp,3*npnt-3);
                XYZ(numUdp,3*npnt+1) = XYZ(numUdp,3*npnt-2);
                XYZ(numUdp,3*npnt+2) = XYZ(numUdp,3*npnt-1);

                npnt++;
                nedge++;
            }
        }
    }
#ifdef DEBUG
    printf("nedge=%d\n", nedge);
#endif

    if (fp != NULL) {
        fclose(fp);
    }

    /* cannot have a wraparound geometry that only has one Edge */
    if (fabs(XYZ(numUdp,3*npnt-3)-XYZ(numUdp,0)) < EPS06 &&
        fabs(XYZ(numUdp,3*npnt-2)-XYZ(numUdp,1)) < EPS06 &&
        fabs(XYZ(numUdp,3*npnt-1)-XYZ(numUdp,2)) < EPS06   ) {
        wraparound = 1;
    } else {
        wraparound = 0;
    }
#ifdef DEBUG
    printf("wraparound=%d\n", wraparound);
#endif

    if (wraparound == 1 && nedge == 1) {
        snprintf(message, 100, "wraparound geometry with only one Edge\n");
        status = EGADS_DEGEN;
        goto cleanup;
    }

    /* fit a Bspline to the data */
    bitflag = ORDERED(numUdp) + 2 * PERIODIC(numUdp);
    status = EG_fitBspline(context, npnt, bitflag, udps[numUdp].arg[8].val,
                           NCP(numUdp), &ecurve, &rms);
    CHECK_STATUS(EG_fitBspline);

#ifdef GRAFIC
    /* plot the fit */
    status = plotCurve(npnt, ecurve);
    CHECK_STATUS(plotCurve);
#endif

    /* get storage of the Nodes and Edges */
    MALLOC(enodes, ego, (nedge+1));
    MALLOC(eedges, ego, (nedge  ));
    MALLOC(senses, int, (nedge  ));

    /* make the Nodes and the Edges */
    status = EG_getRange(ecurve, range, &idum);
    CHECK_STATUS(EG_getRange);

    range_save = range[1];

    status = EG_evaluate(ecurve, &(range[0]), data);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL, &(enodes[0]));
    CHECK_STATUS(EG_makeTopology);

    status = EG_evaluate(ecurve, &(range[1]), data);
    CHECK_STATUS(EG_evaluate);

    status = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL, &(enodes[nedge]));
    CHECK_STATUS(EG_makeTopology);

    jpnt = 1;
    for (iedge = 0; iedge < nedge-1; iedge++) {
        for (ipnt = jpnt; ipnt < npnt; ipnt++) {
            if (fabs(XYZ(numUdp,3*ipnt  )-XYZ(numUdp,3*ipnt-3)) < EPS06 &&
                fabs(XYZ(numUdp,3*ipnt+1)-XYZ(numUdp,3*ipnt-2)) < EPS06 &&
                fabs(XYZ(numUdp,3*ipnt+2)-XYZ(numUdp,3*ipnt-1)) < EPS06   ) {

                status = EG_invEvaluate(ecurve, &(XYZ(numUdp,3*ipnt)), uv_out, xyz_out);
                CHECK_STATUS(EG_invEvaluate);

                range[1] = uv_out[0];

                status = EG_makeTopology(context, NULL, NODE, 0, xyz_out, 0, NULL, NULL, &(enodes[iedge+1]));
                CHECK_STATUS(EG_makeTopology);

                status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                         range, 2, &(enodes[iedge]), NULL, &(eedges[iedge]));
                CHECK_STATUS(EG_makeTopology);

                senses[iedge] = SFORWARD;

                range[0] = range[1];
                jpnt     = ipnt + 1;

                break;
            }
        }
    }

    range[1] = range_save;

    /* re-use first Node if wrap-around */
    if (wraparound == 1) {
        status = EG_deleteObject(enodes[nedge]);
        CHECK_STATUS(EG_deleteObject);

        enodes[nedge] = enodes[0];
    }

    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                             range, 2, &(enodes[nedge-1]), NULL, &(eedges[nedge-1]));
    CHECK_STATUS(EG_makeTopology);

    senses[nedge-1] = SFORWARD;

    /* make a Loop */
    if (wraparound == 0) {
        status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                 NULL, nedge, eedges, senses, &eloop);
        CHECK_STATUS(EG_makeTopology);
    } else {
        status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                 NULL, nedge, eedges, senses, &eloop);
        CHECK_STATUS(EG_makeTopology);
    }

    /* make a WireBody or a FaceBody */
    if (wraparound == 0) {
        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);
    } else {

        /* make Face from the loop */
        status = EG_makeFace(eloop, SFORWARD, NULL, &eface);

        /* if Face could not be made, make a WIREBODY instead */
        if (status != EGADS_SUCCESS) {
            status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                     NULL, 1, &eloop, NULL, ebody);
            CHECK_STATUS(EG_makeTopology);

        } else {

            /* find the direction of the Face normal */
            status = EG_getRange(eface, range, &periodic);
            CHECK_STATUS(EG_getRange);

            range[0] = (range[0] + range[1]) / 2;
            range[1] = (range[2] + range[3]) / 2;

            status = EG_evaluate(eface, range, data);
            CHECK_STATUS(EG_evaluate);

            norm[0] = data[4] * data[8] - data[5] * data[7];
            norm[1] = data[5] * data[6] - data[3] * data[8];
            norm[2] = data[3] * data[7] - data[4] * data[6];

            /* if the normal is not positive, flip the Face */
            if (norm[2] < 0) {
                status = EG_flipObject(eface, &enew);
                CHECK_STATUS(EG_flipObject);

                eface = enew;
            }

            /* create the FaceBody (which will be returned) */
            senses[0] = SFORWARD;
            status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                                     NULL, 1, &eface, senses, ebody);
            CHECK_STATUS(EG_makeTopology);
        }
    }

    /* set the output value(s) */
    NPNT(0) = npnt;
    RMS( 0) = rms;

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

cleanup:
    if (cpdata != NULL) free(cpdata);
    if (enodes != NULL) free(enodes);
    if (eedges != NULL) free(eedges);
    if (senses != NULL) free(senses);

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
 *   EG_fitBspline - fit a degree=3 Bspline curve to a set of points    *
 *                                                                      *
 ************************************************************************
 */

static int
EG_fitBspline(ego    context,           /* (in)  EGADS context */
              int    npnt,              /* (in)  number of points */
              int    bitflag,           /* (in)  1=ordered, 2=periodic */
              double xyz[],             /* (in)  array of points (xyzxyz...) */
              int    ncp,               /* (in)  number of control points */
              ego    *ecurve,           /* (out) Bspline curve (degree=3) */
              double *rms)              /* (out) RMS distance of points from curve */
{
    int    status = EGADS_SUCCESS;

    int    nknot, ndata, idata, j, header[4];
    double *cpdata=NULL;

    ROUTINE(EG_fitBspline);

    /* --------------------------------------------------------------- */

    /* default returns */
    *ecurve = NULL;
    *rms    = 0;

    /* check the inputs */
    if (context == NULL) {
        status = EGADS_NULLOBJ;
        goto cleanup;
    } else if (npnt < 2) {
        status = EGADS_NODATA;
        goto cleanup;
    } else if (xyz == NULL) {
        status = EGADS_NODATA;
        goto cleanup;
    } else if (ncp < 3) {
        status = EGADS_NODATA;
        goto cleanup;
    }

    /* set up arrays needed to define Bspline */
    nknot = ncp + 4;
    ndata = nknot + 3 * ncp;

    header[0] = 0;            // bitflag
    header[1] = 3;            // degree
    header[2] = ncp;          // number of control points
    header[3] = nknot;        // number of knots

    MALLOC(cpdata, double, ndata);
    ndata = 0;

    /* knot vector */
    cpdata[ndata++] = 0;
    cpdata[ndata++] = 0;
    cpdata[ndata++] = 0;
    cpdata[ndata++] = 0;

    for (j = 1; j < ncp-3; j++) {
        cpdata[ndata++] = j;
    }

    cpdata[ndata++] = ncp - 3;
    cpdata[ndata++] = ncp - 3;
    cpdata[ndata++] = ncp - 3;
    cpdata[ndata++] = ncp - 3;

    /* control points at the two ends */
    idata = ndata;
    cpdata[idata  ] = xyz[0];
    cpdata[idata+1] = xyz[1];
    cpdata[idata+2] = xyz[2];

    idata = ndata + 3 * (ncp-1);
    cpdata[idata  ] = xyz[3*npnt-3];
    cpdata[idata+1] = xyz[3*npnt-2];
    cpdata[idata+2] = xyz[3*npnt-1];

    /* perform the fitting (which updates the interior control points) */
    status = fit1dCloud(npnt, bitflag, xyz,
                        ncp,  &(cpdata[ndata]), rms);
    CHECK_STATUS(fit1dCloud);

    /* make the geometry */
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                             header, cpdata, ecurve);

cleanup:
    if (cpdata != NULL) free(cpdata);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit1dCloud - find spline that best-fits the cloud of points        *
 *                                                                      *
 ************************************************************************
 */

static int
fit1dCloud(int    m,                    /* (in)  number of points in cloud */
           int    bitflag,              /* (in)  1=ordered, 2=periodic */
           double XYZcloud[],           /* (in)  array  of points in cloud */
           int    n,                    /* (in)  number of control points */
           double cp[],                 /* (in)  array  of control points (first and last are set) */
                                        /* (out) array  of control points (all set) */
           double *normf)               /* (out) RMS of distances between cloud and fit */
{
    int    status = EGADS_SUCCESS;      /* (out)  return status */

    int    ordered=0, periodic=0, np, nvar, ivar, jvar, nobj, iobj, i, j, k, next;
    int    niter, iter, count, itol, maxiter;
    double xmin, xmax, ymin, ymax, zmin, zmax, scale, xcent, ycent, zcent;
    double frac, toler, lambda, normdelta, normfnew, delta0, delta1, delta2;
    double XYZ[3], dXYZdT[3], omega, dmax, dist1, dist2, dist;
    double xx, xb, xe, yy, yb, ye, zz, zb, ze, tt, dd, dmin;

    int    *MMi=NULL;
    double *MMd=NULL, errmax;
    double *XYZcopy=NULL, *dXYZdP=NULL, *cpnew=NULL;
    double *beta=NULL, *delta=NULL, *betanew=NULL;
    double *f=NULL,  *fnew=NULL;
    double *aa=NULL, *bb=NULL, *cc=NULL, *rhs=NULL;

    ROUTINE(fit1dCloud);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("enter fit1dCloud(m=%d, ordered=%d, n=%d)\n", m, ordered, n);
#endif

    assert(m > 1);                      // needed to avoid clang warning
    assert(n > 2);                      // needed to avoid clang warning

    /* default return */
    *normf  = 1e-12;

    /* extract ordered and periodic flags */
    if (bitflag == 1 || bitflag == 3) ordered  = 1;
    if (bitflag == 2 || bitflag == 3) periodic = 1;

    /* number of design variables and objectives */
    np   = 3 * n - 6;
    nvar = m + np;
    nobj = 3 * m;

    /* if m < n, then assume that the linear spline is the best fit */
    if (m < n) {
        for (j = 1; j < n-1; j++) {
            frac = (double)(j) / (double)(n-1);

            cp[3*j  ] = (1-frac) * cp[0] + frac * cp[3*n-3];
            cp[3*j+1] = (1-frac) * cp[1] + frac * cp[3*n-2];
            cp[3*j+2] = (1-frac) * cp[2] + frac * cp[3*n-1];
        }

#ifdef DEBUG
        printf("making linear fit because not enough points in cloud\n");
#endif
        goto cleanup;
    }

    /* allocate all temporary arrays */
    MALLOC(XYZcopy, double, 3*m);
    MALLOC(dXYZdP,  double,   n);
    MALLOC(cpnew,   double, 3*n);

    MALLOC(beta,    double, nvar);
    MALLOC(delta,   double, nvar);
    MALLOC(betanew, double, nvar);

    MALLOC(f,       double, nobj);
    MALLOC(fnew,    double, nobj);

    MALLOC(aa,      double, m    );
    MALLOC(bb,      double, m *np);
    MALLOC(cc,      double, np*np);
    MALLOC(rhs,     double, nvar );

#define AA(I)       aa[(I)]
#define BB(I,J)     bb[(J)+np*(I)]
#define CC(I,J)     cc[(J)+np*(I)]

    /* transform inputs so that they are centered at origin and
       unit length */
    xmin = XYZcloud[0];
    xmax = XYZcloud[0];
    ymin = XYZcloud[1];
    ymax = XYZcloud[1];
    zmin = XYZcloud[2];
    zmax = XYZcloud[2];

    for (k = 1; k < m; k++) {
        if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
        if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
        if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
        if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
        if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
        if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
    }

    scale = 1.0 / MAX(MAX(xmax-xmin, ymax-ymin), zmax-zmin);
    xcent = scale * (xmin + xmax) / 2;
    ycent = scale * (ymin + ymax) / 2;
    zcent = scale * (zmin + zmax) / 2;

    for (k = 0; k < m; k++) {
        XYZcopy[3*k  ] = scale * (XYZcloud[3*k  ] - xcent);
        XYZcopy[3*k+1] = scale * (XYZcloud[3*k+1] - ycent);
        XYZcopy[3*k+2] = scale * (XYZcloud[3*k+2] - zcent);
    }
    for (j = 0; j < n; j++) {
        cp[3*j  ] = scale * (cp[3*j  ] - xcent);
        cp[3*j+1] = scale * (cp[3*j+1] - ycent);
        cp[3*j+2] = scale * (cp[3*j+2] - zcent);
    }

    /* set up the initial values for the interior control
       points and the initial values of "t" */

    /* XYZcopy is ordered */
    if (ordered == 1) {

        /* set the initial control point locations by picking up evenly
           spaced points (based upon point number) from the cloud */
        for (j = 1; j < n-1; j++) {
            i = (j * (m-1)) / (n-1);

            cp[3*j  ] = XYZcopy[3*i  ];
            cp[3*j+1] = XYZcopy[3*i+1];
            cp[3*j+2] = XYZcopy[3*i+2];
        }

        /* for each point in the cloud, assign the value of "t"
           (which is stored in the first m betas) based upon it
           local pseudo-arc-length */
        beta[0] = 0;
        for (k = 1; k < m; k++) {
            beta[k] = beta[k-1] + sqrt(SQR(XYZcopy[3*k  ]-XYZcopy[3*k-3])
                                      +SQR(XYZcopy[3*k+1]-XYZcopy[3*k-2])
                                      +SQR(XYZcopy[3*k+2]-XYZcopy[3*k-1]));
        }

        for (k = 0; k < m; k++) {
            beta[k] = (n-3) * beta[k] / beta[m-1];
        }

    /* XYZcopy is unordered */
    } else {

        /* set the "center" control point to coincide with the point
           in the cloud that is furthest away from the first and
           last control points */
        dmax = 0;
        for (k = 1; k < m-1; k++) {
            dist1 = SQR(XYZcopy[3*k  ]-cp[    0])
                  + SQR(XYZcopy[3*k+1]-cp[    1])
                  + SQR(XYZcopy[3*k+2]-cp[    2]);
            dist2 = SQR(XYZcopy[3*k  ]-cp[3*n-3])
                  + SQR(XYZcopy[3*k+1]-cp[3*n-2])
                  + SQR(XYZcopy[3*k+2]-cp[3*n-1]);
            dist  = MIN(dist1, dist2);

            if (dist > dmax) {
                dmax = dist;
                cp[3*(n/2)  ] = XYZcopy[3*k  ];
                cp[3*(n/2)+1] = XYZcopy[3*k+1];
                cp[3*(n/2)+2] = XYZcopy[3*k+2];
            }
        }

        /* fill in the other control points */
        for (j = 1; j < (n/2); j++) {
            frac = (double)(j) / (double)(n/2);

            cp[3*j  ] = (1-frac) * cp[0] + frac * cp[3*(n/2)  ];
            cp[3*j+1] = (1-frac) * cp[1] + frac * cp[3*(n/2)+1];
            cp[3*j+2] = (1-frac) * cp[2] + frac * cp[3*(n/2)+2];
        }

        for (j = (n/2)+1; j < n; j++) {
            frac = (double)(j-(n/2)) / (double)(n-1-(n/2));

            cp[3*j  ] = (1-frac) * cp[3*(n/2)  ] + frac * cp[3*n-3];
            cp[3*j+1] = (1-frac) * cp[3*(n/2)+1] + frac * cp[3*n-2];
            cp[3*j+2] = (1-frac) * cp[3*(n/2)+2] + frac * cp[3*n-1];
        }

        /* for each point in the cloud, assign the value of "t"
           (which is stored in the first m betas) as the closest
           point to the control polygon */
        for (k = 0; k < m; k++) {
            xx = XYZcopy[3*k  ];
            yy = XYZcopy[3*k+1];
            zz = XYZcopy[3*k+2];

            dmin = HUGEQ;
            for (j = 1; j < n; j++) {
                xb = cp[3*j-3];   yb = cp[3*j-2];   zb = cp[3*j-1];
                xe = cp[3*j  ];   ye = cp[3*j+1];   ze = cp[3*j+2];

                tt = ((xe-xb) * (xx-xb) + (ye-yb) * (yy-yb) + (ze-zb) * (zz-zb))
                   / ((xe-xb) * (xe-xb) + (ye-yb) * (ye-yb) + (ze-zb) * (ze-zb));
                tt = MIN(MAX(0, tt), 1);

                dd = SQR((1-tt) * xb + tt * xe - xx)
                   + SQR((1-tt) * yb + tt * ye - yy)
                   + SQR((1-tt) * zb + tt * ze - zz);

                if (dd < dmin) {
                    dmin    = dd;
                    beta[k] = ((j-1) + tt) * (double)(n - 3) / (double)(n - 1);
                }
            }
        }
    }

#ifdef DEBUG
    printf("Initialization\n");
    for (j = 0; j < n; j++) {
        printf("%3d: %12.6f %12.6f %12.6f\n", j, cp[3*j], cp[3*j+1], cp[3*j+2]);
    }
    for (k = 0; k < m; k++) {
        printf("%3d: %12.6f\n", k, beta[k]);
    }
#endif

    /* set the relaxation parameter for control points */
    omega = 0.25;

    /* insert the interior control points into the design variables */
    next = m;
    for (j = 1; j < n-1; j++) {
        beta[next++] = cp[3*j  ];
        beta[next++] = cp[3*j+1];
        beta[next++] = cp[3*j+2];
    }

    /* compute the initial objective function */
    for (k = 0; k < m; k++) {
        status = eval1dBspline(beta[k], n, cp, XYZ, NULL, NULL);
        CHECK_STATUS(eval1dBspline);

        f[3*k  ] = XYZcopy[3*k  ] - XYZ[0];
        f[3*k+1] = XYZcopy[3*k+1] - XYZ[1];
        f[3*k+2] = XYZcopy[3*k+2] - XYZ[2];
    }
    *normf = L2norm(f, nobj) / m;
#ifdef DEBUG
    printf("initial   norm(f)=%11.4e\n", *normf);
#endif

    /* initialize the Levenberg-Marquardt algorithm */
    niter  = 501;
    toler  = 1.0e-6;
    lambda = 1;

    /* LM iterations */
    for (iter = 0; iter < niter; iter++) {

        /* initialize [AA  BB]
                      [      ] =  transpose(J) * J + lambda * diag(transpose(J) * J)
                      [BB' CC]

           and        rhs  = -transpose(J) * f
        */
        for (jvar = 0; jvar < np; jvar++) {
            for (ivar = 0; ivar < np; ivar++) {
                CC(ivar,jvar) = 0;
            }
            CC(jvar,jvar) = 1e-6;
        }

        for (jvar = 0; jvar < nvar; jvar++) {
            rhs[jvar] = 0;
        }

        /* accumulate AA, BB, CC, and rhs by looping over points in cloud */
        for (k = 0; k < m; k++) {
            status = eval1dBspline(beta[k], n, cp, XYZ, dXYZdT, dXYZdP);
            CHECK_STATUS(eval1dBspline);

            AA(k) = dXYZdT[0] * dXYZdT[0] + dXYZdT[1] * dXYZdT[1] + dXYZdT[2] * dXYZdT[2];

            for (ivar = 1; ivar < n-1; ivar++) {
                BB(k, 3*ivar-3) = dXYZdT[0] * dXYZdP[ivar];
                BB(k, 3*ivar-2) = dXYZdT[1] * dXYZdP[ivar];
                BB(k, 3*ivar-1) = dXYZdT[2] * dXYZdP[ivar];

                for (jvar = 1; jvar < n-1; jvar++) {
#ifndef __clang_analyzer__
                    CC(3*ivar-3, 3*jvar-3) += dXYZdP[ivar] * dXYZdP[jvar];
                    CC(3*ivar-2, 3*jvar-2) += dXYZdP[ivar] * dXYZdP[jvar];
                    CC(3*ivar-1, 3*jvar-1) += dXYZdP[ivar] * dXYZdP[jvar];
#endif
                }
            }

            rhs[k] = dXYZdT[0] * f[3*k] + dXYZdT[1] * f[3*k+1] + dXYZdT[2] * f[3*k+2];

            for (ivar = 1; ivar < n-1; ivar++) {
#ifndef __clang_analyzer__
                rhs[m+3*ivar-3] += dXYZdP[ivar] * f[3*k  ];
                rhs[m+3*ivar-2] += dXYZdP[ivar] * f[3*k+1];
                rhs[m+3*ivar-1] += dXYZdP[ivar] * f[3*k+2];
#endif
            }
        }

        /* set up sparse-matrix arrays */
        count = m + 2 * m * np + np * np + 1;

        MALLOC(MMd, double, count);
        MALLOC(MMi, int,    count);

        /* store diagonal values (multiplied by (1+lambda)) */
        for (k = 0; k < m; k++) {
            MMd[k] = AA(k) * (1 + lambda);
        }
        for (ivar = 0; ivar < np; ivar++) {
            MMd[m+ivar] = CC(ivar,ivar) * (1 + lambda);
        }

        /* set up off-diagonal elements, including indices */
        MMi[0] = nvar + 1;
        count  = nvar;

        /* BB to the right of AA */
        for (k = 0; k < m; k++) {
            for (jvar = 0; jvar < np; jvar++) {
                count++;
                MMd[count] = BB(k,jvar);
                MMi[count] = m + jvar;
            }
            MMi[k+1] = count + 1;
        }

        for (ivar = 0; ivar < np; ivar++) {
            /* transpose(BB) below A */
            for (k = 0; k < m; k++) {
                count++;
                MMd[count] = BB(k,ivar);
                MMi[count] = k;
            }

            /* CC in bottom-right corner */
            for (jvar = 0; jvar < np; jvar++) {
                if (ivar != jvar) {
                    count++;
                    MMd[count] = CC(ivar,jvar);
                    MMi[count] = m + jvar;
                }
            }
            MMi[m+ivar+1] = count + 1;
        }

        /* arbitrary value (not used) */
        MMd[nvar] = 0;

        /* set up for sparse matrix solve (via biconjugate gradient technique) */
        itol    = 1;
        errmax  = 1.0e-12;
        maxiter = 2 * nvar;
        for (ivar = 0; ivar < nvar; ivar++) {
            delta[ivar] = 0;
        }

        status = solveSparse(MMd, MMi, rhs, delta, itol, &errmax, &maxiter);
        CHECK_STATUS(solveSparse);

        FREE(MMd);
        FREE(MMi);

        /* check for convergence on delta (which corresponds to a small
           change in beta) */
        normdelta = L2norm(delta, nvar);

        if (normdelta < toler) {
#ifdef DEBUG
            printf("converged with norm(delta)=%11.4e\n", normdelta);
#endif
            break;
        }

        /* find the temporary new beta */
        for (ivar = 0; ivar < nvar; ivar++) {

            /* beta associated with Tcloud */
            if (ivar < m) {
                betanew[ivar] = beta[ivar] + delta[ivar];

                if (betanew[ivar] < 0  ) betanew[ivar] = 0;
                if (betanew[ivar] > n-3) betanew[ivar] = n-3;

                /* beta associated with control points */
            } else {
                betanew[ivar] = beta[ivar] + omega * delta[ivar];
            }
        }

        /* gradually increase omega */
        omega = MIN(1.01*omega, 1.0);

        /* extract the temporary control points from betanew */
        next = m;
        for (j = 0; j < n; j++) {
            if (j == 0 || j == n-1) {
                cpnew[3*j  ] = cp[3*j  ];
                cpnew[3*j+1] = cp[3*j+1];
                cpnew[3*j+2] = cp[3*j+2];
            } else {
                cpnew[3*j  ] = betanew[next++];
                cpnew[3*j+1] = betanew[next++];
                cpnew[3*j+2] = betanew[next++];
            }
        }

        /* apply periodicity condition by making sure first and last
           intervals are the same */
        if (periodic == 1) {
            delta0 = (2*cpnew[0] - cpnew[3] - cpnew[3*n-6]) / 2;
            delta1 = (2*cpnew[1] - cpnew[4] - cpnew[3*n-5]) / 2;
            delta2 = (2*cpnew[2] - cpnew[5] - cpnew[3*n-4]) / 2;

            cpnew[    3] += delta0;
            cpnew[    4] += delta1;
            cpnew[    5] += delta2;

            cpnew[3*n-6] += delta0;
            cpnew[3*n-5] += delta1;
            cpnew[3*n-4] += delta2;
        }

        /* compute the objective function based upon the new beta */
        for (k = 0; k < m; k++) {
            status = eval1dBspline(betanew[k], n, cp, XYZ, NULL, NULL);
            CHECK_STATUS(eval1dBspline);

            fnew[3*k  ] = XYZcopy[3*k  ] - XYZ[0];
            fnew[3*k+1] = XYZcopy[3*k+1] - XYZ[1];
            fnew[3*k+2] = XYZcopy[3*k+2] - XYZ[2];
        }
        normfnew = L2norm(fnew, nobj) / m;
#ifdef DEBUG
        if (iter%10 == 0) {
            printf("iter=%4d: norm(delta)=%11.4e, norm(f)=%11.4e  ",
                   iter, normdelta, normfnew);
        }
#endif

        /* if this was a better step, accept it and decrease
           lambda (making it more Newton-like) */
        if (normfnew < *normf) {
            lambda /= 2;
#ifdef DEBUG
            if (iter%10 == 0) {
                printf("ACCEPTED,  lambda=%11.4e, omega=%10.5f\n", lambda, omega);
            }
#endif

            /* save new design variables, control points, and
               objective function */
            for (ivar = 0; ivar < nvar; ivar++) {
                beta[ivar] = betanew[ivar];
            }
            for (j = 0; j < n; j++) {
                cp[3*j  ] = cpnew[3*j  ];
                cp[3*j+1] = cpnew[3*j+1];
                cp[3*j+2] = cpnew[3*j+2];
            }
            for (iobj = 0; iobj < nobj; iobj++) {
                f[iobj] = fnew[iobj];
            }
            *normf = normfnew;

        /* otherwise do not take the step and increase lambda (making it
           more steepest-descent-like) */
        } else {
            lambda *= 2;
#ifdef DEBUG
            if (iter %10 == 0) {
                printf("rejected,  lambda=%11.4e, omega=%10.5f\n", lambda, omega);
            }
#endif
        }

        /* check for convergence (based upon a small value of
           objective function) */
        if (*normf < toler) {
#ifdef DEBUG
            printf("converged with norm(f)=%11.4e\n", *normf);
#endif
            break;
        }
    }

    /* transform control points back to their original scale */
    for (j = 0; j < n; j++) {
        cp[3*j  ] = xcent + cp[3*j  ] / scale;
        cp[3*j+1] = ycent + cp[3*j+1] / scale;
        cp[3*j+2] = zcent + cp[3*j+2] / scale;
    }

    *normf /= scale;

#ifdef DEBUG
    printf("final control points\n");
    for (j = 0; j < n; j++) {
        printf("%3d: %12.6f %12.6f %12.6f\n", j, cp[3*j], cp[3*j+1], cp[3*j+2]);
    }
    printf("*normf: %12.4e\n", *normf);
#endif

cleanup:
    if (XYZcopy != NULL) free(XYZcopy);
    if (dXYZdP  != NULL) free(dXYZdP );
    if (cpnew   != NULL) free(cpnew  );

    if (beta    != NULL) free(beta   );
    if (delta   != NULL) free(delta  );
    if (betanew != NULL) free(betanew);

    if (f       != NULL) free(f      );
    if (fnew    != NULL) free(fnew   );

    if (aa      != NULL) free(aa     );
    if (bb      != NULL) free(bb     );
    if (cc      != NULL) free(cc     );
    if (rhs     != NULL) free(rhs    );

    if (MMi     != NULL) free(MMi    );
    if (MMd     != NULL) free(MMd    );

#undef AA
#undef BB
#undef CC

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   eval1dBspline - evaluate cubic Bspline and its derivatives         *
 *                                                                      *
 ************************************************************************
 */

static int
eval1dBspline(double T,                 /* (in)  independent variable */
              int    n,                 /* (in)  number of control points */
              double cp[],              /* (in)  array  of control points */
              double XYZ[],             /* (out) dependent variables */
    /*@null@*/double dXYZdT[],          /* (out) derivative wrt T (or NULL) */
    /*@null@*/double dXYZdP[])          /* (out) derivative wrt P (or NULL) */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    i, span;
    double N[4], dN[4];

    ROUTINE(eval1dBspline);

    /* --------------------------------------------------------------- */

    assert (n > 3);

    XYZ[0] = 0;
    XYZ[1] = 0;
    XYZ[2] = 0;

    /* set up the Bspline bases */
    status = cubicBsplineBases(n, T, N, dN);
    CHECK_STATUS(cubicBsplineBases);

    span = MIN(floor(T), n-4);

    /* find the dependent variable */
    for (i = 0; i < 4; i++) {
        XYZ[0] += N[i] * cp[3*(i+span)  ];
        XYZ[1] += N[i] * cp[3*(i+span)+1];
        XYZ[2] += N[i] * cp[3*(i+span)+2];
    }

    /* find the deriviative wrt T */
    if (dXYZdT != NULL) {
        dXYZdT[0] = 0;
        dXYZdT[1] = 0;
        dXYZdT[2] = 0;

        for (i = 0; i < 4; i++) {
            dXYZdT[0] += dN[i] * cp[3*(i+span)  ];
            dXYZdT[1] += dN[i] * cp[3*(i+span)+1];
            dXYZdT[2] += dN[i] * cp[3*(i+span)+2];
        }
    }

    /* find the derivative wrt P */
    if (dXYZdP != NULL) {
        for (i = 0; i < n; i++) {
            dXYZdP[i] = 0;
        }
        for (i = 0; i < 4; i++) {
            dXYZdP[i+span] += N[i];
        }
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   cubicBsplineBases - basis function values for cubic Bspline        *
 *                                                                      *
 ************************************************************************
 */

static int
cubicBsplineBases(int    ncp,           /* (in)  number of control points */
                  double T,             /* (in)  independent variable (0<=T<=(ncp-3) */
                  double N[],           /* (out) bases */
                  double dN[])          /* (out) d(bases)/d(T) */
{
    int       status = EGADS_SUCCESS;   /* (out) return status */

    int      i, r, span;
    double   saved, dsaved, num, dnum, den, dden, temp, dtemp;
    double   left[4], dleft[4], rite[4], drite[4];

    ROUTINE(cubicBsplineBases);

    /* --------------------------------------------------------------- */

    span = MIN(floor(T)+3, ncp-1);

    N[ 0] = 1.0;
    dN[0] = 0;

    for (i = 1; i <= 3; i++) {
        left[ i] = T - MAX(0, span-2-i);
        dleft[i] = 1;

        rite[ i] = MIN(ncp-3,span-3+i) - T;
        drite[i] =                     - 1;

        saved  = 0;
        dsaved = 0;

        for (r = 0; r < i; r++) {
            num   = N[ r];
            dnum  = dN[r];

            den   = rite[ r+1] + left[ i-r];
            dden  = drite[r+1] + dleft[i-r];

            temp  = num / den;
            dtemp = (dnum * den - dden * num) / den / den;

            N[ r] = saved  + rite[ r+1] * temp;
            dN[r] = dsaved + drite[r+1] * temp + rite[r+1] * dtemp;

            saved  = left[ i-r] * temp;
            dsaved = dleft[i-r] * temp + left[i-r] * dtemp;
        }

        N[ i] = saved;
        dN[i] = dsaved;
    }

//cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   solveSparse - solve: A * x = b  using biconjugate gradient method  *
 *                                                                      *
 ************************************************************************
 */

static int
solveSparse(double SAv[],               /* (in)  sparse array values */
            int    SAi[],               /* (in)  sparse array indices */
            double b[],                 /* (in)  rhs vector */
            double x[],                 /* (in)  guessed result vector */
                                        /* (out) result vector */
            int    itol,                /* (in)  stopping criterion */
            double *errmax,             /* (in)  convergence tolerance */
                                        /* (out) estimated error at convergence */
            int    *iter)               /* (in)  maximum number of iterations */
                                        /* (out) number of iterations taken */
{
    int    status = EGADS_SUCCESS;      /* (out) return status */

    int    n, i, j, k, itmax;
    double tol, ak, akden, bk, bknum, bkden=1, bnorm, dxnorm, xnorm, znorm_old, znorm=0;
    double *p=NULL, *pp=NULL, *r=NULL, *rr=NULL, *z=NULL, *zz=NULL;

    ROUTINE(solveSparse);

    /* --------------------------------------------------------------- */

    tol   = *errmax;
    itmax = *iter;

    n = SAi[0] - 1;

    MALLOC(p,  double, n);
    MALLOC(pp, double, n);
    MALLOC(r,  double, n);
    MALLOC(rr, double, n);
    MALLOC(z,  double, n);
    MALLOC(zz, double, n);

    /* make sure none of the diagonals are very small */
    for (i = 0; i < n; i++) {
        if (fabs(SAv[i]) < 1.0e-14) {
            printf(" solveSparse: cannot solve since SAv[%d]=%11.4e\n", i, SAv[i]);
            status = EGADS_DEGEN;
            goto cleanup;
        }
    }

    /* calculate initial residual */
    *iter = 0;

    /* r = A * x  */
    for (i = 0; i < n; i++) {
        r[i] = SAv[i] * x[i];

        for (k = SAi[i]; k < SAi[i+1]; k++) {
            r[i] += SAv[k] * x[SAi[k]];
        }
    }

    for (j = 0; j < n; j++) {
        r[ j] = b[j] - r[j];
        rr[j] = r[j];
    }

    if (itol == 1) {
        bnorm = L2norm(b, n);

        for (j = 0; j < n; j++) {
            z[j] = r[j] / SAv[j];
        }
    } else if (itol == 2) {
        for (j = 0; j < n; j++) {
            z[j] = b[j] / SAv[j];
        }

        bnorm = L2norm(z, n);

        for (j = 0; j < n; j++) {
            z[j] = r[j] / SAv[j];
        }
    } else {
        for (j = 0; j < n; j++) {
            z[j] = b[j] / SAv[j];
        }

        bnorm = L2norm(z, n);

        for (j = 0; j < n; j++) {
            z[j] = r[j] / SAv[j];
        }

        znorm = L2norm(z, n);
    }

    /* main iteration loop */
    for (*iter = 0; *iter < itmax; (*iter)++) {

        for (j = 0; j < n; j++) {
            zz[j] = rr[j] / SAv[j];
        }

        /* calculate coefficient bk and direction vectors p and pp */
        bknum = 0;
        for (j = 0; j < n; j++) {
            bknum += z[j] * rr[j];
        }

        if (*iter == 0) {
            for (j = 0; j < n; j++) {
                p[ j] = z[ j];
                pp[j] = zz[j];
            }
        } else {
            bk = bknum / bkden;

            for (j = 0; j < n; j++) {
                p[ j] = bk * p[ j] + z[ j];
                pp[j] = bk * pp[j] + zz[j];
            }
        }

        /* calculate coefficient ak, new iterate x, and new residuals r and rr */
        bkden = bknum;

        /* z = A * p  */
        for (i = 0; i < n; i++) {
            z[i] = SAv[i] * p[i];

            for (k = SAi[i]; k < SAi[i+1]; k++) {
                z[i] += SAv[k] * p[SAi[k]];
            }
        }

        akden = 0;
        for (j = 0; j < n; j++) {
            akden += z[j] * pp[j];
        }

        ak = bknum / akden;

        /* zz = transpose(A) * pp  */
        for (i = 0; i < n; i++) {
            zz[i] = SAv[i] * pp[i];
        }

        for (i = 0; i < n; i++) {
            for (k = SAi[i]; k < SAi[i+1]; k++) {
                j = SAi[k];
                zz[j] += SAv[k] * pp[i];
            }
        }

        for (j = 0; j < n; j++) {
            x[ j] += ak * p[ j];
            r[ j] -= ak * z[ j];
            rr[j] -= ak * zz[j];
        }

        /* solve Abar * z = r */
        for (j = 0; j < n; j++) {
            z[j] = r[j] / SAv[j];
        }

        /* compute and check stopping criterion */
        if (itol == 1) {
            *errmax = L2norm(r, n) / bnorm;
        } else if (itol == 2) {
            *errmax = L2norm(z, n) / bnorm;
        } else {
            znorm_old = znorm;
            znorm = L2norm(z, n);
            if (fabs(znorm_old-znorm) > (1.0e-14)*znorm) {
                dxnorm = fabs(ak) * L2norm(p, n);
                *errmax  = znorm / fabs(znorm_old-znorm) * dxnorm;
            } else {
                *errmax = znorm / bnorm;
                continue;
            }

            xnorm = L2norm(x, n);
            if (*errmax <= xnorm/2) {
                *errmax /= xnorm;
            } else {
                *errmax = znorm / bnorm;
                continue;
            }
        }

        /* exit if converged */
        if (*errmax <= tol) break;
    }

cleanup:
    if (p  != NULL) free(p );
    if (pp != NULL) free(pp);
    if (r  != NULL) free(r );
    if (rr != NULL) free(rr);
    if (z  != NULL) free(z );
    if (zz != NULL) free(zz);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   L2norm - L2-norm of vector                                         *
 *                                                                      *
 ************************************************************************
 */

static double
L2norm(double f[],                      /* (in)  vector */
       int    n)                        /* (in)  length of vector */
{
    double L2norm;                      /* (out) L2-norm */

    int    i;

    ROUTINE(L2norm);

    /* --------------------------------------------------------------- */

    /* L2-norm */
    L2norm = 0;

    for (i = 0; i < n; i++) {
        L2norm += f[i] * f[i];
    }

    L2norm = sqrt(L2norm);

//cleanup:
    return L2norm;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotCurve - plot geometry, control polygon, and fit                *
 *                                                                      *
 ************************************************************************
 */

#ifdef GRAFIC
static int
plotCurve(int    npnt,                  /* (in)  number of points in cloud */
          ego    ecurve)                /* (in)  Curve */
{
    int    status = EGADS_SUCCESS;

    int    io_kbd=5, io_scr=6, nline=0, nplot=0;
    int    indgr=1+2+4+16+64, periodic, oclass, mtype, *tempIlist, icp, ipnt;
    int    ilin[3], isym[3], nper[3];
    int    ncp=NCP(numUdp);
    float  *xplot=NULL, *yplot=NULL, *zplot=NULL;
    double xmin, xmax, ymin, ymax, zmin, zmax, trange[4], frac, tt, data[18];
    double *tempRlist;
    ego    eref;

    ROUTINE(plotCurve);

    /* --------------------------------------------------------------- */

    MALLOC(xplot, float, (npnt+1000+ncp));
    MALLOC(yplot, float, (npnt+1000+ncp));
    MALLOC(zplot, float, (npnt+1000+ncp));

    xmin = XYZ(numUdp,0);
    xmax = XYZ(numUdp,0);
    ymin = XYZ(numUdp,1);
    ymax = XYZ(numUdp,1);
    zmin = XYZ(numUdp,2);
    zmax = XYZ(numUdp,2);

    /* build plot arrays for data points */
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        xplot[nplot] = XYZ(numUdp,3*ipnt  );
        yplot[nplot] = XYZ(numUdp,3*ipnt+1);
        zplot[nplot] = XYZ(numUdp,3*ipnt+2);
        nplot++;

        if (XYZ(numUdp,3*ipnt  ) < xmin) xmin = XYZ(numUdp,3*ipnt  );
        if (XYZ(numUdp,3*ipnt  ) > xmax) xmax = XYZ(numUdp,3*ipnt  );
        if (XYZ(numUdp,3*ipnt+1) < ymin) ymin = XYZ(numUdp,3*ipnt+1);
        if (XYZ(numUdp,3*ipnt+1) > ymax) ymax = XYZ(numUdp,3*ipnt+1);
        if (XYZ(numUdp,3*ipnt+2) < zmin) zmin = XYZ(numUdp,3*ipnt+2);
        if (XYZ(numUdp,3*ipnt+2) > zmax) zmax = XYZ(numUdp,3*ipnt+2);
    }

    ilin[nline] = -GR_DASHED;
    isym[nline] = +GR_CIRCLE;
    nper[nline] = npnt;
    nline++;

    /* build plot arrays for fit */
    status = EG_getRange(ecurve, trange, &periodic);
    CHECK_STATUS(EG_getRange);

    for (ipnt = 0; ipnt < 1000; ipnt++) {
        frac = (double)(ipnt) / (double)(1000-1);
        tt   = (1-frac) * trange[0] + frac * trange[1];

        status = EG_evaluate(ecurve, &tt, data);
        CHECK_STATUS(EG_evaluate);

        xplot[nplot] = data[0];
        yplot[nplot] = data[1];
        zplot[nplot] = data[2];
        nplot++;
    }

    ilin[nline] = +GR_SOLID;
    isym[nline] = -GR_PLUS;
    nper[nline] = 1000;
    nline++;

    /* build plot arrays for control points */
    status = EG_getGeometry(ecurve, &oclass, &mtype, &eref, &tempIlist, &tempRlist);
    CHECK_STATUS(EG_getGeometry);

    if (oclass != CURVE        ) goto cleanup;
    if (mtype  != BSPLINE      ) goto cleanup;

    for (icp = 0; icp < ncp; icp++) {
        xplot[nplot] = tempRlist[tempIlist[3]+3*icp  ];
        yplot[nplot] = tempRlist[tempIlist[3]+3*icp+1];
        zplot[nplot] = tempRlist[tempIlist[3]+3*icp+2];
        nplot++;
    }

    ilin[nline] = +GR_DOTTED;
    isym[nline] = +GR_SQUARE;
    nper[nline] = ncp;
    nline++;

    /* generate plot */
    if        (xmin == xmax) {
        grinit_(&io_kbd, &io_scr, "udpFitcurve", strlen("udpFitcurve"));

        grline_(ilin, isym, &nline,                "~y~z~O=inputs, --=fit, ...=cp",
                &indgr, yplot, zplot, nper, strlen("~y~z~O=inputs, --=fit, ...=cp"));

    } else if (ymin == ymax) {
        grinit_(&io_kbd, &io_scr, "udpFitcurve", strlen("udpFitcurve"));

        grline_(ilin, isym, &nline,                "~z~x~O=inputs, --=fit, ...=cp",
                &indgr, zplot, xplot, nper, strlen("~z~x~O=inputs, --=fit, ...=cp"));

    } else {
        grinit_(&io_kbd, &io_scr, "udpFitcurve", strlen("udpFitcurve"));

        grline_(ilin, isym, &nline,                "~x~y~O=inputs, --=fit, ...=cp",
                &indgr, xplot, yplot, nper, strlen("~x~y~O=inputs, --=fit, ...=cp"));
    }

cleanup:
    if (xplot != NULL) free(xplot);
    if (yplot != NULL) free(yplot);
    if (zplot != NULL) free(zplot);

    return status;
}
#endif
