/*
 ************************************************************************
 *                                                                      *
 * Fitter -- best-fit cubic Bspline to cloud of points                  *
 *                                                                      *
 *           Written by John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 * Algorithms (without smoothing) documented in:                        *
 *    "The Creation of a Static BRep Model Given a Cloud of Points"     *
 *    John F. Dannenhoffer, III                                         *
 *    AIAA-2017-0138                                                    *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2012/2025  John F. Dannenhoffer, III (Syracuse University)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "Fitter.h"

#include "egads.h"

#define  EPS06      1.0e-06
#define  EPS10      1.0e-10
#define  EPS12      1.0e-12
#define  EPS14      1.0e-14
#define  MIN(A,B)   (((A) < (B)) ? (A) : (B))
#define  MAX(A,B)   (((A) < (B)) ? (B) : (A))
#define  SQR(A)     ((A) * (A))

#ifndef  DEBUG
    #define  DEBUG 0
#endif

#ifdef GRAFIC
    #include "grafic.h"
#endif

#define MAKE_NORMALS_PLOTFILE   0

/*
 ************************************************************************
 *                                                                      *
 *   structures                                                         *
 *                                                                      *
 ************************************************************************
 */

typedef struct {
    int       iedge;          /* Edge index */
    int       bitflag;        /* 1=ordered, 2=uPeriodic, 8=intGiven */

    int       m;              /* number of points in cloud */
    double    *XYZcloud;      /* array  of points in cloud (normalized) */
    double    *Tcloud;        /* (t)    of points in cloud */

    int       n;              /* number of control points */
    double    *cp;            /* array  of control points (normalized) */
    double    *srat;          /* stretching ratio between initial cp's */
    double    *norm;          /* normal direction */
    double    *f;             /* objective function components:
                                 3*(m)   distances between cloud and surface
                                +3*(n-2) relative smoothness measure */

    int       rtype;          /* restriction type (1-x, 2=y, 3=z) */
    int       nsmth;          /* number of post-smoothing passes */
    int       iter;           /* number of iterations so far */
    double    lambda;         /* LM switching parameter */

    double    scale;          /* scale factor of inputs */
    double    xavg;           /* average x    of inputs */
    double    yavg;           /* average y    of inputs */
    double    zavg;           /* average z    of inputs */

    FILE      *fp;            /* file pointer for outputs (or NULL) */
} fit1d_T;

typedef struct {
    int       iface;          /* Face index */
    int       bitflag;        /* 2=uPeriodic, 4=vPeriodic, 8=intGiven */

    int       m;              /* number of points in cloud */
    double    *XYZcloud;      /* array  of points in cloud */
    double    *UVcloud;       /* (u,v)  if points in cloud */

    int       nu;             /* number of control points in u direction */
    int       nv;             /* number of control points in v direction */
    double    *cp;            /* array  of control points */
    double    *norm;          /* array  of initial normals */
    double    *f;             /* objective function components:
                                 3*(m)           distances between cloud and surface
                                +3*(nu-2)*(nv-2) relative smoothness measure */

    int       rtype;          /* restriction type (1=x, 2=y, 3=z) */
    int       nsmth;          /* number of post-smoothing passes */
    int       iter;           /* number of iterations so far */
    double    lambda;         /* LM switching parameter */

    double    scale;          /* scale factor of inputs */
    double    xavg;           /* average x    of inputs */
    double    yavg;           /* average y    of inputs */
    double    zavg;           /* average z    of inputs */

    int       *mask;          /* mask for smoothing vectors */

    FILE      *fp;            /* file pointer for outputs (or NULL) */
} fit2d_T;

/*
 ************************************************************************
 *                                                                      *
 *   macros                                                             *
 *                                                                      *
 ************************************************************************
 */

/* macros for error checking */
#define ROUTINE(NAME) char routine[] = #NAME ;                          \
    if (routine[0] == '\0') printf("bad routine(%s)\n", routine);

#define CHECK_STATUS(X)                                                 \
    if (status < FIT_SUCCESS) {                                         \
        printf( "ERROR:: BAD STATUS = %d from %s (called from %s:%d)\n", status, #X, routine, __LINE__); \
        goto cleanup;                                                   \
    }

#define MALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR != NULL) {                                                  \
        printf("ERROR:: MALLOC overwrites for %s=%llx (called from %s:%d)\n", #PTR, (long long)PTR, routine, __LINE__); \
        status = FIT_MALLOC;                                            \
        goto cleanup;                                                   \
    }                                                                   \
    PTR = (TYPE *) malloc((SIZE) * sizeof(TYPE));                       \
    if (PTR == NULL) {                                                  \
        printf("ERROR:: MALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
        status = FIT_MALLOC;                                            \
        goto cleanup;                                                   \
    }

#define RALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR == NULL) {                                                  \
        MALLOC(PTR,TYPE,SIZE);                                          \
    } else {                                                            \
       realloc_temp = realloc(PTR, (SIZE) * sizeof(TYPE));              \
       if (PTR == NULL) {                                               \
           printf("ERROR:: RALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
           status = FIT_MALLOC;                                         \
           goto cleanup;                                                \
       } else {                                                         \
           PTR = (TYPE *)realloc_temp;                                  \
       }                                                                \
    }

#define FREE(PTR)                                               \
    if (PTR != NULL) {                                          \
        free(PTR);                                              \
    }                                                           \
    PTR = NULL;

//static void *realloc_temp=NULL;            /* used by RALLOC macro */

/*
 ************************************************************************
 *                                                                      *
 *   declarations                                                       *
 *                                                                      *
 ************************************************************************
 */

static int    fit1d_objf(fit1d_T *fit1d, double smooth, double Tcloud[],
                         double cp[], double srat[], double f[]);
static int    fit2d_objf(fit2d_T *fit2d, double smooth, double UVcloud[],
                         double cp[], double f[]);
static int    eval1dBspline(double T, int n, double cp[], double XYZ[],
                            /*@null@*/double dXYZdT[], /*@null@*/double dXYZdP[]);
static int    eval2dBspline(double U, double V, int nu, int nv, double P[], double XYZ[],
                            /*@null@*/double dXYZdU[], /*@null@*/double dXYZdV[], /*@null@*/double dXYZdP[]);
static int    cubicBsplineBases(int ncp, double t, double N[], double dN[]);
static int    interp1d(double T, int ntab, double Ttab[], double XYZtab[], double XYZ[]);
static double L2norm(double f[], int n);
static double Linorm(double f[], int n);
static int    fitPlane(double xyz[], int n, double *a, double *b, double *c, double *d);
static int    matsol(double A[], double b[], int n, double x[]);

#ifdef GRAFIC
    static void   plotCurve_image(int*, void*, void*, void*, void*, void*,
                                  void*, void*, void*, void*, void*, float*, char*, int);
    static void   plotSurface_image(int*, void*, void*, void*, void*, void*,
                                    void*, void*, void*, void*, void*, float*, char*, int);
#endif


/*
 ************************************************************************
 *                                                                      *
 *   fit1dCloud - find spline that best-fits the cloud of points        *
 *                                                                      *
 ************************************************************************
 */
int
fit1dCloud(int    iedge,                /* (in)  Edge index */
           int    m,                    /* (in)  number of points in cloud */
           int    bitflag,              /* (in)  1=ordered, 2=uPeriodic, 8=intGiven */
           double XYZcloud[],           /* (in)  array  of points in cloud (x,y,z,x,... ) */
           int    n,                    /* (in)  number of control points */
           double cp[],                 /* (in)  array  of control points (first and last are set) */
                                        /* (out) array  of control points (all set) */
           double smooth,               /* (in)  initial control net smoothing */
           int    rtype,                /* (in)  restriction type (1=x, 2=y, 3=z) */
           int    nsmth,                /* (in)  number of post-smoothing passes */
           double Tcloud[],             /* (out) T-parameters of points in cloud */
           double *normf,               /* (out) RMS of distances between cloud and fit */
           double *maxf,                /* (out) maximum distance between cloud and fit */
           double *dotmin,              /* (in)  minimum allowable normalized dot product of control polygon */
                                        /* (out) minimum actual    normalized dot product of control polygon */
 /*@null@*/int    *nmin,                /* (out) minimum number of cloud points in any interval */
           int    *numiter,             /* (in)  if >=0, number of iterations allowed */
                                        /* (out) number of iterations executed */
 /*@null@*/FILE   *fp)                  /* (in)  file for progress outputs (or NULL) */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    j, iter, accept, niter=100, ismth;
    double old_normf, old_maxf, dotallow, toler=EPS06;
    double xavg, yavg, zavg, *cpsmth=NULL;

    void    *myFit1d;
    fit1d_T *fit1d;

    ROUTINE(fit1dCloud);

    /* --------------------------------------------------------------- */

    if (fp != NULL) {
        fprintf(fp, "enter fit1dCloud(iedge=%d, bitflag=%d, m=%d, n=%d, smooth=%f, rtype=%d, nsmth=%d)\n",
                iedge, bitflag, m, n, smooth, rtype, nsmth);
    }

    if (m <= 1 || XYZcloud == NULL) {
        status = FIT_EMPTYCLOUD;
        goto cleanup;
    } else if (n <= 2 || cp == NULL) {
        status = FIT_NCP;
        goto cleanup;
    } else if (smooth < 0 || smooth > 1) {
        status = FIT_SMOOTH;
        goto cleanup;
    }

    if (*numiter >= 0) {
        niter = *numiter;
    }

    /* default returns */
    *normf   = EPS12;
    *maxf    = 0;
    if (nmin   != NULL) *nmin   = -1;
    *numiter = 0;

    /* minimum allowable normalized dot product */
    if (m >= 3) {
        dotallow = *dotmin;
    } else {
        dotallow = 0;
        *dotmin  = 1;
    }

    /* initialize */
    status = fit1d_init(m, bitflag, smooth, XYZcloud, n, cp, rtype, fp, normf, maxf, &myFit1d);
    CHECK_STATUS(fit1d_init);

    fit1d = (fit1d_T *)myFit1d;
    fit1d->iedge = iedge;
    fit1d->nsmth = nsmth;

    /* Levenberg-Marquardt iterations (if not satisfied with initial guess) */
    if (*normf > toler) {
        for (iter = 0; iter < niter; iter++) {
            (*numiter)++;

            old_normf = *normf;
            old_maxf  = *maxf;
            *dotmin   = dotallow;

            status = fit1d_step(myFit1d, smooth, normf, maxf, dotmin, &accept);
            CHECK_STATUS(fit1d_step);

            /* check for convergence */
            if (accept == 1 && fabs(*normf-old_normf) < toler && fabs(*maxf-old_maxf) < toler) {
                if (fp != NULL) {
                    fprintf(fp, "converged in %d itertions\n", iter);
                    fprintf(fp, "final     normf=%10.4e, maxf=%10.4e\n", *normf, *maxf);
                }
                break;
            } else if (accept == 0) {
                *normf = old_normf;
                *maxf  = old_maxf;
            }

            /* reduce smoothing for next step */
            smooth *= 0.99;
        }
    }

    if (DEBUG) {
        FILE *fp_plot;
        fp_plot = fopen("slugs.plot", "a");
        if (fp_plot != NULL) {
            fprintf(fp_plot, " %9d %9d edge_%d|c\n", n, 1, iedge);
            for (j = 0; j < n; j++) {
                fprintf(fp_plot, "%12.6f %12.6f %12.6f\n",
                        fit1d->scale * fit1d->cp[3*j  ] + fit1d->xavg,
                        fit1d->scale * fit1d->cp[3*j+1] + fit1d->yavg,
                        fit1d->scale * fit1d->cp[3*j+2] + fit1d->zavg);
            }

            fprintf(fp_plot, " %9d %9d norm_%d|g\n", 2, 1, iedge);

            xavg = 0;
            yavg = 0;
            zavg = 0;

            for (j = 0; j < m; j++) {
                xavg += XYZcloud[3*j  ];
                yavg += XYZcloud[3*j+1];
                zavg += XYZcloud[3*j+2];
            }

            xavg = (xavg/m);
            yavg = (yavg/m);
            zavg = (zavg/m);

            fprintf(fp_plot, "%12.6f %12.6f %12.6f  %12.6f %12.6f %12.6f\n",
                    xavg, yavg, zavg, xavg+fit1d->norm[0], yavg+fit1d->norm[1], zavg+fit1d->norm[2]);
            fclose(fp_plot);
        }
    }

    /* compute outputs, statistics, and clean up */
    status = fit1d_done(myFit1d, Tcloud, cp, normf, maxf, nmin);
    CHECK_STATUS(fit1d_done);

    /* smooth the control polygon */
    MALLOC(cpsmth, double, 3*n);

    for (ismth = 0; ismth < nsmth; ismth++) {
        printf("                 iedge %3d smoothing pass %2d\n", iedge, ismth);
        for (j = 1; j < n-1; j++) {
            cpsmth[3*(j)  ] = ( cp[3*(j-1)  ]
                              + cp[3*(j+1)  ]) / 2;
            cpsmth[3*(j)+1] = ( cp[3*(j-1)+1]
                              + cp[3*(j+1)+1]) / 2;
            cpsmth[3*(j)+2] = ( cp[3*(j-1)+2]
                              + cp[3*(j+1)+2]) / 2;
        }
        for (j = 1; j < n-1; j++) {
            cp[3*(j)  ] = ( cp[    3*(j)  ]
                          + cpsmth[3*(j)  ]) / 2;
            cp[3*(j)+1] = ( cp[    3*(j)+1]
                          + cpsmth[3*(j)+1]) / 2;
            cp[3*(j)+2] = ( cp[    3*(j)+2]
                          + cpsmth[3*(j)+2]) / 2;
        }
    }

    FREE(cpsmth);

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit1d_init - initialize the B-spline curve fitter                  *
 *                                                                      *
 ************************************************************************
 */
int
fit1d_init(int     m,                   /* (in)  number of points in cloud */
           int     bitflag,             /* (in)  1=ordered, 2=uPeriodic, 8=intGiven */
           double  smooth,              /* (in)  control net smoothing parameter */
           double  XYZcloud[],          /* (in)  array  of points in cloud (x,y,z,x,...) */
           int     n,                   /* (in)  number of control points */
           double  cp[],                /* (in)  array  of control points (first and last set) */
           int     rtype,               /* (in)  restriction type (1=x, 2=y, 3=z) */
 /*@null@*/FILE    *fp,                 /* (in)  file pointer for progress prints */
           double  *normf,              /* (out) initil RMS between cloud and fit */
           double  *maxf,               /* (out) maximum distance between cloud and fit */
           void    **context)           /* (out) pointer to context */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    ordered=0, uPeriodic=0, intGiven=0, nchange=1;
    int    nobj, j, k;
    double frac, xmin, xmax, ymin, ymax, zmin, zmax, del1, del2;
    double dbest, dtest, v0[3], v1[3], dot, a, b, c, d, worst;
    double xa, ya, za, xb, yb, zb, tt, xx, yy, zz;
    double dx0, dy0, dz0, dx1, dy1, dz1, ddotn, len0, len1, dotallow=0;

    fit1d_T *fit1d=NULL;

    ROUTINE(fit1d_init);

    /* --------------------------------------------------------------- */

    if ((bitflag & 1) != 0) ordered   = 1;
    if ((bitflag & 2) != 0) uPeriodic = 1;
    if ((bitflag & 8) != 0) intGiven  = 1;

    /* periodic is not implemented (yet) */
    if (uPeriodic == 1) {
        status = FIT_BITFLAG;
        goto cleanup;
    } else if (m <= 1 || XYZcloud == NULL) {
        status = FIT_EMPTYCLOUD;
        goto cleanup;
    } else if (n <= 2 || cp == NULL) {
        status = FIT_NCP;
        goto cleanup;
    }

    /* get a fit1d structure */
    MALLOC(fit1d, fit1d_T, 1);

    /* initialize the fit1d structure */
    fit1d->m        = m;
    fit1d->XYZcloud = NULL;
    fit1d->Tcloud   = NULL;

    fit1d->bitflag  = bitflag;

    fit1d->n        = n;
    fit1d->cp       = NULL;
    fit1d->srat     = NULL;
    fit1d->norm     = NULL;
    fit1d->f        = NULL;

    MALLOC(fit1d->XYZcloud, double, 3*fit1d->m);
    MALLOC(fit1d->Tcloud,   double,   fit1d->m);

    MALLOC(fit1d->cp,   double,            3*fit1d->n  );
    MALLOC(fit1d->srat, double,              fit1d->n  );
    MALLOC(fit1d->norm, double,            3           );
    MALLOC(fit1d->f,    double, 3*fit1d->m+3*fit1d->n-6);

    fit1d->rtype    = rtype;
    fit1d->nsmth    = 0;
    fit1d->iter     = 0;
    fit1d->lambda   = 1;

    fit1d->scale    = 1;
    fit1d->xavg     = 0;
    fit1d->yavg     = 0;
    fit1d->zavg     = 0;

    fit1d->fp       = fp;

    /* find extrema of input data */
    xmin = XYZcloud[0];
    xmax = XYZcloud[0];
    ymin = XYZcloud[1];
    ymax = XYZcloud[1];
    zmin = XYZcloud[2];
    zmax = XYZcloud[2];

    for (k = 0; k < fit1d->m; k++) {
        if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
        if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
        if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
        if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
        if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
        if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
    }

    if (cp[    0] < xmin) xmin = cp[    0];
    if (cp[    0] > xmax) xmax = cp[    0];
    if (cp[    1] < ymin) ymin = cp[    1];
    if (cp[    1] > ymax) ymax = cp[    1];
    if (cp[    2] < zmin) zmin = cp[    2];
    if (cp[    2] > zmax) zmax = cp[    2];

    if (cp[3*n-3] < xmin) xmin = cp[3*n-3];
    if (cp[3*n-3] > xmax) xmax = cp[3*n-3];
    if (cp[3*n-2] < ymin) ymin = cp[3*n-2];
    if (cp[3*n-2] > ymax) ymax = cp[3*n-2];
    if (cp[3*n-1] < zmin) zmin = cp[3*n-1];
    if (cp[3*n-1] > zmax) zmax = cp[3*n-1];

    if (intGiven == 1) {
        for (j = 1; j < fit1d->n-1; j++) {
            if (cp[3*j  ] < xmin) xmin = cp[3*j  ];
            if (cp[3*j  ] > xmax) xmax = cp[3*j  ];
            if (cp[3*j+1] < ymin) ymin = cp[3*j+1];
            if (cp[3*j+1] > ymax) ymax = cp[3*j+1];
            if (cp[3*j+2] < zmin) zmin = cp[3*j+2];
            if (cp[3*j+2] > zmax) zmax = cp[3*j+2];
        }
    }

    /* normalize input data */
                                  fit1d->scale = xmax - xmin;
    if (ymax-ymin > fit1d->scale) fit1d->scale = ymax - ymin;
    if (zmax-zmin > fit1d->scale) fit1d->scale = zmax - zmin;

    if (fit1d->scale < EPS12) {
        printf("%20.12e %20.12e %20.12e\n", xmin, ymin, zmin);
        printf("%20.12e %20.12e %20.12e\n", xmax, ymax, zmax);
        status = FIT_DEGENERATE;
        goto cleanup;
    }

    fit1d->xavg = (xmax + xmin) / 2;
    fit1d->yavg = (ymax + ymin) / 2;
    fit1d->zavg = (zmax + zmin) / 2;

    for (k = 0; k < fit1d->m; k++) {
        fit1d->XYZcloud[3*k  ] = (XYZcloud[3*k  ] - fit1d->xavg) / fit1d->scale;
        fit1d->XYZcloud[3*k+1] = (XYZcloud[3*k+1] - fit1d->yavg) / fit1d->scale;
        fit1d->XYZcloud[3*k+2] = (XYZcloud[3*k+2] - fit1d->zavg) / fit1d->scale;
    }

    fit1d->cp[    0] = (cp[    0] - fit1d->xavg) / fit1d->scale;
    fit1d->cp[    1] = (cp[    1] - fit1d->yavg) / fit1d->scale;
    fit1d->cp[    2] = (cp[    2] - fit1d->zavg) / fit1d->scale;

    fit1d->cp[3*n-3] = (cp[3*n-3] - fit1d->xavg) / fit1d->scale;
    fit1d->cp[3*n-2] = (cp[3*n-2] - fit1d->yavg) / fit1d->scale;
    fit1d->cp[3*n-1] = (cp[3*n-1] - fit1d->zavg) / fit1d->scale;

    if (intGiven == 1) {
        for (j = 1; j < fit1d->n-1; j++) {
            fit1d->cp[3*j  ] = (cp[3*j  ] - fit1d->xavg) / fit1d->scale;
            fit1d->cp[3*j+1] = (cp[3*j+1] - fit1d->yavg) / fit1d->scale;
            fit1d->cp[3*j+2] = (cp[3*j+2] - fit1d->zavg) / fit1d->scale;
        }
    }

    /* number of design variables, objectives, and interior control points*3 */
    nobj = 3 * fit1d->m + 3 * (fit1d->n - 2);

    /* initialize the spacing ratios */
    for (j = 0; j < fit1d->n; j++) {
        fit1d->srat[j] = 0;
    }

    if (intGiven == 1) {
        /* determine the spacing ratio, which will be used to hold cp
           spacing if interior points are given */
        for (j = 1; j < fit1d->n-1; j++) {
            del1 = SQR(cp[3*j-3] - cp[3*j  ]) + SQR(cp[3*j-2] - cp[3*j+1]) + SQR(cp[3*j-1] - cp[3*j+2]);
            del2 = SQR(cp[3*j+3] - cp[3*j  ]) + SQR(cp[3*j+4] - cp[3*j+1]) + SQR(cp[3*j+5] - cp[3*j+2]);
            fit1d->srat[j] = 0.5 - del1 / (del1 + del2);
        }

        /* for each point in the cloud, assign the value of t
           that is associated with the closest control point */
        for (k = 0; k < fit1d->m; k++) {
            fit1d->Tcloud[k] = 0;
            dbest = 1e20;

            for (j = 0; j < fit1d->n; j++) {
                dtest = SQR(XYZcloud[3*k  ] - cp[3*j  ])
                      + SQR(XYZcloud[3*k+1] - cp[3*j+1])
                      + SQR(XYZcloud[3*k+2] - cp[3*j+2]);

                if (dtest < dbest) {
                    frac = (double)(j) / (double)(fit1d->n-1);
                    v0[0] = XYZcloud[3*k  ] - cp[3*j  ];
                    v0[1] = XYZcloud[3*k+1] - cp[3*j+1];
                    v0[2] = XYZcloud[3*k+2] - cp[3*j+2];

                    if (j > 0) {
                        v1[0] = cp[3*j-3] - cp[3*j  ];
                        v1[1] = cp[3*j-2] - cp[3*j+1];
                        v1[2] = cp[3*j-1] - cp[3*j+2];
                        dot   = v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
                        if (dot > 0) {
                            frac -= MIN(0.5, dot);
                        }
                    }
                    if (j < fit1d->n-1) {
                        v1[0] = cp[3*j+4] - cp[3*j  ];
                        v1[1] = cp[3*j+5] - cp[3*j+1];
                        v1[2] = cp[3*j+6] - cp[3*j+2];
                        dot   = v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
                        if (dot > 0) {
                            frac += MIN(0.5, dot);
                        }
                    }
                    fit1d->Tcloud[k] = frac * (double)(fit1d->n - 3);

                    dbest = dtest;
                }
            }
        }

    /* if fit1d->m < 3, then assume that the linear spline is the best fit */
    } else if (fit1d->m < 3) {
        for (j = 1; j < fit1d->n-1; j++) {
            frac = (double)(j) / (double)(fit1d->n-1);

            fit1d->cp[3*j  ] = (1-frac) * fit1d->cp[0] + frac * fit1d->cp[3*fit1d->n-3];
            fit1d->cp[3*j+1] = (1-frac) * fit1d->cp[1] + frac * fit1d->cp[3*fit1d->n-2];
            fit1d->cp[3*j+2] = (1-frac) * fit1d->cp[2] + frac * fit1d->cp[3*fit1d->n-1];
        }

        xa = fit1d->cp[           0];
        ya = fit1d->cp[           1];
        za = fit1d->cp[           2];
        xb = fit1d->cp[3*fit1d->n-3];
        yb = fit1d->cp[3*fit1d->n-2];
        zb = fit1d->cp[3*fit1d->n-1];

        for (k = 0; k < fit1d->m; k++) {
            xx = fit1d->XYZcloud[3*k  ];
            yy = fit1d->XYZcloud[3*k+1];
            zz = fit1d->XYZcloud[3*k+2];

            tt = ((xx-xa) * (xb-xa) + (yy-ya) * (yb-ya) + (zz-za) * (zb-za))
               / ((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya) + (zb-za) * (zb-za));

            fit1d->Tcloud[k] = tt * (fit1d->n-3);
        }

        if (fit1d->fp != NULL) {
            fprintf(fit1d->fp, "making linear fit because not enough points in cloud\n");
        }

    } else {

        /* ordered == 0 -> linear interpolation of control points from boundries */
        if (ordered  == 0) {
            for (j = 1; j < fit1d->n-1; j++) {
                frac = (double)(j) / (double)(fit1d->n-1);

                fit1d->cp[3*j  ] = (1-frac) * fit1d->XYZcloud[0] + frac * fit1d->XYZcloud[3*fit1d->m-3];
                fit1d->cp[3*j+1] = (1-frac) * fit1d->XYZcloud[1] + frac * fit1d->XYZcloud[3*fit1d->m-2];
                fit1d->cp[3*j+2] = (1-frac) * fit1d->XYZcloud[2] + frac * fit1d->XYZcloud[3*fit1d->m-1];
            }

            for (k = 0; k < fit1d->m; k++) {
                fit1d->Tcloud[k] = (double)(k) / (double)(fit1d->m - 1) * (double)(fit1d->n - 3);
            }

        /* ordered == 1 -> equi-arclength spacing of control points */
        } else {

            /* compute total arclength */
            fit1d->Tcloud[0] = 0;
            for (k = 1; k < fit1d->m; k++) {
                fit1d->Tcloud[k] = fit1d->Tcloud[k-1]
                    + sqrt( SQR(fit1d->XYZcloud[3*k  ] - fit1d->XYZcloud[3*k-3])
                           +SQR(fit1d->XYZcloud[3*k+1] - fit1d->XYZcloud[3*k-2])
                           +SQR(fit1d->XYZcloud[3*k+2] - fit1d->XYZcloud[3*k-1]));
            }

            /* scale arc-length to range 0 to fit1d->n-3 */
            for (k = 0; k < fit1d->m; k++) {
                fit1d->Tcloud[k] *= (fit1d->n - 3) / fit1d->Tcloud[fit1d->m-1];
            }

            /* now interpolate for the interior control point locations */
            for (j = 1; j < fit1d->n-1; j++) {
                frac = (double)(j) / (double)(fit1d->n - 1) * (double)(fit1d->n - 3);
                status = interp1d(frac, fit1d->m, fit1d->Tcloud, fit1d->XYZcloud, &(fit1d->cp[3*j]));
                CHECK_STATUS(interp1d);
            }
        }
    }

    /* find best-fit plane to the data */
    if (fit1d->m >= 3) {
        status = fitPlane(fit1d->XYZcloud, fit1d->m, &a, &b, &c, &d);
        CHECK_STATUS(fit1d);

        if        (fabs(a) >= fabs(b) && fabs(a) >= fabs(c)) {
            fit1d->norm[0] = 1;
            fit1d->norm[1] = b / a;
            fit1d->norm[2] = c / a;
        } else if (fabs(b) >= fabs(c) && fabs(b) >= fabs(a)) {
            fit1d->norm[0] = a / b;
            fit1d->norm[1] = 1;
            fit1d->norm[2] = c / b;
        } else {
            fit1d->norm[0] = a / c;
            fit1d->norm[1] = b / c;
            fit1d->norm[2] = 1;
        }
    } else {
        fit1d->norm[0] = 0;
        fit1d->norm[1] = 0;
        fit1d->norm[2] = 1;
    }

    if (DEBUG) {
        printf("norm=%12.6f, %12.6f, %12.6f\n",
               fit1d->norm[0],
               fit1d->norm[1],
               fit1d->norm[2]);
    }

    /* make sure initial control net does not have any very acute angles */
    while (nchange > 0) {
        nchange = 0;
        worst   = 1;

        for (j = 1; j < fit1d->n-1; j++) {
            dx0 = fit1d->cp[3*j  ] - fit1d->cp[3*j-3];
            dy0 = fit1d->cp[3*j+1] - fit1d->cp[3*j-2];
            dz0 = fit1d->cp[3*j+2] - fit1d->cp[3*j-1];

            ddotn = dx0 * fit1d->norm[0] + dy0 * fit1d->norm[1] + dz0 * fit1d->norm[2];
            dx0  -= ddotn * fit1d->norm[0];
            dy0  -= ddotn * fit1d->norm[1];
            dz0  -= ddotn * fit1d->norm[2];

            dx1 = fit1d->cp[3*j+3] - fit1d->cp[3*j  ];
            dy1 = fit1d->cp[3*j+4] - fit1d->cp[3*j+1];
            dz1 = fit1d->cp[3*j+5] - fit1d->cp[3*j+2];

            ddotn = dx1 * fit1d->norm[0] + dy1 * fit1d->norm[1] + dz1 * fit1d->norm[2];
            dx1  -= ddotn * fit1d->norm[0];
            dy1  -= ddotn * fit1d->norm[1];
            dz1  -= ddotn * fit1d->norm[2];

            len0 = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
            len1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

            if (len0 < EPS12 || len1 < EPS12) {
                printf("j=%2d, len0=%12.4e, len1=%12.4e\n", j, len0, len1);
            } else {
                dot = (dx0 * dx1 + dy0 * dy1 + dz0 * dz1) / len0 / len1;

                worst = MIN(worst, dot);
                if (dot < dotallow) {
                    fit1d->cp[3*j  ] = (fit1d->cp[3*j-3] + 2 * fit1d->cp[3*j  ] + fit1d->cp[3*j+3]) / 4;
                    fit1d->cp[3*j+1] = (fit1d->cp[3*j-2] + 2 * fit1d->cp[3*j+1] + fit1d->cp[3*j+4]) / 4;
                    fit1d->cp[3*j+2] = (fit1d->cp[3*j-1] + 2 * fit1d->cp[3*j+2] + fit1d->cp[3*j+5]) / 4;
                    nchange++;
                }
            }
        }
    }

    /* compute the initial objective function */
    status = fit1d_objf(fit1d, smooth, fit1d->Tcloud, fit1d->cp, fit1d->srat, fit1d->f);
    CHECK_STATUS(fit1d_objf);

    /* compute and report the norm of the objective function */
    *maxf  = Linorm(fit1d->f, 3*fit1d->m);
    *normf = L2norm(fit1d->f, nobj) / sqrt(nobj);
    if (fit1d->fp != NULL) {
        fprintf(fit1d->fp, "initial   normf=%10.4e, maxf=%10.4e\n", *normf, *maxf);
    }

    /* return the fit1d structure (as the context) */
    *context = (void *)fit1d;

cleanup:
    if (status != FIT_SUCCESS && fit1d != NULL) {
        FREE(fit1d->XYZcloud);
        FREE(fit1d->Tcloud  );
        FREE(fit1d->cp      );
        FREE(fit1d->srat    );
        FREE(fit1d->norm    );
        FREE(fit1d->f       );

        FREE(fit1d);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit1d_step - take one step of the B-spline curve fitter            *
 *                                                                      *
 ************************************************************************
 */
int
fit1d_step(void    *context,            /* (in)  context */
           double  smooth,              /* (in)  control net smoothing parameter */
           double  *normf,              /* (out) RMS of distances between cloud and fit */
           double  *maxf,               /* (out) maximum distance between cloud and fit */
           double  *dotmin,             /* (in)  minimum allowable normalized dot product of control polygon */
                                        /* (out) minimum actual    normalized dot product of control polygon */
           int     *accept)             /* (out) =1 if step was accepted */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    nvar, ivar, jvar, nobj, iobj, i, j, k, next, uPeriodic=0;

    double normfnew, maxfnew, tempw3=0, tempw2=0, tempw1=0, tempc, tempe1, tempe2, tempe3;
    double delta0, delta1, delta2, dotallow, dx0, dy0, dz0, dx1, dy1, dz1, ddotn, len0, len1, dot;
    double XYZ[3], dXYZdT[3], *dXYZdP=NULL;
    double *cpnew=NULL, *beta=NULL, *betanew=NULL, *delta=NULL;
    double *fnew=NULL;

    fit1d_T *fit1d = (fit1d_T *) context;

    int    nn = 3 * (fit1d->n - 2);
    double *AA=NULL, *BB=NULL, *CC=NULL, *DD=NULL, *EE=NULL;

#define A(K)      AA[(K)]
#define B(K,I)    BB[(K)*nn+(I)]
#define C(I,J)    CC[(I)*nn+(J)]
#define D(K)      DD[(K)]
#define E(K)      EE[(K)]

    ROUTINE(fit1d_step);

    /* --------------------------------------------------------------- */

    *accept = 0;

    /* note: first 5 iterations only adjust U */
    (fit1d->iter)++;

    /* number of design variables, objectives, and interior control points*3 */
    nvar =     fit1d->m + 3 * (fit1d->n - 2);
    nobj = 3 * fit1d->m + 3 * (fit1d->n - 2);

    /* allocate all temporary arrays */
    MALLOC(dXYZdP,  double,   fit1d->n);
    MALLOC(cpnew,   double, 3*fit1d->n);

    MALLOC(beta,    double, nvar);
    MALLOC(delta,   double, nvar);
    MALLOC(betanew, double, nvar);

    MALLOC(fnew,    double, nobj);

    /* store the t values in the first m betas */
    next = 0;
    for (k = 0; k < fit1d->m; k++) {
        beta[next++] = fit1d->Tcloud[k];
    }

    /* store the interior control points into beta */
    for (j = 1; j < fit1d->n-1; j++) {
        beta[next++] = fit1d->cp[3*j  ];
        beta[next++] = fit1d->cp[3*j+1];
        beta[next++] = fit1d->cp[3*j+2];
    }
    assert(next == nvar);

    /* minimum allowable dot product */
    dotallow = *dotmin;
    *dotmin  = 1;

    /*
          trans(J) * J =  [   A      B ]      trans(J) * Q =  [ D ]
                          [            ]                      [   ]
                          [trans(B)  C ]                      [ E ]
    */
    MALLOC(AA, double, fit1d->m   );
    MALLOC(BB, double, fit1d->m*nn);
    MALLOC(CC, double, nn      *nn);
    MALLOC(DD, double, fit1d->m   );
    MALLOC(EE, double,          nn);

    for (i = 1; i < fit1d->n-1; i++) {
        for (j = 1; j < fit1d->n-1; j++) {
            C(3*i-3,3*j-3) = 0;
            C(3*i-2,3*j-3) = 0;
            C(3*i-1,3*j-3) = 0;

            C(3*i-3,3*j-2) = 0;
            C(3*i-2,3*j-2) = 0;
            C(3*i-1,3*j-2) = 0;

            C(3*i-3,3*j-1) = 0;
            C(3*i-2,3*j-1) = 0;
            C(3*i-1,3*j-1) = 0;
        }

        E(3*i-3) = 0;
        E(3*i-2) = 0;
        E(3*i-1) = 0;
    }

    /* create top-left (A) and top-right (B) parts of JtJ and top part (D) of JtQ */
    for (k = 0; k < fit1d->m; k++) {
        status = eval1dBspline(beta[k], fit1d->n, fit1d->cp, XYZ, dXYZdT, dXYZdP);
        CHECK_STATUS(eval1dBspline);

        /* top-left (A) part of JtJ */
        A(k) = dXYZdT[0] * dXYZdT[0] + dXYZdT[1] * dXYZdT[1] + dXYZdT[2] * dXYZdT[2];

        /* top-right (B) part of JtJ */
        for (j = 1; j < fit1d->n-1; j++) {
            B(k,3*j-3) = dXYZdT[0] * dXYZdP[j];
            B(k,3*j-2) = dXYZdT[1] * dXYZdP[j];
            B(k,3*j-1) = dXYZdT[2] * dXYZdP[j];
        }

        /* top part (D) of JtQ (negative needed since f = (XYZ_spline - XYZ_cloud) */
        D(k) = - dXYZdT[0] * fit1d->f[3*k] - dXYZdT[1] * fit1d->f[3*k+1] - dXYZdT[2] * fit1d->f[3*k+2];

        for (ivar = 1; ivar < fit1d->n-1; ivar++) {
            E(3*ivar-3) -= dXYZdP[ivar] * fit1d->f[3*k  ];
            E(3*ivar-2) -= dXYZdP[ivar] * fit1d->f[3*k+1];
            E(3*ivar-1) -= dXYZdP[ivar] * fit1d->f[3*k+2];
        }
    }

    /* create bottom-right (C) part of JtJ */
    for (j = 1; j < fit1d->n-1; j++) {
        jvar = 3 * (j - 1);

        if        (j == 1) {
            tempc  =  5;
            tempe1 = -4;
            tempe2 =  1;
            tempe3 =  0;
        } else if (j == 2) {
            tempw1 = -4;
            tempc  =  6;
            tempe1 = -4;
            tempe2 =  1;
            tempe3 =  0;
        } else if (j == fit1d->n-2) {
            tempw3  = 0;
            tempw2 =  1;
            tempw1 = -4;
            tempc  =  5;
        } else if (j == fit1d->n-3) {
            tempw3  = 0;
            tempw2 =  1;
            tempw1 = -4;
            tempc  =  6;
            tempe1 = -4;
        } else {
            tempw3 =  0;
            tempw2 =  1;
            tempw1 = -4;
            tempc  =  6;
            tempe1 = -4;
            tempe2 =  1;
            tempe3 =  0;
        }

        for (k = 0; k < fit1d->m; k++) {
            status = eval1dBspline(beta[k], fit1d->n, fit1d->cp, XYZ, dXYZdT, dXYZdP);
            CHECK_STATUS(eval1dBspline);

            if (j > 3         ) tempw3 += dXYZdP[j] * dXYZdP[j-3];
            if (j > 2         ) tempw2 += dXYZdP[j] * dXYZdP[j-2];
            if (j > 1         ) tempw1 += dXYZdP[j] * dXYZdP[j-1];
            if (1             ) tempc  += dXYZdP[j] * dXYZdP[j  ];
            if (j < fit1d->n-2) tempe1 += dXYZdP[j] * dXYZdP[j+1];
            if (j < fit1d->n-3) tempe2 += dXYZdP[j] * dXYZdP[j+2];
            if (j < fit1d->n-4) tempe3 += dXYZdP[j] * dXYZdP[j+3];
        }

        if (j > 3   ) {
            C(jvar,  jvar- 9) = tempw3;
            C(jvar+1,jvar- 8) = tempw3;
            C(jvar+2,jvar- 7) = tempw3;
        }

        if (j > 2   ) {
            C(jvar,  jvar- 6) = tempw2;
            C(jvar+1,jvar- 5) = tempw2;
            C(jvar+2,jvar- 4) = tempw2;
        }

        if (j > 1   ) {
            C(jvar,  jvar- 3) = tempw1;
            C(jvar+1,jvar- 2) = tempw1;
            C(jvar+2,jvar- 1) = tempw1;
        }

        if (1) {
            C(jvar,  jvar   ) = tempc ;
            C(jvar+1,jvar+ 1) = tempc ;
            C(jvar+2,jvar+ 2) = tempc ;
        }

        if (j < fit1d->n-2) {
            C(jvar  ,jvar+ 3) = tempe1;
            C(jvar+1,jvar+ 4) = tempe1;
            C(jvar+2,jvar+ 5) = tempe1;
        }

        if (j < fit1d->n-3) {
            C(jvar,  jvar+ 6) = tempe2;
            C(jvar+1,jvar+ 7) = tempe2;
            C(jvar+2,jvar+ 8) = tempe2;
        }

        if (j < fit1d->n-4) {
            C(jvar,  jvar+ 9) = tempe3;
            C(jvar+1,jvar+10) = tempe3;
            C(jvar+2,jvar+11) = tempe3;
        }
    }

    /* create bottom (E) part of JtQ */
    for (j = 1; j < fit1d->n-1; j++) {
        ivar = 3 * fit1d->m + 3 * j - 3;

        E(3*j-3) -= 2 * smooth * fit1d->f[ivar  ];
        E(3*j-2) -= 2 * smooth * fit1d->f[ivar+1];
        E(3*j-1) -= 2 * smooth * fit1d->f[ivar+2];

        if (j > 1) {
            E(3*j-3) += smooth * fit1d->f[ivar-3];
            E(3*j-2) += smooth * fit1d->f[ivar-2];
            E(3*j-1) += smooth * fit1d->f[ivar-1];
        }

        if (j < fit1d->n-2) {
            E(3*j-3) += smooth * fit1d->f[ivar+3];
            E(3*j-2) += smooth * fit1d->f[ivar+4];
            E(3*j-1) += smooth * fit1d->f[ivar+5];
        }
    }

    /* multiply diagonals of JtJ by (1 + lambda)  */
    for (k = 0; k < fit1d->m; k++) {
        A(k) *= (1 + fit1d->lambda);
    }

    for (j = 1; j < fit1d->n-1; j++) {
        C(3*j-3,3*j-3) *= (1 + fit1d->lambda);
        C(3*j-2,3*j-2) *= (1 + fit1d->lambda);
        C(3*j-1,3*j-1) *= (1 + fit1d->lambda);
    }

    /* update:  C = C - trans(B) * inv(A) * B
       update:  E = E - trans(B) * inv(A) * D  */
    for (k = 0; k < fit1d->m; k++) {
        for (i = 1; i < fit1d->n-1; i++) {
            for (j = 1; j < fit1d->n-1; j++) {
                C(3*i-3,3*j-3) -= B(k,3*i-3) / A(k) * B(k,3*j-3);
                C(3*i-2,3*j-3) -= B(k,3*i-2) / A(k) * B(k,3*j-3);
                C(3*i-1,3*j-3) -= B(k,3*i-1) / A(k) * B(k,3*j-3);

                C(3*i-3,3*j-2) -= B(k,3*i-3) / A(k) * B(k,3*j-2);
                C(3*i-2,3*j-2) -= B(k,3*i-2) / A(k) * B(k,3*j-2);
                C(3*i-1,3*j-2) -= B(k,3*i-1) / A(k) * B(k,3*j-2);

                C(3*i-3,3*j-1) -= B(k,3*i-3) / A(k) * B(k,3*j-1);
                C(3*i-2,3*j-1) -= B(k,3*i-2) / A(k) * B(k,3*j-1);
                C(3*i-1,3*j-1) -= B(k,3*i-1) / A(k) * B(k,3*j-1);
            }

            E(3*i-3) -= B(k,3*i-3) / A(k) * D(k);
            E(3*i-2) -= B(k,3*i-2) / A(k) * D(k);
            E(3*i-1) -= B(k,3*i-1) / A(k) * D(k);
        }
    }

    /* solve for the second part of beta (the control points) */
    status = matsol(CC, EE, nn, &(delta[fit1d->m]));
    CHECK_STATUS(matsol);

    /* solve for the first part of beta (the parametric coordinates) */
    for (k = 0; k < fit1d->m; k++) {
        delta[k] = D(k);

        for (i = 1; i < fit1d->n-1;i++) {
            delta[k] -= B(k,3*i-3) * delta[fit1d->m+3*i-3]
                      + B(k,3*i-2) * delta[fit1d->m+3*i-2]
                      + B(k,3*i-1) * delta[fit1d->m+3*i-1];
        }

        delta[k] /= A(k);
    }

    /* find the temporary new beta (and clip the Tclouds) */
    for (ivar = 0; ivar < nvar; ivar++) {
        betanew[ivar] = beta[ivar] + delta[ivar];
        if (ivar < fit1d->m) {
            if (betanew[ivar] < 0         ) betanew[ivar] = 0;
            if (betanew[ivar] > fit1d->n-3) betanew[ivar] = fit1d->n-3;
        }
    }

    /* extract the temporary control points from betanew */
    next = fit1d->m;
    for (j = 0; j < fit1d->n; j++) {
        if (j == 0 || j == fit1d->n-1) {
            cpnew[3*j  ] = fit1d->cp[3*j  ];
            cpnew[3*j+1] = fit1d->cp[3*j+1];
            cpnew[3*j+2] = fit1d->cp[3*j+2];
        } else if (fit1d->iter <= 5) {
            cpnew[3*j  ] = fit1d->cp[3*j  ];
            cpnew[3*j+1] = fit1d->cp[3*j+1];
            cpnew[3*j+2] = fit1d->cp[3*j+2];
            next += 3;
        } else {
            cpnew[3*j  ] = betanew[next++];
            cpnew[3*j+1] = betanew[next++];
            cpnew[3*j+2] = betanew[next++];
        }
    }
    assert (next == nvar);

    /* reject changes if it exceeds the allowable minimum
       normalized dot product or if control points get too
       close to each other */
    for (j = 0; j < fit1d->n; j++) {
        if (j == 0 || j == fit1d->n-1) {
            continue;
        } else if (fit1d->iter <= 5) {
            continue;
        } else {
            dx0 = cpnew[3*j  ] - cpnew[3*j-3];
            dy0 = cpnew[3*j+1] - cpnew[3*j-2];
            dz0 = cpnew[3*j+2] - cpnew[3*j-1];

            ddotn = dx0 * fit1d->norm[0] + dy0 * fit1d->norm[1] + dz0 * fit1d->norm[2];
            dx0  -= ddotn * fit1d->norm[0];
            dy0  -= ddotn * fit1d->norm[1];
            dz0  -= ddotn * fit1d->norm[2];

            dx1 = cpnew[3*j+3] - cpnew[3*j  ];
            dy1 = cpnew[3*j+4] - cpnew[3*j+1];
            dz1 = cpnew[3*j+5] - cpnew[3*j+2];

            ddotn = dx1 * fit1d->norm[0] + dy1 * fit1d->norm[1] + dz1 * fit1d->norm[2];
            dx1  -= ddotn * fit1d->norm[0];
            dy1  -= ddotn * fit1d->norm[1];
            dz1  -= ddotn * fit1d->norm[2];

            len0 = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
            len1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);

            if (len0 < EPS12 || len1 < EPS12) {
                cpnew[3*j  ] = fit1d->cp[3*j  ];
                cpnew[3*j+1] = fit1d->cp[3*j+1];
                cpnew[3*j+2] = fit1d->cp[3*j+2];
                continue;
            }

            dot = (dx0 * dx1 + dy0 * dy1 + dz0 * dz1) / len0 / len1;
            if (dot < dotallow) {
                cpnew[3*j  ] = fit1d->cp[3*j  ];
                cpnew[3*j+1] = fit1d->cp[3*j+1];
                cpnew[3*j+2] = fit1d->cp[3*j+2];
            } else {
                *dotmin = MIN(*dotmin, dot);
            }
        }
    }

    /* apply periodicity condition by making sure first and last
       intervals are the same */
#ifndef __clang_analyzer__
    if ((fit1d->bitflag & 2) != 0) uPeriodic = 1;

    if (uPeriodic == 1) {
        delta0 = (2*cpnew[0] - cpnew[3] - cpnew[3*fit1d->n-6]) / 2;
        delta1 = (2*cpnew[1] - cpnew[4] - cpnew[3*fit1d->n-5]) / 2;
        delta2 = (2*cpnew[2] - cpnew[5] - cpnew[3*fit1d->n-4]) / 2;

        cpnew[           3] += delta0;
        cpnew[           4] += delta1;
        cpnew[           5] += delta2;

        cpnew[3*fit1d->n-6] += delta0;
        cpnew[3*fit1d->n-5] += delta1;
        cpnew[3*fit1d->n-4] += delta2;
    }
#endif

    /* compute the objective function based upon the new beta */
    status = fit1d_objf(fit1d, smooth, betanew, cpnew, fit1d->srat, fnew);
    CHECK_STATUS(fit1d_obj);

    maxfnew  = Linorm(fnew, 3*fit1d->m);
    normfnew = L2norm(fnew, nobj) / sqrt(nobj);
    if (fit1d->iter%10 == 0 && fit1d->fp != NULL) {
        fprintf(fit1d->fp, "iter=%4d normf=%10.4e, maxf=%10.4e, dotmin=%7.4f  ", fit1d->iter, normfnew, maxfnew, *dotmin);
    }

    /* if this was a better step, accept it and decrease
       lambda (making it more Newton-like) */
    if (normfnew < *normf) {
        fit1d->lambda = MAX(fit1d->lambda/2, EPS10);
        if (fit1d->iter%10 == 0 && fit1d->fp != NULL) {
            fprintf(fit1d->fp, "ACCEPTED,  lambda=%10.3e,  smooth=%10.3e\n", fit1d->lambda, smooth);
        }

        /* save new design variables, control points, and objective function */
        for (k = 0; k < fit1d->m; k++) {
            fit1d->Tcloud[k] = betanew[k];
        }
        for (j = 0; j < fit1d->n; j++) {
            fit1d->cp[3*j  ] = cpnew[3*j  ];
            fit1d->cp[3*j+1] = cpnew[3*j+1];
            fit1d->cp[3*j+2] = cpnew[3*j+2];
        }
#ifndef __clang_analyzer__
        for (iobj = 0; iobj < nobj; iobj++) {
            fit1d->f[iobj] = fnew[iobj];
        }
#endif
        *normf  = normfnew;
        *maxf   = maxfnew;
        *accept = 1;

    /* otherwise do not take the step and increase lambda (making it
       more steepest-descent-like) */
    } else {
        fit1d->lambda = MIN(fit1d->lambda*2, 1.0e+10);
        if (fit1d->iter%10 == 0 && fit1d->fp != NULL) {
            fprintf(fit1d->fp, "rejected,  lambda=%10.3e,  smooth=%10.3e\n", fit1d->lambda, smooth);
        }
    }

cleanup:
    FREE(fnew   );
    FREE(betanew);
    FREE(delta  );
    FREE(beta   );
    FREE(cpnew  );
    FREE(dXYZdP );

    FREE(AA);
    FREE(BB);
    FREE(CC);
    FREE(DD);
    FREE(EE);
#undef A
#undef B
#undef C
#undef D
#undef E

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit1d_done - clean up after B-spline curve fitter                  *
 *                                                                      *
 ************************************************************************
 */
int
fit1d_done(void    *context,            /* (in)  context */
           double  Tcloud[],            /* (out) T-parameters of points in cloud */
           double  cp[],                /* (out) array of control points */
           double  *normf,              /* (out) RMS of distances between cloud and fit */
           double  *maxf,               /* (out) maximum distance between cloud and fit */
 /*@null@*/int     *nmin)               /* (out) minimum number of cloud points in any interval (or NULL) */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    j, k;
    int    *nper=NULL;

    fit1d_T *fit1d = (fit1d_T *) context;

    ROUTINE(fit1d_done);

    /* --------------------------------------------------------------- */

    /* allocate temporry array */
    MALLOC(nper, int, fit1d->n);

    /* find the minimum number of cloud points in each spline interval */
    if (nmin != NULL) {
        *nmin = fit1d->m;

        for (j = 0; j < fit1d->n-3; j++) {
            nper[j] = 0;
        }

        for (k = 0; k < fit1d->m; k++) {
            j = MIN(MAX(0,floor(fit1d->Tcloud[k])), fit1d->n-4);
            nper[j]++;
        }

        for (j = 0; j < fit1d->n-3; j++) {
            if (nper[j] < *nmin) {
                *nmin = nper[j];
            }
        }
    }

    /* extract the T parameters from the structure */
    for (k = 0; k < fit1d->m; k++) {
        Tcloud[k] = fit1d->Tcloud[k];
    }

    /* extract the control points from the structure (and un-normalized them) */
    for (j = 0; j < fit1d->n; j++) {
        cp[3*j  ] = fit1d->scale * fit1d->cp[3*j  ] + fit1d->xavg;
        cp[3*j+1] = fit1d->scale * fit1d->cp[3*j+1] + fit1d->yavg;
        cp[3*j+2] = fit1d->scale * fit1d->cp[3*j+2] + fit1d->zavg;
    }

    /* print the dot products on the interior of the control net */

    /* return the norm of the error (un-normlized) */
    *normf *= fit1d->scale;
    *maxf  *= fit1d->scale;

    /* remove temporary storage from fit1d structure */
    FREE(fit1d->XYZcloud);
    FREE(fit1d->Tcloud  );
    FREE(fit1d->cp      );
    FREE(fit1d->srat    );
    FREE(fit1d->norm    );
    FREE(fit1d->f       );

    FREE(fit1d);

cleanup:
    FREE(nper);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit1d_objf - compute objective function                            *
 *                                                                      *
 ************************************************************************
 */
static int
fit1d_objf(fit1d_T *fit1d,              /* (in)  pointer to fit1d structure */
           double  smooth,              /* (in)  smoothing parameter */
           double  Tcloud[],            /* (in)  current T parameters */
           double  cp[],                /* (in)  current control points */
           double  srat[],              /* (in)  relative spacing ratio of inputs */
           double  f[])                 /* (out) objective functions */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    next, j, k;
    double XYZ[3];

    ROUTINE(fit1d_obj);

    /* --------------------------------------------------------------- */

    next = 0;

    for (k = 0; k < fit1d->m; k++) {
        status = eval1dBspline(Tcloud[k], fit1d->n, cp, XYZ, NULL, NULL);
        CHECK_STATUS(eval1dBspline);

        f[next++] = XYZ[0] - fit1d->XYZcloud[3*k  ];
        f[next++] = XYZ[1] - fit1d->XYZcloud[3*k+1];
        f[next++] = XYZ[2] - fit1d->XYZcloud[3*k+2];
    }

    for (j = 1; j < fit1d->n-1; j++) {
        f[next++] = smooth * (2 * cp[3*j  ] - (1+srat[j]) * cp[3*j-3] - (1-srat[j]) * cp[3*j+3]);
        f[next++] = smooth * (2 * cp[3*j+1] - (1+srat[j]) * cp[3*j-2] - (1-srat[j]) * cp[3*j+4]);
        f[next++] = smooth * (2 * cp[3*j+2] - (1+srat[j]) * cp[3*j-1] - (1-srat[j]) * cp[3*j+5]);
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit2dCloud - find spline that best-fits the cloud of points        *
 *                                                                      *
 ************************************************************************
 */
int
fit2dCloud(int    iface,                /* (in)  Face index */
           int    m,                    /* (in)  number of points in cloud */
           int    bitflag,              /* (in)  2=uPeriodic, 4=vPeriodic, 8=intGiven */
           double XYZcloud[],           /* (in)  array  of points in cloud (x,y,z,x,...) */
           int    nu,                   /* (in)  number of control points in U direction */
           int    nv,                   /* (in)  number of control points in V direction */
           double cp[],                 /* (in)  array  of control points (boundaries set) */
                                        /* (out) array  of control points (all set) */
           double smooth,               /* (in)  initial control net smoothing parameter */
           int    rtype,                /* (in)  restriction type (1=x, 2=y, 3=z) */
           int    nsmth,                /* (in)  number of post-smoothing passes */
           double UVcloud[],            /* (out) UV-parameters of points in cloud (u,v,u,...) */
           double *normf,               /* (out) RMS of distances between cloud and fit */
           double *maxf,                /* (out) maximum distance between cloud and fit */
           double *dotmin,              /* (in)  minimum allowable normalized dot product of control polygon */
                                        /* (out) minimum actual    normalized dot product of control polygon */
 /*@null@*/int    *nmin,                /* (out) minimum number of cloud points in any interval */
           int    *numiter,             /* (in)  if >=0, number of iterations allowed */
                                        /* (out) number of iterations executed */
 /*@null@*/FILE   *fp)                  /* (in)  file for progress outputs (or NULL) */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    i, j, iter, accept, niter=100, ismth;
    double old_normf, old_maxf, dotallow, toler=EPS06;
    double *cpsmth=NULL;

    fit2d_T *fit2d;
    void    *myFit2d;

#define IJ(I,J,XYZ) 3*((I)+(J)*fit2d->nu)+(XYZ)

    ROUTINE(fit2dCloud);

    /* --------------------------------------------------------------- */

    if (fp != NULL) {
        fprintf(fp, "enter fit2dCloud(iface=%d, m=%d, bitflag=%d, nu=%d, nv=%d, smooth=%f, rtype=%d, nsmth=%d)\n",
                iface, m, bitflag, nu, nv, smooth, rtype, nsmth);
    }

    if (m <= 1 || XYZcloud == NULL) {
        status = FIT_EMPTYCLOUD;
        goto cleanup;
    } else if (nu <= 2 || nv <= 2 || cp == NULL) {
        status = FIT_NCP;
        goto cleanup;
    } else if (smooth < 0 || smooth > 1) {
        status = FIT_SMOOTH;
        goto cleanup;
    }

    if (*numiter >= 0) {
        niter = *numiter;
    }

    /* default returns */
    *normf   = EPS12;
    *maxf    = 0;
    if (nmin != NULL) *nmin = m;
    *numiter = 0;

    /* minimum allowable normalized dot product */
    if (m >= 3) {
        dotallow = *dotmin;
    } else {
        dotallow = 0;
        *dotmin  = 1;
    }

    /* initialize */
    status = fit2d_init(m, bitflag, smooth, XYZcloud, nu, nv, cp, rtype, fp, normf, maxf, &myFit2d);
    CHECK_STATUS(fit2d_init);

    fit2d = (fit2d_T *) myFit2d;
    fit2d->iface = iface;
    fit2d->nsmth = nsmth;

    /* create a plotfile for the original control net and the normals */
    if (MAKE_NORMALS_PLOTFILE) {
        char   filename[80];
        FILE   *fp_plot2;

        snprintf(filename, 79, "normals_%d.plot", fit2d->iface);
        fp_plot2 = fopen(filename, "w");

        /* control net */
        fprintf(fp_plot2, " %9d %9d init_%d|b", fit2d->nu, fit2d->nv, fit2d->iface);

        for (j = 0; j < fit2d->nv; j++) {
            for (i = 0; i < fit2d->nu; i++) {
                fprintf(fp_plot2, "%15.7f %15.7f %15.7f\n",
                        fit2d->scale * fit2d->cp[IJ(i,j,0)] + fit2d->xavg,
                        fit2d->scale * fit2d->cp[IJ(i,j,1)] + fit2d->yavg,
                        fit2d->scale * fit2d->cp[IJ(i,j,2)] + fit2d->zavg);
            }
        }

        /* normals */
        fprintf(fp_plot2, " %9d %9d tuft_%d|r\n", (fit2d->nu)*(fit2d->nv), -1, fit2d->iface);

        for (j = 0; j < fit2d->nv; j++) {
            for (i = 0; i < fit2d->nu; i++) {
                fprintf(fp_plot2, "%15.7f %15.7f %15.7f  ",
                        fit2d->scale * fit2d->cp[IJ(i,j,0)] + fit2d->xavg,
                        fit2d->scale * fit2d->cp[IJ(i,j,1)] + fit2d->yavg,
                        fit2d->scale * fit2d->cp[IJ(i,j,2)] + fit2d->zavg);
                fprintf(fp_plot2, "%15.7f %15.7f %15.7f\n",
                        fit2d->scale * (fit2d->cp[IJ(i,j,0)] + 0.01 * fit2d->norm[IJ(i,j,0)]) + fit2d->xavg,
                        fit2d->scale * (fit2d->cp[IJ(i,j,1)] + 0.01 * fit2d->norm[IJ(i,j,1)]) + fit2d->yavg,
                        fit2d->scale * (fit2d->cp[IJ(i,j,2)] + 0.01 * fit2d->norm[IJ(i,j,2)]) + fit2d->zavg);
            }
        }

        fprintf(fp_plot2, "%9d %9d end\n", 0, 0);

        fclose(fp_plot2);
    }

    if (DEBUG) {
        int  ip0, ip1, ip2, ip3;
        FILE *fp_plot;
        fp_plot = fopen("slugs.plot", "a");
        if (fp_plot != NULL) {
            fprintf(fp_plot, " %9d %9d init_%d\n", (nu-1)*(nv-1), -4, iface);

            for (j = 1; j < nv; j++) {
                for (i = 1; i < nu; i++) {
                    ip0 = 3 * ((i-1) + nu * (j-1));
                    ip1 = 3 * ((i  ) + nu * (j-1));
                    ip2 = 3 * ((i  ) + nu * (j  ));
                    ip3 = 3 * ((i-1) + nu * (j  ));
                    fprintf(fp_plot, "%12.6f %12.6f %12.6f %4.0f  ", fit2d->scale * fit2d->cp[ip0  ] + fit2d->xavg,
                                                                     fit2d->scale * fit2d->cp[ip0+1] + fit2d->yavg,
                                                                     fit2d->scale * fit2d->cp[ip0+2] + fit2d->zavg, -1.);
                    fprintf(fp_plot, "%12.6f %12.6f %12.6f %4.0f  ", fit2d->scale * fit2d->cp[ip1  ] + fit2d->xavg,
                                                                     fit2d->scale * fit2d->cp[ip1+1] + fit2d->yavg,
                                                                     fit2d->scale * fit2d->cp[ip1+2] + fit2d->zavg,  0.);
                    fprintf(fp_plot, "%12.6f %12.6f %12.6f %4.0f  ", fit2d->scale * fit2d->cp[ip2  ] + fit2d->xavg,
                                                                     fit2d->scale * fit2d->cp[ip2+1] + fit2d->yavg,
                                                                     fit2d->scale * fit2d->cp[ip2+2] + fit2d->zavg,  1.);
                    fprintf(fp_plot, "%12.6f %12.6f %12.6f %4.0f\n", fit2d->scale * fit2d->cp[ip3  ] + fit2d->xavg,
                                                                     fit2d->scale * fit2d->cp[ip3+1] + fit2d->yavg,
                                                                     fit2d->scale * fit2d->cp[ip3+2] + fit2d->zavg,  0.);
                }
            }
            fclose(fp_plot);
        }
    }

    /* Levenberg-Marquardt iterations (if not satisfied with initial guess) */
    if (*normf > toler) {
        for (iter = 0; iter < niter; iter++) {
            (*numiter)++;

            old_normf = *normf;
            old_maxf  = *maxf;
            *dotmin   = dotallow;

            status = fit2d_step(myFit2d, smooth, normf, maxf, dotmin, &accept);
            CHECK_STATUS(fit2d_step);

            /* check for convergence */
            if (accept == 1 && fabs(*normf-old_normf) < toler && fabs(*maxf-old_maxf) < toler) {
                if (fp != NULL) {
                    fprintf(fp, "converged in %d itertions\n", iter);
                    fprintf(fp, "final     normf=%10.4e, maxf=%10.4e\n", *normf, *maxf);
                }
                break;
            } else if (accept == 0) {
                *normf = old_normf;
                *maxf  = old_maxf;
            }

            /* reduce smoothing for next step */
            smooth *= 0.99;
        }
    }

    /* if rtype>0, update the uv parametric coordinates by inverse
       evaluations (which may be slow) */
    if (fit2d->rtype > 0) {
        /*@-shadow@*/
        int    ipnt;
        int    header[7], ndata, ij, ncp=nu, periodic;
        double *cpdata=NULL, data[18], uvrange[4];
        double dist_old, dist_new, uvbest[2];
        ego    context, esurf, eface;
        int nkeep=0;
        double data2[18];

        /* new EGADS context */
        status = EG_open(&context);
        CHECK_STATUS(EG_open);

        /* create the Surface */
        header[0] = 0;            // bitflag
        header[1] = 3;            // u degree
        header[2] = ncp;          // u number of control points
        header[3] = ncp + 4;      // u number of knots
        header[4] = 3;            // v degree
        header[5] = ncp;          // v number of control points
        header[6] = ncp + 4;      // v number of knots

        ndata = 2 * (ncp+4) + 3 * ncp * ncp;
        MALLOC(cpdata, double, ndata);

        ndata = 0;

        /* knot vectors for u and v */
        for (i = 0; i < 2; i++) {
            cpdata[ndata++] = 0;
            cpdata[ndata++] = 0;
            cpdata[ndata++] = 0;
            cpdata[ndata++] = 0;

            for (j = 1; j < ncp-3; j++) {
                cpdata[ndata++] = j;
            }

            cpdata[ndata++] = ncp-3;
            cpdata[ndata++] = ncp-3;
            cpdata[ndata++] = ncp-3;
            cpdata[ndata++] = ncp-3;
        }

        /* control points */
        for (ij = 0; ij < ncp*ncp; ij++) {
            cpdata[ndata++] = fit2d->scale * fit2d->cp[3*ij  ] + fit2d->xavg;
            cpdata[ndata++] = fit2d->scale * fit2d->cp[3*ij+1] + fit2d->yavg;
            cpdata[ndata++] = fit2d->scale * fit2d->cp[3*ij+2] + fit2d->zavg;
        }

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                                 header, cpdata, &esurf);
        FREE(cpdata);
        CHECK_STATUS(EG_makeGeometry);

        (void) EG_getRange(esurf, uvrange, &periodic);

        status = EG_makeFace(esurf, SFORWARD, uvrange, &eface);
        CHECK_STATUS(EG_makeFace);

        /* reassign UVcloud by inverse evaluations */
        printf("start inverse evaluations, iface=%d\n", fit2d->iface);
        for (ipnt = 0; ipnt < m; ipnt++) {
            status = EG_evaluate(eface, &(fit2d->UVcloud[2*ipnt]), data2);
            CHECK_STATUS(EG_evaluate);

            dist_old = (XYZcloud[3*ipnt  ]-data2[0]) * (XYZcloud[3*ipnt  ]-data2[0])
                     + (XYZcloud[3*ipnt+1]-data2[1]) * (XYZcloud[3*ipnt+1]-data2[1])
                     + (XYZcloud[3*ipnt+2]-data2[2]) * (XYZcloud[3*ipnt+2]-data2[2]);

            status = EG_invEvaluate(eface, &(XYZcloud[3*ipnt]), uvbest, data);
            CHECK_STATUS(EG_invEvaluate);

            dist_new = (XYZcloud[3*ipnt  ]-data[0]) * (XYZcloud[3*ipnt  ]-data[0])
                     + (XYZcloud[3*ipnt+1]-data[1]) * (XYZcloud[3*ipnt+1]-data[1])
                     + (XYZcloud[3*ipnt+2]-data[2]) * (XYZcloud[3*ipnt+2]-data[2]);

            if (dist_new < dist_old) {
                fit2d->UVcloud[2*ipnt  ] = uvbest[0];
                fit2d->UVcloud[2*ipnt+1] = uvbest[1];
            } else {
                nkeep++;
            }
        }
        printf("end inverse evaluations, iface=%d, nkeep=%d\n", fit2d->iface, nkeep);
    }

    /* compute outputs, statistics, and clean up */
    status = fit2d_done(myFit2d, UVcloud, cp, normf, maxf, nmin);
    CHECK_STATUS(fit2d_done);

    /* smooth the control net */
    MALLOC(cpsmth, double, 3*nu*nv);

    for (ismth = 0; ismth < nsmth; ismth++) {
        printf("                 iface %3d smoothing pass %2d\n", iface, ismth);
        for (j = 1; j < nv-1; j++) {
            for (i = 1; i < nu-1; i++) {
                cpsmth[3*((i)+nu*(j))  ] = ( cp[3*((i-1)+nu*(j  ))  ]
                                           + cp[3*((i+1)+nu*(j  ))  ]
                                           + cp[3*((i  )+nu*(j-1))  ]
                                           + cp[3*((i  )+nu*(j+1))  ]) / 4;
                cpsmth[3*((i)+nu*(j))+1] = ( cp[3*((i-1)+nu*(j  ))+1]
                                           + cp[3*((i+1)+nu*(j  ))+1]
                                           + cp[3*((i  )+nu*(j-1))+1]
                                           + cp[3*((i  )+nu*(j+1))+1]) / 4;
                cpsmth[3*((i)+nu*(j))+2] = ( cp[3*((i-1)+nu*(j  ))+2]
                                           + cp[3*((i+1)+nu*(j  ))+2]
                                           + cp[3*((i  )+nu*(j-1))+2]
                                           + cp[3*((i  )+nu*(j+1))+2]) / 4;
            }
        }
        for (j = 1; j < nv-1; j++) {
            for (i = 1; i < nu-1; i++) {
                cp[3*((i)+nu*(j))  ] = ( cp[    3*((i)+nu*(j))  ]
                                       + cpsmth[3*((i)+nu*(j))  ]) / 2;
                cp[3*((i)+nu*(j))+1] = ( cp[    3*((i)+nu*(j))+1]
                                       + cpsmth[3*((i)+nu*(j))+1]) / 2;
                cp[3*((i)+nu*(j))+2] = ( cp[    3*((i)+nu*(j))+2]
                                       + cpsmth[3*((i)+nu*(j))+2]) / 2;
            }
        }
    }
    FREE(cpsmth);

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit2d_init - initialize the B-spline surface fitter                *
 *                                                                      *
 ************************************************************************
 */
int
fit2d_init(int     m,                   /* (in)  number of points in cloud */
           int     bitflag,             /* (in)  2=uPeriodic, 4=vPeriodic, 8=intGiven */
           double  smooth,              /* (in)  control net smoothing parameter */
           double  XYZcloud[],          /* (in)  array  of points in cloud (x,y,z,x,...) */
           int     nu,                  /* (in)  number of control points in u direction */
           int     nv,                  /* (in)  number of control points in v direction */
           double  cp[],                /* (in)  array  of control points (with outline set) */
           int     rtype,               /* (in)  restriction type (1=x, 2=y, 3=z) */
 /*@null@*/FILE    *fp,                 /* (in)  file pointer for progress prints */
           double  *normf,              /* (out) initial RMS between cloud and fit */
           double  *maxf,               /* (out) maximum distance between cloud and fit */
           void    **context)           /* (out) pointer to context */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    uPeriodic=0, vPeriodic=0, intGiven=0;
    int    nobj, nmask, i, j, k, ivar, jvar, nchange;
    double xavg, yavg, zavg, dx, dy, dz, dotprod;
    double fraci, fracj, dbest, dtest;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double vec1x, vec1y, vec1z, vec2x, vec2y, vec2z, vec3x, vec3y, vec3z, vec3;
    double *temp=NULL, *normsmth=NULL;

    fit2d_T *fit2d=NULL;

#define MASK(I,J) fit2d->mask[(I)+(J)*nmask]

    ROUTINE(fit2d_init);

    /* --------------------------------------------------------------- */

    if ((bitflag & 2) != 0) uPeriodic = 1;
    if ((bitflag & 4) != 0) vPeriodic = 1;
#ifndef __clang_analyzer__
    if ((bitflag & 8) != 0) intGiven  = 1;
#endif

    /* periodic is not implemented (yet) */
    if (uPeriodic == 1 || vPeriodic == 1) {
        status = FIT_BITFLAG;
        goto cleanup;
    }

    /* cannot be both U and V periodic */
    if (uPeriodic == 1 && vPeriodic == 1) {
        status = FIT_BITFLAG;
        goto cleanup;
    }

    /* get a fit2d structure */
    MALLOC(fit2d, fit2d_T, 1);

    /* initialize the fit2d structure */
    fit2d->m        = m;
    fit2d->XYZcloud = NULL;
    fit2d->UVcloud  = NULL;

    fit2d->nu       = nu;
    fit2d->nv       = nv;
    fit2d->cp       = NULL;
    fit2d->norm     = NULL;
    fit2d->f        = NULL;

    MALLOC(fit2d->XYZcloud, double, 3*fit2d->m);
    MALLOC(fit2d->UVcloud,  double, 2*fit2d->m);

    MALLOC(fit2d->cp,    double,            3*fit2d->nu*fit2d->nv);
    MALLOC(fit2d->norm,  double,            3*fit2d->nu*fit2d->nv);
    MALLOC(fit2d->f,     double, 3*fit2d->m+3*fit2d->nu*fit2d->nv);    /* a little too big */

    fit2d->bitflag  = bitflag;
    fit2d->rtype    = rtype;
    fit2d->nsmth    = 0;
    fit2d->iter     = 0;
    fit2d->lambda   = 1;

    fit2d->scale    = 1;
    fit2d->xavg     = 0;
    fit2d->yavg     = 0;
    fit2d->zavg     = 0;

    fit2d->mask     = NULL;

    fit2d->fp       = fp;

    /* find extrema of input data */
    xmin = XYZcloud[0];
    xmax = XYZcloud[0];
    ymin = XYZcloud[1];
    ymax = XYZcloud[1];
    zmin = XYZcloud[2];
    zmax = XYZcloud[2];

    for (k = 0; k < fit2d->m; k++) {
        if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
        if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
        if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
        if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
        if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
        if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
    }

    /* normalize input data */
                                  fit2d->scale = xmax - xmin;
    if (ymax-ymin > fit2d->scale) fit2d->scale = ymax - ymin;
    if (zmax-zmin > fit2d->scale) fit2d->scale = zmax - zmin;

    if (fit2d->scale < EPS12) {
        printf("%20.12e %20.12e %20.12e\n", xmin, ymin, zmin);
        printf("%20.12e %20.12e %20.12e\n", xmax, ymax, zmax);
        status = FIT_DEGENERATE;
        goto cleanup;
    }

    fit2d->xavg = (xmax + xmin) / 2;
    fit2d->yavg = (ymax + ymin) / 2;
    fit2d->zavg = (zmax + zmin) / 2;

    for (k = 0; k < fit2d->m; k++) {
        fit2d->XYZcloud[3*k  ] = (XYZcloud[3*k  ] - fit2d->xavg) / fit2d->scale;
        fit2d->XYZcloud[3*k+1] = (XYZcloud[3*k+1] - fit2d->yavg) / fit2d->scale;
        fit2d->XYZcloud[3*k+2] = (XYZcloud[3*k+2] - fit2d->zavg) / fit2d->scale;
    }

    for (j = 0;  j < fit2d->nv; j++) {
        for (i = 0; i < fit2d->nu; i++) {
            fit2d->cp[IJ(i,j,0)] = (cp[IJ(i,j,0)] - fit2d->xavg) / fit2d->scale;
            fit2d->cp[IJ(i,j,1)] = (cp[IJ(i,j,1)] - fit2d->yavg) / fit2d->scale;
            fit2d->cp[IJ(i,j,2)] = (cp[IJ(i,j,2)] - fit2d->zavg) / fit2d->scale;
        }
    }

    /* number of design variables and objectives */
    nobj  = 3 * fit2d->m + 3 * (fit2d->nu - 2) * (fit2d->nv - 2);
    nmask =                    (fit2d->nu - 2) * (fit2d->nv - 2);

    assert (nmask > 0);

    /* bilinear interpolate control points from boundaries */
#ifndef __clang_analyzer__
    if (intGiven == 0) {
        for (j = 1; j < fit2d->nv-1; j++) {
            for (i =  1; i < fit2d->nu-1; i++) {
                fraci = (double)(i) / (double)(fit2d->nu-1);
                fracj = (double)(j) / (double)(fit2d->nv-1);

                fit2d->cp[IJ(i,j,0)]
                    = (1-fraci)             * fit2d->cp[IJ(0,          j,          0)]
                    + (  fraci)             * fit2d->cp[IJ(fit2d->nu-1,j,          0)]
                    +             (1-fracj) * fit2d->cp[IJ(i,          0,          0)]
                    +             (  fracj) * fit2d->cp[IJ(i,          fit2d->nv-1,0)]
                    - (1-fraci) * (1-fracj) * fit2d->cp[IJ(0,          0,          0)]
                    - (  fraci) * (1-fracj) * fit2d->cp[IJ(fit2d->nu-1,0,          0)]
                    - (1-fraci) * (  fracj) * fit2d->cp[IJ(0,          fit2d->nv-1,0)]
                    - (  fraci) * (  fracj) * fit2d->cp[IJ(fit2d->nu-1,fit2d->nv-1,0)];
                fit2d->cp[IJ(i,j,1)]
                    = (1-fraci)             * fit2d->cp[IJ(0,          j,          1)]
                    + (  fraci)             * fit2d->cp[IJ(fit2d->nu-1,j,          1)]
                    +             (1-fracj) * fit2d->cp[IJ(i,          0,          1)]
                    +             (  fracj) * fit2d->cp[IJ(i,          fit2d->nv-1,1)]
                    - (1-fraci) * (1-fracj) * fit2d->cp[IJ(0,          0,          1)]
                    - (  fraci) * (1-fracj) * fit2d->cp[IJ(fit2d->nu-1,0,          1)]
                    - (1-fraci) * (  fracj) * fit2d->cp[IJ(0,          fit2d->nv-1,1)]
                    - (  fraci) * (  fracj) * fit2d->cp[IJ(fit2d->nu-1,fit2d->nv-1,1)];
                fit2d->cp[IJ(i,j,2)]
                    = (1-fraci)             * fit2d->cp[IJ(0,          j,          2)]
                    + (  fraci)             * fit2d->cp[IJ(fit2d->nu-1,j,          2)]
                    +             (1-fracj) * fit2d->cp[IJ(i,          0,          2)]
                    +             (  fracj) * fit2d->cp[IJ(i,          fit2d->nv-1,2)]
                    - (1-fraci) * (1-fracj) * fit2d->cp[IJ(0,          0,          2)]
                    - (  fraci) * (1-fracj) * fit2d->cp[IJ(fit2d->nu-1,0,          2)]
                    - (1-fraci) * (  fracj) * fit2d->cp[IJ(0,          fit2d->nv-1,2)]
                    - (  fraci) * (  fracj) * fit2d->cp[IJ(fit2d->nu-1,fit2d->nv-1,2)];
            }
        }
    }
#endif

    /* debug prints */
    if (DEBUG) {
        for (j = 0; j < fit2d->nv; j++) {
            printf("X %3d: ", j);
            for (i = 0; i < fit2d->nu; i++) {
                printf("%10.5f ", fit2d->cp[IJ(i,j,0)]);
            }
            printf("\n");
        }
        printf("\n");
        for (j = 0; j < fit2d->nv; j++) {
            printf("Y %3d: ", j);
            for (i = 0; i < fit2d->nu; i++) {
                printf("%10.5f ", fit2d->cp[IJ(i,j,1)]);
            }
            printf("\n");
        }
        printf("\n");
        for (j = 0; j < fit2d->nv; j++) {
            printf("Z %3d: ", j);
            for (i = 0; i < fit2d->nu; i++) {
                printf("%10.5f ", fit2d->cp[IJ(i,j,2)]);
            }
        printf("\n");
        }
        printf("\n");
    }

    /* compute the normals to the original control net */
    nchange = 1;
    while (nchange > 0) {
        nchange = 0;

#ifndef __clang_analyzer__
        for (j = 0; j < fit2d->nv; j++) {
            for (i = 0; i < fit2d->nu; i++) {
                if        (i == 0) {
                    vec1x  = fit2d->cp[IJ(i+1,j,0)] - fit2d->cp[IJ(i  ,j,0)];
                    vec1y  = fit2d->cp[IJ(i+1,j,1)] - fit2d->cp[IJ(i  ,j,1)];
                    vec1z  = fit2d->cp[IJ(i+1,j,2)] - fit2d->cp[IJ(i  ,j,2)];
                } else if (i == fit2d->nu-1) {
                    vec1x  = fit2d->cp[IJ(i  ,j,0)] - fit2d->cp[IJ(i-1,j,0)];
                    vec1y  = fit2d->cp[IJ(i  ,j,1)] - fit2d->cp[IJ(i-1,j,1)];
                    vec1z  = fit2d->cp[IJ(i  ,j,2)] - fit2d->cp[IJ(i-1,j,2)];
                } else {
                    vec1x  = fit2d->cp[IJ(i+1,j,0)] - fit2d->cp[IJ(i-1,j,0)];
                    vec1y  = fit2d->cp[IJ(i+1,j,1)] - fit2d->cp[IJ(i-1,j,1)];
                    vec1z  = fit2d->cp[IJ(i+1,j,2)] - fit2d->cp[IJ(i-1,j,2)];
                }

                if        (j == 0) {
                    vec2x  = fit2d->cp[IJ(i,j+1,0)] - fit2d->cp[IJ(i,j  ,0)];
                    vec2y  = fit2d->cp[IJ(i,j+1,1)] - fit2d->cp[IJ(i,j  ,1)];
                    vec2z  = fit2d->cp[IJ(i,j+1,2)] - fit2d->cp[IJ(i,j  ,2)];
                } else if (j == fit2d->nv-1) {
                    vec2x  = fit2d->cp[IJ(i,j  ,0)] - fit2d->cp[IJ(i,j-1,0)];
                    vec2y  = fit2d->cp[IJ(i,j  ,1)] - fit2d->cp[IJ(i,j-1,1)];
                    vec2z  = fit2d->cp[IJ(i,j  ,2)] - fit2d->cp[IJ(i,j-1,2)];
                } else {
                    vec2x  = fit2d->cp[IJ(i,j+1,0)] - fit2d->cp[IJ(i,j-1,0)];
                    vec2y  = fit2d->cp[IJ(i,j+1,1)] - fit2d->cp[IJ(i,j-1,1)];
                    vec2z  = fit2d->cp[IJ(i,j+1,2)] - fit2d->cp[IJ(i,j-1,2)];
                }

                vec3x = vec1y * vec2z - vec1z * vec2y;
                vec3y = vec1z * vec2x - vec1x * vec2z;
                vec3z = vec1x * vec2y - vec1y * vec2x;
                vec3  = sqrt(vec3x*vec3x + vec3y*vec3y + vec3z*vec3z);

                fit2d->norm[IJ(i,j,0)] = vec3x / vec3;
                fit2d->norm[IJ(i,j,1)] = vec3y / vec3;
                fit2d->norm[IJ(i,j,2)] = vec3z / vec3;
            }
        }
#endif

        /* smooth the normals */
#ifndef __clang_analyzer__
        MALLOC(normsmth, double, 3*fit2d->nu*fit2d->nv);

        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                normsmth[IJ(i,j,0)] = ( fit2d->norm[IJ(i-1,j  ,0)]
                                      + fit2d->norm[IJ(i+1,j  ,0)]
                                      + fit2d->norm[IJ(i  ,j-1,0)]
                                      + fit2d->norm[IJ(i  ,j+1,0)]) / 4;
                normsmth[IJ(i,j,1)] = ( fit2d->norm[IJ(i-1,j  ,1)]
                                      + fit2d->norm[IJ(i+1,j  ,1)]
                                      + fit2d->norm[IJ(i  ,j-1,1)]
                                      + fit2d->norm[IJ(i  ,j+1,1)]) / 4;
                normsmth[IJ(i,j,2)] = ( fit2d->norm[IJ(i-1,j  ,2)]
                                      + fit2d->norm[IJ(i+1,j  ,2)]
                                      + fit2d->norm[IJ(i  ,j-1,2)]
                                      + fit2d->norm[IJ(i  ,j+1,2)]) / 4;
            }
        }

        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                fit2d->norm[IJ(i,j,0)] = (fit2d->norm[IJ(i,j,0)] + normsmth[IJ(i,j,0)]) / 2;
                fit2d->norm[IJ(i,j,1)] = (fit2d->norm[IJ(i,j,1)] + normsmth[IJ(i,j,1)]) / 2;
                fit2d->norm[IJ(i,j,2)] = (fit2d->norm[IJ(i,j,2)] + normsmth[IJ(i,j,2)]) / 2;
            }
        }

        FREE(normsmth);
#endif

        /* bump out any control points to relieve any local concavities */
#ifndef __clang_analyzer__
        if (fit2d->rtype > 0) {
            for (j = 1; j < fit2d->nv-1; j++) {
                for (i = 1; i < fit2d->nu-1; i++) {
                    xavg = (fit2d->cp[IJ(i-1,j,0)] + fit2d->cp[IJ(i+1,j,0)]) / 2;
                    yavg = (fit2d->cp[IJ(i-1,j,1)] + fit2d->cp[IJ(i+1,j,1)]) / 2;
                    zavg = (fit2d->cp[IJ(i-1,j,2)] + fit2d->cp[IJ(i+1,j,2)]) / 2;

                    dx = fit2d->cp[IJ(i,j,0)] - xavg;
                    dy = fit2d->cp[IJ(i,j,1)] - yavg;
                    dz = fit2d->cp[IJ(i,j,2)] - zavg;

                    dotprod = dx * fit2d->norm[IJ(i,j,0)]
                            + dy * fit2d->norm[IJ(i,j,1)]
                            + dz * fit2d->norm[IJ(i,j,2)];

                    if (dotprod < -EPS06) {
                        fit2d->cp[IJ(i,j,0)] = xavg;
                        fit2d->cp[IJ(i,j,1)] = yavg;
                        fit2d->cp[IJ(i,j,2)] = zavg;

                        nchange++;
                    }

                    xavg = (fit2d->cp[IJ(i,j-1,0)] + fit2d->cp[IJ(i,j+1,0)]) / 2;
                    yavg = (fit2d->cp[IJ(i,j-1,1)] + fit2d->cp[IJ(i,j+1,1)]) / 2;
                    zavg = (fit2d->cp[IJ(i,j-1,2)] + fit2d->cp[IJ(i,j+1,2)]) / 2;

                    dx = fit2d->cp[IJ(i,j,0)] - xavg;
                    dy = fit2d->cp[IJ(i,j,1)] - yavg;
                    dz = fit2d->cp[IJ(i,j,2)] - zavg;

                    dotprod = dx * fit2d->norm[IJ(i,j,0)]
                            + dy * fit2d->norm[IJ(i,j,1)]
                            + dz * fit2d->norm[IJ(i,j,2)];

                    if (dotprod < -EPS06) {
                        fit2d->cp[IJ(i,j,0)] = xavg;
                        fit2d->cp[IJ(i,j,1)] = yavg;
                        fit2d->cp[IJ(i,j,2)] = zavg;

                        nchange++;
                    }
                }
            }
        }
#endif
    }

    /* for each point in the cloud, assign the values of "u" and "v"
       that are associated with the closest interior control point.
       note: only interior control points are used since corner points
       can cause problems */
#ifndef __clang_analyzer__
    for (k = 0; k < fit2d->m; k++) {
        fit2d->UVcloud[2*k  ] = 0;
        fit2d->UVcloud[2*k+1] = 0;
        dbest       = 1e20;

        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                if (intGiven == 1) {
                    dtest = SQR(       XYZcloud[3*k  ] -        cp[IJ(i,j,0)])
                          + SQR(       XYZcloud[3*k+1] -        cp[IJ(i,j,1)])
                          + SQR(       XYZcloud[3*k+2] -        cp[IJ(i,j,2)]);
                } else {
                    dtest = SQR(fit2d->XYZcloud[3*k  ] - fit2d->cp[IJ(i,j,0)])
                          + SQR(fit2d->XYZcloud[3*k+1] - fit2d->cp[IJ(i,j,1)])
                          + SQR(fit2d->XYZcloud[3*k+2] - fit2d->cp[IJ(i,j,2)]);
                }
                if (dtest < dbest) {
                    fit2d->UVcloud[2*k  ] = (double)(i) / (double)(fit2d->nu-1) * (double)(fit2d->nu-3);
                    fit2d->UVcloud[2*k+1] = (double)(j) / (double)(fit2d->nv-1) * (double)(fit2d->nv-3);
                    dbest       = dtest;
                }
            }
        }
    }
#endif

    /* compute the initial objective function */
    status = fit2d_objf(fit2d, smooth, fit2d->UVcloud, fit2d->cp, fit2d->f);
    CHECK_STATUS(fit2d_objf);

    /* compute and report the norm of the objective function */
    *maxf  = Linorm(fit2d->f, 3*fit2d->m);
    *normf = L2norm(fit2d->f, nobj) / sqrt(nobj);
    if (fit2d->fp != NULL) {
        fprintf(fit2d->fp, "initial   normf=%10.4e, maxf=%10.4e\n", *normf, *maxf);
    }

    /* create MASK array, which contains the smoothing vectors */
    MALLOC(fit2d->mask, int, nmask*nmask);

    i = 0;
    j = 0;
    for (ivar = 0; ivar < nmask; ivar++) {
        for (jvar = 0; jvar < nmask; jvar++) {
            MASK(ivar,jvar) = 0;
        }

        if (              1)     MASK(ivar,ivar            ) = +4;
        if (i >           0)     MASK(ivar,ivar-1          ) = -2;
        if (i < fit2d->nu-3)     MASK(ivar,ivar+1          ) = -2;

        if (j > 0) {
            if (              1) MASK(ivar,ivar-fit2d->nu+2) = -2;
            if (i >           0) MASK(ivar,ivar-fit2d->nu+1) = +1;
            if (i < fit2d->nu-3) MASK(ivar,ivar-fit2d->nu+3) = +1;
        }

        if (j < fit2d->nv-3) {
            if (              1) MASK(ivar,ivar+fit2d->nu-2) = -2;
            if (i >           0) MASK(ivar,ivar+fit2d->nu-3) = +1;
            if (i < fit2d->nu-3) MASK(ivar,ivar+fit2d->nu-1) = +1;
        }

        i++;
        if (i >= fit2d->nu-2) {
            i = 0;
            j++;
        }
    }

    /* debug prints */
    if (DEBUG) {
        for (ivar = 0; ivar < nmask; ivar++) {
            for (jvar = 0; jvar < nmask; jvar++) {
                if (MASK(ivar,jvar) == 0) {
                    printf("  . ");
                } else {
                    printf("%3d ", MASK(ivar,jvar));
                }
                if (jvar%(fit2d->nv-2) == fit2d->nv-3) printf("| ");
            }
            printf("\n");
            if (ivar%(fit2d->nu-2) == fit2d->nu-3) printf("\n");
        }
    }

    /* return the fit2d structure (as the context) */
    *context = (void*)fit2d;

cleanup:
    if (status != FIT_SUCCESS && fit2d != NULL) {
        FREE(fit2d->XYZcloud);
        FREE(fit2d->UVcloud );
        FREE(fit2d->cp      );
        FREE(fit2d->norm    );
        FREE(fit2d->f       );

        FREE(fit2d);
    }

    FREE(temp);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit2d_step - take one step of the B-spline surface fitter          *
 *                                                                      *
 ************************************************************************
 */
int
fit2d_step(void    *context,            /* (in)  context */
           double  smooth,              /* (in)  control net smoothing parameter */
           double  *normf,              /* (out) RMS of distances between cloud and fit */
           double  *maxf,               /* (out) maximum distance between cloud and fit */
           double  *dotmin,             /* (in)  minimum allowable normalized dot product of control polygon */
                                        /* (out) minimum actual    normalized dot product of control polygon */
           int     *accept)             /* (out) =1 if step was accepted */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    nvar, ivar, jvar, nobj, iobj, nmask, i, j, k, next;
    int    ii, jj, ij;
    double normfnew, maxfnew, sum, fact, sum0, sum1;
    double XYZ[3], dXYZdU[3], dXYZdV[3], *dXYZdP=NULL;
    double *cpnew=NULL, *beta=NULL, *betanew=NULL, *delta=NULL;
    double *fnew=NULL;

    fit2d_T *fit2d = (fit2d_T *) context;

    int   nn = 3 * (fit2d->nu - 2) * (fit2d->nv - 2);
    double *AA=NULL, *BB=NULL, *CC=NULL, *DD=NULL, *EE=NULL;

#define A(K,I)    AA[2*(K)+(I)]
#define B(K,IJ)   BB[(K)*nn+(IJ)]
#define C(I,J)    CC[(I)*nn+(J)]
#define D(K)      DD[(K)]
#define E(K)      EE[(K)]

    ROUTINE(fit2d_step);

    /* --------------------------------------------------------------- */

    *accept = 0;

    /* note: first 5 iterations only adjust U and V */
    (fit2d->iter)++;

    /* number of design variables and objectives */
    nvar  = 2 * fit2d->m + 3 * (fit2d->nu - 2) * (fit2d->nv - 2);
    nobj  = 3 * fit2d->m + 3 * (fit2d->nu - 2) * (fit2d->nv - 2);
#ifndef __clang_analyzer__
    nmask =                    (fit2d->nu - 2) * (fit2d->nv - 2);
#endif

    /* allocate all temporary arrays */
    MALLOC(dXYZdP,  double,   fit2d->nu*fit2d->nv);
    MALLOC(cpnew,   double, 3*fit2d->nu*fit2d->nv);

    MALLOC(beta,    double, nvar);
    MALLOC(delta,   double, nvar);
    MALLOC(betanew, double, nvar);

    MALLOC(fnew,    double, nobj);

    /* store the UV values in the first 2*m betas */
    next = 0;
    for (k = 0; k < fit2d->m; k++) {
        beta[next++] = fit2d->UVcloud[2*k  ];
        beta[next++] = fit2d->UVcloud[2*k+1];
    }

    /* store the interior control points into beta */
    for (j = 1; j < fit2d->nv-1; j++) {
        for (i = 1; i < fit2d->nu-1; i++) {
            beta[next++] = fit2d->cp[IJ(i,j,0)];
            beta[next++] = fit2d->cp[IJ(i,j,1)];
            beta[next++] = fit2d->cp[IJ(i,j,2)];
        }
    }
    assert(next == nvar);

    /* minimum allowable dot product */
    *dotmin  = 1;

    /*
          trans(J) * J =  [   A      B ]      trans(J) * Q =  [ D ]
                          [            ]                      [   ]
                          [trans(B)  C ]                      [ E ]
    */
    MALLOC(AA, double, (2*fit2d->m)*(2 ));
    MALLOC(BB, double, (2*fit2d->m)*(nn));
    MALLOC(CC, double, (  nn      )*(nn));
    MALLOC(DD, double, (2*fit2d->m)     );
    MALLOC(EE, double,              (nn));

    for (i = 0; i < nn; i++) {
        for (j = 0; j < nn; j++) {
            C(i,j) = 0;
        }
        E(i) = 0;
    }

    /* create top-left (A) and top-right (B) parts of JtJ and top part (D) of JtQ */
    for (k = 0; k < fit2d->m; k++) {
        status = eval2dBspline(beta[2*k], beta[2*k+1], fit2d->nu, fit2d->nv, fit2d->cp,
                               XYZ, dXYZdU, dXYZdV, dXYZdP);
        CHECK_STATUS(eval2dBspline);

        /* top-left (A) part of JtJ */
        A(2*k  ,0) = dXYZdU[0] * dXYZdU[0] + dXYZdU[1] * dXYZdU[1] + dXYZdU[2] * dXYZdU[2];
        A(2*k  ,1) = dXYZdU[0] * dXYZdV[0] + dXYZdU[1] * dXYZdV[1] + dXYZdU[2] * dXYZdV[2];
        A(2*k+1,0) = dXYZdV[0] * dXYZdV[0] + dXYZdV[1] * dXYZdV[1] + dXYZdV[2] * dXYZdV[2];
        A(2*k+1,1) = dXYZdV[0] * dXYZdU[0] + dXYZdV[1] * dXYZdU[1] + dXYZdV[2] * dXYZdU[2];

        /* top-right (B) part of JtJ */
        jvar = 0;
        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                ij = i + j * fit2d->nu;

                B(2*k  ,jvar  ) = dXYZdU[0] * dXYZdP[ij];
                B(2*k  ,jvar+1) = dXYZdU[1] * dXYZdP[ij];
                B(2*k  ,jvar+2) = dXYZdU[2] * dXYZdP[ij];

                B(2*k+1,jvar  ) = dXYZdV[0] * dXYZdP[ij];
                B(2*k+1,jvar+1) = dXYZdV[1] * dXYZdP[ij];
                B(2*k+1,jvar+2) = dXYZdV[2] * dXYZdP[ij];

                jvar += 3;
            }
        }

        /* top part (D) of JtQ (negative needed since f = (XYZ_spline - XYZ_cloud) */
        D(2*k  ) = - dXYZdU[0] * fit2d->f[3*k] - dXYZdU[1] * fit2d->f[3*k+1] - dXYZdU[2] * fit2d->f[3*k+2];
        D(2*k+1) = - dXYZdV[0] * fit2d->f[3*k] - dXYZdV[1] * fit2d->f[3*k+1] - dXYZdV[2] * fit2d->f[3*k+2];

#ifndef __clang_analyzer__
        jvar = 0;
        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                ij = i + j * fit2d->nu;

                E(jvar  ) -= dXYZdP[ij] * fit2d->f[3*k  ];
                E(jvar+1) -= dXYZdP[ij] * fit2d->f[3*k+1];
                E(jvar+2) -= dXYZdP[ij] * fit2d->f[3*k+2];

                jvar += 3;
            }
        }
#endif
    }

    /* create bottom-right (C) part of JtJ */
    for (k = 0; k < fit2d->m; k++) {
        status = eval2dBspline(beta[2*k], beta[2*k+1], fit2d->nu, fit2d->nv,
                               fit2d->cp, XYZ, dXYZdU, dXYZdV, dXYZdP);
        CHECK_STATUS(eval2dBspline);

#ifndef __clang_analyzer__
        ivar = 0;
        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                jvar = 0;
                for (jj = 1; jj < fit2d->nv-1; jj++) {
                    for (ii = 1; ii < fit2d->nu-1; ii++) {
                        C(ivar  ,jvar  ) += dXYZdP[i+j*fit2d->nu] * dXYZdP[ii+jj*fit2d->nu];
                        C(ivar+1,jvar+1) += dXYZdP[i+j*fit2d->nu] * dXYZdP[ii+jj*fit2d->nu];
                        C(ivar+2,jvar+2) += dXYZdP[i+j*fit2d->nu] * dXYZdP[ii+jj*fit2d->nu];

                        jvar += 3;
                    }
                }
                ivar += 3;
            }
        }
#endif
    }

#ifndef __clang_analyzer__
    for (jvar = 2*fit2d->m; jvar < nvar; jvar+=3) {
        for (ivar = 2*fit2d->m; ivar < nvar; ivar+=3) {

            i = (ivar - 2 * fit2d->m) / 3;
            j = (jvar - 2 * fit2d->m) / 3;

            sum = 0;
            for (k = 0; k < nmask; k++) {
                sum += MASK(i,k) * MASK(j,k);
            }

            C(ivar-2*fit2d->m  ,jvar-2*fit2d->m  ) += smooth * smooth * sum;
            C(ivar-2*fit2d->m+1,jvar-2*fit2d->m+1) += smooth * smooth * sum;
            C(ivar-2*fit2d->m+2,jvar-2*fit2d->m+2) += smooth * smooth * sum;
        }
    }

    /* create bottom (E) part of JtQ */
    for (iobj = 3*fit2d->m; iobj < nobj; iobj+=3) {
        for (ivar = 0; ivar < nn; ivar+=3) {

            i = (ivar               ) / 3;
            j = (iobj - 3 * fit2d->m) / 3;

            sum = 0;
            for (k = 0; k < nmask; k++) {
                sum += MASK(i,k) * MASK(j,k);
            }

            E(ivar  ) -= smooth * MASK(i,j) * fit2d->f[iobj  ];
            E(ivar+1) -= smooth * MASK(i,j) * fit2d->f[iobj+1];
            E(ivar+2) -= smooth * MASK(i,j) * fit2d->f[iobj+2];
        }
    }
#endif

    /* multiply diagonals of JtJ by (1 + lambda) */
    for (k = 0; k < fit2d->m; k++) {
        A(2*k  ,0) *= (1 + fit2d->lambda);
        A(2*k+1,0) *= (1 + fit2d->lambda);
    }

    for (i = 0; i < nn; i++) {
        C(i,i) *= (1 + fit2d->lambda);
    }

    /* update:  C = C - trans(B) * inv(A) * B
       update:  E = E - trans(B) * inv(A) * D  */
#ifndef __clang_analyzer__
    for (k = 0; k < fit2d->m; k++) {
        fact = 1 / (A(2*k,1) * A(2*k,1) - A(2*k,0) * A(2*k+1,0));

        for (i = 0; i < nn; i++) {
            for (j = 0; j < nn; j++) {
                C(i,j) -= (B(2*k  ,j) * (A(2*k,1) * B(2*k+1,i) - A(2*k+1,0) * B(2*k  ,i))
                         + B(2*k+1,j) * (A(2*k,1) * B(2*k  ,i) - A(2*k  ,0) * B(2*k+1,i))) * fact;
            }

            E(i)       -= (D(2*k  )   * (A(2*k,1) * B(2*k+1,i) - A(2*k+1,0) * B(2*k  ,i))
                         + D(2*k+1)   * (A(2*k,1) * B(2*k  ,i) - A(2*k  ,0) * B(2*k+1,i))) * fact;
        }
    }
#endif

    /* solve for the second part of beta (the control points) */
    status = matsol(CC, EE, nn, &(delta[2*fit2d->m]));
    CHECK_STATUS(matsol);

    /* solve for the first part of beta (the parametric coordinates) */
    for (k = 0; k < fit2d->m; k++) {
        fact = 1 / (A(2*k,1) * A(2*k,1) - A(2*k,0) * A(2*k+1,0));

        sum0 = - D(2*k  );
        sum1 = - D(2*k+1);

#ifndef __clang_analyzer__
        for (i = 0; i < nn; i++) {
            sum0 += B(2*k  ,i) * delta[2*fit2d->m+i];
            sum1 += B(2*k+1,i) * delta[2*fit2d->m+i];
        }
#endif

        delta[2*k  ] = (A(2*k+1,0) * sum0 - A(2*k,1) * sum1) * fact;
        delta[2*k+1] = (A(2*k  ,0) * sum1 - A(2*k,1) * sum0) * fact;
    }

    /* find the temporary new beta (and clip the UVclouds) */
#ifndef __clang_analyzer__
    for (ivar = 0; ivar < nvar; ivar++) {
        betanew[ivar] = beta[ivar] + delta[ivar];
        if (ivar < 2*fit2d->m) {
            if (ivar%2 == 0 && betanew[ivar] < 0          ) betanew[ivar] = 0;
            if (ivar%2 == 0 && betanew[ivar] > fit2d->nu-3) betanew[ivar] = fit2d->nu-3;
            if (ivar%2 == 1 && betanew[ivar] < 0          ) betanew[ivar] = 0;
            if (ivar%2 == 1 && betanew[ivar] > fit2d->nv-3) betanew[ivar] = fit2d->nv-3;
        }
    }
#endif

    /* extract the temporary control points from betanew */
    next = 2 * fit2d->m;
    for (j = 0; j < fit2d->nv; j++) {
        for (i = 0; i < fit2d->nu; i++) {
            if (i == 0 || i == fit2d->nu-1 || j == 0 || j == fit2d->nv-1) {
                cpnew[IJ(i,j,0)] = fit2d->cp[IJ(i,j,0)];
                cpnew[IJ(i,j,1)] = fit2d->cp[IJ(i,j,1)];
                cpnew[IJ(i,j,2)] = fit2d->cp[IJ(i,j,2)];
            } else if (fit2d->iter <= 5) {
                cpnew[IJ(i,j,0)] = fit2d->cp[IJ(i,j,0)];
                cpnew[IJ(i,j,1)] = fit2d->cp[IJ(i,j,1)];
                cpnew[IJ(i,j,2)] = fit2d->cp[IJ(i,j,2)];
                next += 3;
            } else {
                cpnew[IJ(i,j,0)] = betanew[next++];
                cpnew[IJ(i,j,1)] = betanew[next++];
                cpnew[IJ(i,j,2)] = betanew[next++];
            }
        }
    }
    assert (next == nvar);

    /* for restricted Faces, adjust or reject the proposed control point movement */
    if (fit2d->rtype > 0) {
        for (j = 1; j < fit2d->nv-1; j++) {
            for (i = 1; i < fit2d->nu-1; i++) {
                double dotprod;

                dotprod = fit2d->norm[IJ(i,j,0)] * (cpnew[IJ(i,j,0)] - fit2d->cp[IJ(i,j,0)])
                        + fit2d->norm[IJ(i,j,1)] * (cpnew[IJ(i,j,1)] - fit2d->cp[IJ(i,j,1)])
                        + fit2d->norm[IJ(i,j,2)] * (cpnew[IJ(i,j,2)] - fit2d->cp[IJ(i,j,2)]);

                /* reject changes if it causes the control polygon to become concave */
                if (dotprod < 0) {
                    cpnew[IJ(i,j,0)] = fit2d->cp[IJ(i,j,0)];
                    cpnew[IJ(i,j,1)] = fit2d->cp[IJ(i,j,1)];
                    cpnew[IJ(i,j,2)] = fit2d->cp[IJ(i,j,2)];

                /* otherwise, make sure control poit movement is along local normal */
                } else {
                    cpnew[IJ(i,j,0)] = fit2d->cp[IJ(i,j,0)] + fit2d->norm[IJ(i,j,0)] * dotprod;
                    cpnew[IJ(i,j,1)] = fit2d->cp[IJ(i,j,1)] + fit2d->norm[IJ(i,j,1)] * dotprod;
                    cpnew[IJ(i,j,2)] = fit2d->cp[IJ(i,j,2)] + fit2d->norm[IJ(i,j,2)] * dotprod;
                }
            }
        }
    }

    /* compute the objective function based upon the new beta */
    status = fit2d_objf(fit2d, smooth, betanew, cpnew, fnew);
    CHECK_STATUS(fit2d_objf);

    maxfnew  = Linorm(fnew, 3*fit2d->m);
    normfnew = L2norm(fnew, nobj) / sqrt(nobj);
    if (fit2d->iter%10 == 0 && fit2d->fp != NULL) {
        fprintf(fit2d->fp, "iter=%4d normf=%10.4e, maxf=%10.4e, dotmin=%7.4f  ", fit2d->iter, normfnew, maxfnew, *dotmin);
    }

    /* if this was a better step, accept it and decrease
       lambda (making it more Newton-like) */
    if (normfnew < *normf) {
        fit2d->lambda = MAX(fit2d->lambda/2, EPS10);
        if (fit2d->iter%10 == 0 && fit2d->fp != NULL) {
        fprintf(fit2d->fp, "ACCEPTED,  lambda=%10.3e,  smooth=%10.3e\n", fit2d->lambda, smooth);
        }

        /* save new design variables, control points, and objective function */
#ifndef __clang_analyzer__
        for (k = 0; k < fit2d->m; k++) {
            fit2d->UVcloud[2*k  ] = betanew[2*k  ];
            fit2d->UVcloud[2*k+1] = betanew[2*k+1];
        }
#endif
        for (j = 0; j < fit2d->nv; j++) {
            for (i = 0; i < fit2d->nu; i++) {
                fit2d->cp[IJ(i,j,0)] = cpnew[IJ(i,j,0)];
                fit2d->cp[IJ(i,j,1)] = cpnew[IJ(i,j,1)];
                fit2d->cp[IJ(i,j,2)] = cpnew[IJ(i,j,2)];
            }
        }
#ifndef __clang_analyzer__
        for (iobj = 0; iobj < nobj; iobj++) {
            fit2d->f[iobj] = fnew[iobj];
        }
#endif
        *normf  = normfnew;
        *maxf   = maxfnew;
        *accept = 1;

    /* otherwise do not take the step and increase lambda (making it
       more steepest-descent-like) */
    } else {
        fit2d->lambda = MIN(fit2d->lambda*2, 1.0e+10);
        if (fit2d->iter%10 == 0 && fit2d->fp != NULL) {
        fprintf(fit2d->fp, "rejected,  lambda=%10.3e,  smooth=%10.3e\n", fit2d->lambda, smooth);
        }
    }

cleanup:
    FREE(fnew   );
    FREE(betanew);
    FREE(delta  );
    FREE(beta   );
    FREE(cpnew  );
    FREE(dXYZdP );

#undef MASK

    FREE(AA);
    FREE(BB);
    FREE(CC);
    FREE(DD);
    FREE(EE);
#undef A
#undef B
#undef C
#undef D
#undef E

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit2d_done - clean up after B-spline surface fitter                *
 *                                                                      *
 ************************************************************************
 */
int
fit2d_done(void    *context,            /* (in)  context */
           double  UVcloud[],           /* (out) UV-parameters of points in cloud */
           double  cp[],                /* (out) array of control points */
           double  *normf,              /* (out) RMS of distances between cloud and fit */
           double  *maxf,               /* (out) maximum distance between cloud and fit */
 /*@null@*/int     *nmin)               /* (out) minimum number of cloud points in any interval (or NULL) */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    i, j, k;
    int    *nper=NULL;

    fit2d_T *fit2d = (fit2d_T *) context;

    ROUTINE(fit2d_done);

    /* --------------------------------------------------------------- */

    /* allocate temporary array */
    MALLOC(nper, int, fit2d->nu*fit2d->nv);

    /* find the minimum number of cloud points in each spline interval */
    if (nmin != NULL) {
        *nmin = fit2d->m;

        for (j = 0; j < fit2d->nv-3; j++) {
            for (i = 0; i < fit2d->nu-3; i++) {
                nper[i+fit2d->nu*j] = 0;
            }
        }

        for (k = 0; k < fit2d->m; k++) {
            i = MIN(MAX(0,floor(fit2d->UVcloud[2*k+1])), fit2d->nu-4);
            j = MIN(MAX(0,floor(fit2d->UVcloud[2*k  ])), fit2d->nv-4);
            nper[i+fit2d->nu*j]++;
        }

        for (j = 0; j < fit2d->nv-3; j++) {
            for ( i = 0; i < fit2d->nu-3; i++) {
                if (nper[i+fit2d->nu*j] < *nmin) {
                    *nmin = nper[i+fit2d->nu*j];
                }
            }
        }
    }

    /* extract the UV parameters from the structure */
    for (k = 0; k < fit2d->m; k++) {
        UVcloud[2*k  ] = fit2d->UVcloud[2*k  ];
        UVcloud[2*k+1] = fit2d->UVcloud[2*k+1];
    }

    /* extract the control points from the structure (and un-normalize them) */
    for (j = 0; j < fit2d->nv; j++) {
        for (i = 0; i < fit2d->nu; i++) {
            cp[3*(i+j*fit2d->nu)  ] = fit2d->scale * fit2d->cp[IJ(i,j,0)] + fit2d->xavg;
            cp[3*(i+j*fit2d->nu)+1] = fit2d->scale * fit2d->cp[IJ(i,j,1)] + fit2d->yavg;
            cp[3*(i+j*fit2d->nu)+2] = fit2d->scale * fit2d->cp[IJ(i,j,2)] + fit2d->zavg;
        }
    }

    /* return the norm of the error (un-normalized) */
    *normf *= fit2d->scale;
    *maxf  *= fit2d->scale;

    /* remove temporary storage from fit2d structure */
    FREE(fit2d->XYZcloud);
    FREE(fit2d->UVcloud );
    FREE(fit2d->cp      );
    FREE(fit2d->norm    );
    FREE(fit2d->f       );
    FREE(fit2d->mask    );

    FREE(fit2d);

cleanup:
    FREE(nper);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   fit2d_objf - compute objective function                            *
 *                                                                      *
 ************************************************************************
 */
static int
fit2d_objf(fit2d_T *fit2d,              /* (in)  pointer to fit2d structure */
           double  smooth,              /* (in)  smoothing parameter */
           double  UVcloud[],           /* (in)  current UV parameters */
           double  cp[],                /* (in)  current control points */
           double  f[])                 /* (out) objective functions */
{
    int    status = FIT_SUCCESS;        /* (out)  return status */

    int    next, i, j, k;
    double XYZ[3];

    ROUTINE(fit2d_objf);

    /* --------------------------------------------------------------- */

#ifndef __clang_analyzer__
    next = 0;

    for (k = 0; k < fit2d->m; k++) {
        status = eval2dBspline(UVcloud[2*k], UVcloud[2*k+1], fit2d->nu, fit2d->nv, cp, XYZ, NULL, NULL, NULL);
        CHECK_STATUS(eval2dBspline);

        f[next++] = XYZ[0] - fit2d->XYZcloud[3*k  ];
        f[next++] = XYZ[1] - fit2d->XYZcloud[3*k+1];
        f[next++] = XYZ[2] - fit2d->XYZcloud[3*k+2];
    }

    for (j = 1; j < fit2d->nv-1; j++) {
        for (i = 1; i < fit2d->nu-1; i++) {
            f[next++] = smooth * (4 * cp[IJ(i  ,j  ,0)] - 2 * cp[IJ(i-1,j  ,0)] - 2 * cp[IJ(i+1,j  ,0)]
                                                        - 2 * cp[IJ(i  ,j-1,0)] - 2 * cp[IJ(i  ,j+1,0)]
                                                        +     cp[IJ(i-1,j-1,0)] +     cp[IJ(i+1,j-1,0)]
                                                        +     cp[IJ(i-1,j+1,0)] +     cp[IJ(i+1,j+1,0)]);
            f[next++] = smooth * (4 * cp[IJ(i  ,j  ,1)] - 2 * cp[IJ(i-1,j  ,1)] - 2 * cp[IJ(i+1,j  ,1)]
                                                        - 2 * cp[IJ(i  ,j-1,1)] - 2 * cp[IJ(i  ,j+1,1)]
                                                        +     cp[IJ(i-1,j-1,1)] +     cp[IJ(i+1,j-1,1)]
                                                        +     cp[IJ(i-1,j+1,1)] +     cp[IJ(i+1,j+1,1)]);
            f[next++] = smooth * (4 * cp[IJ(i  ,j  ,2)] - 2 * cp[IJ(i-1,j  ,2)] - 2 * cp[IJ(i+1,j  ,2)]
                                                        - 2 * cp[IJ(i  ,j-1,2)] - 2 * cp[IJ(i  ,j+1,2)]
                                                        +     cp[IJ(i-1,j-1,2)] +     cp[IJ(i+1,j-1,2)]
                                                        +     cp[IJ(i-1,j+1,2)] +     cp[IJ(i+1,j+1,2)]);
        }
    }
#endif

cleanup:
    return status;
}

#undef IJ


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
    int    status = FIT_SUCCESS;        /* (out) return status */

    int    i, span;
    double B[4], dB[4];

    ROUTINE(eval1dBspline);

    /* --------------------------------------------------------------- */

    assert (n > 3);

    XYZ[0] = 0;
    XYZ[1] = 0;
    XYZ[2] = 0;

    /* set up the Bspline bases */
    status = cubicBsplineBases(n, T, B, dB);
    CHECK_STATUS(cubicBsplineBases);

    span = MIN(floor(T), n-4);

    /* find the dependent variable */
    for (i = 0; i < 4; i++) {
        XYZ[0] += B[i] * cp[3*(i+span)  ];
        XYZ[1] += B[i] * cp[3*(i+span)+1];
        XYZ[2] += B[i] * cp[3*(i+span)+2];
    }

    /* find the deriviative wrt T */
    if (dXYZdT != NULL) {
        dXYZdT[0] = 0;
        dXYZdT[1] = 0;
        dXYZdT[2] = 0;

        for (i = 0; i < 4; i++) {
            dXYZdT[0] += dB[i] * cp[3*(i+span)  ];
            dXYZdT[1] += dB[i] * cp[3*(i+span)+1];
            dXYZdT[2] += dB[i] * cp[3*(i+span)+2];
        }
    }

    /* find the derivative wrt P */
    if (dXYZdP != NULL) {
        for (i = 0; i < n; i++) {
            dXYZdP[i] = 0;
        }
        for (i = 0; i < 4; i++) {
            dXYZdP[i+span] = B[i];
        }
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   eval2dBspline - evaluate Spline and its derivative                 *
 *                                                                      *
 ************************************************************************
 */
static int
eval2dBspline(double U,                 /* (in)  first  independent variable */
              double V,                 /* (in)  second independent variable */
              int    nu,                /* (in)  number of U control points */
              int    nv,                /* (in)  number of V control points */
              double cp[],              /* (in)  array  of control points */
              double XYZ[],             /* (out) dependent variable */
    /*@null@*/double dXYZdU[],          /* (out) derivative wrt U (or NULL) */
    /*@null@*/double dXYZdV[],          /* (out) derivative wrt V (or NULL) */
    /*@null@*/double dXYZdP[])          /* (out) derivative wrt P (or NULL) */
{
    int    status = FIT_SUCCESS;         /* (out) return status */

    int    i, j, spanu, spanv;
    double Bu[4], dBu[4], Bv[4], dBv[4];

    ROUTINE(eval2dBspline);

    /* --------------------------------------------------------------- */

    assert(nu > 3);
    assert(nv > 3);

    XYZ[0] = 0;
    XYZ[1] = 0;
    XYZ[2] = 0;

    /* set up the Bspline bases */
    status = cubicBsplineBases(nu, U, Bu, dBu);
    CHECK_STATUS(cubicBsplineBases);

    status = cubicBsplineBases(nv, V, Bv, dBv);
    CHECK_STATUS(cubicBsplineBases);

    spanu = MIN(floor(U), nu-4);
    spanv = MIN(floor(V), nv-4);

    /* find the dependent variable */
    for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) {
            XYZ[0] += Bu[i] * Bv[j] * cp[3*((i+spanu)+nu*(j+spanv))  ];
            XYZ[1] += Bu[i] * Bv[j] * cp[3*((i+spanu)+nu*(j+spanv))+1];
            XYZ[2] += Bu[i] * Bv[j] * cp[3*((i+spanu)+nu*(j+spanv))+2];
        }
    }

    /* find the derivative wrt U */
    if (dXYZdU != NULL) {
        dXYZdU[0] = 0;
        dXYZdU[1] = 0;
        dXYZdU[2] = 0;

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                dXYZdU[0] += dBu[i] * Bv[j] * cp[3*((i+spanu)+nu*(j+spanv))  ];
                dXYZdU[1] += dBu[i] * Bv[j] * cp[3*((i+spanu)+nu*(j+spanv))+1];
                dXYZdU[2] += dBu[i] * Bv[j] * cp[3*((i+spanu)+nu*(j+spanv))+2];
            }
        }
    }

    /* find the derivative wrt V */
    if (dXYZdV != NULL) {
        dXYZdV[0] = 0;
        dXYZdV[1] = 0;
        dXYZdV[2] = 0;

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                dXYZdV[0] += Bu[i] * dBv[j] * cp[3*((i+spanu)+nu*(j+spanv))  ];
                dXYZdV[1] += Bu[i] * dBv[j] * cp[3*((i+spanu)+nu*(j+spanv))+1];
                dXYZdV[2] += Bu[i] * dBv[j] * cp[3*((i+spanu)+nu*(j+spanv))+2];
            }
        }
    }

    /* find the derivative wrt P */
    if (dXYZdP != NULL) {
        for (j = 0; j < nv; j++) {
            for (i = 0; i < nu; i++) {
                dXYZdP[i+nu*j] = 0;
            }
        }

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                dXYZdP[(i+spanu)+nu*(j+spanv)] = Bu[i] * Bv[j];
            }
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
                  double B[],           /* (out) bases */
                  double dB[])          /* (out) d(bases)/d(T) */
{
    int       status = FIT_SUCCESS;     /* (out) return status */

    int      i, r, span;
    double   saved, dsaved, num, dnum, den, dden, temp, dtemp;
    double   left[4], dleft[4], rite[4], drite[4];

    ROUTINE(cubicBsplineBases);

    /* --------------------------------------------------------------- */

    span = MIN(floor(T)+3, ncp-1);

    B[ 0] = 1.0;
    dB[0] = 0;

    for (i = 1; i <= 3; i++) {
        left[ i] = T - MAX(0, span-2-i);
        dleft[i] = 1;

        rite[ i] = MIN(ncp-3,span-3+i) - T;
        drite[i] =                     - 1;

        saved  = 0;
        dsaved = 0;

        for (r = 0; r < i; r++) {
            num   = B[ r];
            dnum  = dB[r];

            den   = rite[ r+1] + left[ i-r];
            dden  = drite[r+1] + dleft[i-r];

            temp  = num / den;
            dtemp = (dnum * den - dden * num) / den / den;

            B[ r] = saved  + rite[ r+1] * temp;
            dB[r] = dsaved + drite[r+1] * temp + rite[r+1] * dtemp;

            saved  = left[ i-r] * temp;
            dsaved = dleft[i-r] * temp + left[i-r] * dtemp;
        }

        B[ i] = saved;
        dB[i] = dsaved;
    }

//cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   interp1d - 1D linear interpolation                                 *
 *                                                                      *
 ************************************************************************
 */
static int
interp1d(double T,                      /* (in)  independent variable */
         int    ntab,                   /* (in)  number of   points   in table */
         double Ttab[],                 /* (in)  array  of   ind vars in table */
         double XYZtab[],               /* (in)  array  of 3 dep vars in table */
         double XYZ[])                  /* (out) 3 dep vars at T */
{
    int    status = FIT_SUCCESS;        /* (out) return status */

    int    ileft, irite, imid;
    double frac;

    ROUTINE(interp1d);

    /* --------------------------------------------------------------- */

    /* binary search to find interval */
    ileft = 0;
    irite = ntab - 1;

    while (irite > ileft + 1) {
        imid = (ileft + irite) / 2;

        if (T < Ttab[imid]) {
            irite = imid;
        } else {
            ileft = imid;
        }
    }

    /* make sure we have ascending Ttab */
    if (Ttab[irite]-Ttab[ileft] < EPS12) {
        status = -999;
        goto cleanup;
    }

    /* linear interpolation */
    frac = (T - Ttab[ileft]) / (Ttab[irite] - Ttab[ileft]);

    XYZ[0] = (1-frac) * XYZtab[3*ileft  ] + frac * XYZtab[3*irite  ];
    XYZ[1] = (1-frac) * XYZtab[3*ileft+1] + frac * XYZtab[3*irite+1];
    XYZ[2] = (1-frac) * XYZtab[3*ileft+2] + frac * XYZtab[3*irite+2];

cleanup:
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

#ifndef __clang_analyzer__
    for (i = 0; i < n; i++) {
        L2norm += f[i] * f[i];
    }
#endif

    L2norm = sqrt(L2norm);

//cleanup:
    return L2norm;
}


/*
 ************************************************************************
 *                                                                      *
 *   Linorm - Linfinity-norm of vector                                  *
 *                                                                      *
 ************************************************************************
 */
static double
Linorm(double f[],                      /* (in)  vector */
       int    n)                        /* (in)  length of vector */
{
    double Linorm;                      /* (out) Linfinity-norm */

    int    i;
    double temp;

    ROUTINE(Linorm);

    /* --------------------------------------------------------------- */

    /* Linfinity-norm */
    Linorm = 0;

#ifndef __clang_analyzer__
    for (i = 0; i < n; i+=3) {
        temp = f[i] * f[i] + f[i+1] * f[i+1] + f[i+2] * f[i+2];
        if (temp > Linorm) {
            Linorm = temp;
        }
    }
#endif

    Linorm = sqrt(Linorm);

//cleanup:
    return Linorm;
}


/*
 ************************************************************************
 *                                                                      *
 *   fitPlane - least-square fit a*x + b*y + c*z + d to  points         *
 *                                                                      *
 ************************************************************************
 */

static int
fitPlane(double xyz[],                  /* (in)  array  of points */
         int    n,                      /* (in)  number of points */
         double *a,                     /* (out) x-coefficient */
         double *b,                     /* (out) y-coefficient */
         double *c,                     /* (out) z-coefficient */
         double *d)                     /* (out) constant */
{
    int       status = FIT_SUCCESS;     /* (out) return status */

    int       i;
    double    xcent, ycent, zcent, xdet, ydet, zdet;
    double    xx, xy, xz, yy, yz, zz;

    ROUTINE(fitPlane);

    /* --------------------------------------------------------------- */

    /* default returns */
    *a = 1;
    *b = 1;
    *c = 1;
    *d = 0;

    /* make sure we have enough points */
    if (n < 3) {
        printf("not enough points\n");
        status = -999;
        goto cleanup;
    }

    /* find the centroid of the points */
    xcent = 0;
    ycent = 0;
    zcent = 0;

    for (i = 0; i < n; i++) {
        xcent += xyz[3*i  ];
        ycent += xyz[3*i+1];
        zcent += xyz[3*i+2];
    }

    xcent /= n;
    ycent /= n;
    zcent /= n;

    /* compute the covarience matrix (relative to the controid) */
    xx = 0;
    xy = 0;
    xz = 0;
    yy = 0;
    yz = 0;
    zz = 0;

    for (i = 0; i < n; i++) {
        xx += (xyz[3*i  ] - xcent) * (xyz[3*i  ] - xcent);
        xy += (xyz[3*i  ] - xcent) * (xyz[3*i+1] - ycent);
        xz += (xyz[3*i  ] - xcent) * (xyz[3*i+2] - zcent);
        yy += (xyz[3*i+1] - ycent) * (xyz[3*i+1] - ycent);
        yz += (xyz[3*i+1] - ycent) * (xyz[3*i+2] - zcent);
        zz += (xyz[3*i+2] - zcent) * (xyz[3*i+2] - zcent);
    }

    /* find the deteminants associated with assuming normx=1,
       normy=1, and normz=1 */
    xdet = yy * zz - yz * yz;
    ydet = zz * xx - xz * xz;
    zdet = xx * yy - xy * xy;

    /* use the determinant with the best conditioning */
    if        (fabs(xdet) >= fabs(ydet) && fabs(xdet) >= fabs(zdet)) {
        *a = 1;
        *b = (xz * yz - xy * zz) / xdet;
        *c = (xy * yz - xz * yy) / xdet;
        *d = -xcent              / xdet;
    } else if (fabs(ydet) >= fabs(zdet) && fabs(ydet) >= fabs(xdet)) {
        *a = (yz * xz - xy * zz) / ydet;
        *b = 1;
        *c = (xy * yz - yz * xx) / ydet;
        *d = -ycent              / ydet;
    } else if (fabs(zdet) >= 0) {
        *a = (yz * xy - xz * yy) / zdet;
        *b = (xz * xy - yz * xx) / zdet;
        *c = 1;
        *d = -zcent              / zdet;
    } else {
        printf("cannot find big determinant\n");
        printf("xx=%e, xy=%e, xz=%e, yy=%e, yz=%e, zz=%e\n", xx, xy, xz, yy, yz, zz);
        printf("xdet=%e, ydet=%e, zdet=%e\n", xdet, ydet, zdet);
        status = -999;
        goto cleanup;
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   matsol - Gaussian elimination with partial pivoting                *
 *                                                                      *
 ************************************************************************
 */

static int
matsol(double    A[],                   /* (in)  matrix to be solved (stored rowwise) */
                                        /* (out) upper-triangular form of matrix */
       double    b[],                   /* (in)  right hand side */
                                        /* (out) right-hand side after swapping */
       int       n,                     /* (in)  size of matrix */
       double    x[])                   /* (out) solution of A*x=b */
{
    int       status = FIT_SUCCESS;     /* (out) return status */

    int       ir, jc, kc, imax;
    double    amax, swap, fact;

    ROUTINE(matsol);

    /* --------------------------------------------------------------- */

    /* reduce each column of A */
    for (kc = 0; kc < n; kc++) {

        /* find pivot element */
        imax = kc;
        amax = fabs(A[kc*n+kc]);

        for (ir = kc+1; ir < n; ir++) {
            if (fabs(A[ir*n+kc]) > amax) {
                imax = ir;
                amax = fabs(A[ir*n+kc]);
            }
        }

        /* check for possibly-singular matrix (ie, near-zero pivot) */
        if (amax < EPS12) {
            status = FIT_SINGULAR;
            goto cleanup;
        }

        /* if diagonal is not pivot, swap rows in A and b */
        if (imax != kc) {
            for (jc = 0; jc < n; jc++) {
                swap         = A[kc  *n+jc];
                A[kc  *n+jc] = A[imax*n+jc];
                A[imax*n+jc] = swap;
            }

            swap    = b[kc  ];
            b[kc  ] = b[imax];
            b[imax] = swap;
        }

        /* row-reduce part of matrix to the bottom of and right of [kc,kc] */
        for (ir = kc+1; ir < n; ir++) {
            fact = A[ir*n+kc] / A[kc*n+kc];

            for (jc = kc+1; jc < n; jc++) {
                A[ir*n+jc] -= fact * A[kc*n+jc];
            }

            b[ir] -= fact * b[kc];

            A[ir*n+kc] = 0;
        }
    }

    /* back-substitution pass */
    x[n-1] = b[n-1] / A[(n-1)*n+(n-1)];

    for (jc = n-2; jc >= 0; jc--) {
        x[jc] = b[jc];
        for (kc = jc+1; kc < n; kc++) {
            x[jc] -= A[jc*n+kc] * x[kc];
        }
        x[jc] /= A[jc*n+jc];
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotCurve - plot curve data                                        *
 *                                                                      *
 ************************************************************************
 */
#ifdef GRAFIC
int
plotCurve(int    iedge,                 /* (in)  Edge index */
          int    m,                     /* (in)  number of points in cloud */
          double XYZcloud[],            /* (in)  array  of points in cloud */
/*@null@*/double Tcloud[],              /* (in)  T-parameters of points in cloud (or NULL) */
          int    n,                     /* (in)  number of control points */
          double cp[],                  /* (in)  array  of control points */
          double normf,                 /* (in)  RMS of distances between cloud and fit */
          double dotmin,                /* (in)  minimum normalized dot product of control polygon */
          int    nmin)                  /* (in)  minimum number of control points in any interval */
{
    int    status = FIT_SUCCESS;        /* (out) return status */

    int    indgr=1+2+4+16+64+1024, ixyz=0;
    char   pltitl[255];

    ROUTINE(plotCurve);

    /* --------------------------------------------------------------- */

    snprintf(pltitl, 254, "~x~y~ iedge=%d, m=%d, n=%d, normf=%.7f, dotmin=%7.4f, nmin=%d",
             iedge, m, n, normf, dotmin, nmin);

    grctrl_(plotCurve_image, &indgr, pltitl,
            (void*)(&ixyz),
            (void*)(&m),
            (void*)(XYZcloud),
            (void*)(Tcloud),
            (void*)(&n),
            (void*)(cp),
            (void*)NULL,
            (void*)NULL,
            (void*)NULL,
            (void*)NULL,
            strlen(pltitl));

//cleanup:
    return status;
}
#endif


/*
 ************************************************************************
 *                                                                      *
 *   plotSurface - plot surface data                                    *
 *                                                                      *
 ************************************************************************
 */
#ifdef GRAFIC
int
plotSurface(int    iface,               /* (in)  Face index */
            int    m,                   /* (in)  number of points in cloud */
            double XYZcloud[],          /* (in)  array  of points in cloud */
  /*@null@*/double UVcloud[],           /* (in)  UV-parameters of points in cloud (or NULL) */
            int    n,                   /* (in)  number of control points in each direction */
            double cp[],                /* (in)  array  of control points */
            double normf,               /* (in)  RMS of distances between cloud and fit */
            double dotmin,              /* (in)  minimum normalized dot product of control polygon */
            int    nmin)                /* (in)  minimum number of cloud points in any interval */
{
    int    status = FIT_SUCCESS;        /* (out) return status */

    int    indgr=1+2+4+16+64+1024, ixyz=0, kcp=0, ijcp=0;
    char   pltitl[255];

    ROUTINE(plotSurface);

    /* --------------------------------------------------------------- */

    snprintf(pltitl, 254, "~x~y~ iface=%d, i=0", iface);

    grctrl_(plotSurface_image, &indgr, pltitl,
            (void*)(&iface),
            (void*)(&ixyz),
            (void*)(&kcp),
            (void*)(&ijcp),
            (void*)(&m),
            (void*)(XYZcloud),
            (void*)(UVcloud),
            (void*)(&n),
            (void*)(cp),
            (void*)NULL,
            strlen(pltitl));

//cleanup:
    return status;
}
#endif


/*
 ************************************************************************
 *                                                                      *
 *   plotCurve_image - plot curve data (level 3)                        *
 *                                                                      *
 ************************************************************************
 */
#ifdef GRAFIC
static void
plotCurve_image(int   *ifunct,
                void  *ixyzP,
                void  *mP,
                void  *XYZcloudP,
                void  *TcloudP,
                void  *nP,
                void  *cpP,
    /*@unused@*/void  *a6,
    /*@unused@*/void  *a7,
    /*@unused@*/void  *a8,
    /*@unused@*/void  *a9,
                float *scale,
                char  *text,
                int   textlen)
{
    int    *ixyz     = (int    *)ixyzP;
    int    *m        = (int    *)mP;
    double *XYZcloud = (double *)XYZcloudP;
    double *Tcloud   = (double *)TcloudP;
    int    *n        = (int    *)nP;
    double *cp       = (double *)cpP;

    int    i, k, status;
    float  x4, y4, z4, dum;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double TT, XYZ[3];
    char   anot[24], lablgr[81];

    int    icircle = GR_CIRCLE;
    int    iplus   = GR_PLUS;
    int    isolid  = GR_SOLID;
    int    idotted = GR_DOTTED;
    int    igreen  = GR_GREEN;
    int    iblue   = GR_BLUE;
    int    ired    = GR_RED;
    int    iblack  = GR_BLACK;

    int    izero   = 0;
    int    ione    = 1;

    ROUTINE(plot1dBspine_image);

    /* --------------------------------------------------------------- */

    if (*ifunct == 0) {
        if (*ixyz == 0) {               // y vs. x
            xmin = XYZcloud[0];
            xmax = XYZcloud[0];
            ymin = XYZcloud[1];
            ymax = XYZcloud[1];

            for (k = 1; k < *m; k++) {
                if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
                if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
                if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
                if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
            }

            scale[0] = xmin - EPS06;
            scale[1] = xmax + EPS06;
            scale[2] = ymin - EPS06;
            scale[3] = ymax + EPS06;
        } else if (*ixyz == 1) {        // z vs. y
            ymin = XYZcloud[1];
            ymax = XYZcloud[1];
            zmin = XYZcloud[2];
            zmax = XYZcloud[2];

            for (k = 1; k < *m; k++) {
                if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
                if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
                if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
                if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
            }

            scale[0] = ymin - EPS06;
            scale[1] = ymax + EPS06;
            scale[2] = zmin - EPS06;
            scale[3] = zmax + EPS06;
        } else {                        // x vs. z
            xmin = XYZcloud[0];
            xmax = XYZcloud[0];
            zmin = XYZcloud[2];
            zmax = XYZcloud[2];

            for (k = 1; k < *m; k++) {
                if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
                if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
                if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
                if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
            }

            scale[0] = zmin - EPS06;
            scale[1] = zmax + EPS06;
            scale[2] = xmin - EPS06;
            scale[3] = xmax + EPS06;
        }

        strcpy(text, "ChgXYZ");

    } else if (*ifunct == 1) {

        /* cloud of points */
        grcolr_(&igreen);

        for (k = 0; k < *m; k++) {
            x4 = XYZcloud[3*k  ];
            y4 = XYZcloud[3*k+1];
            z4 = XYZcloud[3*k+2];
            if (*ixyz == 0) {
                grmov2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grmov2_(&y4, &z4);
            } else {
                grmov2_(&z4, &x4);
            }
            grsymb_(&icircle);
        }

        /* control points */
        grcolr_(&iblue);
        grdash_(&idotted);

        x4 = cp[0];
        y4 = cp[1];
        z4 = cp[2];
        if (*ixyz == 0) {
            grmov2_(&x4, &y4);
        } else if (*ixyz == 1) {
            grmov2_(&y4, &z4);
        } else {
            grmov2_(&z4, &x4);
        }
        grsymb_(&iplus);
        snprintf(anot, 23, "~~%d", 0);
        granot_(anot, strlen(anot));

        for (i = 1; i < *n; i++) {
            x4 = cp[3*i  ];
            y4 = cp[3*i+1];
            z4 = cp[3*i+2];
            if (*ixyz == 0) {
                grdrw2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grdrw2_(&y4, &z4);
            } else {
                grdrw2_(&z4, &x4);
            }
            grsymb_(&iplus);
            snprintf(anot, 23, "~~%d", i);
            granot_(anot, strlen(anot));
        }

        /* Bspline curve */
        grcolr_(&iblack);
        grdash_(&isolid);

        x4 = cp[0];
        y4 = cp[1];
        z4 = cp[2];
        if (*ixyz == 0) {
            grmov2_(&x4, &y4);
        } else if (*ixyz == 1) {
            grmov2_(&y4, &z4);
        } else {
            grmov2_(&z4, &x4);
        }

        for (i = 1; i < 201; i++) {
            TT = (*n-3) * (double)(i) / (double)(201-1);

            status = eval1dBspline(TT, *n, cp, XYZ, NULL, NULL);
            if (status < FIT_SUCCESS) {
                printf("ERROR:: eval1dBspline -> status=%d\n", status);
            }

            x4 = XYZ[0];
            y4 = XYZ[1];
            z4 = XYZ[2];
            if (*ixyz == 0) {
                grdrw2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grdrw2_(&y4, &z4);
            } else {
                grdrw2_(&z4, &x4);
            }
        }

        /* distance from point to curve */
        if (Tcloud != NULL) {
            grcolr_(&ired);

            for (k = 0; k < *m; k++) {
                x4 = XYZcloud[3*k  ];
                y4 = XYZcloud[3*k+1];
                z4 = XYZcloud[3*k+2];
                if (*ixyz == 0) {
                    grmov2_(&x4, &y4);
                } else if (*ixyz == 1) {
                    grmov2_(&y4, &z4);
                } else {
                    grmov2_(&z4, &x4);
                }

                status = eval1dBspline(Tcloud[k], *n, cp, XYZ, NULL, NULL);
                if (status < FIT_SUCCESS) {
                    printf("ERROR:: eval1dBspline -> status=%d\n", status);
                }

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                if (*ixyz == 0) {
                    grdrw2_(&x4, &y4);
                } else if (*ixyz == 1) {
                    grdrw2_(&y4, &z4);
                } else {
                    grdrw2_(&z4, &x4);
                }
            }
        }

        grcolr_(&iblack);

    /* C option (change orientation) */
    } else if (*ifunct == -3) {
        *ixyz = (*ixyz + 1) % 3;

        grvalu_("LABLGR", &izero, &dum, lablgr, strlen("LABLGR"), 80);

        if (*ixyz == 0) {
            lablgr[1] = 'x';
            lablgr[3] = 'y';
        } else if (*ixyz == 1) {
            lablgr[1] = 'y';
            lablgr[3] = 'z';
        } else {
            lablgr[1] = 'z';
            lablgr[3] = 'x';
        }
        grvalu_("LABLGR", &ione, &dum, lablgr, strlen("LABLGR"), strlen(lablgr));
        grscpt_(&ione, "O", strlen("O"));

    } else {
        printf("ERROR:: illegal option\n");
    }
}
#endif


/*
 ************************************************************************
 *                                                                      *
 *   plotSurface_image - plot surface data (level 3)                    *
 *                                                                      *
 ************************************************************************
 */
#ifdef GRAFIC
static void
plotSurface_image(int    *ifunct,
                  void   *ifaceP,
                  void   *ixyzP,
                  void   *kcpP,
                  void   *ijcpP,
      /*@unused@*/void   *mP,
      /*@unused@*/void   *XYZcloudP,
      /*@unused@*/void   *UVcloudP,
                  void   *nP,
                  void   *cpP,
      /*@unused@*/void   *a9,
                  float  *scale,
                  char   *text,
                  int    textlen)
{
    int      *iface    = (int    *) ifaceP;
    int      *ixyz     = (int    *) ixyzP;
    int      *kcp      = (int    *) kcpP;
    int      *ijcp     = (int    *) ijcpP;
//$$$    int      *m        = (int    *) mP;
//$$$    double   *XYZcloud = (double *) XYZcloudP;
//$$$    double   *UVcloud  = (double *) UVcloudP;
    int      *n        = (int    *) nP;
    double   *cp       = (double *) cpP;

    int      status = FIT_SUCCESS;

    int      i, j;
    float    x4, y4, z4, dum;
    double   xmin, xmax, ymin, ymax, zmin, zmax, UV[2], XYZ[3];
    char     anot[24], lablgr[81];


    int      iblue    = GR_BLUE;
    int      iblack   = GR_BLACK;

    int      isolid   = GR_SOLID;
    int      idotted  = GR_DOTTED;

    int      ione     = 1;

    ROUTINE(plotSurface_image);

#undef  IJ
#define IJ(I,J,XYZ) 3*((I)+(J)*(*n))+(XYZ)

    /* --------------------------------------------------------------- */

    /* return scales */
    if (*ifunct == 0) {
        if (*ixyz == 0) {               // y vs. x
            if (*kcp == 0) {
                j = *ijcp;

                xmin = cp[IJ(0,j,0)];
                xmax = cp[IJ(0,j,0)];
                ymin = cp[IJ(0,j,1)];
                ymax = cp[IJ(0,j,1)];

                for (i = 1; i < *n; i++) {
                    xmin = MIN(xmin, cp[IJ(i,j,0)]);
                    xmax = MAX(xmax, cp[IJ(i,j,0)]);
                    ymin = MIN(ymin, cp[IJ(i,j,1)]);
                    ymax = MAX(ymax, cp[IJ(i,j,1)]);
                }
            } else {
                i = *ijcp;

                xmin = cp[IJ(i,0,0)];
                xmax = cp[IJ(i,0,0)];
                ymin = cp[IJ(i,0,1)];
                ymax = cp[IJ(i,0,1)];

                for (j = 1; j < *n; j++) {
                    xmin = MIN(xmin, cp[IJ(i,j,0)]);
                    xmax = MAX(xmax, cp[IJ(i,j,0)]);
                    ymin = MIN(ymin, cp[IJ(i,j,1)]);
                    ymax = MAX(ymax, cp[IJ(i,j,1)]);
                }
            }

            scale[0] = xmin;
            scale[1] = xmax;
            scale[2] = ymin;
            scale[3] = ymax;
        } else if (*ixyz == 1) {        // z vs. y
            if (*kcp == 0) {
                j = *ijcp;

                ymin = cp[IJ(0,j,1)];
                ymax = cp[IJ(0,j,1)];
                zmin = cp[IJ(0,j,2)];
                zmax = cp[IJ(0,j,2)];

                for (i = 1; i < *n; i++) {
                    ymin = MIN(ymin, cp[IJ(i,j,1)]);
                    ymax = MAX(ymax, cp[IJ(i,j,1)]);
                    zmin = MIN(zmin, cp[IJ(i,j,2)]);
                    zmax = MAX(zmax, cp[IJ(i,j,2)]);
                }
            } else {
                i = *ijcp;

                ymin = cp[IJ(i,0,1)];
                ymax = cp[IJ(i,0,1)];
                zmin = cp[IJ(i,0,2)];
                zmax = cp[IJ(i,0,2)];

                for (j = 1; j < *n; j++) {
                    ymin = MIN(ymin, cp[IJ(i,j,1)]);
                    ymax = MAX(ymax, cp[IJ(i,j,1)]);
                    zmin = MIN(zmin, cp[IJ(i,j,2)]);
                    zmax = MAX(zmax, cp[IJ(i,j,2)]);
                }
            }

            scale[0] = ymin;
            scale[1] = ymax;
            scale[2] = zmin;
            scale[3] = zmax;
        } else {                        // x vs. z
            if (*kcp == 0) {
                j = *ijcp;

                zmin = cp[IJ(0,j,2)];
                zmax = cp[IJ(0,j,2)];
                xmin = cp[IJ(0,j,0)];
                xmax = cp[IJ(0,j,0)];

                for (i = 1; i < *n; i++) {
                    zmin = MIN(zmin, cp[IJ(i,j,2)]);
                    zmax = MAX(zmax, cp[IJ(i,j,2)]);
                    xmin = MIN(xmin, cp[IJ(i,j,0)]);
                    xmax = MAX(xmax, cp[IJ(i,j,0)]);
                }
            } else {
                i = *ijcp;

                zmin = cp[IJ(i,0,2)];
                zmax = cp[IJ(i,0,2)];
                xmin = cp[IJ(i,0,0)];
                xmax = cp[IJ(i,0,0)];

                for (j = 1; j < *n; j++) {
                    zmin = MIN(zmin, cp[IJ(i,j,2)]);
                    zmax = MAX(zmax, cp[IJ(i,j,2)]);
                    xmin = MIN(xmin, cp[IJ(i,j,0)]);
                    xmax = MAX(xmax, cp[IJ(i,j,0)]);
                }
            }

            scale[0] = zmin;
            scale[1] = zmax;
            scale[2] = xmin;
            scale[3] = xmax;
        }

        strcpy(text, "ChgXYZ SwitchIJ NextIJ");

        goto cleanup;

    /* plot image */
    } else if (*ifunct == 1) {

        /* draw the control net */
        grcolr_(&iblue);
        grdash_(&idotted);

        if (*kcp == 0) {                // constant i
            j = *ijcp;

            x4 = cp[3*((0)+(*n)*(j))  ];
            y4 = cp[3*((0)+(*n)*(j))+1];
            z4 = cp[3*((0)+(*n)*(j))+2];
            if (*ixyz == 0) {
                grmov2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grmov2_(&y4, &z4);
            } else {
                grmov2_(&z4, &x4);
            }
            snprintf(anot, 23, "~~%d,%d", 0, j);
            granot_(anot, strlen(anot));

            for (i = 1; i < *n; i++) {
                x4 = cp[3*((i)+(*n)*(j))  ];
                y4 = cp[3*((i)+(*n)*(j))+1];
                z4 = cp[3*((i)+(*n)*(j))+2];
                if (*ixyz == 0) {
                    grdrw2_(&x4, &y4);
                } else if (*ixyz == 1) {
                    grdrw2_(&y4, &z4);
                } else {
                    grdrw2_(&z4, &x4);
                }
                snprintf(anot, 23, "~~%d,%d", i, j);
                granot_(anot, strlen(anot));
            }
        } else {                        // constant j
            i = *ijcp;

            x4 = cp[3*((i)+(*n)*(0))  ];
            y4 = cp[3*((i)+(*n)*(0))+1];
            z4 = cp[3*((i)+(*n)*(0))+2];
            if (*ixyz == 0) {
                grmov2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grmov2_(&y4, &z4);
            } else {
                grmov2_(&z4, &x4);
            }
            snprintf(anot, 23, "~~%d,%d", i, 0);
            granot_(anot, strlen(anot));

            for (j = 1; j < *n; j++) {
                x4 = cp[3*((i)+(*n)*(j))  ];
                y4 = cp[3*((i)+(*n)*(j))+1];
                z4 = cp[3*((i)+(*n)*(j))+2];
                if (*ixyz == 0) {
                    grdrw2_(&x4, &y4);
                } else if (*ixyz == 1) {
                    grdrw2_(&y4, &z4);
                } else {
                    grdrw2_(&z4, &x4);
                }
                snprintf(anot, 23, "~~%d,%d", i, j);
                granot_(anot, strlen(anot));
            }
        }

        /* Bspline iso-cut */
        grcolr_(&iblack);
        grdash_(&isolid);

        if (*kcp == 0) {                // constant i
            j = *ijcp;

            UV[0] = (double)(0);
            UV[1] = (*n-3) * (double)(j) / (double)(*n-1);

            status= eval2dBspline(UV[0], UV[1], *n, *n, cp, XYZ, NULL, NULL, NULL);
            if (status < FIT_SUCCESS) {
                printf("ERROR:: eval2dBspline -> status=%d\n", status);
            }

            x4 = XYZ[0];
            y4 = XYZ[1];
            z4 = XYZ[2];
            if (*ixyz == 0) {
                grmov2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grmov2_(&y4, &z4);
            } else {
                grmov2_(&z4, &x4);
            }

            for (i = 1; i < 201; i++) {
                UV[0] = (*n-3) * (double)(i) / (double)(201-1);

                status= eval2dBspline(UV[0], UV[1], *n, *n, cp, XYZ, NULL, NULL, NULL);
                if (status < FIT_SUCCESS) {
                    printf("ERROR:: eval2dBspline -> status=%d\n", status);
                }

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                if (*ixyz == 0) {
                    grdrw2_(&x4, &y4);
                } else if (*ixyz == 1) {
                    grdrw2_(&y4, &z4);
                } else {
                    grdrw2_(&z4, &x4);
                }
            }
        } else {                        // constant j
            i = *ijcp;

            UV[0] = (*n-3) * (double)(i) / (double)(*n-1);
            UV[1] = (double)(0);

            status= eval2dBspline(UV[0], UV[1], *n, *n, cp, XYZ, NULL, NULL, NULL);
            if (status < FIT_SUCCESS) {
                printf("ERROR:: eval2dBspline -> status=%d\n", status);
            }

            x4 = XYZ[0];
            y4 = XYZ[1];
            z4 = XYZ[2];
            if (*ixyz == 0) {
                grmov2_(&x4, &y4);
            } else if (*ixyz == 1) {
                grmov2_(&y4, &z4);
            } else {
                grmov2_(&z4, &x4);
            }

            for (j = 1; j < 201; j++) {
                UV[1] = (*n-3) * (double)(j) / (double)(201-1);

                status= eval2dBspline(UV[0], UV[1], *n, *n, cp, XYZ, NULL, NULL, NULL);
                if (status < FIT_SUCCESS) {
                    printf("ERROR:: eval2dBspline -> status=%d\n", status);
                }

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                if (*ixyz == 0) {
                    grdrw2_(&x4, &y4);
                } else if (*ixyz == 1) {
                    grdrw2_(&y4, &z4);
                } else {
                    grdrw2_(&z4, &x4);
                }
            }
        }

        grcolr_(&iblack);

        goto cleanup;

    /* C option (change orientation) */
    } else if (*ifunct == -3) {
        *ixyz = (*ixyz + 1) % 3;

    /* N option (next icp or jcp) */
    } else if (*ifunct == -14) {
        *ijcp = *ijcp + 1;

        if        (*kcp == 0 && *ijcp >= *n-1) {
            *ijcp = 0;
        } else if (*kcp == 1 && *ijcp >= *n-1) {
            *ijcp = 0;
        }

        grscpt_(&ione, "O", strlen("O"));

    /* S option (change i or j) */
    } else if (*ifunct == -19) {
        *kcp  = (*kcp + 1) % 2;
        *ijcp = 0;

        grscpt_(&ione, "O", strlen("O"));

    } else {
        printf("ERROR:: illegal option\n");
    }

    if (*ixyz == 0) {
        if (*kcp == 0) {
            snprintf(lablgr, 79, "~x~y~ iface=%d, i=%d", *iface, *ijcp);
        } else {
            snprintf(lablgr, 79, "~x~y~ iface=%d, j=%d", *iface, *ijcp);
        }
    } else if (*ixyz == 1) {
        if (*kcp == 0) {
            snprintf(lablgr, 79, "~y~z~ iface=%d, i=%d", *iface, *ijcp);
        } else {
            snprintf(lablgr, 79, "~y~z~ iface=%d, j=%d", *iface, *ijcp);
        }
    } else {
        if (*kcp == 0) {
            snprintf(lablgr, 79, "~z~x~ iface=%d, i=%d", *iface, *ijcp);
        } else {
            snprintf(lablgr, 79, "~z~x~ iface=%d, j=%d", *iface, *ijcp);
        }
    }

    grvalu_("LABLGR", &ione, &dum, lablgr, strlen("LABLGR"), strlen(lablgr));
    grscpt_(&ione, "O", strlen("O"));

cleanup:
    return;
}
#endif
