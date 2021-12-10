
#include "egads.h"
#include "egads_dot.h"
#include <math.h>

#define TWOPI 6.2831853071795862319959269

int
pingBodies(ego tess1, ego tess2, double dtime, int iparam, const char *shape, double ftol, double etol, double ntol);

/*
 **********************************************************************
 *                                                                    *
 *   makeNacaBody - Splined NACA airfoil FACEBODY                     *
 *                                                                    *
 **********************************************************************
 */
/* number of points for the NACA spline and spline fit tolerance */
#define NUMPNTS    101
#define DXYTOL     1.0e-8

/* set KNOTS to 0 for arc-length knots, and -1 for equally spaced knots */
#define KNOTS      0

/* NACA 4 series coefficients */
static const double A =  0.2969;
static const double B = -0.1260;
static const double C = -0.3516;
static const double D =  0.2843;
static const double Eb= -0.1015; /* blunt TE */
static const double Es= -0.1036; /* sharp TE */

int makeNacaBody( ego context,       /* (in) EGADS context     */
                  const int sharpte, /* (in) sharp or blunt TE */
                  const double m,    /* (in) camber            */
                  const double p,    /* (in) maxloc            */
                  const double t,    /* (in) thickness         */
                  ego *ebody )      /* (out) result pointer    */
{
  int     status = EGADS_SUCCESS;
  int     ipnt, *header=NULL, sizes[2], sense[3], nedge, oclass, mtype;
  double  *pnts = NULL, *rvec = NULL;
  double  data[18], tdata[2], tle;
  double  x, y, zeta, s, yt, yc, theta, ycm, dycm;
  ego     eref, enodes[4], eedges[3], ecurve, eline, eloop, eplane, eface;

  /* mallocs required by Windows compiler */
  pnts = (double*)EG_alloc((3*NUMPNTS)*sizeof(double));
  if (pnts == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }

  /* points around airfoil (upper and then lower) */
  for (ipnt = 0; ipnt < NUMPNTS; ipnt++) {
    zeta = TWOPI * ipnt / (NUMPNTS-1);
    s    = (1 + cos(zeta)) / 2;

    if (sharpte == 0) {
      yt = 5.*t * (A * sqrt(s) + s * (B + s * (C + s * (D + s * Eb))));
    } else {
      yt = 5.*t * (A * sqrt(s) + s * (B + s * (C + s * (D + s * Es))));
    }

    if (s < p) {
      ycm  = (s * (2*p -   s)) / (p*p);
      dycm = (    (2*p - 2*s)) / (p*p);
    } else {
      ycm  = ((1-2*p) + s * (2*p -   s)) / pow(1-p,2);
      dycm = (              (2*p - 2*s)) / pow(1-p,2);
    }
    yc    = m * ycm;
    theta = atan(m * dycm);

    if (ipnt < NUMPNTS/2) {
      x = s  - yt * sin(theta);
      y = yc + yt * cos(theta);
    } else if (ipnt == NUMPNTS/2) {
      x = 0.;
      y = 0.;
    } else {
      x = s  + yt * sin(theta);
      y = yc - yt * cos(theta);
    }

    pnts[3*ipnt  ] = x;
    pnts[3*ipnt+1] = y;
    pnts[3*ipnt+2] = 0.;
  }

  /* create spline curve from upper TE, to LE, to lower TE
   *
   * finite difference must use knots equally spaced (sizes[1] == -1)
   * arc-length based knots (sizes[1] == 0) causes the t-space to change.
   */
  sizes[0] = NUMPNTS;
  sizes[1] = KNOTS;
  status = EG_approximate(context, 0, DXYTOL, sizes, pnts, &ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Node at upper trailing edge */
  ipnt = 0;
  data[0] = pnts[3*ipnt  ];
  data[1] = pnts[3*ipnt+1];
  data[2] = pnts[3*ipnt+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* node at leading edge as a function of the spline */
  status = EG_getGeometry(ecurve, &oclass, &mtype, &eref, &header, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  ipnt = (NUMPNTS - 1) / 2 + 3; /* index, with knot offset of 3 (cubic)*/
  tle = rvec[ipnt];             /* t-value (should be very close to (0,0,0) */

  status = EG_evaluate(ecurve, &tle, data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (sharpte == 0) {
    /* create Node at lower trailing edge */
    ipnt = NUMPNTS - 1;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    status = EG_makeTopology(context, NULL, NODE, 0,
                             data, 0, NULL, NULL, &enodes[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    enodes[3] = enodes[0];
  } else {
    enodes[2] = enodes[0];
  }

  /* make Edge for upper surface */
  tdata[0] = 0;   /* t-value at lower TE */
  tdata[1] = tle;

  /* construct the upper Edge */
  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                           tdata, 2, &enodes[0], NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make Edge for lower surface */
  tdata[0] = tdata[1]; /* t-value at leading edge */
  tdata[1] = 1;        /* t value at upper TE */

  /* construct the lower Edge */
  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                           tdata, 2, &enodes[1], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (sharpte == 0) {
    nedge = 3;

    /* create line segment at trailing edge */
    ipnt = NUMPNTS - 1;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    data[3] = pnts[0] - data[0];
    data[4] = pnts[1] - data[1];
    data[5] = pnts[2] - data[2];

    status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, data, &eline);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make Edge for this line */
    tdata[0] = 0;
    tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

    status = EG_makeTopology(context, eline, EDGE, TWONODE,
                             tdata, 2, &enodes[2], NULL, &eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;
  } else {
    nedge = 2;
  }

  /* create loop of the Edges */
  sense[0] = SFORWARD;
  sense[1] = SFORWARD;
  sense[2] = SFORWARD;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, nedge, eedges, sense, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create a plane for the loop */
  data[0] = 0.;
  data[1] = 0.;
  data[2] = 0.;
  data[3] = 1.; data[4] = 0.; data[5] = 0.;
  data[6] = 0.; data[7] = 1.; data[8] = 0.;

  status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, data, &eplane);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Face from the plane and the Loop */
  status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                           NULL, 1, &eloop, sense, &eface);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the FaceBody */
  status = EG_makeTopology(context, eplane, BODY, FACEBODY,
                           NULL, 1, &eface, sense, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in makeNacaBody\n", status);
  }
  EG_free(header);
  EG_free(rvec);
  EG_free(pnts);

  return status;
}

/*
 **********************************************************************
 *                                                                    *
 *   setNacaBody_dot - Set sensitivities on a NACA spline             *
 *                       created with makeNacaBody                    *
 *                                                                    *
 **********************************************************************
 */

int setNacaBody_dot( ego eobj,  /* (in/out) body with sensitivities    */
                     const int sharpte,  /* (in) sharp or blunt TE     */
                     const double m,     /* (in) camber                */
                     const double m_dot, /* (in) camber sensitivity    */
                     const double p,     /* (in) maxloc                */
                     const double p_dot, /* (in) maxloc sensitivity    */
                     const double t,     /* (in) thickness             */
                     const double t_dot) /* (in) thickness sensitivity */

{
  int     status = EGADS_SUCCESS;
  int     ipnt, nedge, oclass, mtype, nchild, nloop, nface, *senses, sizes[2];
  double  data[18], data_dot[18], tdata[2], tdata_dot[2];
  double  zeta, s, *pnts=NULL, *pnts_dot=NULL;
  double  yt, yt_dot, yc, yc_dot, theta, theta_dot;
  double  ycm, ycm_dot, dycm, dycm_dot, tle, tle_dot;
  double  x, x_dot, y, y_dot, *rvec=NULL, *rvec_dot=NULL;;
  ego     enodes[3], *efaces, *eedges, ecurve, eline, *eloops, eplane, eref, *echildren;

  /* get the Face from the FaceBody */
  status = EG_getTopology(eobj, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Plane and Loop from the face */
  status = EG_getTopology(efaces[0], &eplane, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the edges */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype, data, &nedge, &eedges,
                          &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the nodes and the curve from the first edge */
  status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype, data, &nchild, &echildren,
                          &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = echildren[0]; // upper trailing edge
  enodes[1] = echildren[1]; // leading edge

  /* mallocs required by Windows compiler */
  pnts     = (double*)EG_alloc((3*NUMPNTS)*sizeof(double));
  pnts_dot = (double*)EG_alloc((3*NUMPNTS)*sizeof(double));
  if (pnts == NULL || pnts_dot == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }

  /* points around airfoil (upper and then lower) */
  for (ipnt = 0; ipnt < NUMPNTS; ipnt++) {
    zeta = TWOPI * ipnt / (NUMPNTS-1);
    s    = (1 + cos(zeta)) / 2;

    if (sharpte == 0) {
      yt     = 5.*t     * (A * sqrt(s) + s * (B + s * (C + s * (D + s * Eb))));
      yt_dot = 5.*t_dot * (A * sqrt(s) + s * (B + s * (C + s * (D + s * Eb))));
    } else {
      yt     = 5.*t     * (A * sqrt(s) + s * (B + s * (C + s * (D + s * Es))));
      yt_dot = 5.*t_dot * (A * sqrt(s) + s * (B + s * (C + s * (D + s * Es))));
    }

    if (s < p) {
      ycm      =         ( s * (2*p -   s)) / (p*p);
      ycm_dot  = p_dot * (-2 * s * (p - s)) / (p*p*p);
      dycm     =         (     (2*p - 2*s)) / (p*p);
      dycm_dot = p_dot * (-2 * (p   - 2*s)) / (p*p*p);
    } else {
      ycm      =         ( (1-2*p) + s * (2*p -   s)) / pow(1-p,2);
      ycm_dot  = p_dot * (         2*(s - p)*(s - 1)) / pow(p-1,3);
      dycm     =         (               (2*p - 2*s)) / pow(1-p,2);
      dycm_dot = p_dot * (          -2*(1 + p - 2*s)) / pow(p-1,3);
    }
    yc        = m * ycm;
    yc_dot    = m_dot * ycm + m * ycm_dot;
    theta     = atan(m * dycm);
    theta_dot = (m_dot * dycm + m * dycm_dot) / (1 + m*m*dycm*dycm);

    if (ipnt < NUMPNTS/2) {
      x     = s  - yt * sin(theta);
      y     = yc + yt * cos(theta);
      x_dot =        - yt_dot * sin(theta) - theta_dot * yt * cos(theta);
      y_dot = yc_dot + yt_dot * cos(theta) - theta_dot * yt * sin(theta);
    } else if (ipnt == NUMPNTS/2) {
      x     = 0;
      y     = 0;
      x_dot = 0;
      y_dot = 0;
    } else {
      x     = s  + yt * sin(theta);
      y     = yc - yt * cos(theta);
      x_dot =        + yt_dot * sin(theta) + theta_dot * yt * cos(theta);
      y_dot = yc_dot - yt_dot * cos(theta) + theta_dot * yt * sin(theta);
    }

    pnts[3*ipnt  ] = x;
    pnts[3*ipnt+1] = y;
    pnts[3*ipnt+2] = 0.;

    pnts_dot[3*ipnt  ] = x_dot;
    pnts_dot[3*ipnt+1] = y_dot;
    pnts_dot[3*ipnt+2] = 0.;
  }

  /* populate spline curve sensitivities */
  sizes[0] = NUMPNTS;
  sizes[1] = KNOTS;
  status = EG_approximate_dot(ecurve, 0, DXYTOL, sizes, pnts, pnts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Node at trailing edge */
  ipnt = 0;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL,
                              &pnts[3*ipnt], &pnts_dot[3*ipnt]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Node at leading edge */
  status = EG_getGeometry_dot(ecurve, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  ipnt = (NUMPNTS - 1) / 2 + 3; /* index, with knot offset of 3 (cubic)*/
  tle     = rvec[ipnt];         /* t-value (should be very close to (0,0,0) */
  tle_dot = rvec_dot[ipnt];     /* t-value sensitivity */

  status = EG_evaluate_dot(ecurve, &tle, &tle_dot, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Edge t-range sensitivity for upper surface */
  tdata[0]     = 0;
  tdata[1]     = tle;
  tdata_dot[0] = 0;
  tdata_dot[1] = tle_dot;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Edge t-range sensitivity for lower surface */
  tdata[0]     = tdata[1];
  tdata[1]     = 1;
  tdata_dot[0] = tdata_dot[1];
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  if (sharpte == 0) {
    /* trailing edge line and lower trailing edge node from the 3rd edge */
    status = EG_getTopology(eedges[2], &eline, &oclass, &mtype, data,
                            &nchild, &echildren, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    enodes[2] = echildren[0]; /* lower trailing edge */

    /* set the sensitivity of the Node at lower trailing edge */
    ipnt = NUMPNTS - 1;
    status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL,
                                &pnts[3*ipnt], &pnts_dot[3*ipnt]);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity of the line segment at trailing edge */
    ipnt = NUMPNTS - 1;
    data[0] = pnts[3*ipnt  ];
    data[1] = pnts[3*ipnt+1];
    data[2] = pnts[3*ipnt+2];
    data[3] = pnts[0] - data[0];
    data[4] = pnts[1] - data[1];
    data[5] = pnts[2] - data[2];

    data_dot[0] = pnts_dot[3*ipnt  ];
    data_dot[1] = pnts_dot[3*ipnt+1];
    data_dot[2] = pnts_dot[3*ipnt+2];
    data_dot[3] = pnts_dot[0] - data_dot[0];
    data_dot[4] = pnts_dot[1] - data_dot[1];
    data_dot[5] = pnts_dot[2] - data_dot[2];

    status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set Edge t-range sensitivity */
    tdata[0] = 0;
    tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

    tdata_dot[0] = 0;
    tdata_dot[1] = (data[3]*data_dot[3] + data[4]*data_dot[4] + data[5]*data_dot[5])/tdata[1];

    status = EG_setRange_dot(eedges[2], EDGE, tdata, tdata_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  /* plane data */
  data[0] = 0.;
  data[1] = 0.;
  data[2] = 0.;
  data[3] = 1.; data[4] = 0.; data[5] = 0.;
  data[6] = 0.; data[7] = 1.; data[8] = 0.;

  /* set the sensitivity of the plane */
  data_dot[0] = 0.;
  data_dot[1] = 0.;
  data_dot[2] = 0.;
  data_dot[3] = 0.; data_dot[4] = 0.; data_dot[5] = 0.;
  data_dot[6] = 0.; data_dot[7] = 0.; data_dot[8] = 0.;

  status = EG_setGeometry_dot(eplane, SURFACE, PLANE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in setNacaBody_dot\n", status);
  }
  EG_free(rvec);
  EG_free(rvec_dot);
  EG_free(pnts);
  EG_free(pnts_dot);

  return status;
}


/*********************************************************************/
/*                                                                   */
/*   main - main program                                             */
/*                                                                   */
/*********************************************************************/

int
main(int       argc,                    /* (in)  number of arguments */
     char      *argv[])                 /* (in)  array of arguments */
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  int    sharpte = 0;
  double x[3], x_dot[3], params[3], dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    context, ebody1, ebody2, tess1, tess2;
  enum naca { im, ip, it };

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  status = EG_open(&context);
  if (status != EGADS_SUCCESS) goto cleanup;

  for (sharpte = 0; sharpte < 2; sharpte++) {
    printf("\n sharpte = %d\n", sharpte);

    /* make a NACA body */
    status = makeNacaBody( context, sharpte, x[im], x[ip], x[it], &ebody1 );
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the tessellation */
    params[0] =  0.1;
    params[1] =  0.01;
    params[2] = 20.0;
    status = EG_makeTessBody(ebody1, params, &tess1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Faces from the Body */
    status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping NACA Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping NACA Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out sensitivities */
    for (iparam = 0; iparam < 3; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 3; iparam++) {

      /* set the analytic sensitivity of the body */
      x_dot[iparam] = 1.0;
      status = setNacaBody_dot(ebody1,
                               sharpte,
                               x[im], x_dot[im],
                               x[ip], x_dot[ip],
                               x[it], x_dot[it]);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make a perturbed Circle for finite difference */
      x[iparam] += dtime;
      status = makeNacaBody(context, sharpte, x[im], x[ip], x[it], &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam,
                          "Naca", 1e-7, 5e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;

      EG_deleteObject(tess2);
      EG_deleteObject(ebody2);
    }

    EG_deleteObject(tess1);
    EG_deleteObject(ebody1);
  }

cleanup:

  EG_close(context);

  /* these statements are in case we used an error return to go to cleanup */
  if (status != EGADS_SUCCESS) {
    printf(" Overall Failure %d\n", status);
    status = EXIT_FAILURE;
  } else {
    printf(" EGADS_SUCCESS!\n");
    status = EXIT_SUCCESS;
  }

  return status;
}
