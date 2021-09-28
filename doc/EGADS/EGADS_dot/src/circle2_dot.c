
#include "egads.h"
#include "egads_dot.h"
#include <math.h>

#define TWOPI 6.2831853071795862319959269

int
pingBodies(ego tess1, ego tess2, double dtime, int iparam, const char *shape, double ftol, double etol, double ntol);

/*
 **********************************************************************
 *                                                                    *
 *   makeCircleBody - Creates a Circle WIREBODY                       *
 *                                                                    *
 **********************************************************************
 */

int
makeCircleBody( ego context,         /* (in)  EGADS context    */
                const double *xcent, /* (in)  Center           */
                const double *xaxis, /* (in)  x-axis           */
                const double *yaxis, /* (in)  y-axis           */
                const double r,      /* (in)  radius           */
                ego *ebody )         /* (out) Circle wire body */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD}, oclass, mtype, *ivec=NULL;
  double data[10], tdata[2], dx[3], *rvec=NULL;
  ego    ecircle, eedge, enode, eloop, eref;

  /* create the Circle */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xaxis[0]; /* x-axis */
  data[4] = xaxis[1];
  data[5] = xaxis[2];
  data[6] = yaxis[0]; /* y-axis */
  data[7] = yaxis[1];
  data[8] = yaxis[2];
  data[9] = r;        /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecircle, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the orthonormal x-axis from the circle */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];

  /* create the Node for the Edge */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle, EDGE, ONENODE,
                           tdata, 1, &enode, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, &eedge, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                           NULL, 1, &eloop, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in makeCircleBody\n", status);
  }
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}

/*
 **********************************************************************
 *                                                                    *
 *   setCircleBody_dot - Set sensitivities on a Circle                *
 *                       by building a 2nd circle                     *
 *                                                                    *
 **********************************************************************
 */

int
setCircleBody_dot( ego ebody,       /* (in/out) body with sensitivities */
                   const double *xcent,     /* (in)  Center             */
                   const double *xcent_dot, /* (in)  Center sensitivity */
                   const double *xaxis,     /* (in)  x-axis             */
                   const double *xaxis_dot, /* (in)  x-axis sensitivity */
                   const double *yaxis,     /* (in)  y-axis             */
                   const double *yaxis_dot, /* (in)  y-axis sensitivity */
                   const double r,          /* (in)  radius             */
                   const double r_dot)      /* (in)  radius sensitivity */
{
  int    status = EGADS_SUCCESS;
  int    oclass, mtype, senses[1] = {SFORWARD}, *ivec=NULL;
  double data[10], data_dot[10], dx[3], dx_dot[3];
  double tdata[2], tdata_dot[2], *rvec=NULL, *rvec_dot=NULL;
  ego    context, ecircle, enode, eloop, eedge, ebody2, eref;

  status = EG_getContext(ebody, &context);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* the Circle data and sensitivity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xaxis[0]; /* x-axis */
  data[4] = xaxis[1];
  data[5] = xaxis[2];
  data[6] = yaxis[0]; /* y-axis */
  data[7] = yaxis[1];
  data[8] = yaxis[2];
  data[9] = r;        /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = xaxis_dot[0]; /* x-axis */
  data_dot[4] = xaxis_dot[1];
  data_dot[5] = xaxis_dot[2];
  data_dot[6] = yaxis_dot[0]; /* y-axis */
  data_dot[7] = yaxis_dot[1];
  data_dot[8] = yaxis_dot[2];
  data_dot[9] = r_dot;        /* radius */

  /* create the Circle */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Circle sensitivity */
  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecircle, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(rvec);

  status = EG_getGeometry_dot(ecircle, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the orthonormal x-axis from the circle */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];

  /* the Node and it's sensitvities */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;

  /* create the Node for the Edge */
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitvity of the Node */
  status = EG_setGeometry_dot(enode, NODE, 0, NULL , data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle, EDGE, ONENODE,
                           tdata, 1, &enode, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedge, EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, &eedge, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                           NULL, 1, &eloop, NULL, &ebody2);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_copyGeometry_dot(ebody2, NULL, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(ebody2);

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in makeCircleBody_dot\n", status);
  }
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;
  EG_free(rvec_dot); rvec_dot = NULL;

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
  int    iparam, np1;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const double *t1, *x1;
  ego    context, ebody1, ebody2, tess1, tess2;

  status = EG_open(&context);
  if (status != EGADS_SUCCESS) goto cleanup;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* make the Circle body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  x[9] = 1.0;
  status = makeCircleBody(context, xcent, xax, yax, x[9], &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Circle np1 = %d\n", np1);

  /* zero out sensitivities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the analytic sensitivity of the body */
    x_dot[iparam] = 1.0;
    status = setCircleBody_dot(ebody1,
                               xcent , xcent_dot,
                               xax   , xax_dot,
                               yax   , yax_dot,
                               x[9]  , x_dot[9]);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeCircleBody(context, xcent, xax, yax, x[9], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam,
                        "Circle", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);
  }

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);

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
