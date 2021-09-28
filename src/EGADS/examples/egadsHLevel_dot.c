#include <math.h>
#include "egads.h"
#include "egads_dot.h"

#include "../src/egadsStack.h"

#if defined(_MSC_VER) && (_MSC_VER < 1900)
#define __func__  __FUNCTION__
#endif

#define TWOPI 6.2831853071795862319959269
#define PI    (TWOPI/2.0)
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define DOT(a,b)          (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define CROSS_DOT(a_dot,b,b_dot,c,c_dot) a_dot[0] = (b_dot[1]*c[2]) + (b[1]*c_dot[2]) - (b_dot[2]*c[1]) - (b[2]*c_dot[1]);\
                                         a_dot[1] = (b_dot[2]*c[0]) + (b[2]*c_dot[0]) - (b_dot[0]*c[2]) - (b[0]*c_dot[2]);\
                                         a_dot[2] = (b_dot[0]*c[1]) + (b[0]*c_dot[1]) - (b_dot[1]*c[0]) - (b[1]*c_dot[0])

/*****************************************************************************/
/*                                                                           */
/*  pingBodies                                                               */
/*                                                                           */
/*****************************************************************************/

int
pingBodies(ego tess1, ego tess2, double dtime, int iparam, const char *shape, double ftol, double etol, double ntol)
{
  int    status = EGADS_SUCCESS;
  int    n, d, np1, np2, nt1, nt2, periodic, nerr=0;
  int    nface, nedge, nnode, iface, iedge, inode, oclass, mtype;
  double p1_dot[18], p1[18], p2[18], fd_dot[3], range1[4], range2[4], range_dot[4];
  const int    *pt1, *pi1, *pt2, *pi2, *ts1, *tc1, *ts2, *tc2;
  const double *t1, *t2, *x1, *x2, *uv1, *uv2;
  ego    ebody1, ebody2;
  ego    *efaces1=NULL, *efaces2=NULL, *eedges1=NULL, *eedges2=NULL, *enodes1=NULL, *enodes2=NULL;
  ego    top, prev, next;

  status = EG_statusTessBody( tess1, &ebody1, &np1, &np2 );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_statusTessBody( tess2, &ebody2, &np1, &np2 );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Faces from the Body 1 */
  status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, &efaces1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body 1 */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, &eedges1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Body 1 */
  status = EG_getBodyTopos(ebody1, NULL, NODE, &nnode, &enodes1);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* get the Faces from the Body 2 */
  status = EG_getBodyTopos(ebody2, NULL, FACE, &nface, &efaces2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body 2 */
  status = EG_getBodyTopos(ebody2, NULL, EDGE, &nedge, &eedges2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Body 2 */
  status = EG_getBodyTopos(ebody2, NULL, NODE, &nnode, &enodes2);
  if (status != EGADS_SUCCESS) goto cleanup;


  for (iface = 0; iface < nface; iface++) {

    /* extract the face tessellation */
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_getTessFace(tess2, iface+1, &np2, &x2, &uv2, &pt2, &pi2,
                                            &nt2, &ts2, &tc2);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (n = 0; n < np1; n++) {

      /* evaluate original edge and velocities*/
      status = EG_evaluate_dot(efaces1[iface], &uv1[2*n], NULL, p1, p1_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* evaluate perturbed edge */
      status = EG_evaluate(efaces2[iface], &uv2[2*n], p2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* compute the configuration velocity based on finite difference */
      fd_dot[0] = (p2[0] - p1[0])/dtime - p1[3]*(uv2[2*n] - uv1[2*n])/dtime - p1[6]*(uv2[2*n+1] - uv1[2*n+1])/dtime;
      fd_dot[1] = (p2[1] - p1[1])/dtime - p1[4]*(uv2[2*n] - uv1[2*n])/dtime - p1[7]*(uv2[2*n+1] - uv1[2*n+1])/dtime;
      fd_dot[2] = (p2[2] - p1[2])/dtime - p1[5]*(uv2[2*n] - uv1[2*n])/dtime - p1[8]*(uv2[2*n+1] - uv1[2*n+1])/dtime;

      for (d = 0; d < 3; d++) {
        if (fabs(p1_dot[d] - fd_dot[d]) > ftol) {
          printf("%s Face %d iparam=%d, p1[%d]=%+le fabs(%+le - %+le) = %+le > %e\n",
                 shape, iface+1, iparam, d, p1[d], p1_dot[d], fd_dot[d], fabs(p1_dot[d] - fd_dot[d]), ftol);
          nerr++;
        }
      }

      //printf("p1_dot = (%+f, %+f, %+f)\n", p1_dot[0], p1_dot[1], p1_dot[2]);
      //printf("fd_dot = (%+f, %+f, %+f)\n", fd_dot[0], fd_dot[1], fd_dot[2]);
      //printf("\n");
    }
  }

  for (iedge = 0; iedge < nedge; iedge++) {

    status = EG_getInfo(eedges1[iedge], &oclass, &mtype, &top, &prev, &next);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (mtype == DEGENERATE) continue;

    /* extract the tessellation from the original edge */
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the tessellation from the perturbed edge */
    status = EG_getTessEdge(tess2, iedge+1, &np2, &x2, &t2);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (n = 0; n < np1; n++) {

      /* evaluate original edge and velocities*/
      status = EG_evaluate_dot(eedges1[iedge], &t1[n], NULL, p1, p1_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* evaluate perturbed edge */
      status = EG_evaluate(eedges2[iedge], &t2[n], p2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* compute the configuration velocity based on finite difference */
      fd_dot[0] = (p2[0] - p1[0])/dtime - p1[3]*(t2[n] - t1[n])/dtime;
      fd_dot[1] = (p2[1] - p1[1])/dtime - p1[4]*(t2[n] - t1[n])/dtime;
      fd_dot[2] = (p2[2] - p1[2])/dtime - p1[5]*(t2[n] - t1[n])/dtime;

      for (d = 0; d < 3; d++) {
        if (fabs(p1_dot[d] - fd_dot[d]) > etol) {
          printf("%s Edge %d iparam=%d, p1[%d]=%+le fabs(%+le - %+le) = %+le > %e\n",
                 shape, iedge+1, iparam, d, p1[d], p1_dot[d], fd_dot[d], fabs(p1_dot[d] - fd_dot[d]), etol);
          nerr++;
        }
      }

      //printf("p1_dot = (%+f, %+f, %+f)\n", p1_dot[0], p1_dot[1], p1_dot[2]);
      //printf("fd_dot = (%+f, %+f, %+f)\n", fd_dot[0], fd_dot[1], fd_dot[2]);
      //printf("\n");


      /* check t-range sensitivity */
      status = EG_getRange_dot( eedges1[iedge], range1, range_dot, &periodic );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_getRange( eedges2[iedge], range2, &periodic );
      if (status != EGADS_SUCCESS) goto cleanup;

      fd_dot[0] = (range2[0] - range1[0])/dtime;
      fd_dot[1] = (range2[1] - range1[1])/dtime;

      for (d = 0; d < 2; d++) {
        if (fabs(range_dot[d] - fd_dot[d]) > etol) {
          printf("%s Edge %d iparam=%d, trng[%d]=%+le fabs(%+le - %+le) = %+le > %e\n",
                 shape, iedge+1, iparam, d, range1[d], range_dot[d], fd_dot[d], fabs(range_dot[d] - fd_dot[d]), etol);
          nerr++;
        }
      }
    }
  }

  for (inode = 0; inode < nnode; inode++) {

    /* evaluate original node and velocities*/
    status = EG_evaluate_dot(enodes1[inode], NULL, NULL, p1, p1_dot);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* evaluate perturbed edge */
    status = EG_evaluate(enodes2[inode], NULL, p2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* compute the configuration velocity based on finite difference */
    fd_dot[0] = (p2[0] - p1[0])/dtime;
    fd_dot[1] = (p2[1] - p1[1])/dtime;
    fd_dot[2] = (p2[2] - p1[2])/dtime;

    for (d = 0; d < 3; d++) {
      if (fabs(p1_dot[d] - fd_dot[d]) > etol) {
        printf("%s Node %d iparam=%d, p1[%d]=%+le fabs(%+le - %+le) = %+le > %e\n",
               shape, inode+1, iparam, d, p1[d], p1_dot[d], fd_dot[d], fabs(p1_dot[d] - fd_dot[d]), etol);
        nerr++;
      }
    }

    //printf("p1_dot = (%+f, %+f, %+f)\n", p1_dot[0], p1_dot[1], p1_dot[2]);
    //printf("fd_dot = (%+f, %+f, %+f)\n", fd_dot[0], fd_dot[1], fd_dot[2]);
    //printf("\n");
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(efaces1);
  EG_free(eedges1);
  EG_free(enodes1);

  EG_free(efaces2);
  EG_free(eedges2);
  EG_free(enodes2);

  return status + nerr;
}

/*****************************************************************************/
/*                                                                           */
/*  Re-make Topology from getTopology                                        */
/*                                                                           */
/*****************************************************************************/

int
remakeTopology(ego etopo)
{
  int    status = EGADS_SUCCESS;
  int    i, oclass, mtype, *senses, nchild, *ivec=NULL;
  double data[4], *rvec=NULL, tol, tolNew;
  ego    context, eref, egeom, eNewTopo=NULL, eNewGeom=NULL, *echild;

  status = EG_getContext(etopo, &context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getTopology(etopo, &egeom, &oclass, &mtype,
                          data, &nchild, &echild, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, egeom, oclass, mtype,
                          data, nchild, echild, senses, &eNewTopo);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_isEquivalent(etopo, eNewTopo);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getTolerance(etopo, &tol);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getTolerance(eNewTopo, &tolNew);
  if (status != EGADS_SUCCESS) goto cleanup;
  if (tolNew > 1.001*tol) {
    printf("Tolerance missmatch!! %le %le\n", tol, tolNew);
    status = EGADS_BADSCALE;
    goto cleanup;
  }

  if (egeom != NULL) {
    status = EG_getGeometry(egeom, &oclass, &mtype, &eref, &ivec, &rvec);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeGeometry(context, oclass, mtype, eref, ivec,
                             rvec, &eNewGeom);
    if (status != EGADS_SUCCESS) goto cleanup;
    EG_deleteObject(eNewGeom);
  }

  for (i = 0; i < nchild; i++) {
    status = remakeTopology(echild[i]);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

cleanup:
  EG_deleteObject(eNewTopo);

  EG_free(ivec);
  EG_free(rvec);

  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}
/*****************************************************************************/
/*                                                                           */
/*  Extrude                                                                  */
/*                                                                           */
/*****************************************************************************/

int
pingExtrudeParam(ego src, double dist, double *dir, double *params, double dtime, const char *shape, double ftol, double etol, double ntol)
{
  int    status = EGADS_SUCCESS;
  int    i, np1, nt1, iedge, nedge, iface, nface, sgn;
  double vec[4], vec_dot[4];
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  /* test extrude in both dkrections */
  for (sgn = 1; sgn >= -1; sgn -=2) {

    vec[0] = sgn*dist;
    vec[1] = dir[0];
    vec[2] = dir[1];
    vec[3] = dir[2];

    /* make the extruded body */
    status = EG_extrude(src, vec[0], &vec[1], &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Faces from the Body */
    status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Extrude %s Edge %d np1 = %d\n", shape, iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Extrude %s Face %d np1 = %d\n", shape, iface+1, np1);
    }

    /* zero out velocities */
    for (i = 0; i < 4; i++) vec_dot[i] = 0;

    for (i = 0; i < 4; i++) {

      /* set the velocity of the extruded body */
      vec_dot[i] = 1.0;
      status = EG_extrude_dot(ebody1, src, vec[0], vec_dot[0], &vec[1], &vec_dot[1]);
      if (status != EGADS_SUCCESS) goto cleanup;
      vec_dot[i] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make a perturbed body for finite difference */
      vec[i] += dtime;
      status = EG_extrude(src, vec[0], &vec[1], &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      vec[i] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, i, shape, ftol, etol, ntol);
      if (status != EGADS_SUCCESS) goto cleanup;

      EG_deleteObject(tess2);
      EG_deleteObject(ebody2);
    }

    EG_deleteObject(tess1);
    EG_deleteObject(ebody1);

  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


int
pingExtrude(ego src, double dist, double *dir, double *params, double dtime, const char *shape, double ftol, double etol, double ntol)
{
  int    status = EGADS_SUCCESS;
  int    nchld, oclass, mtype, *senses;
  double data[18];
  ego    eref, *echld;

  /* ping with the body */
  status = pingExtrudeParam(src, dist, dir, params, dtime, shape, ftol, etol, ntol);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* ping with the underlying Loop/Face directly */
  status = EG_getTopology(src, &eref, &oclass, &mtype,
                          data, &nchld, &echld, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingExtrudeParam(echld[0], dist, dir, params, dtime, shape, ftol, etol, ntol);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Line                                                                     */
/*                                                                           */
/*****************************************************************************/

int
makeLineBody( ego context,      /* (in)  EGADS context                      */
              objStack *stack,  /* (in)  EGADS obj stack                    */
              const double *x0, /* (in)  coordinates of the first point     */
              const double *x1, /* (in)  coordinates of the second point    */
              ego *ebody )      /* (out) Line wire body created from points */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[6], tdata[2];
  ego    eline, eedge, enodes[2], eloop;

  /* create Nodes for the Edge */
  data[0] = x0[0];
  data[1] = x0[1];
  data[2] = x0[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = x1[0];
  data[1] = x1[1];
  data[2] = x1[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Line (point and direction) */
  data[0] = x0[0];
  data[1] = x0[1];
  data[2] = x0[2];
  data[3] = x1[0] - x0[0];
  data[4] = x1[1] - x0[1];
  data[5] = x1[2] - x0[2];

  status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &eline);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Line */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  status = EG_makeTopology(context, eline, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, OPEN,
                           NULL, 1, &eedge, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                           NULL, 1, &eloop, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}

int
setLineBody_dot( const double *x0,     /* (in)  coordinates of the first point  */
                 const double *x0_dot, /* (in)  velocity of the first point     */
                 const double *x1,     /* (in)  coordinates of the second point */
                 const double *x1_dot, /* (in)  velocity of the second point    */
                 ego ebody )           /* (in/out) Line body with velocities    */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses;
  double data[6], data_dot[6], tdata[2], tdata_dot[2];
  ego    eline, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Line from the Edge */
  status = EG_getTopology(eedges[0], &eline, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Nodes */
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, x0, x0_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, x1, x1_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Compute the Line data and velocity */
  data[0] = x0[0];
  data[1] = x0[1];
  data[2] = x0[2];
  data[3] = x1[0] - x0[0];
  data[4] = x1[1] - x0[1];
  data[5] = x1[2] - x0[2];

  data_dot[0] = x0_dot[0];
  data_dot[1] = x0_dot[1];
  data_dot[2] = x0_dot[2];
  data_dot[3] = x1_dot[0] - x0_dot[0];
  data_dot[4] = x1_dot[1] - x0_dot[1];
  data_dot[5] = x1_dot[2] - x0_dot[2];

  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  tdata_dot[0] = 0;
  tdata_dot[1] = (data[3]*data_dot[3] + data[4]*data_dot[4] + data[5]*data_dot[5])/tdata[1];

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingLineExtrude(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[6], x_dot[6], *p1, *p2, *p1_dot, *p2_dot, params[3], dtime = 1e-7;
  double dist, dist_dot=0, dir[3], dir_dot[3]={0,0,0};
  const double *t1, *x1;
  ego    src1, src2, ebody1, ebody2, tess1, tess2;

  p1 = x;
  p2 = x+3;

  p1_dot = x_dot;
  p2_dot = x_dot+3;

  dist = 2.;
  dir[0] = 0.;
  dir[1] = 0.;
  dir[2] = 1.;

  /* make the Line body */
  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 0.50; p2[1] = 0.75; p2[2] = 1.00;
  status = makeLineBody(context, stack, p1, p2, &src1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* tessellation parameters */
  params[0] =  0.05;
  params[1] =  0.001;
  params[2] = 12.0;

  /* zero out velocities */
  for (iparam = 0; iparam < 6; iparam++) x_dot[iparam] = 0;

  /* zero out sensitivities */
  status = setLineBody_dot(p1, p1_dot, p2, p2_dot, src1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check extrude sensitivities */
  status = pingExtrude(src1, dist, dir, params, dtime, "Line", 5e-7, 5e-7, 1e-7);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the extruded body */
  status = EG_extrude(src1, dist, dir, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Line np1 = %d\n", np1);

  for (iparam = 0; iparam < 6; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setLineBody_dot(p1, p1_dot, p2, p2_dot, src1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude_dot(ebody1, src1, dist, dist_dot, dir, dir_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;


    /* make a perturbed Line for finite difference */
    x[iparam] += dtime;
    status = makeLineBody(context, stack, p1, p2, &src2);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude(src2, dist, dir, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Line", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);
  }

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Circle                                                                   */
/*                                                                           */
/*****************************************************************************/

int
makeCircleBody( ego context,         /* (in)  EGADS context        */
                objStack *stack,     /* (in)  EGADS obj stack      */
                const int btype,     /* (in)  WIREBODY or FACEBODY */
                const double *xcent, /* (in)  Center               */
                const double *xax,   /* (in)  x-axis               */
                const double *yax,   /* (in)  y-axis               */
                const double r,      /* (in)  radius               */
                ego *ebody )         /* (out) Circle body          */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD}, oclass, mtype, *ivec=NULL;
  double data[10], tdata[2], dx[3], dy[3], *rvec=NULL;
  ego    ecircle, eedge, enode, eloop, eplane, eface, eref;

  /* create the Circle */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];  /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  data[9] = r;        /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecircle, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];

  /* create the Node for the Edge */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enode);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle, EDGE, ONENODE,
                           tdata, 1, &enode, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, &eedge, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (btype == WIREBODY) {

    status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                             NULL, 1, &eloop, NULL, ebody);
    if (status != EGADS_SUCCESS) goto cleanup;

  } else {

    /* create the Plane */
    data[0] = xcent[0]; /* center */
    data[1] = xcent[1];
    data[2] = xcent[2];
    data[3] = dx[0];   /* x-axis */
    data[4] = dx[1];
    data[5] = dx[2];
    data[6] = dy[0];   /* y-axis */
    data[7] = dy[1];
    data[8] = dy[2];
    status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL,
                             data, &eplane);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eplane);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                             NULL, 1, &eloop, senses, &eface);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eface);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                             NULL, 1, &eface, NULL, ebody);
    if (status != EGADS_SUCCESS) goto cleanup;

  }
  status = EG_stackPush(stack, *ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setCircleBody_dot( const int btype,         /* (in)  WIREBODY or FACEBODY */
                   const double *xcent,     /* (in)  Center          */
                   const double *xcent_dot, /* (in)  Center velocity */
                   const double *xax,       /* (in)  x-axis          */
                   const double *xax_dot,   /* (in)  x-axis velocity */
                   const double *yax,       /* (in)  y-axis          */
                   const double *yax_dot,   /* (in)  y-axis velocity */
                   const double r,          /* (in)  radius          */
                   const double r_dot,      /* (in)  radius velocity */
                   ego ebody )              /* (in/out) Circle body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, tdata[2], tdata_dot[2];
  ego    ecircle, eplane, *efaces, *enodes, *eloops, *eedges, eref;

  if (btype == WIREBODY) {
    /* get the Loop from the Body */
    status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                            data, &nloop, &eloops, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
  } else {
    /* get the Face from the Body */
    status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                            data, &nface, &efaces, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Loop and Plane from the Face */
    status = EG_getTopology(efaces[0], &eplane, &oclass, &mtype,
                            data, &nloop, &eloops, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Circle data and velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  data[9] = r;        /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = xax_dot[0];   /* x-axis */
  data_dot[4] = xax_dot[1];
  data_dot[5] = xax_dot[2];
  data_dot[6] = yax_dot[0];   /* y-axis */
  data_dot[7] = yax_dot[1];
  data_dot[8] = yax_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(ecircle, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];
  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];
  dy_dot[0] = rvec_dot[6];
  dy_dot[1] = rvec_dot[7];
  dy_dot[2] = rvec_dot[8];

  /* set the sensitivity of the Node */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (btype == FACEBODY) {
    data[0] = xcent[0]; /* center */
    data[1] = xcent[1];
    data[2] = xcent[2];
    data[3] = dx[0];    /* x-axis */
    data[4] = dx[1];
    data[5] = dx[2];
    data[6] = dy[0];    /* y-axis */
    data[7] = dy[1];
    data[8] = dy[2];

    data_dot[0] = xcent_dot[0]; /* center */
    data_dot[1] = xcent_dot[1];
    data_dot[2] = xcent_dot[2];
    data_dot[3] = dx_dot[0];    /* x-axis */
    data_dot[4] = dx_dot[1];
    data_dot[5] = dx_dot[2];
    data_dot[6] = dy_dot[0];    /* y-axis */
    data_dot[7] = dy_dot[1];
    data_dot[8] = dy_dot[2];

    status = EG_setGeometry_dot(eplane, SURFACE, PLANE, NULL, data, data_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  status = EGADS_SUCCESS;

cleanup:
  EG_free(rvec); rvec = NULL;
  EG_free(rvec_dot); rvec_dot = NULL;

  return status;
}


int
pingCircleExtrude(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, btype[2] = {FACEBODY, WIREBODY}, i;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double dist, dist_dot=0, dir[3], dir_dot[3]={0,0,0};
  const double *t1, *x1;
  ego    src1, src2, ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  dist = 2.;
  dir[0] = 0.;
  dir[1] = 0.;
  dir[2] = 1.;

  for (i = 0; i < 2; i++) {
    /* make the Circle body */
    xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
    xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
    yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
    x[9] = 1.0;
    status = makeCircleBody(context, stack, btype[i], xcent, xax, yax, x[9], &src1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* tessellation parameters */
    params[0] =  0.1;
    params[1] =  0.1;
    params[2] = 20.0;

    /* zero out velocities */
    for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;


    /* zero out sensitivities */
    status = setCircleBody_dot(btype[i],
                               xcent, xcent_dot,
                               xax, xax_dot,
                               yax, yax_dot,
                               x[9], x_dot[9], src1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* check extrude sensitivities */
    status = pingExtrude(src1, dist, dir, params, dtime, "Circle", 5e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make the extruded body */
    status = EG_extrude(src1, dist, dir, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

#if 0
      {
        char filename[42];
        sprintf(filename, "circle.egads");
        EG_saveModel(ebody1, filename);
      }
#endif

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* extract the tessellation from the edge */
    status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;

    printf(" Circle np1 = %d\n", np1);

    for (iparam = 0; iparam < 10; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setCircleBody_dot(btype[i],
                                 xcent, xcent_dot,
                                 xax, xax_dot,
                                 yax, yax_dot,
                                 x[9], x_dot[9], src1);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_extrude_dot(ebody1, src1, dist, dist_dot, dir, dir_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      /* make a perturbed Circle for finite difference */
      x[iparam] += dtime;
      status = makeCircleBody(context, stack, btype[i], xcent, xax, yax, x[9], &src2);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_extrude(src2, dist, dir, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Circle", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;

      EG_deleteObject(tess2);
      EG_deleteObject(ebody2);
    }

    EG_deleteObject(tess1);
    EG_deleteObject(ebody1);
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Arc                                                                      */
/*                                                                           */
/*****************************************************************************/

int
makeArcBody( ego context,         /* (in)  EGADS context    */
             objStack *stack,     /* (in)  EGADS obj stack  */
             const double *xcent, /* (in)  Center           */
             const double *xax,   /* (in)  x-axis           */
             const double *yax,   /* (in)  y-axis           */
             const double r,      /* (in)  radius           */
             ego *ebody )         /* (out) Arc wire body    */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD}, oclass, mtype, *ivec=NULL;
  double data[10], tdata[2], dx[3], dy[3], *rvec=NULL;
  ego    ecircle, eedge, enodes[2], eloop, eref;

  /* create the Circle geometry*/
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];  /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  data[9] = r;        /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecircle, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];

  /* create the Nodes for the Edge */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dy[0]*r;
  data[1] = xcent[1] + dy[1]*r;
  data[2] = xcent[2] + dy[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = 0;
  tdata[1] = PI/2.;

  status = EG_makeTopology(context, ecircle, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, OPEN,
                           NULL, 1, &eedge, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                           NULL, 1, &eloop, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setArcBody_dot( const double *xcent,     /* (in)  Center          */
                const double *xcent_dot, /* (in)  Center velocity */
                const double *xax,       /* (in)  x-axis          */
                const double *xax_dot,   /* (in)  x-axis velocity */
                const double *yax,       /* (in)  y-axis          */
                const double *yax_dot,   /* (in)  y-axis velocity */
                const double r,          /* (in)  radius          */
                const double r_dot,      /* (in)  radius velocity */
                ego ebody )              /* (in/out) Arc body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, tdata[2], tdata_dot[2];
  ego    ecircle, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = PI/2.;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Circle data and velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  data[9] = r;        /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = xax_dot[0];   /* x-axis */
  data_dot[4] = xax_dot[1];
  data_dot[5] = xax_dot[2];
  data_dot[6] = yax_dot[0];   /* y-axis */
  data_dot[7] = yax_dot[1];
  data_dot[8] = yax_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(ecircle, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];
  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];
  dy_dot[0] = rvec_dot[6];
  dy_dot[1] = rvec_dot[7];
  dy_dot[2] = rvec_dot[8];

  /* set the sensitivity of the Nodes */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dy[0]*r;
  data[1] = xcent[1] + dy[1]*r;
  data[2] = xcent[2] + dy[2]*r;
  data_dot[0] = xcent_dot[0] + dy_dot[0]*r + dy[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dy_dot[1]*r + dy[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dy_dot[2]*r + dy[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(rvec); rvec = NULL;
  EG_free(rvec_dot); rvec_dot = NULL;

  return status;
}


int
pingArcExtrude(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double dist, dist_dot=0, dir[3], dir_dot[3]={0,0,0};
  const double *t1, *x1;
  ego    src1, src2, ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* direction in the plane of the arc! */
  dist = 2.;
  dir[0] = cos(45.*PI/180.);
  dir[1] = sin(45.*PI/180.);
  dir[2] = 0.4;

  /* make the Circle body */
  xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
  xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
  yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
  x[9] = 1.0;
  status = makeArcBody(context, stack, xcent, xax, yax, x[9], &src1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* tessellation parameters */
  params[0] =  0.1;
  params[1] =  0.1;
  params[2] = 20.0;

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;


  /* zero out sensitivities */
  status = setArcBody_dot(xcent, xcent_dot,
                          xax, xax_dot,
                          yax, yax_dot,
                          x[9], x_dot[9], src1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check extrude sensitivities */
  status = pingExtrude(src1, dist, dir, params, dtime, "Arc", 5e-7, 5e-7, 1e-7);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the extruded body */
  status = EG_extrude(src1, dist, dir, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Arc np1 = %d\n", np1);

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setArcBody_dot(xcent, xcent_dot,
                            xax, xax_dot,
                            yax, yax_dot,
                            x[9], x_dot[9], src1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude_dot(ebody1, src1, dist, dist_dot, dir, dir_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    /* make a perturbed Arc for finite difference */
    x[iparam] += dtime;
    status = makeArcBody(context, stack, xcent, xax, yax, x[9], &src2);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude(src2, dist, dir, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Arc", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);
  }

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Plane                                                                    */
/*                                                                           */
/*****************************************************************************/

int
makeLineEdge( ego context,      /* (in)  EGADS context           */
              objStack *stack,  /* (in)  EGADS obj stack         */
              ego n1,           /* (in)  first node              */
              ego n2,           /* (in)  second node             */
              ego *eedge )      /* (out) Edge created from nodes */
{
  int    status = EGADS_SUCCESS;
  double data[6], tdata[2], x1[3], x2[3];
  ego    eline, enodes[2];

  status = EG_evaluate(n1, NULL, x1);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate(n2, NULL, x2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Line (point and direction) */
  data[0] = x1[0];
  data[1] = x1[1];
  data[2] = x1[2];
  data[3] = x2[0] - x1[0];
  data[4] = x2[1] - x1[1];
  data[5] = x2[2] - x1[2];
  status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &eline);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Line */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  enodes[0] = n1;
  enodes[1] = n2;

  status = EG_makeTopology(context, eline, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, eedge);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
setLineEdge_dot( ego eedge )      /* (in/out) Edge with velocity */
{
  int    status = EGADS_SUCCESS;
  int    nnode, oclass, mtype, *senses;
  double data[6], data_dot[6], x1[3], x1_dot[3], x2[3], x2_dot[3], tdata[2], tdata_dot[2];
  ego    eline, *enodes;

  /* get the Nodes and the Line from the Edge */
  status = EG_getTopology(eedge, &eline, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the velocity of the nodes */
  status = EG_evaluate_dot(enodes[0], NULL, NULL, x1, x1_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate_dot(enodes[1], NULL, NULL, x2, x2_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Line data and velocity */
  data[0] = x1[0];
  data[1] = x1[1];
  data[2] = x1[2];
  data[3] = x2[0] - x1[0];
  data[4] = x2[1] - x1[1];
  data[5] = x2[2] - x1[2];

  data_dot[0] = x1_dot[0];
  data_dot[1] = x1_dot[1];
  data_dot[2] = x1_dot[2];
  data_dot[3] = x2_dot[0] - x1_dot[0];
  data_dot[4] = x2_dot[1] - x1_dot[1];
  data_dot[5] = x2_dot[2] - x1_dot[2];

  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  tdata_dot[0] = 0;
  tdata_dot[1] = (data[3]*data_dot[3] + data[4]*data_dot[4] + data[5]*data_dot[5])/tdata[1];

  status = EG_setRange_dot(eedge, EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
makePlaneBody( ego context,         /* (in)  EGADS context    */
               objStack *stack,     /* (in)  EGADS obj stack  */
               const double *xcent, /* (in)  Center           */
               const double *xax,   /* (in)  x-axis           */
               const double *yax,   /* (in)  y-axis           */
               ego *ebody )         /* (out) Plane body       */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[9], dx[3], dy[3], *rvec=NULL;
  ego    eplane, eedges[4], enodes[4], eloop, eface, eref;

  /* create the Plane */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL,
                           data, &eplane);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eplane);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(eplane, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the sphere */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];


  /* create the Nodes for the Edges */
  data[0] = xcent[0] - dx[0] - dy[0];
  data[1] = xcent[1] - dx[1] - dy[1];
  data[2] = xcent[2] - dx[2] - dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dx[0] - dy[0];
  data[1] = xcent[1] + dx[1] - dy[1];
  data[2] = xcent[2] + dx[2] - dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dx[0] + dy[0];
  data[1] = xcent[1] + dx[1] + dy[1];
  data[2] = xcent[2] + dx[2] + dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] - dx[0] + dy[0];
  data[1] = xcent[1] - dx[1] + dy[1];
  data[2] = xcent[2] - dx[2] + dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Edges */
  status =  makeLineEdge(context, stack, enodes[0], enodes[3], &eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[0], enodes[1], &eedges[1] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[1], enodes[2], &eedges[2] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[3], enodes[2], &eedges[3] );
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and the face body */
  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                           NULL, 1, &eloop, NULL, &eface);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eface);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                           NULL, 1, &eface, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec);
  EG_free(rvec);

  return status;
}


int
setPlaneBody_dot( const double *xcent,     /* (in)  Center          */
                  const double *xcent_dot, /* (in)  Center velocity */
                  const double *xax,       /* (in)  x-axis          */
                  const double *xax_dot,   /* (in)  x-axis velocity */
                  const double *yax,       /* (in)  y-axis          */
                  const double *yax_dot,   /* (in)  y-axis velocity */
                  ego ebody )              /* (in/out) Plane body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3], *rvec=NULL, *rvec_dot=NULL;
  ego    eplane, *efaces, *eloops, *eedges, enodes[4], *enode, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Plane from the Face */
  status = EG_getTopology(efaces[0], &eplane, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Edges */
  status = EG_getTopology(eedges[0], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[3] = enode[1];

  status = EG_getTopology(eedges[1], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[1] = enode[1];

  status = EG_getTopology(eedges[2], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[1] = enode[0];
  enodes[2] = enode[1];

  status = EG_getTopology(eedges[3], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[3] = enode[0];
  enodes[2] = enode[1];


  /* set the Plane data and velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = xax_dot[0];   /* x-axis */
  data_dot[4] = xax_dot[1];
  data_dot[5] = xax_dot[2];
  data_dot[6] = yax_dot[0];   /* y-axis */
  data_dot[7] = yax_dot[1];
  data_dot[8] = yax_dot[2];

  status = EG_setGeometry_dot(eplane, SURFACE, PLANE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the plane */
  status = EG_getGeometry_dot(eplane, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];

  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];
  dy_dot[0] = rvec_dot[6];
  dy_dot[1] = rvec_dot[7];
  dy_dot[2] = rvec_dot[8];

  /* create the Nodes for the Edges */
  data[0]     = xcent[0]     - dx[0]     - dy[0];
  data[1]     = xcent[1]     - dx[1]     - dy[1];
  data[2]     = xcent[2]     - dx[2]     - dy[2];
  data_dot[0] = xcent_dot[0] - dx_dot[0] - dy_dot[0];
  data_dot[1] = xcent_dot[1] - dx_dot[1] - dy_dot[1];
  data_dot[2] = xcent_dot[2] - dx_dot[2] - dy_dot[2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     + dx[0]     - dy[0];
  data[1]     = xcent[1]     + dx[1]     - dy[1];
  data[2]     = xcent[2]     + dx[2]     - dy[2];
  data_dot[0] = xcent_dot[0] + dx_dot[0] - dy_dot[0];
  data_dot[1] = xcent_dot[1] + dx_dot[1] - dy_dot[1];
  data_dot[2] = xcent_dot[2] + dx_dot[2] - dy_dot[2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     + dx[0]     + dy[0];
  data[1]     = xcent[1]     + dx[1]     + dy[1];
  data[2]     = xcent[2]     + dx[2]     + dy[2];
  data_dot[0] = xcent_dot[0] + dx_dot[0] + dy_dot[0];
  data_dot[1] = xcent_dot[1] + dx_dot[1] + dy_dot[1];
  data_dot[2] = xcent_dot[2] + dx_dot[2] + dy_dot[2];
  status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     - dx[0]     + dy[0];
  data[1]     = xcent[1]     - dx[1]     + dy[1];
  data[2]     = xcent[2]     - dx[2]     + dy[2];
  data_dot[0] = xcent_dot[0] - dx_dot[0] + dy_dot[0];
  data_dot[1] = xcent_dot[1] - dx_dot[1] + dy_dot[1];
  data_dot[2] = xcent_dot[2] - dx_dot[2] + dy_dot[2];
  status = EG_setGeometry_dot(enodes[3], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Line Edge velocities */
  status = setLineEdge_dot( eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[1] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[2] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[3] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingPlaneExtrude(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double dist, dist_dot=0, dir[3], dir_dot[3]={0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    src1, src2, ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  dist = 2.;
  dir[0] = 0.;
  dir[1] = 0.;
  dir[2] = 1.;

  /* make the Plane Face body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  status = makePlaneBody(context, stack, xcent, xax, yax, &src1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* tessellation parameters */
  params[0] =  0.5;
  params[1] =  0.1;
  params[2] = 20.0;


  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  /* zero out sensitivities */
  status = setPlaneBody_dot(xcent, xcent_dot,
                            xax, xax_dot,
                            yax, yax_dot, src1);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* check extrude sensitivities */
  status = pingExtrude(src1, dist, dir, params, dtime, "Plane", 5e-7, 5e-7, 1e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the extruded body */
  status = EG_extrude(src1, dist, dir, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the tessellation */
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* get the Faces from the Body */
  status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Plane Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Plane Face %d np1 = %d\n", iface+1, np1);
  }

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setPlaneBody_dot(xcent, xcent_dot,
                              xax, xax_dot,
                              yax, yax_dot, src1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude_dot(ebody1, src1, dist, dist_dot, dir, dir_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makePlaneBody(context, stack, xcent, xax, yax, &src2);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude(src2, dist, dir, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Plane", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);
  }

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Plane with Holes                                                         */
/*                                                                           */
/*****************************************************************************/

int
makeHolyPlaneBody( ego context,         /* (in)  EGADS context    */
                   objStack *stack,     /* (in)  EGADS obj stack  */
                   const double *xcent, /* (in)  Center           */
                   const double *xax,   /* (in)  x-axis           */
                   const double *yax,   /* (in)  y-axis           */
                   ego *ebody )         /* (out) Plane body       */
{
  int    status = EGADS_SUCCESS;
  int    oclass, mtype, *ivec=NULL;
  int    psens[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE};
  int    csens[4][2] = {{SREVERSE,0}, {SFORWARD,0}, {SREVERSE,SREVERSE}, {SFORWARD,SFORWARD}};
  int    lsens[5] = {SFORWARD, SREVERSE, SREVERSE, SREVERSE, SREVERSE};
  double data[10], tdata[2], dx[3], dy[3], r, *rvec=NULL;
  ego    eplane, ecircles[5], eedges[10], enodes[10], eloops[5], nodes[2], eface, eref;

  /* create the Plane */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL,
                           data, &eplane);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eplane);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(eplane, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the plane to make a face and circles */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];


  /* create the Nodes for the Edges */
  data[0] = xcent[0] - dx[0] - dy[0];
  data[1] = xcent[1] - dx[1] - dy[1];
  data[2] = xcent[2] - dx[2] - dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dx[0] - dy[0];
  data[1] = xcent[1] + dx[1] - dy[1];
  data[2] = xcent[2] + dx[2] - dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dx[0] + dy[0];
  data[1] = xcent[1] + dx[1] + dy[1];
  data[2] = xcent[2] + dx[2] + dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] - dx[0] + dy[0];
  data[1] = xcent[1] - dx[1] + dy[1];
  data[2] = xcent[2] - dx[2] + dy[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Edges */
  status =  makeLineEdge(context, stack, enodes[0], enodes[3], &eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[0], enodes[1], &eedges[1] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[1], enodes[2], &eedges[2] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[3], enodes[2], &eedges[3] );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the outer loop */
  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 4, eedges, psens, &eloops[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloops[0]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* radius of the Circles */
  r = 0.25;

  /* create a forward Circle */
  data[0] = xcent[0] + 2*r*dx[0] + 2*r*dy[0]; /* center */
  data[1] = xcent[1] + 2*r*dx[1] + 2*r*dy[1];
  data[2] = xcent[2] + 2*r*dx[2] + 2*r*dy[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];  /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;  /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircles[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircles[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Node for the Edge */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[4]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[4]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create a reversed Circle */
  data[0] = xcent[0] - 2*r*dx[0] + 2*r*dy[0]; /* center */
  data[1] = xcent[1] - 2*r*dx[1] + 2*r*dy[1];
  data[2] = xcent[2] - 2*r*dx[2] + 2*r*dy[2];
  data[3] =  dx[0];   /* x-axis */
  data[4] =  dx[1];
  data[5] =  dx[2];
  data[6] = -dy[0];  /* y-axis */
  data[7] = -dy[1];
  data[8] = -dy[2];
  data[9] = r;  /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircles[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircles[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Node for the Edge */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[5]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[5]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the Edges and Loops on the Circles */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircles[0], EDGE, ONENODE,
                           tdata, 1, &enodes[4], NULL, &eedges[4]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, &eedges[4], csens[0], &eloops[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloops[1]);
  if (status != EGADS_SUCCESS) goto cleanup;


  status = EG_makeTopology(context, ecircles[1], EDGE, ONENODE,
                           tdata, 1, &enodes[5], NULL, &eedges[5]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[5]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, &eedges[5], csens[1], &eloops[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloops[2]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create a forward split Circle */
  data[0] = xcent[0] - 2*r*dx[0] - 2*r*dy[0]; /* center */
  data[1] = xcent[1] - 2*r*dx[1] - 2*r*dy[1];
  data[2] = xcent[2] - 2*r*dx[2] - 2*r*dy[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];  /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;  /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircles[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircles[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Nodes for the Edges */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[6]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[6]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = data[0] - 2*dx[0]*r;
  data[1] = data[1] - 2*dx[1]*r;
  data[2] = data[2] - 2*dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[7]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[7]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Edges and Loop */
  tdata[0] = 0;
  tdata[1] = TWOPI/2.;

  nodes[0] = enodes[6];
  nodes[1] = enodes[7];
  status = EG_makeTopology(context, ecircles[3], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[6]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[6]);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = TWOPI/2.;
  tdata[1] = TWOPI;

  nodes[0] = enodes[7];
  nodes[1] = enodes[6];
  status = EG_makeTopology(context, ecircles[3], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[7]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[7]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 2, &eedges[6], csens[2], &eloops[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloops[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create a reversed split Circle */
  data[0] = xcent[0] + 2*r*dx[0] - 2*r*dy[0]; /* center */
  data[1] = xcent[1] + 2*r*dx[1] - 2*r*dy[1];
  data[2] = xcent[2] + 2*r*dx[2] - 2*r*dy[2];
  data[3] =  dx[0];   /* x-axis */
  data[4] =  dx[1];
  data[5] =  dx[2];
  data[6] = -dy[0];   /* y-axis */
  data[7] = -dy[1];
  data[8] = -dy[2];
  data[9] = r;  /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircles[4]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircles[4]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Nodes for the Edges */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[8]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[8]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = data[0] - 2*dx[0]*r;
  data[1] = data[1] - 2*dx[1]*r;
  data[2] = data[2] - 2*dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[9]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[9]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Edges and Loop */
  tdata[0] = 0;
  tdata[1] = TWOPI/2.;

  nodes[0] = enodes[8];
  nodes[1] = enodes[9];
  status = EG_makeTopology(context, ecircles[4], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[8]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[8]);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = TWOPI/2.;
  tdata[1] = TWOPI;

  nodes[0] = enodes[9];
  nodes[1] = enodes[8];
  status = EG_makeTopology(context, ecircles[4], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[9]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[9]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 2, &eedges[8], csens[3], &eloops[4]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloops[4]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the face body */
  status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                           NULL, 5, eloops, lsens, &eface);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eface);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                           NULL, 1, &eface, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec);
  EG_free(rvec);

  return status;
}


int
setHolyPlaneBody_dot( const double *xcent,     /* (in)  Center          */
                      const double *xcent_dot, /* (in)  Center velocity */
                      const double *xax,       /* (in)  x-axis          */
                      const double *xax_dot,   /* (in)  x-axis velocity */
                      const double *yax,       /* (in)  y-axis          */
                      const double *yax_dot,   /* (in)  y-axis velocity */
                      ego ebody )              /* (in/out) Plane body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double r, data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, tdata[2], tdata_dot[2];
  ego    eplane, *efaces, *eloops, *eedges, enodes[4], *enode, ecircle, eref;

  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Plane from the Face */
  status = EG_getTopology(efaces[0], &eplane, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop for the Plane */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Edges */
  status = EG_getTopology(eedges[0], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[3] = enode[1];

  status = EG_getTopology(eedges[1], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[1] = enode[1];

  status = EG_getTopology(eedges[2], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[1] = enode[0];
  enodes[2] = enode[1];

  status = EG_getTopology(eedges[3], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[3] = enode[0];
  enodes[2] = enode[1];


  /* set the Plane data and velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = xax_dot[0];   /* x-axis */
  data_dot[4] = xax_dot[1];
  data_dot[5] = xax_dot[2];
  data_dot[6] = yax_dot[0];   /* y-axis */
  data_dot[7] = yax_dot[1];
  data_dot[8] = yax_dot[2];

  status = EG_setGeometry_dot(eplane, SURFACE, PLANE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the plane */
  status = EG_getGeometry_dot(eplane, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];

  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];
  dy_dot[0] = rvec_dot[6];
  dy_dot[1] = rvec_dot[7];
  dy_dot[2] = rvec_dot[8];

  /* create the Nodes for the Edges */
  data[0]     = xcent[0]     - dx[0]     - dy[0];
  data[1]     = xcent[1]     - dx[1]     - dy[1];
  data[2]     = xcent[2]     - dx[2]     - dy[2];
  data_dot[0] = xcent_dot[0] - dx_dot[0] - dy_dot[0];
  data_dot[1] = xcent_dot[1] - dx_dot[1] - dy_dot[1];
  data_dot[2] = xcent_dot[2] - dx_dot[2] - dy_dot[2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     + dx[0]     - dy[0];
  data[1]     = xcent[1]     + dx[1]     - dy[1];
  data[2]     = xcent[2]     + dx[2]     - dy[2];
  data_dot[0] = xcent_dot[0] + dx_dot[0] - dy_dot[0];
  data_dot[1] = xcent_dot[1] + dx_dot[1] - dy_dot[1];
  data_dot[2] = xcent_dot[2] + dx_dot[2] - dy_dot[2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     + dx[0]     + dy[0];
  data[1]     = xcent[1]     + dx[1]     + dy[1];
  data[2]     = xcent[2]     + dx[2]     + dy[2];
  data_dot[0] = xcent_dot[0] + dx_dot[0] + dy_dot[0];
  data_dot[1] = xcent_dot[1] + dx_dot[1] + dy_dot[1];
  data_dot[2] = xcent_dot[2] + dx_dot[2] + dy_dot[2];
  status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     - dx[0]     + dy[0];
  data[1]     = xcent[1]     - dx[1]     + dy[1];
  data[2]     = xcent[2]     - dx[2]     + dy[2];
  data_dot[0] = xcent_dot[0] - dx_dot[0] + dy_dot[0];
  data_dot[1] = xcent_dot[1] - dx_dot[1] + dy_dot[1];
  data_dot[2] = xcent_dot[2] - dx_dot[2] + dy_dot[2];
  status = EG_setGeometry_dot(enodes[3], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Line Edge velocities */
  status = setLineEdge_dot( eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[1] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[2] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[3] );
  if (status != EGADS_SUCCESS) goto cleanup;



  /* get the Edge from the Loop for the first Circle */
  status = EG_getTopology(eloops[1], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];

  /* radius of the Circles */
  r = 0.25;

  /* forward Circle sensitivity */
  data[0] = xcent[0] + 2*r*dx[0] + 2*r*dy[0]; /* center */
  data[1] = xcent[1] + 2*r*dx[1] + 2*r*dy[1];
  data[2] = xcent[2] + 2*r*dx[2] + 2*r*dy[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */

  data_dot[0] = xcent_dot[0] + 2*r*dx_dot[0] + 2*r*dy_dot[0]; /* center */
  data_dot[1] = xcent_dot[1] + 2*r*dx_dot[1] + 2*r*dy_dot[1];
  data_dot[2] = xcent_dot[2] + 2*r*dx_dot[2] + 2*r*dy_dot[2];
  data_dot[3] = dx_dot[0];   /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];   /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = 0;           /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Node sensitivity */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;

  data_dot[0] = data_dot[0] + dx_dot[0]*r;
  data_dot[1] = data_dot[1] + dx_dot[1]*r;
  data_dot[2] = data_dot[2] + dx_dot[2]*r;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* get the Edge from the Loop for the reversed Circle */
  status = EG_getTopology(eloops[2], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];

  /* reverse Circle sensitivity */
  data[0] = xcent[0] - 2*r*dx[0] + 2*r*dy[0]; /* center */
  data[1] = xcent[1] - 2*r*dx[1] + 2*r*dy[1];
  data[2] = xcent[2] - 2*r*dx[2] + 2*r*dy[2];
  data[3] =  dx[0];   /* x-axis */
  data[4] =  dx[1];
  data[5] =  dx[2];
  data[6] = -dy[0];   /* y-axis */
  data[7] = -dy[1];
  data[8] = -dy[2];
  data[9] = r;       /* radius */

  data_dot[0] = xcent_dot[0] - 2*r*dx_dot[0] + 2*r*dy_dot[0]; /* center */
  data_dot[1] = xcent_dot[1] - 2*r*dx_dot[1] + 2*r*dy_dot[1];
  data_dot[2] = xcent_dot[2] - 2*r*dx_dot[2] + 2*r*dy_dot[2];
  data_dot[3] =  dx_dot[0];   /* x-axis */
  data_dot[4] =  dx_dot[1];
  data_dot[5] =  dx_dot[2];
  data_dot[6] = -dy_dot[0];   /* y-axis */
  data_dot[7] = -dy_dot[1];
  data_dot[8] = -dy_dot[2];
  data_dot[9] =  0;           /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Node sensitivity */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;

  data_dot[0] = data_dot[0] + dx_dot[0]*r;
  data_dot[1] = data_dot[1] + dx_dot[1]*r;
  data_dot[2] = data_dot[2] + dx_dot[2]*r;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* get the Edge from the Loop for the forward split Circle */
  status = EG_getTopology(eloops[3], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = TWOPI/2;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = TWOPI/2;
  tdata[1] = TWOPI;

  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[1] = enode[1];

  /* forward splie Circle sensitivity */
  data[0] = xcent[0] - 2*r*dx[0] - 2*r*dy[0]; /* center */
  data[1] = xcent[1] - 2*r*dx[1] - 2*r*dy[1];
  data[2] = xcent[2] - 2*r*dx[2] - 2*r*dy[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */

  data_dot[0] = xcent_dot[0] - 2*r*dx_dot[0] - 2*r*dy_dot[0]; /* center */
  data_dot[1] = xcent_dot[1] - 2*r*dx_dot[1] - 2*r*dy_dot[1];
  data_dot[2] = xcent_dot[2] - 2*r*dx_dot[2] - 2*r*dy_dot[2];
  data_dot[3] = dx_dot[0];   /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];   /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = 0;           /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Node sensitivity */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;

  data_dot[0] = data_dot[0] + dx_dot[0]*r;
  data_dot[1] = data_dot[1] + dx_dot[1]*r;
  data_dot[2] = data_dot[2] + dx_dot[2]*r;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = data[0] - 2*dx[0]*r;
  data[1] = data[1] - 2*dx[1]*r;
  data[2] = data[2] - 2*dx[2]*r;

  data_dot[0] = data_dot[0] - 2*dx_dot[0]*r;
  data_dot[1] = data_dot[1] - 2*dx_dot[1]*r;
  data_dot[2] = data_dot[2] - 2*dx_dot[2]*r;

  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;



  /* get the Edge from the Loop for the reversed split Circle */
  status = EG_getTopology(eloops[4], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = TWOPI/2;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = TWOPI/2;
  tdata[1] = TWOPI;

  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[1] = enode[1];

  /* forward splie Circle sensitivity */
  data[0] = xcent[0] + 2*r*dx[0] - 2*r*dy[0]; /* center */
  data[1] = xcent[1] + 2*r*dx[1] - 2*r*dy[1];
  data[2] = xcent[2] + 2*r*dx[2] - 2*r*dy[2];
  data[3] =  dx[0];   /* x-axis */
  data[4] =  dx[1];
  data[5] =  dx[2];
  data[6] = -dy[0];   /* y-axis */
  data[7] = -dy[1];
  data[8] = -dy[2];
  data[9] =  r;       /* radius */

  data_dot[0] = xcent_dot[0] + 2*r*dx_dot[0] - 2*r*dy_dot[0]; /* center */
  data_dot[1] = xcent_dot[1] + 2*r*dx_dot[1] - 2*r*dy_dot[1];
  data_dot[2] = xcent_dot[2] + 2*r*dx_dot[2] - 2*r*dy_dot[2];
  data_dot[3] =  dx_dot[0];   /* x-axis */
  data_dot[4] =  dx_dot[1];
  data_dot[5] =  dx_dot[2];
  data_dot[6] = -dy_dot[0];   /* y-axis */
  data_dot[7] = -dy_dot[1];
  data_dot[8] = -dy_dot[2];
  data_dot[9] =  0;           /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Node sensitivity */
  data[0] = data[0] + dx[0]*r;
  data[1] = data[1] + dx[1]*r;
  data[2] = data[2] + dx[2]*r;

  data_dot[0] = data_dot[0] + dx_dot[0]*r;
  data_dot[1] = data_dot[1] + dx_dot[1]*r;
  data_dot[2] = data_dot[2] + dx_dot[2]*r;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = data[0] - 2*dx[0]*r;
  data[1] = data[1] - 2*dx[1]*r;
  data[2] = data[2] - 2*dx[2]*r;

  data_dot[0] = data_dot[0] - 2*dx_dot[0]*r;
  data_dot[1] = data_dot[1] - 2*dx_dot[1]*r;
  data_dot[2] = data_dot[2] - 2*dx_dot[2]*r;

  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  status = EGADS_SUCCESS;

cleanup:
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingHolyPlaneExtrude(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double dist, dist_dot=0, dir[3], dir_dot[3]={0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    src1, src2, ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  dist = 2.;
  dir[0] = 0.;
  dir[1] = 0.;
  dir[2] = 1.;

  /* make the Plane Face body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
//  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
//  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
  yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
  status = makeHolyPlaneBody(context, stack, xcent, xax, yax, &src1);
  if (status != EGADS_SUCCESS) goto cleanup;


//  ego model;
//  status = EG_makeTopology(context, NULL, MODEL, 0,
//                           NULL, 1, &src1, NULL, &model);
//  EG_saveModel(model, "holyPlane.egads");

  /* make the extruded body */

//  status = EG_extrude(src1, dist, dir, &ebody1);
//  if (status != EGADS_SUCCESS) goto cleanup;
//
//
//  ego model;
//  status = EG_makeTopology(context, NULL, MODEL, 0,
//                           NULL, 1, &ebody1, NULL, &model);
//  EG_saveModel(model, "holyPlane.egads");



  /* tessellation parameters */
  params[0] =  0.5;
  params[1] =  0.1;
  params[2] = 20.0;


  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  /* zero out sensitivities */
  status = setHolyPlaneBody_dot(xcent, xcent_dot,
                                xax, xax_dot,
                                yax, yax_dot, src1);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* check extrude sensitivities */
  status = pingExtrude(src1, dist, dir, params, dtime, "HolyPlane", 5e-7, 5e-7, 1e-7);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the extruded body */
  status = EG_extrude(src1, dist, dir, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the tessellation */
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* get the Faces from the Body */
  status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" HolyPlane Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" HolyPlane Face %d np1 = %d\n", iface+1, np1);
  }

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setHolyPlaneBody_dot(xcent, xcent_dot,
                                  xax, xax_dot,
                                  yax, yax_dot, src1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude_dot(ebody1, src1, dist, dist_dot, dir, dir_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeHolyPlaneBody(context, stack, xcent, xax, yax, &src2);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_extrude(src2, dist, dir, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "HolyPlane", 5e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);
  }

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


int main(int argc, char *argv[])
{
  int status, i, oclass, mtype;
  ego context, ref, prev, next;
  objStack stack;

  /* create an EGADS context */
  status = EG_open(&context);
  if (status != EGADS_SUCCESS) {
    printf(" EG_open return = %d\n", status);
    return EXIT_FAILURE;
  }

  /* create stack for gracefully cleaning up objects */
  status  = EG_stackInit(&stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingLineExtrude(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingCircleExtrude(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingArcExtrude(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingPlaneExtrude(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingHolyPlaneExtrude(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if (i != EGADS_SUCCESS)
      printf(" EGADS Internal: EG_deleteObject = %d!\n", i);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  /* check to make sure the context is clean */
  EG_getInfo(context, &oclass, &mtype, &ref, &prev, &next);
  if (next != NULL) {
    status = EGADS_CONSTERR;
    printf("Context is not properly clean!\n");
  }

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
