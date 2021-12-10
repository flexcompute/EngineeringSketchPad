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
  int    n, d, np1, np2, nt1, nt2, nerr=0;
  int    nface, nedge, nnode, iface, iedge, inode, oclass, mtype, periodic;
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
    }

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
  double data[4], *rvec=NULL;
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
  EG_deleteObject(eNewTopo);

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
  EG_free(ivec);
  EG_free(rvec);

  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Box                                                                      */
/*                                                                           */
/*****************************************************************************/

int
pingBox(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double data[6], data_dot[6], params[3], dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  data[0] = 4; /* base coordinate */
  data[1] = 5;
  data[2] = 6;
  data[3] = 1; /* dx, dy, dz sizes */
  data[4] = 2;
  data[5] = 3;

  /* make the body */
  status = EG_makeSolidBody(context, BOX, data, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Faces from the Body */
  status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.2;
  params[1] =  0.01;
  params[2] = 12.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" BOX Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" BOX Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 6; iparam++) data_dot[iparam] = 0;

  for (iparam = 0; iparam < 6; iparam++) {

    /* set the velocity of the original body */
    data_dot[iparam] = 1.0;
    status = EG_makeSolidBody_dot(ebody1, BOX, data, data_dot);
    if (status != EGADS_SUCCESS) goto cleanup;
    data_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    data[iparam] += dtime;
    status = EG_makeSolidBody(context, BOX, data, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    data[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "BOX", 1e-7, 1e-7, 1e-7);
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
/*  Sphere                                                                   */
/*                                                                           */
/*****************************************************************************/

int
pingSphere(ego context)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, sgn;
  double data[4], data_dot[4], params[3], dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  data[0] = 4; /* center */
  data[1] = 5;
  data[2] = 6;
  data[3] = 2; /* radius */

  for (sgn = -1; sgn <= 1; sgn += 2) {
    /* make the body */
    status = EG_makeSolidBody(context, sgn*SPHERE, data, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Faces from the Body */
    status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the tessellation */
    params[0] =  0.4;
    params[1] =  0.2;
    params[2] = 20.0;
    status = EG_makeTessBody(ebody1, params, &tess1);

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" SPHERE Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" SPHERE Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 4; iparam++) data_dot[iparam] = 0;

    for (iparam = 0; iparam < 4; iparam++) {

      /* set the velocity of the original body */
      data_dot[iparam] = 1.0;
      status = EG_makeSolidBody_dot(ebody1, sgn*SPHERE, data, data_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      data_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      /* make a perturbed body for finite difference */
      data[iparam] += dtime;
      status = EG_makeSolidBody(context, sgn*SPHERE, data, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      data[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "SPHERE", 1e-7, 1e-7, 1e-7);
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
/*  Cone                                                                     */
/*                                                                           */
/*****************************************************************************/

int
pingCone(ego context)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, sgn;
  double data[7], data_dot[7], params[3], dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  data[0] = 5; /* vertex */
  data[1] = 4;
  data[2] = 6;
  data[3] = 1; /* base */
  data[4] = 2;
  data[5] = 3;
  data[6] = 4; /* radius */

  for (sgn = -1; sgn <= 1; sgn += 2) {
    /* make the body */
    status = EG_makeSolidBody(context, sgn*CONE, data, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Faces from the Body */
    status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the tessellation */
    params[0] =  0.5;
    params[1] =  0.3;
    params[2] = 20.0;
    status = EG_makeTessBody(ebody1, params, &tess1);

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" CONE Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" CONE Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 7; iparam++) data_dot[iparam] = 0;

    for (iparam = 0; iparam < 7; iparam++) {

      /* set the velocity of the original body */
      data_dot[iparam] = 1.0;
      status = EG_makeSolidBody_dot(ebody1, sgn*CONE, data, data_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      data_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make a perturbed body for finite difference */
      data[iparam] += dtime;
      status = EG_makeSolidBody(context, sgn*CONE, data, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      data[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "CONE", 1e-7, 1e-7, 1e-7);
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
/*  Cylinder                                                                 */
/*                                                                           */
/*****************************************************************************/

int
pingCylinder(ego context)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, sgn;
  double data[7], data_dot[7], params[3], dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  data[0] = 5; /* vertex */
  data[1] = 4;
  data[2] = 6;
  data[3] = 1; /* base */
  data[4] = 2;
  data[5] = 3;
  data[6] = 4; /* radius */

  for (sgn = -1; sgn <= 1; sgn += 2) {
    /* make the body */
    status = EG_makeSolidBody(context, sgn*CYLINDER, data, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Faces from the Body */
    status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the tessellation */
    params[0] =  0.4;
    params[1] =  0.2;
    params[2] = 20.0;
    status = EG_makeTessBody(ebody1, params, &tess1);

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" CYLINDER Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" CYLINDER Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 7; iparam++) data_dot[iparam] = 0;

    for (iparam = 0; iparam < 7; iparam++) {

      /* set the velocity of the original body */
      data_dot[iparam] = 1.0;
      status = EG_makeSolidBody_dot(ebody1, sgn*CYLINDER, data, data_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      data_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make a perturbed body for finite difference */
      data[iparam] += dtime;
      status = EG_makeSolidBody(context, sgn*CYLINDER, data, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      data[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "CYLINDER", 1e-7, 1e-7, 1e-7);
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
/*  Torus                                                                    */
/*                                                                           */
/*****************************************************************************/

int
pingTorus(ego context)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, sgn;
  double data[8], data_dot[8], params[3], dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  data[0] = 5; /* center */
  data[1] = 4;
  data[2] = 6;
  data[3] = 1; /* axis */
  data[4] = 2;
  data[5] = 3;
  data[6] = 4;   /* major radius */
  data[7] = 0.5; /* minor radius */

  /* sgn =1 because EG_mapTessBody does not work for a torus with 1 node */
  for (sgn = 1; sgn <= 1; sgn += 2) {

    /* make the body */
    status = EG_makeSolidBody(context, sgn*TORUS, data, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Faces from the Body */
    status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make the tessellation */
    params[0] =  0.5;
    params[1] =  0.3;
    params[2] = 20.0;
    status = EG_makeTessBody(ebody1, params, &tess1);

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" TORUS Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" TORUS Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 8; iparam++) data_dot[iparam] = 0;

    for (iparam = 0; iparam < 8; iparam++) {

      /* set the velocity of the original body */
      data_dot[iparam] = 1.0;
      status = EG_makeSolidBody_dot(ebody1, sgn*TORUS, data, data_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      data_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make a perturbed body for finite difference */
      data[iparam] += dtime;
      status = EG_makeSolidBody(context, sgn*TORUS, data, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      data[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "TORUS", 1e-7, 1e-7, 1e-7);
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
/*  makeFace                                                                 */
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
  double data[6], data_dot[6], x1[3], x1_dot[3], x2[3], x2_dot[3];
  double tdata[2], tdata_dot[2];
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

  /* set the Edge t-range sensitivity */
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
makePlaneLoop( ego context,         /* (in)  EGADS context    */
               objStack *stack,     /* (in)  EGADS obj stack  */
               const double *x0,    /* (in)  Node 0 coord     */
               const double *x1,    /* (in)  Node 1 coord     */
               const double *x2,    /* (in)  Node 2 coord     */
               ego *eloop )         /* (out) Loop             */
{
  int    status = EGADS_SUCCESS;
  int    senses[3] = {SFORWARD, SFORWARD, SFORWARD};
  double data[3];
  ego    eedges[3], enodes[3];

  /* create the Nodes for the Edges */
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

  data[0] = x2[0];
  data[1] = x2[1];
  data[2] = x2[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Edges */
  status =  makeLineEdge(context, stack, enodes[0], enodes[1], &eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[1], enodes[2], &eedges[1] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status =  makeLineEdge(context, stack, enodes[2], enodes[0], &eedges[2] );
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loopy */
  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 3, eedges, senses, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:

  return status;
}


int
setPlaneLoop_dot( const double *x0,     /* (in)  Node 0 coord          */
                  const double *x0_dot, /* (in)  Node 0 coord velocity */
                  const double *x1,     /* (in)  Node 1 coord          */
                  const double *x1_dot, /* (in)  Node 1 coord velocity */
                  const double *x2,     /* (in)  Node 2 coord          */
                  const double *x2_dot, /* (in)  Node 2 coord velocity */
                  ego eloop )           /* (in/out) Loop with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, oclass, mtype, *lsens, *senses;
  double data[3], data_dot[3];
  ego    *eedges, enodes[3]={NULL,NULL,NULL}, *enode, eref;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloop, &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &lsens);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Edges */
  status = EG_getTopology(eedges[0], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = enode[0];
  enodes[1] = enode[1];

  status = EG_getTopology(eedges[1], &eref, &oclass, &mtype,
                          data, &nnode, &enode, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[1] = enode[0];
  enodes[2] = enode[1];

  /* create the Nodes for the Edges */
  data[0]     = x0[0];
  data[1]     = x0[1];
  data[2]     = x0[2];
  data_dot[0] = x0_dot[0];
  data_dot[1] = x0_dot[1];
  data_dot[2] = x0_dot[2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = x1[0];
  data[1]     = x1[1];
  data[2]     = x1[2];
  data_dot[0] = x1_dot[0];
  data_dot[1] = x1_dot[1];
  data_dot[2] = x1_dot[2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = x2[0];
  data[1]     = x2[1];
  data[2]     = x2[2];
  data_dot[0] = x2_dot[0];
  data_dot[1] = x2_dot[1];
  data_dot[2] = x2_dot[2];
  status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Line Edge velocities */
  status = setLineEdge_dot( eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[1] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = setLineEdge_dot( eedges[2] );
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:

  return status;
}


int
pingMakeFace(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[9], x_dot[9], params[3], dtime = 1e-8;
  double *x0, *x1, *x2, *x0_dot, *x1_dot, *x2_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *xyz1, *uv1;
  ego    eloop1, eloop2, eface1, eface2, ebody, ebody1, ebody2, tess1, tess2;

  x0 = x;
  x1 = x+3;
  x2 = x+6;

  x0_dot = x_dot;
  x1_dot = x_dot+3;
  x2_dot = x_dot+6;

  /* make the Plane Face body */
  x0[0] = 0.00; x0[1] = 0.00; x0[2] = 0.00;
  x1[0] = 1.10; x1[1] = 0.10; x1[2] = 0.05;
  x2[0] = 0.05; x2[1] = 1.20; x2[2] = 0.10;
  status = makePlaneLoop(context, stack, x0, x1, x2, &eloop1);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeFace(eloop1, SFORWARD, NULL, &eface1);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                           NULL, 1, &eface1, NULL, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Faces from the Body */
  status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.5;
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);
  if (status != EGADS_SUCCESS) goto cleanup;


  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &xyz1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" makeFace Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &xyz1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" makeFace Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 9; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 9; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setPlaneLoop_dot(x0, x0_dot,
                              x1, x1_dot,
                              x2, x2_dot, eloop1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_makeFace_dot(eface1, eloop1, NULL, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                             NULL, 1, &eface1, NULL, &ebody);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_copyGeometry_dot(ebody, NULL, NULL, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    EG_deleteObject(ebody);

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makePlaneLoop(context, stack, x0, x1, x2, &eloop2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    status = EG_makeFace(eloop2, SFORWARD, NULL, &eface2);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                             NULL, 1, &eface2, NULL, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "makeFace", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);
    EG_deleteObject(eface2);
  }

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);
  EG_deleteObject(eface1);

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

  status = pingBox(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;
#if 0
  status = pingSphere(context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingCone(context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingCylinder(context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingTorus(context);
  if (status != EGADS_SUCCESS) goto cleanup;
#endif
  status = pingMakeFace(context, &stack);
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
