#include <math.h>
#include "egads.h"
#include "egads_dot.h"
#include "egadsSplineVels.h"

#include "../src/egadsStack.h"

int
EG_isoCurve(const int *header2d, const double *data2d,
            const int ik, const int jk, int *header, double **data);

int
EG_isoCurve_dot(const int *header2d, const double *data2d,
                const double *data2d_dot, const int ik, const int jk,
                int *header, double **data, double **data_dot);


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

    /* note that a negative index means the .fMap Attribute is taken care of in
       EG_getTessFace() */
    status = EG_getTessFace(tess2,-iface-1, &np2, &x2, &uv2, &pt2, &pi2,
                                            &nt2, &ts2, &tc2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* check dot information */
    status = EG_hasGeometry_dot(efaces1[iface]);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (n = 0; n < np1; n++) {

      /* evaluate original face and velocities*/
      status = EG_evaluate_dot(efaces1[iface], &uv1[2*n], NULL, p1, p1_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* evaluate perturbed face */
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
    status = EG_getTessEdge(tess2,-iedge-1, &np2, &x2, &t2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* check dot information */
    status = EG_hasGeometry_dot(eedges1[iedge]);
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

    //printf("range_dot = (%+f, %+f)\n", range_dot[0], range_dot[1]);
    //printf("   fd_dot = (%+f, %+f)\n", fd_dot[0], fd_dot[1]);
    //printf("\n");
  }

  for (inode = 0; inode < nnode; inode++) {

    /* check dot information */
    status = EG_hasGeometry_dot(enodes1[inode]);
    if (status != EGADS_SUCCESS) goto cleanup;

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
/*  equivDotVels                                                             */
/*                                                                           */
/*****************************************************************************/

int velocityOfRange( void* usrData,        /* (in)  blind pointer to user data */
                     const ego *sections,  /* (in)  array of sections */
                     int isec,             /* (in)  current section (bias-0) */
                     ego eedge,            /* (in)  the Edge */
                     double *trange,       /* (out) t-range of the eedge */
                     double *trange_dot )  /* (out) t-range velocity of eedge */
{
  int status = EGADS_SUCCESS;
  int periodic;

  /* get the velocity of the edge t-range */
  status = EG_getRange_dot(eedge, trange, trange_dot, &periodic);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

int
velocityOfNode(void*     usrData,   /* (in)  blind pointer to user data */
               const ego secs[],    /* (in)  array of sections */
               int       isec,      /* (in)  current section (bias-0) */
               ego       enode,     /* (in)  Node for which sensitivity is desired */
               ego       eedge,     /* (in)  Edge attached to the node */
               double    xyz[],     /* (out) coordinates of enode */
               double    xyz_dot[]) /* (out) velocity    of enode */
{
  int status = EGADS_SUCCESS;

  /* evaluate the point and sensitivity */
  status = EG_evaluate_dot(enode, NULL, NULL, xyz, xyz_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

int
velocityOfEdge(void*     usrData,        /* (in)  blind pointer to user data */
               const ego secs[],         /* (in)  array of sections */
               int       isec,           /* (in)  currnt section (bias-0) */
               ego       eedge,          /* (in)  Edge for which sensitivity is desired */
         const int       npnt,           /* (in)  number   of t for evaluation */
         const double    ts[],           /* (in)  array    of t for evaluation */
         const double    ts_dot[],       /* (in)  velocity of t for evaluation */
               double    xyz[],          /* (out) coordinates at ts */
               double    xyz_dot[],      /* (out) velocity    at ts */
               double    dxdt_beg[],     /* (out) tangent vector at beg of eedge */
               double    dxdt_beg_dot[], /* (out) velocity of dxdt_beg */
               double    dxdt_end[],     /* (out) tangent vector ay end of eedge */
               double    dxdt_end_dot[]) /* (out) velocity of dxdt_end */
{
  int status = EGADS_SUCCESS;
  int ipnt;
  double x[18], x_dot[18];

  for (ipnt = 0; ipnt < npnt; ipnt++) {
    /* evaluate the points and sensitivity */
    status = EG_evaluate_dot(eedge, &ts[ipnt], &ts_dot[ipnt], x, x_dot);
    if (status != EGADS_SUCCESS) goto cleanup;

    xyz[3*ipnt+0] = x[0];
    xyz[3*ipnt+1] = x[1];
    xyz[3*ipnt+2] = x[2];

    xyz_dot[3*ipnt+0] = x_dot[0];
    xyz_dot[3*ipnt+1] = x_dot[1];
    xyz_dot[3*ipnt+2] = x_dot[2];

    /* set the sensitivity of the tangent at the beginning and end */
    if (ipnt == 0) {
      dxdt_beg[0] = x[3];
      dxdt_beg[1] = x[4];
      dxdt_beg[2] = x[5];

      dxdt_beg_dot[0] = x_dot[3];
      dxdt_beg_dot[1] = x_dot[4];
      dxdt_beg_dot[2] = x_dot[5];
    }

    if (ipnt == npnt-1) {
      dxdt_end[0] = x[3];
      dxdt_end[1] = x[4];
      dxdt_end[2] = x[5];

      dxdt_end_dot[0] = x_dot[3];
      dxdt_end_dot[1] = x_dot[4];
      dxdt_end_dot[2] = x_dot[5];
    }
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


int
velocityOfBspline( void* usrData,        /* (in)  blind pointer to user data */
                   const ego *secs,      /* (in)  array of sections */
                   int isec,             /* (in)  current section (bias-0) */
                   ego eedge,            /* (in)  Edge for which sensitivity is desired */
                   ego egeom,            /* (in)  B-spline CURVE for the Edge */
                   int **ivec,           /* (out) integer data for the B-spline */
                   double **rvec,        /* (out) real    data for the B-spline */
                   double **rvec_dot )   /* (out) velocity of B-spline real data */
{
  int status = EGADS_SUCCESS;
  int oclass, mtype;
  ego ref;

  /* get the b-spline data sensitivity */
  status = EG_getGeometry(egeom, &oclass, &mtype, &ref, ivec, rvec);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(*rvec);

  /* get the b-spline data sensitivity */
  status = EG_getGeometry_dot(egeom, rvec, rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


int
equivDotVels(ego tess1, ego tess2, int iparam, const char *shape, double ftol, double etol, double ntol)
{
  int    status = EGADS_SUCCESS;
  int    n, d, np1, np2, nt1, nt2, nerr=0;
  int    nface, nedge, nnode, iface, iedge, inode, oclass, mtype, periodic;
  double p1[18], p1_dot[18], p2[18], p2_dot[18], range1[2], range1_dot[2], range2[2], range2_dot[2];
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

      /* evaluate face velocities with _dot functions */
      status = EG_evaluate_dot(efaces1[iface], &uv1[2*n], NULL, p1, p1_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* evaluate face velocities with _vels functions */
      status = EG_evaluate_dot(efaces2[iface], &uv2[2*n], NULL, p2, p2_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      for (d = 0; d < 3; d++) {
        if (fabs(p1_dot[d] - p2_dot[d]) > ftol) {
          printf("%s Face %d iparam=%d, p1[%d]=%+le diff fabs(%+le - %+le) = %+le > %e\n",
                 shape, iface+1, iparam, d, p1[d], p1_dot[d], p2_dot[d], fabs(p1_dot[d] - p2_dot[d]), ftol);
          nerr++;
        }
      }

      //printf("p1_dot = (%+f, %+f, %+f)\n", p1_dot[0], p1_dot[1], p1_dot[2]);
      //printf("p2_dot = (%+f, %+f, %+f)\n", p2_dot[0], p2_dot[1], p2_dot[2]);
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

    /* check dot information */
    status = EG_hasGeometry_dot(eedges1[iedge]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_hasGeometry_dot(eedges2[iedge]);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (n = 0; n < np1; n++) {

      /* evaluate edge velocities with _dot functions */
      status = EG_evaluate_dot(eedges1[iedge], &t1[n], NULL, p1, p1_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* evaluate face velocities with _vels functions */
      status = EG_evaluate_dot(eedges2[iedge], &t2[n], NULL, p2, p2_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      for (d = 0; d < 3; d++) {
        if (fabs(p1_dot[d] - p2_dot[d]) > etol) {
          printf("%s Edge %d iparam=%d, p1[%d]=%+le diff fabs(%+le - %+le) = %+le > %e\n",
                 shape, iedge+1, iparam, d, p1[d], p1_dot[d], p2_dot[d], fabs(p1_dot[d] - p2_dot[d]), etol);
          nerr++;
        }
      }

      //printf("p1_dot = (%+f, %+f, %+f)\n", p1_dot[0], p1_dot[1], p1_dot[2]);
      //printf("p2_dot = (%+f, %+f, %+f)\n", p2_dot[0], p2_dot[1], p2_dot[2]);
      //printf("\n");
    }

    /* check t-range sensitivity */
    status = EG_getRange_dot( eedges1[iedge], range1, range1_dot, &periodic );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_getRange_dot( eedges2[iedge], range2, range2_dot, &periodic );
    if (status != EGADS_SUCCESS) goto cleanup;

    for (d = 0; d < 2; d++) {
      if (fabs(range1_dot[d] - range2_dot[d]) > etol) {
        printf("%s Edge %d iparam=%d, trng[%d]=%+le fabs(%+le - %+le) = %+le > %e\n",
               shape, iedge+1, iparam, d, range1[d], range1_dot[d], range2_dot[d], fabs(range1_dot[d] - range2_dot[d]), etol);
        nerr++;
      }
    }

    //printf("range_dot = (%+f, %+f)\n", range_dot[0], range_dot[1]);
    //printf("   fd_dot = (%+f, %+f)\n", fd_dot[0], fd_dot[1]);
    //printf("\n");
  }

  for (inode = 0; inode < nnode; inode++) {

    /* check dot information */
    status = EG_hasGeometry_dot(enodes1[inode]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_hasGeometry_dot(enodes2[inode]);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* evaluate node velocities with _dot functions */
    status = EG_evaluate_dot(enodes1[inode], NULL, NULL, p1, p1_dot);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* evaluate node velocities with _vels functions */
    status = EG_evaluate_dot(enodes2[inode], NULL, NULL, p2, p2_dot);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (d = 0; d < 3; d++) {
      if (fabs(p1_dot[d] - p2_dot[d]) > etol) {
        printf("%s Node %d iparam=%d, p1[%d]=%+le diff fabs(%+le - %+le) = %+le > %e\n",
               shape, inode+1, iparam, d, p1[d], p1_dot[d], p2_dot[d], fabs(p1_dot[d] - p2_dot[d]), etol);
        nerr++;
      }
    }

    //printf("p1_dot = (%+f, %+f, %+f)\n", p1_dot[0], p1_dot[1], p1_dot[2]);
    //printf("p2_dot = (%+f, %+f, %+f)\n", p2_dot[0], p2_dot[1], p2_dot[2]);
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
/*  Transform                                                                */
/*                                                                           */
/*****************************************************************************/

int
makeTransform( objStack *stack, /* (in)  EGADS obj stack             */
               ego eobj,        /* (in) the object to be transformed */
               double *xforms,  /* (in) sequence of transformations  */
               ego *result )    /* (out) transformed object          */
{
  int    status = EGADS_SUCCESS;
  double scale, offset[3];
  ego    context, exform;

  double mat[12] = {1.00, 0.00, 0.00, 0.0,
                    0.00, 1.00, 0.00, 0.0,
                    0.00, 0.00, 1.00, 0.0};

  status = EG_getContext(eobj, &context);
  if (status != EGADS_SUCCESS) goto cleanup;

  scale     = xforms[0];
  offset[0] = xforms[1];
  offset[1] = xforms[2];
  offset[2] = xforms[3];

  mat[ 0] = scale; mat[ 1] =    0.; mat[ 2] =    0.; mat[ 3] = offset[0];
  mat[ 4] =    0.; mat[ 5] = scale; mat[ 6] =    0.; mat[ 7] = offset[1];
  mat[ 8] =    0.; mat[ 9] =    0.; mat[10] = scale; mat[11] = offset[2];

  status = EG_makeTransform(context, mat, &exform); if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_copyObject(eobj, exform, result);     if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_deleteObject(exform);                 if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_stackPush(stack, *result);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
setTransform_dot(ego eobj,           /* (in) the object to be transformed           */
                 double *xforms,     /* (in) sequence of transformations            */
                 double *xforms_dot, /* (in) velocity of transformations            */
                 ego result)         /* (in/out) transformed object with velocities */
{
  int status = EGADS_SUCCESS;
  double scale, scale_dot, offset[3], offset_dot[3];

  double mat[12] = {1.00, 0.00, 0.00, 0.0,
                    0.00, 1.00, 0.00, 0.0,
                    0.00, 0.00, 1.00, 0.0};
  double mat_dot[12];

  scale     = xforms[0];
  offset[0] = xforms[1];
  offset[1] = xforms[2];
  offset[2] = xforms[3];

  scale_dot     = xforms_dot[0];
  offset_dot[0] = xforms_dot[1];
  offset_dot[1] = xforms_dot[2];
  offset_dot[2] = xforms_dot[3];

  mat[ 0] = scale; mat[ 1] =    0.; mat[ 2] =    0.; mat[ 3] = offset[0];
  mat[ 4] =    0.; mat[ 5] = scale; mat[ 6] =    0.; mat[ 7] = offset[1];
  mat[ 8] =    0.; mat[ 9] =    0.; mat[10] = scale; mat[11] = offset[2];

  mat_dot[ 0] = scale_dot; mat_dot[ 1] =        0.; mat_dot[ 2] =        0.; mat_dot[ 3] = offset_dot[0];
  mat_dot[ 4] =        0.; mat_dot[ 5] = scale_dot; mat_dot[ 6] =        0.; mat_dot[ 7] = offset_dot[1];
  mat_dot[ 8] =        0.; mat_dot[ 9] =        0.; mat_dot[10] = scale_dot; mat_dot[11] = offset_dot[2];

  status = EG_copyGeometry_dot(eobj, mat, mat_dot, result);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Node                                                                     */
/*                                                                           */
/*****************************************************************************/

int
makeNode( ego context,      /* (in)  EGADS context                   */
          objStack *stack,  /* (in)  EGADS obj stack                 */
          const double *x0, /* (in)  coordinates of the point        */
          ego *enode )      /* (out) Node loop created from points   */
{
  int    status = EGADS_SUCCESS;
  double data[3];

  /* create Node */
  data[0] = x0[0];
  data[1] = x0[1];
  data[2] = x0[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, enode);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
setNode_dot( const double *x0,     /* (in)  coordinates of the point     */
             const double *x0_dot, /* (in)  velocity of the point        */
             ego enode )           /* (in/out) Line loop with velocities */
{
  int    status = EGADS_SUCCESS;

  /* set the sensitivity of the Nodes */
  status = EG_setGeometry_dot(enode, NODE, 0, NULL, x0, x0_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingNodeRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, iedge, nedge, nsec = 4;
  double x[3*5], x_dot[3*5], *p1, *p2, *p3, *p4, *p1_dot, *p2_dot, *p3_dot, *p4_dot, params[3], dtime = 1e-7;
  const double *t1, *x1;
  ego    secs1[4], secs2[4], ebody1, ebody2, tess1, tess2;

  p1 = x;
  p2 = x+3;
  p3 = x+6;
  p4 = x+9;

  p1_dot = x_dot;
  p2_dot = x_dot+3;
  p3_dot = x_dot+6;
  p4_dot = x_dot+9;

  printf(" ---------------------------------\n");
  printf(" Ping Ruled Node\n");

  /* make the ruled body */
  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
  p3[0] = 1.00; p3[1] = 1.20; p3[2] = 0.10;
  p4[0] = 1.00; p4[1] = 1.20; p4[2] = 1.10;
  status = makeNode(context, stack, p1, &secs1[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = makeNode(context, stack, p2, &secs1[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = makeNode(context, stack, p3, &secs1[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = makeNode(context, stack, p4, &secs1[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_ruled(nsec, secs1, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* tessellation parameters */
  params[0] =  0.4;
  params[1] =  0.01;
  params[2] = 12.0;

  /* make the tessellation */
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Ping Ruled Node Edge %d np1 = %d\n", iedge+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 12; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 12; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setNode_dot(p1, p1_dot, secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = setNode_dot(p2, p2_dot, secs1[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = setNode_dot(p3, p3_dot, secs1[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = setNode_dot(p4, p4_dot, secs1[3]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled_dot(ebody1, nsec, secs1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeNode(context, stack, p1, &secs2[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p2, &secs2[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p3, &secs2[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p4, &secs2[3]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled(nsec, secs2, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Ping Ruled Node", 1e-7, 1e-7, 1e-7);
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


int
pingNodeBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    j, ci, iparam, np1, iedge, nedge, nsec = 5;
  double x[3*5], x_dot[3*5], *p1, *p2, *p3, *p4, *p5, *p1_dot, *p2_dot, *p3_dot, *p4_dot, *p5_dot, params[3], dtime = 1e-7;
  const double *t1, *x1;
  ego    secs1[7]={NULL}, secs2[7], ebody1, ebody2, tess1, tess2;

  p1 = x;
  p2 = x+3;
  p3 = x+6;
  p4 = x+9;
  p5 = x+12;

  p1_dot = x_dot;
  p2_dot = x_dot+3;
  p3_dot = x_dot+6;
  p4_dot = x_dot+9;
  p5_dot = x_dot+12;

  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
  p3[0] = 1.00; p3[1] = 1.20; p3[2] = 0.10;
  p4[0] = 1.00; p4[1] = 1.20; p4[2] = 1.10;
  p5[0] = 2.00; p5[1] = 1.20; p5[2] = 1.10;

  for (ci = 0; ci < 3; ci++) {
    printf(" ---------------------------------\n");
    printf(" Ping Blend Node C%d\n", 2-ci);

    nsec = 5+ci;

    /* make the blend body */
    status = makeNode(context, stack, p1, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p2, &secs1[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p3, &secs1[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (j = 0; j <= ci; j++) secs1[2+j] = secs1[2];

    status = makeNode(context, stack, p4, &secs1[3+ci]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p5, &secs1[4+ci]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs1, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Blend Node Edge %d np1 = %d\n", iedge+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 3*5; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 3*5; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setNode_dot(p1, p1_dot, secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setNode_dot(p2, p2_dot, secs1[1]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setNode_dot(p3, p3_dot, secs1[2]);
      if (status != EGADS_SUCCESS) goto cleanup;

      for (j = 0; j <= ci; j++) secs1[2+j] = secs1[2];

      status = setNode_dot(p4, p4_dot, secs1[3+ci]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setNode_dot(p5, p5_dot, secs1[4+ci]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs1, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeNode(context, stack, p1, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNode(context, stack, p2, &secs2[1]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNode(context, stack, p3, &secs2[2]);
      if (status != EGADS_SUCCESS) goto cleanup;

      for (j = 0; j <= ci; j++) secs2[2+j] = secs2[2];

      status = makeNode(context, stack, p4, &secs2[3+ci]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNode(context, stack, p5, &secs2[4+ci]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs2, NULL, NULL, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Blend Node", 5e-7, 1e-7, 1e-7);
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
equivNodeRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, iedge, nedge, nsec = 4;
  double x[3*4], x_dot[3*4], *p1, *p2, *p3, *p4, *p1_dot, *p2_dot, *p3_dot, *p4_dot, params[3];
  const double *t1, *x1;
  ego    secs[4], ebody1, ebody2, tess1, tess2;

  egadsSplineVels vels;
  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  p1 = x;
  p2 = x+3;
  p3 = x+6;
  p4 = x+9;

  p1_dot = x_dot;
  p2_dot = x_dot+3;
  p3_dot = x_dot+6;
  p4_dot = x_dot+9;

  printf(" ---------------------------------\n");
  printf(" Equiv Ruled Node\n");

  /* make the body for _dot sensitivities */
  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
  p3[0] = 1.00; p3[1] = 1.20; p3[2] = 0.10;
  p4[0] = 1.00; p4[1] = 1.20; p4[2] = 1.10;
  status = makeNode(context, stack, p1, &secs[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = makeNode(context, stack, p2, &secs[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = makeNode(context, stack, p3, &secs[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = makeNode(context, stack, p4, &secs[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_ruled(nsec, secs, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the body for _vels sensitivities */
  status = EG_ruled(nsec, secs, &ebody2);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* tessellation parameters */
  params[0] =  0.4;
  params[1] =  0.01;
  params[2] = 12.0;

  /* make the tessellation */
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* map the tessellation */
  status = EG_mapTessBody(tess1, ebody2, &tess2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Equiv Ruled Line Edge %d np1 = %d\n", iedge+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 12; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 12; iparam++) {

    /* set the velocity with _dot functions */
    x_dot[iparam] = 1.0;
    status = setNode_dot(p1, p1_dot, secs[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = setNode_dot(p2, p2_dot, secs[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = setNode_dot(p3, p3_dot, secs[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = setNode_dot(p4, p4_dot, secs[3]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled_dot(ebody1, nsec, secs);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* set the velocity with _vels functions */
    status = EG_ruled_vels(nsec, secs, &vels, ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* ping the bodies */
    status = equivDotVels(tess1, tess2, iparam, "Equiv Ruled Node", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

  }
  EG_deleteObject(tess2);
  EG_deleteObject(ebody2);

  EG_deleteObject(tess1);
  EG_deleteObject(ebody1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
equivNodeBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    j, ci, iparam, np1, iedge, nedge, nsec = 5;
  double x[3*5], x_dot[3*5], *p1, *p2, *p3, *p4, *p5, *p1_dot, *p2_dot, *p3_dot, *p4_dot, *p5_dot, params[3];
  const double *t1, *x1;
  ego    secs[7], ebody1, ebody2, tess1, tess2;

  egadsSplineVels vels;
  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  p1 = x;
  p2 = x+3;
  p3 = x+6;
  p4 = x+9;
  p5 = x+12;

  p1_dot = x_dot;
  p2_dot = x_dot+3;
  p3_dot = x_dot+6;
  p4_dot = x_dot+9;
  p5_dot = x_dot+12;

  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
  p3[0] = 1.00; p3[1] = 1.20; p3[2] = 0.10;
  p4[0] = 1.00; p4[1] = 1.20; p4[2] = 1.10;
  p5[0] = 2.00; p5[1] = 1.20; p5[2] = 1.10;

  for (ci = 0; ci < 3; ci++) {
    printf(" ---------------------------------\n");
    printf(" Equiv Blend Node C%d\n", 2-ci);

    nsec = 5+ci;

    /* make the blend body */
    status = makeNode(context, stack, p1, &secs[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p2, &secs[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p3, &secs[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (j = 0; j <= ci; j++) secs[2+j] = secs[2];

    status = makeNode(context, stack, p4, &secs[3+ci]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeNode(context, stack, p5, &secs[4+ci]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make the body for _vels sensitivities */
    status = EG_blend(nsec, secs, NULL, NULL, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the Edges from the Body */
    status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Equiv Blend Node Edge %d np1 = %d\n", iedge+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 3*5; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 3*5; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setNode_dot(p1, p1_dot, secs[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setNode_dot(p2, p2_dot, secs[1]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setNode_dot(p3, p3_dot, secs[2]);
      if (status != EGADS_SUCCESS) goto cleanup;

      for (j = 0; j <= ci; j++) secs[2+j] = secs[2];

      status = setNode_dot(p4, p4_dot, secs[3+ci]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setNode_dot(p5, p5_dot, secs[4+ci]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* set the velocity with _vels functions */
      status = EG_blend_vels(nsec, secs, NULL, NULL, NULL, NULL, &vels, ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* check sensitivity equivalence */
      status = equivDotVels(tess1, tess2, iparam, "Equiv Blend Node", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;
    }
    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);

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
/*  Line                                                                     */
/*                                                                           */
/*****************************************************************************/

int
makeLineLoop( ego context,      /* (in)  EGADS context                   */
              objStack *stack,  /* (in)  EGADS obj stack                 */
              const double *x0, /* (in)  coordinates of the first point  */
              const double *x1, /* (in)  coordinates of the second point */
              ego *eloop )      /* (out) Line loop created from points   */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[6], tdata[2];
  ego    eline, eedge, enodes[2];

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
                           NULL, 1, &eedge, senses, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
setLineLoop_dot( const double *x0,     /* (in)  coordinates of the first point  */
                 const double *x0_dot, /* (in)  velocity of the first point     */
                 const double *x1,     /* (in)  coordinates of the second point */
                 const double *x1_dot, /* (in)  velocity of the second point    */
                 ego eloop )           /* (in/out) Line loop with velocities    */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, oclass, mtype, *senses;
  double data[6], data_dot[6], tdata[2], tdata_dot[2];
  ego    eline, *enodes, *eedges, eref;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloop, &eref, &oclass, &mtype,
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
pingLineRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[6], x_dot[6], *p1, *p2, *p1_dot, *p2_dot, params[3], dtime = 1e-7;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2;

  p1 = x;
  p2 = x+3;

  p1_dot = x_dot;
  p2_dot = x_dot+3;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Ruled Line dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 0.75;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the ruled body */
    p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
    p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
    status = makeLineLoop(context, stack, p1, p2, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled(nsec, secs1, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

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
      printf(" Ping Ruled Line Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Ruled Line Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 6; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 6; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setLineLoop_dot(p1, p1_dot, p2, p2_dot, secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled_dot(ebody1, nsec, secs1);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeLineLoop(context, stack, p1, p2, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled(nsec, secs2, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Ruled Line", 1e-7, 1e-7, 1e-7);
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
pingLineBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[6], x_dot[6], *p1, *p2, *p1_dot, *p2_dot, params[3], dtime = 1e-7;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2;

  p1 = x;
  p2 = x+3;

  p1_dot = x_dot;
  p2_dot = x_dot+3;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Blend Line dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 0.75;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the ruled body */
    p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
    p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
    status = makeLineLoop(context, stack, p1, p2, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs1, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

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
      printf(" Ping Blend Line Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Blend Line Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 6; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 6; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setLineLoop_dot(p1, p1_dot, p2, p2_dot, secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs1, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeLineLoop(context, stack, p1, p2, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs2, NULL, NULL, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Blend Line", 5e-7, 1e-7, 1e-7);
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
equivLineRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[6], x_dot[6], *p1, *p2, *p1_dot, *p2_dot, params[3];
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs[3], ebody1, ebody2, tess1, tess2;

  egadsSplineVels vels;
  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  p1 = x;
  p2 = x+3;

  p1_dot = x_dot;
  p2_dot = x_dot+3;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Equiv Ruled Line dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 0.75;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the body for _dot sensitivities */
    p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
    p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
    status = makeLineLoop(context, stack, p1, p2, &secs[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs[0], xform1, &secs[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs[0], xform2, &secs[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled(nsec, secs, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make the body for _vels sensitivities */
    status = EG_ruled(nsec, secs, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
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
      printf(" Equiv Ruled Line Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Equiv Ruled Line Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 6; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 6; iparam++) {

      /* set the velocity with _dot functions */
      x_dot[iparam] = 1.0;
      status = setLineLoop_dot(p1, p1_dot, p2, p2_dot, secs[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs[0], xform1, xform_dot, secs[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs[0], xform2, xform_dot, secs[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled_dot(ebody1, nsec, secs);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* set the velocity with _vels functions */
      status = EG_ruled_vels(nsec, secs, &vels, ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* ping the bodies */
      status = equivDotVels(tess1, tess2, iparam, "Equiv Ruled Line", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;

    }
    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);

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
equivLineBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[6], x_dot[6], *p1, *p2, *p1_dot, *p2_dot, params[3];
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs[3], ebody1, ebody2, tess1, tess2;

  egadsSplineVels vels;

  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  p1 = x;
  p2 = x+3;

  p1_dot = x_dot;
  p2_dot = x_dot+3;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Equiv Blend Line dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 0.75;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the body for _dot sensitivities */
    p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
    p2[0] = 1.00; p2[1] = 0.20; p2[2] = 0.10;
    status = makeLineLoop(context, stack, p1, p2, &secs[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs[0], xform1, &secs[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs[0], xform2, &secs[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make the body for _vels sensitivities */
    status = EG_blend(nsec, secs, NULL, NULL, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
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
      printf(" Equiv Blend Line Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Equiv Blend Line Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 6; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 6; iparam++) {

      /* set the velocity with _dot functions */
      x_dot[iparam] = 1.0;
      status = setLineLoop_dot(p1, p1_dot, p2, p2_dot, secs[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs[0], xform1, xform_dot, secs[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs[0], xform2, xform_dot, secs[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* set the velocity with _vels functions */
      status = EG_blend_vels(nsec, secs, NULL, NULL, NULL, NULL, &vels, ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* check sensitivity equivalence */
      status = equivDotVels(tess1, tess2, iparam, "Equiv Blend Line", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;

    }
    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);

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
/*  Line2                                                                    */
/*                                                                           */
/*****************************************************************************/

int
makeLine2Loop( ego context,      /* (in)  EGADS context                   */
               objStack *stack,  /* (in)  EGADS obj stack                 */
               const double *x0, /* (in)  the point                       */
               const double *v0, /* (in)  the vector                      */
               double *ts,       /* (in)  t values for end nodes          */
               ego *eloop )      /* (out) Line loop created from points   */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[18];
  ego    eline, eedge, enodes[2];

  /* create the Line (point and direction) */
  data[0] = x0[0];
  data[1] = x0[1];
  data[2] = x0[2];
  data[3] = v0[0];
  data[4] = v0[1];
  data[5] = v0[2];

  status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &eline);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Nodes for the Edge */
  status = EG_evaluate(eline, &ts[0], data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate(eline, &ts[1], data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the Edge on the Line */
  status = EG_makeTopology(context, eline, EDGE, TWONODE,
                           ts, 2, enodes, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, OPEN,
                           NULL, 1, &eedge, senses, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
setLine2Loop_dot( const double *x0,     /* (in)  the point                       */
                  const double *x0_dot, /* (in)  velocity of the point           */
                  const double *v0,     /* (in)  the vector                      */
                  const double *v0_dot, /* (in)  velocity of the vector          */
                  const double *ts,     /* (in)  t values for end nodes          */
                  const double *ts_dot, /* (in)  velocity of t-values            */
                  ego eloop )           /* (in/out) Line loop with velocities    */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, oclass, mtype, *senses;
  double data[18], data_dot[18];
  ego    eline, *enodes, *eedges, eref;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloop, &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Line from the Edge */
  status = EG_getTopology(eedges[0], &eline, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* Compute the Line data and velocity */
  data[0] = x0[0];
  data[1] = x0[1];
  data[2] = x0[2];
  data[3] = v0[0];
  data[4] = v0[1];
  data[5] = v0[2];

  data_dot[0] = x0_dot[0];
  data_dot[1] = x0_dot[1];
  data_dot[2] = x0_dot[2];
  data_dot[3] = v0_dot[0];
  data_dot[4] = v0_dot[1];
  data_dot[5] = v0_dot[2];

  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the sensitivity of the Nodes */
  status = EG_evaluate_dot(eline, &ts[0], &ts_dot[0], data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate_dot(eline, &ts[1], &ts_dot[1], data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  status = EG_setRange_dot(eedges[0], EDGE, ts, ts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingLine2Ruled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[8], x_dot[8], params[3], dtime = 1e-7;
  double *x0, *v0, *x0_dot, *v0_dot, *ts, *ts_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2;

  x0 = x;
  v0 = x+3;
  ts = x+6;

  x0_dot = x_dot;
  v0_dot = x_dot+3;
  ts_dot = x_dot+6;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Ruled Line2 dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 1.0;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the ruled body */
    x0[0] = 0.00; x0[1] = 0.00; x0[2] = 0.00;
    v0[0] = 1.00; v0[1] = 0.00; v0[2] = 0.00;
    ts[0] = -1; ts[1] = 1;
    status = makeLine2Loop(context, stack, x0, v0, ts, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled(nsec, secs1, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

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
      printf(" Ping Ruled Line2 Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Ruled Line2 Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 8; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 8; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setLine2Loop_dot(x0, x0_dot, v0, v0_dot, ts, ts_dot, secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled_dot(ebody1, nsec, secs1);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeLine2Loop(context, stack, x0, v0, ts, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled(nsec, secs2, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Ruled Line2", 1e-7, 1e-7, 1e-7);
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
pingLine2Blend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[8], x_dot[8], params[3], dtime = 1e-7;
  double *x0, *v0, *x0_dot, *v0_dot, *ts, *ts_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2;

  x0 = x;
  v0 = x+3;
  ts = x+6;

  x0_dot = x_dot;
  v0_dot = x_dot+3;
  ts_dot = x_dot+6;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Blend Line2 dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 1.0;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the ruled body */
    x0[0] = 0.00; x0[1] = 0.00; x0[2] = 0.00;
    v0[0] = 1.00; v0[1] = 0.00; v0[2] = 0.00;
    ts[0] = -1; ts[1] = 1;
    status = makeLine2Loop(context, stack, x0, v0, ts, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs1, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.01;
    params[2] = 12.0;

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
      printf(" Ping Blend Line2 Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Blend Line2 Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 8; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 8; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setLine2Loop_dot(x0, x0_dot, v0, v0_dot, ts, ts_dot, secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs1, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeLine2Loop(context, stack, x0, v0, ts, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs2, NULL, NULL, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Ruled Line2", 1e-7, 1e-7, 1e-7);
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
/*  Circle                                                                   */
/*                                                                           */
/*****************************************************************************/

int
makeCircle( ego context,         /* (in)  EGADS context        */
            objStack *stack,     /* (in)  EGADS obj stack      */
            const int btype,     /* (in)  LOOP or FACE         */
            const double *xcent, /* (in)  Center               */
            const double *xax,   /* (in)  x-axis               */
            const double *yax,   /* (in)  y-axis               */
            const double r,      /* (in)  radius               */
            ego *eobj )          /* (out) Circle Face/Loop     */
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

  if (btype == LOOP) {

    *eobj = eloop;

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

    *eobj = eface;
  }

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setCircle_dot( const double *xcent,     /* (in)  Center          */
               const double *xcent_dot, /* (in)  Center velocity */
               const double *xax,       /* (in)  x-axis          */
               const double *xax_dot,   /* (in)  x-axis velocity */
               const double *yax,       /* (in)  y-axis          */
               const double *yax_dot,   /* (in)  y-axis velocity */
               const double r,          /* (in)  radius          */
               const double r_dot,      /* (in)  radius velocity */
               ego eobj )               /* (in/out) Circle with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses, btype;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3];
  double tdata[2], tdata_dot[2];
  double *rvec=NULL, *rvec_dot=NULL;
  ego    ecircle, eplane, *enodes, *eloops, *eedges, eref;

  /* get the type */
  status = EG_getTopology(eobj, &eplane, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (oclass == LOOP) {
    nloop = 1;
    eloops = &eobj;
    btype = LOOP;
  } else {
    btype = FACE;
  }

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
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

  if (btype == FACE) {
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
pingCircleRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[10], x_dot[10], params[3], dtime = 1e-7;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2, eloop1, eloop2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Ruled Circle dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 1.0;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the ruled body */
    xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
    xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
    yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
    x[9] = 1.0;
    status = makeCircle(context, stack, FACE, xcent, xax, yax, x[9], &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled(nsec, secs1, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.2;
    params[2] = 20.0;

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
      printf(" Ping Ruled Circle Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Ruled Circle Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 10; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], eloop1);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled_dot(ebody1, nsec, secs1);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeCircle(context, stack, FACE, xcent, xax, yax, x[9], &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop2);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled(nsec, secs2, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Ruled Circle", 1e-7, 1e-7, 1e-7);
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
pingCircleBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[10], x_dot[10], params[3], dtime = 1e-7;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2, eloop1, eloop2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Blend Circle dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 1.0;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the ruled body */
    xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
    xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
    yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
    x[9] = 1.0;
    status = makeCircle(context, stack, FACE, xcent, xax, yax, x[9], &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs1, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.2;
    params[2] = 20.0;

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
      printf(" Ping Blend Circle Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Blend Circle Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 10; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], secs1[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], eloop1);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs1, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = makeCircle(context, stack, FACE, xcent, xax, yax, x[9], &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop2);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs2, NULL, NULL, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Blend Circle", 1e-7, 1e-7, 1e-7);
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
equivCircleRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[10], x_dot[10], params[3];
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs[3], ebody1, ebody2, tess1, tess2, eloop;

  egadsSplineVels vels;
  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Equiv Ruled Circle dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 0.75;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the body for _dot sensitivities */
    xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
    xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
    yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
    x[9] = 1.0;
    status = makeCircle(context, stack, FACE, xcent, xax, yax, x[9], &secs[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop, xform1, &secs[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs[0], xform2, &secs[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_ruled(nsec, secs, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make the body for _vels sensitivities */
    status = EG_ruled(nsec, secs, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.2;
    params[2] = 20.0;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
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
      printf(" Equiv Ruled Circle Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Equiv Ruled Circle Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 10; iparam++) {

      /* set the velocity with _dot functions */
      x_dot[iparam] = 1.0;
      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], secs[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], eloop);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop, xform1, xform_dot, secs[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs[0], xform2, xform_dot, secs[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled_dot(ebody1, nsec, secs);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* set the velocity with _vels functions */
      status = EG_ruled_vels(nsec, secs, &vels, ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      /* ping the bodies */
      status = equivDotVels(tess1, tess2, iparam, "Equiv Ruled Circle", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;
    }

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);

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
equivCircleBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  double x[10], x_dot[10], params[3];
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs[3], ebody1, ebody2, tess1, tess2, eloop;

  egadsSplineVels vels;
  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Equiv Blend Circle dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 0.75;
    xform2[1] = 0.;
    xform2[2] = 0.;
    xform2[3] = 2.*dir;

    /* make the body for _dot sensitivities */
    xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
    xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
    yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
    x[9] = 1.0;
    status = makeCircle(context, stack, FACE, xcent, xax, yax, x[9], &secs[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop, xform1, &secs[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs[0], xform2, &secs[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs, NULL, NULL, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make the body for _vels sensitivities */
    status = EG_blend(nsec, secs, NULL, NULL, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  0.4;
    params[1] =  0.2;
    params[2] = 20.0;

    /* make the tessellation */
    status = EG_makeTessBody(ebody1, params, &tess1);

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
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
      printf(" Equiv Blend Circle Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Equiv Blend Circle Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;


    for (iparam = 0; iparam < 10; iparam++) {

      /* set the velocity with _dot functions */
      x_dot[iparam] = 1.0;
      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], secs[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], eloop);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop, xform1, xform_dot, secs[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs[0], xform2, xform_dot, secs[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs, NULL, NULL, NULL, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* set the velocity with _vels functions */
      status = EG_blend_vels(nsec, secs, NULL, NULL, NULL, NULL, &vels, ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;


      /* ping the bodies */
      status = equivDotVels(tess1, tess2, iparam, "Equiv Blend Circle", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;
    }

    EG_deleteObject(tess2);
    EG_deleteObject(ebody2);

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
/*  Nose treatment                                                           */
/*                                                                           */
/*****************************************************************************/

int
makeCircle2( ego context,         /* (in)  EGADS context        */
             objStack *stack,     /* (in)  EGADS obj stack      */
             const int btype,     /* (in)  LOOP or FACE         */
             const double *xcent, /* (in)  Center               */
             const double *xax,   /* (in)  x-axis               */
             const double *yax,   /* (in)  y-axis               */
             const double r,      /* (in)  radius               */
             ego *eobj )          /* (out) Circle Face/Loop     */
{
  int    status = EGADS_SUCCESS;
  int    senses[2] = {SFORWARD,SFORWARD}, oclass, mtype, *ivec=NULL;
  double data[10], tdata[2], dx[3], dy[3], *rvec=NULL;
  ego    ecircle, eedges[2], enodes[3], eloop, eplane, eface, eref;

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

  /* create the Nodes for the Edge */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] - dx[0]*r;
  data[1] = xcent[1] - dx[1]*r;
  data[2] = xcent[2] - dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  enodes[2] = enodes[0];

  /* make the Edges on the Circle */
  tdata[0] = 0;
  tdata[1] = PI;

  status = EG_makeTopology(context, ecircle, EDGE, TWONODE,
                           tdata, 2, &enodes[0], NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = PI;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle, EDGE, TWONODE,
                           tdata, 2, &enodes[1], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;


  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 2, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (btype == LOOP) {

    *eobj = eloop;

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

    *eobj = eface;
  }

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setCircle2_dot( const double *xcent,     /* (in)  Center          */
                const double *xcent_dot, /* (in)  Center velocity */
                const double *xax,       /* (in)  x-axis          */
                const double *xax_dot,   /* (in)  x-axis velocity */
                const double *yax,       /* (in)  y-axis          */
                const double *yax_dot,   /* (in)  y-axis velocity */
                const double r,          /* (in)  radius          */
                const double r_dot,      /* (in)  radius velocity */
                ego eobj )               /* (in/out) Circle with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses, btype;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3];
  double tdata[2], tdata_dot[2];
  double *rvec=NULL, *rvec_dot=NULL;
  ego    ecircle, eplane, *enodes, *eloops, *eedges, eref;

  /* get the type */
  status = EG_getTopology(eobj, &eplane, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (oclass == LOOP) {
    nloop = 1;
    eloops = &eobj;
    btype = LOOP;
  } else {
    btype = FACE;
  }

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = PI;

  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);

  tdata[0] = PI;
  tdata[1] = TWOPI;

  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);

  /* get the Nodes and the Circle from the Edge */
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

  data[0] = xcent[0] - dx[0]*r;
  data[1] = xcent[1] - dx[1]*r;
  data[2] = xcent[2] - dx[2]*r;
  data_dot[0] = xcent_dot[0] - dx_dot[0]*r - dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] - dx_dot[1]*r - dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] - dx_dot[2]*r - dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (btype == FACE) {
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
pingNoseCircleBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 5;
  double x[26], x_dot[26], params[3], dtime = 1e-7;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot, *rc1, *rc1_dot, *rcN, *rcN_dot;
  double xform1[4], xform2[4], xform3[4], xform4[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[5], secs2[5], ebody1, ebody2, tess1, tess2, eloop1, eloop2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;
  rc1   = x+10;
  rcN   = x+18;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;
  rc1_dot   = x_dot+10;
  rcN_dot   = x_dot+18;

  rc1[0] = 0.2; /* radius of curvature */
  rc1[1] = 1.;  /* first axis */
  rc1[2] = 0.;
  rc1[3] = 0.;

  rc1[4] = 0.1; /* radius of curvature */
  rc1[5] = 0.;  /* first axis */
  rc1[6] = 1.;
  rc1[7] = 0.;


  rcN[0] = 0.1; /* radius of curvature */
  rcN[1] = 1.;  /* first axis */
  rcN[2] = 0.;
  rcN[3] = 0.;

  rcN[4] = 0.2; /* radius of curvature */
  rcN[5] = 0.;  /* first axis */
  rcN[6] = 1.;
  rcN[7] = 0.;

  /* circle parameters */
  xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
  xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
  yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
  x[9] = 1.0;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Nose Circle dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 1.0;
    xform2[1] = -0.1;
    xform2[2] = -0.2;
    xform2[3] = 2.*dir;

    xform3[0] = 1.0;
    xform3[1] = 0.;
    xform3[2] = 0.;
    xform3[3] = 3.*dir;

    xform4[0] = 1.0;
    xform4[1] = 0.;
    xform4[2] = 0.;
    xform4[3] = 4.*dir;

    /* make the body */
    status = EG_makeTopology(context, NULL, NODE, 0,
                             xcent, 0, NULL, NULL, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform3, &secs1[3] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform4, &secs1[4] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs1, rc1, rcN, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  1.0;
    params[1] =  1.0;
    params[2] = 30.0;

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
      printf(" Ping Nose Circle Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Nose Circle Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 26; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 26; iparam++) {
      if (iparam == 10+1 ||
          iparam == 10+2 ||
          iparam == 10+3 ||
          iparam == 10+5 ||
          iparam == 10+6 ||
          iparam == 10+7 ||

          iparam == 18+1 ||
          iparam == 18+2 ||
          iparam == 18+3 ||
          iparam == 18+5 ||
          iparam == 18+6 ||
          iparam == 18+7 ) continue; /* skip axis as they must be orthogonal */

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = EG_setGeometry_dot(secs1[0], NODE, 0, NULL, xcent, xcent_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setCircle_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], eloop1);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform3, xform_dot, secs1[3] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform4, xform_dot, secs1[4] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs1, rc1, rc1_dot, rcN, rcN_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = EG_makeTopology(context, NULL, NODE, 0,
                               xcent, 0, NULL, NULL, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;
      status = EG_stackPush(stack, secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeCircle(context, stack, LOOP, xcent, xax, yax, x[9], &eloop2);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform3, &secs2[3] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform4, &secs2[4] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs2, rc1, rcN, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Nose Circle", 5e-7, 5e-7, 1e-7);
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
pingNoseCircle2Blend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 5;
  double x[26], x_dot[26], params[3], dtime = 1e-7;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot, *rc1, *rc1_dot, *rcN, *rcN_dot;
  double xform1[4], xform2[4], xform3[4], xform4[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[5], secs2[5], ebody1, ebody2, tess1, tess2, eloop1, eloop2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;
  rc1   = x+10;
  rcN   = x+18;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;
  rc1_dot   = x_dot+10;
  rcN_dot   = x_dot+18;

  rc1[0] = 0.2; /* radius of curvature */
  rc1[1] = 1.;  /* first axis */
  rc1[2] = 0.;
  rc1[3] = 0.;

  rc1[4] = 0.1; /* radius of curvature */
  rc1[5] = 0.;  /* first axis */
  rc1[6] = 1.;
  rc1[7] = 0.;


  rcN[0] = 0.1; /* radius of curvature */
  rcN[1] = 1.;  /* first axis */
  rcN[2] = 0.;
  rcN[3] = 0.;

  rcN[4] = 0.2; /* radius of curvature */
  rcN[5] = 0.;  /* first axis */
  rcN[6] = 1.;
  rcN[7] = 0.;

  /* circle parameters */
  xcent[0] = 0.0; xcent[1] = 0.0; xcent[2] = 0.0;
  xax[0]   = 1.0; xax[1]   = 0.0; xax[2]   = 0.0;
  yax[0]   = 0.0; yax[1]   = 1.0; yax[2]   = 0.0;
  x[9] = 1.0;

  for (dir = -1; dir <= 1; dir += 2) {

    printf(" ---------------------------------\n");
    printf(" Ping Nose Circle2 dir %+d\n", dir);

    xform1[0] = 1.0;
    xform1[1] = 0.1;
    xform1[2] = 0.2;
    xform1[3] = 1.*dir;

    xform2[0] = 1.0;
    xform2[1] = -0.1;
    xform2[2] = -0.2;
    xform2[3] = 2.*dir;

    xform3[0] = 1.0;
    xform3[1] = 0.;
    xform3[2] = 0.;
    xform3[3] = 3.*dir;

    xform4[0] = 1.0;
    xform4[1] = 0.;
    xform4[2] = 0.;
    xform4[3] = 4.*dir;

    /* make the body */
    status = EG_makeTopology(context, NULL, NODE, 0,
                             xcent, 0, NULL, NULL, &secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, secs1[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeCircle2(context, stack, LOOP, xcent, xax, yax, x[9], &eloop1);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform1, &secs1[1] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform2, &secs1[2] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, eloop1, xform3, &secs1[3] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = makeTransform(stack, secs1[0], xform4, &secs1[4] );
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_blend(nsec, secs1, rc1, rcN, &ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* test re-making the topology */
    status = remakeTopology(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* tessellation parameters */
    params[0] =  1.0;
    params[1] =  1.0;
    params[2] = 30.0;

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
      printf(" Ping Nose Circle2 Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Ping Nose Circle2 Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 26; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 26; iparam++) {
      if (iparam == 10+1 ||
          iparam == 10+2 ||
          iparam == 10+3 ||
          iparam == 10+5 ||
          iparam == 10+6 ||
          iparam == 10+7 ||

          iparam == 18+1 ||
          iparam == 18+2 ||
          iparam == 18+3 ||
          iparam == 18+5 ||
          iparam == 18+6 ||
          iparam == 18+7 ) continue; /* skip axis as they must be orthogonal */

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = EG_setGeometry_dot(secs1[0], NODE, 0, NULL, xcent, xcent_dot);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setCircle2_dot(xcent, xcent_dot,
                              xax, xax_dot,
                              yax, yax_dot,
                              x[9], x_dot[9], eloop1);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform1, xform_dot, secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform2, xform_dot, secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( eloop1, xform3, xform_dot, secs1[3] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = setTransform_dot( secs1[0], xform4, xform_dot, secs1[4] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend_dot(ebody1, nsec, secs1, rc1, rc1_dot, rcN, rcN_dot);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;


      /* make a perturbed body for finite difference */
      x[iparam] += dtime;
      status = EG_makeTopology(context, NULL, NODE, 0,
                               xcent, 0, NULL, NULL, &secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;
      status = EG_stackPush(stack, secs2[0]);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeCircle2(context, stack, LOOP, xcent, xax, yax, x[9], &eloop2);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform1, &secs2[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform2, &secs2[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop2, xform3, &secs2[3] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs2[0], xform4, &secs2[4] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs2, rc1, rcN, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Ping Nose Circle2", 5e-7, 5e-7, 1e-7);
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
/*  NACA Airfoil                                                             */
/*                                                                           */
/*****************************************************************************/

#define NUMPNTS    101
#define DXYTOL     1.0e-8

/* set KNOTS to 0 for arc-length knots, and -1 for equally spaced knots
 */
#define KNOTS      0

int makeNaca( ego context,     /* (in) EGADS context                         */
              objStack *stack, /* (in)  EGADS obj stack                      */
              int btype,       /* (in) result type (LOOP, FACE, or FACEBODY) */
              int sharpte,     /* (in) sharp or blunt TE                     */
              double m,        /* (in) camber                                */
              double p,        /* (in) maxloc                                */
              double t,        /* (in) thickness                             */
              ego *eobj )      /* (out) result pointer                       */
{
  int status = EGADS_SUCCESS;

  int     ipnt, *header=NULL, sizes[2], sense[3], nedge, oclass, mtype, nPos;
  double  *pnts = NULL, *rdata = NULL;
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
      yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
    } else {
      yt   = t/0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
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
  status = EG_stackPush(stack, ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (btype == CURVE) {
    /* return the loop */
    *eobj = ecurve;
    goto cleanup;
  }

  /* get the B-spline CURVE data */
  status = EG_getGeometry(ecurve, &oclass, &mtype, &eref, &header, &rdata);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Node at trailing edge */
  ipnt = 0;
  data[0] = pnts[3*ipnt  ];
  data[1] = pnts[3*ipnt+1];
  data[2] = pnts[3*ipnt+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* node at leading edge as a function of the spline */
  ipnt = (NUMPNTS - 1) / 2 + 3; /* leading edge index, with knot offset of 3 (cubic)*/
  tle = rdata[ipnt];            /* leading edge t-value (should be very close to (0,0,0) */

  status = EG_evaluate(ecurve, &tle, data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
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
    status = EG_stackPush(stack, enodes[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    enodes[3] = enodes[0];
  } else {
    enodes[2] = enodes[0];
  }

  /* make Edge for lower surface */
  tdata[0] = 0;   /* t value at lower TE because EG_spline1dFit was used */
  tdata[1] = tle;

  /* construct the upper Edge */
  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                           tdata, 2, &enodes[0], NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make Edge for upper surface */
  tdata[0] = tdata[1]; /* t-value at leading edge */
  tdata[1] = 1;        /* t value at upper TE because EG_spline1dFit was used */

  /* construct the lower Edge */
  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                           tdata, 2, &enodes[1], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
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
    status = EG_stackPush(stack, eline);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make Edge for this line */
    tdata[0] = 0;
    tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

    status = EG_makeTopology(context, eline, EDGE, TWONODE,
                             tdata, 2, &enodes[2], NULL, &eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make sure there are vertexes on the trailing edge for finite differencing */
    nPos = 5;
    status = EG_attributeAdd(eedges[2], ".nPos", ATTRINT, 1, &nPos, NULL, NULL);
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
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (btype == FACE || btype == FACEBODY) {
    /* create a plane for the loop */
    data[0] = 0.;
    data[1] = 0.;
    data[2] = 0.;
    data[3] = 1.; data[4] = 0.; data[5] = 0.;
    data[6] = 0.; data[7] = 1.; data[8] = 0.;

    status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, data, &eplane);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eplane);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* create the face from the plane and the loop */
    status = EG_makeTopology(context, eplane, FACE, SFORWARD,
                             NULL, 1, &eloop, sense, &eface);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eface);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (btype == FACE) {
      *eobj = eface;
    } else {
      /* make a face body */
      status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                               NULL, 1, &eface, sense, eobj);
      if (status != EGADS_SUCCESS) goto cleanup;
      status = EG_stackPush(stack, *eobj);
      if (status != EGADS_SUCCESS) goto cleanup;
    }
  } else {
    /* return the loop */
    *eobj = eloop;
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(header);
  EG_free(rdata);
  EG_free(pnts);

  return status;
}


int setNaca_dot( int sharpte,     /* (in) sharp or blunt TE */
                 double m,        /* (in) camber */
                 double m_dot,    /* (in) camber sensitivity */
                 double p,        /* (in) maxloc */
                 double p_dot,    /* (in) maxloc sensitivity */
                 double t,        /* (in) thickness */
                 double t_dot,    /* (in) thickness sensitivity */
                 ego eobj )       /* (out) ego with sensitivities */
{
  int status = EGADS_SUCCESS;

  int     ipnt, nedge, oclass, btype=0, mtype, nchild, nloop, *senses, sizes[2];
  double  data[18], data_dot[18], trange[4], tdata[2], tdata_dot[2];
  double  zeta, s, *pnts=NULL, *pnts_dot=NULL;
  double  yt, yt_dot, yc, yc_dot, theta, theta_dot;
  double  ycm, ycm_dot, dycm, dycm_dot, tle, tle_dot;
  double  x, x_dot, y, y_dot, *rvec=NULL, *rvec_dot=NULL;
  ego     enodes[3]={NULL,NULL,NULL}, *eedges, ecurve, eline, *eloops, eplane, eface, eref, *echildren;

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
      yt     = t     / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
      yt_dot = t_dot / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1015)))));
    } else {
      yt     = t     / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
      yt_dot = t_dot / 0.20 * (0.2969 * sqrt(s) + s * (-0.1260 + s * (-0.3516 + s * ( 0.2843 + s * (-0.1036)))));
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
      x = s  - yt * sin(theta);
      y = yc + yt * cos(theta);
      x_dot =        - yt_dot * sin(theta) - theta_dot * yt * cos(theta);
      y_dot = yc_dot + yt_dot * cos(theta) - theta_dot * yt * sin(theta);
    } else if (ipnt == NUMPNTS/2) {
      x = 0;
      y = 0;
      x_dot = 0;
      y_dot = 0;
    } else {
      x = s  + yt * sin(theta);
      y = yc - yt * cos(theta);
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



  if (eobj->oclass == CURVE) {
    ecurve = eobj;

  } else {

    /* get the type */
    status = EG_getTopology(eobj, &eplane, &oclass, &mtype,
                            data, &nloop, &eloops, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (oclass == LOOP) {
      nloop = 1;
      eloops = &eobj;
      btype = LOOP;
    } else if (oclass == FACE) {
      btype = FACE;
    } else {
      btype = FACE;
      eface = eloops[0];
      status = EG_getTopology(eface, &eplane, &oclass, &mtype,
                              data, &nloop, &eloops, &senses);
      if (status != EGADS_SUCCESS) goto cleanup;
    }

    /* get the edges */
    status = EG_getTopology(eloops[0], &eref, &oclass, &mtype, data, &nedge, &eedges,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* get the nodes and the curve from the first edge */
    status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype, trange, &nchild, &echildren,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    enodes[0] = echildren[0]; // upper trailing edge
    enodes[1] = echildren[1]; // leading edge
  }

  /* populate spline curve sensitivities */
  sizes[0] = NUMPNTS;
  sizes[1] = KNOTS;
  status = EG_approximate_dot(ecurve, 0, DXYTOL, sizes, pnts, pnts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (eobj->oclass == CURVE) goto cleanup;

  /* set the sensitivity of the Node at trailing edge */
  ipnt = 0;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, &pnts[3*ipnt], &pnts_dot[3*ipnt]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the t-sensitivity at the leading edge */
  status = EG_getGeometry_dot(ecurve, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  ipnt = (NUMPNTS - 1) / 2 + 3; /* leading edge index, with knot offset of 3 (cubic)*/
  tle     = rvec[ipnt];         /* leading edge t-value (should be very close to (0,0,0) */
  tle_dot = rvec_dot[ipnt];


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

  /* set the sensitivity of the Node at leading edge */
  status = EG_evaluate_dot(ecurve, &tle, &tle_dot, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  if (sharpte == 0) {
    /* get trailing edge line and the lower trailing edge node from the 3rd edge */
    status = EG_getTopology(eedges[2], &eline, &oclass, &mtype, data, &nchild, &echildren,
                            &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    enodes[2] = echildren[0]; // lower trailing edge

    /* set the sensitivity of the Node at lower trailing edge */
    ipnt = NUMPNTS - 1;
    status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, &pnts[3*ipnt], &pnts_dot[3*ipnt]);
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

  if (btype == FACE) {
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
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(rvec);
  EG_free(rvec_dot);
  EG_free(pnts);
  EG_free(pnts_dot);

  return status;
}


int
pingNaca(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  int    sharpte = 0;
  double x[3], x_dot[3], params[3], dtime = 1e-8;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;
  enum naca { im, ip, it };

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  /* make the Naca body */
  status = makeNaca( context, stack, FACEBODY, sharpte, x[im], x[ip], x[it], &ebody1 );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* tessellation parameters */
  params[0] =  0.05;
  params[1] =  0.01;
  params[2] = 15.0;

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
    printf(" Ping NACA Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Ping NACA Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 3; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 3; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setNaca_dot( sharpte,
                          x[im], x_dot[im],
                          x[ip], x_dot[ip],
                          x[it], x_dot[it],
                          ebody1 );
    if (status != EGADS_SUCCESS) goto cleanup;

    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;


    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeNaca( context, stack, FACEBODY, sharpte, x[im], x[ip], x[it], &ebody2 );
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Ping NACA", 1e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
pingNacaRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  int    sharpte = 0;
  double x[3], x_dot[3], params[3], dtime = 1e-8;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2, eloop1, eloop2;
  enum naca { im, ip, it };

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  for (sharpte = 0; sharpte <= 1; sharpte++) {
    for (dir = -1; dir <= 1; dir += 2) {

      printf(" ---------------------------------\n");
      printf(" Ping Ruled NACA dir %+d sharpte %d\n", dir, sharpte);

      xform1[0] = 1.0;
      xform1[1] = 0.1;
      xform1[2] = 0.2;
      xform1[3] = 1.*dir;

      xform2[0] = 0.75;
      xform2[1] = 0.;
      xform2[2] = 0.;
      xform2[3] = 2.*dir;

      /* make the ruled body */
      status = makeNaca( context, stack, FACEBODY, sharpte, x[im], x[ip], x[it], &secs1[0] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNaca( context, stack, LOOP, sharpte, x[im], x[ip], x[it], &eloop1 );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop1, xform1, &secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled(nsec, secs1, &ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* test re-making the topology */
      status = remakeTopology(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* tessellation parameters */
      params[0] =  0.5;
      params[1] =  4.0;
      params[2] = 35.0;

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
        printf(" Ping Ruled NACA Edge %d np1 = %d\n", iedge+1, np1);
      }

      for (iface = 0; iface < nface; iface++) {
        status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                &nt1, &ts1, &tc1);
        if (status != EGADS_SUCCESS) goto cleanup;
        printf(" Ping Ruled NACA Face %d np1 = %d\n", iface+1, np1);
      }

      /* zero out velocities */
      for (iparam = 0; iparam < 3; iparam++) x_dot[iparam] = 0;

      for (iparam = 0; iparam < 3; iparam++) {

        /* set the velocity of the original body */
        x_dot[iparam] = 1.0;
        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              secs1[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              eloop1 );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( eloop1, xform1, xform_dot, secs1[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_ruled_dot(ebody1, nsec, secs1);
        if (status != EGADS_SUCCESS) goto cleanup;
        x_dot[iparam] = 0.0;

        status = EG_hasGeometry_dot(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;


        /* make a perturbed body for finite difference */
        x[iparam] += dtime;
        status = makeNaca( context, stack, FACE, sharpte, x[im], x[ip], x[it], &secs2[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeNaca( context, stack, LOOP, sharpte, x[im], x[ip], x[it], &eloop2 );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, eloop2, xform1, &secs2[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_ruled(nsec, secs2, &ebody2);
        if (status != EGADS_SUCCESS) goto cleanup;
        x[iparam] -= dtime;

        /* map the tessellation */
        status = EG_mapTessBody(tess1, ebody2, &tess2);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* ping the bodies */
        status = pingBodies(tess1, tess2, dtime, iparam, "Ping Ruled NACA", 5e-7, 5e-7, 1e-7);
        if (status != EGADS_SUCCESS) goto cleanup;

        EG_deleteObject(tess2);
        EG_deleteObject(ebody2);
      }

      EG_deleteObject(tess1);
      EG_deleteObject(ebody1);
    }
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
pingNacaBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  int    sharpte = 0;
  double x[7], x_dot[7], params[3], dtime = 1e-8;
  double *RC1, *RC1_dot, *RCn, *RCn_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[3], secs2[3], ebody1, ebody2, tess1, tess2, eloop1, eloop2;
  enum naca { im, ip, it };

  RC1 = x + 3;
  RCn = x + 5;

  RC1_dot = x_dot + 3;
  RCn_dot = x_dot + 5;

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  RC1[0] = 0;
  RC1[1] = 1;

  RCn[0] = 0;
  RCn[1] = 2;

  for (sharpte = 0; sharpte <= 1; sharpte++) {
    for (dir = -1; dir <= 1; dir += 2) {

      xform1[0] = 1.0;
      xform1[1] = 0.1;
      xform1[2] = 0.2;
      xform1[3] = 1.*dir;

      xform2[0] = 0.75;
      xform2[1] = 0.;
      xform2[2] = 0.;
      xform2[3] = 2.*dir;

      printf(" ---------------------------------\n");
      printf(" Ping Blend NACA dir %+d sharpte %d\n", dir, sharpte);

      /* make the ruled body */
      status = makeNaca( context, stack, FACEBODY, sharpte, x[im], x[ip], x[it], &secs1[0] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNaca( context, stack, LOOP, sharpte, x[im], x[ip], x[it], &eloop1 );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop1, xform1, &secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs1, RC1, RCn, &ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* test re-making the topology */
      status = remakeTopology(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* tessellation parameters */
      params[0] =  0.5;
      params[1] = 20.0;
      params[2] = 35.0;

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
        printf(" Ping Blend NACA Edge %d np1 = %d\n", iedge+1, np1);
      }

      for (iface = 0; iface < nface; iface++) {
        status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                &nt1, &ts1, &tc1);
        if (status != EGADS_SUCCESS) goto cleanup;
        printf(" Ping Blend NACA Face %d np1 = %d\n", iface+1, np1);
      }

      /* zero out velocities */
      for (iparam = 0; iparam < 7; iparam++) x_dot[iparam] = 0;

      for (iparam = 0; iparam < 7; iparam++) {
        if (iparam == 3) continue; // RC1[0] swithc (not a parameter)
        if (iparam == 5) continue; // RCn[0] swithc (not a parameter)

        /* set the velocity of the original body */
        x_dot[iparam] = 1.0;
        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              secs1[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              eloop1 );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( eloop1, xform1, xform_dot, secs1[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_blend_dot(ebody1, nsec, secs1, RC1, RC1_dot, RCn, RCn_dot);
        if (status != EGADS_SUCCESS) goto cleanup;
        x_dot[iparam] = 0.0;

        status = EG_hasGeometry_dot(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;


        /* make a perturbed body for finite difference */
        x[iparam] += dtime;
        status = makeNaca( context, stack, FACE, sharpte, x[im], x[ip], x[it], &secs2[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeNaca( context, stack, LOOP, sharpte, x[im], x[ip], x[it], &eloop2 );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, eloop2, xform1, &secs2[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_blend(nsec, secs2, RC1, RCn, &ebody2);
        if (status != EGADS_SUCCESS) goto cleanup;
        x[iparam] -= dtime;

        /* map the tessellation */
        status = EG_mapTessBody(tess1, ebody2, &tess2);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* ping the bodies */
        status = pingBodies(tess1, tess2, dtime, iparam, "Ping Blend NACA", 5e-7, 5e-7, 1e-7);
        if (status != EGADS_SUCCESS) goto cleanup;

        EG_deleteObject(tess2);
        EG_deleteObject(ebody2);
      }

      EG_deleteObject(tess1);
      EG_deleteObject(ebody1);
    }
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
equivNacaRuled(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  int    sharpte = 0;
  double x[3], x_dot[3], params[3];
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs[3], ebody1, ebody2, tess1, tess2, eloop;
  enum naca { im, ip, it };

  egadsSplineVels vels;

  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  for (sharpte = 0; sharpte <= 1; sharpte++) {
    for (dir = -1; dir <= 1; dir += 2) {

      printf(" ---------------------------------\n");
      printf(" Equiv Ruled NACA dir %+d sharpte %d\n", dir, sharpte);

      xform1[0] = 1.0;
      xform1[1] = 0.1;
      xform1[2] = 0.2;
      xform1[3] = 1.*dir;

      xform2[0] = 0.75;
      xform2[1] = 0.;
      xform2[2] = 0.;
      xform2[3] = 2.*dir;

      /* make the ruled body */
      status = makeNaca( context, stack, FACE, sharpte, x[im], x[ip], x[it], &secs[0] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNaca( context, stack, LOOP, sharpte, x[im], x[ip], x[it], &eloop );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop, xform1, &secs[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs[0], xform2, &secs[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_ruled(nsec, secs, &ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make the body for _vels sensitivities */
      status = EG_ruled(nsec, secs, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* tessellation parameters */
      params[0] =  0.5;
      params[1] =  4.0;
      params[2] = 35.0;

      /* make the tessellation */
      status = EG_makeTessBody(ebody1, params, &tess1);

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
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
        printf(" Equiv Ruled NACA Edge %d np1 = %d\n", iedge+1, np1);
      }

      for (iface = 0; iface < nface; iface++) {
        status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                &nt1, &ts1, &tc1);
        if (status != EGADS_SUCCESS) goto cleanup;
        printf(" Equiv Ruled NACA Face %d np1 = %d\n", iface+1, np1);
      }

      /* zero out velocities */
      for (iparam = 0; iparam < 3; iparam++) x_dot[iparam] = 0;

      for (iparam = 0; iparam < 3; iparam++) {

        /* set the velocity of the original body */
        x_dot[iparam] = 1.0;
        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              secs[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              eloop );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( eloop, xform1, xform_dot, secs[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs[0], xform2, xform_dot, secs[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_ruled_dot(ebody1, nsec, secs);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_hasGeometry_dot(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;


        /* set the velocity with _vels functions */
        status = EG_ruled_vels(nsec, secs, &vels, ebody2);
        if (status != EGADS_SUCCESS) goto cleanup;
        x_dot[iparam] = 0.0;


        /* ping the bodies */
        status = equivDotVels(tess1, tess2, iparam, "Equiv Ruled NACA", 1e-7, 1e-7, 1e-7);
        if (status != EGADS_SUCCESS) goto cleanup;
      }

      EG_deleteObject(tess2);
      EG_deleteObject(ebody2);

      EG_deleteObject(tess1);
      EG_deleteObject(ebody1);
    }
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
equivNacaBlend(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nsec = 3;
  int    sharpte = 0;
  double x[7], x_dot[7], params[3];
  double *RC1, *RC1_dot, *RCn, *RCn_dot;
  double xform1[4], xform2[4], xform_dot[4]={0,0,0,0};
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs[3], ebody1, ebody2, tess1, tess2, eloop;
  enum naca { im, ip, it };

  egadsSplineVels vels;

  vels.usrData = NULL;
  vels.velocityOfRange = &velocityOfRange;
  vels.velocityOfNode = &velocityOfNode;
  vels.velocityOfEdge = &velocityOfEdge;
  vels.velocityOfBspline = &velocityOfBspline;

  RC1 = x + 3;
  RCn = x + 5;

  RC1_dot = x_dot + 3;
  RCn_dot = x_dot + 5;

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  RC1[0] = 0;
  RC1[1] = 1;

  RCn[0] = 0;
  RCn[1] = 2;

  for (sharpte = 0; sharpte <= 1; sharpte++) {
    for (dir = -1; dir <= 1; dir += 2) {

      xform1[0] = 1.0;
      xform1[1] = 0.1;
      xform1[2] = 0.2;
      xform1[3] = 1.*dir;

      xform2[0] = 0.75;
      xform2[1] = 0.;
      xform2[2] = 0.;
      xform2[3] = 2.*dir;

      printf(" ---------------------------------\n");
      printf(" Equiv Blend NACA dir %+d sharpte %d\n", dir, sharpte);

      /* make the ruled body */
      status = makeNaca( context, stack, FACE, sharpte, x[im], x[ip], x[it], &secs[0] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeNaca( context, stack, LOOP, sharpte, x[im], x[ip], x[it], &eloop );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, eloop, xform1, &secs[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs[0], xform2, &secs[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = EG_blend(nsec, secs, RC1, RCn, &ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make the body for _vels sensitivities */
      status = EG_blend(nsec, secs, RC1, RCn, &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* tessellation parameters */
      params[0] =  0.5;
      params[1] = 20.0;
      params[2] = 35.0;

      /* make the tessellation */
      status = EG_makeTessBody(ebody1, params, &tess1);

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
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
        printf(" Equiv Blend NACA Edge %d np1 = %d\n", iedge+1, np1);
      }

      for (iface = 0; iface < nface; iface++) {
        status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                &nt1, &ts1, &tc1);
        if (status != EGADS_SUCCESS) goto cleanup;
        printf(" Equiv Blend NACA Face %d np1 = %d\n", iface+1, np1);
      }

      /* zero out velocities */
      for (iparam = 0; iparam < 7; iparam++) x_dot[iparam] = 0;

      for (iparam = 0; iparam < 7; iparam++) {
        if (iparam == 3) continue; // RC1[0] swithc (not a parameter)
        if (iparam == 5) continue; // RCn[0] swithc (not a parameter)

        /* set the velocity of the original body */
        x_dot[iparam] = 1.0;
        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              secs[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              eloop );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( eloop, xform1, xform_dot, secs[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs[0], xform2, xform_dot, secs[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_blend_dot(ebody1, nsec, secs, RC1, RC1_dot, RCn, RCn_dot);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_hasGeometry_dot(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;


        /* set the velocity with _vels functions */
        status = EG_blend_vels(nsec, secs, RC1, RC1_dot, RCn, RCn_dot, &vels, ebody2);
        if (status != EGADS_SUCCESS) goto cleanup;
        x_dot[iparam] = 0.0;


        /* ping the bodies */
        status = equivDotVels(tess1, tess2, iparam, "Equiv Blend NACA", 1e-7, 1e-7, 1e-7);
        if (status != EGADS_SUCCESS) goto cleanup;

      }
      EG_deleteObject(tess2);
      EG_deleteObject(ebody2);

      EG_deleteObject(tess1);
      EG_deleteObject(ebody1);
    }
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int makeSplineFaceBody( objStack *stack,  /* (in)  EGADS obj stack    */
                        ego esurf,        /* (in) b-spline surface    */
                        int sharpte,      /* (in) sharp trailing edge */
                        ego *ebody )      /* (out) result FaceBody    */
{
  int status = EGADS_SUCCESS;

  int i, oclass, mtype, *ivec=NULL, icurv[4][4], esens[4], lsens[1] = {SFORWARD};
  double *rvec=NULL, *rcurv[4]={NULL,NULL,NULL,NULL};
  double xyz[18], uv[2], tdata[2] = {0,1}, data[4];
  ego context, eref, ecurves[4], eedges[8], enodes[4], echild[2], eloop, eface, epcurvs[4];

  status = EG_getContext(esurf, &context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_stackPush(stack, esurf);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make Nodes */
  uv[0] = 0; uv[1] = 0;
  status = EG_evaluate(esurf, uv, xyz);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           xyz, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  uv[0] = 1; uv[1] = 0;
  status = EG_evaluate(esurf, uv, xyz);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           xyz, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  uv[0] = 1; uv[1] = 1;
  status = EG_evaluate(esurf, uv, xyz);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           xyz, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  uv[0] = 0; uv[1] = 1;
  status = EG_evaluate(esurf, uv, xyz);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           xyz, 0, NULL, NULL, &enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the v = 0 curve */
  status = EG_isoCurve(ivec, rvec,
                       -1, 0,
                       icurv[0], &rcurv[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the u = 1 curve */
  status = EG_isoCurve(ivec, rvec,
                       ivec[2]-1, -1,
                       icurv[1], &rcurv[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the v = 1 curve */
  status = EG_isoCurve(ivec, rvec,
                       -1, ivec[5]-1,
                       icurv[2], &rcurv[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the u = 0 curve */
  status = EG_isoCurve(ivec, rvec,
                       0, -1,
                       icurv[3], &rcurv[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make Curves */
  for (i = 0; i < 4; i++) {
    status = EG_makeGeometry(context, CURVE, BSPLINE,
                             NULL, icurv[i], rcurv[i], &ecurves[i]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, ecurves[i]);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  /* make Edges */
  if (sharpte == 1) {
    echild[0] = enodes[1];
    status = EG_makeTopology(context, ecurves[0], EDGE, ONENODE,
                             tdata, 1, echild, NULL, &eedges[0]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    echild[0] = enodes[1]; echild[1] = enodes[2];
    status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE,
                             tdata, 2, echild, NULL, &eedges[1]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    echild[0] = enodes[2];
    status = EG_makeTopology(context, ecurves[2], EDGE, ONENODE,
                             tdata, 1, echild, NULL, &eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    eedges[3] = eedges[1];

  } else {
    echild[0] = enodes[0]; echild[1] = enodes[1];
    status = EG_makeTopology(context, ecurves[0], EDGE, TWONODE,
                             tdata, 2, echild, NULL, &eedges[0]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[0]);
    if (status != EGADS_SUCCESS) goto cleanup;

    echild[0] = enodes[1]; echild[1] = enodes[2];
    status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE,
                             tdata, 2, echild, NULL, &eedges[1]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[1]);
    if (status != EGADS_SUCCESS) goto cleanup;

    echild[0] = enodes[3]; echild[1] = enodes[2];
    status = EG_makeTopology(context, ecurves[2], EDGE, TWONODE,
                             tdata, 2, echild, NULL, &eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[2]);
    if (status != EGADS_SUCCESS) goto cleanup;

    echild[0] = enodes[0]; echild[1] = enodes[3];
    status = EG_makeTopology(context, ecurves[3], EDGE, TWONODE,
                             tdata, 2, echild, NULL, &eedges[3]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(stack, eedges[3]);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  /* make p-curves */
  data[0] = 0.; data[1] = 0.; /* v == 0 VMIN */
  data[2] = 1.; data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &epcurvs[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, epcurvs[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 1.; data[1] = 0.; /* u == 1 UMAX */
  data[2] = 0.; data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &epcurvs[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, epcurvs[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.; data[1] = 1.; /* v == 1 VMAX */
  data[2] = 1.; data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &epcurvs[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, epcurvs[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.; data[1] = 0.; /* u == 0 UMIN */
  data[2] = 0.; data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &epcurvs[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, epcurvs[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  esens[0] = SFORWARD;
  esens[1] = SFORWARD;
  esens[2] = SREVERSE;
  esens[3] = SREVERSE;

  eedges[0+4] = epcurvs[0];
  eedges[1+4] = epcurvs[1];
  eedges[2+4] = epcurvs[2];
  eedges[3+4] = epcurvs[3];

  /* make the loop from the edges */
  status = EG_makeTopology(context, esurf, LOOP, CLOSED, NULL, 4, eedges, esens,
                         &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the face from the loop */
  status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL, 1, &eloop, lsens,
                           &eface);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eface);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the FaceBody */
  status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                           NULL, 1, &eface, NULL, ebody);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, *ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  EG_free(ivec);
  EG_free(rvec);
  for (i = 0; i < 4; i++) {
    EG_free(rcurv[i]);
  }
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int makeSplineFaceBody_dot( ego ebody )  /* (in/out) FaceBody with sensitvities */
{
  int status = EGADS_SUCCESS;

  int i, oclass, mtype, *ivec=NULL, icurv[4][4];
  int nchild, *senses;
  double *rvec=NULL, *rvec_dot=NULL, *rcurv[4]={NULL,NULL,NULL,NULL}, *rcurv_dot[4]={NULL,NULL,NULL,NULL};
  double xyz[18], xyz_dot[18];
  double uv[2], tdata[2] = {0,1}, tdata_dot[2] = {0,0}, data[4];
  ego *echildren=NULL, eref, esurf, ecurves[4], *eedges=NULL, enodes[4], eloop, eface;

    /* get the Face */
   status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                           data, &nchild, &echildren, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;

   /* get the Loop */
   eface = echildren[0];
   status = EG_getTopology(eface, &esurf, &oclass, &mtype,
                           data, &nchild, &echildren, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;
   eloop = echildren[0];

   /* get the Edges */
   status = EG_getTopology(eloop, &esurf, &oclass, &mtype,
                           data, &nchild, &eedges, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;

   /* get the Curves/Nodes */
   status = EG_getTopology(eedges[0], &ecurves[0], &oclass, &mtype,
                           data, &nchild, &echildren, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;
   enodes[0] = echildren[0]; enodes[1] = echildren[1];

   status = EG_getTopology(eedges[1], &ecurves[1], &oclass, &mtype,
                           data, &nchild, &echildren, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;
   enodes[1] = echildren[0];  enodes[2] = echildren[1];

   status = EG_getTopology(eedges[2], &ecurves[2], &oclass, &mtype,
                           data, &nchild, &echildren, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;
   enodes[3] = echildren[0];  enodes[2] = echildren[1];

   status = EG_getTopology(eedges[3], &ecurves[3], &oclass, &mtype,
                           data, &nchild, &echildren, &senses);
   if (status != EGADS_SUCCESS) goto cleanup;
   enodes[0] = echildren[0];  enodes[3] = echildren[1];



   /* set Node sensitivities */
   uv[0] = 0; uv[1] = 0;
   status = EG_evaluate_dot(esurf, uv, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   uv[0] = 1; uv[1] = 0;
   status = EG_evaluate_dot(esurf, uv, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   uv[0] = 1; uv[1] = 1;
   status = EG_evaluate_dot(esurf, uv, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   uv[0] = 0; uv[1] = 1;
   status = EG_evaluate_dot(esurf, uv, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;

   status = EG_setGeometry_dot(enodes[3], NODE, 0, NULL, xyz, xyz_dot);
   if (status != EGADS_SUCCESS) goto cleanup;


   status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
   if (status != EGADS_SUCCESS) goto cleanup;
   EG_free(rvec);

   status = EG_getGeometry_dot(esurf, &rvec, &rvec_dot);
   if (status != EGADS_SUCCESS) goto cleanup;


   /* get the v = 0 curve */
   status = EG_isoCurve_dot(ivec, rvec, rvec_dot,
                            -1, 0,
                            icurv[0], &rcurv[0], &rcurv_dot[0]);
   if (status != EGADS_SUCCESS) goto cleanup;

   /* get the u = 1 curve */
   status = EG_isoCurve_dot(ivec, rvec, rvec_dot,
                            ivec[2]-1, -1,
                            icurv[1], &rcurv[1], &rcurv_dot[1]);
   if (status != EGADS_SUCCESS) goto cleanup;

   /* get the v = 1 curve */
   status = EG_isoCurve_dot(ivec, rvec, rvec_dot,
                            -1, ivec[5]-1,
                            icurv[2], &rcurv[2], &rcurv_dot[2]);
   if (status != EGADS_SUCCESS) goto cleanup;

   /* get the u = 0 curve */
   status = EG_isoCurve_dot(ivec, rvec, rvec_dot,
                            0, -1,
                            icurv[3], &rcurv[3], &rcurv_dot[3]);
   if (status != EGADS_SUCCESS) goto cleanup;


   /* set Curve/Edge sensitivities */
   for (i = 0; i < 4; i++) {
     status = EG_setGeometry_dot(ecurves[i], CURVE, BSPLINE,
                                 icurv[i], rcurv[i], rcurv_dot[i]);
     if (status != EGADS_SUCCESS) goto cleanup;

     status = EG_setRange_dot(eedges[i], EDGE, tdata, tdata_dot);
     if (status != EGADS_SUCCESS) goto cleanup;
   }


cleanup:
  EG_free(ivec);
  EG_free(rvec);
  EG_free(rvec_dot);
  for (i = 0; i < 4; i++) {
    EG_free(rcurv[i]);
    EG_free(rcurv_dot[i]);
  }
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
pingNacaSkinned(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir, nchild, nsec = 4;
  int    sharpte = 0, oclass, mtype, *senses=NULL;
  double x[3], x_dot[3], params[3], data[4], dtime = 1e-8;
  double xform1[4], xform2[4], xform3[4], xform_dot[4]={0,0,0,0}, massProp[14];
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    secs1[4], secs2[4], ebody1, ebody2, tess1, tess2, esurf1, esurf2, eface, *echildren=NULL;
  enum naca { im, ip, it };

  const double areas[2] = {6.135844537797e+00,6.137250618148e+00};

  x[im] = 0.1;  /* camber */
  x[ip] = 0.4;  /* max loc */
  x[it] = 0.16; /* thickness */

  for (sharpte = 0; sharpte <= 1; sharpte++) {
    for (dir = -1; dir <= 1; dir += 2) {

      printf(" ---------------------------------\n");
      printf(" Ping Skinned NACA dir %+d sharpte %d\n", dir, sharpte);

      xform1[0] = 1.0;
      xform1[1] = 0.1;
      xform1[2] = 0.2;
      xform1[3] = 1.*dir;

      xform2[0] = 0.75;
      xform2[1] = 0.;
      xform2[2] = 0.;
      xform2[3] = 2.*dir;

      xform3[0] = 1.25;
      xform3[1] = 0.;
      xform3[2] = 0.;
      xform3[3] = 3.*dir;

      /* make the ruled body */
      status = makeNaca( context, stack, CURVE, sharpte, x[im], x[ip], x[it], &secs1[0] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs1[0], xform1, &secs1[1] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs1[0], xform2, &secs1[2] );
      if (status != EGADS_SUCCESS) goto cleanup;

      status = makeTransform(stack, secs1[0], xform3, &secs1[3] );
      if (status != EGADS_SUCCESS) goto cleanup;

      /* perform the skinning operation */
      status = EG_skinning(nsec, secs1, 3, &esurf1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make the FaceBody */
      status = makeSplineFaceBody(stack, esurf1, sharpte, &ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* test re-making the topology */
      status = remakeTopology(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* get the surface back */
      status = EG_getTopology(ebody1, &esurf1, &oclass, &mtype,
                              data, &nchild, &echildren, &senses);
      if (status != EGADS_SUCCESS) goto cleanup;

      eface = echildren[0];
      status = EG_getTopology(eface, &esurf1, &oclass, &mtype,
                              data, &nchild, &echildren, &senses);
      if (status != EGADS_SUCCESS) goto cleanup;


      status = EG_getMassProperties(ebody1, massProp);
      if (status != EGADS_SUCCESS) goto cleanup;

      if (fabs(massProp[1] - areas[sharpte]) > 1e-9) {
        printf("Skinning area failure! fabs(%+le - %+le) = %+le > %e\n",
               massProp[1], areas[sharpte], fabs(massProp[1] - areas[sharpte]), 1e-9);
        status = EGADS_GEOMERR;
        goto cleanup;
      }

      /* tessellation parameters */
      params[0] =  0.5;
      params[1] =  4.0;
      params[2] = 35.0;

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
        printf(" Ping Skinned NACA Edge %d np1 = %d\n", iedge+1, np1);
      }

      for (iface = 0; iface < nface; iface++) {
        status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                &nt1, &ts1, &tc1);
        if (status != EGADS_SUCCESS) goto cleanup;
        printf(" Ping Skinned NACA Face %d np1 = %d\n", iface+1, np1);
      }

      /* zero out velocities */
      for (iparam = 0; iparam < 3; iparam++) x_dot[iparam] = 0;

      for (iparam = 0; iparam < 3; iparam++) {

        /* set the velocity of the original body */
        x_dot[iparam] = 1.0;
        status = setNaca_dot( sharpte,
                              x[im], x_dot[im],
                              x[ip], x_dot[ip],
                              x[it], x_dot[it],
                              secs1[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs1[0], xform1, xform_dot, secs1[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs1[0], xform2, xform_dot, secs1[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = setTransform_dot( secs1[0], xform3, xform_dot, secs1[3] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_skinning_dot(esurf1, nsec, secs1);
        if (status != EGADS_SUCCESS) goto cleanup;
        x_dot[iparam] = 0.0;

        status = makeSplineFaceBody_dot(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_hasGeometry_dot(ebody1);
        if (status != EGADS_SUCCESS) goto cleanup;


        /* make a perturbed body for finite difference */
        x[iparam] += dtime;
        status = makeNaca( context, stack, CURVE, sharpte, x[im], x[ip], x[it], &secs2[0] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, secs2[0], xform1, &secs2[1] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, secs2[0], xform2, &secs2[2] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = makeTransform(stack, secs2[0], xform3, &secs2[3] );
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_skinning(nsec, secs2, 3, &esurf2);
        if (status != EGADS_SUCCESS) goto cleanup;
        x[iparam] -= dtime;

        /* make the FaceBody */
        status = makeSplineFaceBody(stack, esurf2, sharpte, &ebody2);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* map the tessellation */
        status = EG_mapTessBody(tess1, ebody2, &tess2);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* ping the bodies */
        status = pingBodies(tess1, tess2, dtime, iparam, "Ping Skinned NACA", 5e-7, 5e-7, 1e-7);
        if (status != EGADS_SUCCESS) goto cleanup;

        EG_deleteObject(tess2);
      }

      EG_deleteObject(tess1);
    }
  }

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

  /*-------*/
  status = pingNodeRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivNodeRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingNodeBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivNodeBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /*-------*/
  status = pingLineRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingLineBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivLineRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivLineBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /*-------*/
  status = pingLine2Ruled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingLine2Blend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /*-------*/
  status = pingCircleRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingCircleBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivCircleRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivCircleBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /*-------*/
  status = pingNoseCircleBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingNoseCircle2Blend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /*-------*/
  status = pingNaca(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingNacaRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingNacaBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivNacaRuled(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = equivNacaBlend(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /*-------*/
  status = pingNacaSkinned(context, &stack);
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
