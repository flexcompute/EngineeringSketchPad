#include <math.h>
#include "egads.h"
#include "egads_dot.h"

#include "../src/egadsStack.h"

#if defined(_MSC_VER) && (_MSC_VER < 1900)
#define __func__  __FUNCTION__
#endif

//#define PCURVE_SENSITIVITY

int
EG_isoCurve(const int *header2d, const double *data2d,
            const int ik, const int jk, int *header, double **data);

int
EG_isoCurve_dot(const int *header2d, const double *data2d, const double *data2d_dot,
                const int ik, const int jk, int *header, double **data, double **data_dot);


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
/*  pingBodiesExtern                                                         */
/*                                                                           */
/*****************************************************************************/

int
pingBodiesExtern(ego tess1, ego ebody2, double dtime, int iparam, const char *shape, double ftol, double etol, double ntol)
{
  int status = EGADS_SUCCESS;
  int iedge, nedge, iface, nface, oclass, mtype;
  int np1, nt1;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego ebody1, tess=NULL, tess2=NULL, *eedges1=NULL, *efaces1=NULL;
  ego top, prev, next;

  status = EG_statusTessBody( tess1, &ebody1, &oclass, &mtype );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebody1, NULL, EDGE, &nedge, &eedges1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Faces from the Body */
  status = EG_getBodyTopos(ebody1, NULL, FACE, &nface, &efaces1);
  if (status != EGADS_SUCCESS) goto cleanup;


  status = EG_initTessBody(ebody1, &tess);
  if (status != EGADS_SUCCESS) goto cleanup;

  for (iedge = 0; iedge < nedge; iedge++) {

    status = EG_getInfo(eedges1[iedge], &oclass, &mtype, &top, &prev, &next);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (mtype == DEGENERATE) continue;

    /* extract the tessellation from the edge */
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the tessellation on the edge */
    status = EG_setTessEdge(tess, iedge+1, np1, x1, t1);
    if (status != EGADS_SUCCESS) goto cleanup;
  }


  for (iface = 0; iface < nface; iface++) {

    /* extract the face tessellation */
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* set the face tessellation */
    status = EG_setTessFace(tess, iface+1, np1, x1, uv1,
                                           nt1, ts1);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  status = EG_statusTessBody(tess, &ebody1, &oclass, &mtype);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* map the tessellation */
  status = EG_mapTessBody(tess, ebody2, &tess2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* ping the bodies */
  status = pingBodies(tess, tess2, dtime, iparam, shape, ftol, etol, ntol);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(eedges1);
  EG_free(efaces1);
  EG_deleteObject(tess);
  EG_deleteObject(tess2);

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Transform                                                                */
/*                                                                           */
/*****************************************************************************/

int
makeTransformBody(ego ebody,      /* (in) the body to be transformed   */
                  double *xforms, /* (in) sequence of transformations  */
                  ego *ebodys)    /* (out) array of transformed bodies */
{
  int    status = EGADS_SUCCESS;
  double scale, offset[3], ax, ay, az, cosa, sina;
  ego    context, exform;

  double mat[12] = {1.00, 0.00, 0.00, 0.0,
                    0.00, 1.00, 0.00, 0.0,
                    0.00, 0.00, 1.00, 0.0};

  status = EG_getContext(ebody, &context);
  if (status != EGADS_SUCCESS) goto cleanup;

  ax = xforms[0];

  cosa = cos(ax);
  sina = sin(ax);

  mat[ 0] =    1.; mat[ 1] =    0.; mat[ 2] =    0.; mat[ 3] = 0.;
  mat[ 4] =    0.; mat[ 5] =  cosa; mat[ 6] = -sina; mat[ 7] = 0.;
  mat[ 8] =    0.; mat[ 9] =  sina; mat[10] =  cosa; mat[11] = 0.;

  status = EG_makeTransform(context, mat, &exform);  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_copyObject(ebody, exform, &ebodys[0]); if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_deleteObject(exform);                  if (status != EGADS_SUCCESS) goto cleanup;

  ay = xforms[1];

  cosa = cos(ay);
  sina = sin(ay);

  mat[ 0] =  cosa; mat[ 1] =    0.; mat[ 2] =  sina; mat[ 3] = 0.;
  mat[ 4] =    0.; mat[ 5] =    1.; mat[ 6] =    0.; mat[ 7] = 0.;
  mat[ 8] = -sina; mat[ 9] =    0.; mat[10] =  cosa; mat[11] = 0.;

  status = EG_makeTransform(context, mat, &exform);      if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_copyObject(ebodys[0], exform, &ebodys[1]); if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_deleteObject(exform);                      if (status != EGADS_SUCCESS) goto cleanup;

  az = xforms[2];

  cosa = cos(az);
  sina = sin(az);

  mat[ 0] =  cosa; mat[ 1] = -sina; mat[ 2] =    0.; mat[ 3] = 0.;
  mat[ 4] =  sina; mat[ 5] =  cosa; mat[ 6] =    0.; mat[ 7] = 0.;
  mat[ 8] =    0.; mat[ 9] =    0.; mat[10] =    1.; mat[11] = 0.;

  status = EG_makeTransform(context, mat, &exform);      if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_copyObject(ebodys[1], exform, &ebodys[2]); if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_deleteObject(exform);                      if (status != EGADS_SUCCESS) goto cleanup;

  scale     = xforms[3];
  offset[0] = xforms[4];
  offset[1] = xforms[5];
  offset[2] = xforms[6];

  mat[ 0] = scale; mat[ 1] =    0.; mat[ 2] =    0.; mat[ 3] = offset[0];
  mat[ 4] =    0.; mat[ 5] = scale; mat[ 6] =    0.; mat[ 7] = offset[1];
  mat[ 8] =    0.; mat[ 9] =    0.; mat[10] = scale; mat[11] = offset[2];

  status = EG_makeTransform(context, mat, &exform);      if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_copyObject(ebodys[2], exform, &ebodys[3]); if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_deleteObject(exform);                      if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
setTransformBody_dot(ego ebody,          /* (in) the body to be transformed   */
                     double *xforms,     /* (in) sequence of transformations  */
                     double *xforms_dot, /* (in) velocity of transformations  */
                     ego *ebodys)        /* (out) transformed bodies with velocities */
{
  int status = EGADS_SUCCESS;
  double scale, scale_dot, offset[3], offset_dot[3], ax, ax_dot, ay, ay_dot, az, az_dot;
  double cosa, sina, cosa_dot, sina_dot;

  double mat[12] = {1.00, 0.00, 0.00, 0.0,
                    0.00, 1.00, 0.00, 0.0,
                    0.00, 0.00, 1.00, 0.0};
  double mat_dot[12];

  ax     = xforms[0];
  ax_dot = xforms_dot[0];

  cosa = cos(ax);
  sina = sin(ax);

  cosa_dot = -sin(ax) * ax_dot;
  sina_dot =  cos(ax) * ax_dot;

  mat[ 0] =    1.; mat[ 1] =    0.; mat[ 2] =    0.; mat[ 3] = 0.;
  mat[ 4] =    0.; mat[ 5] =  cosa; mat[ 6] = -sina; mat[ 7] = 0.;
  mat[ 8] =    0.; mat[ 9] =  sina; mat[10] =  cosa; mat[11] = 0.;

  mat_dot[ 0] =    0.; mat_dot[ 1] =        0.; mat_dot[ 2] =        0.; mat_dot[ 3] = 0.;
  mat_dot[ 4] =    0.; mat_dot[ 5] =  cosa_dot; mat_dot[ 6] = -sina_dot; mat_dot[ 7] = 0.;
  mat_dot[ 8] =    0.; mat_dot[ 9] =  sina_dot; mat_dot[10] =  cosa_dot; mat_dot[11] = 0.;

  status = EG_copyGeometry_dot(ebody, mat, mat_dot, ebodys[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  ay     = xforms[1];
  ay_dot = xforms_dot[1];

  cosa = cos(ay);
  sina = sin(ay);

  cosa_dot = -sin(ay) * ay_dot;
  sina_dot =  cos(ay) * ay_dot;

  mat[ 0] =  cosa; mat[ 1] =    0.; mat[ 2] =  sina; mat[ 3] = 0.;
  mat[ 4] =    0.; mat[ 5] =    1.; mat[ 6] =    0.; mat[ 7] = 0.;
  mat[ 8] = -sina; mat[ 9] =    0.; mat[10] =  cosa; mat[11] = 0.;

  mat_dot[ 0] =  cosa_dot; mat_dot[ 1] =        0.; mat_dot[ 2] =  sina_dot; mat_dot[ 3] = 0.;
  mat_dot[ 4] =        0.; mat_dot[ 5] =        0.; mat_dot[ 6] =        0.; mat_dot[ 7] = 0.;
  mat_dot[ 8] = -sina_dot; mat_dot[ 9] =        0.; mat_dot[10] =  cosa_dot; mat_dot[11] = 0.;

  status = EG_copyGeometry_dot(ebodys[0], mat, mat_dot, ebodys[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  az     = xforms[2];
  az_dot = xforms_dot[2];

  cosa = cos(az);
  sina = sin(az);

  cosa_dot = -sin(az) * az_dot;
  sina_dot =  cos(az) * az_dot;

  mat[ 0] =  cosa; mat[ 1] = -sina; mat[ 2] =    0.; mat[ 3] = 0.;
  mat[ 4] =  sina; mat[ 5] =  cosa; mat[ 6] =    0.; mat[ 7] = 0.;
  mat[ 8] =    0.; mat[ 9] =    0.; mat[10] =    1.; mat[11] = 0.;

  mat_dot[ 0] =  cosa_dot; mat_dot[ 1] = -sina_dot; mat_dot[ 2] =        0.; mat_dot[ 3] = 0.;
  mat_dot[ 4] =  sina_dot; mat_dot[ 5] =  cosa_dot; mat_dot[ 6] =        0.; mat_dot[ 7] = 0.;
  mat_dot[ 8] =        0.; mat_dot[ 9] =        0.; mat_dot[10] =        0.; mat_dot[11] = 0.;

  status = EG_copyGeometry_dot(ebodys[1], mat, mat_dot, ebodys[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  scale     = xforms[3];
  offset[0] = xforms[4];
  offset[1] = xforms[5];
  offset[2] = xforms[6];

  scale_dot     = xforms_dot[3];
  offset_dot[0] = xforms_dot[4];
  offset_dot[1] = xforms_dot[5];
  offset_dot[2] = xforms_dot[6];

  mat[ 0] = scale; mat[ 1] =    0.; mat[ 2] =    0.; mat[ 3] = offset[0];
  mat[ 4] =    0.; mat[ 5] = scale; mat[ 6] =    0.; mat[ 7] = offset[1];
  mat[ 8] =    0.; mat[ 9] =    0.; mat[10] = scale; mat[11] = offset[2];

  mat_dot[ 0] = scale_dot; mat_dot[ 1] =        0.; mat_dot[ 2] =        0.; mat_dot[ 3] = offset_dot[0];
  mat_dot[ 4] =        0.; mat_dot[ 5] = scale_dot; mat_dot[ 6] =        0.; mat_dot[ 7] = offset_dot[1];
  mat_dot[ 8] =        0.; mat_dot[ 9] =        0.; mat_dot[10] = scale_dot; mat_dot[11] = offset_dot[2];

  status = EG_copyGeometry_dot(ebodys[2], mat, mat_dot, ebodys[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


int
pingTransform(ego ebody, double *params, const char *shape, double ftol, double etol, double ntol)
{
  int    status = EGADS_SUCCESS;
  int    i, j, np1, nt1, iedge, nedge, iface, nface;
  double dtime = 1e-7;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebodys1[4], ebodys2[4], tess1, tess2;

  double xforms[7] = {45.*PI/180., 30.*PI/180., 10.*PI/180.,
                      1.25, 1.00, 2.00, 3.0};
  double xforms_dot[7];

  /* make the transformed body */
  status = makeTransformBody(ebody, xforms, ebodys1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Faces from the Body */
  status = EG_getBodyTopos(ebodys1[3], NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Body */
  status = EG_getBodyTopos(ebodys1[3], NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  status = EG_makeTessBody(ebodys1[3], params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Trsf %s Edge %d np1 = %d\n", shape, iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Trsf %s Face %d np1 = %d\n", shape, iface+1, np1);
  }

  /* zero out velocities */
  for (i = 0; i < 7; i++) xforms_dot[i] = 0;

  for (i = 0; i < 7; i++) {

    /* set the velocity of the original body */
    xforms_dot[i] = 1.0;
    status = setTransformBody_dot(ebody, xforms, xforms_dot, ebodys1);
    if (status != EGADS_SUCCESS) goto cleanup;
    xforms_dot[i] = 0.0;

    status = EG_hasGeometry_dot(ebodys1[3]);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    xforms[i] += dtime;
    status = makeTransformBody(ebody, xforms, ebodys2);
    if (status != EGADS_SUCCESS) goto cleanup;
    xforms[i] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebodys2[3], &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, i, shape, ftol, etol, ntol);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
    for (j = 0; j < 4; j++)
      EG_deleteObject(ebodys2[j]);
  }

  EG_deleteObject(tess1);
  for (j = 0; j < 4; j++)
    EG_deleteObject(ebodys1[j]);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
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
pingLine(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[6], x_dot[6], *p1, *p2, *p1_dot, *p2_dot, params[3], dtime = 1e-7;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

  p1 = x;
  p2 = x+3;

  p1_dot = x_dot;
  p2_dot = x_dot+3;

  /* make the Line body */
  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 0.50; p2[1] = 0.75; p2[2] = 1.00;
  status = makeLineBody(context, stack, p1, p2, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.05;
  params[1] =  0.001;
  params[2] = 12.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Line np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 6; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 6; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setLineBody_dot(p1, p1_dot, p2, p2_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Line for finite difference */
    x[iparam] += dtime;
    status = makeLineBody(context, stack, p1, p2, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Line", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Line", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setLineBody_dot(p1, p1_dot, p2, p2_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Line", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;


  EG_deleteObject(tess1);

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
makeCircleBody( ego context,         /* (in)  EGADS context    */
                objStack *stack,     /* (in)  EGADS obj stack  */
                const double *xcent, /* (in)  Center           */
                const double *xax,   /* (in)  x-axis           */
                const double *yax,   /* (in)  y-axis           */
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
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
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
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];

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
setCircleBody_dot( const double *xcent,     /* (in)  Center          */
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
  int    nnode, nedge, nloop, oclass, mtype, *senses;
  double data[10], data_dot[10], dx[3], dx_dot[3], tdata[2], tdata_dot[2];
  double *rvec=NULL, *rvec_dot=NULL;
  ego    ecircle, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
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
  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];

  /* set the sensitivity of the Node */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL , data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(rvec); rvec = NULL;
  EG_free(rvec_dot); rvec_dot = NULL;

  return status;
}


int
pingCircle(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

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
  status = makeCircleBody(context, stack, xcent, xax, yax, x[9], &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Circle np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setCircleBody_dot(xcent, xcent_dot,
                               xax, xax_dot,
                               yax, yax_dot,
                               x[9], x_dot[9], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeCircleBody(context, stack, xcent, xax, yax, x[9], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Circle", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Circle", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setCircleBody_dot(xcent, xcent_dot,
                             xax, xax_dot,
                             yax, yax_dot,
                             x[9], x_dot[9], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Circle", 1e-7, 5e-7, 1e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Ellipse                                                                  */
/*                                                                           */
/*****************************************************************************/

int
makeEllipseBody( ego context,         /* (in)  EGADS context     */
                 objStack *stack,     /* (in)  EGADS obj stack   */
                 const double *xcent, /* (in)  Center            */
                 const double *xax,   /* (in)  x-axis            */
                 const double *yax,   /* (in)  y-axis            */
                 const double majr,   /* (in)  major radius      */
                 const double minr,   /* (in)  minor radius      */
                 ego *ebody )         /* (out) Ellipse wire body */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[11], tdata[2];
  ego    eellipse, eedge, enode, eloop;

  /* create the Ellipse */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = majr;     /* major radius */
  data[10] = minr;     /* minor radius */
  status = EG_makeGeometry(context, CURVE, ELLIPSE, NULL, NULL,
                           data, &eellipse);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eellipse);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Node for the Edge */
  data[0] = xcent[0] + xax[0]*majr;
  data[1] = xcent[1] + xax[1]*majr;
  data[2] = xcent[2] + xax[2]*majr;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enode);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, eellipse, EDGE, ONENODE,
                           tdata, 1, &enode, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, CLOSED,
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
setEllipseBody_dot( const double *xcent,     /* (in)  Center                */
                    const double *xcent_dot, /* (in)  Center velocity       */
                    const double *xax,       /* (in)  x-axis                */
                    const double *xax_dot,   /* (in)  x-axis velocity       */
                    const double *yax,       /* (in)  y-axis                */
                    const double *yax_dot,   /* (in)  y-axis velocity       */
                    const double majr,       /* (in)  major radius          */
                    const double majr_dot,   /* (in)  major radius velocity */
                    const double minr,       /* (in)  minor radius          */
                    const double minr_dot,   /* (in)  minor radius velocity */
                    ego ebody )              /* (in/out) Circle body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses;
  double data[11], data_dot[11], tdata[2], tdata_dot[2];
  ego    eellipse, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &eellipse, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Node */
  data[0] = xcent[0] + xax[0]*majr;
  data[1] = xcent[1] + xax[1]*majr;
  data[2] = xcent[2] + xax[2]*majr;
  data_dot[0] = xcent_dot[0] + xax_dot[0]*majr + xax[0]*majr_dot;
  data_dot[1] = xcent_dot[1] + xax_dot[1]*majr + xax[1]*majr_dot;
  data_dot[2] = xcent_dot[2] + xax_dot[2]*majr + xax[2]*majr_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Circle data and velocity */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = majr;     /* major radius */
  data[10] = minr;     /* minor radius */

  data_dot[ 0] = xcent_dot[0]; /* center */
  data_dot[ 1] = xcent_dot[1];
  data_dot[ 2] = xcent_dot[2];
  data_dot[ 3] = xax_dot[0];   /* x-axis */
  data_dot[ 4] = xax_dot[1];
  data_dot[ 5] = xax_dot[2];
  data_dot[ 6] = yax_dot[0];   /* y-axis */
  data_dot[ 7] = yax_dot[1];
  data_dot[ 8] = yax_dot[2];
  data_dot[ 9] = majr_dot;     /* major radius */
  data_dot[10] = minr_dot;     /* minor radius */

  status = EG_setGeometry_dot(eellipse, CURVE, ELLIPSE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingEllipse(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[11], x_dot[11], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* make the Circle body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.20; xax[1]   = 0.05; xax[2]   = 0.10;
  yax[0]   = 0.10; yax[1]   = 1.10; yax[2]   = 0.05;
  x[ 9] = 2.0;
  x[10] = 1.0;
  status = makeEllipseBody(context, stack, xcent, xax, yax, x[9], x[10], &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Ellipse np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 11; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 11; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setEllipseBody_dot(xcent, xcent_dot,
                                xax, xax_dot,
                                yax, yax_dot,
                                x[ 9], x_dot[ 9],
                                x[10], x_dot[10], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeEllipseBody(context, stack, xcent, xax, yax, x[9], x[10], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Ellipse", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Ellipse", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setEllipseBody_dot(xcent, xcent_dot,
                              xax, xax_dot,
                              yax, yax_dot,
                              x[ 9], x_dot[ 9],
                              x[10], x_dot[10], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Ellipse", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Parabola                                                                  */
/*                                                                           */
/*****************************************************************************/

int
makeParabolaBody( ego context,         /* (in)  EGADS context    */
                  objStack *stack,     /* (in)  EGADS obj stack  */
                  const double *xcent, /* (in)  Center           */
                  const double *xax,   /* (in)  x-axis           */
                  const double *yax,   /* (in)  y-axis           */
                  const double focus,  /* (in)  focus            */
                  ego *ebody )         /* (out) Circle wire body */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[10], tdata[2];
  ego    eparabola, eedge, enodes[2], eloop;

  /* create the Parabola */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = xax[0];   /* x-axis */
  data[4] = xax[1];
  data[5] = xax[2];
  data[6] = yax[0];   /* y-axis */
  data[7] = yax[1];
  data[8] = yax[2];
  data[9] = focus;    /* focus */
  status = EG_makeGeometry(context, CURVE, PARABOLA, NULL, NULL,
                           data, &eparabola);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eparabola);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Parabola */
  tdata[0] = -1;
  tdata[1] =  1;

  /* create the Nodes */
  status = EG_evaluate(eparabola, &tdata[0], data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate(eparabola, &tdata[1], data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge */
  status = EG_makeTopology(context, eparabola, EDGE, TWONODE,
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
setParabolaBody_dot( const double *xcent,   /* (in)  Center          */
                   const double *xcent_dot, /* (in)  Center velocity */
                   const double *xax,       /* (in)  x-axis          */
                   const double *xax_dot,   /* (in)  x-axis velocity */
                   const double *yax,       /* (in)  y-axis          */
                   const double *yax_dot,   /* (in)  y-axis velocity */
                   const double focus,      /* (in)  focus           */
                   const double focus_dot,  /* (in)  focus velocity  */
                   ego ebody )              /* (in/out) Circle body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses;
  double data[18], data_dot[18], tdata[2], tdata_dot[2];
  ego    eparabola, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &eparabola, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = -1;
  tdata[1]     =  1;
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
  data[9] = focus;    /* focus */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = xax_dot[0];   /* x-axis */
  data_dot[4] = xax_dot[1];
  data_dot[5] = xax_dot[2];
  data_dot[6] = yax_dot[0];   /* y-axis */
  data_dot[7] = yax_dot[1];
  data_dot[8] = yax_dot[2];
  data_dot[9] = focus_dot;    /* focus */

  status = EG_setGeometry_dot(eparabola, CURVE, PARABOLA, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Nodes */
  tdata[0] = -1;
  tdata[1] =  1;

  status = EG_evaluate_dot(eparabola, &tdata[0], NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate_dot(eparabola, &tdata[1], NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingParabola(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

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
  status = makeParabolaBody(context, stack, xcent, xax, yax, x[9], &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Parabola np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setParabolaBody_dot(xcent, xcent_dot,
                                 xax, xax_dot,
                                 yax, yax_dot,
                                 x[9], x_dot[9], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeParabolaBody(context, stack, xcent, xax, yax, x[9], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Parabola", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Parabola", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setParabolaBody_dot(xcent, xcent_dot,
                               xax, xax_dot,
                               yax, yax_dot,
                               x[9], x_dot[9], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Parabola", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Hyperbola                                                                */
/*                                                                           */
/*****************************************************************************/

int
makeHyperbolaBody( ego context,         /* (in)  EGADS context       */
                   objStack *stack,     /* (in)  EGADS obj stack     */
                   const double *xcent, /* (in)  Center              */
                   const double *xax,   /* (in)  x-axis              */
                   const double *yax,   /* (in)  y-axis              */
                   const double majr,   /* (in)  major radius        */
                   const double minr,   /* (in)  minor radius        */
                   ego *ebody )         /* (out) Hyperbola wire body */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[11], tdata[2];
  ego    ehyperbola, eedge, enodes[2], eloop;

  /* create the Hyperbola */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = majr;     /* major radius */
  data[10] = minr;     /* minor radius */
  status = EG_makeGeometry(context, CURVE, HYPERBOLA, NULL, NULL,
                           data, &ehyperbola);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ehyperbola);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Parabola */
  tdata[0] = -1;
  tdata[1] =  1;

  /* create the Nodes */
  status = EG_evaluate(ehyperbola, &tdata[0], data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate(ehyperbola, &tdata[1], data);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge */
  status = EG_makeTopology(context, ehyperbola, EDGE, TWONODE,
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
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}

int
setHyperbolaBody_dot( const double *xcent,     /* (in)  Center                */
                      const double *xcent_dot, /* (in)  Center velocity       */
                      const double *xax,       /* (in)  x-axis                */
                      const double *xax_dot,   /* (in)  x-axis velocity       */
                      const double *yax,       /* (in)  y-axis                */
                      const double *yax_dot,   /* (in)  y-axis velocity       */
                      const double majr,       /* (in)  major radius          */
                      const double majr_dot,   /* (in)  major radius velocity */
                      const double minr,       /* (in)  minor radius          */
                      const double minr_dot,   /* (in)  minor radius velocity */
                      ego ebody )              /* (in/out) Circle body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses;
  double data[11], data_dot[11], tdata[2], tdata_dot[2];
  ego    ehyperbola, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Node and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ehyperbola, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = -1;
  tdata[1]     =  1;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Circle data and velocity */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = majr;     /* major radius */
  data[10] = minr;     /* minor radius */

  data_dot[ 0] = xcent_dot[0]; /* center */
  data_dot[ 1] = xcent_dot[1];
  data_dot[ 2] = xcent_dot[2];
  data_dot[ 3] = xax_dot[0];   /* x-axis */
  data_dot[ 4] = xax_dot[1];
  data_dot[ 5] = xax_dot[2];
  data_dot[ 6] = yax_dot[0];   /* y-axis */
  data_dot[ 7] = yax_dot[1];
  data_dot[ 8] = yax_dot[2];
  data_dot[ 9] = majr_dot;     /* major radius */
  data_dot[10] = minr_dot;     /* minor radius */

  status = EG_setGeometry_dot(ehyperbola, CURVE, HYPERBOLA, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Nodes */
  tdata[0] = -1;
  tdata[1] =  1;

  status = EG_evaluate_dot(ehyperbola, &tdata[0], NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_evaluate_dot(ehyperbola, &tdata[1], NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingHyperbola(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[11], x_dot[11], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* make the Hyperbola body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.20; xax[1]   = 0.05; xax[2]   = 0.10;
  yax[0]   = 0.10; yax[1]   = 1.10; yax[2]   = 0.05;
  x[ 9] = 2.0;
  x[10] = 1.0;
  status = makeHyperbolaBody(context, stack, xcent, xax, yax, x[9], x[10], &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Hyperbola np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 11; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 11; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setHyperbolaBody_dot(xcent, xcent_dot,
                                  xax, xax_dot,
                                  yax, yax_dot,
                                  x[ 9], x_dot[ 9],
                                  x[10], x_dot[10], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeHyperbolaBody(context, stack, xcent, xax, yax, x[9], x[10], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Hyperbola", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Hyperbola", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setHyperbolaBody_dot(xcent, xcent_dot,
                                xax, xax_dot,
                                yax, yax_dot,
                                x[ 9], x_dot[ 9],
                                x[10], x_dot[10], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Hyperbola", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Offset Curve                                                             */
/*                                                                           */
/*****************************************************************************/

int
makeOffsetCurveBody( ego context,         /* (in)  EGADS context                   */
                     objStack *stack,     /* (in)  EGADS obj stack                 */
                     const double *x0,    /* (in)  coordinates of the first point  */
                     const double *x1,    /* (in)  coordinates of the second point */
                     const double *vec,   /* (in)  offset vector                   */
                     const double offset, /* (in)  offset magnitude                */
                     ego *ebody )         /* (out) Wire body created from points   */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD};
  double data[6], tdata[2];
  ego    eline, ecurve, eedge, enodes[2], eloop;

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

  /* t-data for the edge */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  /* create the Offset curve */
  data[0] = vec[0];
  data[1] = vec[1];
  data[2] = vec[2];
  data[3] = offset;

  status = EG_makeGeometry(context, CURVE, OFFSET, eline, NULL,
                           data, &ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Offset curve */
  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
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
setOffsetCurveBody_dot( const double *x0,        /* (in)  coordinates of the first point  */
                        const double *x0_dot,    /* (in)  velocity of the first point     */
                        const double *x1,        /* (in)  coordinates of the second point */
                        const double *x1_dot,    /* (in)  velocity of the second point    */
                        const double *vec,       /* (in)  offset vector                   */
                        const double *vec_dot,   /* (in)  offset vector velocity          */
                        const double offset,     /* (in)  offset magnitude                */
                        const double offset_dot, /* (in)  offset magnitude   velocity     */
                        ego ebody )              /* (in/out) Line body with velocities    */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses, *ivec=NULL;
  double data[6], data_dot[6], *rvec=NULL, tdata[2], tdata_dot[2];
  ego    ecurve, eline, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Curve from the Edge */
  status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Line from the Trimmed curve */
  status = EG_getGeometry(ecurve, &oclass, &mtype, &eline, &ivec, &rvec);
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

  /* set the Edge t-range sensitivity */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  tdata_dot[0] = 0;
  tdata_dot[1] = (data[3]*data_dot[3] + data[4]*data_dot[4] + data[5]*data_dot[5])/tdata[1];

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the offset curve velocity */
  data[0] = vec[0];
  data[1] = vec[1];
  data[2] = vec[2];
  data[3] = offset;

  data_dot[0] = vec_dot[0];
  data_dot[1] = vec_dot[1];
  data_dot[2] = vec_dot[2];
  data_dot[3] = offset_dot;

  status = EG_setGeometry_dot(ecurve, CURVE, OFFSET, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  EG_free(ivec);
  EG_free(rvec);

  return status;
}

int
pingOffsetCurve(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double x[10], x_dot[10], *p1, *p2, *vec, *p1_dot, *p2_dot, *vec_dot;
  double params[3], dtime = 1e-8;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

  p1  = x;
  p2  = x+3;
  vec = x+6;

  p1_dot  = x_dot;
  p2_dot  = x_dot+3;
  vec_dot = x_dot+6;

  /* make the Line body */
  p1[0] = 0.00; p1[1] = 0.00; p1[2] = 0.00;
  p2[0] = 0.50; p2[1] = 0.75; p2[2] = 1.00;
  vec[0] = 1.0; vec[1] = 2.0; vec[2] = 3.0;
  x[9] = 1.1;
  status = makeOffsetCurveBody(context, stack, p1, p2, vec, x[9], &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.05;
  params[1] =  1.0;
  params[2] = 12.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Offset Curve np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setOffsetCurveBody_dot(p1, p1_dot,
                                    p2, p2_dot,
                                    vec, vec_dot,
                                    x[9], x_dot[9], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Line for finite difference */
    x[iparam] += dtime;
    status = makeOffsetCurveBody(context, stack, p1, p2, vec, x[9], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Offset Curve", 1e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Offset Curve", 1e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setOffsetCurveBody_dot(p1, p1_dot,
                                  p2, p2_dot,
                                  vec, vec_dot,
                                  x[9], x_dot[9], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Offset Curve", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Bezier Curve                                                             */
/*                                                                           */
/*****************************************************************************/

int
makeBezierCurveBody( ego context,        /* (in)  EGADS context                        */
                     objStack *stack,    /* (in)  EGADS obj stack                      */
                     const int npts,     /* (in)  number of points                     */
                     const double *pts,  /* (in)  coordinates                          */
                     ego *ebody )        /* (out) Bezier wire body created from points */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD}, header[3];
  double data[6], tdata[2];
  ego    ecurve, eedge, enodes[2], eloop;

  /* make the bezier curve */
  header[0] = 0;
  header[1] = npts-1;
  header[2] = npts;

  /* create the Bezier curve */
  status = EG_makeGeometry(context, CURVE, BEZIER, NULL,
                           header, pts, &ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Nodes for the Edge */
  data[0] = pts[0];
  data[1] = pts[1];
  data[2] = pts[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = pts[3*(npts-1)+0];
  data[1] = pts[3*(npts-1)+1];
  data[2] = pts[3*(npts-1)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Curve */
  tdata[0] = 0;
  tdata[1] = 1;

  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
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
setBezierCurveBody_dot( const int npts,        /* (in)  number of points             */
                        const double *pts,     /* (in)  coordinates                  */
                        const double *pts_dot, /* (in)  velocity of coordinates      */
                        ego ebody )            /* (in/out) Line body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses, header[3];
  double data[6], data_dot[6], tdata[2], tdata_dot[2];
  ego    ecurve, *enodes, *eloops, *eedges, eref;

  /* the bezier curve header */
  header[0] = 0;
  header[1] = npts-1;
  header[2] = npts;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Curve from the Edge */
  status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = 1;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Nodes */
  data[0] = pts[0];
  data[1] = pts[1];
  data[2] = pts[2];
  data_dot[0] = pts_dot[0];
  data_dot[1] = pts_dot[1];
  data_dot[2] = pts_dot[2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = pts[3*(npts-1)+0];
  data[1] = pts[3*(npts-1)+1];
  data[2] = pts[3*(npts-1)+2];
  data_dot[0] = pts_dot[3*(npts-1)+0];
  data_dot[1] = pts_dot[3*(npts-1)+1];
  data_dot[2] = pts_dot[3*(npts-1)+2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Set the Bezier velocity */
  status = EG_setGeometry_dot(ecurve, CURVE, BEZIER, header, pts, pts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}

int
pingBezierCurve(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1;
  double params[3], dtime = 1e-7;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

  const int npts = 4;
  double pts[12] = {0.00, 0.00, 0.00,
                    1.00, 0.00, 0.10,
                    1.50, 1.00, 0.70,
                    0.25, 0.75, 0.60};
  double pts_dot[12];

  /* make the Bezier body */
  status = makeBezierCurveBody(context, stack, npts, pts, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.01;
  params[2] = 12.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" Bezier np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 3*npts; iparam++) pts_dot[iparam] = 0;

  for (iparam = 0; iparam < 3*npts; iparam++) {

    /* set the velocity of the original body */
    pts_dot[iparam] = 1.0;
    status = setBezierCurveBody_dot(npts, pts, pts_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    pts[iparam] += dtime;
    status = makeBezierCurveBody(context, stack, npts, pts, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Bezier", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Bezier", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setBezierCurveBody_dot(npts, pts, pts_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Bezier", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  B-spline Curve                                                           */
/*                                                                           */
/*****************************************************************************/

/* set KNOTS to 0 for arc-length knots, and -1 for equally spaced knots
 */
#define           KNOTS           0

/* tolerance for the spline fit */
#define           DXYTOL          1.0e-8

int
makeBsplineCurveBody( ego context,        /* (in)  EGADS context                        */
                      objStack *stack,    /* (in)  EGADS obj stack                      */
                      const int npts,     /* (in)  number of points                     */
                      const double *pts,  /* (in)  coordinates                          */
                      ego *ebody )        /* (out) Spline wire body created from points */
{
  int    status = EGADS_SUCCESS;
  int    senses[1] = {SFORWARD}, sizes[2];
  double data[6], tdata[2];
  ego    ecurve, eedge, enodes[2], eloop;

  /* create the B-spline curve */
  sizes[0] = npts;
  sizes[1] = KNOTS;
  status = EG_approximate(context, 0, DXYTOL, sizes, pts, &ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurve);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Nodes for the Edge */
  data[0] = pts[0];
  data[1] = pts[1];
  data[2] = pts[2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = pts[3*(npts-1)+0];
  data[1] = pts[3*(npts-1)+1];
  data[2] = pts[3*(npts-1)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Curve */
  tdata[0] = 0;
  tdata[1] = 1;

  status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
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
setBsplineCurveBody_dot( const int npts,        /* (in)  number of points             */
                         const double *pts,     /* (in)  coordinates                  */
                         const double *pts_dot, /* (in)  velocity of coordinates      */
                         ego ebody )            /* (in/out) Line body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, oclass, mtype, *senses, sizes[2];
  double data[6], data_dot[6], tdata[2], tdata_dot[2];
  ego    ecurve, *enodes, *eloops, *eedges, eref;

  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Curve from the Edge */
  status = EG_getTopology(eedges[0], &ecurve, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = 1;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the sensitivity of the Nodes */
  data[0] = pts[0];
  data[1] = pts[1];
  data[2] = pts[2];
  data_dot[0] = pts_dot[0];
  data_dot[1] = pts_dot[1];
  data_dot[2] = pts_dot[2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = pts[3*(npts-1)+0];
  data[1] = pts[3*(npts-1)+1];
  data[2] = pts[3*(npts-1)+2];
  data_dot[0] = pts_dot[3*(npts-1)+0];
  data_dot[1] = pts_dot[3*(npts-1)+1];
  data_dot[2] = pts_dot[3*(npts-1)+2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* Compute the B-spline velocity */
  sizes[0] = npts;
  sizes[1] = KNOTS;
  status = EG_approximate_dot(ecurve, 0, DXYTOL, sizes, pts, pts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  return status;
}


int
pingBsplineCurve(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    i, iparam, np1;
  double params[3], dtime = 1e-7;
  const double *t1, *x1;
  ego    ebody1, ebody2, tess1, tess2;

  double p[] = {0, 0.1, 0.2, 0.5, 0.8, 0.9, 1.0};
#define NP (sizeof(p)/sizeof(double))

  const int npts = NP;
  double pts[3*NP];
  double pts_dot[3*NP];
#undef NP

  /* create half circle points for the fit */
  for (i = 0; i < npts; i++) {
    pts[3*i  ] = cos(PI*p[i]);
    pts[3*i+1] = sin(PI*p[i]);
    pts[3*i+2] = 0.0;
  }

  /* make the B-spline body */
  status = makeBsplineCurveBody(context, stack, npts, pts, &ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* test re-making the topology */
  status = remakeTopology(ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the tessellation */
  params[0] =  0.1;
  params[1] =  0.01;
  params[2] = 12.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  /* extract the tessellation from the edge */
  status = EG_getTessEdge(tess1, 1, &np1, &x1, &t1);
  if (status != EGADS_SUCCESS) goto cleanup;

  printf(" B-spline np1 = %d\n", np1);

  /* zero out velocities */
  for (iparam = 0; iparam < 3*npts; iparam++) pts_dot[iparam] = 0;

  for (iparam = 0; iparam < 3*npts; iparam++) {

    /* set the velocity of the original body */
    pts_dot[iparam] = 1.0;
    status = setBsplineCurveBody_dot(npts, pts, pts_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Line for finite difference */
    pts[iparam] += dtime;
    status = makeBsplineCurveBody(context, stack, npts, pts, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "B-spline Curve", 1e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "B-spline Curve", 1e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setBsplineCurveBody_dot(npts, pts, pts_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "B-spline Curve", 1e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

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

  /* get the axes for the plane */
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
pingPlane(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* make the Plane Face body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  status = makePlaneBody(context, stack, xcent, xax, yax, &ebody1);
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

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setPlaneBody_dot(xcent, xcent_dot,
                              xax, xax_dot,
                              yax, yax_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makePlaneBody(context, stack, xcent, xax, yax, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Plane", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Plane", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setPlaneBody_dot(xcent, xcent_dot,
                            xax, xax_dot,
                            yax, yax_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Plane", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Spherical                                                                */
/*                                                                           */
/*****************************************************************************/

int
makeSphericalBody( ego context,         /* (in)  EGADS context   */
                   objStack *stack,     /* (in)  EGADS obj stack */
                   const double *xcent, /* (in)  Center          */
                   const double *xax,   /* (in)  x-axis          */
                   const double *yax,   /* (in)  y-axis          */
                   double r,            /* (in)  radius          */
                   ego *ebody )         /* (out) Spherical body  */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[10], tdata[2], dx[3], dy[3], dz[3], *rvec=NULL;
  ego    esphere, ecircle, eedges[8], enodes[2], eloop, eface, eref;

  /* create the Spherical */
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
  status = EG_makeGeometry(context, SURFACE, SPHERICAL, NULL, NULL,
                           data, &esphere);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, esphere);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(esphere, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the sphere */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];
  r     = rvec[9];

  /* compute the axis of the sphere between the poles */
  CROSS(dz, dx, dy);

  /* create the Circle curve for the edge */
  data[0] = rvec[0]; /* center */
  data[1] = rvec[1];
  data[2] = rvec[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dz[0];   /* y-axis */
  data[7] = dz[1];
  data[8] = dz[2];
  data[9] = r;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Nodes for the Edge */
  data[0] = xcent[0] - dz[0]*r;
  data[1] = xcent[1] - dz[1]*r;
  data[2] = xcent[2] - dz[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dz[0]*r;
  data[1] = xcent[1] + dz[1]*r;
  data[2] = xcent[2] + dz[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = -PI/2.;
  tdata[1] =  PI/2.;

  status = EG_makeTopology(context, ecircle, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Degenerate Edges on the Nodes */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &enodes[0], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[2] = eedges[0]; /* repeat the circle edge */

  status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &enodes[1], NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = -PI/2; /* v == -PI/2 VMIN */
  data[2] = 1.;    data[3] =  0.  ;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] =  PI/2.; /* v == PI/2 VMAX */
  data[2] = 1.;    data[3] =  0.   ;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the Loop and Face body */
  status = EG_makeTopology(context, esphere, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, esphere, FACE, SFORWARD,
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
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setSphericalBody_dot( const double *xcent,     /* (in)  Center          */
                      const double *xcent_dot, /* (in)  Center velocity */
                      const double *xax,       /* (in)  x-axis          */
                      const double *xax_dot,   /* (in)  x-axis velocity */
                      const double *yax,       /* (in)  y-axis          */
                      const double *yax_dot,   /* (in)  y-axis velocity */
                      double r,                /* (in)  radius          */
                      double r_dot,            /* (in)  radius velocity */
                      ego ebody )              /* (in/out) Circle body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3], dz[3], dz_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, tdata[2], tdata_dot[2];
  ego    esphere, ecircle, *enodes, *eloops, *eedges, *efaces, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Sphere from the Face */
  status = EG_getTopology(efaces[0], &esphere, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ecircle, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Edge t-range sensitivity */
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  tdata[0]     = -PI/2;
  tdata[1]     =  PI/2;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Sphere data and velocity */
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

  status = EG_setGeometry_dot(esphere, SURFACE, SPHERICAL, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(esphere, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the sphere */
  dx[0] = rvec[3];
  dx[1] = rvec[4];
  dx[2] = rvec[5];
  dy[0] = rvec[6];
  dy[1] = rvec[7];
  dy[2] = rvec[8];
  r     = rvec[9];

  dx_dot[0] = rvec_dot[3];
  dx_dot[1] = rvec_dot[4];
  dx_dot[2] = rvec_dot[5];
  dy_dot[0] = rvec_dot[6];
  dy_dot[1] = rvec_dot[7];
  dy_dot[2] = rvec_dot[8];
  r_dot     = rvec_dot[9];

  /* compute the axis of the sphere between the poles */
  CROSS(dz, dx, dy);
  CROSS_DOT(dz_dot, dx, dx_dot, dy, dy_dot);

  /* set the Circle curve velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dz[0];    /* y-axis */
  data[7] = dz[1];
  data[8] = dz[2];
  data[9] = r;         /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dz_dot[0];    /* y-axis */
  data_dot[7] = dz_dot[1];
  data_dot[8] = dz_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the sensitivity of the Nodes */
  data[0] = xcent[0] - dz[0]*r;
  data[1] = xcent[1] - dz[1]*r;
  data[2] = xcent[2] - dz[2]*r;
  data_dot[0] = xcent_dot[0] - dz_dot[0]*r - dz[0]*r_dot;
  data_dot[1] = xcent_dot[1] - dz_dot[1]*r - dz[1]*r_dot;
  data_dot[2] = xcent_dot[2] - dz_dot[2]*r - dz[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dz[0]*r;
  data[1] = xcent[1] + dz[1]*r;
  data[2] = xcent[2] + dz[2]*r;
  data_dot[0] = xcent_dot[0] + dz_dot[0]*r + dz[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dz_dot[1]*r + dz[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dz_dot[2]*r + dz[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
  /* set P-curve sensitivities */
  data_dot[0] = data_dot[1] = data_dot[2] = data_dot[3] = 0;

  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_setGeometry_dot(eedges[4+0], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = -PI/2; /* v == -PI/2 VMIN */
  data[2] = 1.;    data[3] =  0.  ;
  status = EG_setGeometry_dot(eedges[4+1], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_setGeometry_dot(eedges[4+2], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] =  PI/2.; /* v == PI/2 VMAX */
  data[2] = 1.;    data[3] =  0.   ;
  status = EG_setGeometry_dot(eedges[4+3], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
#endif

  status = EGADS_SUCCESS;

cleanup:
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingSpherical(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface, dir;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* check both directions */
  for (dir = -1; dir < 2; dir+=2) {

    /* make the Spherical body */
    xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
    xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
    yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
    //xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
    //yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
    x[9] = 1.0*dir;
    status = makeSphericalBody(context, stack, xcent, xax, yax, x[9], &ebody1);
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
    params[1] =  0.1;
    params[2] = 20.0;
    status = EG_makeTessBody(ebody1, params, &tess1);

    for (iedge = 0; iedge < nedge; iedge++) {
      status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Spherical Edge %d np1 = %d\n", iedge+1, np1);
    }

    for (iface = 0; iface < nface; iface++) {
      status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                              &nt1, &ts1, &tc1);
      if (status != EGADS_SUCCESS) goto cleanup;
      printf(" Spherical Face %d np1 = %d\n", iface+1, np1);
    }

    /* zero out velocities */
    for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

    for (iparam = 0; iparam < 10; iparam++) {

      /* set the velocity of the original body */
      x_dot[iparam] = 1.0;
      status = setSphericalBody_dot(xcent, xcent_dot,
                                    xax, xax_dot,
                                    yax, yax_dot,
                                    x[9], x_dot[9], ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;
      x_dot[iparam] = 0.0;

      status = EG_hasGeometry_dot(ebody1);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* make a perturbed Circle for finite difference */
      x[iparam] += dtime;
      status = makeSphericalBody(context, stack, xcent, xax, yax, x[9], &ebody2);
      if (status != EGADS_SUCCESS) goto cleanup;
      x[iparam] -= dtime;

      /* map the tessellation */
      status = EG_mapTessBody(tess1, ebody2, &tess2);
      if (status != EGADS_SUCCESS) goto cleanup;

      /* ping the bodies */
      status = pingBodies(tess1, tess2, dtime, iparam, "Spherical", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;

      status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Spherical", 1e-7, 1e-7, 1e-7);
      if (status != EGADS_SUCCESS) goto cleanup;

      EG_deleteObject(tess2);
    }

    /* zero out sensitivities */
    status = setSphericalBody_dot(xcent, xcent_dot,
                                  xax, xax_dot,
                                  yax, yax_dot,
                                  x[9], x_dot[9], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* check transformations */
    status = pingTransform(ebody1, params, "Spherical", 5e-7, 5e-7, 5e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess1);
  }

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Conical                                                                  */
/*                                                                           */
/*****************************************************************************/

int
makeConicalBody( ego context,         /* (in)  EGADS context   */
                 objStack *stack,     /* (in)  EGADS obj stack */
                 const double *xcent, /* (in)  Center          */
                 const double *xax,   /* (in)  x-axis          */
                 const double *yax,   /* (in)  y-axis          */
                 const double *zax,   /* (in)  z-axis          */
                 const double angle,  /* (in)  angle           */
                 const double r,      /* (in)  radius          */
                 ego *ebody )         /* (out) Spherical body  */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[14], tdata[2], dx[3], dy[3], dz[3], x1[3], x2[3], *rvec=NULL, vmin, h;
  ego    econe, ecircle, eline, eedges[8], enodes[2], eloop, eface, eref;

  /* create the Conical */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = zax[0];   /* z-axis */
  data[10] = zax[1];
  data[11] = zax[2];
  data[12] = angle;    /* angle  */
  data[13] = r;        /* radius */
  status = EG_makeGeometry(context, SURFACE, CONICAL, NULL, NULL,
                           data, &econe);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, econe);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(econe, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the cone */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];
  dz[0] = rvec[ 9];
  dz[1] = rvec[10];
  dz[2] = rvec[11];

  /* create the Circle curve for the base */
  data[0] = rvec[0]; /* center */
  data[1] = rvec[1];
  data[2] = rvec[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;

  vmin = -r/sin(angle);
  h = vmin*cos(angle);

  /* create the Node on the Tip */
  x1[0] = xcent[0] + dz[0]*h;
  x1[1] = xcent[1] + dz[1]*h;
  x1[2] = xcent[2] + dz[2]*h;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           x1, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Node on the Circle */
  x2[0] = xcent[0] + dx[0]*r;
  x2[1] = xcent[1] + dx[1]*r;
  x2[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           x2, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Line (point and direction) */
  data[0] = x2[0];
  data[1] = x2[1];
  data[2] = x2[2];
  data[3] = x2[0] - x1[0];
  data[4] = x2[1] - x1[1];
  data[5] = x2[2] - x1[2];
  status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &eline);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Line */
  tdata[0] = vmin;
  tdata[1] = 0;

  status = EG_makeTopology(context, eline, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the Degenerate Edges on the tip Node */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &enodes[0], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[2] = eedges[0]; /* repeat the line edge */

  /* make the Edge on the Circle */
  tdata[0] = 0.;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle, EDGE, ONENODE,
                           tdata, 1, &enodes[1], NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = vmin;  /* v == vmin VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0.;    /* v == 0 VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, econe, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, econe, FACE, SFORWARD,
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
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setConicalBody_dot( const double *xcent,     /* (in)  Center          */
                    const double *xcent_dot, /* (in)  Center velocity */
                    const double *xax,       /* (in)  x-axis          */
                    const double *xax_dot,   /* (in)  x-axis velocity */
                    const double *yax,       /* (in)  y-axis          */
                    const double *yax_dot,   /* (in)  y-axis velocity */
                    const double *zax,       /* (in)  z-axis          */
                    const double *zax_dot,   /* (in)  z-axis velocity */
                    const double angle,      /* (in)  angle           */
                    const double angle_dot,  /* (in)  angle velocity  */
                    const double r,          /* (in)  radius          */
                    const double r_dot,      /* (in)  radius velocity */
                    ego ebody )              /* (in/out) Circle body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double data[14], data_dot[14], dx[3], dx_dot[3], dy[3], dy_dot[3], dz[3], dz_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, vmin, vmin_dot, h, h_dot, x1[3], x1_dot[3], x2[3], x2_dot[3];
  double tdata[2], tdata_dot[2];
  ego    econe, ecircle, eline, *enodes, *eloops, *eedges, *efaces, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Sphere from the Face */
  status = EG_getTopology(efaces[0], &econe, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[3], &ecircle, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Line from the Edge */
  status = EG_getTopology(eedges[0], &eline, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Cone data and velocity */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = zax[0];   /* z-axis */
  data[10] = zax[1];
  data[11] = zax[2];
  data[12] = angle;    /* angle  */
  data[13] = r;        /* radius */

  data_dot[ 0] = xcent_dot[0]; /* center */
  data_dot[ 1] = xcent_dot[1];
  data_dot[ 2] = xcent_dot[2];
  data_dot[ 3] = xax_dot[0];   /* x-axis */
  data_dot[ 4] = xax_dot[1];
  data_dot[ 5] = xax_dot[2];
  data_dot[ 6] = yax_dot[0];   /* y-axis */
  data_dot[ 7] = yax_dot[1];
  data_dot[ 8] = yax_dot[2];
  data_dot[ 9] = zax_dot[0];   /* z-axis */
  data_dot[10] = zax_dot[1];
  data_dot[11] = zax_dot[2];
  data_dot[12] = angle_dot;    /* angle  */
  data_dot[13] = r_dot;        /* radius */

  status = EG_setGeometry_dot(econe, SURFACE, CONICAL, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(econe, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the cone */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];
  dz[0] = rvec[ 9];
  dz[1] = rvec[10];
  dz[2] = rvec[11];

  dx_dot[0] = rvec_dot[ 3];
  dx_dot[1] = rvec_dot[ 4];
  dx_dot[2] = rvec_dot[ 5];
  dy_dot[0] = rvec_dot[ 6];
  dy_dot[1] = rvec_dot[ 7];
  dy_dot[2] = rvec_dot[ 8];
  dz_dot[0] = rvec_dot[ 9];
  dz_dot[1] = rvec_dot[10];
  dz_dot[2] = rvec_dot[11];

  /* set the Circle curve velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];    /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;         /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];    /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  vmin = -r/sin(angle);
  vmin_dot = -r_dot/sin(angle) + angle_dot/tan(angle) * r/sin(angle) ;

  h = vmin*cos(angle);
  h_dot = vmin_dot*cos(angle) - vmin*sin(angle)*angle_dot;


  /* sensitivity of Node on the Tip*/
  x1[0] = xcent[0] + dz[0]*h;
  x1[1] = xcent[1] + dz[1]*h;
  x1[2] = xcent[2] + dz[2]*h;
  x1_dot[0] = xcent_dot[0] + dz_dot[0]*h + dz[0]*h_dot;
  x1_dot[1] = xcent_dot[1] + dz_dot[1]*h + dz[1]*h_dot;
  x1_dot[2] = xcent_dot[2] + dz_dot[2]*h + dz[2]*h_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, x1, x1_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* sensitivity of Node on the Circle*/
  x2[0] = xcent[0] + dx[0]*r;
  x2[1] = xcent[1] + dx[1]*r;
  x2[2] = xcent[2] + dx[2]*r;
  x2_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  x2_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  x2_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, x2, x2_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Line data and velocity */
  data[0] = x2[0];
  data[1] = x2[1];
  data[2] = x2[2];
  data[3] = x2[0] - x1[0];
  data[4] = x2[1] - x1[1];
  data[5] = x2[2] - x1[2];

  data_dot[0] = x2_dot[0];
  data_dot[1] = x2_dot[1];
  data_dot[2] = x2_dot[2];
  data_dot[3] = x2_dot[0] - x1_dot[0];
  data_dot[4] = x2_dot[1] - x1_dot[1];
  data_dot[5] = x2_dot[2] - x1_dot[2];

  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = vmin;
  tdata[1]     = 0;
  tdata_dot[0] = vmin_dot;
  tdata_dot[1] = 0;

  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
  /* set P-curve sensitivities */
  data_dot[0] = data_dot[1] = data_dot[2] = data_dot[3] = 0;

  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_setGeometry_dot(eedges[4+0], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = vmin;  /* v == vmin VMIN */
  data[2] = 1.;    data[3] = 0.  ;
                   data_dot[1] = vmin_dot;
  status = EG_setGeometry_dot(eedges[4+1], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  data_dot[1] = 0;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_setGeometry_dot(eedges[4+2], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] =  0.;   /* v == 0. VMAX */
  data[2] = 1.;    data[3] =  0.;
  status = EG_setGeometry_dot(eedges[4+3], PCURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
#endif

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingConical(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[14], x_dot[14], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *zax, *xcent_dot, *xax_dot, *yax_dot, *zax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;
  zax   = x+9;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;
  zax_dot   = x_dot+9;

  /* make the Circle body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  zax[0]   = 0.20; zax[1]   = 0.05; zax[2]   = 1.15;
//  xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
//  yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
//  zax[0]   = 0.; zax[1]   = 0.; zax[2]   = 1.;
  x[12] = 45.0*PI/180.;
  x[13] = 2.0;
  status = makeConicalBody(context, stack, xcent, xax, yax, zax, x[12], x[13], &ebody1);
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
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Conical Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Conical Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 14; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 14; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setConicalBody_dot(xcent, xcent_dot,
                                xax, xax_dot,
                                yax, yax_dot,
                                zax, zax_dot,
                                x[12], x_dot[12],
                                x[13], x_dot[13], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeConicalBody(context, stack, xcent, xax, yax, zax, x[12], x[13], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Conical", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Conical", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setConicalBody_dot(xcent, xcent_dot,
                              xax, xax_dot,
                              yax, yax_dot,
                              zax, zax_dot,
                              x[12], x_dot[12],
                              x[13], x_dot[13], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Conical", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Cylindrical                                                              */
/*                                                                           */
/*****************************************************************************/

int
makeCylindricalBody( ego context,         /* (in)  EGADS context    */
                     objStack *stack,     /* (in)  EGADS obj stack  */
                     const double *xcent, /* (in)  Center           */
                     const double *xax,   /* (in)  x-axis           */
                     const double *yax,   /* (in)  y-axis           */
                     const double *zax,   /* (in)  z-axis           */
                     const double r,      /* (in)  radius           */
                     ego *ebody )         /* (out) Cylindrical body */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[13], tdata[2], dx[3], dy[3], dz[3], x1[3], x2[3], *rvec=NULL;
  ego    ecylinder, ecircle[2], eedges[8], enodes[2], eloop, eface, eref;

  /* create the Cylindrical */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = zax[0];   /* z-axis */
  data[10] = zax[1];
  data[11] = zax[2];
  data[12] = r;        /* radius */
  status = EG_makeGeometry(context, SURFACE, CYLINDRICAL, NULL, NULL,
                           data, &ecylinder);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecylinder);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecylinder, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the cylinder */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];
  dz[0] = rvec[ 9];
  dz[1] = rvec[10];
  dz[2] = rvec[11];

  /* create the Circle curve for the base */
  data[0] = rvec[0]; /* center */
  data[1] = rvec[1];
  data[2] = rvec[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Circle curve for the top */
  data[0] = rvec[0] + dz[0]*r; /* center */
  data[1] = rvec[1] + dz[1]*r;
  data[2] = rvec[2] + dz[2]*r;
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Node on the Base */
  x1[0] = xcent[0] + dx[0]*r;
  x1[1] = xcent[1] + dx[1]*r;
  x1[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           x1, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Node on the Top */
  x2[0] = xcent[0] + dx[0]*r + dz[0]*r;
  x2[1] = xcent[1] + dx[1]*r + dz[1]*r;
  x2[2] = xcent[2] + dx[2]*r + dz[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           x2, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Edges */
  status =  makeLineEdge(context, stack, enodes[0], enodes[1], &eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle Base */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle[0], EDGE, ONENODE,
                           tdata, 1, &enodes[0], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[2] = eedges[0]; /* repeat the line edge */

  status = EG_makeTopology(context, ecircle[1], EDGE, ONENODE,
                           tdata, 1, &enodes[1], NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0;    /* v == 0 VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = r;    /* v == r VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, ecylinder, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, ecylinder, FACE, SFORWARD,
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
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setCylindricalBody_dot( const double *xcent,     /* (in)  Center          */
                        const double *xcent_dot, /* (in)  Center velocity */
                        const double *xax,       /* (in)  x-axis          */
                        const double *xax_dot,   /* (in)  x-axis velocity */
                        const double *yax,       /* (in)  y-axis          */
                        const double *yax_dot,   /* (in)  y-axis velocity */
                        const double *zax,       /* (in)  z-axis          */
                        const double *zax_dot,   /* (in)  z-axis velocity */
                        const double r,          /* (in)  radius          */
                        const double r_dot,      /* (in)  radius velocity */
                        ego ebody )              /* (in/out) body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double data[13], data_dot[13], dx[3], dx_dot[3], dy[3], dy_dot[3], dz[3], dz_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, x1[3], x1_dot[3], x2[3], x2_dot[3];
  double tdata[2], tdata_dot[2];
  ego    ecylinder, ecircle[2], *enodes, *eloops, *eedges, *efaces, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Sphere from the Face */
  status = EG_getTopology(efaces[0], &ecylinder, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle Base from the Edge */
  status = EG_getTopology(eedges[1], &ecircle[0], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle Top from the Edge */
  status = EG_getTopology(eedges[3], &ecircle[1], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Line Edge */
  status = EG_getTopology(eedges[0], &eref, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Cone data and velocity */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = zax[0];   /* z-axis */
  data[10] = zax[1];
  data[11] = zax[2];
  data[12] = r;        /* radius */

  data_dot[ 0] = xcent_dot[0]; /* center */
  data_dot[ 1] = xcent_dot[1];
  data_dot[ 2] = xcent_dot[2];
  data_dot[ 3] = xax_dot[0];   /* x-axis */
  data_dot[ 4] = xax_dot[1];
  data_dot[ 5] = xax_dot[2];
  data_dot[ 6] = yax_dot[0];   /* y-axis */
  data_dot[ 7] = yax_dot[1];
  data_dot[ 8] = yax_dot[2];
  data_dot[ 9] = zax_dot[0];   /* z-axis */
  data_dot[10] = zax_dot[1];
  data_dot[11] = zax_dot[2];
  data_dot[12] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecylinder, SURFACE, CYLINDRICAL, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(ecylinder, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the cone */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];
  dz[0] = rvec[ 9];
  dz[1] = rvec[10];
  dz[2] = rvec[11];

  dx_dot[0] = rvec_dot[ 3];
  dx_dot[1] = rvec_dot[ 4];
  dx_dot[2] = rvec_dot[ 5];
  dy_dot[0] = rvec_dot[ 6];
  dy_dot[1] = rvec_dot[ 7];
  dy_dot[2] = rvec_dot[ 8];
  dz_dot[0] = rvec_dot[ 9];
  dz_dot[1] = rvec_dot[10];
  dz_dot[2] = rvec_dot[11];

  /* set the Base Circle curve velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];    /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;         /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];    /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle[0], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Top Circle curve velocity */
  data[0] = xcent[0] + dz[0]*r; /* center */
  data[1] = xcent[1] + dz[1]*r;
  data[2] = xcent[2] + dz[2]*r;
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];    /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;         /* radius */

  data_dot[0] = xcent_dot[0] + dz_dot[0]*r + dz[0]*r_dot; /* center */
  data_dot[1] = xcent_dot[1] + dz_dot[1]*r + dz[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dz_dot[2]*r + dz[2]*r_dot;
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];    /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle[1], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* sensitivity of Node on the Base*/
  x1[0] = xcent[0] + dx[0]*r;
  x1[1] = xcent[1] + dx[1]*r;
  x1[2] = xcent[2] + dx[2]*r;
  x1_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  x1_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  x1_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, x1, x1_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* sensitivity of Node on the Circle*/
  x2[0] = xcent[0] + dx[0]*r + dz[0]*r;
  x2[1] = xcent[1] + dx[1]*r + dz[1]*r;
  x2[2] = xcent[2] + dx[2]*r + dz[2]*r;
  x2_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot + dz_dot[0]*r + dz[0]*r_dot;
  x2_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot + dz_dot[1]*r + dz[1]*r_dot;
  x2_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot + dz_dot[2]*r + dz[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, x2, x2_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Line Edge velocities */
  status = setLineEdge_dot( eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingCylindrical(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[13], x_dot[13], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *zax, *xcent_dot, *xax_dot, *yax_dot, *zax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;
  zax   = x+9;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;
  zax_dot   = x_dot+9;

  /* make the Cylindrical body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  zax[0]   = 0.20; zax[1]   = 0.05; zax[2]   = 1.15;
//  xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
//  yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
//  zax[0]   = 0.; zax[1]   = 0.; zax[2]   = 1.;
  x[12] = 1.0;
  status = makeCylindricalBody(context, stack, xcent, xax, yax, zax, x[12], &ebody1);
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
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Cylindrical Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Cylindrical Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 13; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 13; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setCylindricalBody_dot(xcent, xcent_dot,
                                    xax, xax_dot,
                                    yax, yax_dot,
                                    zax, zax_dot,
                                    x[12], x_dot[12], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeCylindricalBody(context, stack, xcent, xax, yax, zax, x[12], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Cylindrical", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Cylindrical", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setCylindricalBody_dot(xcent, xcent_dot,
                                  xax, xax_dot,
                                  yax, yax_dot,
                                  zax, zax_dot,
                                  x[12], x_dot[12], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Cylindrical", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Toroidal                                                                 */
/*                                                                           */
/*****************************************************************************/

int
makeToroidalBody( ego context,         /* (in)  EGADS context   */
                  objStack *stack,     /* (in)  EGADS obj stack */
                  const double *xcent, /* (in)  Center          */
                  const double *xax,   /* (in)  x-axis          */
                  const double *yax,   /* (in)  y-axis          */
                  const double *zax,   /* (in)  z-axis          */
                  const double majr,   /* (in)  major radius    */
                  const double minr,   /* (in)  minor radius    */
                  ego *ebody )         /* (out) Toroidal body   */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[14], tdata[2], dx[3], dy[3], dz[3], R, *rvec=NULL;
  ego    etorus, ecircle[3], eedges[8], enodes[2], eloop, eface, eref;

  R = minr + majr;

  /* create the Toroidal */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = zax[0];   /* z-axis */
  data[10] = zax[1];
  data[11] = zax[2];
  data[12] = majr;     /* major radius */
  data[13] = minr;     /* minor radius */
  status = EG_makeGeometry(context, SURFACE, TOROIDAL, NULL, NULL,
                           data, &etorus);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, etorus);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(etorus, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the torus */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];
  dz[0] = rvec[ 9];
  dz[1] = rvec[10];
  dz[2] = rvec[11];

  /* create Circle curve */
  data[0] = xcent[0] + dx[0]*majr; /* center */
  data[1] = xcent[1] + dx[1]*majr;
  data[2] = xcent[2] + dx[2]*majr;
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dz[0];   /* y-axis */
  data[7] = dz[1];
  data[8] = dz[2];
  data[9] = minr;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the outer Circle curve */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = R;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Circle curve */
  data[0] = xcent[0] - dx[0]*majr; /* center */
  data[1] = xcent[1] - dx[1]*majr;
  data[2] = xcent[2] - dx[2]*majr;
  data[3] = -dx[0];   /* x-axis */
  data[4] = -dx[1];
  data[5] = -dx[2];
  data[6] = dz[0];   /* y-axis */
  data[7] = dz[1];
  data[8] = dz[2];
  data[9] = minr;    /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[2]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Nodes */
  data[0] = xcent[0] + dx[0]*R;
  data[1] = xcent[1] + dx[1]*R;
  data[2] = xcent[2] + dx[2]*R;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] - dx[0]*R;
  data[1] = xcent[1] - dx[1]*R;
  data[2] = xcent[2] - dx[2]*R;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edges on the Circles */
  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle[0], EDGE, ONENODE,
                           tdata, 1, &enodes[0], NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = 0;
  tdata[1] = PI;

  status = EG_makeTopology(context, ecircle[1], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle[2], EDGE, ONENODE,
                           tdata, 1, &enodes[1], NULL, &eedges[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[3] = eedges[1]; /* repeat the edge */


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0.;    /* v == 0 VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = PI;    data[1] = 0.;    /* u == PI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = TWOPI; /* v == TWOPI VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, etorus, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, etorus, FACE, SFORWARD,
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
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setToroidalBody_dot( const double *xcent,     /* (in)  Center                  */
                     const double *xcent_dot, /* (in)  Center velocity         */
                     const double *xax,       /* (in)  x-axis                  */
                     const double *xax_dot,   /* (in)  x-axis velocity         */
                     const double *yax,       /* (in)  y-axis                  */
                     const double *yax_dot,   /* (in)  y-axis velocity         */
                     const double *zax,       /* (in)  z-axis                  */
                     const double *zax_dot,   /* (in)  z-axis velocity         */
                     const double majr,       /* (in)  major radius            */
                     const double majr_dot,   /* (in)  major radius velocity   */
                     const double minr,       /* (in)  minor radius            */
                     const double minr_dot,   /* (in)  minor radius velocity   */
                     ego ebody )              /* (in/out) body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses;
  double data[14], data_dot[14], dx[3], dx_dot[3], dy[3], dy_dot[3], dz[3], dz_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, R, R_dot, tdata[2], tdata_dot[2];
  ego    etorus, ecircle[3], *enodes, *eloops, *eedges, *efaces, eref;

  R     = minr     + majr;
  R_dot = minr_dot + majr_dot;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Sphere from the Face */
  status = EG_getTopology(efaces[0], &etorus, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[0], &ecircle[0], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[2], &ecircle[2], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle and Nodes from the Edge */
  status = EG_getTopology(eedges[1], &ecircle[1], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[2], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0]     = 0;
  tdata[1]     = PI;

  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Cone data and velocity */
  data[ 0] = xcent[0]; /* center */
  data[ 1] = xcent[1];
  data[ 2] = xcent[2];
  data[ 3] = xax[0];   /* x-axis */
  data[ 4] = xax[1];
  data[ 5] = xax[2];
  data[ 6] = yax[0];   /* y-axis */
  data[ 7] = yax[1];
  data[ 8] = yax[2];
  data[ 9] = zax[0];   /* z-axis */
  data[10] = zax[1];
  data[11] = zax[2];
  data[12] = majr;     /* major radius */
  data[13] = minr;     /* minor radius */

  data_dot[ 0] = xcent_dot[0]; /* center */
  data_dot[ 1] = xcent_dot[1];
  data_dot[ 2] = xcent_dot[2];
  data_dot[ 3] = xax_dot[0];   /* x-axis */
  data_dot[ 4] = xax_dot[1];
  data_dot[ 5] = xax_dot[2];
  data_dot[ 6] = yax_dot[0];   /* y-axis */
  data_dot[ 7] = yax_dot[1];
  data_dot[ 8] = yax_dot[2];
  data_dot[ 9] = zax_dot[0];   /* z-axis */
  data_dot[10] = zax_dot[1];
  data_dot[11] = zax_dot[2];
  data_dot[12] = majr_dot;     /* major radius */
  data_dot[13] = minr_dot;     /* minor radius */

  status = EG_setGeometry_dot(etorus, SURFACE, TOROIDAL, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(etorus, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the torus */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];
  dz[0] = rvec[ 9];
  dz[1] = rvec[10];
  dz[2] = rvec[11];

  dx_dot[0] = rvec_dot[ 3];
  dx_dot[1] = rvec_dot[ 4];
  dx_dot[2] = rvec_dot[ 5];
  dy_dot[0] = rvec_dot[ 6];
  dy_dot[1] = rvec_dot[ 7];
  dy_dot[2] = rvec_dot[ 8];
  dz_dot[0] = rvec_dot[ 9];
  dz_dot[1] = rvec_dot[10];
  dz_dot[2] = rvec_dot[11];

  /* set the Circle curve velocity */
  data[0] = xcent[0] + dx[0]*majr; /* center */
  data[1] = xcent[1] + dx[1]*majr;
  data[2] = xcent[2] + dx[2]*majr;
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dz[0];    /* y-axis */
  data[7] = dz[1];
  data[8] = dz[2];
  data[9] = minr;     /* radius */

  data_dot[0] = xcent_dot[0] + dx_dot[0]*majr + dx[0]*majr_dot; /* center */
  data_dot[1] = xcent_dot[1] + dx_dot[1]*majr + dx[1]*majr_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*majr + dx[2]*majr_dot;
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dz_dot[0];    /* y-axis */
  data_dot[7] = dz_dot[1];
  data_dot[8] = dz_dot[2];
  data_dot[9] = minr_dot;     /* radius */

  status = EG_setGeometry_dot(ecircle[0], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the outer Circle curve velocity */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];    /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = R;         /* radius */

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];    /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = R_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle[1], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Circle curve velocity */
  data[0] = xcent[0] - dx[0]*majr; /* center */
  data[1] = xcent[1] - dx[1]*majr;
  data[2] = xcent[2] - dx[2]*majr;
  data[3] = -dx[0];    /* x-axis */
  data[4] = -dx[1];
  data[5] = -dx[2];
  data[6] = dz[0];    /* y-axis */
  data[7] = dz[1];
  data[8] = dz[2];
  data[9] = minr;     /* radius */

  data_dot[0] = xcent_dot[0] - dx_dot[0]*majr - dx[0]*majr_dot; /* center */
  data_dot[1] = xcent_dot[1] - dx_dot[1]*majr - dx[1]*majr_dot;
  data_dot[2] = xcent_dot[2] - dx_dot[2]*majr - dx[2]*majr_dot;
  data_dot[3] = -dx_dot[0];    /* x-axis */
  data_dot[4] = -dx_dot[1];
  data_dot[5] = -dx_dot[2];
  data_dot[6] = dz_dot[0];    /* y-axis */
  data_dot[7] = dz_dot[1];
  data_dot[8] = dz_dot[2];
  data_dot[9] = minr_dot;     /* radius */

  status = EG_setGeometry_dot(ecircle[2], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* sensitivity of Nodes */
  data[0] = xcent[0] + dx[0]*R;
  data[1] = xcent[1] + dx[1]*R;
  data[2] = xcent[2] + dx[2]*R;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*R + dx[0]*R_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*R + dx[1]*R_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*R + dx[2]*R_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] - dx[0]*R;
  data[1] = xcent[1] - dx[1]*R;
  data[2] = xcent[2] - dx[2]*R;
  data_dot[0] = xcent_dot[0] - dx_dot[0]*R - dx[0]*R_dot;
  data_dot[1] = xcent_dot[1] - dx_dot[1]*R - dx[1]*R_dot;
  data_dot[2] = xcent_dot[2] - dx_dot[2]*R - dx[2]*R_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingToroidal(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[14], x_dot[14], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *zax, *xcent_dot, *xax_dot, *yax_dot, *zax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;
  zax   = x+9;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;
  zax_dot   = x_dot+9;

  /* make the Toroidal body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  zax[0]   = 0.20; zax[1]   = 0.05; zax[2]   = 1.15;
//  xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
//  yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
//  zax[0]   = 0.; zax[1]   = 0.; zax[2]   = 1.;
  x[12] = 2.5;
  x[13] = 1.0;
  status = makeToroidalBody(context, stack, xcent, xax, yax, zax, x[12], x[13], &ebody1);
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
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Toroidal Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Toroidal Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 14; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 14; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setToroidalBody_dot(xcent, xcent_dot,
                                 xax, xax_dot,
                                 yax, yax_dot,
                                 zax, zax_dot,
                                 x[12], x_dot[12],
                                 x[13], x_dot[13], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeToroidalBody(context, stack, xcent, xax, yax, zax, x[12], x[13], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Toroidal", 5e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Toroidal", 5e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setToroidalBody_dot(xcent, xcent_dot,
                               xax, xax_dot,
                               yax, yax_dot,
                               zax, zax_dot,
                               x[12], x_dot[12],
                               x[13], x_dot[13], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Toroidal", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Revolution                                                               */
/*                                                                           */
/*****************************************************************************/

int
makeRevolutionBody( ego context,         /* (in)  EGADS context    */
                    objStack *stack,     /* (in)  EGADS obj stack  */
                    const double *xcent, /* (in)  Center           */
                    const double *xax,   /* (in)  x-axis           */
                    const double *yax,   /* (in)  y-axis           */
                    const double r,      /* (in)  radius           */
                    ego *ebody )         /* (out) Revolved body    */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, nnode, oclass, mtype, *ivec=NULL, *sens=NULL;
  double data[14], tdata[2], dx[3], dy[3], dz[3], *rvec=NULL;
  ego    esurf, ecircle[2], eline, eedges[8], enodes[2], eloop, eface, eref, *echld;

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
                           data, &ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecircle[0], &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];

  /* compute the axis of the revolution */
  CROSS(dz, dx, dy);

  /* create Circle curve */
  data[0] = xcent[0] + dz[0]*2*r; /* center */
  data[1] = xcent[1] + dz[1]*2*r;
  data[2] = xcent[2] + dz[2]*2*r;
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Nodes */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dz[0]*2*r + dx[0]*r;
  data[1] = xcent[1] + dz[1]*2*r + dx[1]*r;
  data[2] = xcent[2] + dz[2]*2*r + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the Edge on a Line */
  status =  makeLineEdge(context, stack, enodes[0], enodes[1], &eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Line from the Edge */
  status = EG_getTopology(eedges[0], &eline, &oclass, &mtype,
                          data, &nnode, &echld, &sens);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Revolution surface */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dz[0]; /* direction */
  data[4] = dz[1];
  data[5] = dz[2];
  status = EG_makeGeometry(context, SURFACE, REVOLUTION, eline, NULL,
                           data, &esurf);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, esurf);
  if (status != EGADS_SUCCESS) goto cleanup;


  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle[0], EDGE, ONENODE,
                           tdata, 1, &enodes[0], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[2] = eedges[0]; /* repeat the line edge */

  status = EG_makeTopology(context, ecircle[1], EDGE, ONENODE,
                           tdata, 1, &enodes[1], NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0.;    /* v == 0 VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 2*r;   /* v == 2*r VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, esurf, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, esurf, FACE, SFORWARD,
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
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setRevolutionBody_dot( const double *xcent,     /* (in)  Center                    */
                       const double *xcent_dot, /* (in)  Center velocity           */
                       const double *xax,       /* (in)  x-axis                    */
                       const double *xax_dot,   /* (in)  x-axis velocity           */
                       const double *yax,       /* (in)  y-axis                    */
                       const double *yax_dot,   /* (in)  y-axis velocity           */
                       const double r,          /* (in)  radius                    */
                       const double r_dot,      /* (in)  radius velocity           */
                       ego ebody )              /* (in/out) body with velocities   */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses, *ivec=NULL;
  double data[14], data_dot[14], dx[3], dx_dot[3], dy[3], dy_dot[3], dz[3], dz_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, tdata[2], tdata_dot[2];
  ego    esurf, ecircle[2], eline, *enodes, *eloops, *eedges, *efaces, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Surface from the Face */
  status = EG_getTopology(efaces[0], &esurf, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[1], &ecircle[0], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[3], &ecircle[1], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Line Edge */
  status = EG_getTopology(eedges[0], &eline, &oclass, &mtype,
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

  status = EG_setGeometry_dot(ecircle[0], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry_dot(ecircle[0], &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];

  dx_dot[0] = rvec_dot[ 3];
  dx_dot[1] = rvec_dot[ 4];
  dx_dot[2] = rvec_dot[ 5];
  dy_dot[0] = rvec_dot[ 6];
  dy_dot[1] = rvec_dot[ 7];
  dy_dot[2] = rvec_dot[ 8];

  /* compute the axis of the sphere between the poles */
  CROSS(dz, dx, dy);
  CROSS_DOT(dz_dot, dx, dx_dot, dy, dy_dot);

  /* set the Circle curve velocity */
  data[0] = xcent[0] + dz[0]*2*r; /* center */
  data[1] = xcent[1] + dz[1]*2*r;
  data[2] = xcent[2] + dz[2]*2*r;
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];    /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */

  data_dot[0] = xcent_dot[0] + dz_dot[0]*2*r + dz[0]*2*r_dot; /* center */
  data_dot[1] = xcent_dot[1] + dz_dot[1]*2*r + dz[1]*2*r_dot;
  data_dot[2] = xcent_dot[2] + dz_dot[2]*2*r + dz[2]*2*r_dot;
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];    /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle[1], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* sensitivity of Nodes */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dz[0]*2*r + dx[0]*r;
  data[1] = xcent[1] + dz[1]*2*r + dx[1]*r;
  data[2] = xcent[2] + dz[2]*2*r + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dz_dot[0]*2*r + dz[0]*2*r_dot + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dz_dot[1]*2*r + dz[1]*2*r_dot + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dz_dot[2]*2*r + dz[2]*2*r_dot + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Line Edge velocities */
  status = setLineEdge_dot( eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the velocity on the surface */
  data[0] = xcent[0]; /* center */
  data[1] = xcent[1];
  data[2] = xcent[2];
  data[3] = dz[0]; /* direction */
  data[4] = dz[1];
  data[5] = dz[2];

  data_dot[0] = xcent_dot[0]; /* center */
  data_dot[1] = xcent_dot[1];
  data_dot[2] = xcent_dot[2];
  data_dot[3] = dz_dot[0]; /* direction */
  data_dot[4] = dz_dot[1];
  data_dot[5] = dz_dot[2];

  status = EG_setGeometry_dot(esurf, SURFACE, REVOLUTION, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the revolution reference geometry velocity */
  EG_free(rvec); rvec = NULL;
  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_free(rvec); rvec = NULL;
  EG_free(rvec_dot); rvec_dot = NULL;
  status = EG_getGeometry_dot(eline, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(eref, CURVE, LINE, NULL, rvec, rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(ivec);
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingRevolution(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* make the Revolution body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
//  xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
//  yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
  x[9] = 1.0;
  status = makeRevolutionBody(context, stack, xcent, xax, yax, x[9], &ebody1);
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
  params[1] =  0.2;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Revolution Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Revolution Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setRevolutionBody_dot(xcent, xcent_dot,
                                   xax, xax_dot,
                                   yax, yax_dot,
                                   x[9], x_dot[9], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeRevolutionBody(context, stack, xcent, xax, yax, x[9], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Revolution", 5e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Revolution", 5e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setRevolutionBody_dot(xcent, xcent_dot,
                                 xax, xax_dot,
                                 yax, yax_dot,
                                 x[9], x_dot[9], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Revolution", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Extrusion                                                                */
/*                                                                           */
/*****************************************************************************/

int
makeExtrusionBody( ego context,         /* (in)  EGADS context    */
                   objStack *stack,     /* (in)  EGADS obj stack  */
                   const double *xcent, /* (in)  Center           */
                   const double *xax,   /* (in)  x-axis           */
                   const double *yax,   /* (in)  y-axis           */
                   const double r,      /* (in)  radius           */
                   const double *vec,   /* (in)  extrusion vector */
                   ego *ebody )         /* (out) Extruded body    */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[14], tdata[2], dx[3], dy[3], *rvec=NULL, vmag;
  ego    esurf, ecircle[2], eedges[8], enodes[2], eloop, eface, eref;

  vmag = sqrt(DOT(vec, vec));

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
                           data, &ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(ecircle[0], &oclass, &mtype, &eref, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];

  /* create Circle extruded curve */
  data[0] = xcent[0] + vec[0]; /* center */
  data[1] = xcent[1] + vec[1];
  data[2] = xcent[2] + vec[2];
  data[3] = dx[0];   /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];   /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */
  status = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           data, &ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecircle[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Extrusion surface */
  data[0] = vec[0]; /* direction */
  data[1] = vec[1];
  data[2] = vec[2];
  status = EG_makeGeometry(context, SURFACE, EXTRUSION, ecircle[0], NULL,
                           data, &esurf);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, esurf);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Nodes */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + vec[0] + dx[0]*r;
  data[1] = xcent[1] + vec[1] + dx[1]*r;
  data[2] = xcent[2] + vec[2] + dx[2]*r;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the Edges on the Circles */
  status =  makeLineEdge(context, stack, enodes[0], enodes[1], &eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  tdata[0] = 0;
  tdata[1] = TWOPI;

  status = EG_makeTopology(context, ecircle[0], EDGE, ONENODE,
                           tdata, 1, &enodes[0], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[2] = eedges[0]; /* repeat the line edge */

  status = EG_makeTopology(context, ecircle[1], EDGE, ONENODE,
                           tdata, 1, &enodes[1], NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0.;    /* v == 0 VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = vmag;  /* v == vmag VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, esurf, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, esurf, FACE, SFORWARD,
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
  EG_free(ivec); ivec = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setExtrusionBody_dot( const double *xcent,     /* (in)  Center                    */
                      const double *xcent_dot, /* (in)  Center velocity           */
                      const double *xax,       /* (in)  x-axis                    */
                      const double *xax_dot,   /* (in)  x-axis velocity           */
                      const double *yax,       /* (in)  y-axis                    */
                      const double *yax_dot,   /* (in)  y-axis velocity           */
                      const double r,          /* (in)  radius                    */
                      const double r_dot,      /* (in)  radius velocity           */
                      const double *vec,       /* (in)  extrusion vector          */
                      const double *vec_dot,   /* (in)  extrusion vector velocity */
                      ego ebody )              /* (in/out) body with velocities   */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses, *ivec=NULL;
  double data[14], data_dot[14], dx[3], dx_dot[3], dy[3], dy_dot[3];
  double *rvec=NULL, *rvec_dot=NULL, tdata[2], tdata_dot[2];
  ego    esurf, ecircle[2], eline, *enodes, *eloops, *eedges, *efaces, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Surface from the Face */
  status = EG_getTopology(efaces[0], &esurf, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edges from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[1], &ecircle[0], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Circle from the Edge */
  status = EG_getTopology(eedges[3], &ecircle[1], &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Line Edge */
  status = EG_getTopology(eedges[0], &eline, &oclass, &mtype,
                          data, &nnode, &enodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Extruded surface velocity */
  data[0] = vec[0]; /* direction */
  data[1] = vec[1];
  data[2] = vec[2];
  data_dot[0] = vec_dot[0];
  data_dot[1] = vec_dot[1];
  data_dot[2] = vec_dot[2];
  status = EG_setGeometry_dot(esurf, SURFACE, EXTRUSION, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
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

  status = EG_setGeometry_dot(ecircle[0], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the extrude reference geometry velocity */
  status = EG_setGeometry_dot(eref, CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_free(rvec); rvec = NULL;
  status = EG_getGeometry_dot(ecircle[0], &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the axes for the circle */
  dx[0] = rvec[ 3];
  dx[1] = rvec[ 4];
  dx[2] = rvec[ 5];
  dy[0] = rvec[ 6];
  dy[1] = rvec[ 7];
  dy[2] = rvec[ 8];

  dx_dot[0] = rvec_dot[ 3];
  dx_dot[1] = rvec_dot[ 4];
  dx_dot[2] = rvec_dot[ 5];
  dy_dot[0] = rvec_dot[ 6];
  dy_dot[1] = rvec_dot[ 7];
  dy_dot[2] = rvec_dot[ 8];

  /* set the Circle curve velocity */
  data[0] = xcent[0] + vec[0]; /* center */
  data[1] = xcent[1] + vec[1];
  data[2] = xcent[2] + vec[2];
  data[3] = dx[0];    /* x-axis */
  data[4] = dx[1];
  data[5] = dx[2];
  data[6] = dy[0];    /* y-axis */
  data[7] = dy[1];
  data[8] = dy[2];
  data[9] = r;       /* radius */

  data_dot[0] = xcent_dot[0] + vec_dot[0]; /* center */
  data_dot[1] = xcent_dot[1] + vec_dot[1];
  data_dot[2] = xcent_dot[2] + vec_dot[2];
  data_dot[3] = dx_dot[0];    /* x-axis */
  data_dot[4] = dx_dot[1];
  data_dot[5] = dx_dot[2];
  data_dot[6] = dy_dot[0];    /* y-axis */
  data_dot[7] = dy_dot[1];
  data_dot[8] = dy_dot[2];
  data_dot[9] = r_dot;        /* radius */

  status = EG_setGeometry_dot(ecircle[1], CURVE, CIRCLE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;



  /* sensitivity of Nodes */
  data[0] = xcent[0] + dx[0]*r;
  data[1] = xcent[1] + dx[1]*r;
  data[2] = xcent[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + vec[0] + dx[0]*r;
  data[1] = xcent[1] + vec[1] + dx[1]*r;
  data[2] = xcent[2] + vec[2] + dx[2]*r;
  data_dot[0] = xcent_dot[0] + vec_dot[0] + dx_dot[0]*r + dx[0]*r_dot;
  data_dot[1] = xcent_dot[1] + vec_dot[1] + dx_dot[1]*r + dx[1]*r_dot;
  data_dot[2] = xcent_dot[2] + vec_dot[2] + dx_dot[2]*r + dx[2]*r_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set Line Edge velocities */
  status = setLineEdge_dot( eedges[0] );
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = TWOPI;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(ivec);
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingExtrusion(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[13], x_dot[13], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *vec, *xcent_dot, *xax_dot, *yax_dot, *vec_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;
  vec   = x+10;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;
  vec_dot   = x_dot+10;

  /* make the Extruded body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
//  xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
//  yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
  x[9] = 1.0;
  vec[0] = 0.0; vec[1] = 5.0; vec[2] = 5.0;
  status = makeExtrusionBody(context, stack, xcent, xax, yax, x[9], vec, &ebody1);
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
  params[1] =  0.2;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Extrusion Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Extrusion Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 13; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 13; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setExtrusionBody_dot(xcent, xcent_dot,
                                  xax, xax_dot,
                                  yax, yax_dot,
                                  x[9], x_dot[9],
                                  vec, vec_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    x[iparam] += dtime;
    status = makeExtrusionBody(context, stack, xcent, xax, yax, x[9], vec, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Extrusion", 5e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Extrusion", 5e-7, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setExtrusionBody_dot(xcent, xcent_dot,
                                xax, xax_dot,
                                yax, yax_dot,
                                x[9], x_dot[9],
                                vec, vec_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Extrusion", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}

/*****************************************************************************/
/*                                                                           */
/*  Bezier Surface                                                           */
/*                                                                           */
/*****************************************************************************/

int
makeBezierSurfaceBody( ego context,        /* (in)  EGADS context         */
                       objStack *stack,    /* (in)  EGADS obj stack       */
                       const int nCPu,     /* (in)  number of points in u */
                       const int nCPv,     /* (in)  number of points in v */
                       const double *pts,  /* (in)  coordinates           */
                       ego *ebody )        /* (out) Bezier Surface body   */
{
  int    status = EGADS_SUCCESS;
  int    i, j, senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, header[5];
  double data[4], tdata[2], *lpts=NULL;
  ego    esurf, ecurves[4], eedges[8], enodes[4], eloop, eface, nodes[2];

  lpts = (double*)EG_alloc(3*MAX(nCPu,nCPv)*sizeof(double));

  header[0] = 0;
  header[1] = nCPu-1;
  header[2] = nCPu;
  header[3] = nCPv-1;
  header[4] = nCPv;

  /* create the Bezier surface */
  status = EG_makeGeometry(context, SURFACE, BEZIER, NULL,
                           header, pts, &esurf);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, esurf);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create Nodes for the Edges */
  i = 0;
  j = 0;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = 0;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = nCPv-1;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = 0;
  j = nCPv-1;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the curves for the Edges */
  header[0] = 0;
  header[1] = nCPv-1;
  header[2] = nCPv;

  i = 0;
  for (j = 0; j < nCPv; j++) {
    lpts[3*j+0] = pts[3*(i+j*nCPu)+0];
    lpts[3*j+1] = pts[3*(i+j*nCPu)+1];
    lpts[3*j+2] = pts[3*(i+j*nCPu)+2];
  }
  status = EG_makeGeometry(context, CURVE, BEZIER, NULL, header,
                           lpts, &ecurves[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurves[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  header[0] = 0;
  header[1] = nCPu-1;
  header[2] = nCPu;

  j = 0;
  for (i = 0; i < nCPu; i++) {
    lpts[3*i+0] = pts[3*(i+j*nCPu)+0];
    lpts[3*i+1] = pts[3*(i+j*nCPu)+1];
    lpts[3*i+2] = pts[3*(i+j*nCPu)+2];
  }
  status = EG_makeGeometry(context, CURVE, BEZIER, NULL, header,
                           lpts, &ecurves[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurves[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  header[0] = 0;
  header[1] = nCPv-1;
  header[2] = nCPv;

  i = nCPu-1;
  for (j = 0; j < nCPv; j++) {
    lpts[3*j+0] = pts[3*(i+j*nCPu)+0];
    lpts[3*j+1] = pts[3*(i+j*nCPu)+1];
    lpts[3*j+2] = pts[3*(i+j*nCPu)+2];
  }
  status = EG_makeGeometry(context, CURVE, BEZIER, NULL, header,
                           lpts, &ecurves[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurves[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  header[0] = 0;
  header[1] = nCPu-1;
  header[2] = nCPu;

  j = nCPv-1;
  for (i = 0; i < nCPu; i++) {
    lpts[3*i+0] = pts[3*(i+j*nCPu)+0];
    lpts[3*i+1] = pts[3*(i+j*nCPu)+1];
    lpts[3*i+2] = pts[3*(i+j*nCPu)+2];
  }
  status = EG_makeGeometry(context, CURVE, BEZIER, NULL, header,
                           lpts, &ecurves[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, ecurves[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Curves */
  tdata[0] = 0;
  tdata[1] = 1;

  nodes[0] = enodes[0];
  nodes[1] = enodes[3];
  status = EG_makeTopology(context, ecurves[0], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  nodes[0] = enodes[0];
  nodes[1] = enodes[1];
  status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  nodes[0] = enodes[1];
  nodes[1] = enodes[2];
  status = EG_makeTopology(context, ecurves[2], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  nodes[0] = enodes[3];
  nodes[1] = enodes[2];
  status = EG_makeTopology(context, ecurves[3], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0.;    /* v == 0 VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 1;     data[1] = 0.;    /* u == 1 UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 1.;    /* v == 1 VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, esurf, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, esurf, FACE, SFORWARD,
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
  EG_free(lpts); lpts = NULL;

  return status;
}


int
setBezierSurfaceBody_dot( const int nCPu,        /* (in)  number of points in u   */
                          const int nCPv,        /* (in)  number of points in v   */
                          const double *pts,     /* (in)  coordinates             */
                          const double *pts_dot, /* (in)  velocity of coordinates */
                          ego ebody )            /* (in/out) body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    i, j, nnode, nedge, nloop, nface, oclass, mtype, *senses, header[5];
  double data[18], data_dot[18], tdata[2], tdata_dot[2];
  double *lpts=NULL, *lpts_dot=NULL;
  ego    esurf, ecurves[4], enodes[4], *nodes, *eloops, *eedges, *efaces, eref;


  lpts     = (double*)EG_alloc(3*MAX(nCPu,nCPv)*sizeof(double));
  lpts_dot = (double*)EG_alloc(3*MAX(nCPu,nCPv)*sizeof(double));

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Sphere from the Face */
  status = EG_getTopology(efaces[0], &esurf, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Curves and Nodes from the Edges */
  status = EG_getTopology(eedges[0], &ecurves[0], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = nodes[0];
  enodes[3] = nodes[1];

  status = EG_getTopology(eedges[1], &ecurves[1], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = nodes[0];
  enodes[1] = nodes[1];

  status = EG_getTopology(eedges[2], &ecurves[2], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[1] = nodes[0];
  enodes[2] = nodes[1];

  status = EG_getTopology(eedges[3], &ecurves[3], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[3] = nodes[0];
  enodes[2] = nodes[1];


  header[0] = 0;
  header[1] = nCPu-1;
  header[2] = nCPu;
  header[3] = nCPv-1;
  header[4] = nCPv;

  /* set the surface sensitivity */
  status = EG_setGeometry_dot(esurf, SURFACE, BEZIER, header, pts, pts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Node velocities */
  i = 0;
  j = 0;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = 0;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = nCPv-1;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = 0;
  j = nCPv-1;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[3], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Curve velocities */
  header[0] = 0;
  header[1] = nCPv-1;
  header[2] = nCPv;

  i = 0;
  for (j = 0; j < nCPv; j++) {
    lpts[3*j+0]     = pts[3*(i+j*nCPu)+0];
    lpts[3*j+1]     = pts[3*(i+j*nCPu)+1];
    lpts[3*j+2]     = pts[3*(i+j*nCPu)+2];
    lpts_dot[3*j+0] = pts_dot[3*(i+j*nCPu)+0];
    lpts_dot[3*j+1] = pts_dot[3*(i+j*nCPu)+1];
    lpts_dot[3*j+2] = pts_dot[3*(i+j*nCPu)+2];
  }
  status = EG_setGeometry_dot(ecurves[0], CURVE, BEZIER, header, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  header[0] = 0;
  header[1] = nCPu-1;
  header[2] = nCPu;

  j = 0;
  for (i = 0; i < nCPu; i++) {
    lpts[3*i+0]     = pts[3*(i+j*nCPu)+0];
    lpts[3*i+1]     = pts[3*(i+j*nCPu)+1];
    lpts[3*i+2]     = pts[3*(i+j*nCPu)+2];
    lpts_dot[3*i+0] = pts_dot[3*(i+j*nCPu)+0];
    lpts_dot[3*i+1] = pts_dot[3*(i+j*nCPu)+1];
    lpts_dot[3*i+2] = pts_dot[3*(i+j*nCPu)+2];
  }
  status = EG_setGeometry_dot(ecurves[1], CURVE, BEZIER, header, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  header[0] = 0;
  header[1] = nCPv-1;
  header[2] = nCPv;

  i = nCPu-1;
  for (j = 0; j < nCPv; j++) {
    lpts[3*j+0]     = pts[3*(i+j*nCPu)+0];
    lpts[3*j+1]     = pts[3*(i+j*nCPu)+1];
    lpts[3*j+2]     = pts[3*(i+j*nCPu)+2];
    lpts_dot[3*j+0] = pts_dot[3*(i+j*nCPu)+0];
    lpts_dot[3*j+1] = pts_dot[3*(i+j*nCPu)+1];
    lpts_dot[3*j+2] = pts_dot[3*(i+j*nCPu)+2];
  }
  status = EG_setGeometry_dot(ecurves[2], CURVE, BEZIER, header, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  header[0] = 0;
  header[1] = nCPu-1;
  header[2] = nCPu;

  j = nCPv-1;
  for (i = 0; i < nCPu; i++) {
    lpts[3*i+0]     = pts[3*(i+j*nCPu)+0];
    lpts[3*i+1]     = pts[3*(i+j*nCPu)+1];
    lpts[3*i+2]     = pts[3*(i+j*nCPu)+2];
    lpts_dot[3*i+0] = pts_dot[3*(i+j*nCPu)+0];
    lpts_dot[3*i+1] = pts_dot[3*(i+j*nCPu)+1];
    lpts_dot[3*i+2] = pts_dot[3*(i+j*nCPu)+2];
  }
  status = EG_setGeometry_dot(ecurves[3], CURVE, BEZIER, header, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = 1;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[2], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(lpts);
  EG_free(lpts_dot);

  return status;
}


int
pingBezierSurface(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    i, np1, nt1, iedge, nedge, iface, nface;
  double params[3], dtime = 1e-8;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  const int nCPu = 4;
  const int nCPv = 4;
  double pts[3*4*4] = {0.00, 0.00, 0.00,
                       1.00, 0.00, 0.10,
                       1.50, 1.00, 0.70,
                       0.25, 0.75, 0.60,

                       0.00, 0.00, 1.00,
                       1.00, 0.00, 1.10,
                       1.50, 1.00, 1.70,
                       0.25, 0.75, 1.60,

                       0.00, 0.00, 2.00,
                       1.00, 0.00, 2.10,
                       1.50, 1.00, 2.70,
                       0.25, 0.75, 2.60,

                       0.00, 0.00, 3.00,
                       1.00, 0.00, 3.10,
                       1.50, 1.00, 3.70,
                       0.25, 0.75, 3.60};
  double pts_dot[3*4*4];

  status = makeBezierSurfaceBody(context, stack, nCPu, nCPv, pts, &ebody1);
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
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Bezier Surface Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Bezier Surface Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (i = 0; i < 3*nCPu*nCPv; i++) pts_dot[i] = 0;

  for (i = 0; i < 3*nCPu*nCPv; i++) {

    /* set the velocity of the original body */
    pts_dot[i] = 1.0;
    status = setBezierSurfaceBody_dot(nCPu, nCPv, pts, pts_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts_dot[i] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    pts[i] += dtime;
    status = makeBezierSurfaceBody(context, stack, nCPu, nCPv, pts, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts[i] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, i, "Bezier Surface", 1e-6, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, i, "Bezier Surface", 1e-6, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setBezierSurfaceBody_dot(nCPu, nCPv, pts, pts_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Bezier Surface", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Offset Surface                                                           */
/*                                                                           */
/*****************************************************************************/

int
makeOffsetSurfaceBody( ego context,         /* (in)  EGADS context    */
                       objStack *stack,     /* (in)  EGADS obj stack  */
                       const double *xcent, /* (in)  Center           */
                       const double *xax,   /* (in)  x-axis           */
                       const double *yax,   /* (in)  y-axis           */
                       const double offset, /* (in)  offset           */
                       ego *ebody )         /* (out) Plane body       */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, oclass, mtype, *ivec=NULL;
  double data[9], dx[3], dy[3], dz[3], *rvec=NULL;
  ego    eplane, esurf, eedges[8], enodes[4], eloop, eface, eref;

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

  /* compute the normal of the plane */
  CROSS(dz, dx, dy);

  /* create the Offset surface */
  data[0] = offset;
  status = EG_makeGeometry(context, SURFACE, OFFSET, eplane, NULL,
                           data, &esurf);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, esurf);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create the Nodes for the Edges */
  data[0] = xcent[0] - dx[0] - dy[0] + dz[0]*offset;
  data[1] = xcent[1] - dx[1] - dy[1] + dz[1]*offset;
  data[2] = xcent[2] - dx[2] - dy[2] + dz[2]*offset;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dx[0] - dy[0] + dz[0]*offset;
  data[1] = xcent[1] + dx[1] - dy[1] + dz[1]*offset;
  data[2] = xcent[2] + dx[2] - dy[2] + dz[2]*offset;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] + dx[0] + dy[0] + dz[0]*offset;
  data[1] = xcent[1] + dx[1] + dy[1] + dz[1]*offset;
  data[2] = xcent[2] + dx[2] + dy[2] + dz[2]*offset;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = xcent[0] - dx[0] + dy[0] + dz[0]*offset;
  data[1] = xcent[1] - dx[1] + dy[1] + dz[1]*offset;
  data[2] = xcent[2] - dx[2] + dy[2] + dz[2]*offset;
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


  /* create P-curves */
  data[0] = -1.;    data[1] = -1.;    /* u == -1 UMIN       */
  data[2] =  0.;    data[3] =  1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = -1.;    data[1] = -1.;    /* v == -1 VMIN */
  data[2] =  1.;    data[3] =  0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 1.;    data[1] = -1.;     /* u ==  1 UMAX   */
  data[2] = 0.;    data[3] =  1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = -1.;    data[1] = 1.;     /* v ==  1 VMAX */
  data[2] =  1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and the face body */
  status = EG_makeTopology(context, esurf, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, esurf, FACE, SFORWARD,
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
setOffsetSurfaceBody_dot( const double *xcent,     /* (in)  Center          */
                          const double *xcent_dot, /* (in)  Center velocity */
                          const double *xax,       /* (in)  x-axis          */
                          const double *xax_dot,   /* (in)  x-axis velocity */
                          const double *yax,       /* (in)  y-axis          */
                          const double *yax_dot,   /* (in)  y-axis velocity */
                          const double offset,     /* (in)  offset          */
                          const double offset_dot, /* (in)  offset_dot      */
                          ego ebody )              /* (in/out) body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses, *ivec=NULL;
  double data[10], data_dot[10], dx[3], dx_dot[3], dy[3], dy_dot[3], dz[3], dz_dot[3], *rvec=NULL, *rvec_dot=NULL;
  ego    eplane, esurf, *efaces, *eloops, *eedges, enodes[4], *enode, eref;

  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Offset surface from the Face */
  status = EG_getTopology(efaces[0], &esurf, &oclass, &mtype,
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

  /* get the Plane from the offset surface */
  status = EG_getGeometry(esurf, &oclass, &mtype, &eplane, &ivec, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(ivec); ivec=NULL;
  EG_free(rvec); rvec=NULL;


  /* set velocity of the Offset surface */
  data[0] = offset;
  data_dot[0] = offset_dot;
  status = EG_setGeometry_dot(esurf, SURFACE, OFFSET, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

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

  /* compute the axis of the sphere between the poles */
  CROSS(dz, dx, dy);
  CROSS_DOT(dz_dot, dx, dx_dot, dy, dy_dot);

  /* create the Nodes for the Edges */
  data[0]     = xcent[0]     - dx[0]     - dy[0]     + dz[0]*offset;
  data[1]     = xcent[1]     - dx[1]     - dy[1]     + dz[1]*offset;
  data[2]     = xcent[2]     - dx[2]     - dy[2]     + dz[2]*offset;
  data_dot[0] = xcent_dot[0] - dx_dot[0] - dy_dot[0] + dz_dot[0]*offset + dz[0]*offset_dot;
  data_dot[1] = xcent_dot[1] - dx_dot[1] - dy_dot[1] + dz_dot[1]*offset + dz[1]*offset_dot;
  data_dot[2] = xcent_dot[2] - dx_dot[2] - dy_dot[2] + dz_dot[2]*offset + dz[2]*offset_dot;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     + dx[0]     - dy[0]     + dz[0]*offset;
  data[1]     = xcent[1]     + dx[1]     - dy[1]     + dz[1]*offset;
  data[2]     = xcent[2]     + dx[2]     - dy[2]     + dz[2]*offset;
  data_dot[0] = xcent_dot[0] + dx_dot[0] - dy_dot[0] + dz_dot[0]*offset + dz[0]*offset_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1] - dy_dot[1] + dz_dot[1]*offset + dz[1]*offset_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2] - dy_dot[2] + dz_dot[2]*offset + dz[2]*offset_dot;
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     + dx[0]     + dy[0]     + dz[0]*offset;
  data[1]     = xcent[1]     + dx[1]     + dy[1]     + dz[1]*offset;
  data[2]     = xcent[2]     + dx[2]     + dy[2]     + dz[2]*offset;
  data_dot[0] = xcent_dot[0] + dx_dot[0] + dy_dot[0] + dz_dot[0]*offset + dz[0]*offset_dot;
  data_dot[1] = xcent_dot[1] + dx_dot[1] + dy_dot[1] + dz_dot[1]*offset + dz[1]*offset_dot;
  data_dot[2] = xcent_dot[2] + dx_dot[2] + dy_dot[2] + dz_dot[2]*offset + dz[2]*offset_dot;
  status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0]     = xcent[0]     - dx[0]     + dy[0]     + dz[0]*offset;
  data[1]     = xcent[1]     - dx[1]     + dy[1]     + dz[1]*offset;
  data[2]     = xcent[2]     - dx[2]     + dy[2]     + dz[2]*offset;
  data_dot[0] = xcent_dot[0] - dx_dot[0] + dy_dot[0] + dz_dot[0]*offset + dz[0]*offset_dot;
  data_dot[1] = xcent_dot[1] - dx_dot[1] + dy_dot[1] + dz_dot[1]*offset + dz[1]*offset_dot;
  data_dot[2] = xcent_dot[2] - dx_dot[2] + dy_dot[2] + dz_dot[2]*offset + dz[2]*offset_dot;
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
  EG_free(ivec); ivec=NULL;
  EG_free(rvec); rvec=NULL;
  EG_free(rvec_dot); rvec_dot=NULL;

  return status;
}


int
pingOffsetSurface(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    iparam, np1, nt1, iedge, nedge, iface, nface;
  double x[10], x_dot[10], params[3], dtime = 1e-8;
  double *xcent, *xax, *yax, *xcent_dot, *xax_dot, *yax_dot;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  xcent = x;
  xax   = x+3;
  yax   = x+6;

  xcent_dot = x_dot;
  xax_dot   = x_dot+3;
  yax_dot   = x_dot+6;

  /* make the Plane Face body */
  xcent[0] = 0.00; xcent[1] = 0.00; xcent[2] = 0.00;
  xax[0]   = 1.10; xax[1]   = 0.10; xax[2]   = 0.05;
  yax[0]   = 0.05; yax[1]   = 1.20; yax[2]   = 0.10;
  //xax[0]   = 1.; xax[1]   = 0.; xax[2]   = 0.;
  //yax[0]   = 0.; yax[1]   = 1.; yax[2]   = 0.;
  x[9] = 1.1;
  status = makeOffsetSurfaceBody(context, stack, xcent, xax, yax, x[9], &ebody1);
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


  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Offset Surface Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" Offset Surface Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (iparam = 0; iparam < 10; iparam++) x_dot[iparam] = 0;

  for (iparam = 0; iparam < 10; iparam++) {

    /* set the velocity of the original body */
    x_dot[iparam] = 1.0;
    status = setOffsetSurfaceBody_dot(xcent, xcent_dot,
                                      xax, xax_dot,
                                      yax, yax_dot,
                                      x[9], x_dot[9], ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    x_dot[iparam] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed Circle for finite difference */
    x[iparam] += dtime;
    status = makeOffsetSurfaceBody(context, stack, xcent, xax, yax, x[9], &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    x[iparam] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, iparam, "Offset Surface", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, iparam, "Offset Surface", 1e-7, 1e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setOffsetSurfaceBody_dot(xcent, xcent_dot,
                                    xax, xax_dot,
                                    yax, yax_dot,
                                    x[9], x_dot[9], ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "Offset Surface", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  B-spline Surface                                                         */
/*                                                                           */
/*****************************************************************************/

int
makeBsplineSurfaceBody( ego context,        /* (in)  EGADS context         */
                        objStack *stack,    /* (in)  EGADS obj stack       */
                        const int nCPu,     /* (in)  number of points in u */
                        const int nCPv,     /* (in)  number of points in v */
                        const double *pts,  /* (in)  coordinates           */
                        ego *ebody )        /* (out) Bezier Surface body   */
{
  int    status = EGADS_SUCCESS;
  int    senses[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE};
  int    i, j, sizes[2], *header=NULL, lheader[4], oclass, mtype;
  double data[4], tdata[2], *lpts=NULL, *rvec=NULL;
  ego    esurf, ecurves[4], eedges[8], enodes[4], eloop, eface, nodes[2], eref;

  /* create the B-spline surface */
  sizes[0] = nCPu;
  sizes[1] = nCPv;
  status = EG_approximate(context, 0, DXYTOL, sizes, pts, &esurf);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, esurf);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &header, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create Nodes for the Edges */
  i = 0;
  j = 0;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = 0;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = nCPv-1;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = 0;
  j = nCPv-1;
  data[0] = pts[3*(i+j*nCPu)+0];
  data[1] = pts[3*(i+j*nCPu)+1];
  data[2] = pts[3*(i+j*nCPu)+2];
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, enodes[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the curves for the Edges */
  /* umin */
  status = EG_isoCurve(header, rvec, 0, -1, lheader, &lpts);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, lheader,
                           lpts, &ecurves[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  status = EG_stackPush(stack, ecurves[0]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* vmin */
  status = EG_isoCurve(header, rvec, -1, 0, lheader, &lpts);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, lheader,
                           lpts, &ecurves[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  status = EG_stackPush(stack, ecurves[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* umax */
  status = EG_isoCurve(header, rvec, header[2]-1, -1, lheader, &lpts);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, lheader,
                           lpts, &ecurves[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  status = EG_stackPush(stack, ecurves[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* vmax */
  status = EG_isoCurve(header, rvec, -1, header[5]-1, lheader, &lpts);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, lheader,
                           lpts, &ecurves[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  status = EG_stackPush(stack, ecurves[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Curves */
  tdata[0] = 0;
  tdata[1] = 1;

  nodes[0] = enodes[0];
  nodes[1] = enodes[3];
  status = EG_makeTopology(context, ecurves[0], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  nodes[0] = enodes[0];
  nodes[1] = enodes[1];
  status = EG_makeTopology(context, ecurves[1], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  nodes[0] = enodes[1];
  nodes[1] = enodes[2];
  status = EG_makeTopology(context, ecurves[2], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  nodes[0] = enodes[3];
  nodes[1] = enodes[2];
  status = EG_makeTopology(context, ecurves[3], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 0.;    /* v == 0 VMIN */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 1;     data[1] = 0.;    /* u == 1 UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = 1.;    /* v == 1 VMAX */
  data[2] = 1.;    data[3] = 0.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the loop and Face body */
  status = EG_makeTopology(context, esurf, LOOP, CLOSED,
                           NULL, 4, eedges, senses, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_stackPush(stack, eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, esurf, FACE, SFORWARD,
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
  EG_free(header); header = NULL;
  EG_free(lpts); lpts = NULL;
  EG_free(rvec); rvec = NULL;

  return status;
}


int
setBsplineSurfaceBody_dot( const int nCPu,        /* (in)  number of points in u   */
                          const int nCPv,        /* (in)  number of points in v   */
                          const double *pts,     /* (in)  coordinates             */
                          const double *pts_dot, /* (in)  velocity of coordinates */
                          ego ebody )            /* (in/out) body with velocities */
{
  int    status = EGADS_SUCCESS;
  int    i, j, nnode, nedge, nloop, nface, oclass, mtype, *senses, sizes[2], *header=NULL, lheader[4];
  double data[18], data_dot[18], tdata[2], tdata_dot[2];
  double *lpts=NULL, *lpts_dot=NULL, *rvec=NULL, *rvec_dot=NULL;
  ego    esurf, ecurves[4], enodes[4], *nodes, *eloops, *eedges, *efaces, eref;


  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &efaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop and Sphere from the Face */
  status = EG_getTopology(efaces[0], &esurf, &oclass, &mtype,
                          data, &nloop, &eloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &eedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Curves and Nodes from the Edges */
  status = EG_getTopology(eedges[0], &ecurves[0], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = nodes[0];
  enodes[3] = nodes[1];

  status = EG_getTopology(eedges[1], &ecurves[1], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[0] = nodes[0];
  enodes[1] = nodes[1];

  status = EG_getTopology(eedges[2], &ecurves[2], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[1] = nodes[0];
  enodes[2] = nodes[1];

  status = EG_getTopology(eedges[3], &ecurves[3], &oclass, &mtype,
                          data, &nnode, &nodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;
  enodes[3] = nodes[0];
  enodes[2] = nodes[1];

  /* set the surface sensitivity */
  sizes[0] = nCPu;
  sizes[1] = nCPv;
  status = EG_approximate_dot(esurf, 0, DXYTOL, sizes, pts, pts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &header, &rvec);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(rvec);

  status = EG_getGeometry_dot(esurf, &rvec, &rvec_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Node velocities */
  i = 0;
  j = 0;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = 0;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = nCPu-1;
  j = nCPv-1;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[2], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  i = 0;
  j = nCPv-1;
  data[0]     = pts[3*(i+j*nCPu)+0];
  data[1]     = pts[3*(i+j*nCPu)+1];
  data[2]     = pts[3*(i+j*nCPu)+2];
  data_dot[0] = pts_dot[3*(i+j*nCPu)+0];
  data_dot[1] = pts_dot[3*(i+j*nCPu)+1];
  data_dot[2] = pts_dot[3*(i+j*nCPu)+2];
  status = EG_setGeometry_dot(enodes[3], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* set the Curve velocities */
  /* umin */
  status = EG_isoCurve_dot(header, rvec, rvec_dot, 0, -1, lheader, &lpts, &lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(ecurves[0], CURVE, BSPLINE, lheader, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  EG_free(lpts_dot); lpts_dot = NULL;

  /* vmin */
  status = EG_isoCurve_dot(header, rvec, rvec_dot, -1, 0, lheader, &lpts, &lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(ecurves[1], CURVE, BSPLINE, lheader, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  EG_free(lpts_dot); lpts_dot = NULL;

  /* umax */
  status = EG_isoCurve_dot(header, rvec, rvec_dot, header[2]-1, -1, lheader, &lpts, &lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(ecurves[2], CURVE, BSPLINE, lheader, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  EG_free(lpts_dot); lpts_dot = NULL;

  /* vmax */
  status = EG_isoCurve_dot(header, rvec, rvec_dot, -1, header[5]-1, lheader, &lpts, &lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_setGeometry_dot(ecurves[3], CURVE, BSPLINE, lheader, lpts, lpts_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_free(lpts); lpts = NULL;
  EG_free(lpts_dot); lpts_dot = NULL;

  /* set the Edge t-range sensitivity */
  tdata[0]     = 0;
  tdata[1]     = 1;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;
  status = EG_setRange_dot(eedges[0], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[1], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[2], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_setRange_dot(eedges[3], EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_free(header);
  EG_free(lpts);
  EG_free(lpts_dot);
  EG_free(rvec);
  EG_free(rvec_dot);

  return status;
}


int
pingBsplineSurface(ego context, objStack *stack)
{
  int    status = EGADS_SUCCESS;
  int    i, np1, nt1, iedge, nedge, iface, nface;
  double params[3], dtime = 1e-8;
  const int    *pt1, *pi1, *ts1, *tc1;
  const double *t1, *x1, *uv1;
  ego    ebody1, ebody2, tess1, tess2;

  const int nCPu = 4;
  const int nCPv = 4;
  double pts[3*4*4] = {0.00, 0.00, 0.00,
                       1.00, 0.00, 0.10,
                       1.50, 1.00, 0.70,
                       0.25, 0.75, 0.60,

                       0.00, 0.00, 1.00,
                       1.00, 0.00, 1.10,
                       1.50, 1.00, 1.70,
                       0.25, 0.75, 1.60,

                       0.00, 0.00, 2.00,
                       1.00, 0.00, 2.10,
                       1.50, 1.00, 2.70,
                       0.25, 0.75, 2.60,

                       0.00, 0.00, 3.00,
                       1.00, 0.00, 3.10,
                       1.50, 1.00, 3.70,
                       0.25, 0.75, 3.60};
  double pts_dot[3*4*4];

  status = makeBsplineSurfaceBody(context, stack, nCPu, nCPv, pts, &ebody1);
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
  params[1] =  0.1;
  params[2] = 20.0;
  status = EG_makeTessBody(ebody1, params, &tess1);

  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess1, iedge+1, &np1, &x1, &t1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" B-spline Surface Edge %d np1 = %d\n", iedge+1, np1);
  }

  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess1, iface+1, &np1, &x1, &uv1, &pt1, &pi1,
                                            &nt1, &ts1, &tc1);
    if (status != EGADS_SUCCESS) goto cleanup;
    printf(" B-spline Surface Face %d np1 = %d\n", iface+1, np1);
  }

  /* zero out velocities */
  for (i = 0; i < 3*nCPu*nCPv; i++) pts_dot[i] = 0;

  for (i = 0; i < 3*nCPu*nCPv; i++) {

    /* set the velocity of the original body */
    pts_dot[i] = 1.0;
    status = setBsplineSurfaceBody_dot(nCPu, nCPv, pts, pts_dot, ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts_dot[i] = 0.0;

    status = EG_hasGeometry_dot(ebody1);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* make a perturbed body for finite difference */
    pts[i] += dtime;
    status = makeBsplineSurfaceBody(context, stack, nCPu, nCPv, pts, &ebody2);
    if (status != EGADS_SUCCESS) goto cleanup;
    pts[i] -= dtime;

    /* map the tessellation */
    status = EG_mapTessBody(tess1, ebody2, &tess2);
    if (status != EGADS_SUCCESS) goto cleanup;

    /* ping the bodies */
    status = pingBodies(tess1, tess2, dtime, i, "B-spline Surface", 1e-6, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = pingBodiesExtern(tess1, ebody2, dtime, i, "B-spline Surface", 1e-6, 5e-7, 1e-7);
    if (status != EGADS_SUCCESS) goto cleanup;

    EG_deleteObject(tess2);
  }

  /* zero out sensitivities */
  status = setBsplineSurfaceBody_dot(nCPu, nCPv, pts, pts_dot, ebody1);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check transformations */
  status = pingTransform(ebody1, params, "B-spline Surface", 5e-7, 5e-7, 5e-7);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(tess1);

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  return status;
}


/*****************************************************************************/
/*                                                                           */
/*  Check functionality of _dot functions                                    */
/*                                                                           */
/*****************************************************************************/

int
checkNode_dot(ego context)
{
  int    status = EGADS_SUCCESS;
  double data[18], data_dot[18];
  ego    enode;

  /* make a node */
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  status = EG_makeTopology(context, NULL, NODE, 0,
                           data, 0, NULL, NULL, &enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it does not have sensitivities */
  status = EG_hasGeometry_dot(enode);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* attempt to set sensitivity with incorrect data */
  printf("Check error handling...\n");
  data[0] = 2;
  data[1] = 0;
  data[2] = 1;
  data_dot[0] = 3;
  data_dot[1] = 4;
  data_dot[2] = 5;
  status = EG_setGeometry_dot(enode, NODE, 0, NULL, data, data_dot);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* add sensitivity */
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  status = EG_setGeometry_dot(enode, NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it has sensitivities */
  status = EG_hasGeometry_dot(enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* attempt to remove sensitives with incorrect oclass */
  printf("Check error handling...\n");
  status = EG_setGeometry_dot(enode, LOOP, 0, NULL, NULL, NULL);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* remove sensitives */
  status = EG_setGeometry_dot(enode, NODE, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that the node does not have sensitivities */
  status = EG_hasGeometry_dot(enode);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* remove sensitives when none exist */
  status = EG_setGeometry_dot(enode, 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  EG_deleteObject(enode); enode = NULL;

  /* create a node directly with sensitivities */
  data[0] = 0;
  data[1] = 1;
  data[2] = 2;
  data_dot[0] = 3;
  data_dot[1] = 4;
  data_dot[2] = 5;
  status = EG_makeTopology_dot(context, NULL, NODE, 0,
                               data, data_dot, 0, NULL, NULL, &enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it has sensitivities */
  status = EG_hasGeometry_dot(enode);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }
  EG_deleteObject(enode);

  return status;
}


int
checkCurve_dot(ego context)
{
  int    status = EGADS_SUCCESS;
  int    esens[1] = {SFORWARD}, *senses, oclass, mtype, nloop, nedge, nnode;
  double data[18], data_dot[18], tdata[2], tdata_dot[2];
  ego    eline=NULL, enodes[2]={NULL,NULL}, eedge=NULL, eloop=NULL, ebody=NULL;
  ego    bline, *bnodes, *bedges, *bloops, eref, eloop2=NULL, ebody2=NULL;

  /* create the Line (point and direction) */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  data[3] = 2;
  data[4] = 1;
  data[5] = 3;
  status = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it does not have sensitivities */
  status = EG_hasGeometry_dot(eline);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }


  /* attempt to set sensitivity with incorrect data */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  data[3] = 1;
  data[4] = 2;
  data[5] = 3;
  data_dot[0] = 3;
  data_dot[1] = 4;
  data_dot[2] = 5;
  data_dot[3] = 6;
  data_dot[4] = 7;
  data_dot[5] = 8;
  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* add sensitivity */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  data[3] = 2;
  data[4] = 1;
  data[5] = 3;
  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it has sensitivities */
  status = EG_hasGeometry_dot(eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* attempt to remove sensitives with incorrect oclass */
  printf("Check error handling...\n");
  status = EG_setGeometry_dot(eline, SURFACE, LINE, NULL, NULL, NULL);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* attempt to remove sensitives with incorrect mtype */
  printf("Check error handling...\n");
  status = EG_setGeometry_dot(eline, CURVE, CIRCLE, NULL, NULL, NULL);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* remove sensitives */
  status = EG_setGeometry_dot(eline, CURVE, LINE, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that the node does not have sensitivities */
  status = EG_hasGeometry_dot(eline);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* remove sensitives when none exist */
  status = EG_setGeometry_dot(eline, 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;


  EG_deleteObject(eline); eline = NULL;

  /* crate a new line with sensitivities */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  data[3] = 2;
  data[4] = 1;
  data[5] = 3;
  data_dot[0] = 3;
  data_dot[1] = 4;
  data_dot[2] = 5;
  data_dot[3] = 6;
  data_dot[4] = 7;
  data_dot[5] = 8;
  status = EG_makeGeometry_dot(context, CURVE, LINE, NULL, NULL,
                               data, data_dot, &eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it has sensitivities */
  status = EG_hasGeometry_dot(eline);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Nodes for the Edge with sensitvities */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  status = EG_makeTopology_dot(context, NULL, NODE, 0,
                               data, data_dot, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = data[3];
  data[1] = data[4];
  data[2] = data[5];
  status = EG_makeTopology_dot(context, NULL, NODE, 0,
                               data, data_dot, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /******************************/
  /* Edge                       */
  /******************************/

  /* make the Edge on the Line */
  tdata[0] = 0;
  tdata[1] = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);

  tdata_dot[0] = 0;
  tdata_dot[1] = (data[3]*data_dot[3] + data[4]*data_dot[4] + data[5]*data_dot[5])/tdata[1];

  status = EG_makeTopology_dot(context, eline, EDGE, TWONODE,
                               tdata, tdata_dot, 2, enodes, NULL, &eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot checks the Edge t-range, curve and nodes for an edge */
  status = EG_hasGeometry_dot(eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* remove the t-range sensitivity */
  status = EG_setRange_dot(eedge, EDGE, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot checks the Edge t-range, curve and nodes for an edge */
  status = EG_hasGeometry_dot(eedge);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  status = EG_setRange_dot(eedge, EDGE, tdata, tdata_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_hasGeometry_dot(eedge);
  if (status != EGADS_SUCCESS) goto cleanup;

  /******************************/
  /* Loop                       */
  /******************************/

  status = EG_makeTopology_dot(context, NULL, LOOP, OPEN,
                               NULL, NULL, 1, &eedge, esens, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the loop has sensitivities */
  status = EG_hasGeometry_dot(eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = EG_makeTopology(context, NULL, LOOP, OPEN,
                           NULL, 1, &eedge, esens, &eloop2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the loop has sensitivities */
  status = EG_hasGeometry_dot(eloop2);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* remove sensitives from one node */
  status = EG_setGeometry_dot(enodes[0], 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot should now be NOTFOUND */
  status = EG_hasGeometry_dot(eloop);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* hasGeometry_dot should now be NOTFOUND */
  status = EG_hasGeometry_dot(eloop2);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }
  EG_deleteObject(eloop2); eloop2=NULL;

  /* add the sensitivity back */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  status = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /***********************************/
  /* Body                            */
  /*                                 */
  /* Bodies create new ego instances */
  /***********************************/

  /* make the WIREBODY */
  status = EG_makeTopology_dot(context, NULL, BODY, WIREBODY,
                               NULL, NULL, 1, &eloop, NULL, &ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the body has sensitivities */
  status = EG_hasGeometry_dot(ebody);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* make the WIREBODY */
  status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                           NULL, 1, &eloop, NULL, &ebody2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the body has sensitivities */
  status = EG_hasGeometry_dot(ebody2);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_deleteObject(ebody2); ebody2=NULL;


  /* get the Loop from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nloop, &bloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(bloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &bedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes and the Line from the Edge */
  status = EG_getTopology(bedges[0], &bline, &oclass, &mtype,
                          data, &nnode, &bnodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* remove sensitives from one node */
  status = EG_setGeometry_dot(bnodes[0], 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot should now be NOTFOUND */
  status = EG_hasGeometry_dot(ebody);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* add the sensitivity back */
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  status = EG_setGeometry_dot(bnodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the body has sensitivities again */
  status = EG_hasGeometry_dot(ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* remove all sensitivities from all geometry */
  status = EG_setGeometry_dot(ebody, BODY, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* redundant check */
  status = EG_setGeometry_dot(ebody, 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot should now be NOTFOUND for all geometry */
  status = EG_hasGeometry_dot(bline);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  status = EG_hasGeometry_dot(bnodes[0]);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  status = EG_hasGeometry_dot(bnodes[1]);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  EG_deleteObject(ebody);
  EG_deleteObject(eloop);
  EG_deleteObject(eedge);
  EG_deleteObject(enodes[0]);
  EG_deleteObject(enodes[1]);
  EG_deleteObject(eline);

  return status;
}


int
checkSurface_dot(ego context)
{
  int    status = EGADS_SUCCESS;
  int    i, esens[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE}, *senses;
  int    oclass, mtype, nface, nloop, nedge, nnode;
  double data[18], data_dot[18], tdata[2], tdata_dot[2];
  ego    esphere=NULL, ecircle=NULL, enodes[2]={NULL,NULL};
  ego    eedges[8]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL}, eloop=NULL, eface=NULL, ebody=NULL;
  ego    bcircle, *bnodes, *bedges, *bloops, *bfaces, eref, eface2=NULL;

  /* create the Sphere */
  data[0] = 0;   /* center */
  data[1] = 0;
  data[2] = 0;
  data[3] = 1;   /* x-axis */
  data[4] = 0;
  data[5] = 0;
  data[6] = 0;   /* y-axis */
  data[7] = 1;
  data[8] = 0;
  data[9] = 2;   /* radius */
  status = EG_makeGeometry(context, SURFACE, SPHERICAL, NULL, NULL,
                           data, &esphere);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it does not have sensitivities */
  status = EG_hasGeometry_dot(esphere);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }


  /* attempt to set sensitivity with incorrect data */
  printf("Check error handling...\n");
  data[0] = 0;   /* center */
  data[1] = 0;
  data[2] = 0;
  data[3] = 1;   /* x-axis */
  data[4] = 0;
  data[5] = 0;
  data[6] = 0;   /* y-axis */
  data[7] = 1;
  data[8] = 0;
  data[9] = 1.4;   /* radius */
  data_dot[0] = 3;
  data_dot[1] = 4;
  data_dot[2] = 5;
  data_dot[3] = 6;
  data_dot[4] = 7;
  data_dot[5] = 8;
  data_dot[6] = 9;
  data_dot[7] = 10;
  data_dot[8] = 11;
  data_dot[9] = 12;
  status = EG_setGeometry_dot(esphere, SURFACE, SPHERICAL, NULL, data, data_dot);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* add sensitivity */
  data[0] = 0;   /* center */
  data[1] = 0;
  data[2] = 0;
  data[3] = 1;   /* x-axis */
  data[4] = 0;
  data[5] = 0;
  data[6] = 0;   /* y-axis */
  data[7] = 1;
  data[8] = 0;
  data[9] = 2;   /* radius */
  status = EG_setGeometry_dot(esphere, SURFACE, SPHERICAL, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that it has sensitivities */
  status = EG_hasGeometry_dot(esphere);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* attempt to remove sensitives with incorrect oclass */
  printf("Check error handling...\n");
  status = EG_setGeometry_dot(esphere, CURVE, SPHERICAL, NULL, NULL, NULL);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* attempt to remove sensitives with incorrect mtype */
  printf("Check error handling...\n");
  status = EG_setGeometry_dot(esphere, SURFACE, CONICAL, NULL, NULL, NULL);
  if (status != EGADS_GEOMERR) {
    status = EGADS_GEOMERR;
    goto cleanup;
  }

  /* remove sensitives */
  status = EG_setGeometry_dot(esphere, SURFACE, SPHERICAL, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check that the node does not have sensitivities */
  status = EG_hasGeometry_dot(esphere);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* remove sensitives when none exist */
  status = EG_setGeometry_dot(esphere, 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* add sensitivity to the sphere again */
  status = EG_setGeometry_dot(esphere, SURFACE, SPHERICAL, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* create circle for the Edge */
  data[0] = 0;   /* center */
  data[1] = 0;
  data[2] = 0;
  data[3] = 1;   /* x-axis */
  data[4] = 0;
  data[5] = 0;
  data[6] = 0;   /* y-axis */
  data[7] = 0;
  data[8] = 1;
  data[9] = 2;   /* radius */
  status = EG_makeGeometry_dot(context, CURVE, CIRCLE, NULL, NULL,
                               data, data_dot, &ecircle);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create Nodes with sensitivities for the Edge */
  data[0] =  0;
  data[1] =  0;
  data[2] = -2;
  status = EG_makeTopology_dot(context, NULL, NODE, 0,
                               data, data_dot, 0, NULL, NULL, &enodes[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0;
  data[1] = 0;
  data[2] = 2;
  status = EG_makeTopology_dot(context, NULL, NODE, 0,
                               data, data_dot, 0, NULL, NULL, &enodes[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Edge on the Circle */
  tdata[0] = -PI/2.;
  tdata[1] =  PI/2.;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_makeTopology_dot(context, ecircle, EDGE, TWONODE,
                               tdata, tdata_dot, 2, enodes, NULL, &eedges[0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Degenerate Edges on the Nodes */
  tdata[0] = 0;
  tdata[1] = TWOPI;
  tdata_dot[0] = 0;
  tdata_dot[1] = 0;

  status = EG_makeTopology_dot(context, NULL, EDGE, DEGENERATE,
                               tdata, tdata_dot, 1, &enodes[0], NULL, &eedges[1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  eedges[2] = eedges[0]; /* repeat the circle edge */

  status = EG_makeTopology_dot(context, NULL, EDGE, DEGENERATE,
                               tdata, tdata_dot, 1, &enodes[1], NULL, &eedges[3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create P-curves */
  data[0] = 0.;    data[1] = 0.;    /* u == 0 UMIN       */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+0]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] = -PI/2; /* v == -PI/2 VMIN */
  data[2] = 1.;    data[3] =  0.  ;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+1]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = TWOPI; data[1] = 0.;    /* u == TWOPI UMAX   */
  data[2] = 0.;    data[3] = 1.;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+2]);
  if (status != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.;    data[1] =  PI/2.; /* v == PI/2 VMAX */
  data[2] = 1.;    data[3] =  0.   ;
  status = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &eedges[4+3]);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make the Loop */
  status = EG_makeTopology_dot(context, esphere, LOOP, CLOSED,
                               NULL, NULL, 4, eedges, esens, &eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the loop has sensitivities */
  status = EG_hasGeometry_dot(eloop);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make a Face */
  status = EG_makeTopology_dot(context, esphere, FACE, SFORWARD,
                               NULL, NULL, 1, &eloop, NULL, &eface);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the Face has sensitivities */
  status = EG_hasGeometry_dot(eface);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* make a Face2 */
  status = EG_makeTopology(context, esphere, FACE, SFORWARD,
                           NULL, 1, &eloop, NULL, &eface2);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the Face has sensitivities */
  status = EG_hasGeometry_dot(eface2);
  if (status != EGADS_SUCCESS) goto cleanup;
  EG_deleteObject(eface2); eface2 = NULL;


  /* make the Body */
  status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                           NULL, 1, &eface, NULL, &ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the Body has sensitivities */
  status = EG_hasGeometry_dot(ebody);
  if (status != EGADS_SUCCESS) goto cleanup;


  /* get the Face from the Body */
  status = EG_getTopology(ebody, &eref, &oclass, &mtype,
                          data, &nface, &bfaces, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Loop from the Face */
  status = EG_getTopology(bfaces[0], &eref, &oclass, &mtype,
                          data, &nloop, &bloops, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Edge from the Loop */
  status = EG_getTopology(bloops[0], &eref, &oclass, &mtype,
                          data, &nedge, &bedges, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* get the Nodes from the Edge */
  status = EG_getTopology(bedges[0], &bcircle, &oclass, &mtype,
                          data, &nnode, &bnodes, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* remove sensitives from one node */
  status = EG_setGeometry_dot(bnodes[0], 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot should now be NOTFOUND */
  status = EG_hasGeometry_dot(ebody);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  /* add the sensitivity back */
  data[0] =  0;
  data[1] =  0;
  data[2] = -2;
  status = EG_setGeometry_dot(bnodes[0], NODE, 0, NULL, data, data_dot);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* check all geometry in the body has sensitivities again */
  status = EG_hasGeometry_dot(ebody);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* remove all sensitivities from all geometry */
  status = EG_setGeometry_dot(ebody, BODY, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* redundant check */
  status = EG_setGeometry_dot(ebody, 0, 0, NULL, NULL, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* hasGeometry_dot should now be NOTFOUND for all geometry */
  status = EG_hasGeometry_dot(bcircle);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  status = EG_hasGeometry_dot(bnodes[0]);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }

  status = EG_hasGeometry_dot(bnodes[1]);
  if (status != EGADS_NOTFOUND) {
    status = EGADS_NODATA;
    goto cleanup;
  }


  status = EGADS_SUCCESS;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf(" Failure %d in %s\n", status, __func__);
  }

  EG_deleteObject(ebody);
  EG_deleteObject(eface);
  EG_deleteObject(eloop);
  for (i = 1; i < 8; i++) { /* skip the first repeated Edge */
    EG_deleteObject(eedges[i]);
  }
  EG_deleteObject(enodes[0]);
  EG_deleteObject(enodes[1]);
  EG_deleteObject(ecircle);
  EG_deleteObject(esphere);

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

  status = pingLine(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingCircle(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingEllipse(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingParabola(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingHyperbola(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* The trimmed curve disappears with OCC */
  //status = pingTrimmedCurve(context);
  //if (status != EGADS_SUCCESS) goto cleanup;

  status = pingOffsetCurve(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingBezierCurve(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingBsplineCurve(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingPlane(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingSpherical(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingConical(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingCylindrical(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingToroidal(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingRevolution(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingExtrusion(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingBezierSurface(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* The trimmed surface disappears with OCC */
//  status = pingTrimmedSurface(context);
//  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingOffsetSurface(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = pingBsplineSurface(context, &stack);
  if (status != EGADS_SUCCESS) goto cleanup;

  /******************************/
  /* Check _dot functions       */
  /******************************/

  status = checkNode_dot(context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = checkCurve_dot(context);
  if (status != EGADS_SUCCESS) goto cleanup;

  status = checkSurface_dot(context);
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
