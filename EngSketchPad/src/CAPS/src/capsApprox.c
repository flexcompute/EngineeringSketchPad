/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             CAPS Spline Approximate functions
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>          /* Needed in some systems for DBL_MAX definition */
#include <math.h>

#include "prm.h"
#include "egads.h"
#include "egadsTris.h"
#include "capsTypes.h"


#define TOLCOPO         1.e-8   /* Tolerance for coincident points in
                                   normalized coordinates */
#define NOTFILLED      -1

extern void EG_makeConnect(int k1, int k2, int *tri, int *kedge, int *ntable,
                           connect *etable, int face);


static int caps_spline(int natural, int nu, double *r, double *aux, double *t)
{
    int    i;
    double rv, bet, *diag, *u, *z;

   /* Finds the tangents at the POINTs defining a ferguson spline */

   /* natural == 0  for fixed (prescribed) slopes at ends
              == 1  for zero second derivatives at ends
              == 2  for periodic ends (r[0] == r[nu-1]) & t must be 4*nu in len
   */

    if (nu <= 1) {
        printf("caps_spline: cannot interpolate a spline with %d points\n", nu);
        return -1;
    } else if (nu == 2) {      /* Override natural condition */
        t[0] = r[1] - r[0]; 
        t[1] = t[0];
    } else if (nu == 3) {
        if (natural == 0) {
            t[1] = 0.25 * (3.0 * (r[2] - r[0]) - t[0] - t[2]); 
        } else {
            t[0] = -1.25 * r[0] + 1.5 * r[1] - 0.25 * r[2]; 
            t[1] = -0.50 * r[0]              + 0.50 * r[2]; 
            t[2] =  0.25 * r[0] - 1.5 * r[1] + 1.25 * r[2]; 
        }

    } else if (natural == 0) {

        rv   = 3.0 * (r[2] - r[0]) - t[0];
        bet  = 4.0;
        t[1] = rv / bet;

        for (i = 2; i < nu-2; i++) {
            aux[i] = 1.0 / bet; 
            rv     = 3.0 * (r[i+1] - r[i-1]);
            bet    = 4.0 - aux[i];
            t[i]   = (rv - t[i-1]) / bet;
        }

        aux[nu-2] = 1.0 / bet;
        rv        = 3.0 * (r[nu-1] - r[nu-3]) - t[nu-1];
        bet       = 4.0 - aux[nu-2];
        t[nu-2]   = (rv - t[nu-3]) / bet;

        for (i = nu-3; i > 0; i--) {
            t[i] -= aux[i+1] * t[i+1];
        }

    } else if (natural == 1) {

        rv   = 3.0 * (r[2] - r[0]) - 1.5 * (r[1] - r[0]);
        bet  = 3.5;
        t[1] = rv / bet;

        for (i = 2; i < nu-2; i++) {
            aux[i] = 1.0 / bet; 
            rv     = 3.0 * (r[i+1] - r[i-1]);
            bet    = 4.0 - aux[i];
            t[i]   = (rv - t[i-1]) / bet;
        }

        aux[nu-2] = 1.0 / bet;
        rv        = 3.0 * (r[nu-1] - r[nu-3]) - 1.5 * (r[nu-1] - r[nu-2]);
        bet       = 3.5 - aux[nu-2];
        t[nu-2]   = (rv - t[nu-3]) / bet;

        for (i = nu-3; i > 0; i--) {
            t[i] -= aux[i+1] * t[i+1];
        }

        t[0   ] = 1.5 * (r[1   ] - r[0   ]) - 0.5 * t[1   ];
        t[nu-1] = 1.5 * (r[nu-1] - r[nu-2]) - 0.5 * t[nu-2];

    } else {        /* natural == 2 */

        diag = aux  + nu;
        u    = diag + nu;
        z    = u    + nu;

        for (i = 0; i < nu-1; i++) diag[i] = 4.0;
        diag[0   ] += 4.0;
        diag[nu-2] += .25;

        bet  = diag[0];
        t[0] = 3.0 * (r[1] - r[nu-2]) / bet;

        for (i = 1; i < nu-1; i++) {
            aux[i] = 1.0 / bet;
            bet    = diag[i] - aux[i];
            t[i]   = (3.0 * (r[i+1] - r[i-1]) - t[i-1]) / bet;
        }

        for (i = nu-3; i >= 0; i--) t[i] -= aux[i+1] * t[i+1];

        for (i = 1; i < nu-1; i++) u[i] = 0;
        u[0   ] = -4.0;
        u[nu-2] =  1.0;
        
        bet  = diag[0];
        z[0] = u[0] / bet;

        for (i = 1; i < nu-1; i++) {
            aux[i] = 1.0 / bet;
            bet    = diag[i] - aux[i];
            z[i]   = (u[i] - z[i-1]) / bet;
        }

        for (i = nu-3; i >= 0; i--) z[i] -= aux[i+1] * z[i+1];

        bet = (t[0] - t[nu-2] / 4.0) / (1.0 + z[0] - z[nu-2] / 4.0);

        for (i = 0; i < nu-1; i++) t[i] -= bet * z[i];

        t[nu-1] = t[0];      /* enforce periodicity */
    }

    return 0;
} 


static int caps_fillCoeff1D(int nrank, int ntx, double *fit, double *coeff, 
                            double *r)
{
  /* NOTE: r must be at least 6*nt in length */
  int    i, j, nt, per = 1;
  double *t, *aux;

  nt  = ntx;
  if (nt < 0) {
    nt  = -nt;
    per = 2;
  }
  t   = r + nt;
  aux = t + nt;

  for (i = 0; i < nrank; i++) {
    for (j = 0; j < nt; j++) r[j] = fit[nrank*j+i];
    if (caps_spline(per, nt, r, aux, t)) return 1;
    for (j = 0; j < nt; j++) {
      coeff[2*nrank*j      +i] = r[j];
      coeff[2*nrank*j+nrank+i] = t[j];
    }
  }

  return 0;
}


static void caps_eval1D(int nrank, int nt, double *coeff, double t, double *sv,
                        /*@null@*/ double *dt1, /*@null@*/ double *dt2)
{
  int    i, l0, l1;
  double a, b, c;

  l0    = (int) t;
  if (l0 <     0) l0 = 0;
  if (l0 >= nt-1) l0 = nt-2;
  t     = t - (double) l0;
  l1    = l0 + 1;

  for (i = 0; i < nrank; i++) {
    c = coeff[2*nrank*l1+i] - coeff[2*nrank*l0+i];
    a =  3.0*c - 2.0*coeff[2*nrank*l0+nrank+i] - coeff[2*nrank*l1+nrank+i];
    b = -2.0*c +     coeff[2*nrank*l0+nrank+i] + coeff[2*nrank*l1+nrank+i];
    sv[i] =          coeff[2*nrank*l0+i      ] +
                  t*(coeff[2*nrank*l0+nrank+i] + t*(a + t*b));
    if (dt1 != NULL) dt1[i] = coeff[2*nrank*l0+nrank+i] + t*(2.0*a + 3.0*t*b);
    if (dt2 != NULL) dt2[i] = 2.0*a + 6.0*t*b;
  }
}


static void caps_invEval1D(int nrank, int ntx, double *coeff, double *sv, 
                           double *t)
{
  int    i, j, k, nt, l0, l1, per, newton, again = 0, count = 0;
  double dis, dis0, disn, tmin, tmax, a, b, c, d, t1, t2, r;
  double cs, kcs, jcs, step, tf, tm, a00, b0;

  per   = 0;
  nt    = ntx;
  if (ntx < 0) {
    nt  = -ntx;
    per = 1;
  }
  tmin  = cs = 0.0;
  tmax  = nt - 1;
  dis0  = DBL_MAX;

  for (i = 0; i < nt; i++) {
    dis = 0.0;
    for (j = 0; j < nrank; j++)
      dis += (coeff[2*nrank*i+j]-sv[j])*(coeff[2*nrank*i+j]-sv[j]);
    if (dis < dis0) {
       dis0 = dis;
       cs   = (double) i;
    }
  }

  for (k = 0, step = 0.5; k < 20; k++, step *= 0.5) {
    kcs = cs;
    for (i = 0, jcs = kcs-3.0*step; i < 4; i++, jcs += 2.0*step) {
      if ((jcs >= tmin) && (jcs <= tmax)) {
        l0  = (int) jcs;
        if (l0 <     0) l0 = 0;
        if (l0 >= nt-1) l0 = nt-2;
        tf  = jcs - (double) l0;
        l1  = l0 + 1;
        dis = 0.0;
        for (j = 0; j < nrank; j++) {
          c = coeff[2*nrank*l1+j] - coeff[2*nrank*l0+j];
          a =  3.0*c-2.0*coeff[2*nrank*l0+nrank+j] - coeff[2*nrank*l1+nrank+j];
          b = -2.0*c+    coeff[2*nrank*l0+nrank+j] + coeff[2*nrank*l1+nrank+j];
          d =            coeff[2*nrank*l0+j      ] +
                     tf*(coeff[2*nrank*l0+nrank+j] + tf*(a + tf*b));
          dis += (d-sv[j])*(d-sv[j]);
        }
        if (dis < dis0) {
          dis0 = dis;
          cs   = jcs;
        }
      }
    }
    if (cs < tmin || cs > tmax) {
      printf(" caps_invEval1D Info: cs = %10.5e %10.5e %10.5e\n",
             cs, tmin, tmax);
      *t = cs;
      return;
    }

    /* perform newton-raphson on close location */
    tm     = cs;
    disn   = dis0;
    newton = 1;
    do {
      if ((tm >= tmin) && (tm <= tmax)) {
        l0  = (int) tm;
        if (l0 <     0) l0 = 0;
        if (l0 >= nt-1) l0 = nt-2;
        tf  = tm - (double) l0;
        l1  = l0 + 1;
        dis = a00 = b0 = r = 0.0;
        for (j = 0; j < nrank; j++) {
          c  = coeff[2*nrank*l1+j] - coeff[2*nrank*l0+j];
          a  =  3.0*c-2.0*coeff[2*nrank*l0+nrank+j] - coeff[2*nrank*l1+nrank+j];
          b  = -2.0*c+    coeff[2*nrank*l0+nrank+j] + coeff[2*nrank*l1+nrank+j];
          d  =            coeff[2*nrank*l0+j      ] +
                      tf*(coeff[2*nrank*l0+nrank+j] + tf*(a + tf*b));
          t1 =           (coeff[2*nrank*l0+nrank+j] + tf*(2.0*a + 3.0*tf*b));
          t2 = 2.0*a + 6.0*tf*b;
          r   += t1*t1;
          dis += (d-sv[j])*(d-sv[j]);
          b0  -= (d-sv[j])*t1;
          a00 += (d-sv[j])*t2 + t1*t1;
        }
        if (dis > disn) {
          newton = 0;
          break;
        }
        disn = dis;
        c    = b0/a00;
        a00  = r*c*c;
        tm  += c;
        if (per == 1) {
          if (tm < tmin) {
            if (again == 1) {
              newton = 0;
              break;
            } else {
              again = 1;
              tm += (tmax - tmin);
            }
          }
          if (tm > tmax) {
            if (again == -1) {
              newton = 0;
              break;
            } else {
              again = -1;
              tm -= (tmax - tmin);
            }
          }
        }
        if (a00 < TOLCOPO*TOLCOPO) break;
      } else {
        newton = 0;
      }

      count++;
      if (count > 100) newton = 0;
    } while (newton);

    if (newton == 1) {
      *t = tm;
      return;
    }
  }

  *t = cs;
}


int caps_fillCoeff2D(int nrank, int nux, int nvx, double *fit, double *coeff,
                     double *r)
{
  /* NOTE: r must be at least 3*max(nux,nvx) in length -- 6 if any periodics */
  int    i, j, k, nu, nv, maxsize, peru = 1, perv = 1;
  double *t, *aux;

  nu = nux;
  nv = nvx;
  if (nu < 0) {
    nu   = -nu;
    peru = 2;
  }
  if (nv < 0) {
    nv   = -nv;
    perv = 2;
  }
  if (nv > nu) {
    maxsize = nv;
  } else {
    maxsize = nu;
  }
  t   = r + maxsize;
  aux = t + maxsize;

  for (i = 0; i < nrank; i++) {

    for (j = 0; j < nv; j++)            /* sv */
      for (k = 0; k < nu; k++)
        coeff[4*nrank*(j*nu+k)+i] = fit[nrank*(j*nu+k)+i];

    for (j = 0; j < nv; j++) {          /* du */
      for (k = 0; k < nu; k++) r[k] = fit[nrank*(j*nu+k)+i];
      if (caps_spline(peru, nu, r, aux, t)) return 1;
      for (k = 0; k < nu; k++) coeff[4*nrank*(j*nu+k)+  nrank+i] = t[k];
    }

    for (k = 0; k < nu; k++) {          /* dv */
      for (j = 0; j < nv; j++) r[j] = fit[nrank*(j*nu+k)+i];
      if (caps_spline(perv, nv, r, aux, t)) return 1;
      for (j = 0; j < nv; j++) coeff[4*nrank*(j*nu+k)+2*nrank+i] = t[j];
    }

    for (j = 0; j < nv; j++) r[j] = coeff[4*nrank*j*nu+nrank+i];
    if (caps_spline(perv, nv, r, aux, t)) return 1;
    for (j = 0; j < nv; j++) coeff[4*nrank*j*nu+3*nrank+i] = t[j];
    for (j = 0; j < nv; j++) r[j] = coeff[4*nrank*((j+1)*nu-1)+nrank+i];
    if (caps_spline(perv, nv, r, aux, t)) return 1;
    for (j = 0; j < nv; j++) coeff[4*nrank*((j+1)*nu-1)+3*nrank+i] = t[j];

    for (j = 0; j < nv; j++) {          /* duv */
      for (k = 0; k < nu; k++) r[k] = coeff[4*nrank*(j*nu+k)+2*nrank+i];
      if (caps_spline(peru, nu, r, aux, t)) return 1;
      for (k = 0; k < nu; k++) coeff[4*nrank*(j*nu+k)+3*nrank+i] = t[k];
    }
  }

  return 0;
}


static void caps_eval2D(int nrank, int nu, int nv, double *coeff, 
                        double *uv, double *sv, /*@null@*/ double *du,  
                        /*@null@*/ double *dv,  /*@null@*/ double *duu, 
                        /*@null@*/ double *duv, /*@null@*/ double *dvv)
{
  int    i, l0, l1, l2, l3, nNULL = 0;
  double u, v, s0, s1, s2, s3, t0, t1, t2, t3;
  double a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34;
  double a41, a42, a43, a44, s10, s20, s30, s40;

  if (du  == NULL) nNULL++;
  if (dv  == NULL) nNULL++;
  if (duu == NULL) nNULL++;
  if (duv == NULL) nNULL++;
  if (dvv == NULL) nNULL++;

  l0 = (int) uv[0];
  if (l0 <     0) l0 = 0;
  if (l0 >= nu-1) l0 = nu-2;
  u  = uv[0] - (double) l0;
  l1 = (int) uv[1];
  if (l1 <     0) l1 = 0;
  if (l1 >= nv-1) l1 = nv-2;
  v  = uv[1] - (double) l1;
  l0 = l0 + nu*l1;
  l1 = l0 + 1;
  l2 = l0 + nu;
  l3 = l2 + 1;

  for (i = 0; i < nrank; i++) {
    s0 = -3.0*coeff[4*nrank*l0        +i]  + 3.0*coeff[4*nrank*l2        +i] -
          2.0*coeff[4*nrank*l0+2*nrank+i]  -     coeff[4*nrank*l2+2*nrank+i];
    s1 = -3.0*coeff[4*nrank*l1        +i]  + 3.0*coeff[4*nrank*l3        +i] -
          2.0*coeff[4*nrank*l1+2*nrank+i]  -     coeff[4*nrank*l3+2*nrank+i];
    s2 = -3.0*coeff[4*nrank*l0+  nrank+i]  + 3.0*coeff[4*nrank*l2+  nrank+i] -
          2.0*coeff[4*nrank*l0+3*nrank+i]  -     coeff[4*nrank*l2+3*nrank+i];
    s3 = -3.0*coeff[4*nrank*l1+  nrank+i]  + 3.0*coeff[4*nrank*l3+  nrank+i] -
          2.0*coeff[4*nrank*l1+3*nrank+i]  -     coeff[4*nrank*l3+3*nrank+i];
    t0 =  2.0*coeff[4*nrank*l0        +i]  - 2.0*coeff[4*nrank*l2        +i] +
              coeff[4*nrank*l0+2*nrank+i]  +     coeff[4*nrank*l2+2*nrank+i];
    t1 =  2.0*coeff[4*nrank*l1        +i]  - 2.0*coeff[4*nrank*l3        +i] +
              coeff[4*nrank*l1+2*nrank+i]  +     coeff[4*nrank*l3+2*nrank+i];
    t2 =  2.0*coeff[4*nrank*l0+  nrank+i]  - 2.0*coeff[4*nrank*l2+  nrank+i] +
              coeff[4*nrank*l0+3*nrank+i]  +     coeff[4*nrank*l2+3*nrank+i];
    t3 =  2.0*coeff[4*nrank*l1+  nrank+i]  - 2.0*coeff[4*nrank*l3+  nrank+i] +
              coeff[4*nrank*l1+3*nrank+i]  +     coeff[4*nrank*l3+3*nrank+i];

    a11 = coeff[4*nrank*l0        +i];
    a12 = coeff[4*nrank*l0+2*nrank+i];
    a13 = s0;
    a14 = t0;
    a21 = coeff[4*nrank*l0+  nrank+i];
    a22 = coeff[4*nrank*l0+3*nrank+i];
    a23 = s2;
    a24 = t2;
    a31 = -3.0*coeff[4*nrank*l0        +i]  + 3.0*coeff[4*nrank*l1        +i] -
           2.0*coeff[4*nrank*l0+  nrank+i]  -     coeff[4*nrank*l1+  nrank+i];
    a32 = -3.0*coeff[4*nrank*l0+2*nrank+i]  + 3.0*coeff[4*nrank*l1+2*nrank+i] -
           2.0*coeff[4*nrank*l0+3*nrank+i]  -     coeff[4*nrank*l1+3*nrank+i];
    a33 = -3.0*s0 + 3.0*s1 - 2.0*s2 - s3;
    a34 = -3.0*t0 + 3.0*t1 - 2.0*t2 - t3;
    a41 =  2.0*coeff[4*nrank*l0        +i]  - 2.0*coeff[4*nrank*l1        +i] +
               coeff[4*nrank*l0+  nrank+i]  +     coeff[4*nrank*l1+  nrank+i];
    a42 =  2.0*coeff[4*nrank*l0+2*nrank+i]  - 2.0*coeff[4*nrank*l1+2*nrank+i] +
               coeff[4*nrank*l0+3*nrank+i]  +     coeff[4*nrank*l1+3*nrank+i];
    a43 =  2.0*s0 - 2.0*s1 + s2 + s3;
    a44 =  2.0*t0 - 2.0*t1 + t2 + t3;

    s10 = a11 + v*(a12 + v*(a13 + v*a14));
    s20 = a21 + v*(a22 + v*(a23 + v*a24));
    s30 = a31 + v*(a32 + v*(a33 + v*a34));
    s40 = a41 + v*(a42 + v*(a43 + v*a44));

    sv[i] = s10 + u*(s20 + u*(s30 + u*s40));
    if (nNULL == 5) continue;

    if (du  != NULL) du[i]  = s20 + u*(2.0*s30 + 3.0*u*s40);
    if (duu != NULL) duu[i] =          2.0*s30 + 6.0*u*s40;

    s10 = a12 + v*(2.0*a13 + 3.0*v*a14);
    s20 = a22 + v*(2.0*a23 + 3.0*v*a24);
    s30 = a32 + v*(2.0*a33 + 3.0*v*a34);
    s40 = a42 + v*(2.0*a43 + 3.0*v*a44);
    if (dv  != NULL) dv[i]  = s10 + u*(s20 + u*(s30 +u*s40));
    if (duv != NULL) duv[i] = s20 + u*(2.0*s30 + 3.0*u*s40);

    s10 = 2.0*a13 + 6.0*v*a14;
    s20 = 2.0*a23 + 6.0*v*a24;
    s30 = 2.0*a33 + 6.0*v*a34;
    s40 = 2.0*a43 + 6.0*v*a44;
    if (dvv != NULL) dvv[i] = s10 + u*(s20 + u*(s30 + u*s40));
  }
}


static int caps_newton2D(int nrank, int nu, int nv, double *coeff, double *sv, 
                         double *uv, double *tmp)
{
  int    k;
  double dis, dis0 = DBL_MAX, told = TOLCOPO*TOLCOPO;
  double cu, cv, a00, a10, a11, b0, b1, det;

  while ((uv[0] >= 0.0) && (uv[0] <= (double)(nu-1)) &&
         (uv[1] >= 0.0) && (uv[1] <= (double)(nv-1))) {

    caps_eval2D(nrank, nu, nv, coeff, uv, &tmp[0], &tmp[nrank], &tmp[2*nrank], 
                &tmp[3*nrank], &tmp[4*nrank], &tmp[5*nrank]);
    for (dis = 0.0, k = 0; k < nrank; k++) dis += (tmp[k]-sv[k])*(tmp[k]-sv[k]);
    if (dis > dis0) return 1;
    dis0 = dis;
    a00  = a10 = a11 = 0.0;
    b0   = b1 = 0;
    for (k = 0; k < nrank; k++) {
      a00 += tmp[  nrank+k]*tmp[  nrank+k] + (tmp[k]-sv[k])*tmp[3*nrank+k];
      a10 += tmp[  nrank+k]*tmp[2*nrank+k] + (tmp[k]-sv[k])*tmp[4*nrank+k];
      a11 += tmp[2*nrank+k]*tmp[2*nrank+k] + (tmp[k]-sv[k])*tmp[5*nrank+k];
      b0  -= (tmp[k]-sv[k])*tmp[  nrank+k];
      b1  -= (tmp[k]-sv[k])*tmp[2*nrank+k];
    }
    det  = (a00*a11 - a10*a10);
    if (det == 0.0) return 1;
    det = 1.0/det;
    cu  = det*(b0*a11 - b1*a10);
    cv  = det*(b1*a00 - b0*a10);
    for (a00 = 0.0, k = 0; k < nrank; k++)
      a00 += (tmp[  nrank+k]*cu + tmp[2*nrank+k]*cv)*
             (tmp[  nrank+k]*cu + tmp[2*nrank+k]*cv);
    uv[0] += cu;
    uv[1] += cv;
    if (a00 < told) break;
  }

  if ((uv[0] >= 0.0) && (uv[0] <= (double)(nu-1)) &&
      (uv[1] >= 0.0) && (uv[1] <= (double)(nv-1))) return 0;

  return 1;
}


static void caps_invEval2D(int nrank, int nux, int nvx, double *coeff, 
                           double *sv, double *uv, double *tmp)
{
  /* NOTE: tmp must be at least 6*nrank in length */
  int    i, j, k, l, ik, jk, nu, nv;
  double cu, cv, kcu, kcv, jcu, jcv, stepu, stepv, dis, dis0;

  nu   = nux;
  nv   = nvx;
  if (nu < 0) nu = -nu;
  if (nv < 0) nv = -nv;
  cu   = cv = 0.0;
  dis0 = DBL_MAX;

  for (j = 0; j < nv; j++)
    for (i = 0; i < nu; i++) {
      dis = 0.0;
      for (k = 0; k < nrank; k++)
        dis += (coeff[4*nrank*(j*nu+i)+k]-sv[k])*
               (coeff[4*nrank*(j*nu+i)+k]-sv[k]);
      if (dis < dis0) {
         dis0 = dis;
         cu   = (double) i;
         cv   = (double) j;
      }
    }

  stepu = stepv = 1.0;
  for (l = 0; l < 20; l++) {
    uv[0] = cu;
    uv[1] = cv;
    caps_eval2D(nrank, nu, nv, coeff, uv, &tmp[0], &tmp[nrank], &tmp[2*nrank], 
                NULL, NULL, NULL);
    jcu = jcv = 0.0;
    for (k = 0; k < nrank; k++) {
      jcu += tmp[  nrank+k]*tmp[  nrank+k];
      jcv += tmp[2*nrank+k]*tmp[2*nrank+k];
    }
    jcu = sqrt(jcu)*stepu;
    jcv = sqrt(jcv)*stepv;
    if (jcu > 2.0*jcv) {
       stepu = 0.5*stepu;
    } else if (jcv > 2.0*jcu) {
       stepv = 0.5*stepv;
    } else {
       stepu = 0.5*stepu;
       stepv = 0.5*stepv;
    }

again:
    kcu = cu;
    kcv = cv;

    for (j = 0, jcv = kcv-4.0*stepv, jk = -1; j < 9; j++, jcv += stepv) {
      for (i = 0, jcu = kcu-4.0*stepu, ik = -1; i < 9; i++, jcu += stepu) {
        /* Avoid already computed locations */
        if (i%2 == 0 && j%2 == 0) continue;
        if ((jcu >= 0.0) && (jcu <= (double) (nu-1)) &&
            (jcv >= 0.0) && (jcv <= (double) (nv-1))) {
          uv[0] = jcu;
          uv[1] = jcv;
          caps_eval2D(nrank, nu, nv, coeff, uv, &tmp[0], 
                      NULL, NULL, NULL, NULL, NULL);
          dis = 0.0;
          for (k = 0; k < nrank; k++) dis += (tmp[k]-sv[k])*(tmp[k]-sv[k]);
          if (dis < dis0) {
            dis0 = dis;
            cu   = jcu;
            cv   = jcv;
            ik   = i;
            jk   = j;
          }
        }
      }
    }
    uv[0] = cu;
    uv[1] = cv;
    if (caps_newton2D(nrank, nu, nv, coeff, sv, uv, tmp) == 0) {
      for (dis = 0.0, k = 0; k < nrank; k++) 
        dis += (tmp[k]-sv[k])*(tmp[k]-sv[k]);
      if (dis < dis0) return;
    }

    if (ik > -1 && (ik < 2 || ik > 6 || jk < 2 || jk > 6)) goto again;
  }
  
  uv[0] = cu;
  uv[1] = cv;
}


static int caps_Aprx1DFit(int nrank, int tr0, int tr1, int npts,
                          /*@null@*/ double *tx, double *values, double tol, 
                          capsAprx1D *interp)
{
  int    i, j, stat, nt, ntx, ntm = 0, periodic = 11;
  double sq, rmserr, maxerr, *r, *coeff, *ts;
  double *fit = NULL, *tfit = NULL, *tmap = NULL;

  if (npts <= 1) return EGADS_INDEXERR;
  ts = (double *) EG_alloc(npts*sizeof(double));
  if (ts == NULL) return EGADS_MALLOC;

  if (npts == 2) {
    nt  = npts;
    fit = (double *) EG_alloc(npts*nrank*sizeof(double));
    if (fit == NULL) {
      EG_free(ts);
      return EGADS_MALLOC;
    }
    if (tx != NULL) {
      ntm  = npts;
      tfit = (double *) EG_alloc(npts*sizeof(double));
      if (tfit == NULL) {
        EG_free(fit);
        EG_free(ts);
        return EGADS_MALLOC;
      }
      for (i = 0; i < npts; i++) tfit[i] = tx[i];
    }
    for (i = 0; i < npts*nrank; i++)  fit[i] = values[i];
    ts[0] = 0.0;
    ts[1] = 1.0;
    goto fill;
  }

  /* are we monotonic in tx? */

  if (tx != NULL) {
    maxerr = tx[1] - tx[0];
    for (i = 2; i < npts; i++) {
      sq = tx[i] - tx[i-1];
      if (sq*maxerr <= 0.0) {
        printf(" caps_Interp1DFit: Ts are not Monotonic!\n");
        EG_free(ts);
        return EGADS_CONSTERR;
      }
      maxerr = sq;
    }
  }

  /* make the fit */

  ts[0] = 0.0;
  for (i = 1; i < npts; i++) {
    for (sq = 0.0, j = tr0; j < tr1; j++)
      sq += (values[nrank*i+j] - values[nrank*(i-1)+j])*
            (values[nrank*i+j] - values[nrank*(i-1)+j]);
    ts[i] = ts[i-1] + sqrt(sq);
  }
  stat = prm_NormalizeU(0.0, periodic, npts, ts);
  if (stat != EGADS_SUCCESS) {
    printf(" caps_Interp1DFit: prm_NormalizeU = %d!\n", stat);
    EG_free(ts);
    return EGADS_NOTFOUND;
  }
  nt   = npts;
  stat = prm_BestCfit(npts, nrank, ts, values, tol, periodic, 
                      &nt, &fit, &rmserr, &maxerr);
#ifdef DEBUG
  printf(" caps_Interp1DFit: prm_BestCfit Vals errors = %d %d  %lf %lf\n",
         stat, nt, rmserr, maxerr);
#endif
  if ((stat != EGADS_SUCCESS) && (stat != PRM_TOLERANCEUNMET)) {
    printf(" caps_Interp1DFit: npts = %d, prm_BestCfit = %d\n", npts, stat);
    if (fit != NULL) EG_free(fit);
    EG_free(ts);
    return EGADS_NOTFOUND;
  }
  if (fit == NULL) {
    printf(" caps_Interp1DFit: prm_BestCfit returns NULL\n");
    EG_free(ts);
    return EGADS_NULLOBJ;
  }
  if (tx != NULL) {
    ntm  = 1.5*npts;
    stat = prm_BestCfit(npts, 1, ts, tx, tol, periodic, 
                        &ntm, &tfit, &rmserr, &maxerr);
#ifdef DEBUG
    printf(" caps_Interp1DFit: prm_BestCfit Ts   errors = %d %d  %lf %lf\n",
           stat, ntm, rmserr, maxerr);
#endif
    if ((stat != EGADS_SUCCESS) && (stat != PRM_TOLERANCEUNMET)) {
      printf(" caps_Interp1DFit: Ts npts = %d, prm_BestCfit = %d", npts, stat);
      if (tfit != NULL) EG_free(tfit);
      EG_free(fit);
      EG_free(ts);
      return EGADS_NOTFOUND;
    }
    if (tfit == NULL) {
      printf(" caps_Interp1DFit: prm_BestCfit returns NULL for Ts");
      EG_free(fit);
      EG_free(ts);
      return EGADS_NOTFOUND;
    }
  }

  /* fill in the structure */

fill:
  if (nt > ntm) {
    i = nt;
  } else {
    i = ntm;
  }
  r = (double *) EG_alloc(6*i*sizeof(double));
  if (r == NULL) {
    if (tfit != NULL) EG_free(tfit);
    EG_free(fit);
    EG_free(ts);
    return EGADS_MALLOC;
  }
  coeff = (double *) EG_alloc(nrank*2*nt*sizeof(double));
  if (coeff == NULL) {
    EG_free(r);
    if (tfit != NULL) EG_free(tfit);
    EG_free(fit);
    EG_free(ts);
    return EGADS_MALLOC;
  }
  if (tfit != NULL) {
    tmap = (double *) EG_alloc(2*ntm*sizeof(double));
    if (tmap == NULL) {
      EG_free(tmap);
      EG_free(r);
      EG_free(coeff);
      EG_free(tfit);
      EG_free(fit);
      EG_free(ts);
      return EGADS_MALLOC;
    }
  }
  ntx = nt;
  if (periodic == 1) ntx = -nt;
  if (caps_fillCoeff1D(nrank, ntx, fit, coeff, r) == 1) {
    EG_free(tmap);
    EG_free(r);
    EG_free(coeff);
    if (tfit != NULL) EG_free(tfit);
    EG_free(fit);
    EG_free(ts);
    return EGADS_NULLOBJ;
  }
  if ((tfit != NULL) && (tmap != NULL))
    if (caps_fillCoeff1D(1, ntm, tfit, tmap, r) == 1) {
      EG_free(tmap);
      EG_free(r);
      EG_free(coeff);
      EG_free(tfit);
      EG_free(fit);
      EG_free(ts);
      return EGADS_NULLOBJ;
   }
  EG_free(r);
  if (tfit != NULL) EG_free(tfit);
  EG_free(fit);

  interp->nrank    = nrank;
  interp->periodic = periodic;
  interp->nts      = nt;
  interp->interp   = coeff;
  interp->ntm      = ntm;
  interp->tmap     = tmap;
  if (tx == NULL) {
    interp->trange[0] = 0.0;
    interp->trange[1] = nt-1;
  } else {
    interp->trange[0] = tx[0];
    interp->trange[1] = tx[npts-1];
  }
  EG_free(ts);

  return EGADS_SUCCESS;
}


int caps_Interp1DFit(int nrank, int npts, /*@null@*/ double *tx, double *values,
                     double tol, capsAprx1D *interp)
{
  return caps_Aprx1DFit(nrank, 0, nrank, npts, tx, values, tol, interp);
}


static int
caps_triFill(int npts, int ntris, int *tris, prmTri *vtris)
{
  int     i, j, *vtab, nside;
  connect *etab;
  
  vtab = (int *) EG_alloc(npts*sizeof(int));
  if (vtab == NULL) {
    printf(" CAPS Error: Vert Table Malloc (caps_triFill)!\n");
    return EGADS_MALLOC;    
  }
  etab = (connect *) EG_alloc(ntris*3*sizeof(connect));
  if (etab == NULL) {
    printf(" CAPS Error: Side Table Malloc (caps_triFill)!\n");
    EG_free(vtab);
    return EGADS_MALLOC;    
  }

  for (i = 0; i < ntris; i++) {
    vtris[i].indices[0] = tris[3*i  ];
    vtris[i].indices[1] = tris[3*i+1];
    vtris[i].indices[2] = tris[3*i+2];
    vtris[i].neigh[0]   = i+1;
    vtris[i].neigh[1]   = i+1;
    vtris[i].neigh[2]   = i+1;
    vtris[i].own        = 1;
  }
  nside = NOTFILLED;
  for (j = 0; j < npts; j++) vtab[j] = NOTFILLED;
  for (i = 0; i < ntris;  i++) {
    EG_makeConnect( vtris[i].indices[1], vtris[i].indices[2],
                   &vtris[i].neigh[0], &nside, vtab, etab, 0);
    EG_makeConnect( vtris[i].indices[0], vtris[i].indices[2],
                   &vtris[i].neigh[1], &nside, vtab, etab, 0);
    EG_makeConnect( vtris[i].indices[0], vtris[i].indices[1],
                   &vtris[i].neigh[2], &nside, vtab, etab, 0);
  }
                   
  /* find any unconnected triangle sides */
  for (j = 0; j <= nside; j++) {
    if (etab[j].tri == NULL) continue;
/*  printf(" CAPS Info: Unconnected Side %d %d = %d\n",
           etab[j].node1+1, etab[j].node2+1, *etab[j].tri); */
    *etab[j].tri = 0;
  }

  EG_free(etab);
  EG_free(vtab);
  return EGADS_SUCCESS;
}


int caps_Interp2DFit(int nrank, int npts, double *uvx, double *values, 
                     int ntris, int *tris, double tol, capsAprx2D *interp)
{
  int    i, j, k, nneg, maxsize, stat, nux, nvx, nu, num, nv = 0, nvm = 0;
  int    iv0, iv1, iv2, periodic;
  double rmserr, maxerr, dotmin, area, *r, *coeff, *uvmap;
  double ll[2], ur[2], *fit = NULL, *uvfit = NULL;
  prmTri *vtris;
  prmUV  *uvs;
 
  if ((npts == 0) || (ntris == 0)) return EGADS_RANGERR;
  
  uvs   = (prmUV *)  EG_alloc( npts*sizeof(prmUV));
  vtris = (prmTri *) EG_alloc(ntris*sizeof(prmTri));
  if ((uvs == NULL) || (vtris == NULL)) {
    if (vtris != NULL) EG_free(vtris);
    if (uvs   != NULL) EG_free(uvs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < npts; i++) {
    uvs[i].u = uvx[2*i  ];
    uvs[i].v = uvx[2*i+1];
  }
  stat = caps_triFill(npts, ntris, tris, vtris);
  if (stat != EGADS_SUCCESS) {
    printf(" caps_Interp2DFit: caps_triFill = %d!\n", stat);
    EG_free(vtris);
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }
  for (nneg = i = 0; i < ntris; i++) {
    iv0  = vtris[i].indices[0] - 1;
    iv1  = vtris[i].indices[1] - 1;
    iv2  = vtris[i].indices[2] - 1;
    area = (uvs[iv1].u - uvs[iv0].u)*(uvs[iv2].v - uvs[iv0].v) -
           (uvs[iv1].v - uvs[iv0].v)*(uvs[iv2].u - uvs[iv0].u);
    if (area < 0.0) nneg++;
  }
  if (nneg == ntris) {
    for (i = 0; i < ntris; i++) {
      iv0                 = vtris[i].indices[0];
      vtris[i].indices[0] = vtris[i].indices[1];
      vtris[i].indices[1] = iv0;
      iv0                 = vtris[i].neigh[0];
      vtris[i].neigh[0]   = vtris[i].neigh[1];
      vtris[i].neigh[1]   = iv0;
    }
  } else if (nneg != 0) {
    printf(" caps_Interp2DFit: Input has %d Negative (UV) tris (of %d)!\n",
           nneg, ntris);
    EG_free(vtris);
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }

  /* make the fit */

  periodic = 0;
  stat = prm_SmoothUV(2, periodic, NULL, ntris, vtris, npts, nrank, uvs, values);
  EG_free(vtris);
  if ((stat != EGADS_SUCCESS) && (stat != PRM_NOTCONVERGED)) {
    printf(" caps_Interp2DFit: prm_SmoothUV = %d!\n", stat);
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }
//#ifdef DEBUG
  if (stat == PRM_NOTCONVERGED)
    printf(" caps_Interp2DFit: prm_SmoothUV Not Converged!\n");
//#endif

  stat = prm_NormalizeUV(0.0, periodic, npts, uvs);
  if (stat != EGADS_SUCCESS) {
    printf(" caps_Interp2DFit: prm_NormalizeUV = %d!\n", stat);
    EG_free(uvs);
    return EGADS_NOTFOUND;
  }
  nu   = npts;
  stat = prm_BestGrid(npts, nrank, uvs, values, 0, NULL, tol, periodic, NULL,
                      &nu, &nv, &fit, &rmserr, &maxerr, &dotmin);
#ifdef DEBUG
  printf(" caps_Interp2DFit: prm_BestGrid Values errors = %d  %d %d  %lf %lf\n",
         stat, nu, nv, rmserr, maxerr);
#endif
  if ((stat != EGADS_SUCCESS) && (stat != PRM_TOLERANCEUNMET)) {
    printf(" caps_Interp2DFit: prm_BestGrid = %d!\n", stat);
    if (fit != NULL) EG_free(fit);
    EG_free(uvs);
    return EGADS_NOTFOUND;
  }
  if (fit == NULL) {
    printf(" caps_Interp2DFit: prm_BestGrid returns NULL!\n");
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }
  num  = 1.5*npts;
  stat = prm_BestGrid(npts, 2, uvs, uvx, 0, NULL, tol, periodic,  NULL,
                      &num, &nvm, &uvfit, &rmserr, &maxerr, &dotmin);
#ifdef DEBUG
  printf(" caps_Interp2DFit: prm_BestGrid UVs    errors = %d  %d %d  %lf %lf\n",
         stat, num, nvm, rmserr, maxerr);
#endif
  if ((stat != EGADS_SUCCESS) && (stat != PRM_TOLERANCEUNMET)) {
    printf(" caps_Interp2DFit: prm_BestGrid UV = %d!\n", stat);
    if (uvfit != NULL) EG_free(uvfit);
    EG_free(fit);
    EG_free(uvs);
    return EGADS_NOTFOUND;
  }
  if (uvfit == NULL) {
    printf(" caps_Interp2DFit: prm_BestGrid UV returns NULL!\n");
    EG_free(fit);
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }
  ll[0] = ur[0] = uvfit[0];
  ll[1] = ur[1] = uvfit[1];
  for (k = j = 0; j < nvm; j++)
    for (i = 0; i < num; i++, k++) {
      if (uvfit[2*k  ] < ll[0]) ll[0] = uvfit[2*k  ];
      if (uvfit[2*k  ] > ur[0]) ur[0] = uvfit[2*k  ];
      if (uvfit[2*k+1] < ll[1]) ll[1] = uvfit[2*k+1];
      if (uvfit[2*k+1] > ur[1]) ur[1] = uvfit[2*k+1];
/*	printf(" %d,%d: %lf %lf\n", j, i, uvfit[2*k], uvfit[2*k+1]);  */
    }

  /* fill in the structure */

  maxsize = nu;
  if (maxsize < nv)  maxsize = nv;
  if (maxsize < num) maxsize = num;
  if (maxsize < nvm) maxsize = nvm;
  coeff = (double *) EG_alloc(nrank*4*nu*nv*sizeof(double));
  if (coeff == NULL) {
    EG_free(uvfit);
    EG_free(fit);
    EG_free(uvs);
    return EGADS_MALLOC;
  }
  uvmap = (double *) EG_alloc(2*4*num*nvm*sizeof(double));
  if (uvmap == NULL) {
    EG_free(coeff);
    EG_free(uvfit);
    EG_free(fit);
    EG_free(uvs);
    return EGADS_MALLOC;
  }
  r = (double *) EG_alloc(6*maxsize*sizeof(double));
  if (r == NULL) {
    EG_free(uvmap);
    EG_free(coeff);
    EG_free(uvfit);
    EG_free(fit);
    EG_free(uvs);
    return EGADS_MALLOC;
  }
  nux = nu;
  nvx = nv;
  if ((periodic & 1) != 0) nux = -nu;
  if ((periodic & 2) != 0) nvx = -nv;
  if (caps_fillCoeff2D(nrank, nux, nvx, fit, coeff, r) == 1) {
    EG_free(r);
    EG_free(uvmap);
    EG_free(coeff);
    EG_free(uvfit);
    EG_free(fit);
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }
  if (caps_fillCoeff2D(2, num, nvm, uvfit, uvmap, r) == 1) {
    EG_free(r);
    EG_free(uvmap);
    EG_free(coeff);
    EG_free(uvfit);
    EG_free(fit);
    EG_free(uvs);
    return EGADS_NULLOBJ;
  }
  EG_free(r);
  EG_free(uvfit);
  EG_free(fit);

  interp->nrank     = nrank;
  interp->periodic  = periodic;
  interp->nus       = nu;
  interp->nvs       = nv;
  interp->interp    = coeff;
  interp->urange[0] = ll[0];
  interp->urange[1] = ur[0];
  interp->vrange[0] = ll[1];
  interp->vrange[1] = ur[1];
  interp->num       = num;
  interp->nvm       = nvm;
  interp->uvmap     = uvmap;
  
  EG_free(uvs);
  return EGADS_SUCCESS;
}


int caps_Interpolate1D(capsAprx1D *interp, double tx, double *sv,
                       /*@null@*/ double *dt1, /*@null@*/ double *dt2)
{
  int    i, nrank;
  double t, mt0, mt1, r = 1.0;

  nrank = interp->nrank;
  t     = tx;
  if (interp->tmap != NULL) {
    caps_invEval1D(1, interp->ntm, interp->tmap, &tx, &t);
    if ((dt1 != NULL) || (dt2 != NULL)) {
      caps_eval1D(1, interp->ntm, interp->tmap, t, &mt0, &mt1, NULL);
      if (mt1 == 0.0) mt1 = 1.0;
      r  = (interp->nts-1)/mt1;
      r /=  interp->ntm-1;
    }
    t *= interp->nts-1;
    t /= interp->ntm-1;
  }

  caps_eval1D(nrank, interp->nts, interp->interp, t, sv, dt1, dt2);
  if (dt1 != NULL) for (i = 0; i < nrank; i++) dt1[i] *= r;
  if (dt2 != NULL) for (i = 0; i < nrank; i++) dt2[i] *= r*r;

  return EGADS_SUCCESS;
}


int caps_Interpolate2D(capsAprx2D *interp, double *uvx, double *sv,
                       /*@null@*/ double *du,  /*@null@*/ double *dv,
                       /*@null@*/ double *duu, /*@null@*/ double *duv,
                       /*@null@*/ double *dvv)
{
  int    i, nrank;
  double uvn[2], uv[2], tmp[12], *store;
#ifdef DIFFERENCE
  double step = 1.e-5;
#else
  double mp[2], mu[2], mv[2], sav, det;
#endif

  nrank = interp->nrank;
  if (interp->uvmap == NULL) {
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uvx, sv,
                du, dv, duu, duv, dvv);
    return EGADS_SUCCESS;
  }

  caps_invEval2D(2, interp->num, interp->nvm, interp->uvmap, uvx, uvn, tmp);
  uv[0] = (interp->nus-1)*uvn[0]/(interp->num-1);
  uv[1] = (interp->nvs-1)*uvn[1]/(interp->nvm-1);
  if ((du == NULL) && (dv == NULL) && (duu == NULL) && 
                                      (duv == NULL) && (dvv == NULL)) {
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv, sv,
               NULL, NULL, NULL, NULL, NULL);
    return EGADS_SUCCESS;
  }

  store = (double *) EG_alloc(nrank*5*sizeof(double));
  if (store == NULL) return EGADS_MALLOC;

#ifndef DIFFERENCE

  caps_eval2D(2, interp->num, interp->nvm, interp->uvmap, uvn, mp, mu, mv,
            NULL, NULL, NULL);

  /* invert the mapping matrix */
  mu[0] *= (interp->num-1);
  mu[1] *= (interp->num-1);
  mv[0] *= (interp->nvm-1);
  mv[1] *= (interp->nvm-1);
  det = mu[0]*mv[1] - mu[1]*mv[0];
  if (det != 0.0) det = 1.0/det;
  sav    =  mu[0];
  mu[0]  =  det*mv[1]*(interp->nus-1);
  mu[1] *= -det      *(interp->nus-1)/(interp->nvs-1);
  mv[0] *= -det      *(interp->nvs-1)/(interp->nus-1);
  mv[1]  =  det*sav  *(interp->nvs-1);

  caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv, sv, 
              &store[0], &store[nrank], &store[2*nrank], &store[3*nrank], 
              &store[4*nrank]);

  if (du != NULL)
   for (i = 0; i < nrank; i++)
     du[i] = store[i]*mu[0] + store[nrank+i]*mu[1];

  if (dv != NULL)
   for (i = 0; i < nrank; i++) 
     dv[i] = store[i]*mv[0] + store[nrank+i]*mv[1];

  if (duu != NULL)
   for (i = 0; i < nrank; i++) 
     duu[i] = store[2*nrank+i]*mu[0]*mu[0] + store[3*nrank+i]*mu[0]*mu[1] +
                                             store[4*nrank+i]*mu[1]*mu[1];

  if (duv != NULL)
   for (i = 0; i < nrank; i++) 
     duv[i] = store[2*nrank+i]*mu[0]*mv[0] + store[3*nrank+i]*mu[0]*mv[1] +
              store[3*nrank+i]*mu[1]*mv[0] + store[4*nrank+i]*mu[1]*mv[1];

  if (dvv != NULL)
   for (i = 0; i < nrank; i++) 
     dvv[i] = store[2*nrank+i]*mv[0]*mv[0] + store[3*nrank+i]*mv[0]*mv[1] + 
                                             store[4*nrank+i]*mv[1]*mv[1];

#else

  caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv, sv,
              NULL, NULL, NULL, NULL, NULL);

  if ((du != NULL) || (duu != NULL) || (duv != NULL)) {
    uv[0] = uvx[0] - step;
    uv[1] = uvx[1];
    caps_invEval2D(2, interp->num, interp->nvm, interp->uvmap, uv, uvn, tmp);
    uv[0] = (interp->nus-1)*uvn[0]/(interp->num-1);
    uv[1] = (interp->nvs-1)*uvn[1]/(interp->nvm-1);
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv, &store[0],
                NULL, NULL, NULL, NULL, NULL);
    uv[0] = uvx[0] + step;
    uv[1] = uvx[1];
    caps_invEval2D(2, interp->num, interp->nvm, interp->uvmap, uv, uvn, tmp);
    uv[0] = (interp->nus-1)*uvn[0]/(interp->num-1);
    uv[1] = (interp->nvs-1)*uvn[1]/(interp->nvm-1);
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv,
                &store[nrank], NULL, NULL, NULL, NULL, NULL);
    if (du != NULL)
      for (i = 0; i < nrank; i++) 
        du[i] = (store[nrank+i]-store[i]) / (2.0*step);
    if (duu != NULL)
      for (i = 0; i < nrank; i++) {
        store[      i] = (sv[i]         -store[i]) / step;
        store[nrank+i] = (store[nrank+i]-   sv[i]) / step;
        duu[i]         = (store[nrank+i]-store[i]) / step;
      }
  }

  if ((dv != NULL) || (duv != NULL) || (dvv != NULL)) {
    uv[0] = uvx[0];
    uv[1] = uvx[1] - step;
    caps_invEval2D(2, interp->num, interp->nvm, interp->uvmap, uv, uvn, tmp);
    uv[0] = (interp->nus-1)*uvn[0]/(interp->num-1);
    uv[1] = (interp->nvs-1)*uvn[1]/(interp->nvm-1);
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv,
                &store[2*nrank], NULL, NULL, NULL, NULL, NULL);
    uv[0] = uvx[0];
    uv[1] = uvx[1] + step;
    caps_invEval2D(2, interp->num, interp->nvm, interp->uvmap, uv, uvn, tmp);
    uv[0] = (interp->nus-1)*uvn[0]/(interp->num-1);
    uv[1] = (interp->nvs-1)*uvn[1]/(interp->nvm-1);
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv,
                &store[3*nrank], NULL, NULL, NULL, NULL, NULL);
    if (dv != NULL)
      for (i = 0; i < nrank; i++) 
        dv[i] = (store[3*nrank+i]-store[2*nrank+i]) / (2.0*step);
    if (dvv != NULL)
      for (i = 0; i < nrank; i++) {
        store[2*nrank+i] = (sv[i]           -store[2*nrank+i]) / step;
        store[3*nrank+i] = (store[3*nrank+i]-           sv[i]) / step;
        dvv[i]           = (store[3*nrank+i]-store[2*nrank+i]) / step;
      }
  }

  if (duv != NULL) {
    uv[0] = uvx[0] + step;
    uv[1] = uvx[1] + step;
    caps_invEval2D(2, interp->num, interp->nvm, interp->uvmap, uv, uvn, tmp);
    uv[0] = (interp->nus-1)*uvn[0]/(interp->num-1);
    uv[1] = (interp->nvs-1)*uvn[1]/(interp->nvm-1);
    caps_eval2D(nrank, interp->nus, interp->nvs, interp->interp, uv,
                &store[4*nrank], NULL, NULL, NULL, NULL, NULL);
    for (i = 0; i < nrank; i++)
      duv[i] = ((store[4*nrank+i]-sv[i])/step - store[  nrank+i] -
                                                store[3*nrank+i]) / step;
  }
#endif

  EG_free(store);
  return EGADS_SUCCESS;
}


int caps_Aprx1DFree(/*@only@*/ capsAprx1D *approx)
{
  if (approx == NULL) return EGADS_SUCCESS;
  if (approx->interp != NULL) EG_free(approx->interp);
  if (approx->tmap   != NULL) EG_free(approx->tmap);
  approx->interp = NULL;
  approx->tmap   = NULL;
  EG_free(approx);

  return EGADS_SUCCESS;
}


int caps_Aprx2DFree(/*@only@*/ capsAprx2D *approx)
{
  if (approx == NULL) return EGADS_SUCCESS;
  if (approx->interp != NULL) EG_free(approx->interp);
  if (approx->uvmap  != NULL) EG_free(approx->uvmap);
  approx->interp = NULL;
  approx->uvmap  = NULL;
  EG_free(approx);

  return EGADS_SUCCESS;
}


int caps_invInterpolate1D(capsAprx1D *interp, double *sv, double *t)
{
  double mt, ntx;

  ntx = interp->nts;
  if (interp->periodic != 0) ntx = -ntx;
  caps_invEval1D(interp->nrank, ntx, interp->interp, sv, t);
  if (interp->tmap != NULL) {
    mt  = *t*(interp->ntm-1);
    mt /=     interp->nts-1;
    caps_eval1D(1, interp->ntm, interp->tmap, mt, t, NULL, NULL);
  }
  return caps_Interpolate1D(interp, *t, sv, NULL, NULL);
}


int caps_invInterpolate2D(capsAprx2D *interp, double *sv, double *uv)
{
  int    nrank, nux, nvx;
  double uvx[2], *tmp;

  nrank = interp->nrank;
  tmp   = (double *) EG_alloc(6*nrank*sizeof(double));
  if (tmp == NULL) return EGADS_MALLOC;
  nux   = interp->nus;
  nvx   = interp->nvs;
  if ((interp->periodic & 1) != 0) nux = -nux;
  if ((interp->periodic & 2) != 0) nvx = -nvx;
  caps_invEval2D(nrank, nux, nvx, interp->interp, sv, uv, tmp);
  EG_free(tmp);
  if (interp->uvmap != NULL) {
    uvx[0]  = uv[0]*(interp->num-1);
    uvx[0] /=        interp->nus-1;
    uvx[1]  = uv[1]*(interp->nvm-1);
    uvx[1] /=        interp->nvs-1;
    caps_eval2D(2, interp->num, interp->nvm, interp->uvmap, uvx, uv, NULL, NULL,
                NULL, NULL, NULL);
  }
  return caps_Interpolate2D(interp, uv, sv, NULL, NULL, NULL, NULL, NULL);
}
