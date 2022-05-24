/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             (unconstrained) conjugate gradient optimization
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egadsErrors.h"

#define MAX(A,B)     (((A) < (B)) ? (B) : (A))


/*
 * golden - golden section search in the grad direction
 */

static int
golden(int    (*func)(int, double[], void *, double *, /*@null@*/ double[]),
                                        /* (in)  name of objective function */
       void   *blind,                   /* (in)  blind pointer to structure */
       int    n,                        /* (in)  number of design variables */
       double xbase[],                  /* (in)  base design variables */
       double dx[],                     /* (in)  change in design variables */
       double xtol,                     /* (in)  convergence tolerance in x */
       double *xmin,                    /* (out) size of step */
       double *fmin)                    /* (out) objective function after step */
{
  int    status = EGADS_SUCCESS;        /* (out) return status */
  
  int    j;
  double x0, f0, x1, f1, x2, f2, x3, f3, u, fu;
  double r, q, ulim, swap;
  double *temp;
  double GOLD   = 1.618034;
  double R      = GOLD - 1;
  double GLIMIT = 100;
  double TINY   = 1e-20;
  
  /* default returns */
  *xmin = 0;
  *fmin = 0;
  
  /* allocate storage */
  temp = (double *) malloc(n*sizeof(double));
  if (temp == NULL) {status = EGADS_MALLOC; goto cleanup;}
  
  /* find (x0,x1,x3) that bracket minimum.  that is, f1 < min(f0,f3) */
  x0 = 0.0;
  for (j = 0; j < n; j++) {
    temp[j] = xbase[j] + x0 * dx[j];
  }
  status = (*func)(n, temp, blind, &f0, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  
  x1 = 0.1;
  for (j = 0; j < n; j++) {
    temp[j] = xbase[j] + x1 * dx[j];
  }
  status = (*func)(n, temp, blind, &f1, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  
  /* switch x0 and x1 so that we go downhill from x0 to x1 */
  if (f1 > f0) {
    swap = x0;  x0  = x1;  x1  = swap;
    swap = f1;  f1  = f0;  f0  = swap;
  }
  
  /* first guess for x3 */
  x3 = x1 + GOLD * (x1-x0);
  
  for (j = 0; j < n; j++) {
    temp[j] = xbase[j] + x3 * dx[j];
  }
  status = (*func)(n, temp, blind, &f3, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  
  /* loop until f1 < min(f0,f3) */
  while (f1 >= f3) {
    
    /* compute u by parabolic extrapolation */
    r = (x1 - x0) * (f1 - f3);
    q = (x1 - x3) * (f1 - f0);
    if (q >= r) {
      u = x1 - ((x1-x3) * q - (x1-x0)*r) / (2.0 * MAX(q-r, TINY));
    } else {
      u = x1 + ((x1-x3) * q - (x1-x0)*r) / (2.0 * MAX(r-q, TINY));
    }
    ulim = x1 + GLIMIT * (x3-x1);
    
    /* parabolic u is between x1 anc x3 */
    if ((x1-u)*(u-x3) > 0) {
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + u * dx[j];
      }
      status = (*func)(n, temp, blind, &fu, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      /* got a minimum between x1 and x3 */
      if (fu < f3) {
        x0 = x1;
        x1 = u;
        break;
        
        /* got a minimum between x0 and u */
      } else if (fu > f1) {
        x3 = u;
        break;
      }
      
      /* parabolic fit was no use, so use default magnification */
      u  = x3 + GOLD * (x3-x1);
      
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + u * dx[j];
      }
      status = (*func)(n, temp, blind, &fu, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      /* parabolic fit is between x3 and its allowed limit */
    } else if ((x3-u)*(u-ulim) > 0) {
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + u * dx[j];
      }
      status = (*func)(n, temp, blind, &fu, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      if (fu < f3) {
        x1 = x3;   f1 = f3;
        x3 = u;    f3 = fu;
        u  = x3 + GOLD * (x3-x1);
        
        for ( j = 0; j < n; j++) {
          temp[j] = xbase[j] + u * dx[j];
        }
        status = (*func)(n, temp, blind, &fu, NULL);
        if (status != EGADS_SUCCESS) goto cleanup;
      }
      
      /* limit parabolic u to maximum allowed value */
    } else if ((u-ulim)*(ulim-x3) >= 0) {
      u  = ulim;
      
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + u * dx[j];
      }
      status = (*func)(n, temp, blind, &fu, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      /* reject parabolic u, use default magnification */
    } else {
      u = x3 + GOLD * (x3-x1);
      
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + u * dx[j];
      }
      status = (*func)(n, temp, blind, &fu, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
    }
    
    /* eliminate oldest point and continue */
    x0 = x1;  f0 = f1;
    x1 = x3;  f1 = f3;
    x3 = u;   f3 = fu;
  }
  
  /* create a point x2 such that x0 to x1 is the smaller segment */
  if (fabs(x3-x1) > (x1-x0)) {
    x2 = x1 + (1-R) * (x1-x0);
    
    for (j = 0; j < n; j++) {
      temp[j] = xbase[j] + x2 * dx[j];
    }
    status = (*func)(n, temp, blind, &f2, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;
  } else {
    x2 = x1;   f2 = f1;
    
    x1 = x1 - (1-R) * (x1-x0);
    
    for (j = 0; j < n; j++) {
      temp[j] = xbase[j] + x1 * dx[j];
    }
    status = (*func)(n, temp, blind, &f1, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;
  }
  
  /* keep shrinking the intervals until the minimum is found */
  while (fabs(x3-x0) > xtol*(fabs(x1)+fabs(x2))) {
    
    /* f2<f1, so move x0 and x1 down to make room for new x2 */
    if (f2 < f1) {
      x0 = x1;
      x1 = x2;   f1 = f2;
      x2 = R * x1 + (1-R) * x3;
      
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + x2 * dx[j];
      }
      status = (*func)(n, temp, blind, &f2, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      /* otherwise, move x1 and x2 up to make room for new x1 */
    } else {
      x3 = x2;
      x2 = x1;  f2 = f1;
      x1 = R * x2 + (1-R) * x0;
      
      for (j = 0; j < n; j++) {
        temp[j] = xbase[j] + x1 * dx[j];
      }
      status = (*func)(n, temp, blind, &f1, NULL);
      if (status != EGADS_SUCCESS) goto cleanup;
    }
  }
  
  /* output the smaller of f1 and f2 */
  if (f1 < f2) {
    *xmin = x1;
    *fmin = f1;
  } else {
    *xmin = x2;
    *fmin = f2;
  }
  
cleanup:
  if (temp != NULL) free(temp);
  
  return(status);
}


int
caps_conjGrad(int (*func)(int, double[], void *, double *, /*@null@*/ double[]),
                                       /* (in)  name of objective function */
              void   *blind,           /* (in)  blind pointer to structure */
              int    n,                /* (in)  number of variables */
              double x[],              /* (in)  initial set of variables */
                                       /* (out) optimized variables */
              double ftol,             /* (in)  convergence tolerance on objective
                                                function */
              /*@null@*/
              FILE   *fp,              /* (in)  unit to write iteration history */
              double *fopt)            /* (out) optimized objective function */
{
    int    status = EGADS_SUCCESS;

    int    j, iter;
    double fx, gg, hh, gamma, xmin;

    double *g    = NULL;
    double *h    = NULL;
    double *grad = NULL;

    int    ITMAX = 200;
    double EPS   = 1e-10;

    /* default returns */
    *fopt = 0;

    /* allocate storage */
    g    = (double*) malloc(n*sizeof(double));
    if (g    == NULL) {status = EGADS_MALLOC; goto cleanup;}

    h    = (double*) malloc(n*sizeof(double));
    if (h    == NULL) {status = EGADS_MALLOC; goto cleanup;}

    grad = (double*) malloc(n*sizeof(double));
    if (grad == NULL) {status = EGADS_MALLOC; goto cleanup;}

    /* initializations */
    status = (*func)(n, x, blind, &fx, grad);
    if (status != EGADS_SUCCESS) goto cleanup;

    *fopt = fx;

    if (fp != NULL) {
        fprintf(fp, "conjGrad Info: iter %3d, fopt=%12.4e\n", 0, *fopt);
    }
#ifdef DEBUG
    printf(" conjGrad Info: iter %3d, fopt=%12.4e\n", 0, *fopt);
#endif

    for (j = 0; j < n; j++) {
        g[   j] = -grad[j];
        h[   j] =     g[j];
        grad[j] =     h[j];
    }

    /* main optimization loop */
    for (iter = 1; iter <= ITMAX; iter++) {

        /* find the minimum along a line in the grad direction */
        status = golden(func, blind, n, x, grad, ftol, &xmin, fopt);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* update point at end of step */
        for (j = 0; j < n; j++) {
            grad[j] *= xmin;
            x[j]    += grad[j];
        }

        if (fp != NULL) {
            fprintf(fp, "conjGrad Info: iter %3d, fopt=%12.4e\n", iter, *fopt);
        }
#ifdef DEBUG
        printf(" conjGrad Info: iter %3d, fopt=%12.4e\n", iter, *fopt);
#endif

        /* if function change was very small, we are done */
        if (fabs(*fopt-fx) <= ftol/2*fabs(*fopt+fabs(fx)+EPS)) {
            if (fp != NULL) {
                fprintf(fp, "conjGrad Info: small function change!\n");
            }
#ifdef DEBUG
            printf(" conjGrad Info: small function change!\n");
#endif
            goto cleanup;
        }

        /* evaluate solution at the "line-minimum" location */
        status = (*func)(n, x, blind, &fx, grad);
        if (status != EGADS_SUCCESS) goto cleanup;

        /* Polak-Ribiere conjugate gradient method */
        gg = 0;
        hh = 0;
        for (j = 0; j < n; j++) {
            gg += g[j] * g[j];
            hh += (grad[j] + g[j]) * grad[j];
        }

        /* if gradient is zero, then we are at the (local) optimum */
        if (gg == 0) {
            if (fp != NULL) {
                fprintf(fp, "conjGrad Info: zero gradient!\n");
            }
#ifdef DEBUG
            printf(" conjGrad Info: zero gradient!\n");
#endif
            goto cleanup;
        }

        /* next direction is conjugate to those we hav already taken */
        gamma = hh / gg;
        for (j = 0; j < n; j++) {
            g[   j] = -grad[j];
            h[   j] =     g[j] + gamma * h[j];
            grad[j] =     h[j];
        }
    }

    if (fp != NULL) {
        fprintf(fp, "conjGrad Info: exceeded maxiter!\n");
    }
#ifdef DEBUG
    printf(" conjGrad Info: exceeded maxiter!\n");
#endif

cleanup:
    if (g    != NULL) free(g   );
    if (h    != NULL) free(h   );
    if (grad != NULL) free(grad);

    return(status);
}
