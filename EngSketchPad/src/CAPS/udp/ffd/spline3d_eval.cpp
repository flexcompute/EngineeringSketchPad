#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "spline3d_eval.hpp"

#ifdef WIN32
#define DllExport   __declspec( dllexport )
#else
#define DllExport
#endif

#define MIN(a,b)        (((a) < (b)) ? (a) : (b))


/* 3D BSplines are stored in a integer header and a double block of memory

   header 0 -- degree U
          1 -- # CP   U
          2 -- # Knot U
          3 -- degree V
          4 -- # CP   V
          5 -- # Knot V
          6 -- degree W
          7 -- # CP   W
          8 -- # Knot W

   data:  knots U, knots V, knots W (#KnotU + #KnotV + #KnotW in length)
          CPs -- 3*(#CPU * #CPV * #CPW) in length

   restrictions: No NURBS only BSplines */

template<class T, class T2>
static int FindSpan(int nKnots, int degree, T2 u, const T *U)
{
  int n, low, mid, high;

  if (u <= U[degree]) return degree;
  n = nKnots - degree - 1;
  if (u >= U[n]) return n-1;

  low  = degree;
  high = n;
  mid  = (low+high)/2;
  while ((u < U[mid]) || (u >= U[mid+1])) {
    if (u < U[mid]) {
      high = mid;
    } else {
      low  = mid;
    }
    mid = (low+high)/2;
  }

  return mid;
}

template<class T, class T2>
static void DersBasisFuns(int i, int p, T2 u, const T *knot, int der, T **ders)
{
  int    j, k, j1, j2, r, s1, s2, rk, pk;
  T d, saved, temp, ndu[SPLINE3D_MAXDEG+1][SPLINE3D_MAXDEG+1];
  T a[2][SPLINE3D_MAXDEG+1], left[SPLINE3D_MAXDEG+1], right[SPLINE3D_MAXDEG+1];

  ndu[0][0] = 1.0;
  for (j = 1; j <= p; j++) {
    left[j]  = u - knot[i+1-j];
    right[j] = knot[i+j] - u;
    saved = 0.0;
    for (r = 0; r < j; r++) {
      ndu[j][r] = right[r+1] + left[j-r];
      temp      = ndu[r][j-1]/ndu[j][r];
      ndu[r][j] = saved + right[r+1]*temp;
      saved     = left[j-r]*temp;
    }
    ndu[j][j] = saved;
  }

  for (j = 0; j <= p; j++) ders[0][j] = ndu[j][p];  /* basis function */

  /* compute derivatives */
  for (r = 0; r <= p; r++ ) {
    s1      = 0;
    s2      = 1;
    a[0][0] = 1.0;
    /* compute k'th derivative */
    for (k = 1; k <= der; k++ ) {
      d  = 0.0;
      rk = r - k;
      pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0]/ndu[pk+1][rk];
        d        = a[s2][0]*ndu[rk][pk];
      }
      j1 = rk >= -1  ? 1   : -rk;
      j2 = (r-1<=pk) ? k-1 : p-r;
      for (j = j1; j <= j2; j++) {
        a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
        d       +=  a[s2][j]*ndu[rk+j][pk];
      }
      if (r <= pk) {
        a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
        d       +=  a[s2][k]*ndu[r][pk];
      }
      ders[k][r] = d;
      /* switch rows */
      j  = s1;
      s1 = s2;
      s2 = j;
    }
  }

  r = p;
  for (k = 1; k <= der; k++) {
    for (j = 0; j <= p; j++ ) ders[k][j] *= r;
    r *= p - k;
  }

}


/*
 * evaluates the 3D BSpline, 1st and 2nd derivatives
 *       where: ivec  - integer BSpline data
 *              data  - double BSpline data
 *              uvw   - BSpline coordinates for the evalation (3 in len)
 *              deriv - returned position and derivatives
 *                       0 -  2  position
 *                       3 -  5  Du
 *                       6 -  8  Du2
 *                       9 - 11  Dv
 *                      12 - 14  DuDv
 *                      15 - 17  Dv2
 *                      18 - 20  Dw
 *                      21 - 23  DuDw
 *                      24 - 26  DvDw
 *                      27 - 29  Dw2
 *
 *              returns 0 - OK, -1 degree is too big
 *
 */
template<class T, class T2>
static int
spline3dEval_impl( const int *ivec, const T *data, const T2 *uvw, T *deriv)
{
  int    der = 2;
  int    degu, degv, degw, nKu, nKv, nKw, nCPu, nCPv;
  int    i, j, k, l, m, n, s, t, spanu, spanv, spanw, du, dv, dw;
  T *NderU[SPLINE3D_MAXDEG+1], *NderV[SPLINE3D_MAXDEG+1], *NderW[SPLINE3D_MAXDEG+1];
  T Nu[SPLINE3D_MAXDEG+1][SPLINE3D_MAXDEG+1], Nv[SPLINE3D_MAXDEG+1][SPLINE3D_MAXDEG+1];
  T Nw[SPLINE3D_MAXDEG+1][SPLINE3D_MAXDEG+1], temp[3*SPLINE3D_MAXDEG], tmp[3*SPLINE3D_MAXDEG];
  const T *Ku, *Kv, *Kw, *CP;

  degu = ivec[0];
  nCPu = ivec[1];
  nKu  = ivec[2];
  degv = ivec[3];
  nCPv = ivec[4];
  nKv  = ivec[5];
  degw = ivec[6];
//nCPw = ivec[7];
  nKw  = ivec[8];
  Ku   = data;
  Kv   = data + nKu;
  Kw   = data + nKu + nKv;
  CP   = data + nKu + nKv + nKw;
  du   = MIN(der, degu);
  dv   = MIN(der, degv);
  dw   = MIN(der, degw);
  for (m = l = 0; l <= der; l++)
    for (k = 0; k <= der-l; k++)
      for (n = 0; n <= der-l-k; n++, m++)
        deriv[3*m  ] = deriv[3*m+1] = deriv[3*m+2] = 0.0;
  if (degu >= SPLINE3D_MAXDEG) {
    printf(" spline3dEVal: degreeU %d >= %d!\n", degu, SPLINE3D_MAXDEG);
    return -1;
  }
  if (degv >= SPLINE3D_MAXDEG) {
    printf(" spline3dEval: degreeV %d >= %d!\n", degv, SPLINE3D_MAXDEG);
    return -1;
  }
  if (degw >= SPLINE3D_MAXDEG) {
    printf(" spline3dEval: degreeW %d >= %d!\n", degw, SPLINE3D_MAXDEG);
    return -1;
  }
  for (i = 0; i <= degu; i++) NderU[i] = &Nu[i][0];
  for (i = 0; i <= degv; i++) NderV[i] = &Nv[i][0];
  for (i = 0; i <= degw; i++) NderW[i] = &Nw[i][0];

  spanu = FindSpan(nKu, degu, uvw[0], Ku);
  DersBasisFuns(spanu,  degu, uvw[0], Ku, du, NderU);
  spanv = FindSpan(nKv, degv, uvw[1], Kv);
  DersBasisFuns(spanv,  degv, uvw[1], Kv, dv, NderV);
  spanw = FindSpan(nKw, degw, uvw[2], Kw);
  DersBasisFuns(spanw,  degw, uvw[2], Kw, dw, NderW);

  for (m = l = 0; l <= dw; l++) {
    for (k = 0; k <= der-l; k++) {
      for (n = 0; n <= der-l-k; n++, m++) {
        if (k > dv) continue;
        if (n > du) continue;
        for (t = 0; t <= degw; t++) {
          temp[3*t] = temp[3*t+1] = temp[3*t+2] = 0.0;
          for (s = 0; s <= degv; s++) {
            tmp[3*s] = tmp[3*s+1] = tmp[3*s+2] = 0.0;
            for (j = 0; j <= degu; j++) {
              i            = spanu-degu+j +      nCPu*(spanv-degv+s) +
                                            nCPu*nCPv*(spanw-degw+t);
              tmp[3*s  ] += Nu[n][j]*CP[3*i  ];
              tmp[3*s+1] += Nu[n][j]*CP[3*i+1];
              tmp[3*s+2] += Nu[n][j]*CP[3*i+2];
            }
          }
          for (s = 0; s <= degv; s++) {
            temp[3*t  ] += Nv[k][s]*tmp[3*s  ];
            temp[3*t+1] += Nv[k][s]*tmp[3*s+1];
            temp[3*t+2] += Nv[k][s]*tmp[3*s+2];
          }
        }
        for (t = 0; t <= degw; t++) {
          deriv[3*m  ] += Nw[l][t]*temp[3*t  ];
          deriv[3*m+1] += Nw[l][t]*temp[3*t+1];
          deriv[3*m+2] += Nw[l][t]*temp[3*t+2];
        }
      }
    }
  }
  return 0;
}

extern "C"
int
spline3dEval( const int *ivec, const double *data, const double *uvw, double *deriv)
{
  return spline3dEval_impl(ivec, data, uvw, deriv);
}


template<int N, class T>
int
spline3dEval( const int *ivec, const SurrealS<N> *data, const T *uvw,
              SurrealS<N> *deriv)
{
  return spline3dEval_impl(ivec, data, uvw, deriv);
}


// Create explicit instantiations of the function
template DllExport int
spline3dEval<1, double>( const int *, const SurrealS<1> *, const double *, SurrealS<1> *);

template DllExport int
spline3dEval<1, SurrealS<1> >( const int *, const SurrealS<1> *, const SurrealS<1> *,
                               SurrealS<1> *);



/*
 * inverse evaluates the 3D BSpline
 *       where: ivec  - integer BSpline data
 *              data  - double BSpline data
 *              point - XYZ target (3 in len)
 *              uvw   - on input the seed parameters
 *                      returned UVW position
 *
 *              returns > 0 - OK (the number of iterations used)
 *                      -2 zero determinate, -3 hit 20 iterations
 *
 */
// template<class T, class T2>
// static int
// EG_spline3dInvEval_impl(int *ivec, T *data, T2 *uv, T *point)
int
spline3dInvEval( const int *ivec, const double *data, const double *point,
                 double *uvw)
{
  int    nKu, nKv, nKw, count, stat;
  double result[30], a00, a10, a20, a11, a21, a22, det, dist, dx[3];
  double b0, b1, b2, delta[3];
  const double *Kv, *Kw;

  /* set up the knot sequences */
  nKu = ivec[2];
  nKv = ivec[5];
  nKw = ivec[8];
  Kv  = data + nKu;
  Kw  = data + nKu + nKv;

  /* newton iteration */
  for (count = 0; count < 20; count++) {
    stat = spline3dEval(ivec, data, uvw, result);
    if (stat != 0) return stat;
    dx[0] = result[0] - point[0];
    dx[1] = result[1] - point[1];
    dx[2] = result[2] - point[2];
    dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (dist < 1.e-12) break;

    b0  = -dx[0]*result[ 3] - dx[1]*result[ 4] - dx[2]*result[ 5];
    b1  = -dx[0]*result[ 9] - dx[1]*result[10] - dx[2]*result[11];
    b2  = -dx[0]*result[18] - dx[1]*result[19] - dx[2]*result[20];

    a00 = result[ 3]*result[ 3] + result[ 4]*result[ 4] + result[ 5]*result[ 5] +
              dx[ 0]*result[ 6] +     dx[ 1]*result[ 7] +     dx[ 2]*result[ 8];
    a10 = result[ 3]*result[ 9] + result[ 4]*result[10] + result[ 5]*result[11] +
              dx[ 0]*result[12] +     dx[ 1]*result[13] +     dx[ 2]*result[14];
    a20 = result[ 3]*result[18] + result[ 4]*result[19] + result[ 5]*result[20] +
              dx[ 0]*result[21] +     dx[ 1]*result[22] +     dx[ 2]*result[23];
    a11 = result[ 9]*result[ 9] + result[10]*result[10] + result[11]*result[11] +
              dx[ 0]*result[15] +     dx[ 1]*result[16] +     dx[ 2]*result[17];
    a21 = result[ 9]*result[18] + result[10]*result[19] + result[11]*result[20] +
              dx[ 0]*result[24] +     dx[ 1]*result[25] +     dx[ 2]*result[26];
    a22 = result[18]*result[18] + result[19]*result[19] + result[20]*result[20] +
              dx[ 0]*result[27] +     dx[ 1]*result[28] +     dx[ 2]*result[29];

    det = a00*(a11*a22 - a21*a21) - a10*(a10*a22 - a21*a20) -
          a20*(a10*a21 - a11*a20);
    if (det == 0.0) return -2;
    det      = 1.0/det;
    delta[0] = det*( b0*(a11*a22 - a21*a21) -  b1*(a10*a22 - a21*a20) -
                     b2*(a10*a21 - a11*a20));
    delta[1] = det*(a00*( b1*a22 -  b2*a21) - a10*( b0*a22 -  b2*a20) -
                    a20*( b0*a21 -  b1*a20));
    delta[2] = det*(a00*(a11*b2  - a21*b1)  - a10*(a10*b2  - a21*b0)  -
                    a20*(a10*b1  - a11*b0));
 // printf(" deltas = %le %le %le\n", delta[0], delta[1], delta[2]);

    /* update and limit (if necessary) */
    uvw[0] += delta[0];
    if (uvw[0] < data[0])     uvw[0] = data[0];
    if (uvw[0] > data[nKu-1]) uvw[0] = data[nKu-1];
    uvw[1] += delta[1];
    if (uvw[1] < Kv[0])       uvw[1] = Kv[0];
    if (uvw[1] > Kv[nKv-1])   uvw[1] = Kv[nKv-1];
    uvw[2] += delta[2];
    if (uvw[2] < Kw[0])       uvw[2] = Kw[0];
    if (uvw[2] > Kw[nKw-1])   uvw[2] = Kw[nKw-1];
  }
  if (count == 20) return -3;

  return count;
}


#ifdef STANDALONE
int main(int argc, char *argv[])
{
  int    i, j, k, m, n, stat;
  double data[2048], result[30], uvw[3], val;
/*
  int    header[9] = { 1, 2, 4,  1, 2, 4,  1, 2, 4 };
  double knots[ 4] = { 0.0, 0.0, 1.0, 1.0 };
 */
/*
  int    header[9] = { 1, 5, 7,  1, 5, 7,  1, 5, 7 };
  double knots[ 7] = { 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0 };
 */
/*
  int    header[9] = { 2, 6, 9,  2, 6, 9,  2, 6, 9 };
  double knots[ 9] = { 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0 };
*/
/*
  int    header[9] = { 3, 4, 8,  3, 4, 8,  3, 4, 8 };
  double knots[ 8]  = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
 */
  int    header[9] = { 3, 6, 10,  3, 6, 10,  3, 6, 10 };
  double knots[10]  = { 0.0, 0.0, 0.0, 0.0, 0.3333, 0.66667, 1.0, 1.0, 1.0, 1.0 };

  if (argc != 4) {
    printf(" uasge: spline3d u v w\n");
    return 1;
  }

  if (header[0]+header[1]+1 != header[2]) {
    printf(" Problem with U spline header!\n");
    return 1;
  }
  if (header[3]+header[4]+1 != header[5]) {
    printf(" Problem with V spline header!\n");
    return 1;
  }
  if (header[6]+header[7]+1 != header[8]) {
    printf(" Problem with W spline header!\n");
    return 1;
  }
  if ((header[2] != header[5]) || (header[2] != header[8])) {
    printf(" All Knot sequences must be the same for this tester!\n");
    return 1;
  }

  /* 3 Knot sequences -- mst all be the same */
  for (n = 0; n < header[2]; n++)       data[n] = knots[n];
  for (j = 0; j < header[5]; j++, n++)  data[n] = knots[j];
  for (j = 0; j < header[8]; j++, n++)  data[n] = knots[j];

  /* Control Points -- just an equally spaced XYZ grid */
  for (m = k = 0; k < header[7]; k++)
    for (j = 0; j < header[4]; j++)
      for (i = 0; i < header[1]; i++, m++) {
        val = i;
        data[3*m  +n] = val/(header[1]-1);
        val = j;
        data[3*m+1+n] = val/(header[4]-1);
        val = k;
        data[3*m+2+n] = val/(header[7]-1);
      }

  /* set our parametric position */
  uvw[0] = atof(argv[1]);
  uvw[1] = atof(argv[2]);
  uvw[2] = atof(argv[3]);
  printf(" uvw = %lf %lf %lf\n\n", uvw[0], uvw[1], uvw[2]);

  stat = spline3dEval(header, data, uvw, result);
  printf(" xyz =  0: %10lf  1: %10lf  2: %10lf   %d\n", result[ 0], result[ 1], result[ 2], stat);
  printf(" du  =  3: %10lf  4: %10lf  5: %10lf\n",      result[ 3], result[ 4], result[ 5]);
  printf(" du2 =  6: %10lf  7: %10lf  8: %10lf\n",      result[ 6], result[ 7], result[ 8]);
  printf(" dv  =  9: %10lf 10: %10lf 11: %10lf\n",      result[ 9], result[10], result[11]);
  printf(" duv = 12: %10lf 13: %10lf 14: %10lf\n",      result[12], result[13], result[14]);
  printf(" dv2 = 15: %10lf 16: %10lf 17: %10lf\n",      result[15], result[16], result[17]);
  printf(" dw  = 18: %10lf 19: %10lf 20: %10lf\n",      result[18], result[19], result[20]);
  printf(" duw = 21: %10lf 22: %10lf 23: %10lf\n",      result[21], result[22], result[23]);
  printf(" dvw = 24: %10lf 25: %10lf 26: %10lf\n",      result[24], result[25], result[26]);
  printf(" dw2 = 27: %10lf 27: %10lf 29: %10lf\n\n",    result[27], result[28], result[29]);

  uvw[0] = uvw[1] = uvw[2] = 0.5;
  stat = spline3dInvEval(header, data, result, uvw);
  printf(" uvw = %lf %lf %lf   %d\n", uvw[0], uvw[1], uvw[2], stat);

  return 0;
}
#endif


// clang spline3d.c -DSTANDALONE
// ./a.out 0.02 0.02 0.02
