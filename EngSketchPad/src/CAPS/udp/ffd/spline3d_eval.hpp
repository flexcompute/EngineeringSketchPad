
#ifndef EGADS_SLINE3D_EVAL_H_
#define EGADS_SLINE3D_EVAL_H_

/*
 * Copyright (C) 2011/2025
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

#include "Surreal/SurrealS.h"

#define SPLINE3D_MAXDEG   5

/*
 * evaluates the 3D BSpline, 1st and 2nd derivatives
 *       where: ivec  - integer BSpline data
 *              data  - DOUBLE BSpline data
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
extern "C"
int
spline3dEval( const int *ivec, const double *data, const double *uvw, double *deriv);

template<int N, class T>
int
spline3dEval( const int *ivec,const SurrealS<N> *data, const T *uvw,
              SurrealS<N> *deriv);


/*
 * inverse evaluates the 3D BSpline
 *       where: ivec  - integer BSpline data
 *              data  - DOUBLE BSpline data
 *              point - XYZ target (3 in len)
 *              uvw   - on input the seed parameters
 *                      returned UVW position
 *
 *              returns > 0 - OK (the number of iterations used)
 *                      -2 zero determinate, -3 hit 20 iterations
 *
 */
int
spline3dInvEval(const int *ivec, const double *data, const double *point,
                double *uvw);

#endif // EGADS_SLINE3D_EVAL_H_
