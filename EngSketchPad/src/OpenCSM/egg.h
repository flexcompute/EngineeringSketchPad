/*
 ************************************************************************
 *                                                                      *
 * egg.h -- header for external grid generator (egg)                    *
 *                                                                      *
 *           Written by John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2018/2022  John F. Dannenhoffer, III (Syracuse University)
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

#ifndef _EGG_H_
#define _EGG_H_

#include <stdio.h>


/*
 ************************************************************************
 *                                                                      *
 * Callable routines                                                    *
 *                                                                      *
 ************************************************************************
 */

/* generate a Grid */
int eggGenerate(double uv[],            /* (in)  array of coordinates for Boundary Points */
                                        /* (out) array of coordinates for Grid     Points */
                int    nbnd[],          /* (in)  number of Boundary Points in each loop */
                void   **gridP);        /* (out) pointer to grid (which is a grid_T*) */

/* morph a Grid */
int eggMorph(void   *gridP,             /* (in)  pointer to grid */
             double uvnew[],            /* (in)  array of coordinates for Boundary Points */
                                        /* (out) array of coordinates for Grid     Points */
             void   **newGridP);        /* (out) pointer to new grid (which is a grid_T*) */

/* get info about a Grid */
int eggInfo(void   *gridP,              /* (in)  pointer to Grid */
            int    *npnt,               /* (out) number of Points */
            int    *nbnd,               /* (out) number of boundary Points */
      const double *uv[],               /* (out) parametric coordinates */
      const int    *p[],                /* (out) parent points */
            int    *ntri,               /* (out) number of Triangles */
      const int    *tris[]);            /* (out) indices of Points in Triangles (bias-0) */

/* dump a Grid to a file */
int eggDump(void   *gridP,              /* (in)  pointer to grid */
            FILE   *fp);                /* (in)  pointer to ASCII file */

/* load a Grid from a file */
int eggLoad(FILE   *fp,                 /* (in)  pointer to ASCII file */
            void   **gridP);            /* (out) pointer to *grid (which is a grid_T*) */

/* free a Grid structure */
int eggFree(void   *gridP);             /* (in)  pointer to grid */

/*
 ************************************************************************
 *                                                                      *
 * Return codes                                                         *
 *                                                                      *
 ************************************************************************
 */

#define  SUCCESS                    0

#define  MALLOC_ERROR              -901
#define  BAD_POINT_INDEX           -902
#define  BAD_TRIANGLE_INDEX        -903
#define  CANNOT_SWAP               -904
#define  COULD_NOT_RECOVER_BND     -905
#define  NUMBER_OF_POINT_MISMATCH  -906

#endif   /* _EGG_H_ */
