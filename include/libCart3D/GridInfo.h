/*
 * $Id: GridInfo.h,v 1.4 2004/03/24 18:38:08 smurman Exp $
 */

#ifndef __GRIDINFO_H_
#define __GRIDINFO_H_

#include "c3d_global.h"
#include "basicTypes.h"
#include "geomStructures.h"
#include "IOinfo.h"

#ifndef MAXNBITS
# define MAXNBITS 21
#endif

typedef struct GridInfoStructure tsGinfo;
typedef tsGinfo *p_tsGinfo;

struct GridInfoStructure {
  float      minBound[DIM], maxBound[DIM];/*      Domain extent (Xmin, Xmax) */
  double     fine_spacing[DIM];/*     mesh spacing on finest allowable level */
  int        M[DIM];       /*  max integer dimension in each coord direction */
  int        nBits;           /*     # bits of resolution for each direction */
  int        maxAllowRef[DIM];/*  Max # of refinements supported in each dir */
  int        maxCurrentRef[DIM];
  dpoint3    h[MAXNBITS];           /* mesh size corresponding to each level */
  int        nDiv[DIM];               /* initial background mesh size        */
  int        coarseIntDelta[DIM];     /* coarsest cell size, integer coords  */
  int        finestCellLevel;         /* finest refinement level in the mesh */
  bool       meshInternal;
  bool       isotropic;
  p_tsIOinfo p_fileInfo;
  p_tsTriangulation p_surf;         /*   triangulation this grid is built on */
};

#endif /* __GRIDINFO_H_ */
