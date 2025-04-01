
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: GridInfo.h,v 1.5 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

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
