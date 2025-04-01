
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
 *  $Id: gridStructures.h,v 1.3 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __GRIDSTRUCTURES_H_
#define __GRIDSTRUCTURES_H_

#include "geomStructures.h"
#include "cartCells.h"

/*
 * this is a more heavyweight structure for working with "complete"
 * meshes
 */
struct basicGridStructure {
  p_tsFace    a_Faces;     /*  flow face list for this grid -indexd into a_U */
  p_tsCutFace a_cFaces;    /*   cut face list for this grid -indexd into a_U */
  p_tsTinyHex a_Cells;     /*      all HEXES (volume + cut) in the subdomain */
  p_tsCutCell a_cCells;    /*  cutCell (incl. splitcells) for this subdomain */
  p_tsState   a_U;         /*           state vector for all Control Volumes */
  p_tsState   a_Uo;        /* any other state vector for all Control Volumes */
  tstPolys      tPolys;    /*     hook to triangle poly info struct & arrays */
  int nFaces;        /*                          total # of faces in a_Faces */
  int nFacesXYZ[DIM];/*                  No. of Cart faces in each direction */
  int nCutFaces;     /*                         total # of faces in a_cFaces */
  int nVolHexes;     /*  No. of Cart. Flow field cells not touching geometry */
  int nCutHexes;     /*   the No. of CARTESIAN HEXAHEDRA cut by the boundary */
  int nSplitCells;   /*              # of CntlVols from multiRegion cuthexes */
  int nCells;        /*  total no of control vol in this grid (in partition) */
};

typedef struct basicGridStructure tsBasicGrid;
typedef tsBasicGrid *p_tsBasicGrid;

#endif /* __GRIDSTRUCTURES_H_ */
