
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
 * $Id: cartCells.h,v 1.5 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __CARTCELLS_H_
#define __CARTCELLS_H_

#include "int64.h"
#include "basicTypes.h"

                                        /* -- define tiny hex type flags --  */
typedef enum {UNSET_HEX, FLOW_HEX, CUT_HEX, SPLIT_HEX, SOLID_HEX} tinyHexType;

#define UNSPLIT_INDEX -9999999          /* <- used in tsCutCell.splitIndex   */

typedef struct TinyHexStructure tsTinyHex;          /* ... low storage hexes */
typedef tsTinyHex * p_tsTinyHex;

struct TinyHexStructure {          /* bare-bones version of a cartesian cell */
  INT64       name;
  char        ref[DIM];
  byte        flagByte;
};
/* NOTES:
 * 1. This struct only uses 12 bytes, but occupies 16 bytes because of the
 *    INT64 width setting the struct stride. There's plenty of room.
 */

typedef  struct      CutCellStructure     tsCutCell;
typedef  tsCutCell*  p_tsCutCell;
struct CutCellStructure{    /*  ...Generic Cut Cell ------------   */
  int      *p_IntTriList;   /*     ptr to list of intersect tris   */
  double   *p_area;         /*     ptr to list of tPoly areas      */
  p_dpoint3 p_centroid;     /*     ptr to list of tPoly centroids  */
  dpoint3   normal;         /*     agglomerated weighted normV     */
  dpoint3   surfCentroid;   /*  agglomerated weighted surf centroid*/
  dpoint3   centroid;       /*     volume centroid FLOW polyhedron */
  double    volume;         /*     total volume of FLOW polyhedron */
  int       nIntTri;        /*     # of intersected triangles      */
  int       splitIndex;     /*     Split=indx of 1st kid, else flag*/
  int       bc_id;          /*     id of surf boundary type        */
  char      nMarked;        /* number of Marked and Touched cells  */
  char      nTouched;       /* during BC info restriction, cf. bc.c::BC_initialize() */
};
/* NOTES:
 * 1. This struct was recently reduced from 128 bytes to 120 bytes by
 *    moving the nIntTri field lower down, for optimal packing.
 **/

/* Linearized cut-cell data */
typedef  struct         LinCutCellStructure   tsLinCutCell;
typedef  tsLinCutCell * p_tsLinCutCell;
struct LinCutCellStructure{
  double   *p_LINarea;
  p_dpoint3 p_LINcentroid;
  dpoint3   LINnormal;
  dpoint3   LINcentroid;
  dpoint3   LINsurfCentroid;
};

typedef  struct     CutFaceStructure     tsCutFace;
typedef  tsCutFace * p_tsCutFace;
struct CutFaceStructure {         /* this is "fuller" struc for all          */
  int     adjCell[2];             /* faces attached to at least 1 cut cell   */
  dpoint3 centroid;               /* in X,Y,Z coordinates                    */
  double  area;
  char    dir;                    /* orientation of face X,Y,Z               */
};
/* NOTES:
 * 1. This struct uses 41 bytes, but the struct is padded to 48 bytes
 *
 *   int   int
 *   double
 *   double
 *   double
 *   double
 *   char [7 bytes leftover on this row]
 **/

/* Linearized cut-face data */
typedef  struct         LinCutFaceStructure   tsLinCutFace;
typedef  tsLinCutFace * p_tsLinCutFace;
struct LinCutFaceStructure {
  double  LINarea;
  dpoint3 LINcentroid;
};

typedef  struct     HexFaceStructure     tsFace;
typedef  tsFace * p_tsFace;

struct HexFaceStructure{
  int       adjCell[2];          /* faces attached to at least 1 cut cell    */
  short int faceLoc[2];          /* location of face on interface Cells      */
  byte      dir;                 /* orientation of face X,Y,Z                */
  byte      size0;               /* refinement level of transverse dirs      */
  byte      size1;               /* (these are cyclically ordered)           */
  mtype     mark:2;              /* mark for keeping track                   */
  unsigned  :0;                  /* pad to end of word                       */
};
/* NOTES:
 * 1. This struct is nearly fully packed into 16 bytes:
 *     int
 *     int
 *     short short
 *     byte byte byte 2-bit-mark [6 bits left on this row]
 **/

#endif /* __CARTCELLS_H_ */
