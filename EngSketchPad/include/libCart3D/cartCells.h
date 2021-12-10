/*
 * $Id: cartCells.h,v 1.3 2011/08/29 17:22:16 mnemec Exp $
 *
 */
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

typedef  struct      CutCellStructure     tsCutCell;
typedef  tsCutCell*  p_tsCutCell;
struct CutCellStructure{    /*  ...Generic Cut Cell ------------   */
  int      *p_IntTriList;   /*     ptr to list of intersect tris   */
  double   *p_area;         /*     ptr to list of tPoly areas      */
  p_dpoint3 p_centroid;     /*     ptr to list of tPoly centroids  */
  int       nIntTri;        /*     # of intersected triangles      */
  dpoint3   normal;         /*     agglomerated weighted normV     */
  dpoint3   surfCentroid;   /*  agglomerated weighted surf centroid*/
  dpoint3   centroid;       /*     volume centroid FLOW polyhedron */
  double    volume;         /*     total volume of FLOW polyhedron */
  int       splitIndex;     /*     Split=indx of 1st kid, else flag*/
  int       bc_id;          /*     id of surf boundary type        */
  char      nMarked;        /* number of Marked and Touched cells  */
  char      nTouched;       /* during BC info restriction, cf. bc.c::BC_initialize() */
};

/* Linearized cut-cell data, MN, April 2005, lib Aug 2011 */
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

/* Linearized cut-face data, MN, April 2005, lib Aug 2011 */
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
/*
 ---------
*/

#endif /* __CARTCELLS_H_ */
