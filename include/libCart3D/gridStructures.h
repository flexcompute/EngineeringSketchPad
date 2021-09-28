/*
 *  $Id: gridStructures.h,v 1.1.1.1 2004/03/04 23:22:30 aftosmis Exp $
 *
 */

#ifndef __GRIDSTRUCTURES_H_
#define __GRIDSTRUCTURES_H_

#include "geomStructures.h"
#include "cartCells.h"

/* 
 * this is a more heavyweight structure for working with "complete"
 * meshes 
 */
struct basicGridStructure {
  p_tsFace    a_Faces;      /* flow face list for this grid -indexd into a_U */
  p_tsCutFace a_cFaces;     /*  cut face list for this grid -indexd into a_U */
  p_tsTinyHex a_Cells;      /*     all HEXES (volume + cut) in the subdomain */
  p_tsCutCell a_cCells;     /* cutCell (incl. splitcells) for this subdomain */
  p_tsState   a_U;          /*          state vector for all Control Volumes */
  tstPolys      tPolys;     /*    hook to triangle poly info struct & arrays */
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
