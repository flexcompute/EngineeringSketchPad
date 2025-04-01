//
// Written by Dr. Ryan Durscher AFRL/RQVC
// Place clearance statement here.
//
// AFLR3 interface functions - Modified from functions provided with
// AFLR3_LIB source (aflr3.c) written by David L. Marcum


#ifdef WIN32
#define strcasecmp stricmp
#define strncasecmp _strnicmp
#define strtok_r   strtok_s
#endif

#include <math.h>

#include "aimUtil.h"
#include "meshUtils.h"     // Collection of helper functions for meshing
#include "miscUtils.h"
#include "egads.h"

#define _ENABLE_BL_

#include <aflr3/AFLR3_LIB.h> // AFLR3_API Library include
#include <aflr4/AFLR4_LIB.h> // Bring in AFLR4 API library
#include <aflr43/AFLR43_LIB.h>
#include <egads_aflr4/EGADS_AFLR4_LIB_INC.h>

// UG_IO Library (file I/O) include
// UG_CPP Library (c++ code) include

#include <ug_io/UG_IO_LIB_INC.h>
#include <ug_cpp/UG_CPP_LIB_INC.h>

// UG_GQ Library (grid quality) include
// This is optional and not required for implementation.

#include <ug_gq/UG_GQ_LIB_INC.h>

// Includes for version functions
// This is optional and not required for implementation.

#include <otb/OTB_LIB_INC.h>
#include <rec3/REC3_LIB_INC.h>
#include <ug3/UG3_LIB_INC.h>
#include <dftr3/DFTR3_LIB_INC.h>
#include <ice3/ICE3_LIB_INC.h>

#ifdef _ENABLE_BL_
#include <aflr2c/AFLR2_LIB.h>
#include <anbl3/ANBL3_LIB.h>

/*
#include "bl1/BL1_LIB_INC.h"
#include "dgeom/DGEOM_LIB_INC.h"
#include "egen/EGEN_LIB_INC.h"
#include "qtb/QTB_LIB_INC.h"
#include "rec2/REC2_LIB_INC.h"
#include "ug2/UG2_LIB_INC.h"
#include "dftr2/DFTR2_LIB_INC.h"
#include "ice2/ICE2_LIB_INC.h"
 */
#endif

#include "aflr3_Interface.h"

#ifndef S_SPLINT_S
#define AFLR_STATUS(aimInfo, statys, ...) \
if (status != 0) { status = CAPS_EXECERR; AIM_STATUS(aimInfo, status, ##__VA_ARGS__); }
#else
extern void AFLR_STATUS(void *aimInfo, int status, ...);
#endif

int aflr3_to_MeshStruct( const AFLR_Grid *grid,
                         meshStruct *genUnstrMesh)
{

    int status; // Function return status

    int i, j, elementIndex; // Indexing variable

    int numPoint;
    int defaultVolID = 1; // Defailt volume ID

    meshAnalysisTypeEnum analysisType;

    analysisType = genUnstrMesh->analysisType;

    // Cleanup existing node and elements
    (void) destroy_meshNodes(genUnstrMesh);

    (void) destroy_meshElements(genUnstrMesh);

    (void) destroy_meshQuickRefStruct(&genUnstrMesh->meshQuickRef);
    genUnstrMesh->meshType = VolumeMesh;

    //printf ("Transferring mesh to general unstructured structure\n");

    // Numbers
    genUnstrMesh->numNode = grid->Number_of_Nodes;
    genUnstrMesh->numElement = grid->Number_of_Surf_Trias  +
                               grid->Number_of_Surf_Quads  +
                               grid->Number_of_Vol_Tets    +
                               grid->Number_of_Vol_Pents_5 +
                               grid->Number_of_Vol_Pents_6 +
                               grid->Number_of_Vol_Hexs;

    genUnstrMesh->meshQuickRef.useStartIndex = (int) true;

    genUnstrMesh->meshQuickRef.numTriangle      = grid->Number_of_Surf_Trias;
    genUnstrMesh->meshQuickRef.numQuadrilateral = grid->Number_of_Surf_Quads;

    genUnstrMesh->meshQuickRef.numTetrahedral = grid->Number_of_Vol_Tets;
    genUnstrMesh->meshQuickRef.numPyramid     = grid->Number_of_Vol_Pents_5;
    genUnstrMesh->meshQuickRef.numPrism       = grid->Number_of_Vol_Pents_6;
    genUnstrMesh->meshQuickRef.numHexahedral  = grid->Number_of_Vol_Hexs;

    // Allocation

    // Nodes - allocate
    genUnstrMesh->node = (meshNodeStruct *)
                         EG_alloc(genUnstrMesh->numNode*sizeof(meshNodeStruct));
    if (genUnstrMesh->node == NULL) {
#if !defined(_MSC_VER) || (_MSC_VER >= 1800)
/*@-formatcode@*/
      printf("Failed to allocate %d meshNodeStruct (%zu bytes)\n",
             genUnstrMesh->numNode, genUnstrMesh->numNode*sizeof(meshNodeStruct));
/*@+formatcode@*/
#endif
      return EGADS_MALLOC;
    }

    // Elements - allocate
    genUnstrMesh->element = (meshElementStruct *)
                   EG_alloc(genUnstrMesh->numElement*sizeof(meshElementStruct));
    if (genUnstrMesh->element == NULL) {
#if !defined(_MSC_VER) || (_MSC_VER >= 1800)
/*@-formatcode@*/
        printf("Failed to allocate %d meshElementStruct (%zu bytes)\n",
               genUnstrMesh->numElement,
               genUnstrMesh->numElement*sizeof(meshElementStruct));
/*@+formatcode@*/
#endif
        EG_free(genUnstrMesh->node);
        genUnstrMesh->node = NULL;
        return EGADS_MALLOC;
    }

    // Initialize
    for (i = 0; i < genUnstrMesh->numNode; i++) {
        status = initiate_meshNodeStruct(&genUnstrMesh->node[i], analysisType);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    for (i = 0; i < genUnstrMesh->numElement; i++ ) {
        status = initiate_meshElementStruct(&genUnstrMesh->element[i],
                                            analysisType);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Nodes - set
    for (i = 0; i < genUnstrMesh->numNode; i++) {

        // Copy node data
        genUnstrMesh->node[i].nodeID = i+1;

        genUnstrMesh->node[i].xyz[0] = grid->Coordinates[i+1][0];
        genUnstrMesh->node[i].xyz[1] = grid->Coordinates[i+1][1];
        genUnstrMesh->node[i].xyz[2] = grid->Coordinates[i+1][2];
    }


    // Start of element index
    elementIndex = 0;

    // Elements-Set triangles
    if (grid->Number_of_Surf_Trias > 0)
        genUnstrMesh->meshQuickRef.startIndexTriangle = elementIndex;

    numPoint = 0;
    for (i = 0; i < grid->Number_of_Surf_Trias; i++) {

        genUnstrMesh->element[elementIndex].elementType = Triangle;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;

        genUnstrMesh->element[elementIndex].markerID = grid->Surf_ID_Flag[i+1];

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                grid->Surf_Tria_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set quadrilateral
    if (grid->Number_of_Surf_Quads > 0)
        genUnstrMesh->meshQuickRef.startIndexQuadrilateral = elementIndex;

    for (i = 0; i < grid->Number_of_Surf_Quads; i++) {

        genUnstrMesh->element[elementIndex].elementType = Quadrilateral;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID    =
            grid->Surf_ID_Flag[grid->Number_of_Surf_Trias+i+1];

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                grid->Surf_Quad_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Tetrahedral
    if (grid->Number_of_Vol_Tets > 0)
        genUnstrMesh->meshQuickRef.startIndexTetrahedral = elementIndex;

    for (i = 0; i < grid->Number_of_Vol_Tets; i++) {

        genUnstrMesh->element[elementIndex].elementType = Tetrahedral;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID    = defaultVolID;

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                grid->Vol_Tet_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Pyramid
    if (grid->Number_of_Vol_Pents_5 > 0)
        genUnstrMesh->meshQuickRef.startIndexPyramid = elementIndex;

    for (i = 0; i < grid->Number_of_Vol_Pents_5; i++) {

        genUnstrMesh->element[elementIndex].elementType = Pyramid;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID    = defaultVolID;

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                grid->Vol_Pent_5_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Prism
    if (grid->Number_of_Vol_Pents_6 > 0)
        genUnstrMesh->meshQuickRef.startIndexPrism = elementIndex;

    for (i = 0; i < grid->Number_of_Vol_Pents_6; i++) {

        genUnstrMesh->element[elementIndex].elementType = Prism;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID    = defaultVolID;

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                grid->Vol_Pent_6_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Hexa
    if (grid->Number_of_Vol_Hexs > 0)
        genUnstrMesh->meshQuickRef.startIndexHexahedral = elementIndex;

    for (i = 0; i < grid->Number_of_Vol_Hexs; i++) {

        genUnstrMesh->element[elementIndex].elementType = Hexahedral;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID    = defaultVolID;

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                grid->Vol_Hex_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    status = CAPS_SUCCESS;


cleanup:
    if (status != CAPS_SUCCESS)
        printf("Premature exit in aflr3_to_MeshStruct status = %d\n", status);

    return status;
}


void initialize_AFLR_Grid(AFLR_Grid* grid)
{
  grid->Edge_ID_Flag = NULL;
  grid->Surf_Grid_BC_Flag = NULL;
  grid->Surf_ID_Flag = NULL;
  grid->Surf_Reconnection_Flag = NULL;
  grid->Surf_Edge_Connectivity = NULL;
  grid->Surf_Tria_Connectivity = NULL;
  grid->Surf_Quad_Connectivity = NULL;
  grid->Vol_ID_Flag = NULL;
  grid->Vol_Tet_Connectivity = NULL;
  grid->Vol_Pent_5_Connectivity = NULL;
  grid->Vol_Pent_6_Connectivity = NULL;
  grid->Vol_Hex_Connectivity = NULL;

  grid->Coordinates = NULL;

  grid->BL_Normal_Spacing = NULL;
  grid->BL_Thickness = NULL;

  grid->Number_of_BL_Vol_Tets = 0;
  grid->Number_of_Nodes = 0;
  grid->Number_of_Surf_Edges = 0;
  grid->Number_of_Surf_Quads = 0;
  grid->Number_of_Surf_Trias = 0;
  grid->Number_of_Vol_Hexs = 0;
  grid->Number_of_Vol_Pents_5 = 0;
  grid->Number_of_Vol_Pents_6 = 0;
  grid->Number_of_Vol_Tets = 0;

  initiate_mapAttrToIndexStruct(&grid->groupMap);
}


void destroy_AFLR_Grid(AFLR_Grid* grid)
{
  // Free grid generation and parameter array space in structures.
  // Note that aflr3_grid_generator frees all background data.
  ug_free( grid->Edge_ID_Flag );
  ug_free( grid->Surf_Grid_BC_Flag );
  ug_free( grid->Surf_ID_Flag );
  ug_free( grid->Surf_Reconnection_Flag );
  ug_free( grid->Surf_Edge_Connectivity );
  ug_free( grid->Surf_Quad_Connectivity );
  ug_free( grid->Surf_Tria_Connectivity );
  ug_free( grid->Vol_Hex_Connectivity );
  ug_free( grid->Vol_ID_Flag );
  ug_free( grid->Vol_Pent_5_Connectivity );
  ug_free( grid->Vol_Pent_6_Connectivity );
  ug_free( grid->Vol_Tet_Connectivity );
  ug_free( grid->Coordinates );
  ug_free( grid->BL_Normal_Spacing );
  ug_free( grid->BL_Thickness );

  destroy_mapAttrToIndexStruct(&grid->groupMap);

  initialize_AFLR_Grid(grid);
}

int append_AFLR_Grid(void *aimInfo,
                     AFLR_Grid* domain,
                     int zone,
                     AFLR_Grid* grid)
{
#if 0
  int i, j, ig;
  INT_ ierr = 0;

  INT_ Number_of_BL_Vol_Tets = domain->Number_of_BL_Vol_Tets + grid->Number_of_BL_Vol_Tets;
  INT_ Number_of_Nodes       = domain->Number_of_Nodes       + grid->Number_of_Nodes      ;
  INT_ Number_of_Surf_Edges  = domain->Number_of_Surf_Edges  + grid->Number_of_Surf_Edges ;
  INT_ Number_of_Surf_Trias  = domain->Number_of_Surf_Trias  + grid->Number_of_Surf_Trias ;
  INT_ Number_of_Surf_Quads  = domain->Number_of_Surf_Quads  + grid->Number_of_Surf_Quads ;
  INT_ Number_of_Vol_Tets    = domain->Number_of_Vol_Tets    + grid->Number_of_Vol_Tets   ;
  INT_ Number_of_Vol_Pents_5 = domain->Number_of_Vol_Pents_5 + grid->Number_of_Vol_Pents_5;
  INT_ Number_of_Vol_Pents_6 = domain->Number_of_Vol_Pents_6 + grid->Number_of_Vol_Pents_6;
  INT_ Number_of_Vol_Hexs    = domain->Number_of_Vol_Hexs    + grid->Number_of_Vol_Hexs   ;

  INT_ Number_of_Surf = Number_of_Surf_Trias + Number_of_Surf_Quads;
  INT_ Number_of_Vol = Number_of_Vol_Hexs + Number_of_Vol_Pents_5 + Number_of_Vol_Pents_6 + Number_of_Vol_Tets;

  if (Number_of_Surf_Quads > 0) {
      AIM_ERROR(aimInfo, "Developer error! append_AFLR_Grid currently does not work correctly for Quad elements.");
      return CAPS_NOTIMPLEMENT;
  }

  grid->Edge_ID_Flag            = (INT_1D *) ug_realloc (&ierr, grid->Edge_ID_Flag           , (Number_of_Surf_Edges +1) * sizeof (INT_1D));  AFLR_STATUS(aimInfo, ierr);
  grid->Surf_Grid_BC_Flag       = (INT_1D *) ug_realloc (&ierr, grid->Surf_Grid_BC_Flag      , (Number_of_Surf       +1) * sizeof (INT_1D));  AFLR_STATUS(aimInfo, ierr);
  grid->Surf_ID_Flag            = (INT_1D *) ug_realloc (&ierr, grid->Surf_ID_Flag           , (Number_of_Surf       +1) * sizeof (INT_1D));  AFLR_STATUS(aimInfo, ierr);
  grid->Surf_Reconnection_Flag  = (INT_1D *) ug_realloc (&ierr, grid->Surf_Reconnection_Flag , (Number_of_Surf       +1) * sizeof (INT_1D));  AFLR_STATUS(aimInfo, ierr);
  grid->Surf_Edge_Connectivity  = (INT_2D *) ug_realloc (&ierr, grid->Surf_Edge_Connectivity , (Number_of_Surf_Edges +1) * sizeof (INT_2D));  AFLR_STATUS(aimInfo, ierr);
  grid->Surf_Tria_Connectivity  = (INT_3D *) ug_realloc (&ierr, grid->Surf_Tria_Connectivity , (Number_of_Surf_Trias +1) * sizeof (INT_3D));  AFLR_STATUS(aimInfo, ierr);
  grid->Surf_Quad_Connectivity  = (INT_4D *) ug_realloc (&ierr, grid->Surf_Quad_Connectivity , (Number_of_Surf_Quads +1) * sizeof (INT_4D));  AFLR_STATUS(aimInfo, ierr);

  if (grid->Vol_ID_Flag != NULL) {
    grid->Vol_ID_Flag           = (INT_1D *) ug_realloc (&ierr, grid->Vol_ID_Flag            , (Number_of_Vol        +1) * sizeof (INT_1D));  AFLR_STATUS(aimInfo, ierr);
  }
  grid->Vol_Tet_Connectivity    = (INT_4D *) ug_realloc (&ierr, grid->Vol_Tet_Connectivity   , (Number_of_Vol_Tets   +1) * sizeof (INT_4D));  AFLR_STATUS(aimInfo, ierr);
  grid->Vol_Pent_5_Connectivity = (INT_5D *) ug_realloc (&ierr, grid->Vol_Pent_5_Connectivity, (Number_of_Vol_Pents_5+1) * sizeof (INT_5D));  AFLR_STATUS(aimInfo, ierr);
  grid->Vol_Pent_6_Connectivity = (INT_6D *) ug_realloc (&ierr, grid->Vol_Pent_6_Connectivity, (Number_of_Vol_Pents_6+1) * sizeof (INT_6D));  AFLR_STATUS(aimInfo, ierr);
  grid->Vol_Hex_Connectivity    = (INT_8D *) ug_realloc (&ierr, grid->Vol_Hex_Connectivity   , (Number_of_Vol_Hexs   +1) * sizeof (INT_8D));  AFLR_STATUS(aimInfo, ierr);

  grid->Coordinates       = (DOUBLE_3D *) ug_realloc (&ierr, grid->Coordinates      , (Number_of_Nodes+1) * sizeof (DOUBLE_3D));  AFLR_STATUS(aimInfo, ierr);

  if (grid->BL_Normal_Spacing != NULL) {
      grid->BL_Normal_Spacing = (DOUBLE_1D *) ug_realloc (&ierr, grid->BL_Normal_Spacing, (Number_of_Nodes+1) * sizeof (DOUBLE_1D));  AFLR_STATUS(aimInfo, ierr);
  }
  if (grid->BL_Thickness != NULL) {
      grid->BL_Thickness      = (DOUBLE_1D *) ug_realloc (&ierr, grid->BL_Thickness     , (Number_of_Nodes+1) * sizeof (DOUBLE_1D));  AFLR_STATUS(aimInfo, ierr);
  }

  ig = grid->Number_of_Nodes;
  for (i = 0; i < domain->Number_of_Nodes; i++) {
      grid->Coordinates[ig+i+1][0] = domain->Coordinates[i+1][0];
      grid->Coordinates[ig+i+1][1] = domain->Coordinates[i+1][1];
      grid->Coordinates[ig+i+1][2] = domain->Coordinates[i+1][2];
  }

  ig = grid->Number_of_Surf_Edges;
  for (i = 0; i < domain->Number_of_Surf_Edges; i++) {
      grid->Edge_ID_Flag[ig+i+1] = domain->Edge_ID_Flag[i+1];
  }

  ig = grid->Number_of_Surf_Trias + grid->Number_of_Surf_Quads;
  for (i = 0; i < domain->Number_of_Surf_Trias + domain->Number_of_Surf_Quads; i++) {
      grid->Surf_Grid_BC_Flag[ig+i+1] = domain->Surf_Grid_BC_Flag[i+1];
  }

  ig = grid->Number_of_Surf_Trias + grid->Number_of_Surf_Quads;
  for (i = 0; i < domain->Number_of_Surf_Trias + domain->Number_of_Surf_Quads; i++) {
      grid->Surf_ID_Flag[ig+i+1] = domain->Surf_ID_Flag[i+1];
  }

  ig = grid->Number_of_Surf_Trias + grid->Number_of_Surf_Quads;
  for (i = 0; i < domain->Number_of_Surf_Trias + domain->Number_of_Surf_Quads; i++) {
      grid->Surf_Reconnection_Flag[ig+i+1] = domain->Surf_Reconnection_Flag[i+1];
  }

  ig = grid->Number_of_Surf_Edges;
  for (i = 0; i < domain->Number_of_Surf_Edges; i++) {
      for (j = 0; j < 2; j++)
          grid->Surf_Edge_Connectivity[ig+i+1][j] = domain->Surf_Edge_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }

  ig = grid->Number_of_Surf_Trias;
  for (i = 0; i < domain->Number_of_Surf_Trias; i++) {
      for (j = 0; j < 3; j++)
        grid->Surf_Tria_Connectivity[ig+i+1][j] = domain->Surf_Tria_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }

  ig = grid->Number_of_Surf_Quads;
  for (i = 0; i < domain->Number_of_Surf_Quads; i++) {
      for (j = 0; j < 4; j++)
          grid->Surf_Quad_Connectivity[ig+i+1][j] = domain->Surf_Quad_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }

  if (grid->Vol_ID_Flag != NULL) {
      ig = grid->Number_of_Vol_Hexs + grid->Number_of_Vol_Pents_5 + grid->Number_of_Vol_Pents_6 + grid->Number_of_Vol_Tets;
      for (i = 0; i < Number_of_Vol - ig; i++) {
          grid->Vol_ID_Flag[ig+i+1] = domain->Vol_ID_Flag[i+1];
      }
  }

  ig = grid->Number_of_Vol_Tets;
  for (i = 0; i < domain->Number_of_Vol_Tets; i++) {
      for (j = 0; j < 4; j++)
          grid->Vol_Tet_Connectivity[ig+i+1][j] = domain->Vol_Tet_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }

  ig = grid->Number_of_Vol_Pents_5;
  for (i = 0; i < domain->Number_of_Vol_Pents_5; i++) {
      for (j = 0; j < 5; j++)
          grid->Vol_Pent_5_Connectivity[ig+i+1][j] = domain->Vol_Pent_5_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }

  ig = grid->Number_of_Vol_Pents_6;
  for (i = 0; i < domain->Number_of_Vol_Pents_6; i++) {
      for (j = 0; j < 6; j++)
          grid->Vol_Pent_6_Connectivity[ig+i+1][j] = domain->Vol_Pent_6_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }

  ig = grid->Number_of_Vol_Hexs;
  for (i = 0; i < domain->Number_of_Vol_Hexs; i++) {
      for (j = 0; j < 8; j++)
          grid->Vol_Hex_Connectivity[ig+i+1][j] = domain->Vol_Hex_Connectivity[i+1][j] + grid->Number_of_Nodes;
  }


  if (grid->BL_Normal_Spacing != NULL) {
    ig = grid->Number_of_Nodes;
    for (i = 0; i < domain->Number_of_Nodes; i++) {
      grid->BL_Normal_Spacing[ig+i+1]      = domain->BL_Normal_Spacing[i+1];
    }
  }

  if (grid->BL_Thickness != NULL) {
    ig = grid->Number_of_Nodes;
    for (i = 0; i < domain->Number_of_Nodes; i++) {
      grid->BL_Thickness[ig+i+1]      = domain->BL_Thickness[i+1];
    }
  }

  grid->Number_of_BL_Vol_Tets = Number_of_BL_Vol_Tets;
  grid->Number_of_Nodes       = Number_of_Nodes      ;
  grid->Number_of_Surf_Edges  = Number_of_Surf_Edges ;
  grid->Number_of_Surf_Trias  = Number_of_Surf_Trias ;
  grid->Number_of_Surf_Quads  = Number_of_Surf_Quads ;
  grid->Number_of_Vol_Tets    = Number_of_Vol_Tets   ;
  grid->Number_of_Vol_Pents_5 = Number_of_Vol_Pents_5;
  grid->Number_of_Vol_Pents_6 = Number_of_Vol_Pents_6;
  grid->Number_of_Vol_Hexs    = Number_of_Vol_Hexs   ;
#endif

  int status = CAPS_SUCCESS;
  INT_ final_mesh = 0;

  status = ug3_merge_mesh (final_mesh,
                           &domain->Number_of_Surf_Edges,
                           &domain->Number_of_Surf_Trias,
                           &domain->Number_of_Surf_Quads,
                           &domain->Number_of_BL_Vol_Tets,
                           &domain->Number_of_Vol_Tets,
                           &domain->Number_of_Vol_Pents_5,
                           &domain->Number_of_Vol_Pents_6,
                           &domain->Number_of_Vol_Hexs,
                           &domain->Number_of_Nodes,
                           zone+1,
                           &domain->Edge_ID_Flag,
                           &domain->Surf_Edge_Connectivity,
                           &domain->Surf_Grid_BC_Flag,
                           &domain->Surf_ID_Flag,
                           &domain->Surf_Reconnection_Flag,
                           &domain->Surf_Tria_Connectivity,
                           &domain->Surf_Quad_Connectivity,
                           &domain->Vol_ID_Flag,
                           &domain->Vol_Tet_Connectivity,
                           &domain->Vol_Pent_5_Connectivity,
                           &domain->Vol_Pent_6_Connectivity,
                           &domain->Vol_Hex_Connectivity,
                           &domain->Coordinates,
                           &domain->BL_Normal_Spacing,
                           &domain->BL_Thickness,

                           &grid->Number_of_Surf_Edges,
                           &grid->Number_of_Surf_Trias,
                           &grid->Number_of_Surf_Quads,
                           &grid->Number_of_BL_Vol_Tets,
                           &grid->Number_of_Vol_Tets,
                           &grid->Number_of_Vol_Pents_5,
                           &grid->Number_of_Vol_Pents_6,
                           &grid->Number_of_Vol_Hexs,
                           &grid->Number_of_Nodes,
                           &grid->Edge_ID_Flag,
                           &grid->Surf_Edge_Connectivity,
                           &grid->Surf_Grid_BC_Flag,
                           &grid->Surf_ID_Flag,
                           &grid->Surf_Reconnection_Flag,
                           &grid->Surf_Tria_Connectivity,
                           &grid->Surf_Quad_Connectivity,
                           &grid->Vol_ID_Flag,
                           &grid->Vol_Tet_Connectivity,
                           &grid->Vol_Pent_5_Connectivity,
                           &grid->Vol_Pent_6_Connectivity,
                           &grid->Vol_Hex_Connectivity,
                           &grid->Coordinates,
                           &grid->BL_Normal_Spacing,
                           &grid->BL_Thickness);
  AFLR_STATUS(aimInfo, status);

  status = merge_mapAttrToIndexStruct(&domain->groupMap, &grid->groupMap, &grid->groupMap);
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


int write_AFLR_Grid(void *aimInfo,
                    const char *fileName,
                    AFLR_Grid *grid)
{
    int status = CAPS_SUCCESS;
    int i;
    FILE *fp = NULL;
    char aimFile[PATH_MAX];

    INT_ Message_Flag = 0;

    // Write the mesh to disk

    snprintf(aimFile, PATH_MAX, "%s.lb8.ugrid", fileName);

    status = ug_io_write_grid_file(aimFile,
                                   Message_Flag,
                                   grid->Number_of_BL_Vol_Tets,
                                   grid->Number_of_Nodes,
                                   grid->Number_of_Surf_Quads,
                                   grid->Number_of_Surf_Trias,
                                   grid->Number_of_Vol_Hexs,
                                   grid->Number_of_Vol_Pents_5,
                                   grid->Number_of_Vol_Pents_6,
                                   grid->Number_of_Vol_Tets,
                                   grid->Surf_Grid_BC_Flag,
                                   grid->Surf_ID_Flag,
                                   grid->Surf_Reconnection_Flag,
                                   grid->Surf_Quad_Connectivity,
                                   grid->Surf_Tria_Connectivity,
                                   grid->Vol_Hex_Connectivity,
                                   grid->Vol_ID_Flag,
                                   grid->Vol_Pent_5_Connectivity,
                                   grid->Vol_Pent_6_Connectivity,
                                   grid->Vol_Tet_Connectivity,
                                   grid->Coordinates,
                                   grid->BL_Normal_Spacing,
                                   grid->BL_Thickness);
    AFLR_STATUS(aimInfo, status);

    snprintf(aimFile, PATH_MAX, "%s.mapbc", fileName);
    fp = fopen(aimFile, "w");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open file: %s", aimFile);
      status = CAPS_IOERR;
      goto cleanup;
    }

    fprintf(fp, "%d\n", grid->groupMap.numAttribute);
    for (i = 0; i < grid->groupMap.numAttribute; i++) {
      fprintf(fp, "%d 0 %s\n", grid->groupMap.attributeIndex[i], grid->groupMap.attributeName[i]);
    }

    status = CAPS_SUCCESS;
cleanup:
  /*@-dependenttrans@*/
    if (fp != NULL) fclose(fp);
  /*@+dependenttrans@*/

  return status;
}


int aflr3_Volume_Mesh (void *aimInfo,
                       capsValue *aimInputs,
                       int ibodyOffset,
                       meshInputStruct meshInput,
                       int boundingBoxIndex,
                       int createBL,
                       double globalBLSpacing,
                       double globalBLThickness,
                       double capsMeshLength,
                       const mapAttrToIndexStruct *groupMap,
                       const mapAttrToIndexStruct *meshMap,
                       int numMeshProp,
                       meshSizingStruct *meshProp,
                       int numSurfaceMesh,
                       ego *surfaceMesh,
                       int *skipVolume,
                       AFLR_Grid *grid)
{
    int i, d; // Indexing

    int propIndex;
    int bodyIndex, state, np, ibody;

    // Command line variables
    int  aflr3_argc   = 1;    // Number of arguments
    char **aflr3_argv = NULL; // String arrays
    char *meshInputString = NULL;
    char *rest = NULL, *token = NULL;
    char aimFile[PATH_MAX];

    ego *copy_body_tess=NULL, context, body, model=NULL;
    ego *faces=NULL, *modelFaces=NULL;
    int numFace = 0, numModelFace = 0;
    int *faceBodyIndex = NULL;
    int *faceGroupIndex = NULL;
    int newID;

    const char *keyWord;
    int index;

    int itransp = 0;
    int *transpBody = NULL;

    int iface = 0, nface, meshIndex;
    const char *groupName = NULL;
    int nnode_face, *face_node_map = NULL;

    const char *pstring = NULL;
    const int *pints = NULL;
    const double *preals = NULL;
    int atype, n;

    char bodyNumber[42], attrname[128];
    FILE *fp = NULL;

    const char* bcType = NULL;
    INT_ nbl, nbldiff;

    // Declare AFLR3 grid generation variables.
    INT_ nzone = 0;
    INT_ zone = -1;

    DOUBLE_2D *u = NULL;

    INT_1D *Surf_Error_Flag= NULL;

    INT_4D *BG_Vol_Tet_Neigbors = NULL;
    INT_4D *BG_Vol_Tet_Connectivity = NULL;

    DOUBLE_3D *BG_Coordinates = NULL;
    DOUBLE_1D *BG_Spacing = NULL;
    DOUBLE_6D *BG_Metric = NULL;

    DOUBLE_3D *Source_Coordinates = NULL;
    DOUBLE_1D *Source_Spacing = NULL;
    DOUBLE_6D *Source_Metric = NULL;

    INT_ Number_of_BG_Nodes = 0;
    INT_ Number_of_BG_Vol_Tets = 0;

    INT_ Number_of_Source_Nodes = 0;

    // Declare main program variables.

    INT_ status= 0;
    INT_ Message_Flag = 1;

    INT_ bc_mod = 1; // use optional transparent BC modification
    INT_ nlist = 0;
    INT_2D *bclist = NULL;
    INT_2D *idlist = NULL;

    INT_ final_mesh;

    CHAR_UG_MAX Output_Case_Name;

    DOUBLE_1D *BG_U_Scalars = NULL;
    DOUBLE_6D *BG_U_Metrics = NULL;

    INT_ * bc_ids_vector = NULL;
    DOUBLE_1D *bl_ds_vector = NULL;
    DOUBLE_1D *bl_del_vector = NULL;

    void *ext_cad_data = NULL;
    egads_struct *ptr = NULL;

    UG_Param_Struct *AFLR3_Param_Struct_Ptr = NULL;
    CHAR_UG_MAX Output_Grid_File_Name;
    char analysisPath[PATH_MAX];
    char currentPath[PATH_MAX];

    AFLR_Grid grid_c;

    initialize_AFLR_Grid(&grid_c);

    // Set and register program parameter functions.

    ug_set_prog_param_code (-3);

    ug_set_prog_param_function1 (ug_initialize_aflr_param);
    ug_set_prog_param_function1 (ug_gq_initialize_param); // optional
    ug_set_prog_param_function2 (aflr3_initialize_param);
    ug_set_prog_param_function2 (aflr3_anbl3_initialize_param);
    ug_set_prog_param_function2 (ice3_initialize_param);
    ug_set_prog_param_function2 (ug3_qchk_initialize_param); // optional

    // Register routines for BL mode
#ifdef _ENABLE_BL_
    aflr3_anbl3_register_grid_generator (anbl3_grid_generator);
    aflr3_anbl3_register_initialize_param (anbl3_initialize_param);
    aflr3_anbl3_register_be_set_surf_edge_data (anbl3_be_set_surf_edge_data);
    aflr3_anbl3_register_be_get_surf_edge_data (anbl3_be_get_surf_edge_data);
    aflr3_anbl3_register_be_set_surf_edge_data_def (anbl3_be_set_surf_edge_data_def);
    aflr3_anbl3_register_be_free_data (anbl3_be_free_data);
#endif

    // Register external routines for evaluation of the distribution function,
    // metrics and transformation vectors.

    // If external routines are used then they must be registered prior to calling
    // aflr3_grid_generator. The external evaluation routines must be registered
    // using dftr3_register_eval, either dftr3_register_eval_inl or
    // dftr3_register_eval_inl_flags, and possibly dftr3_register_eval_free.
    // Either the full initialization routine must be registered using
    // dftr3_register_eval_inl or the flag intialization routine must be
    // registered using dftr3_register_eval_inl_flag_data. Do not specify both.
    // If your initialization routine allocates and retains memory for subsequent
    // evaluations of the sizing function then a routine that frees that memory at
    // completion of grid generation or when a fatal error occurs must be
    // registered using dftr3_register_eval_free.

    dftr3_register_eval (dftr3_test_eval);
    dftr3_register_eval_inl (dftr3_test_eval_inl);

    // Register AFLR4-EGADS routines for CAD related setup & cleanup,
    // cad evaluation, cad bounds and generating boundary edge grids.

    // these calls are in aflr4_main_register - if that changes then these
    // need to change
    aflr4_register_auto_cad_geom_setup (egads_auto_cad_geom_setup);
    aflr4_register_cad_geom_create_tess (egads_aflr4_create_tess);
    aflr4_register_cad_geom_data_cleanup (egads_cad_geom_data_cleanup);
    aflr4_register_cad_geom_file_read (egads_cad_geom_file_read);
    aflr4_register_cad_geom_file_write (egads_cad_geom_file_write);
    aflr4_register_cad_geom_reset_attr (egads_cad_geom_reset_attr);
    aflr4_register_cad_geom_setup (egads_cad_geom_setup);
    aflr4_register_cad_tess_to_dgeom (egads_aflr4_tess_to_dgeom);
    aflr4_register_get_mzone_bedge_data (egads_get_mzone_bedge_data);
    aflr4_register_merge_mzone_data (egads_merge_mzone_data);
    aflr4_register_set_ext_cad_data (egads_set_ext_cad_data);
    aflr4_register_set_mzone_data (egads_set_mzone_data);
    aflr4_register_set_mzone_edge_data (egads_set_mzone_edge_data);

    dgeom_register_cad_eval_curv_at_uv (egads_eval_curv_at_uv);
    dgeom_register_cad_eval_edge_arclen (egads_eval_edge_arclen);
    dgeom_register_cad_eval_uv_bounds (egads_eval_uv_bounds);
    dgeom_register_cad_eval_uv_at_t (egads_eval_uv_at_t);
    dgeom_register_cad_eval_uv_at_xyz (egads_eval_uv_at_xyz);
    dgeom_register_cad_eval_xyz_at_t (egads_eval_xyz_at_u);
    dgeom_register_cad_eval_xyz_at_uv (egads_eval_xyz_at_uv);
    dgeom_register_discrete_eval_xyz_at_t (surfgen_discrete_eval_xyz_at_t);


    status = ug_add_new_arg (&aflr3_argv, (char*)"allocate_and_initialize_argv");
    AFLR_STATUS(aimInfo, status);

    // Set other command options
    if (createBL == (int) true) {
        status = ug_add_flag_arg ((char*)"mbl=1", &aflr3_argc, &aflr3_argv);
        AFLR_STATUS(aimInfo, status);

        if (aimInputs[BL_Max_Layers-1].nullVal == NotNull) {
          nbl = aimInputs[BL_Max_Layers-1].vals.integer;
          status = ug_add_flag_arg ("nbl", &aflr3_argc, &aflr3_argv); AFLR_STATUS(aimInfo, status);
          status = ug_add_int_arg (  nbl , &aflr3_argc, &aflr3_argv); AFLR_STATUS(aimInfo, status);
        }

        if (aimInputs[BL_Max_Layer_Diff-1].nullVal == NotNull) {
          nbldiff = aimInputs[BL_Max_Layer_Diff-1].vals.integer;
          status = ug_add_flag_arg ("nbldiff", &aflr3_argc, &aflr3_argv); AFLR_STATUS(aimInfo, status);
          status = ug_add_int_arg (  nbldiff , &aflr3_argc, &aflr3_argv); AFLR_STATUS(aimInfo, status);
        }

        status = ug_add_flag_arg ((char*)"mblelc=2", &aflr3_argc, &aflr3_argv);
        AFLR_STATUS(aimInfo, status);

    } else {
        status = ug_add_flag_arg ((char*)"mbl=0", &aflr3_argc, &aflr3_argv);
        AFLR_STATUS(aimInfo, status);
    }

    status = ug_add_flag_arg ((char*)"mrecm=3" , &aflr3_argc, &aflr3_argv);
    AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ((char*)"mrecqm=3", &aflr3_argc, &aflr3_argv);
    AFLR_STATUS(aimInfo, status);

    // Parse input string
    if (meshInput.aflr3Input.meshInputString != NULL) {

        rest = meshInputString = EG_strdup(meshInput.aflr3Input.meshInputString);
        while ((token = strtok_r(rest, " ", &rest))) {
            status = ug_add_flag_arg (token, &aflr3_argc, &aflr3_argv);
            if (status != 0) {
                printf("Error: Failed to parse input string: %s\n", token);
                if (meshInputString != NULL)
                    printf("Complete input string: %s\n", meshInputString);
                goto cleanup;
            }
        }
    }

    // setup input parameter structure

    status = aflr3_setup_param2 (3, 0, aflr3_argc, aflr3_argv,
                                 &AFLR3_Param_Struct_Ptr);
    AIM_STATUS(aimInfo, status);

    // Set meshInputs
    if (meshInput.quiet == 1) Message_Flag = 0;
    else Message_Flag = 1;

    // find all bodies with all AFLR_GBC TRANSP SRC/INTRNL

    AIM_ALLOC(transpBody, numSurfaceMesh, int, aimInfo, status);
    for (i = 0; i < numSurfaceMesh; i++) transpBody[i] = 0;

    for (bodyIndex = 0; bodyIndex < numSurfaceMesh; bodyIndex++) {

      status = EG_statusTessBody(surfaceMesh[bodyIndex], &body, &state, &np);
      AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(body, NULL, FACE, &numFace, &faces);
      AIM_STATUS(aimInfo, status);
      AIM_NOTNULL(faces, aimInfo, status);

      for (iface = 0; iface < numFace; iface++) {

        bcType = NULL;

        // check to see if AFLR_GBC is already set
        status = EG_attributeRet(faces[iface], "AFLR_GBC", &atype, &n,
                                 &pints, &preals, &pstring);
        if (status == CAPS_SUCCESS) {
          if (atype != ATTRSTRING) {
            AIM_ERROR(aimInfo, "AFLR_GBC on Body %d Face %d must be a string!", bodyIndex+1, iface+1);
            status = CAPS_BADVALUE;
            goto cleanup;
          }
          bcType = pstring;
        }

        status = retrieve_CAPSMeshAttr(faces[iface], &groupName);
        if (status == CAPS_SUCCESS) {

          status = get_mapAttrToIndexIndex(meshMap, groupName, &meshIndex);
          AIM_STATUS(aimInfo, status);

          for (propIndex = 0; propIndex < numMeshProp; propIndex++) {
            if (meshIndex != meshProp[propIndex].attrIndex) continue;

            // If bcType specified in meshProp
            if ((meshProp[propIndex].bcType != NULL)) {
              bcType = meshProp[propIndex].bcType;
            }
            break;
          }
        }

        // check to see if all faces on a body are TRANSP SRC/INTRNL
        if (bcType != NULL &&
            strncasecmp(bcType, "TRANSP_SRC_UG3_GBC", 18) == 0) {
          if (transpBody[bodyIndex] == -1) {
            AIM_ERROR(aimInfo, "Body %d has mixture of TRANSP_INTRNL_UG3_GBC/TRANSP_SRC_UG3_GBC and other BCs!");
            status = CAPS_BADTYPE;
            goto cleanup;
          }
          transpBody[bodyIndex] = 1;
        } else {
          if (transpBody[bodyIndex] == 1) {
            AIM_ERROR(aimInfo, "Body %d has mixture of TRANSP_INTRNL_UG3_GBC/TRANSP_SRC_UG3_GBC and other BCs!");
            status = CAPS_BADTYPE;
            goto cleanup;
          }
          transpBody[bodyIndex] = -1;
        }
      }
      AIM_FREE(faces);
    }

    AIM_ALLOC(copy_body_tess, 2*numSurfaceMesh, ego, aimInfo, status);

    // Copy bodies making sure TRANSP bodies are last
    ibody = 0;
    for (itransp = 1; itransp >= -1; itransp -= 2) {
      for (bodyIndex = 0; bodyIndex < numSurfaceMesh; bodyIndex++) {
        if (transpBody[bodyIndex] == itransp) continue;

        status = EG_statusTessBody(surfaceMesh[bodyIndex], &body, &state, &np);
        AIM_STATUS(aimInfo, status);

        status = EG_copyObject(body, NULL, &copy_body_tess[ibody]);
        AIM_STATUS(aimInfo, status);

        status = EG_copyObject(surfaceMesh[bodyIndex], copy_body_tess[ibody],
                               &copy_body_tess[numSurfaceMesh+ibody]);
        AIM_STATUS(aimInfo, status);

        status = EG_getBodyTopos(copy_body_tess[ibody], NULL, FACE, &numFace, &faces);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(faces, aimInfo, status);

        // Build up an array of all faces in the model
        AIM_REALL(modelFaces, numModelFace+numFace, ego, aimInfo, status);
        for (i = 0; i < numFace; i++) modelFaces[numModelFace+i] = faces[i];
        AIM_REALL(faceBodyIndex, numModelFace+numFace, int, aimInfo, status);
        for (i = 0; i < numFace; i++) faceBodyIndex[numModelFace+i] = bodyIndex;

        AIM_REALL(faceGroupIndex, numModelFace+numFace, int, aimInfo, status);
        for (i = 0; i < numFace; i++) {
          status = retrieve_CAPSGroupAttr(faces[i], &groupName);
          if (status == EGADS_NOTFOUND) {
            AIM_ERROR(aimInfo, "No capsGroup found on Face %d of body %d", i+1,bodyIndex+1);
            print_AllAttr(aimInfo, faces[i]);
            goto cleanup;
          } else
            AIM_STATUS(aimInfo, status);

          status = get_mapAttrToIndexIndex(groupMap, groupName, &faceGroupIndex[numModelFace+i]);
          AIM_STATUS(aimInfo, status);
        }

        AIM_FREE(faces);
        numModelFace += numFace;
        ibody++;
      }
    }

    // Map model face ID to property index

    bc_ids_vector = (INT_      *) ug_malloc (&status, (numModelFace) * sizeof (INT_));
    bl_ds_vector  = (DOUBLE_1D *) ug_malloc (&status, (numModelFace) * sizeof (DOUBLE_1D));
    bl_del_vector = (DOUBLE_1D *) ug_malloc (&status, (numModelFace) * sizeof (DOUBLE_1D));

    if (status != 0) {
      AIM_ERROR(aimInfo, "AFLR memory allocation error");
      status = EGADS_MALLOC;
      goto cleanup;
    }

    AIM_NOTNULL(modelFaces, aimInfo, status);
    AIM_NOTNULL(faceBodyIndex, aimInfo, status);

    for (iface = 0; iface < numModelFace; iface++) {

      bc_ids_vector[iface] = iface+1;
      bl_ds_vector[iface]  = globalBLSpacing*capsMeshLength;
      bl_del_vector[iface] = globalBLThickness*capsMeshLength;

      // check to see if AFLR_GBC is already set
      status = EG_attributeRet(modelFaces[iface], "AFLR_GBC", &atype, &n,
                               &pints, &preals, &pstring);
      if (status == EGADS_NOTFOUND) {
        bcType = (bl_ds_vector[iface] != 0 && bl_del_vector[iface] != 0 &&
                 faceBodyIndex[iface] != boundingBoxIndex) ? "-STD_UG3_GBC" : "STD_UG3_GBC";

        status = EG_attributeAdd(modelFaces[iface], "AFLR_GBC", ATTRSTRING, 0,
                                 NULL, NULL, bcType);
        AIM_STATUS(aimInfo, status);
      }

      status = retrieve_CAPSMeshAttr(modelFaces[iface], &groupName);
      if (status == CAPS_SUCCESS) {

        status = get_mapAttrToIndexIndex(meshMap, groupName, &meshIndex);
        AIM_STATUS(aimInfo, status);

        for (propIndex = 0; propIndex < numMeshProp; propIndex++) {
          if (meshIndex != meshProp[propIndex].attrIndex) continue;

          bl_ds_vector[iface]  = meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
          bl_del_vector[iface] = meshProp[propIndex].boundaryLayerThickness*capsMeshLength;

          status = EG_attributeRet(modelFaces[iface], "AFLR_GBC", &atype, &n,
                                   &pints, &preals, &pstring);
          AIM_STATUS(aimInfo, status);

          bcType = (bl_ds_vector[iface] != 0 && bl_del_vector[iface] != 0 &&
                    faceBodyIndex[iface] != boundingBoxIndex) ? "-STD_UG3_GBC" : pstring;

          // If bcType specified in meshProp
          if ((meshProp[propIndex].bcType != NULL)) {

            if      (strncasecmp(meshProp[propIndex].bcType, "Farfield"             ,  8) == 0 ||
                     strncasecmp(meshProp[propIndex].bcType, "Freestream"           , 10) == 0 ||
                     strncasecmp(meshProp[propIndex].bcType, "FARFIELD_UG3_GBC"     , 16) == 0)
              bcType = "FARFIELD_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "Viscous"              ,  7) == 0 ||
                     strncasecmp(meshProp[propIndex].bcType, "-STD_UG3_GBC"         , 12) == 0 ||
                (meshProp[propIndex].boundaryLayerSpacing > 0 &&
                    meshProp[propIndex].boundaryLayerThickness > 0))
              bcType = "-STD_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "Inviscid"             ,  8) == 0 ||
                     strncasecmp(meshProp[propIndex].bcType, "STD_UG3_GBC"          , 11) == 0)
              bcType = "STD_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "Symmetry"             ,  8) == 0 ||
                     strncasecmp(meshProp[propIndex].bcType, "BoundaryLayerIntersect",22) == 0 ||
                     strncasecmp(meshProp[propIndex].bcType, "BL_INT_UG3_GBC"       , 14) == 0)
              bcType = "BL_INT_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_SRC_UG3_GBC"   , 18) == 0)
              bcType = "TRANSP_SRC_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_BL_INT_UG3_GBC", 21) == 0)
              bcType = "TRANSP_BL_INT_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_UG3_GBC"       , 14) == 0)
              bcType = "TRANSP_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "-TRANSP_UG3_GBC"      , 15) == 0)
              bcType = "-TRANSP_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_INTRNL_UG3_GBC", 20) == 0)
              bcType = "TRANSP_INTRNL_UG3_GBC";
            else if (strncasecmp(meshProp[propIndex].bcType, "FIXED_BL_INT_UG3_GBC" , 19) == 0)
              bcType = "FIXED_BL_INT_UG3_GBC";
          }

          // Set face BC flag on the copy of the body
          if (bcType != pstring) {
            status = EG_attributeAdd(modelFaces[iface], "AFLR_GBC", ATTRSTRING, 0,
                                     NULL, NULL, bcType);
            AIM_STATUS(aimInfo, status);
          }

          break;
        }
      }
    }

    // create the model

    status = EG_getContext(copy_body_tess[0], &context);
    AIM_STATUS(aimInfo, status);

    status = EG_makeTopology(context, NULL, MODEL, 2*numSurfaceMesh, NULL, numSurfaceMesh,
                             copy_body_tess, NULL, &model);
    AIM_STATUS(aimInfo, status);

    if (createBL == (int)true) {

      status = ug_add_flag_arg ("BC_IDs", &aflr3_argc, &aflr3_argv);                             AFLR_STATUS(aimInfo, status);
      status = ug_add_int_vector_arg (numModelFace, bc_ids_vector, &aflr3_argc, &aflr3_argv);    AFLR_STATUS(aimInfo, status);
      status = ug_add_flag_arg ("BL_DS", &aflr3_argc, &aflr3_argv);                              AFLR_STATUS(aimInfo, status);
      status = ug_add_double_vector_arg (numModelFace, bl_ds_vector, &aflr3_argc, &aflr3_argv);  AFLR_STATUS(aimInfo, status);
      status = ug_add_flag_arg ("BL_DEL", &aflr3_argc, &aflr3_argv);                             AFLR_STATUS(aimInfo, status);
      status = ug_add_double_vector_arg (numModelFace, bl_del_vector, &aflr3_argc, &aflr3_argv); AFLR_STATUS(aimInfo, status);

    }

    // Set memory reduction output file format flag parameter, mpfrmt.
    //
    // If mpfrmt=1,2,3,4,5 then the output grid file is created within AFLR3 and
    // only the surface grid is passed back from AFLR3.
    //
    // If mpfrmt=0 then the output grid file is not created within AFLR3 and
    // the complete volume grid is passed back from AFLR3.
    //

    // Set write mesh flag - do not write out the mesh internally
    status = ug_add_flag_arg ((char *) "mpfrmt=0", &aflr3_argc, &aflr3_argv);
    AFLR_STATUS(aimInfo, status);

    status = ug_add_flag_arg ((char *) "mmsg", &aflr3_argc, &aflr3_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_int_arg(Message_Flag, &aflr3_argc, &aflr3_argv);      AFLR_STATUS(aimInfo, status);

    // note that if mpfrmt is not set to 0 and BL mesh generation is on then only
    // the surface mesh is returned

    // check that all the inputs

    status = ug_check_prog_param(aflr3_argv, aflr3_argc, Message_Flag);
    AFLR_STATUS(aimInfo, status);

    // set CAD geometry data structure

    status = aflr43_tess_to_dgeom (Message_Flag, model);
    AFLR_STATUS(aimInfo, status);

#if DUMP_DEBUG
    remove("aflr3t_in_debug.egads");
    EG_saveModel(egads_get_model (0), "aflr3t_in_debug.egads");
#endif


//#define DUMP_TECPLOT_DEBUG_FILE
#ifdef DUMP_TECPLOT_DEBUG_FILE
    {
      int surf = 0;
      int numEdge = 0;
      int numSurface = 0;
      int numTriFace = 0;
      int numNodes = 0;
      int numQuadFace = 0;
      int glueId; // Id of glued composite surface

      // AFRL4 output arrays
      INT_1D *bcFlag = NULL;
      INT_1D *ieFlag = NULL;
      INT_1D *idFlag = NULL;
      INT_2D *edgeCon = NULL;
      INT_3D *triCon = NULL;
      INT_4D *quadCon = NULL;
      DOUBLE_2D *uv = NULL;
      DOUBLE_3D *xyz = NULL;

      // Get output id index (glue-only composite)
      dgeom_def_get_idef (0, &glueId);

      FILE *fp = aim_fopen(aimInfo, "aflr4_debug.tec", "w");
      fprintf(fp, "VARIABLES = X, Y, Z, u, v\n");

      numSurface = dgeom_get_ndef(); // Get number of surfaces meshed

      for (surf = 0; surf < numSurface ; surf++) {

        if (surf+1 == glueId) continue;

        status = aflr4_get_def (surf+1,
                                0,
                                &numEdge,
                                &numTriFace,
                                &numNodes,
                                &numQuadFace,
                                &bcFlag,
                                &ieFlag,
                                &idFlag,
                                &edgeCon,
                                &triCon,
                                &quadCon,
                                &uv,
                                &xyz);
        AFLR_STATUS(aimInfo, status);

        fprintf(fp, "ZONE T=\"def %d\" N=%d, E=%d, F=FEPOINT, ET=Quadrilateral, DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n",
                surf+1, numNodes, numTriFace+numQuadFace);
        for (int i = 0; i < numNodes; i++)
          fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e\n",
                  xyz[i+1][0], xyz[i+1][1], xyz[i+1][2], uv[i+1][0], uv[i+1][1]);
        for (int i = 0; i < numTriFace; i++)
          fprintf(fp, "%d %d %d %d\n", triCon[i+1][0], triCon[i+1][1], triCon[i+1][2], triCon[i+1][2]);
        for (int i = 0; i < numQuadFace; i++)
          fprintf(fp, "%d %d %d %d\n", quadCon[i+1][0], quadCon[i+1][1], quadCon[i+1][2], quadCon[i+1][3]);

        ug_free (bcFlag);  bcFlag = NULL;
        ug_free (ieFlag);  ieFlag = NULL;
        ug_free (idFlag);  idFlag = NULL;
        ug_free (edgeCon); edgeCon = NULL;
        ug_free (triCon);  triCon = NULL;
        ug_free (quadCon); quadCon = NULL;
        ug_free (uv);      uv = NULL;
        ug_free (xyz);     xyz = NULL;

      }
      fclose(fp); fp = NULL;

    }
#endif


    //===================== IF THERE ARE MULTIPLE ZONES (BODIES) WITH MATCHING FACES
    nzone = dgeom_get_nzone();

    zone = (nzone) ? 0: -1;
    //========================================================= THEN LOOP OVER ZONES

    //................................................  GET LIST OF GRID BCs AND IDs

    status = aflr43_bc_list (&bc_mod, &nlist, &bclist, &idlist);
    AFLR_STATUS(aimInfo, status);

    do {

    //------------------------------------------------ GET COPY OF SURFACE MESH DATA

      status = aflr43_tess_to_surf (model, Message_Flag, zone,
                                    &grid->Number_of_Surf_Edges,
                                    &grid->Number_of_Surf_Trias,
                                    &grid->Number_of_Surf_Quads,
                                    &grid->Number_of_Nodes,
                                    &grid->Edge_ID_Flag,
                                    &grid->Surf_Edge_Connectivity,
                                    &grid->Surf_Grid_BC_Flag,
                                    &grid->Surf_ID_Flag,
                                    &grid->Surf_Reconnection_Flag,
                                    &grid->Surf_Tria_Connectivity,
                                    &grid->Surf_Quad_Connectivity,
                                    &grid->Coordinates);
      AFLR_STATUS(aimInfo, status);

      //..................................... MODIFY SURFACE MESH TRANSPARENT GRID BCs

      aflr43_bc_mod (bc_mod,
                     grid->Number_of_Surf_Trias,
                     grid->Number_of_Surf_Quads,
                     nlist,
                     grid->Surf_Grid_BC_Flag, grid->Surf_ID_Flag, bclist, idlist);

      // ************************************************** GENERATE AFLR3 VOLUME MESH

      if (skipVolume[zone == -1 ? 0 : zone] == 0) {
          status = aflr3_vol_gen (aflr3_argc, aflr3_argv, -3, Message_Flag,
                                  &grid->Number_of_Surf_Edges,
                                  &grid->Number_of_Surf_Trias,
                                  &grid->Number_of_Surf_Quads,
                                  &grid->Number_of_BL_Vol_Tets,
                                  &grid->Number_of_Vol_Tets,
                                  &grid->Number_of_Vol_Pents_5,
                                  &grid->Number_of_Vol_Pents_6,
                                  &grid->Number_of_Vol_Hexs,
                                  &grid->Number_of_Nodes,
                                  &Number_of_BG_Vol_Tets,
                                  &Number_of_BG_Nodes,
                                  &Number_of_Source_Nodes,
                                  &grid->Edge_ID_Flag,
                                  &grid->Surf_Edge_Connectivity,
                                  &grid->Surf_Grid_BC_Flag,
                                  &grid->Surf_ID_Flag,
                                  &Surf_Error_Flag,
                                  &grid->Surf_Reconnection_Flag,
                                  &grid->Surf_Tria_Connectivity,
                                  &grid->Surf_Quad_Connectivity,
                                  &grid->Vol_ID_Flag,
                                  &grid->Vol_Tet_Connectivity,
                                  &grid->Vol_Pent_5_Connectivity,
                                  &grid->Vol_Pent_6_Connectivity,
                                  &grid->Vol_Hex_Connectivity,
                                  &BG_Vol_Tet_Neigbors,
                                  &BG_Vol_Tet_Connectivity,
                                  &grid->Coordinates,
                                  &grid->BL_Normal_Spacing,
                                  &grid->BL_Thickness,
                                  &BG_Coordinates,
                                  &BG_Spacing,
                                  &BG_Metric,
                                  &Source_Coordinates,
                                  &Source_Spacing,
                                  &Source_Metric);

          if (status != 0) {
              strcpy (Output_Case_Name, "debug");
              ug3_write_surf_grid_error_file (Output_Case_Name,
                                              status,
                                              grid->Number_of_Nodes,
                                              grid->Number_of_Surf_Trias,
                                              Surf_Error_Flag,
                                              grid->Surf_Grid_BC_Flag,
                                              grid->Surf_ID_Flag,
                                              grid->Surf_Tria_Connectivity,
                                              grid->Coordinates);

              strcpy(aimFile, "aflr3_surf_debug.tec");

              fp = fopen(aimFile, "w");
              if (fp == NULL) goto cleanup;
              fprintf(fp, "VARIABLES = X, Y, Z, BC, ID\n");

              fprintf(fp, "ZONE N=%d, E=%d, F=FEBLOCK, ET=Triangle\n",
                      grid->Number_of_Nodes, grid->Number_of_Surf_Trias);
              fprintf(fp, ", VARLOCATION=([1,2,3]=NODAL,[4,5]=CELLCENTERED)\n");
              // write nodal coordinates
              if (grid->Coordinates != NULL)
                  for (d = 0; d < 3; d++) {
                      for (i = 0; i < grid->Number_of_Nodes; i++) {
                          if (i % 5 == 0) fprintf( fp, "\n");
                          fprintf(fp, "%22.15e ", grid->Coordinates[i+1][d]);
                      }
                      fprintf(fp, "\n");
                  }

              if (grid->Surf_Grid_BC_Flag != NULL)
                  for (i = 0; i < grid->Number_of_Surf_Trias; i++) {
                      if (i % 5 == 0) fprintf( fp, "\n");
                      fprintf(fp, "%d ", grid->Surf_Grid_BC_Flag[i+1]);
                  }

              if (grid->Surf_ID_Flag != NULL)
                  for (i = 0; i < grid->Number_of_Surf_Trias; i++) {
                      if (i % 5 == 0) fprintf( fp, "\n");
                      if ((Surf_Error_Flag != NULL) && (Surf_Error_Flag[i+1] < 0))
                          fprintf(fp, "-1 ");
                      else
                          fprintf(fp, "%d ", grid->Surf_ID_Flag[i+1]);
                  }

              // cell connectivity
              if (grid->Surf_Tria_Connectivity != NULL)
                  for (i = 0; i < grid->Number_of_Surf_Trias; i++)
                      fprintf(fp, "%d %d %d\n", grid->Surf_Tria_Connectivity[i+1][0],
                                                grid->Surf_Tria_Connectivity[i+1][1],
                                                grid->Surf_Tria_Connectivity[i+1][2]);
      /*@-dependenttrans@*/
              fclose(fp); fp = NULL;
      /*@+dependenttrans@*/
              status = CAPS_EXECERR;
              AIM_ERROR(aimInfo, "AFLR3 Grid generation error. The input surfaces mesh has been written to: %s", aimFile);
              goto cleanup;
          }
      }

      ug_free(Surf_Error_Flag);
      Surf_Error_Flag=NULL;

      if (nzone) {

        final_mesh = (zone == nzone-1) ? 1: 0;

        status = ug3_merge_mesh (final_mesh,
                                 &grid->Number_of_Surf_Edges,
                                 &grid->Number_of_Surf_Trias,
                                 &grid->Number_of_Surf_Quads,
                                 &grid->Number_of_BL_Vol_Tets,
                                 &grid->Number_of_Vol_Tets,
                                 &grid->Number_of_Vol_Pents_5,
                                 &grid->Number_of_Vol_Pents_6,
                                 &grid->Number_of_Vol_Hexs,
                                 &grid->Number_of_Nodes,
                                 zone+1,
                                 &grid->Edge_ID_Flag,
                                 &grid->Surf_Edge_Connectivity,
                                 &grid->Surf_Grid_BC_Flag,
                                 &grid->Surf_ID_Flag,
                                 &grid->Surf_Reconnection_Flag,
                                 &grid->Surf_Tria_Connectivity,
                                 &grid->Surf_Quad_Connectivity,
                                 &grid->Vol_ID_Flag,
                                 &grid->Vol_Tet_Connectivity,
                                 &grid->Vol_Pent_5_Connectivity,
                                 &grid->Vol_Pent_6_Connectivity,
                                 &grid->Vol_Hex_Connectivity,
                                 &grid->Coordinates,
                                 &grid->BL_Normal_Spacing,
                                 &grid->BL_Thickness,

                                 &grid_c.Number_of_Surf_Edges,
                                 &grid_c.Number_of_Surf_Trias,
                                 &grid_c.Number_of_Surf_Quads,
                                 &grid_c.Number_of_BL_Vol_Tets,
                                 &grid_c.Number_of_Vol_Tets,
                                 &grid_c.Number_of_Vol_Pents_5,
                                 &grid_c.Number_of_Vol_Pents_6,
                                 &grid_c.Number_of_Vol_Hexs,
                                 &grid_c.Number_of_Nodes,
                                 &grid_c.Edge_ID_Flag,
                                 &grid_c.Surf_Edge_Connectivity,
                                 &grid_c.Surf_Grid_BC_Flag,
                                 &grid_c.Surf_ID_Flag,
                                 &grid_c.Surf_Reconnection_Flag,
                                 &grid_c.Surf_Tria_Connectivity,
                                 &grid_c.Surf_Quad_Connectivity,
                                 &grid_c.Vol_ID_Flag,
                                 &grid_c.Vol_Tet_Connectivity,
                                 &grid_c.Vol_Pent_5_Connectivity,
                                 &grid_c.Vol_Pent_6_Connectivity,
                                 &grid_c.Vol_Hex_Connectivity,
                                 &grid_c.Coordinates,
                                 &grid_c.BL_Normal_Spacing,
                                 &grid_c.BL_Thickness);
        AFLR_STATUS(aimInfo, status);
      }

      //======================================================= END OF LOOP OVER ZONES

      zone++;

    } while (zone < nzone);

    //======================================== GENERATE TESS FROM AFLR3 SURFACE MESH

    status = aflr43_tess("",
                         Message_Flag, 0, model,
                         grid->Number_of_Surf_Edges,
                         grid->Number_of_Surf_Trias,
                         grid->Number_of_Surf_Quads,
                         grid->Edge_ID_Flag,
                         grid->Surf_ID_Flag,
                         grid->Surf_Edge_Connectivity,
                         grid->Surf_Tria_Connectivity,
                         grid->Surf_Quad_Connectivity,
                         grid->Coordinates);
    AFLR_STATUS(aimInfo, status);

    //....................................... RESET VOLUME MESH TRANSPARENT GRID BCs

    aflr43_bc_reset (bc_mod,
                     &grid->Number_of_Surf_Trias,
                     &grid->Number_of_Surf_Quads,
                     nlist, bclist, idlist,
                     grid->Surf_Grid_BC_Flag,
                     grid->Surf_ID_Flag,
                     grid->Surf_Reconnection_Flag,
                     grid->Surf_Quad_Connectivity,
                     grid->Surf_Tria_Connectivity);

    //......................................... MERGE DUPLICATE FACES IN VOLUME MESH

    if (nzone) {
      status = aflr43_merge_vol_mesh_faces (aflr3_argc, aflr3_argv, -3,
                                            &grid->Number_of_Surf_Trias,
                                            &grid->Number_of_Surf_Quads,
                                            grid->Number_of_Vol_Tets,
                                            grid->Number_of_Vol_Pents_5,
                                            grid->Number_of_Vol_Pents_6,
                                            grid->Number_of_Vol_Hexs,
                                            &grid->Number_of_Nodes,
                                            nlist,
                                            idlist,
                                            grid->Surf_Grid_BC_Flag,
                                            grid->Surf_ID_Flag,
                                            grid->Surf_Reconnection_Flag,
                                            grid->Surf_Quad_Connectivity,
                                            grid->Surf_Tria_Connectivity,
                                            grid->Vol_Tet_Connectivity,
                                            grid->Vol_Pent_5_Connectivity,
                                            grid->Vol_Pent_6_Connectivity,
                                            grid->Vol_Hex_Connectivity,
                                            grid->BL_Normal_Spacing,
                                            grid->BL_Thickness,
                                            grid->Coordinates);
      AFLR_STATUS(aimInfo, status);
    }


    (void)ug_get_char_param ("Output_Grid_File_Name", Output_Grid_File_Name, AFLR3_Param_Struct_Ptr);

    if (strcmp (Output_Grid_File_Name, "_null_") != 0) {

      // write volume mesh data file
      status = aim_file(aimInfo, Output_Grid_File_Name, aimFile);
      AIM_STATUS(aimInfo, status);

      status = aim_file(aimInfo, "", analysisPath);
      AIM_STATUS(aimInfo, status);

      // Get the current path
      (void) getcwd(currentPath, PATH_MAX);

      // Some AFLR writes are in the working directory, so switch to analysis PATH
      (void) chdir(analysisPath);

      status = ug_io_write_grid_file (aimFile,
                                      ug_abs (Message_Flag),
                                      grid->Number_of_BL_Vol_Tets,
                                      grid->Number_of_Nodes,
                                      grid->Number_of_Surf_Quads,
                                      grid->Number_of_Surf_Trias,
                                      grid->Number_of_Vol_Hexs,
                                      grid->Number_of_Vol_Pents_5,
                                      grid->Number_of_Vol_Pents_6,
                                      grid->Number_of_Vol_Tets,
                                      grid->Surf_Grid_BC_Flag,
                                      grid->Surf_ID_Flag,
                                      grid->Surf_Reconnection_Flag,
                                      grid->Surf_Quad_Connectivity,
                                      grid->Surf_Tria_Connectivity,
                                      grid->Vol_Hex_Connectivity,
                                      grid->Vol_ID_Flag,
                                      grid->Vol_Pent_5_Connectivity,
                                      grid->Vol_Pent_6_Connectivity,
                                      grid->Vol_Tet_Connectivity,
                                      grid->Coordinates,
                                      grid->BL_Normal_Spacing,
                                      grid->BL_Thickness);
      (void) chdir(currentPath); // switch back
      AFLR_STATUS(aimInfo, status);
    }

    status = egads_face_node_map_check (Message_Flag, grid->Coordinates);
    AFLR_STATUS(aimInfo, status);

    ext_cad_data = dgeom_get_ext_cad_data ();
    ptr = (egads_struct *) ext_cad_data;

    iface = 1;

    for (bodyIndex = 0, ibody = 0; bodyIndex < numSurfaceMesh; bodyIndex++) {
        if (transpBody[bodyIndex] == 1) continue;

        // set the file name to write the egads file
        snprintf(bodyNumber, 42, AFLR3TESSFILE, bodyIndex+ibodyOffset);
        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        status = EG_getBodyTopos(ptr->bodies[ibody], NULL, FACE, &nface, NULL);
        AIM_STATUS(aimInfo, status);

        for (i = 0; i < nface; i++, iface++) {

            egads_face_node_map_get_ptr (iface, &nnode_face, &face_node_map);

            // Add the unique indexing of the tessellation
            snprintf(attrname, 128, "face_node_map_%d",i+1);
            status = EG_attributeAdd(ptr->tess[ibody], attrname, ATTRINT,
                                     nnode_face, face_node_map+1, NULL, NULL); // face_node_map is index on [i+1]
            AIM_STATUS(aimInfo, status);
        }

        remove(aimFile);
        status = EG_saveTess(ptr->tess[ibody], aimFile);
        AIM_STATUS(aimInfo, status);
        ibody++;
    }

    // map the face index to the capsGroup index
    AIM_NOTNULL(grid->Surf_ID_Flag, aimInfo, status);
    AIM_NOTNULL(faceGroupIndex, aimInfo, status);

    newID = (int) true;
    for (i = 0; i < grid->Number_of_Surf_Trias + grid->Number_of_Surf_Quads; i++) {
      // Some attributes may be lost to interior domains
      if (newID == (int) true) {
        status = get_mapAttrToIndexKeyword(groupMap, faceGroupIndex[grid->Surf_ID_Flag[i+1]-1], &keyWord);
        AIM_STATUS(aimInfo, status);

        status = increment_mapAttrToIndexStruct(&grid->groupMap, keyWord);
        if (status != CAPS_SUCCESS && status != EGADS_EXISTS)
          AIM_STATUS(aimInfo, status);

        status = get_mapAttrToIndexIndex(&grid->groupMap, keyWord, &index);
        AIM_STATUS(aimInfo, status);
      }

      if (i < grid->Number_of_Surf_Trias + grid->Number_of_Surf_Quads-1)
        newID = grid->Surf_ID_Flag[i+2] != grid->Surf_ID_Flag[i+1] ? (int) true : (int) false;
      else
        newID = (int)true;
      grid->Surf_ID_Flag[i+1] = index;
    }

    // Remove the temporary grid created by AFLR
    remove(".tmp.b8.ugrid");

    status = CAPS_SUCCESS;
cleanup:

    // Free program arguements
    ug_free_argv(aflr3_argv); aflr3_argv = NULL;

    ug_free(Surf_Error_Flag);
    Surf_Error_Flag= NULL;

    ug_free(u);

    ug_free (BG_Vol_Tet_Neigbors);
    ug_free (BG_Vol_Tet_Connectivity);
    ug_free (BG_Coordinates);
    ug_free (BG_Spacing);
    ug_free (BG_Metric);

    BG_Vol_Tet_Neigbors = NULL;
    BG_Vol_Tet_Connectivity = NULL;
    BG_Coordinates = NULL;
    BG_Spacing = NULL;
    BG_Metric = NULL;

    ug_free(BG_U_Scalars);
    ug_free(BG_U_Metrics);

    BG_U_Scalars = NULL;
    BG_U_Metrics = NULL;

    ug_io_free_node(Source_Coordinates, Source_Spacing, Source_Metric);

    Source_Coordinates = NULL;
    Source_Spacing = NULL;
    Source_Metric = NULL;

    ug_free (bclist);
    ug_free (idlist);

    ug_free (bc_ids_vector);
    ug_free (bl_ds_vector);
    ug_free (bl_del_vector);

    ug_free_param (AFLR3_Param_Struct_Ptr);

/*@-mustfreefresh@*/
    bc_ids_vector = NULL;
    bl_ds_vector = NULL;
    bl_del_vector = NULL;
/*@+mustfreefresh@*/

    if (ptr != NULL) {
      EG_free (ptr->bodies);
      EG_deleteObject (ptr->model);
    }

    // cleanup aflr4 structures
    aflr4_free_all(0);
    egads_face_node_map_free();

    // Shut off memory and file status monitors and close output file.
    // This is required for implementation!
    ug_shutdown ();

    AIM_FREE(meshInputString);
    AIM_FREE(copy_body_tess);
    AIM_FREE(modelFaces);
    AIM_FREE(faceBodyIndex);
    AIM_FREE(faceGroupIndex);
    AIM_FREE(transpBody);

    destroy_AFLR_Grid(&grid_c);

    if (fp != NULL) fclose(fp);

    return status;
}
