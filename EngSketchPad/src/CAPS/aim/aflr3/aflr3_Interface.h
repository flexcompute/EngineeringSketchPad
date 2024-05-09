// AFLR3 interface functions - Modified from functions provided with
//    AFLR3_LIB source (aflr3.c) written by David L. Marcum


#include <ug/UG_LIB.h>

enum aimInputs
{
  Proj_Name = 1,               /* index is 1-based */
  Mesh_Quiet_Flag,
  Mesh_Format,
  Mesh_Gen_Input_String,
  Multiple_Mesh,
  Mesh_Sizing,
  BL_Initial_Spacing,
  BL_Thickness,
  BL_Max_Layers,
  BL_Max_Layer_Diff,
  Surface_Mesh,
  NUMINPUT = Surface_Mesh     /* Total number of inputs */
};

#define AFLR3TESSFILE "aflr3_%d.eto"

typedef struct {

  // AFLR3 grid variables.
  INT_1D *Edge_ID_Flag;
  INT_1D *Surf_Grid_BC_Flag;
  INT_1D *Surf_ID_Flag;
  INT_1D *Surf_Reconnection_Flag;
  INT_2D *Surf_Edge_Connectivity;
  INT_3D *Surf_Tria_Connectivity;
  INT_4D *Surf_Quad_Connectivity;
  INT_1D *Vol_ID_Flag;
  INT_4D *Vol_Tet_Connectivity;
  INT_5D *Vol_Pent_5_Connectivity;
  INT_6D *Vol_Pent_6_Connectivity;
  INT_8D *Vol_Hex_Connectivity;

  DOUBLE_3D *Coordinates;

  DOUBLE_1D *BL_Normal_Spacing;
  DOUBLE_1D *BL_Thickness;

  INT_ Number_of_BL_Vol_Tets;
  INT_ Number_of_Nodes;
  INT_ Number_of_Surf_Edges;
  INT_ Number_of_Surf_Trias;
  INT_ Number_of_Surf_Quads;
  INT_ Number_of_Vol_Tets;
  INT_ Number_of_Vol_Pents_5;
  INT_ Number_of_Vol_Pents_6;
  INT_ Number_of_Vol_Hexs;

  mapAttrToIndexStruct groupMap;

} AFLR_Grid;


void initialize_AFLR_Grid(AFLR_Grid *grid);
void destroy_AFLR_Grid(AFLR_Grid *grid);
int append_AFLR_Grid(void *aimInfo,
                     AFLR_Grid* domain,
                     int zone,
                     AFLR_Grid* grid);

int write_AFLR_Grid(void *aimInfo,
                    const char *fileName,
                    AFLR_Grid *grid);

int aflr3_to_MeshStruct( const AFLR_Grid *grid,
                         meshStruct *genUnstrMesh);

extern
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
                       AFLR_Grid *grid);
