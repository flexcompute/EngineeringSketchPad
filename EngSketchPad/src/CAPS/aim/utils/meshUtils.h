// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Mesh related utility functions - Written by Dr. Ryan Durscher AFRL/RQVC

#include "meshTypes.h"  // Bring in mesh structures
#include "capsTypes.h"  // Bring in CAPS types
#include "cfdTypes.h"   // Bring in cfd structures
#include "feaTypes.h"   // Bring in fea structures
#include "miscTypes.h"  // Bring in miscellanous structures

#ifdef __cplusplus
extern "C" {
#endif

// extracts boundary elements and adds them to the surface mesh
int mesh_addTess2Dbc(void *aimInfo, meshStruct *surfaceMesh, mapAttrToIndexStruct *attrMap);

int mesh_bodyTessellation(void *aimInfo, ego tess, mapAttrToIndexStruct *attrMap,
                          int *numNodes, double *xyzCoord[],
                          int *numTriFace, int *triFaceConn[], int *triFaceCompID[], int *triFaceTopoID[],
                          int *numBndEdge, int *bndEdgeConn[], int *bndEdgeCompID[], int *bndEdgeTopoID[],
                          int *numNodeEle, int *nodeEleConn[], int *nodeEleCompID[], int *nodeEleTopoID[],
                          int *twoDMesh,
               /*@null@*/ const int *tessFaceQuadMap,
                          int *numQuadFace, int *quadFaceConn[], int *quadFaceCompID[], int *quadFaceTopoID[]);

// Create a surface mesh in meshStruct format using the EGADS body object
int mesh_surfaceMeshEGADSBody(void *aimInfo, ego body, double refLen, double tessParams[3], int quadMesh, ego *tess);

// Create a surface mesh in meshStruct format using the EGADS body tessellation
int mesh_surfaceMeshEGADSTess(void *aimInfo, meshStruct *surfMesh);

// Modify the EGADS body tessellation based on given inputs
int mesh_modifyBodyTess(int numMeshProp,
             /*@null@*/ meshSizingStruct meshProp[],
                        int minEdgePointGlobal,
                        int maxEdgePointGlobal,
                        int quadMesh,
                        double *refLen,
                        double tessParamGlobal[],
                        mapAttrToIndexStruct attrMap,
                        int numBody,
                        ego bodies[]);

int write_MAPBC(void *aimInfo,
                char *fname,
                int numBnds,
                int *bndIds,
                int *bndVals);

void get_Surface_Norm(double p1[3],
                      double p2[3],
                      double p3[3],
                      double norm[3]);

// Populate bndCondStruct boundary condition information - Boundary condition values get filled with 99
int populate_bndCondStruct_from_bcPropsStruct(cfdBoundaryConditionStruct *bcProps,
                                              bndCondStruct *bndConds);

// Populate bndCondStruct boundary condition information from attribute map - Boundary condition values get filled with 99
int populate_bndCondStruct_from_mapAttrToIndexStruct(mapAttrToIndexStruct *attrMap,
                                                     bndCondStruct *bndConds);

// Initiate (0 out all values and NULL all pointers) a bndCond in the bndCondStruct structure format
int initiate_bndCondStruct(bndCondStruct *bndCond);

// Destroy (0 out all values, free all arrays, and NULL all pointers) a bndCond in the bndCondStruct structure format
int destroy_bndCondStruct(bndCondStruct *bndCond);

// Populate a tetgenRegionsStruct regions data structure
int populate_regions(tetgenRegionsStruct* regions,
                     int length,
                     const capsTuple* tuples);

// Initialize a tetgenRegionsStruct regions data structure
int initiate_regions(tetgenRegionsStruct* regions);

// Destroy a tetgenRegionsStruct regions data structure
int destroy_regions(tetgenRegionsStruct* regions);

// Populate a tetgenHolesStruct holes data structure
int populate_holes(tetgenHolesStruct* holes,
                     int length,
                     const capsTuple* tuples);

// Initialize a tetgenHolesStruct holes data structure
int initiate_holes(tetgenHolesStruct* holes);

// Destroy a tetgenHolesStruct holes data structure
int destroy_holes(tetgenHolesStruct* holes);

// Initiate (0 out all values and NULL all pointers) a meshInput in the meshInputStruct structure format
int initiate_meshInputStruct(meshInputStruct *meshInput);

// Destroy (0 out all values and NULL all pointers) a meshInput in the meshInputStruct structure format
int destroy_meshInputStruct(meshInputStruct *meshInput);

// Initiate (0 out all values and NULL all pointers) a meshProp in the meshSizingStruct structure format
int initiate_meshSizingStruct (meshSizingStruct *meshProp);

// Destroy (0 out all values and NULL all pointers) a meshProp in the meshSizingStruct structure format
int destroy_meshSizingStruct (meshSizingStruct *meshProp);

// Fill meshProps in a meshBCStruct format with mesh boundary condition information from incoming Mesh Sizing Tuple
int mesh_getSizingProp(void *aimInfo,
                       int numTuple,
                       capsTuple meshBCTuple[],
                       mapAttrToIndexStruct *attrMap,
                       int *numMeshProp,
                       meshSizingStruct *meshProps[]);

// Initiate (0 out all values and NULL all pointers) a CFD mesh data in the cfdMeshDataStruct structure format
int initiate_cfdMeshDataStruct (cfdMeshDataStruct *data);

// Destroy (0 out all values and NULL all pointers) a CFD mesh data in the cfdMeshDataStruct structure format
int destroy_cfdMeshDataStruct (cfdMeshDataStruct *data);

// Copy a CFD mesh data in the cfdMeshDataStruct structure format
int copy_cfdMeshDataStruct (cfdMeshDataStruct *dataIn, cfdMeshDataStruct *dataOut);

// Initiate (0 out all values and NULL all pointers) a FEA mesh data in the feaMeshDataStruct structure format
int initiate_feaMeshDataStruct (feaMeshDataStruct *data);

// Destroy (0 out all values and NULL all pointers) a FEA mesh data in the feaMeshDataStruct structure format
int destroy_feaMeshDataStruct (feaMeshDataStruct *data);

// Copy a FEA mesh data in the feaMeshDataStruct structure format
int copy_feaMeshDataStruct (feaMeshDataStruct *dataIn, feaMeshDataStruct *dataOut);

// Initiate  (0 out all values and NULL all pointers) and allocate the analysisData void pointer. Creation selected based on type.
int initiate_analysisData(void **analysisData, meshAnalysisTypeEnum analysisType);

// Destroy (0 out all values and NULL all pointers) and free of the analysisData void pointer. Correct destroy function selected based on type.
int destroy_analysisData(void **analysisData, meshAnalysisTypeEnum analysisType);

// Initiate (0 out all values and NULL all pointers) a node data in the meshGeomDataStruct structure format
int initiate_meshGeomDataStruct(meshGeomDataStruct *geom);

// Destroy (0 out all values and NULL all pointers) a node data in the meshGeomDataStruct structure format
int destroy_meshGeomDataStruct(meshGeomDataStruct *geom);

// Copy geometry mesh data in the meshGeomDataStruct structure format
int copy_meshGeomDataStruct(meshGeomDataStruct *dataIn, meshGeomDataStruct *dataOut);

// Initiate (0 out all values and NULL all pointers) a node data in the meshNode structure format
int initiate_meshNodeStruct(meshNodeStruct *node, meshAnalysisTypeEnum meshAnalysisType);

// Destroy (0 out all values and NULL all pointers) a node data in the meshNode structure format
int destroy_meshNodeStruct(meshNodeStruct *node);

// Update/change the analysis data type in a meshNodeStruct
int change_meshNodeAnalysis(meshNodeStruct *node, meshAnalysisTypeEnum meshAnalysisType);

// Initiate (0 out all values and NULL all pointers) a element data in the meshElement structure format
int initiate_meshElementStruct(meshElementStruct *element, meshAnalysisTypeEnum meshAnalysisType);

// Destroy (0 out all values and NULL all pointers) a element data in the meshElementStruct structure format
int destroy_meshElementStruct(meshElementStruct *element);

// Update/change the analysis data type in a meshElementStruct
int change_meshElementAnalysis(meshElementStruct *element, meshAnalysisTypeEnum meshAnalysisType);

// Initiate (0 out all values and NULL all pointers) the meshQuickRef data in the meshQuickRefStruct structure format
int initiate_meshQuickRefStruct(meshQuickRefStruct *quickRef);

// Destroy (0 out all values and NULL all pointers) the meshQuickRef data in the meshQuickRefStruct structure format
int destroy_meshQuickRefStruct(meshQuickRefStruct *quickRef);

// Destroy (0 out all values and NULL all pointers) all nodes in the mesh
int destroy_meshNodes(meshStruct *mesh);

// Destroy (0 out all values and NULL all pointers) all elements in the mesh
int destroy_meshElements(meshStruct *mesh);

// Initiate (0 out all values and NULL all pointers) a mesh data in the meshStruct structure format
int initiate_meshStruct(meshStruct *mesh);

// Destroy (0 out all values and NULL all pointers) a mesh data in the meshStruct structure format
int destroy_meshStruct(meshStruct *mesh);

// Update/change the analysis data in a meshStruct
int change_meshAnalysis(meshStruct *mesh, meshAnalysisTypeEnum meshAnalysisType);

// Return the number of connectivity points based on type
int mesh_numMeshConnectivity(meshElementTypeEnum elementType);

// Return the number of connectivity points based on type of element provided
int mesh_numMeshElementConnectivity(meshElementStruct *element);

// Allocate mesh element connectivity array based on type
int mesh_allocMeshElementConnectivity(meshElementStruct *element);

// Retrieve the number of mesh element of a given type
int mesh_retrieveNumMeshElements(int numElement,
                                 meshElementStruct element[],
                                 meshElementTypeEnum elementType,
                                 int *numElementType);

// Retrieve the starting index of a given type -assume elements were put in order
int mesh_retrieveStartIndexMeshElements(int numElement,
                                        meshElementStruct element[],
                                        meshElementTypeEnum elementType,
                                        int *numElementType,
                                        int *startIndex);

// Retrieve list of mesh element of a given type - elementTypeList is freeable
int mesh_retrieveMeshElements(int numElement,
                              meshElementStruct element[],
                              meshElementTypeEnum elementType,
                              int *numElementType,
                              int *elementTypeList[]);

// Fill out the QuickRef lists for all element types
int mesh_fillQuickRefList( meshStruct *mesh);

// Make a copy of the analysis Data
int mesh_copyMeshAnalysisData(void *in, meshAnalysisTypeEnum analysisType, void *out);

// Make a copy of an element - may offset the element and connectivity indexing
int mesh_copyMeshElementStruct(meshElementStruct *in, int elementOffSetIndex, int connOffSetIndex, meshElementStruct *out);

// Make a copy of an node - may offset the node indexing
int mesh_copyMeshNodeStruct(meshNodeStruct *in, int nodeOffSetIndex, meshNodeStruct *out);

// Copy a mesh structure
int mesh_copyMeshStruct( meshStruct *in, meshStruct *out );

// Combine mesh structures
int mesh_combineMeshStruct(int numMesh, meshStruct mesh[], meshStruct *combineMesh );

// Write a mesh contained in the mesh structure in AFLR3 format (*.ugrid, *.lb8.ugrid, *.b8.ugrid)
int mesh_writeAFLR3(void *aimInfo,
                    char *fname,
                    int asciiFlag, // 0 for binary, anything else for ascii
                    meshStruct *mesh,
                    double scaleFactor); // Scale factor for coordinates

// Read a mesh into the mesh structure from an AFLR3 format (*.ugrid, *.lb8.ugrid, *.b8.ugrid)
int mesh_readAFLR3(void *aimInfo,
                   char *fname,
                   meshStruct *mesh,
                   double scaleFactor) ;// Scale factor for coordinates

// Write a mesh contained in the mesh structure in VTK format (*.vtk)
int mesh_writeVTK(void *aimInfo,
                  char *fname,
                  int asciiFlag, // 0 for binary, anything else for ascii
                  meshStruct *mesh,
                  double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in SU2 format (*.su2)
int mesh_writeSU2(void *aimInfo,
                  char *fname,
                  int asciiFlag, // 0 for binary, anything else for ascii
                  meshStruct *mesh,
                  int numBnds,
                  int bndID[],
                  double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in NASTRAN format (*.bdf)
int mesh_writeNASTRAN(void *aimInfo,
                      char *fname,
                      int asciiFlag, // 0 for binary, anything else for ascii
                      meshStruct *nasMesh,
                      feaFileTypeEnum gridFileType,
                      double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in Astros format (*.bdf)
int mesh_writeAstros(void *aimInfo,
                     char *fname,
                     int asciiFlag, // 0 for binary, anything else for ascii
                     meshStruct *mesh,
                     feaFileTypeEnum gridFileType,
                     int numDesignVariable, feaDesignVariableStruct feaDesignVariable[],
                     double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in STL format (*.stl)
int mesh_writeSTL(void *aimInfo,
                  char *fname,
                  int asciiFlag, // 0 for binary, anything else for ascii
                  meshStruct *mesh,
                  double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in Tecplot format (*.dat)
int mesh_writeTecplot(void *aimInfo,
                      const char *fname,
                      int asciiFlag,
                      meshStruct *mesh,
                      double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in Airfoil format (boundary edges only [Lines]) (*.af)
//  "Character Name"
//	x[0] y[0] x y coordinates
//	x[1] y[1]
//	 ...  ...
int mesh_writeAirfoil(void *aimInfo,
                      char *fname,
                      int asciiFlag,
                      meshStruct *mesh,
                      double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in FAST mesh format (*.msh)
int mesh_writeFAST(void *aimInfo,
                   char *fname,
                   int asciiFlag,
                   meshStruct *mesh,
                   double scaleFactor); // Scale factor for coordinates

// Write a mesh contained in the mesh structure in Abaqus mesh format (*_Mesh.inp)
int mesh_writeAbaqus(void *aimInfo,
                     char *fname,
                     int asciiFlag,
                     meshStruct *mesh,
                     double scaleFactor); // Scale factor for coordinates

// Extrude a surface mesh a single unit the length of extrusionLength - return a
// volume mesh, cell volume and left-handness is not checked
int extrude_SurfaceMesh(double extrusionLength,
                        int extrusionMarker,
                        meshStruct *surfaceMesh, meshStruct *volumeMesh);

// Retrieve the max valence and valence of each node in the mesh - currently only supports Triangle and Quadrilateral elements
int mesh_retrieveMaxValence(meshStruct *mesh,
                            int *maxValence,     // Max valence
                            int *nodeValence[]); // Array of valences for each node  (freeable)

// Look at the nodeID for each node and check to see if it is being used in the element connectivity; if not it is removed
// Note: that the nodeIDs for the nodes and element connectivity isn't changed, as such if using element connectivity to blindly
// access a given node this could lead to seg-faults!. mesh_nodeID2Array must be used to access the node array index.
int mesh_removeUnusedNodes(meshStruct *mesh);

// Constructs a map that maps from nodeID to the mesh->node array index
int mesh_nodeID2Array(const meshStruct *mesh,
                      int **n2a);

// Create a new mesh with topology tagged with capsIgnore being removed, if capsIgnore isn't found the mesh is simply copied.
int mesh_createIgnoreMesh(meshStruct *mesh, meshStruct *meshIgnore);

// Changes the analysisType of a mesh
int mesh_setAnalysisType(meshAnalysisTypeEnum analysisType, meshStruct *mesh);

// Find meshNodeStructs with `isMatch` function
// Returns array of borrowed pointers
int mesh_findNodes(meshStruct *mesh, 
                   int (*isMatch)(meshNodeStruct *, void *), 
                   void *isMatchArg, 
                   int *numFound, 
                   meshNodeStruct ***foundSet);

// Find meshElementStructs with `isMatch` function
// Returns array of borrowed pointers
int mesh_findElements(meshStruct *mesh, 
                      int (*isMatch)(meshElementStruct *, void *), 
                      void *isMatchArg, 
                      int *numFound, 
                      meshElementStruct ***foundSet);

// Find meshElementStructs with given groupName(s)
// Returns array of borrowed pointers
int mesh_findGroupElements(meshStruct *mesh,
                      mapAttrToIndexStruct *attrMap, 
                      int numGroupName, 
                      char **groupName, 
                      int *numGroupElement, 
                      meshElementStruct ***groupElementSet);

// General routine to do fill up a AIM capsDiscr data structure
int mesh_fillDiscr(char *tname, mapAttrToIndexStruct *groupMap,
                   int numBody, ego *tess, capsDiscr *discr);

#ifdef __cplusplus
}
#endif
