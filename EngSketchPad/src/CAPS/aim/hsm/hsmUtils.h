#include "meshTypes.h"  // Bring in mesh structures
#include "capsTypes.h"  // Bring in CAPS types
#include "feaTypes.h"   // Bring in fea structures

#include "hsmTypes.h" // Bring in hsm structures

// Initiate hsmMemory structure
int initiate_hsmMemoryStruct(hsmMemoryStruct *mem);

// Destroy hsmMemory structure
int destroy_hsmMemoryStruct(hsmMemoryStruct *mem);

// Allocate hsmMemory structure
int allocate_hsmMemoryStruct(int numNode, int numElement, int maxDim,
                             hsmMemoryStruct *mem);

// Initiate hsmTempMemory structure
int initiate_hsmTempMemoryStruct(hsmTempMemoryStruct *mem);

// Destroy hsmTempMemory structure
int destroy_hsmTempMemoryStruct(hsmTempMemoryStruct *mem);

// Allocate hsmTempMemory structure
int allocate_hsmTempMemoryStruct(int numNode, int maxValence, int maxDim,
                                 hsmTempMemoryStruct *mem);

// Convert an EGADS body to a boundary element model - disjointed at edges
int hsm_bodyToBEM(void *aimInfo,
                  ego    ebody,                        // (in)  EGADS Body
                  double paramTess[3],                 // (in)  Tessellation parameters
                  int    edgePointMin,                 // (in)  minimum points along any Edge
                  int    edgePointMax,                 // (in)  maximum points along any Edge
                  int    quadMesh,                     // (in)  only do tris-for faces
                  mapAttrToIndexStruct *attrMap,       // (in)  map from CAPSGroup names to indexes
                  mapAttrToIndexStruct *coordSystemMap,// (in)  map from CoordSystem names to indexes
                  mapAttrToIndexStruct *constraintMap, // (in)  map from CAPSConstraint names to indexes
                  mapAttrToIndexStruct *loadMap,       // (in)  map from CAPSLoad names to indexes
                  mapAttrToIndexStruct *transferMap,   // (in)  map from CAPSTransfer names to indexes
       /*@null@*/ mapAttrToIndexStruct *connectMap,    // (in)  map from CAPSConnect names to indexes
                  meshStruct *feaMesh);                // (out) FEA mesh structure

// Write hsm data to a Tecplot file
int hsm_writeTecplot(void *aimInfo,
                     char *projectName,
                     meshStruct feaMesh,
                     hsmMemoryStruct *hsmMemory,
                     int *permutation);

// Set global parameters in hsmMemory structure
int hsm_setGlobalParameter(feaProblemStruct feaProblem,
                           hsmMemoryStruct *hsmMemory);

// Set parameters in hsmMemory structure
int hsm_setSurfaceParameter(feaProblemStruct feaProblem, int permutation[],
                            hsmMemoryStruct *hsmMemory);
int hsm_setEdgeBCParameter(feaProblemStruct feaProblem, int permutation[],
                           hsmMemoryStruct *hsmMemory);
int hsm_setNodeBCParameter(feaProblemStruct feaProblem, int permutation[],
                           hsmMemoryStruct *hsmMemory);

// Generates the adjacency structure (non-zero matrix pattern) for HSM
extern
int hsm_Adjacency(meshStruct *feaMesh,
                  const int numJoint, // Number of joints
                  const int *kjoint,  // Joint connectivity
                  int *maxAdjacency,  // Max adjacency (columns in a row)
                  int **xadj_out,     // Pointers (indices) into adj (freeable)
                  int **adj_out);     // The adjacency lists (freeable)
