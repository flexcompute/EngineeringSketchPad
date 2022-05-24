// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

// Structures for general meshing - Written by Dr. Ryan Durscher AFRL/RQVC


#ifndef MESHTYPES_H
#define MESHTYPES_H

#include "egads.h"
#include "miscTypes.h"

typedef enum {UnknownDistribution, EvenDistribution, TanhDistribution } edgeDistributionEnum;

typedef enum {UnknownMeshElement, Node, Line, Triangle, Triangle_6, Quadrilateral, Quadrilateral_8, Tetrahedral, Tetrahedral_10, Pyramid, Prism, Hexahedral} meshElementTypeEnum;

typedef enum {UnknownMeshSubElement, ConcentratedMassElement, BarElement, BeamElement, ShellElement, ShearElement, MembraneElement} meshElementSubTypeEnum;

typedef enum {UnknownMeshAnalysis, MeshCFD, MeshStructure, MeshOrigami} meshAnalysisTypeEnum;
//typedef enum {UnknownMeshDimension, TwoDimensional, ThreeDimensional} meshDimensionalityEnum;
typedef enum {UnknownMeshType, Surface2DMesh, SurfaceMesh, VolumeMesh} meshTypeEnum;

// Container for boundary conditions
typedef struct {
    int numBND;
    int *bndID;
    int *bcVal;
} bndCondStruct;

typedef struct {
    int size;
    double* x;
    double* y;
    double* z;
    int* attribute;
    double* volume_constraint;
} tetgenRegionsStruct;

typedef struct {
    int size;
    double* x;
    double* y;
    double* z;
} tetgenHolesStruct;

// Container for Tetgen specific inputs
typedef struct {
    double meshQuality_rad_edge;// Tetgen: maximum radius-edge ratio
    double meshQuality_angle;   // Tetgen: minimum dihedral angle
    char  *meshInputString;     // Tetgen: Input string (optional) if NULL use default values
    int   verbose;              // 0 = False, anything else True - Verbose output from mesh generator
    int   ignoreSurfaceExtract; // 0 = False, anything else True - Don't extract the new surface mesh if
                                //      Steiner points are added.
    double meshTolerance; // Tetgen : mesh tolerance
    tetgenRegionsStruct regions;
    tetgenHolesStruct holes;
} tetgenInputStruct;

// Container for AFLR3 specific inputs
typedef struct {
    char  *meshInputString;    // AFLR3: Input string (optional) if NULL use default values
} aflr3InputStruct;

// Container for AFLR4 specific inputs
typedef struct {
    char  *meshInputString;    // AFLR4: Input string (optional) if NULL use default values
} aflr4InputStruct;

// Container for hoTess specific inputs
typedef struct {
    meshElementTypeEnum meshElementType;

    // inputs to EG_tessHOverts
    // number of vertices local to the elevated element
    int numLocalElevatedVerts;

    // weights of local verts relative to ref element
    // 2*numLocalElevatedVerts in length
    double *weightsLocalElevatedVerts;

    // number of internal elevated tris created per source triangle
    // negative number indicates quads (paired triangles)
    int numLocalElevatedTris;

    // local elevated triangle indices (1-bias)
    // 3*numLocalElevatedTris in length
    int *orderLocalElevatedTris;
} hoTessInputStruct;

// Container for meshing inputs
typedef struct {
    double paramTess[3]; /*a set of 3 parameters that drive the EDGE discretization and the
                           FACE triangulation. The first is the maximum length of an EDGE
                           segment or triangle side (in physical space). A zero is flag that allows
                           for any length. The second is a curvature-based value that looks
                           locally at the deviation between the centroid of the discrete object and
                           the underlying geometry. Any deviation larger than the input value will
                           cause the tessellation to be enhanced in those regions. The third is
                           the maximum interior dihedral angle (in degrees) between triangle
                           facets (or Edge segment tangents for a WIREBODY tessellation),
                           note that a zero ignores this phase*/

    int preserveSurfMesh; // 0 = False , anything else True - Use the body tessellation as the surface mesh

    int quiet;            // 0 = False , anything else True - No output from mesh generator
    char *outputFormat;   // Mesh output formats - AFLR3, TECPLOT, VTK, SU2
    char *outputFileName; // Filename prefix for mesh
    int outputASCIIFlag;  // 0 = Binary output, anything else for ASCII

    bndCondStruct bndConds; // Structure of boundary conditions
    tetgenInputStruct tetgenInput; // Structure of Tetgen specific inputs
    aflr3InputStruct aflr3Input; // Structure of AFLR3 specific inputs
    aflr4InputStruct aflr4Input; // Structure of AFLR4 specific inputs
    hoTessInputStruct hoTessInput; // Structure of hoTess specific inputs

} meshInputStruct;

typedef struct {
    char *name; // Attribute name
    int attrIndex;  // Attribute index

    int numEdgePoints; // Number of points along an edge

    edgeDistributionEnum edgeDistribution; // Distribution function along an edge

    double minSpacing;            // Minimum allowed spacing on EDGE/FACE
    double maxSpacing;            // Maximum allowed spacing on EDGE/FACE
    double avgSpacing;            // Average allowed spacing on EDGE

    double maxAngle;              // Maximum angle on an EDGE to control spacing
    double maxDeviation;          // Maximum deviation on an EDGE/FACE to control spacing
    double boundaryDecay;         // Decay of influence of the boundary spacing on the interior spacing

    double nodeSpacing;           // Node spacing at a NODE or end points of an EDGE
    double initialNodeSpacing[2]; // Initial node spacing along an edge

    int useTessParams; // Trigger to use the specified face tessellation parameters
    double tessParams[3]; // Tessellation parameters on a face.

    double boundaryLayerThickness;  // Boundary layer thickness on a face (3D) or edge (2D)
    double boundaryLayerSpacing;    // Boundary layer spacing on a face (3D) or edge (2D)
    int boundaryLayerMaxLayers;     // Maximum number of layers
    int boundaryLayerFullLayers;    // Number of complete layers
    double boundaryLayerGrowthRate; // Growth rate of the boundary layer

    char *bcType;       // Name of the meshing boundary condition type
    double scaleFactor; // Scaling factor to apply to an individual face
    double edgeWeight;  // Interpolation weight on edge mesh size between faces with large angles

} meshSizingStruct;

// Container to store geometric data
typedef struct {

    // These may need to change in the future

    // Node/Vertex geometry information

    double uv[2];
    double firstDerivative[6];

    int type; // The point type (-) Face local index, (0) Node, (+) Edge local index
    int topoIndex; // The point topological index (1 bias)

} meshGeomDataStruct;

// Container for mesh elements
typedef struct {
    meshElementTypeEnum elementType;

    int elementID;

    int markerID;

    // The point topological index (1 bias)
    // The topological type is dictated by elementType:
    //   Node - NODE
    //   Line - EDGE
    //   Triangle, Triangle_6, Quadrilateral, Quadrilateral_8 - FACE
    //   Tetrahedral, Tetrahedral_10, Pyramid, Prism, Hexahedral - Undefined
    int topoIndex;

    int *connectivity; // size[elementType-specific]

    meshAnalysisTypeEnum analysisType;
    void *analysisData;

} meshElementStruct;

// Container for mesh nodes
typedef struct {

    double xyz[3];

    int nodeID;

    meshAnalysisTypeEnum analysisType;
    void *analysisData;

    // Optional - store away geometric data for the element
    meshGeomDataStruct *geomData; // Must be separately allocated and initiated

} meshNodeStruct;

// Container for quickly referencing the mesh
typedef struct {

    int useStartIndex; // Use the start index reference
    int useListIndex; // Use the list of indexes

    // Number of elements per type
    int numNode;
    int numLine;
    int numTriangle;
    int numTriangle_6;
    int numQuadrilateral;
    int numQuadrilateral_8;
    int numTetrahedral;
    int numTetrahedral_10;
    int numPyramid;
    int numPrism;
    int numHexahedral;

    // If the element types are created in order - starting index in the element array of a particular element in type
    int startIndexNode;
    int startIndexLine;
    int startIndexTriangle;
    int startIndexTriangle_6;
    int startIndexQuadrilateral;
    int startIndexQuadrilateral_8;
    int startIndexTetrahedral;
    int startIndexTetrahedral_10;
    int startIndexPyramid;
    int startIndexPrism;
    int startIndexHexahedral;

    // Array of element indexes containing a specific element type
    int *listIndexNode; // size[numNode]
    int *listIndexLine; // size[numLine]
    int *listIndexTriangle; // size[numTriangle]
    int *listIndexTriangle_6; // size[numTriangle_6]
    int *listIndexQuadrilateral; // size[numQuadrilateral]
    int *listIndexQuadrilateral_8; // size[numQuadrilateral_8]
    int *listIndexTetrahedral; // size[numTetrahedral]
    int *listIndexTetrahedral_10; // size[numTetrahedral_10]
    int *listIndexPyramid; // size[numPyramid]
    int *listIndexPrism; // size[numPrism]
    int *listIndexHexahedral; // size[numHexahedral]

} meshQuickRefStruct;

// Container for an unstructured mesh
typedef struct  meshStruct meshStruct;
struct meshStruct {

    //meshDimensionalityEnum meshDimensionality;
    meshTypeEnum meshType;
    meshAnalysisTypeEnum analysisType;

    int numNode;
    meshNodeStruct *node; // size[numNode]

    int numElement;
    meshElementStruct *element; // size[numElement]

    ego egadsTess; // EGADS body tessellation

    int numReferenceMesh; // Number of reference meshes
    meshStruct *referenceMesh; // Pointers to other meshes, should be freed but not individual references, size[numReferenceMesh]

    meshQuickRefStruct meshQuickRef;

    // Attribute to index map that maps CAPSGroup names to indexes
    mapAttrToIndexStruct groupMap;
};

// Container for mesh data relevant for CFD analysis
typedef struct {
    int bcID;
} cfdMeshDataStruct;


// Container for mesh data relevant for structural/fea analysis
typedef struct {
    int coordID;
    int propertyID;

    int attrIndex;

    int constraintIndex;
    int loadIndex;
    int transferIndex;

    int connectIndex;
    int connectLinkIndex;

    int responseIndex;

    meshElementSubTypeEnum elementSubType;

} feaMeshDataStruct;


// Container for mesh data relevant for origami analysis
typedef struct {
    int propertyID;

    int constraintIndex;
    int loadIndex;
    int transferIndex;

    int neighborNodes[2];
    int foldLine;

} origamiMeshDataStruct;

#endif
