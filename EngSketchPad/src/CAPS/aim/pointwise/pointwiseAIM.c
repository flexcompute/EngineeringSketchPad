// Replace "Pointwise" and "pointwise" the "Name" and "name" of AIM respectfully!

/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Pointwise AIM
 *
 *      Modified from code written by Author, Affiliation
 *
 */


/*!\mainpage Introduction
 *
 * \section overviewPointwise Pointwise AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact with the
 * general grid generator <a href="https://www.pointwise.com/">Pointwise</a>.
 *
 * The Pointwise AIM provides the CAPS users with the ability to generate volume meshes mostly suitable for CFD analysis.
 * This includes both inviscid analysis and viscous analysis with boundary layers using the pointwise T-Rex algorithm.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsPointwise and
 * \ref aimOutputsPointwise and \ref attributePointwise, respectively.
 *
 * Files output:
 *  - caps.egads - Pointwise egads file generated
 *  - capsUserDefaults.glf - Glyph script with parameters set with \ref aimInputsPointwise
 *
 * Pointwise should be executed on Linux/macOS with the command line:
 *
 * \code{.sh}
 * pointwise -b $CAPS_GLYPH/GeomToMesh.glf caps.egads capsUserDefaults.glf
 * \endcode
 *
 * and on Windows with:
 *
 * \code{.sh}
 * %PW_HOME%\win64\bin\tclsh.exe %CAPS_GLYPH%\GeomToMesh.glf caps.egads capsUserDefaults.glf
 * \endcode
 */

 /*! \page attributePointwise AIM Attributes
 * The following list of attributes are available to guide the mesh generation with the Pointwise AIM.
 *
 * | Key                      | Value                                                                  | Geometry Location  |  Description |
 * | :------------------------| :--------------------------------------------------------------------: | :----------------: | :----------- |
 * | PW:NodeSpacing           | > 0.0                                                                  | Node               | Specified connector endpoint spacing for a node.                                                                                    |
 * |                          |                                                                        |                    |                                                                                                                                     |
 * | PW:ConnectorMaxEdge      | > 0.0                                                                  | Edge               | Maximum Edge Length in connector.                                                                                                   |
 * | PW:ConnectorEndSpacing   | > 0.0                                                                  | Edge               | Specified connector endpoint spacing.                                                                                               |
 * | PW:ConnectorDimension    | > 0                                                                    | Edge               | Specify connector dimension.                                                                                                        |
 * | PW:ConnectorAverageDS    | > 0.0                                                                  | Edge               | Specified average delta spacing for connector dimension.                                                                            |
 * | PW:ConnectorMaxAngle     | [ 0, 180 )                                                             | Edge               | Connector Maximum Angle. (0.0 = NOT APPLIED)                                                                                        |
 * | PW:ConnectorMaxDeviation | [ 0, infinity )                                                        | Edge               | Connector Maximum Deviation. (0.0 = NOT APPLIED)                                                                                    |
 * | PW:ConnectorAdaptSource  | $true or $false                                                        | Edge               | Set connector up for adaptation as a source                                                                                         |
 * |                          |                                                                        |                    |                                                                                                                                     |
 * | PW:Name                  | Set by pointwiseAIM to the value of capsGroup                          | Face               | Boundary name for domain or collection of domains.                                                                                  |
 * | PW:QuiltName             | Quilting is not supported with CAPS, but input files are generated     | Face               | Name to give one or more quilts that are assembled into a single quilt. <br>No angle test is performed.                             |
 * | PW:Baffle                | $Baffle or $Intersect                                                  | Face               | Either a true baffle surface or a surface intersected by a baffle.                                                                  |
 * | PW:DomainAlgorithm       | $Delaunay, $AdvancingFront, $AdvancingFrontOrtho                       | Face               | Surface meshing algorithm.                                                                                                          |
 * | PW:DomainIsoType         | $Triangle, $TriangleQuad                                               | Face               | Surface cell type. Global default is Triangle.                                                                                      |
 * | PW:DomainMinEdge         | > 0.0                                                                  | Face               | Cell Minimum Equilateral Edge Length in domain. (0.0 = USE BOUNDARY)                                                                |
 * | PW:DomainMaxEdge         | > 0.0                                                                  | Face               | Cell Maximum Equilateral Edge Length in domain. (0.0 = USE BOUNDARY)                                                                |
 * | PW:DomainMaxAngle        | [ 0, 180 )                                                             | Face               | Cell Maximum Angle in domain (0.0 = NOT APPLIED)                                                                                    |
 * | PW:DomainMaxDeviation    | [ 0, infinity )                                                        | Face               | Cell Maximum Deviation in domain (0.0 = NOT APPLIED)                                                                                |
 * | PW:DomainSwapCells       | $true or $false                                                        | Face               | Swap cells with no interior points.                                                                                                 |
 * | PW:DomainQuadMaxAngle    | ( 90, 180 )                                                            | Face               | Quad Maximum Included Angle in domain.                                                                                              |
 * | PW:DomainQuadMaxWarp     | ( 0, 90 )                                                              | Face               | Cell Maximum Warp Angle in domain.                                                                                                  |
 * | PW:DomainDecay           | [ 0, 1 ]                                                               | Face               | Boundary decay applied on domain.                                                                                                   |
 * | PW:DomainMaxLayers       | [ 0, infinity )                                                        | Face               | Maximum T-Rex layers in domain.                                                                                                     |
 * | PW:DomainFullLayers      | [ 0, infinity )                                                        | Face               | Number of full T-Rex layers in domain. (0 allows multi-normals)                                                                     |
 * | PW:DomainTRexGrowthRate  | [ 1, infinity )                                                        | Face               | T-Rex growth rate in domain.                                                                                                        |
 * | PW:DomainTRexType        | $Triangle, $TriangleQuad                                               | Face               | Cell types in T-Rex layers in domain.                                                                                               |
 * | PW:DomainTRexIsoHeight   | > 0.0                                                                  | Face               | Isotropic height for T-Rex cells in domain. Default is 1.0.                                                                         |
 * | PW:PeriodicTranslate     | "tx; ty; tz"                                                           | Face               | Periodic domain with given translation vector.                                                                                      |
 * | PW:PeriodicRotate        | "px; py; pz; nx; ny; nz; angle"                                        | Face               | Periodic domain with given point, normal and rotation angle.                                                                        |
 * | PW:PeriodicTarget        | $true or $false                                                        | Face               | Target domain of a translate or rotate periodic domain. <br>This domain will be deleted before the creation of the periodic domain. |
 * | PW:DomainAdaptSource     | $true or $false                                                        | Face               | Set domain up for adaptation as a source                                                                                            |
 * | PW:DomainAdaptTarget     | $true or $false                                                        | Face               | Set domain up for adaptation as a target                                                                                            |
 * | PW:DomainShapeConstraint | $DataBase or $Free                                                     | Face               | Set the domain shape constraint                                                                                                     |
 * | PW:WallSpacing           | $Wall or > 0.0                                                         | Face               | Viscous normal spacing for T-Rex extrusion. $Wall uses domParams(WallSpacing)                                                       |
 */

 /* These are set by capsUserDefaults.glf rather than as attributes the model
 * |                          |                                                                        |                    |                                                                                                                                     |
 * | PW:TRexIsoHeight         | > 0.0                                                                  | Model              | Isotropic height for volume T-Rex cells. Default is 1.0.                                                                            |
 * | PW:TRexCollisionBuffer   | > 0.0                                                                  | Model              | T-Rex collision buffer. Default is 0.5.                                                                                             |
 * | PW:TRexMaxSkewAngle      | [ 0, 180 ]                                                             | Model              | T-Rex maximum skew angle. Default 180 (Off).                                                                                        |
 * | PW:TRexGrowthRate        | [ 1, infinity )                                                        | Model              | T-Rex growth rate.                                                                                                                  |
 * | PW:TRexType              | $TetPyramid, $TetPyramidPrismHex, <br>or $AllAndConvertWallDoms        | Model              | T-Rex cell type.                                                                                                                    |
 * | PW:BoundaryDecay         | [ 0, 1 ]                                                               | Model              | Volumetric boundary decay. Default is 0.5.                                                                                          |
 * | PW:EdgeMaxGrowthRate     | [ 1, infinity )                                                        | Model              | Volumetric edge maximum growth rate. Default is 1.8.                                                                                |
 * | PW:MinEdge               | $Boundary or > 0.0                                                     | Model              | Tetrahedral Minimum Equilateral Edge Length in block.                                                                               |
 * | PW:MaxEdge               | $Boundary or > 0.0                                                     | Model              | Tetrahedral Maximum Equilateral Edge Length in block.                                                                               |
 */

// Header functions
#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "aimMesh.h"

#include "meshUtils.h"       // Collection of helper functions for meshing
#include "miscUtils.h"       // Collection of helper functions for miscellaneous analysis
#include "deprecateUtils.h"

#include "egads.h"

// Windows aliasing
#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

#define PI                 3.1415926535897931159979635
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define MAX(A,B)         (((A) < (B)) ? (B) : (A))


enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Mesh_Format,
  Mesh_ASCII_Flag,
  Mesh_Sizing,
  Mesh_Length_Factor,
  Connector_Initial_Dim,
  Connector_Max_Dim,
  Connector_Min_Dim,
  Connector_Turn_Angle,
  Connector_Deviation,
  Connector_Split_Angle,
  Connector_Turn_Angle_Hard,
  Connector_Prox_Growth_Rate,
  Connector_Adapt_Sources,
  Connector_Source_Spacing,
  Domain_Algorithm,
  Domain_Full_Layers,
  Domain_Max_Layers,
  Domain_Growth_Rate,
  Domain_Iso_Type,
  Domain_TRex_Type,
  Domain_TRex_ARLimit,
  Domain_TRex_AngleBC,
  Domain_Decay,
  Domain_Min_Edge,
  Domain_Max_Edge,
  Domain_Adapt,
  Domain_Wall_Spacing,
  Domain_Structure_AR_Convert,
  Block_Algorithm,
  Block_Voxel_Layers,
  Block_Boundary_Decay,
  Block_Collision_Buffer,
  Block_Max_Skew_Angle,
  Block_TRex_Skew_Delay,
  Block_Edge_Max_Growth_Rate,
  Block_Full_Layers,
  Block_Max_Layers,
  Block_Growth_Rate,
  Block_TRexType,
  Gen_Source_Box_Length_Scale,
  Gen_Source_Box_Direction,
  Gen_Source_Box_Angle,
  Gen_Source_Growth_Factor,
  Elevate_Cost_Threshold,
  Elevate_Max_Include_Angle,
  Elevate_Relax,
  Elevate_Smoothing_Passes,
  Elevate_WCN_Weight,
  Elevate_WCN_Mode,
  NUMINPUT = Elevate_WCN_Mode    /* Total number of inputs */
};

enum aimOutputs
{
  Volume_Mesh = 1,             /* index is 1-based */
  NUMOUT = Volume_Mesh         /* Total number of outputs */
};

static char egadsFileName[] = "caps.egads";

typedef struct {
  int    npts;
  double *xyz;
  double *t;
  int    *ivp;   // volume node index
} edgeData;

typedef struct {
  int    npts;
  double *xyz;
  double *uv;
  int    ntri;
  int    nquad;
  int    *tris;
  int    *ivp;   // volume node index
} faceData;

typedef struct {
  double **rvec;
  ego    *surfaces;
  ego    body;
  ego    *faces;
  ego    *edges;
  ego    *nodes;
  int    nfaces;
  int    nedges;
  int    nnodes;

  edgeData *tedges;
  faceData *tfaces;
} bodyData;

typedef struct {
  int ind;         // global index into UGRID file
  int egadsID;     // egadsID encoding type, body, and type-index
  double param[2]; // parametric coordinates of the vertex
} gmaVertex;

//#define DEBUG
//#define CHECKGRID // check that parametric evaluation gives the grid xyz


// Additional storage values for the instance of the AIM
typedef struct {

  // Attribute to index map
  mapAttrToIndexStruct groupMap;

  mapAttrToIndexStruct meshMap;

  aimMeshRef meshRef;

} aimStorage;


static int initiate_bodyData(int numBody, bodyData *bodydata)
{

  int i;

  for (i = 0; i < numBody; i++) {
      bodydata[i].rvec = NULL;
      bodydata[i].surfaces = NULL;
      bodydata[i].faces = NULL;
      bodydata[i].edges = NULL;
      bodydata[i].nodes = NULL;
      bodydata[i].nfaces = 0;
      bodydata[i].nedges = 0;
      bodydata[i].nnodes = 0;
      bodydata[i].tedges = NULL;
      bodydata[i].tfaces = NULL;
  }

  return CAPS_SUCCESS;
}

static int destroy_bodyData(int numBody, bodyData *bodydata)
{

  int i, j;

  if (bodydata == NULL) return CAPS_SUCCESS;

  for (i = 0; i < numBody; i++) {
    for (j = 0; j < bodydata[i].nfaces; j++) {
      if (bodydata[i].surfaces != NULL)
        if (bodydata[i].surfaces[j+bodydata[i].nfaces] != NULL)
          EG_deleteObject(bodydata[i].surfaces[j+bodydata[i].nfaces]);
      if (bodydata[i].rvec != NULL)
        EG_free(bodydata[i].rvec[j]);
    }
    EG_free(bodydata[i].nodes);
    EG_free(bodydata[i].edges);
    EG_free(bodydata[i].faces);
    EG_free(bodydata[i].surfaces);
    EG_free(bodydata[i].rvec);

    if (bodydata[i].tedges != NULL) {
        for (j = 0; j < bodydata[i].nedges; j++) {
            EG_free(bodydata[i].tedges[j].xyz);
            EG_free(bodydata[i].tedges[j].t);
            EG_free(bodydata[i].tedges[j].ivp);
        }
        EG_free(bodydata[i].tedges);
    }

    if (bodydata[i].tfaces != NULL) {
        for (j = 0; j < bodydata[i].nfaces; j++) {
            EG_free(bodydata[i].tfaces[j].xyz);
            EG_free(bodydata[i].tfaces[j].uv);
            EG_free(bodydata[i].tfaces[j].tris);
            EG_free(bodydata[i].tfaces[j].ivp);
        }
        EG_free(bodydata[i].tfaces);
    }
  }

  return CAPS_SUCCESS;
}

typedef enum {
  MODELID,  // 0
  SHELLID,  // 1
  FACEID,   // 2
  LOOPID,   // 3
  EDGEID,   // 4
  COEDGEID, // 5
  NODEID    // 6
} egadsIDEnum;


void decodeEgadsID( int id, int *type, int *bodyID, int *index )
{

    //#define PACK(t, m, i)       ((t)<<28 | (m)<<20 | (i))
    //#define UNPACK(v, t, m, i)  t = v>>28; m = (v>>20)&255; i = v&0xFFFFF;

  *type   =  id >> 28;
  *bodyID = (id >> 20) & 255;
  *index  =  id        & 0xFFFFF;

  // change to 0-based indexing
  *bodyID -= 1;
  *index  -= 1;
}


static int initiate_aimStorage(aimStorage *pointwiseInstance)
{

  int status = CAPS_SUCCESS;

  // Destroy attribute to index map
  status = initiate_mapAttrToIndexStruct(&pointwiseInstance->groupMap);
  if (status != CAPS_SUCCESS) return status;

  status = initiate_mapAttrToIndexStruct(&pointwiseInstance->meshMap);
  if (status != CAPS_SUCCESS) return status;

  // Mesh reference passed to solver
  status = aim_initMeshRef(&pointwiseInstance->meshRef);
  if (status != CAPS_SUCCESS) return status;

  return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *pointwiseInstance)
{
  int status; // Function return status

  // Destroy attribute to index map
  status = destroy_mapAttrToIndexStruct(&pointwiseInstance->groupMap);
  if (status != CAPS_SUCCESS)
    printf("Status = %d, pointwiseAIM attributeMap group cleanup!!!\n", status);

  status = destroy_mapAttrToIndexStruct(&pointwiseInstance->meshMap);
  if (status != CAPS_SUCCESS)
    printf("Status = %d, pointwiseAIM attributeMap mesh cleanup!!!\n", status);

  // Free the meshRef
  aim_freeMeshRef(&pointwiseInstance->meshRef);

  return CAPS_SUCCESS;
}


static int getUGRID(FILE *fp, meshStruct *volumeMesh)
{
  int    status = CAPS_SUCCESS;

  int    numNode, numTriangle, numQuadrilateral;
  int    numTetrahedral, numPyramid, numPrism, numHexahedral;
  int    i, elementIndex, numPoint;
  int    defaultVolID = 1; // Defailt volume ID
  double *coords = NULL;

  /* we get a binary UGRID file from Pointwise */
  status = fread(&numNode,          sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }
  status = fread(&numTriangle,      sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }
  status = fread(&numQuadrilateral, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }
  status = fread(&numTetrahedral,   sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }
  status = fread(&numPyramid,       sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }
  status = fread(&numPrism,         sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }
  status = fread(&numHexahedral,    sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; goto cleanup; }

  /*
  printf("\n Header from UGRID file: %d  %d %d  %d %d %d %d\n", numNode,
         numTriangle, numQuadrilateral, numTetrahedral, numPyramid, numPrism,
         numHexahedral);
   */

  coords = (double *) EG_alloc(3*numNode*sizeof(double));
  if (coords == NULL) return EGADS_MALLOC;

  /* read all of the vertices */
  status = fread(coords, sizeof(double), 3*numNode, fp);
  if (status != 3*numNode) { status = CAPS_IOERR; goto cleanup; }

  // TODO: Should this be something else?
  volumeMesh->analysisType = UnknownMeshAnalysis;

  // Set that this is a volume mesh
  volumeMesh->meshType = VolumeMesh;

  // Numbers
  volumeMesh->numNode = numNode;
  volumeMesh->numElement = numTetrahedral   +
                           numPyramid       +
                           numPrism         +
                           numHexahedral;

  volumeMesh->meshQuickRef.useStartIndex = (int) false;

  volumeMesh->meshQuickRef.numTriangle      = 0;
  volumeMesh->meshQuickRef.numQuadrilateral = 0;

  volumeMesh->meshQuickRef.numTetrahedral = numTetrahedral;
  volumeMesh->meshQuickRef.numPyramid     = numPyramid;
  volumeMesh->meshQuickRef.numPrism       = numPrism;
  volumeMesh->meshQuickRef.numHexahedral  = numHexahedral;

  printf("Volume mesh:\n");
  printf("\tNumber of nodes          = %d\n", numNode);
  printf("\tNumber of elements       = %d\n", volumeMesh->numElement);
  printf("\tNumber of triangles      = %d\n", numTriangle);
  printf("\tNumber of quadrilatarals = %d\n", numQuadrilateral);
  printf("\tNumber of tetrahedrals   = %d\n", numTetrahedral);
  printf("\tNumber of pyramids       = %d\n", numPyramid);
  printf("\tNumber of prisms         = %d\n", numPrism);
  printf("\tNumber of hexahedrals    = %d\n", numHexahedral);

  // Nodes - allocate
  volumeMesh->node = (meshNodeStruct *) EG_alloc(volumeMesh->numNode*
                                                 sizeof(meshNodeStruct));
  if (volumeMesh->node == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }

  // Initialize
  for (i = 0; i < volumeMesh->numNode; i++) {
      status = initiate_meshNodeStruct(&volumeMesh->node[i],
                                        volumeMesh->analysisType);
      if (status != CAPS_SUCCESS) goto cleanup;
  }

  // Nodes - set
  for (i = 0; i < volumeMesh->numNode; i++) {

    // Copy node data
    volumeMesh->node[i].nodeID = i+1;

    volumeMesh->node[i].xyz[0] = coords[3*i+0];
    volumeMesh->node[i].xyz[1] = coords[3*i+1];
    volumeMesh->node[i].xyz[2] = coords[3*i+2];
  }
  EG_free(coords); coords = NULL;

  // Elements - allocate
  volumeMesh->element = (meshElementStruct *) EG_alloc(volumeMesh->numElement*
                                                       sizeof(meshElementStruct));
  if (volumeMesh->element == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }

  // Initialize
  for (i = 0; i < volumeMesh->numElement; i++ ) {
      status = initiate_meshElementStruct(&volumeMesh->element[i],
                                           volumeMesh->analysisType);
      if (status != CAPS_SUCCESS) goto cleanup;
  }

  // Start of element index
  elementIndex = 0;

  // Skip all triangles and quads. These will be retrieved from the gma file instead
  status = fseek(fp, (mesh_numMeshConnectivity(Triangle)     *numTriangle +
                      mesh_numMeshConnectivity(Quadrilateral)*numQuadrilateral)*
                 sizeof(int), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; goto cleanup; }

#if 0
  // Elements -Set triangles
  if (numTriangle > 0) volumeMesh->meshQuickRef.startIndexTriangle = elementIndex;

  numPoint = mesh_numMeshConnectivity(Triangle);
  for (i = 0; i < numTriangle; i++) {

      volumeMesh->element[elementIndex].elementType = Triangle;
      volumeMesh->element[elementIndex].elementID   = elementIndex+1;

      status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(volumeMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; goto cleanup; }

      elementIndex += 1;
  }

  // Elements -Set quadrilateral
  if (numQuadrilateral > 0)
    volumeMesh->meshQuickRef.startIndexQuadrilateral = elementIndex;

  numPoint = mesh_numMeshConnectivity(Quadrilateral);
  for (i = 0; i < numQuadrilateral; i++) {

      volumeMesh->element[elementIndex].elementType = Quadrilateral;
      volumeMesh->element[elementIndex].elementID   = elementIndex+1;

      status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(volumeMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; goto cleanup; }

      elementIndex += 1;
  }
#endif

  // skip face ID section of the file
  // they do not map to the elements on faces
  status = fseek(fp, (numTriangle + numQuadrilateral)*sizeof(int), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; goto cleanup; }

  // Elements -Set Tetrahedral
  if (numTetrahedral > 0)
    volumeMesh->meshQuickRef.startIndexTetrahedral = elementIndex;

  numPoint = mesh_numMeshConnectivity(Tetrahedral);
  for (i = 0; i < numTetrahedral; i++) {

      volumeMesh->element[elementIndex].elementType = Tetrahedral;
      volumeMesh->element[elementIndex].elementID   = elementIndex+1;
      volumeMesh->element[elementIndex].markerID    = defaultVolID;

      status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(volumeMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; goto cleanup; }

      elementIndex += 1;
  }

  // Elements -Set Pyramid
  if (numPyramid > 0) volumeMesh->meshQuickRef.startIndexPyramid = elementIndex;

  numPoint = mesh_numMeshConnectivity(Pyramid);
  for (i = 0; i < numPyramid; i++) {

      volumeMesh->element[elementIndex].elementType = Pyramid;
      volumeMesh->element[elementIndex].elementID   = elementIndex+1;
      volumeMesh->element[elementIndex].markerID    = defaultVolID;

      status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(volumeMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; goto cleanup; }

      elementIndex += 1;
  }

  // Elements -Set Prism
  if (numPrism > 0) volumeMesh->meshQuickRef.startIndexPrism = elementIndex;

  numPoint = mesh_numMeshConnectivity(Prism);
  for (i = 0; i < numPrism; i++) {

      volumeMesh->element[elementIndex].elementType = Prism;
      volumeMesh->element[elementIndex].elementID   = elementIndex+1;
      volumeMesh->element[elementIndex].markerID    = defaultVolID;

      status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(volumeMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; goto cleanup; }

      elementIndex += 1;
  }

  // Elements -Set Hexa
  if (numHexahedral > 0)
    volumeMesh->meshQuickRef.startIndexHexahedral = elementIndex;

  numPoint = mesh_numMeshConnectivity(Hexahedral);
  for (i = 0; i < numHexahedral; i++) {

      volumeMesh->element[elementIndex].elementType = Hexahedral;
      volumeMesh->element[elementIndex].elementID   = elementIndex+1;
      volumeMesh->element[elementIndex].markerID    = defaultVolID;

      status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(volumeMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; goto cleanup; }

      elementIndex += 1;
  }

  status = CAPS_SUCCESS;

  cleanup:
      if (status != CAPS_SUCCESS)
        printf("Premature exit in getUGRID status = %d\n", status);

      EG_free(coords); coords = NULL;

      return status;
}


static int
writeMesh(void *aimInfo, aimStorage *pointwiseInstance,
          int numVolumeMesh, meshStruct *volumeMesh)
{

  int status = CAPS_SUCCESS;

  int surfIndex, volIndex;

  // analysis input values
  capsValue *ProjName;
  capsValue *MeshFormat;
  capsValue *MeshASCII_Flag;

  char *filename = NULL, surfNumber[11];
  const char *outputFileName = NULL;
  const char *outputFormat = NULL;
  int outputASCIIFlag;  // 0 = Binary output, anything else for ASCII

  bndCondStruct bndConds; // Structure of boundary conditions

  initiate_bndCondStruct(&bndConds);

  status = aim_getValue(aimInfo, Proj_Name,       ANALYSISIN, &ProjName);
  AIM_STATUS(aimInfo, status);
  if (ProjName == NULL) return CAPS_NULLVALUE;

  status = aim_getValue(aimInfo, Mesh_Format,     ANALYSISIN, &MeshFormat);
  AIM_STATUS(aimInfo, status);
  if (MeshFormat == NULL) return CAPS_NULLVALUE;

  status = aim_getValue(aimInfo, Mesh_ASCII_Flag, ANALYSISIN, &MeshASCII_Flag);
  AIM_STATUS(aimInfo, status);
  if (MeshASCII_Flag == NULL) return CAPS_NULLVALUE;

  // Project Name
  if (ProjName->nullVal == IsNull) {
    printf("No project name (\"Proj_Name\") provided - A volume mesh will not be written out\n");
    return CAPS_SUCCESS;
  }

  outputFileName  = ProjName->vals.string;
  outputFormat    = MeshFormat->vals.string;
  outputASCIIFlag = MeshASCII_Flag->vals.integer;

  if (strcasecmp(outputFormat, "SU2") != 0) {
    for (volIndex = 0; volIndex < numVolumeMesh; volIndex++) {
      volumeMesh = &volumeMesh[volIndex];
      for (surfIndex = 0; surfIndex < volumeMesh->numReferenceMesh; surfIndex++) {
/*@-bufferoverflowhigh@*/
        sprintf(surfNumber, "%d", surfIndex+1);
/*@+bufferoverflowhigh@*/
        filename = (char *) EG_alloc((strlen(outputFileName) +strlen("_Surf_") +
                                      strlen(surfNumber)+1)*sizeof(char));
        if (filename == NULL) { status = EGADS_MALLOC; goto cleanup; }
/*@-bufferoverflowhigh@*/
        sprintf(filename, "%s_Surf_%s", outputFileName, surfNumber);
/*@+bufferoverflowhigh@*/

        if (strcasecmp(outputFormat, "AFLR3") == 0) {

          status = mesh_writeAFLR3(aimInfo,
                                   filename,
                                   outputASCIIFlag,
                                   &volumeMesh->referenceMesh[surfIndex],
                                   1.0);

        } else if (strcasecmp(outputFormat, "VTK") == 0) {

          status = mesh_writeVTK(aimInfo,
                                 filename,
                                 outputASCIIFlag,
                                 &volumeMesh->referenceMesh[surfIndex],
                                 1.0);

        } else if (strcasecmp(outputFormat, "Tecplot") == 0) {

          status = mesh_writeTecplot(aimInfo,
                                     filename,
                                     outputASCIIFlag,
                                     &volumeMesh->referenceMesh[surfIndex],
                                     1.0);
        } else {
          printf("Unrecognized mesh format, \"%s\", the surface mesh will not be written out\n",
                 outputFormat);
        }

        EG_free(filename); filename = NULL;

        if (status != CAPS_SUCCESS) goto cleanup;
      }
    }
  }

  for (volIndex = 0; volIndex < numVolumeMesh; volIndex++) {

//      if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {
//          sprintf(bodyIndexNumber, "%d", bodyIndex);
//          filename = (char *) EG_alloc((strlen(outputFileName) + 1 +
//                                        strlen("_Vol")+strlen(bodyIndexNumber))*sizeof(char));
//      } else {
          filename = EG_strdup(outputFileName);
//      }

      if (filename == NULL) {
          status = EGADS_MALLOC;
          goto cleanup;
      }

      //if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {
      //    strcat(filename,"_Vol");
      //    strcat(filename,bodyIndexNumber);
      //}

      if (strcasecmp(outputFormat, "AFLR3") == 0) {

          status = mesh_writeAFLR3(aimInfo,
                                   filename,
                                   outputASCIIFlag,
                                   &volumeMesh[volIndex],
                                   1.0);

      } else if (strcasecmp(outputFormat, "VTK") == 0) {

          status = mesh_writeVTK(aimInfo,
                                 filename,
                                 outputASCIIFlag,
                                 &volumeMesh[volIndex],
                                 1.0);

      } else if (strcasecmp(outputFormat, "SU2") == 0) {

          // construct boundary condition information for SU2
          status = populate_bndCondStruct_from_mapAttrToIndexStruct(&pointwiseInstance->groupMap,
                                                                    &bndConds);
          if (status != CAPS_SUCCESS) goto cleanup;

          status = mesh_writeSU2(aimInfo,
                                 filename,
                                 outputASCIIFlag,
                                 &volumeMesh[volIndex],
                                 bndConds.numBND,
                                 bndConds.bndID,
                                 1.0);

      } else if (strcasecmp(outputFormat, "Tecplot") == 0) {

          status = mesh_writeTecplot(aimInfo,
                                     filename,
                                     outputASCIIFlag,
                                     &volumeMesh[volIndex],
                                     1.0);

      } else if (strcasecmp(outputFormat, "Nastran") == 0) {

          status = mesh_writeNASTRAN(aimInfo,
                                     filename,
                                     outputASCIIFlag,
                                     &volumeMesh[volIndex],
                                     LargeField,
                                     1.0);
      } else {
          printf("Unrecognized mesh format, \"%s\", the volume mesh will not be written out\n",
                 outputFormat);
      }

      EG_free(filename); filename = NULL;
      if (status != CAPS_SUCCESS) goto cleanup;
  }

  cleanup:
      if (status != CAPS_SUCCESS) printf("Premature exit in writeMesh status = %d\n",
                                         status);

      destroy_bndCondStruct(&bndConds);

      EG_free(filename);
      return status;
}


static int
writeGlobalGlyph(void *aimInfo, capsValue *aimInputs, double capsMeshLength)
{
    int status;

    FILE *fp = NULL;
    char filename[] = "capsUserDefaults.glf";

    // connector controls
    int         conInitDim           = aimInputs[Connector_Initial_Dim-1].vals.integer;
    int         conMaxDim            = aimInputs[Connector_Max_Dim-1].vals.integer;
    int         conMinDim            = aimInputs[Connector_Min_Dim-1].vals.integer;
    double      conTurnAngle         = aimInputs[Connector_Turn_Angle-1].vals.real;
    double      conDeviation         = aimInputs[Connector_Deviation-1].vals.real;
    double      conSplitAngle        = aimInputs[Connector_Split_Angle-1].vals.real;
    double      conProxGrowthRate    = aimInputs[Connector_Prox_Growth_Rate-1].vals.real;
    int         conAdaptSources      = aimInputs[Connector_Adapt_Sources-1].vals.integer;
    int         conSourceSpacing     = aimInputs[Connector_Source_Spacing-1].vals.integer;
    double      conTurnAngleHard     = aimInputs[Connector_Turn_Angle_Hard-1].vals.real;

    // domain controls
    const char *domAlgorithm         = aimInputs[Domain_Algorithm-1].vals.string;
    int         domFullLayers        = aimInputs[Domain_Full_Layers-1].vals.integer;
    int         domMaxLayers         = aimInputs[Domain_Max_Layers-1].vals.integer;
    double      domGrowthRate        = aimInputs[Domain_Growth_Rate-1].vals.real;
    const char *domIsoType           = aimInputs[Domain_Iso_Type-1].vals.string;
    const char *domTRexType          = aimInputs[Domain_TRex_Type-1].vals.string;
    double      domTRexARLimit       = aimInputs[Domain_TRex_ARLimit-1].vals.real;
    double      domTRexAngleBC       = aimInputs[Domain_TRex_AngleBC-1].vals.real;
    double      domDecay             = aimInputs[Domain_Decay-1].vals.real;
    double      domMinEdge           = aimInputs[Domain_Min_Edge-1].vals.real;
    double      domMaxEdge           = aimInputs[Domain_Max_Edge-1].vals.real;
    int         domAdapt             = aimInputs[Domain_Adapt-1].vals.integer;
    double      domWallSpacing       = aimInputs[Domain_Wall_Spacing-1].vals.real;
    double      domStrDomConvARTrig  = aimInputs[Domain_Structure_AR_Convert-1].vals.real;

    // block controls
    const char *blkAlgorithm         = aimInputs[Block_Algorithm-1].vals.string;
    int         blkVoxelLayers       = aimInputs[Block_Voxel_Layers-1].vals.integer;
    double      blkboundaryDecay     = aimInputs[Block_Boundary_Decay-1].vals.real;
    double      blkcollisionBuffer   = aimInputs[Block_Collision_Buffer-1].vals.real;
    double      blkmaxSkewAngle      = aimInputs[Block_Max_Skew_Angle-1].vals.real;
    int         blkTRexSkewDelay     = aimInputs[Block_TRex_Skew_Delay-1].vals.integer;
    double      blkedgeMaxGrowthRate = aimInputs[Block_Edge_Max_Growth_Rate-1].vals.real;
    int         blkfullLayers        = aimInputs[Block_Full_Layers-1].vals.integer;
    int         blkmaxLayers         = aimInputs[Block_Max_Layers-1].vals.integer;
    double      blkgrowthRate        = aimInputs[Block_Growth_Rate-1].vals.real;
    const char *blkTRexType          = aimInputs[Block_TRexType-1].vals.string;

    // general controls
    double genSourceBoxLengthScale   = aimInputs[Gen_Source_Box_Length_Scale-1].vals.real;
    double *genSourceBoxDirection    = aimInputs[Gen_Source_Box_Direction-1].vals.reals;
    double genSourceBoxAngle         = aimInputs[Gen_Source_Box_Angle-1].vals.real;
    double genSourceGrowthFactor     = aimInputs[Gen_Source_Growth_Factor-1].vals.real;

/*  These parameters are for high-order mesh generation. This should be hooked up in the future.
    // Elevate On Export controls
    double      eoecostThreshold     = aimInputs[Elevate_Cost_Threshold-1].vals.real;
    double      eoemaxIncAngle       = aimInputs[Elevate_Max_Include_Angle-1].vals.real;
    double      eoerelax             = aimInputs[Elevate_Relax-1].vals.real;
    int         eoesmoothingPasses   = aimInputs[Elevate_Smoothing_Passes-1].vals.integer;
    double      eoeWCNWeight         = aimInputs[Elevate_WCN_Weight-1].vals.real;
    const char *eoeWCNMode           = aimInputs[Elevate_WCN_Mode-1].vals.string;
*/

    // Apply the capsMeshLength scaling
    domMinEdge *= capsMeshLength;
    domMaxEdge *= capsMeshLength;
    domWallSpacing *= capsMeshLength;

    // Open the file in the analysis directory
    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        status = CAPS_IOERR;
        goto cleanup;
    }

    fprintf(fp, "# Connector level\n");
    fprintf(fp, "set conParams(InitDim)                          %d; # Initial connector dimension\n", conInitDim);
    fprintf(fp, "set conParams(MaxDim)                           %d; # Maximum connector dimension\n", conMaxDim);
    fprintf(fp, "set conParams(MinDim)                           %d; # Minimum connector dimension\n", conMinDim);
    fprintf(fp, "set conParams(TurnAngle)                       %lf; # Maximum turning angle on connectors for dimensioning (0 - not used)\n", conTurnAngle);
    fprintf(fp, "set conParams(Deviation)                       %lf; # Maximum deviation on connectors for dimensioning (0 - not used)\n", conDeviation);
    fprintf(fp, "set conParams(SplitAngle)                      %lf; # Turning angle on connectors to split (0 - not used)\n", conSplitAngle);
    fprintf(fp, "set conParams(JoinCons)                          0; # Perform joining operation on 2 connectors at one endpoint\n"); // This must be 0 in order to put the surface mesh back on an egads tessellation
    fprintf(fp, "set conParams(ProxGrowthRate)                  %lf; # Connector proximity growth rate\n", conProxGrowthRate);
    fprintf(fp, "set conParams(AdaptSources)                     %d; # Compute sources using connectors (0 - not used) V18.2+ (experimental)\n", conAdaptSources);
    fprintf(fp, "set conParams(SourceSpacing)                    %d; # Use source cloud for adaptive pass on connectors V18.2+\n", conSourceSpacing);
    fprintf(fp, "set conParams(TurnAngleHard)                   %lf; # Hard edge turning angle limit for domain T-Rex (0.0 - not used)\n", conTurnAngleHard);
    fprintf(fp, "\n");
    fprintf(fp, "# Domain level\n");
    fprintf(fp, "set domParams(Algorithm)                    \"%s\"; # Isotropic (Delaunay, AdvancingFront or AdvancingFrontOrtho)\n", domAlgorithm);
    fprintf(fp, "set domParams(FullLayers)                       %d; # Domain full layers (0 for multi-normals, >= 1 for single normal)\n", domFullLayers);
    fprintf(fp, "set domParams(MaxLayers)                        %d; # Domain maximum layers\n", domMaxLayers);
    fprintf(fp, "set domParams(GrowthRate)                      %lf; # Domain growth rate for 2D T-Rex extrusion\n", domGrowthRate);
    fprintf(fp, "set domParams(IsoType)                      \"%s\"; # Domain iso cell type (Triangle or TriangleQuad)\n", domIsoType);
    fprintf(fp, "set domParams(TRexType)                     \"%s\"; # Domain T-Rex cell type (Triangle or TriangleQuad)\n", domTRexType);
    fprintf(fp, "set domParams(TRexARLimit)                     %lf; # Domain T-Rex maximum aspect ratio limit (0 - not used)\n", domTRexARLimit);
    fprintf(fp, "set domParams(TRexAngleBC)                     %lf; # Domain T-Rex spacing from surface curvature\n", domTRexAngleBC);
    fprintf(fp, "set domParams(Decay)                           %lf; # Domain boundary decay\n", domDecay);
    fprintf(fp, "set domParams(MinEdge)                         %lf; # Domain minimum edge length\n", domMinEdge);
    fprintf(fp, "set domParams(MaxEdge)                         %lf; # Domain maximum edge length\n", domMaxEdge);
    fprintf(fp, "set domParams(Adapt)                            %d; # Set up all domains for adaptation (0 - not used) V18.2+ (experimental)\n", domAdapt);
    fprintf(fp, "set domParams(WallSpacing)                     %lf; # defined spacing when geometry attributed with $wall\n", domWallSpacing);
    fprintf(fp, "set domParams(StrDomConvertARTrigger)          %lf; # Aspect ratio to trigger converting domains to structured\n", domStrDomConvARTrig);
    fprintf(fp, "\n");
    fprintf(fp, "# Block level\n");
    fprintf(fp, "set blkParams(Algorithm)                    \"%s\"; # Isotropic (Delaunay, Voxel) (V18.3+)\n", blkAlgorithm);
    fprintf(fp, "set blkParams(VoxelLayers)                      %d; # Number of Voxel transition layers if Algorithm set to Voxel (V18.3+)\n", blkVoxelLayers);
    fprintf(fp, "set blkParams(boundaryDecay)                   %lf; # Volumetric boundary decay\n", blkboundaryDecay);
    fprintf(fp, "set blkParams(collisionBuffer)                 %lf; # Collision buffer for colliding T-Rex fronts\n", blkcollisionBuffer);
    fprintf(fp, "set blkParams(maxSkewAngle)                    %lf; # Maximum skew angle for T-Rex extrusion\n", blkmaxSkewAngle);
    fprintf(fp, "set blkParams(TRexSkewDelay)                    %d; # Number of layers to delay enforcement of skew criteria\n", blkTRexSkewDelay);
    fprintf(fp, "set blkParams(edgeMaxGrowthRate)               %lf; # Volumetric edge ratio\n", blkedgeMaxGrowthRate);
    fprintf(fp, "set blkParams(fullLayers)                       %d; # Full layers (0 for multi-normals, >= 1 for single normal)\n", blkfullLayers);
    fprintf(fp, "set blkParams(maxLayers)                        %d; # Maximum layers\n", blkmaxLayers);
    fprintf(fp, "set blkParams(growthRate)                      %lf; # Growth rate for volume T-Rex extrusion\n", blkgrowthRate);
    fprintf(fp, "set blkParams(TRexType)                     \"%s\"; # T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms)\n", blkTRexType);
    fprintf(fp, "set blkParams(volInitialize)                     1; # Initialize block after setup\n");
    fprintf(fp, "\n");
    fprintf(fp, "# General\n");
    fprintf(fp, "set genParams(SkipMeshing)                       1; # Skip meshing of domains during interim processing (V18.3+)\n");
    fprintf(fp, "set genParams(CAESolver)                 \"UGRID\"; # Selected CAE Solver (Currently support CGNS, Gmsh and UGRID)\n");
    fprintf(fp, "set genParams(outerBoxScale)                     0; # Enclose geometry in box with specified scale (0 - no box)\n");
    fprintf(fp, "set genParams(sourceBoxLengthScale)            %lf; # Length scale of enclosed viscous walls in source box (0 - no box)\n", genSourceBoxLengthScale);
    fprintf(fp, "set genParams(sourceBoxDirection)  { %lf %lf %lf }; # Principal direction vector (i.e. normalized freestream vector)\n", genSourceBoxDirection[0], genSourceBoxDirection[1], genSourceBoxDirection[2]);
    fprintf(fp, "set genParams(sourceBoxAngle)                  %lf; # Angle for widening source box in the assigned direction\n", genSourceBoxAngle);
    fprintf(fp, "set genParams(sourceGrowthFactor)              %lf; # Growth rate for spacing value along box\n", genSourceGrowthFactor);
    fprintf(fp, "set genParams(ModelSize)                         0; # Set model size before CAD import (0 - get from file)\n");
    fprintf(fp, "set genParams(writeGMA)                    \"2.0\"; # Write out geometry-mesh associativity file version (0.0 - none, 1.0 or 2.0)\n");
    fprintf(fp, "set genParams(assembleTolMult)                 2.0; # Multiplier on model assembly tolerance for allowed MinEdge\n");
    fprintf(fp, "set genParams(modelOrientIntoMeshVolume)         1; # Whether the model is oriented so normals point into the mesh\n");
    fprintf(fp, "\n");
/*  These parameters are for high-order mesh generation. This should be hooked up in the future.
    fprintf(fp, "# Elevate On Export V18.2+\n");
    fprintf(fp, "set eoeParams(degree)                           Q1; # Polynomial degree (Q1, Q2, Q3 or Q4) NOTE: ONLY APPLIES TO CGNS AND GMSH\n");
    fprintf(fp, "set eoeParams(costThreshold)                   %lf; # Cost convergence threshold\n", eoecostThreshold);
    fprintf(fp, "set eoeParams(maxIncAngle)                     %lf; # Maximum included angle tolerance\n", eoemaxIncAngle);
    fprintf(fp, "set eoeParams(relax)                           %lf; # Iteration relaxation factor\n", eoerelax);
    fprintf(fp, "set eoeParams(smoothingPasses)                  %d; # Number of smoothing passes\n", eoesmoothingPasses);
    fprintf(fp, "set eoeParams(WCNWeight)                       %lf; # WCN cost component weighting factor\n", eoeWCNWeight);
    fprintf(fp, "set eoeParams(WCNMode)                      \"%s\"; # WCN weight factor method (UseValue or Calculate)\n", eoeWCNMode);
    fprintf(fp, "set eoeParams(writeVTU)                  \"false\"; # Write out ParaView VTU files (true or false)\n");
*/

    status = CAPS_SUCCESS;

cleanup:
    if (fp != NULL) fclose(fp);

    return status;
}


#if 0
static int set_autoQuiltAttr(ego *body) {

    int status; // Function return status

    int i, j; // Indexing

    int numBodyFace, numChildrenLoop;
    ego *bodyFace = NULL, *childrenLoop=NULL;

    ego geomRef, geomRef2;
    int oclass, mtypeFace, *senses;
    double data[18];

    char pwQuiltAttr[] = "PW:QuiltName";

    const char *face1Attr = NULL, *face2Attr=NULL;

    char quiltNamePrefix[] = "capsAutoQuilt_";
    char quiltName[30];
    int quiltIndex = 1;

    status = EG_getBodyTopos(*body, NULL, FACE, &numBodyFace, &bodyFace);
    if (status != EGADS_SUCCESS) goto cleanup;

    // Loop through the faces and find what edges
    for (i = 0; i < numBodyFace; i++) {

        status = EG_getTopology(bodyFace[i], &geomRef, &oclass, &mtypeFace, data,
                                &numChildrenLoop, &childrenLoop, &senses);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = retrieve_stringAttr(bodyFace[i], pwQuiltAttr, &face1Attr);
        if (status == EGADS_NOTFOUND) {
            face1Attr = NULL;
        } else if (status != CAPS_SUCCESS) {
            goto cleanup;
        }

        for (j = 0; j < numBodyFace; j++) {

            // If faces match
            if (i == j) continue;

            status = EG_getTopology(bodyFace[j], &geomRef2, &oclass, &mtypeFace, data,
                                    &numChildrenLoop, &childrenLoop, &senses);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_isSame(geomRef, geomRef2);
            if (status < CAPS_SUCCESS ) goto cleanup;

            if (status == EGADS_OUTSIDE) continue;

            // I have the same piece of geometry - should be continuous - let's quilt
            status = retrieve_stringAttr(bodyFace[j], pwQuiltAttr, &face2Attr);
            if (status == EGADS_NOTFOUND) {
                face2Attr = NULL;
            } else if (status != CAPS_SUCCESS) {
                goto cleanup;
            }

            // If neither face has an attribute
            if (face1Attr == NULL && face2Attr == NULL) {

                sprintf(quiltName, "%s%d", quiltNamePrefix, quiltIndex);

                // Face 1
                status = EG_attributeAdd(bodyFace[i], pwQuiltAttr, ATTRSTRING,
                                         0, NULL, NULL, quiltName);
                if (status != EGADS_SUCCESS) goto cleanup;

                // Face 2
                status = EG_attributeAdd(bodyFace[j], pwQuiltAttr, ATTRSTRING,
                                         0, NULL, NULL, quiltName);
                if (status != EGADS_SUCCESS) goto cleanup;

                quiltIndex += 1;

            } else if (face1Attr == NULL) {

                status = EG_attributeAdd(bodyFace[i], pwQuiltAttr, ATTRSTRING,
                                         0, NULL, NULL, face2Attr);
                if (status != EGADS_SUCCESS) goto cleanup;

            } else if (face2Attr == NULL) {

                status = EG_attributeAdd(bodyFace[i], pwQuiltAttr, ATTRSTRING,
                                         0, NULL, NULL, face1Attr);
                if (status != EGADS_SUCCESS) goto cleanup;
            } else {
                // Both Face 1 and Face 2 already have attributes
                continue;
            }
        } // Second face loop
    } // Face loop

    status = CAPS_SUCCESS;

    cleanup:
        if (status != CAPS_SUCCESS)
          printf("Error: Premature exit in set_autoQuiltAttr, status %d\n",
                 status);

        EG_free(bodyFace);

        return status;
}
#endif


static int setPWAttr(void *aimInfo,
                     ego body,
        /*@unused@*/ mapAttrToIndexStruct *groupMap,
                     mapAttrToIndexStruct *meshMap,
                     int numMeshProp,
                     meshSizingStruct *meshProp,
                     double capsMeshLength,
                     int *quilting)
{

    int status; // Function return status

    int nodeIndex, edgeIndex, faceIndex; // Indexing

    int numNode;
    ego *nodes = NULL;

    int numEdge;
    ego *edges = NULL;

    int numFace;
    ego *faces = NULL;

    const char *meshName = NULL;
    const char *groupName = NULL;
    int attrIndex = 0;
    int propIndex = 0;
    double real = 0;

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char   *string;

    status = EG_getBodyTopos(body, NULL, NODE, &numNode, &nodes);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (nodes == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    // Loop through the nodes and set PW:NodeSpacing attribute
    for (nodeIndex = 0; nodeIndex < numNode; nodeIndex++) {

        status = retrieve_CAPSMeshAttr(nodes[nodeIndex], &meshName);
        if (status != EGADS_SUCCESS) continue;

/*@-nullpass@*/
        status = get_mapAttrToIndexIndex(meshMap, meshName, &attrIndex);
        if (status != CAPS_SUCCESS) {
          printf("Error: Unable to retrieve index from capsMesh %s\n", meshName);
          goto cleanup;
        }
/*@+nullpass@*/

        for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

            // Check if the mesh property applies to the capsGroup of this edge
            if (meshProp[propIndex].attrIndex != attrIndex) continue;

            // Is the attribute set?
            if (meshProp[propIndex].nodeSpacing > 0) {

                real = capsMeshLength*meshProp[propIndex].nodeSpacing;

                // add the attribute
                status = EG_attributeAdd(nodes[nodeIndex], "PW:NodeSpacing",
                                         ATTRREAL, 1, NULL, &real, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }
        }
    }

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (edges == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    // Loop through the edges and set PW:Connector* attributes
    for (edgeIndex = 0; edgeIndex < numEdge; edgeIndex++) {

        status = retrieve_CAPSMeshAttr(edges[edgeIndex], &meshName);
        if (status != EGADS_SUCCESS) continue;

/*@-nullpass@*/
        status = get_mapAttrToIndexIndex(meshMap, meshName, &attrIndex);
        if (status != CAPS_SUCCESS) {
          printf("Error: Unable to retrieve index from capsMesh %s\n", meshName);
          goto cleanup;
        }
/*@+nullpass@*/

        for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

            // Check if the mesh property applies to the capsGroup of this edge
            if (meshProp[propIndex].attrIndex != attrIndex) continue;

            // Is the attribute set?
            if (meshProp[propIndex].maxSpacing > 0) {

                real = capsMeshLength*meshProp[propIndex].maxSpacing;

                // add the attribute
                status = EG_attributeAdd(edges[edgeIndex],
                                         "PW:ConnectorMaxEdge", ATTRREAL, 1,
                                         NULL, &real, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }

            // Is the attribute set?
            if (meshProp[propIndex].nodeSpacing > 0) {

                real = capsMeshLength*meshProp[propIndex].nodeSpacing;

                // add the attribute
                status = EG_attributeAdd(edges[edgeIndex],
                                         "PW:ConnectorEndSpacing", ATTRREAL, 1,
                                         NULL, &real, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }

            // Is the attribute set?
            if (meshProp[propIndex].numEdgePoints > 0) {

                // add the attribute
                status = EG_attributeAdd(edges[edgeIndex],
                                         "PW:ConnectorDimension", ATTRINT, 1,
                                         &meshProp[propIndex].numEdgePoints,
                                         NULL, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }

            // Is the attribute set?
            if (meshProp[propIndex].avgSpacing > 0) {

                real = capsMeshLength*meshProp[propIndex].avgSpacing;

                // add the attribute
                status = EG_attributeAdd(edges[edgeIndex],
                                         "PW:ConnectorAverageDS", ATTRREAL, 1,
                                         NULL, &real, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }

            // Is the attribute set?
            if (meshProp[propIndex].maxAngle > 0) {

                real = meshProp[propIndex].maxAngle;

                // add the attribute
                status = EG_attributeAdd(edges[edgeIndex], "PW:ConnectorMaxAngle",
                                         ATTRREAL, 1, NULL, &real, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }

            // Is the attribute set?
            if (meshProp[propIndex].maxDeviation > 0) {

                real = capsMeshLength*meshProp[propIndex].maxDeviation;

                // add the attribute
                status = EG_attributeAdd(edges[edgeIndex],
                                         "PW:ConnectorMaxDeviation", ATTRREAL, 1,
                                         NULL, &real, NULL);
                if (status != EGADS_SUCCESS) goto cleanup;
            }
        }
    }

    status = EG_getBodyTopos(body, NULL, FACE, &numFace, &faces);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (faces == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    // Loop through the faces and copy capsGroup to PW:Name and set PW:Domain* attributes
    for (faceIndex = 0; faceIndex < numFace; faceIndex++) {

        status = retrieve_CAPSGroupAttr(faces[faceIndex], &groupName);
        if (status == EGADS_SUCCESS) {

            status = EG_attributeAdd(faces[faceIndex], "PW:Name",
                                     ATTRSTRING, 0, NULL, NULL, groupName);
            if (status != EGADS_SUCCESS) goto cleanup;

        } else {
            AIM_ERROR(aimInfo, "No capsGroup attribute found on Face %d",
                   faceIndex+1);
            print_AllAttr( aimInfo, faces[faceIndex] );
            goto cleanup;
        }

        status = retrieve_CAPSMeshAttr(faces[faceIndex], &meshName);
        if (status == EGADS_SUCCESS) {

/*@-nullpass@*/
            status = get_mapAttrToIndexIndex(meshMap, meshName, &attrIndex);
            if (status != CAPS_SUCCESS) {
                printf("Error: Unable to retrieve index from capsGroup %s\n",
                       meshName);
                goto cleanup;
            }
/*@+nullpass@*/

            for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

                // Check if the mesh property applies to the capsGroup of this face
                if (meshProp[propIndex].attrIndex != attrIndex) continue;

                // Is the attribute set?
                if (meshProp[propIndex].minSpacing > 0) {

                    real = capsMeshLength*meshProp[propIndex].minSpacing;

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainMinEdge",
                                             ATTRREAL, 1, NULL, &real, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].maxSpacing > 0) {

                    real = capsMeshLength*meshProp[propIndex].maxSpacing;

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainMaxEdge",
                                             ATTRREAL, 1, NULL, &real, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].maxAngle > 0) {

                    real = meshProp[propIndex].maxAngle;

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainMaxAngle",
                                             ATTRREAL, 1, NULL, &real, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].maxDeviation > 0) {

                    real = capsMeshLength*meshProp[propIndex].maxDeviation;

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainMaxDeviation",
                                             ATTRREAL, 1, NULL, &real, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].boundaryDecay > 0) {

                    real = meshProp[propIndex].boundaryDecay;

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainDecay",
                                             ATTRREAL, 1, NULL, &real, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].boundaryLayerMaxLayers > 0) {

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainMaxLayers",
                                             ATTRINT, 1,
                                             &meshProp[propIndex].boundaryLayerMaxLayers,
                                             NULL, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].boundaryLayerFullLayers > 0) {

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:DomainFullLayers",
                                             ATTRINT, 1,
                                             &meshProp[propIndex].boundaryLayerFullLayers,
                                             NULL, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].boundaryLayerGrowthRate > 0) {

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex],
                                             "PW:DomainTRexGrowthRate",
                                             ATTRREAL, 1, NULL,
                                             &meshProp[propIndex].boundaryLayerGrowthRate,
                                             NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }

                // Is the attribute set?
                if (meshProp[propIndex].boundaryLayerSpacing > 0) {

                    real = capsMeshLength*meshProp[propIndex].boundaryLayerSpacing;

                    // add the attribute
                    status = EG_attributeAdd(faces[faceIndex], "PW:WallSpacing",
                                             ATTRREAL, 1, NULL, &real, NULL);
                    if (status != EGADS_SUCCESS) goto cleanup;
                }
            }
        }

        // Check for quilting on faces
        status = EG_attributeRet(faces[faceIndex], "PW:QuiltName", &atype, &alen,
                                 &ints, &reals, &string);
        if (status == EGADS_SUCCESS) *quilting = (int)true;

    } // Face loop

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in setPWAttr, status %d\n", status);

    EG_free(nodes);
    EG_free(edges);
    EG_free(faces);

    return status;
}


// corrects uv values from pointwise gma files
static int correctUV(ego face, ego surf, ego eBsplSurf, /*@unused@*/ double *rvec,
                     double *offuv)
{
  int    status;
  int    i, oclass, mtype, fper, sper;
  double offU, offV, frange[4], srange[4], perU, perV;
  ego    ref;

  /* get the surface's representation */
  status = EG_getGeometry(surf, &oclass, &mtype, &ref, NULL, NULL);
  if (status != EGADS_SUCCESS) {
      printf(" Error: EG_getGeometry = %d (correctUV)!\n", status);
      return status;
  }
  offuv[0] = offuv[1] = 0;
  if (mtype == BSPLINE) return EGADS_SUCCESS;
  offU = offV = 0.0;

  status = EG_getRange(surf, frange, &fper);
  if (status != EGADS_SUCCESS) {
      printf(" Error: EG_getRange0 = %d (correctUV)!\n", status);
      return status;
  }
  perU = perV = 0.0;
  if ((fper & 1) != 0) perU = frange[1] - frange[0];
  if ((fper & 2) != 0) perV = frange[3] - frange[2];
  status = EG_getRange(face, frange, &i);
  if (status != EGADS_SUCCESS) {
      printf(" Error: EG_getRange1 = %d (correctUV)!\n", status);
      return status;
  }
  status = EG_getRange(eBsplSurf, srange, &sper);
  if (status != EGADS_SUCCESS) {
      printf(" Error: EG_getRange0 = %d (correctUV)!\n", status);
      return status;
  }
  if ((fabs(frange[0] - srange[0]) > 1.e-4) ||
      (fabs(frange[1] - srange[1]) > 1.e-4)) {
      //printf(" Internal: Uspan = %lf %lf   %lf %lf -- %lf  %d\n",
      //    frange[0], frange[1], srange[0], srange[1], perU, fper);
      offU = perU;
      if ((frange[0] - srange[0]) < 0.0) offU = -perU;
      //printf("    uOffs = %lf %lf   %lf\n",
      //    srange[0] + offU, srange[1] + offU, offU);
  }
  if ((fabs(frange[2] - srange[2]) > 1.e-4) ||
      (fabs(frange[3] - srange[3]) > 1.e-4)) {
      //printf(" Internal: Vspan = %lf %lf   %lf %lf -- %lf  %d\n",
      //    frange[2], frange[3], srange[2], srange[3], perV, fper);
      offV = perV;
      if ((frange[2] - srange[2]) < 0.0) offV = -perV;
      //printf("    vOffs = %lf %lf   %lf\n",
      //    srange[2] + offV, srange[3] + offV, offV);
  }

  offuv[0] = offU;
  offuv[1] = offV;

  return EGADS_SUCCESS;

#if 0
  status = EG_getRange(face, limits, &iper);
  if (status != EGADS_SUCCESS) {
    printf(" getRange = %d\n", status);
    return;
  }
  status = EG_evaluate(newSurf, uv, results);
  if (status != EGADS_SUCCESS) {
    printf(" Error: Surface eval status = %d!\n", status);
    return;
  }

  if ((surf->mtype == CYLINDRICAL) || (surf->mtype == CONICAL)) {
    x = (results[0]-rvec[0])*rvec[3] + (results[1]-rvec[1])*rvec[4] +
        (results[2]-rvec[2])*rvec[5];
    y = (results[0]-rvec[0])*rvec[6] + (results[1]-rvec[1])*rvec[7] +
        (results[2]-rvec[2])*rvec[8];
    d = atan2(y, x);
    while (d < limits[0]) d += 2*PI;
    while (d > limits[1]) d -= 2*PI;
    uv[0] = d;
  } else if (surf->mtype == SPHERICAL) {
    d  =  rvec[9];
    x1 = &rvec[3];
    x2 = &rvec[6];
    CROSS(norm, x1, x2);
    di = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (d < 0.0) {
      di = -di;
      d  = -d;
    }
    if (di != 0.0) {
      norm[0] /= di;
      norm[1] /= di;
      norm[2] /= di;
    }
    x = (results[0]-rvec[0])*rvec[3] + (results[1]-rvec[1])*rvec[4] +
        (results[2]-rvec[2])*rvec[5];
    y = (results[0]-rvec[0])*rvec[6] + (results[1]-rvec[1])*rvec[7] +
        (results[2]-rvec[2])*rvec[8];
    z = (results[0]-rvec[0])*norm[0] + (results[1]-rvec[1])*norm[1] +
        (results[2]-rvec[2])*norm[2];
    uv[1] = asin(z/d);
    while (uv[1] < limits[2]) uv[1] += PI;
    while (uv[1] > limits[3]) uv[1] -= PI;
    uv[0] = atan2(y, x);
    while (uv[0] < limits[0]) uv[0] += 2*PI;
    while (uv[0] > limits[1]) uv[0] -= 2*PI;
  } else {
    x = (results[0]-rvec[0])*rvec[3] + (results[1]-rvec[1])*rvec[4] +
        (results[2]-rvec[2])*rvec[5];
    y = (results[0]-rvec[0])*rvec[6] + (results[1]-rvec[1])*rvec[7] +
        (results[2]-rvec[2])*rvec[8];
    z = (results[0]-rvec[0])*rvec[9] + (results[1]-rvec[1])*rvec[10] +
        (results[2]-rvec[2])*rvec[11];
    uv[1] = asin(z/rvec[13]);
    while (uv[1] < limits[2]) uv[1] += 2*PI;
    while (uv[1] > limits[3]) uv[1] -= 2*PI;
    uv[0] = atan2(y, x);
    while (uv[0] < limits[0]) uv[0] += 2*PI;
    while (uv[0] > limits[1]) uv[0] -= 2*PI;
  }
#endif
}


typedef struct {
  int ibody;    // 0-based body index
  int iedge;    // 0-based edge index
} edgeMapData;


static int matchSameEdges(ego *bodies, int numBody, edgeMapData ***edgeMapOut)
{
  int status = CAPS_SUCCESS;
  int i, j, isheet, isolid;
  int numSolid=0, numSheet=0;
  int oclass, mtype, nchild, *senses;
  int nSheetEdge, isheetedge, nSolidEdge, isolidedge, nSheetNode, nSolidNode;
  int mtypeSheet, mtypeSolid, numEdge;
  int *solidBodies=NULL, *sheetBodies=NULL;
  double data[4];
  ego geom, *children, *edges=NULL, *sheetEdges=NULL, *solidEdges=NULL;
  ego *sheetNodes, *solidNodes;
  edgeMapData **edgeMap = NULL;

  solidBodies = (int*)EG_alloc(numBody*sizeof(int));
  if (solidBodies == NULL) { status = EGADS_MALLOC; goto cleanup; }

  sheetBodies = (int*)EG_alloc(numBody*sizeof(int));
  if (sheetBodies == NULL) { status = EGADS_MALLOC; goto cleanup; }

  edgeMap = (edgeMapData**)EG_alloc(numBody*sizeof(edgeMapData*));
  if (edgeMap == NULL) { status = EGADS_MALLOC; goto cleanup; }
  for (i = 0; i < numBody; i++) {
    edgeMap[i] = NULL;
  }

  for (i = 0; i < numBody; i++) {
    status = EG_getTopology(bodies[i], &geom, &oclass, &mtype, data,
                            &nchild, &children, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
         if (mtype == SOLIDBODY) solidBodies[numSolid++] = i;
    else if (mtype == SHEETBODY) sheetBodies[numSheet++] = i;
    else {
      printf(" Error: Unsupported body type!\n");
      status = CAPS_BADTYPE;
      goto cleanup;
    }

    // get the Edges to fill the default identity edgeMap
    status = EG_getBodyTopos(bodies[i], NULL, EDGE, &numEdge, &edges);
    if (status != EGADS_SUCCESS) goto cleanup;

    edgeMap[i] = (edgeMapData*)EG_alloc(numEdge*sizeof(edgeMapData));
    if (edgeMap[i] == NULL) { status = EGADS_MALLOC; goto cleanup; }

    for (j = 0; j < numEdge; j++) {
      edgeMap[i][j].ibody = i;
      edgeMap[i][j].iedge = j;
    }
    EG_free(edges); edges = NULL;
  }

  //Find all sheet body edges that also exist in solid bodies
  for ( isheet = 0; isheet < numSheet; isheet++ ) {

    status = EG_getBodyTopos(bodies[sheetBodies[isheet]], NULL, EDGE,
                             &nSheetEdge, &sheetEdges);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (sheetEdges == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    //Loop over the solid bodies for each sheet body
    for ( isolid = 0; isolid < numSolid; isolid++ ) {

      status = EG_getBodyTopos(bodies[solidBodies[isolid]], NULL, EDGE,
                               &nSolidEdge, &solidEdges);
      if (status != EGADS_SUCCESS) goto cleanup;
      if (solidEdges == NULL) {
          status = CAPS_NULLOBJ;
          goto cleanup;
      }

      //Loop over all edges in the sheet body and solid body
      for ( isheetedge = 0; isheetedge < nSheetEdge; isheetedge++ ) {
        for ( isolidedge = 0; isolidedge < nSolidEdge; isolidedge++ ) {

          // Check if the sheet edge and solid edge are on the same geometry
          if (EG_isSame(sheetEdges[isheetedge], solidEdges[isolidedge]) ==
              EGADS_SUCCESS) {

            // Get the nodes on the two edges
            status = EG_getTopology(solidEdges[isolidedge], &geom, &oclass,
                                    &mtypeSolid, data, &nSolidNode, &solidNodes,
                                    &senses);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_getTopology(sheetEdges[isheetedge], &geom, &oclass,
                                    &mtypeSheet, data, &nSheetNode, &sheetNodes,
                                    &senses);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Check that the nodes match on the edges
            if (mtypeSolid == ONENODE) {

              if ( !(EG_isSame(solidNodes[0], sheetNodes[0]) == EGADS_SUCCESS) ) continue;

            } else {
              if ( !( (EG_isSame(solidNodes[0], sheetNodes[0])==EGADS_SUCCESS &&
                       EG_isSame(solidNodes[1], sheetNodes[1])==EGADS_SUCCESS) ||
                  (EG_isSame(solidNodes[0], sheetNodes[1])==EGADS_SUCCESS &&
                   EG_isSame(solidNodes[1], sheetNodes[0])==EGADS_SUCCESS) ) ) continue;
            }

            //Update the index map to register the equality
            edgeMap[sheetBodies[isheet]][isheetedge].ibody = solidBodies[isolid];
            edgeMap[sheetBodies[isheet]][isheetedge].iedge = isolidedge;

            edgeMap[solidBodies[isolid]][isolidedge].ibody = sheetBodies[isheet];
            edgeMap[solidBodies[isolid]][isolidedge].iedge = isheetedge;
            break;
          }
        }
      }
      EG_free(solidEdges); solidEdges = NULL;
    }
    EG_free(sheetEdges); sheetEdges = NULL;
  }

  (*edgeMapOut) = edgeMap;

cleanup:
  if (status != EGADS_SUCCESS) {
    printf("Error: Premature exit in matchSameEdges status = %d\n", status);
    if (edgeMap != NULL) {
      for (i = 0; i < numBody; i++)
        EG_free(edgeMap[i]);
/*@-dependenttrans@*/
      EG_free(edgeMap);
/*@+dependenttrans@*/
    }
  }

  EG_free(solidBodies);
  EG_free(sheetBodies);
  EG_free(edges);

  EG_free(sheetEdges);
  EG_free(solidEdges);
  return status;
}


#if 0
// count faces that must occur twice in the tessellation of a face
// these are on the bounds of a periodic uv-space
static int getFaceEdgeCount(ego body, ego face, int *nedge, ego **edges,
                            int **edgeCount)
{

  int status = CAPS_SUCCESS;

  int *senses, *count = NULL;
  int oclass, mtype, iedge, nnode;
  double trange[4], uv[2];
  ego *nodes, ref;

  *edgeCount = NULL;

  status = EG_getBodyTopos(body, face, EDGE, nedge, edges);
  if (status != EGADS_SUCCESS) {
      printf("EG_getBodyTopos EDGE = %d\n", status);
      goto cleanup;
  }

  count = (int*)EG_alloc((*nedge)*sizeof(int));
  if (count == NULL) { status = EGADS_MALLOC; goto cleanup; }

  for (iedge = 0; iedge < *nedge; iedge++) {
    count[iedge] = 1;

    status = EG_getTopology((*edges)[iedge], &ref, &oclass, &mtype, trange,
                            &nnode, &nodes, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    // try to get edge UV at t-mid with sense of 0
    status = EG_getEdgeUV(face, (*edges)[iedge], 0, (trange[0]+trange[1])/2, uv);
    if (status == EGADS_TOPOERR) {
      /* sense in Face twice! */
      count[iedge] = 2;
    }
  }

  *edgeCount = count;

  status = CAPS_SUCCESS;
  cleanup:
    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in getFaceEdgeCount status = %d\n", status);

    if (status != CAPS_SUCCESS) EG_free(count);

    return status;
}


// gets a face coordinates and UV values from FACE/EDGE/NODE points
static int getFacePoints(bodyData *bodydata, int ibody, int iface,
                         meshStruct *volumeMesh, int nSurfPts,
                         gmaVertex *surfacedata, int *face_pnt, int *face_ind,
                         double *uv)
{
  int status = CAPS_SUCCESS;
  int it, ib, in, iper, edgeIndex;
  int ifp, isp, ivp; // face point, surface point, volume point
  int oclass, mtype, npts, *senses, *tmp, offset;
  int nloop, iloop, nedge, iedge, nnode, i;
  double limits[4], uvbox[4], trange[4], t, *exyz; //, result[18];
  ego face, geom, *loops, *edges, *nodes, ref;

  face = bodydata->faces[iface];
  status = EG_getTopology(face, &ref, &oclass, &mtype, uvbox, &nloop, &loops,
                          &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  geom = bodydata->surfaces[iface];
  if (geom->mtype == BSPLINE) {
    status = EG_getRange(geom, limits, &iper);
    if (status != EGADS_SUCCESS) goto cleanup;
  } else {
    for (i = 0; i < 4; i++)
      limits[i] = uvbox[i];
  }

  // populate UV values interior to the face
  for (isp = 0; isp < nSurfPts; isp++) {
    if (face_pnt[isp] == 1) {

      decodeEgadsID(surfacedata[isp].egadsID, &it, &ib, &in);
      if (it != FACEID) continue;

      if (ibody != ib) {
        printf(" Error: Inconsistent body index!\n" );
        goto cleanup;
      }

      ivp = surfacedata[isp].ind-1;
      ifp = face_ind[isp];

      uv[0] = limits[0] + surfacedata[isp].param[0]*(limits[1]-limits[0]);
      uv[1] = limits[2] + surfacedata[isp].param[1]*(limits[3]-limits[2]);
      if ((geom->mtype == CYLINDRICAL) || (geom->mtype == CONICAL)  ||
          (geom->mtype == SPHERICAL)   || (geom->mtype == TOROIDAL))
        correctUV(face, geom, bodydata->surfaces[iface+bodydata->nfaces],
                  bodydata->rvec[iface], uv);

      bodydata->faces_xyz[iface][3*ifp  ] = volumeMesh->node[ivp].xyz[0];
      bodydata->faces_xyz[iface][3*ifp+1] = volumeMesh->node[ivp].xyz[1];
      bodydata->faces_xyz[iface][3*ifp+2] = volumeMesh->node[ivp].xyz[2];

      bodydata->faces_uv[iface][2*ifp  ] = uv[0];
      bodydata->faces_uv[iface][2*ifp+1] = uv[1];
    }
  }

  for (iloop = 0; iloop < nloop; iloop++) {
    status = EG_getTopology(loops[iloop], &ref, &oclass, &mtype, limits, &nedge,
                            &edges, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (iedge = 0; iedge < nedge; iedge++) {

      // skip degenerate edges
      status = EG_getTopology(edges[iedge], &ref, &oclass, &mtype, trange,
                              &nnode, &nodes, &tmp);
      if (status != EGADS_SUCCESS) goto cleanup;
      if (mtype == DEGENERATE) continue;

      // try to get edge UV at t-mid with sense of 0
      offset = 0;
      status = EG_getEdgeUV(face, edges[iedge], 0, (trange[0]+trange[1])/2, uv);
      if (status == EGADS_TOPOERR) {
        /* sense in Face twice! */
        offset = nSurfPts;
      }

      edgeIndex = EG_indexBodyTopo(bodydata->body, edges[iedge]);

      exyz = bodydata->edges_xyz[edgeIndex-1];

      npts = bodydata->edges_npts[edgeIndex-1];
      for (i = 0; i < npts; i++) {

        t   = bodydata->edges_t[edgeIndex-1][i];
        isp = bodydata->edges_isp[edgeIndex-1][i] + offset;

        ifp = face_ind[isp];

        status = EG_getEdgeUV(face, edges[iedge], senses[iedge], t, uv);
        if (status != EGADS_SUCCESS) goto cleanup;

        bodydata->faces_xyz[iface][3*ifp  ] = exyz[3*i  ];
        bodydata->faces_xyz[iface][3*ifp+1] = exyz[3*i+1];
        bodydata->faces_xyz[iface][3*ifp+2] = exyz[3*i+2];

        bodydata->faces_uv[iface][2*ifp  ] = uv[0];
        bodydata->faces_uv[iface][2*ifp+1] = uv[1];
      }
    }
  }

  status = CAPS_SUCCESS;
  cleanup:
    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in getFaceUV status = %d\n", status);

    return status;
}
#endif


static void swapd(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

static void swapi(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}


#if 0
// A function to implement bubble sort
static void bubbleSort(int n, double t[], double xyz[], int isp[])
{

  int i, j;
  for (i = 0; i < n-1; i++)
    // Last i elements are already in place
    for (j = 0; j < n-i-1; j++)
      if (t[j] > t[j+1]) {
        swapd(&t[j]      , &t[j+1]        );
        swapd(&xyz[3*j+0], &xyz[3*(j+1)+0]);
        swapd(&xyz[3*j+1], &xyz[3*(j+1)+1]);
        swapd(&xyz[3*j+2], &xyz[3*(j+1)+2]);
        swapi(&isp[j]    , &isp[j+1]      );
      }
}
#endif


// A function to reverse an Edge tessellation if needed
static void orientEdgeTess(edgeData *tedge)
{

  int i;
  int n = tedge->npts;
  double *t = tedge->t;
  double *xyz = tedge->xyz;

  if (t[0] < t[n-1]) return;

  for (i = 0; i < n/2; i++) {
    swapd(&t[i]      , &t[n-1-i]        );
    swapd(&xyz[3*i+0], &xyz[3*(n-1-i)+0]);
    swapd(&xyz[3*i+1], &xyz[3*(n-1-i)+1]);
    swapd(&xyz[3*i+2], &xyz[3*(n-1-i)+2]);
  }
}


// Read and skips any comment lines in a gma file
static int nextHeader(char **line, size_t *nline, int *iline, FILE *fp)
{

    while (getline(line, nline, fp) != -1) {
        (*iline)++;
        if (nline == 0) continue;
        if ((*line)[0] == '#') continue;
        break;
    }
    if (feof(fp)) return CAPS_IOERR;

    return CAPS_SUCCESS;
}


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS;
    aimStorage *pointwiseInstance=NULL;

#ifdef DEBUG
    printf("\n pointwiseAIM/aimInitialize   instance = %d!\n", inst);
#endif

    // specify the number of analysis input and out "parameters"
    *nIn     = NUMINPUT;
    *nOut    = NUMOUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames  = NULL;
    *franks  = NULL;
    *fInOut  = NULL;

    // Allocate pointwiseInstance
    AIM_ALLOC(pointwiseInstance, 1, aimStorage, aimInfo, status);
    *instStore = pointwiseInstance;

    // Set initial values for pointwiseInstance
    initiate_aimStorage(pointwiseInstance);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);
    return status;
}


// Available AIM inputs
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{

    /*! \page aimInputsPointwise AIM Inputs
     * The following list outlines the Pointwise options along with their default value available
     * through the AIM interface.
     */

#ifdef DEBUG
    printf(" pointwiseAIM/aimInputs  index = %d!\n", index);
#endif

    // Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name"); // If NULL a volume grid won't be written by the AIM
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->vals.string  = NULL;
        defval->lfixed       = Change;

        /*! \page aimInputsPointwise
         * - <B> Proj_Name = NULL</B> <br>
         * This corresponds to the output name of the mesh. If left NULL, the mesh is not written to a file.
         */
    } else if (index == Mesh_Format) {
        *ainame               = EG_strdup("Mesh_Format");
        defval->type          = String;
        defval->vals.string   = EG_strdup("VTK"); // TECPLOT, VTK, AFLR3, STL, AF, FAST, NASTRAN
        defval->lfixed        = Change;
        defval->nullVal       = IsNull;

        /*! \page aimInputsPointwise
         * - <B> Mesh_Format = NULL</B> <br>
         * Mesh output format. Available format names include: "AFLR3", "VTK", "TECPLOT", "SU2", "Nastran".
         * This file format is written from CAPS, and is not the CAE solver in Pointwise.
         */
    } else if (index == Mesh_ASCII_Flag) {
        *ainame               = EG_strdup("Mesh_ASCII_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = true;

        /*! \page aimInputsPointwise
         * - <B> Mesh_ASCII_Flag = True</B> <br>
         * Output mesh in ASCII format, otherwise write a binary file, if applicable.
         */
    } else  if (index == Mesh_Sizing) {
        *ainame              = EG_strdup("Mesh_Sizing");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsPointwise
         * - <B>Mesh_Sizing = NULL </B> <br>
         * These parameters are implemented by overriding PW: attributes.
         * See \ref meshSizingProp for additional details.
         */
    } else if (index == Mesh_Length_Factor) {
        *ainame              = EG_strdup("Mesh_Length_Factor");
        defval->type         = Double;
        defval->dim          = Scalar;
        defval->vals.real    = 1;
        defval->nullVal      = NotNull;

        /*! \page aimInputsPointwise
         * - <B> Mesh_Length_Factor = 1</B> <br>
         * Scaling factor to compute a meshing Reference_Length via:<br>
         * <br>
         * Reference_Length = capsMeshLength*Mesh_Length_Factor<br>
         * <br>
         * Reference_Length scales all input parameters with units of length
         */
    } else if (index == Connector_Initial_Dim) {
        *ainame              = EG_strdup("Connector_Initial_Dim");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 11;

        /*! \page aimInputsPointwise
         * - <B> Connector_Initial_Dim = 11</B> <br>
         * Initial connector dimension.
         */
    } else if (index == Connector_Max_Dim) {
        *ainame              = EG_strdup("Connector_Max_Dim");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 1024;

        /*! \page aimInputsPointwise
         * - <B> Connector_Max_Dim = 1024</B> <br>
         * Maximum connector dimension.
         */
    } else if (index == Connector_Min_Dim) {
        *ainame              = EG_strdup("Connector_Min_Dim");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 4;

        /*! \page aimInputsPointwise
         * - <B> Connector_Min_Dim = 4</B> <br>
         * Minimum connector dimension.
         */
    } else if (index == Connector_Turn_Angle) {
        *ainame              = EG_strdup("Connector_Turn_Angle");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Connector_Turn_Angle = 0.0</B> <br>
         * Maximum turning angle [degree] on connectors for dimensioning (0 - not used).
         * Influences connector resolution in high curvature regions. Suggested values from 5 to 20 degrees.
         */
    } else if (index == Connector_Deviation) {
        *ainame              = EG_strdup("Connector_Deviation");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Connector_Deviation = 0.0</B> <br>
         * Maximum deviation on connectors for dimensioning (0 - not used).
         * This is the maximum distance between the center of a segment on the connector to the CAD surface.
         * Influences connector resolution in high curvature regions.
         */
    } else if (index == Connector_Split_Angle) {
        *ainame              = EG_strdup("Connector_Split_Angle");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Connector_Split_Angle = 0.0</B> <br>
         * Turning angle on connectors to split (0 - not used).
         */
    } else if (index == Connector_Turn_Angle_Hard) {
        *ainame              = EG_strdup("Connector_Turn_Angle_Hard");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 70;

        /*! \page aimInputsPointwise
         * - <B> Connector_Turn_Angle_Hard = 70</B> <br>
         * Hard edge turning angle [degree] limit for domain T-Rex (0.0 - not used).
         */
    } else if (index == Connector_Prox_Growth_Rate) {
        *ainame              = EG_strdup("Connector_Prox_Growth_Rate");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 1.3;

        /*! \page aimInputsPointwise
         * - <B> Connector_Prox_Growth_Rate = 1.3</B> <br>
         *  Connector proximity growth rate.
         */
    } else if (index == Connector_Adapt_Sources) {
        *ainame              = EG_strdup("Connector_Adapt_Sources");
        defval->type         = Boolean;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Connector_Adapt_Sources = False</B> <br>
         *  Compute sources using connectors.
         */
    } else if (index == Connector_Source_Spacing) {
        *ainame              = EG_strdup("Connector_Source_Spacing");
        defval->type         = Boolean;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Connector_Source_Spacing = False</B> <br>
         *  Use source cloud for adaptive pass on connectors V18.2+.
         */
    } else if (index == Domain_Algorithm) {
        *ainame              = EG_strdup("Domain_Algorithm");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Delaunay");

        /*! \page aimInputsPointwise
         * - <B> Domain_Algorithm = "Delaunay"</B> <br>
         * Isotropic (Delaunay, AdvancingFront or AdvancingFrontOrtho).
         */
    } else if (index == Domain_Full_Layers) {
        *ainame              = EG_strdup("Domain_Full_Layers");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Full_Layers = 0</B> <br>
         * Domain full layers (0 for multi-normals, >= 1 for single normal).
         */
    } else if (index == Domain_Max_Layers) {
        *ainame              = EG_strdup("Domain_Max_Layers");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Max_Layers = 0</B> <br>
         * Domain maximum layers.
         */
    } else if (index == Domain_Growth_Rate) {
        *ainame              = EG_strdup("Domain_Growth_Rate");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 1.3;

        /*! \page aimInputsPointwise
         * - <B> Domain_Growth_Rate = 1.3</B> <br>
         * Domain growth rate for 2D T-Rex extrusion.
         */
    } else if (index == Domain_Iso_Type) {
        *ainame              = EG_strdup("Domain_Iso_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Triangle");

        /*! \page aimInputsPointwise
         * - <B> Domain_Iso_Type = "Triangle"</B> <br>
         * Domain iso cell type (Triangle or TriangleQuad).
         */
    } else if (index == Domain_TRex_Type) {
        *ainame              = EG_strdup("Domain_TRex_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Triangle");

        /*! \page aimInputsPointwise
         * - <B> Domain_TRex_Type = "Triangle"</B> <br>
         * Domain T-Rex cell type (Triangle or TriangleQuad).
         */
    } else if (index == Domain_TRex_ARLimit) {
        *ainame              = EG_strdup("Domain_TRex_ARLimit");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 200.0;

        /*! \page aimInputsPointwise
         * - <B> Domain_TRex_ARLimit = 200.0</B> <br>
         * Domain T-Rex maximum aspect ratio limit (0 - not used).
         */
    } else if (index == Domain_TRex_AngleBC) {
        *ainame              = EG_strdup("Domain_TRex_AngleBC");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Domain_TRex_AngleBC = 0.0</B> <br>
         * Domain T-Rex spacing from surface curvature.
         */
    } else if (index == Domain_Decay) {
        *ainame              = EG_strdup("Domain_Decay");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.5;

        /*! \page aimInputsPointwise
         * - <B> Domain_Decay = 0.5</B> <br>
         * Domain boundary decay.
         */
    } else if (index == Domain_Min_Edge) {
        *ainame              = EG_strdup("Domain_Min_Edge");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Min_Edge = 0.0</B> <br>
         * Domain minimum edge length (relative to capsMeshLength).
         */
    } else if (index == Domain_Max_Edge) {
        *ainame              = EG_strdup("Domain_Max_Edge");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Max_Edge = 0.0</B> <br>
         * Domain minimum edge length (relative to capsMeshLength).
         */
    } else if (index == Domain_Adapt) {
        *ainame              = EG_strdup("Domain_Adapt");
        defval->type         = Boolean;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Adapt = False</B> <br>
         * Set up all domains for adaptation.
         */
    } else if (index == Domain_Wall_Spacing) {
        *ainame              = EG_strdup("Domain_Wall_Spacing");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Wall_Spacing = 0.0</B> <br>
         * Defined spacing when geometry attributed with PW:WallSpacing $wall  (relative to capsMeshLength)
         */
    } else if (index == Domain_Structure_AR_Convert) {
        *ainame              = EG_strdup("Domain_Structure_AR_Convert");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Domain_Structure_AR_Convert = 0.0</B> <br>
         * Aspect ratio to trigger converting domains to structured.
         */
    } else if (index == Block_Algorithm) {
        *ainame              = EG_strdup("Block_Algorithm");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Delaunay");

        /*! \page aimInputsPointwise
         * - <B> Domain_Algorithm = "Delaunay"</B> <br>
         * Isotropic (Delaunay, Voxel) (V18.3+).
         */
    } else if (index == Block_Voxel_Layers) {
        *ainame              = EG_strdup("Block_Voxel_Layers");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 3;

        /*! \page aimInputsPointwise
         * - <B> Block_Voxel_Layers = 3</B> <br>
         * Number of Voxel transition layers if Algorithm set to Voxel (V18.3+).
         */
    } else if (index == Block_Boundary_Decay) {
        *ainame              = EG_strdup("Block_Boundary_Decay");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.5;

        /*! \page aimInputsPointwise
         * - <B> Block_Boundary_Decay = 0.5</B> <br>
         * Volumetric boundary decay.
         */
    } else if (index == Block_Collision_Buffer) {
        *ainame              = EG_strdup("Block_Collision_Buffer");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.5;

        /*! \page aimInputsPointwise
         * - <B> Block_Collision_Buffer = 0.5</B> <br>
         * Collision buffer for colliding T-Rex fronts.
         */
    } else if (index == Block_Max_Skew_Angle) {
        *ainame              = EG_strdup("Block_Max_Skew_Angle");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 180.0;

        /*! \page aimInputsPointwise
         * - <B> Block_Max_Skew_Angle = 180.0</B> <br>
         * Maximum skew angle for T-Rex extrusion.
         */
    } else if (index == Block_TRex_Skew_Delay) {
        *ainame              = EG_strdup("Block_TRex_Skew_Delay");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Block_TRex_Skew_Delay = 0</B> <br>
         * Number of layers to delay enforcement of skew criteria
         */
    } else if (index == Block_Edge_Max_Growth_Rate) {
        *ainame              = EG_strdup("Block_Edge_Max_Growth_Rate");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 1.8;

        /*! \page aimInputsPointwise
         * - <B> Block_Edge_Max_Growth_Rate = 1.8</B> <br>
         * Volumetric edge ratio.
         */
    } else if (index == Block_Full_Layers) {
        *ainame              = EG_strdup("Block_Full_Layers");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Block_Full_Layers = 0</B> <br>
         * Full layers (0 for multi-normals, >= 1 for single normal).
         */
    } else if (index == Block_Max_Layers) {
        *ainame              = EG_strdup("Block_Max_Layers");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 0;

        /*! \page aimInputsPointwise
         * - <B> Block_Max_Layers = 0</B> <br>
         * Maximum layers.
         */
    } else if (index == Block_Growth_Rate) {
        *ainame              = EG_strdup("Block_Growth_Rate");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 1.3;

        /*! \page aimInputsPointwise
         * - <B> Block_Growth_Rate = 1.3</B> <br>
         * Growth rate for volume T-Rex extrusion.
         */
    } else if (index == Block_TRexType) {
        *ainame              = EG_strdup("Block_TRexType");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("TetPyramid");

        /*! \page aimInputsPointwise
         * - <B> Block_TRexType = "TetPyramid"</B> <br>
         * T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms).
         */
    } else if (index == Gen_Source_Box_Length_Scale) {
        *ainame              = EG_strdup("Gen_Source_Box_Length_Scale");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Gen_Source_Box_Length_Scale = 0.0</B> <br>
         * Length scale of enclosed viscous walls in source box (0 - no box)  (relative to capsMeshLength).
         */
    } else if (index == Gen_Source_Box_Direction) {
        *ainame              = EG_strdup("Gen_Source_Box_Direction");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Vector;
        defval->nrow          = 3;
        defval->ncol          = 1;
        defval->vals.reals    = (double *) EG_alloc(defval->nrow*sizeof(double));
        if (defval->vals.reals != NULL) {
            defval->vals.reals[0] = 1.0;
            defval->vals.reals[1] = 0.0;
            defval->vals.reals[2] = 0.0;
        } else return EGADS_MALLOC;

        /*! \page aimInputsPointwise
         * - <B> Gen_Source_Box_Direction = [1.0, 0.0, 0.0]</B> <br>
         * Principal direction vector (i.e. normalized freestream vector).
         */
    } else if (index == Gen_Source_Box_Angle) {
        *ainame              = EG_strdup("Gen_Source_Box_Angle");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.0;

        /*! \page aimInputsPointwise
         * - <B> Gen_Source_Box_Angle = 0.0</B> <br>
         * Angle for widening source box in the assigned direction.
         */
    } else if (index == Gen_Source_Growth_Factor) {
        *ainame              = EG_strdup("Gen_Source_Growth_Factor");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 10.0;

        /*! \page aimInputsPointwise
         * - <B> Gen_Source_Growth_Factor = 10.0</B> <br>
         * Growth rate for spacing value along box.
         */
    } else if (index == Elevate_Cost_Threshold) {
        *ainame              = EG_strdup("Elevate_Cost_Threshold");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.8;

        /*! \page aimInputsPointwise
         * - <B> Elevate_Cost_Threshold = 0.8</B> <br>
         * Cost convergence threshold.
         */
    } else if (index == Elevate_Max_Include_Angle) {
        *ainame              = EG_strdup("Elevate_Max_Include_Angle");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 175.0;

        /*! \page aimInputsPointwise
         * - <B> Elevate_Max_Include_Angle = 175.0</B> <br>
         * Maximum included angle tolerance.
         */
    } else if (index == Elevate_Relax) {
        *ainame              = EG_strdup("Elevate_Relax");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.05;

        /*! \page aimInputsPointwise
         * - <B> Elevate_Relax = 0.05</B> <br>
         * Iteration relaxation factor.
         */
    } else if (index == Elevate_Smoothing_Passes) {
        *ainame              = EG_strdup("Elevate_Smoothing_Passes");
        defval->type         = Integer;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.integer = 1000;

        /*! \page aimInputsPointwise
         * - <B> Elevate_Smoothing_Passes = 1000</B> <br>
         * Number of smoothing passes.
         */
    } else if (index == Elevate_WCN_Weight) {
        *ainame              = EG_strdup("Elevate_WCN_Weight");
        defval->type         = Double;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.real    = 0.5;

        /*! \page aimInputsPointwise
         * - <B> Elevate_WCN_Weight = 0.5</B> <br>
         * WCN cost component weighting factor.
         */
    } else if (index == Elevate_WCN_Mode) {
        *ainame              = EG_strdup("Elevate_WCN_Mode");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->units        = NULL;
        defval->lfixed       = Fixed;
        defval->dim          = Scalar;
        defval->vals.string  = EG_strdup("Calculate");

        /*! \page aimInputsPointwise
         * - <B> Elevate_WCN_Mode = "Calculate"</B> <br>
         * WCN weight factor method (UseValue or Calculate).
         */
    }

    return CAPS_SUCCESS;
}


// AIM preAnalysis function
int aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{
    int status; // Status return

    int i; // Index

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;

    double capsMeshLength = 0;

    const char ugrid[] = "caps.GeomToMesh.ugrid";
    const char ugrid_lb8[] = "caps.GeomToMesh.lb8.ugrid";

    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects
    ego *bodyCopy = NULL;

    ego context, model = NULL;
    aimStorage *pointwiseInstance;
    char aimEgadsFile[PATH_MAX];
    char aimFile[PATH_MAX];

    int quilting = (int)false;

    status = aim_rmFile(aimInfo, ugrid);
    AIM_STATUS(aimInfo, status);
    status = aim_rmFile(aimInfo, ugrid_lb8);
    AIM_STATUS(aimInfo, status);

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" pointwiseAIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif

    if (numBody <= 0 || bodies == NULL) {
#ifdef DEBUG
        printf(" pointwiseAIM/aimPreAnalysis No Bodies!\n");
#endif
        return CAPS_SOURCEERR;
    }
    if (aimInputs == NULL) return CAPS_NULLOBJ;
  
    pointwiseInstance = (aimStorage *) instStore;

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(pointwiseInstance);
    AIM_STATUS(aimInfo, status, "pointwiseAIM aimStorage cleanup!!!");

    // set the filename without extensions where the grid is written for solvers
    status = aim_file(aimInfo, "caps.GeomToMesh", aimFile);
    AIM_STATUS(aimInfo, status);
    AIM_STRDUP(pointwiseInstance->meshRef.fileName, aimFile, aimInfo, status);

    // remove previous mesh
    status = aim_deleteMeshes(aimInfo, &pointwiseInstance->meshRef);
    AIM_STATUS(aimInfo, status);

    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            2, // Only search down to the face level of the EGADS bodyIndex
                                            &pointwiseInstance->groupMap);
    AIM_STATUS(aimInfo, status);

    status = create_CAPSMeshAttrToIndexMap(numBody,
                                            bodies,
                                            3, // Only search down to the face level of the EGADS bodyIndex
                                            &pointwiseInstance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Get mesh sizing parameters
    if (aimInputs[Mesh_Sizing-1].nullVal != IsNull) {

        status = deprecate_SizingAttr(aimInfo,
                                      aimInputs[Mesh_Sizing-1].length,
                                      aimInputs[Mesh_Sizing-1].vals.tuple,
                                      &pointwiseInstance->meshMap,
                                      &pointwiseInstance->groupMap);
        AIM_STATUS(aimInfo, status);

        status = mesh_getSizingProp(aimInfo,
                                    aimInputs[Mesh_Sizing-1].length,
                                    aimInputs[Mesh_Sizing-1].vals.tuple,
                                    &pointwiseInstance->meshMap,
                                    &numMeshProp,
                                    &meshProp);
        AIM_STATUS(aimInfo, status);
    }

    // Get the capsMeshLenght if boundary layer meshing has been requested
    status = check_CAPSMeshLength(numBody, bodies, &capsMeshLength);

    if (capsMeshLength <= 0 || status != CAPS_SUCCESS) {
        if (status != CAPS_SUCCESS) {
          AIM_ERROR(aimInfo, "capsMeshLength is not set on any body.");
        } else {
          AIM_ERROR(aimInfo, "capsMeshLength: %f", capsMeshLength);
        }
        AIM_ADDLINE(aimInfo, "\nThe capsMeshLength attribute must present on at least one body.\n"
                             "\n"
                             "capsMeshLength should be a a positive value representative\n"
                             "of a characteristic length of the geometry,\n"
                             "e.g. the MAC of a wing or diameter of a fuselage.\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Scale the reference length
    capsMeshLength *= aimInputs[Mesh_Length_Factor-1].vals.real;

    bodyCopy = (ego *) EG_alloc(numBody*sizeof(ego));
    if (bodyCopy == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < numBody; i++) bodyCopy[i] = NULL;

    // Get context
    status = EG_getContext(bodies[0], &context);
    if (status != EGADS_SUCCESS) goto cleanup;

    // Make a copy of the bodies
    for (i = 0; i < numBody; i++) {
        status = EG_copyObject(bodies[i], NULL, &bodyCopy[i]);
        if (status != EGADS_SUCCESS) goto cleanup;
/*@-nullpass@*/
        status = setPWAttr(aimInfo, bodyCopy[i],
                           &pointwiseInstance->groupMap,
                           &pointwiseInstance->meshMap,
                           numMeshProp, meshProp, capsMeshLength, &quilting);
/*@+nullpass@*/
        if (status != EGADS_SUCCESS) goto cleanup;
    }

    // Auto quilt faces
    /*
    if (aimInputs[Auto_Quilt_Flag-1].vals.integer == (int) true) {
        printf("Automatically quilting faces...\n");
        for (i = 0; i < numBody; i++) {
            status = set_autoQuiltAttr(&bodyCopy[i]);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }
     */

    // Create a model from the copied bodies
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, numBody, bodyCopy,
                             NULL, &model);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (model == NULL) {
      status = CAPS_NULLOBJ;
      goto cleanup;
    }

    printf("Writing global Glyph inputs...\n");
    status = writeGlobalGlyph(aimInfo, aimInputs, capsMeshLength);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = aim_file(aimInfo, egadsFileName, aimEgadsFile);
    AIM_STATUS(aimInfo, status);

    printf("Writing egads file '%s'....\n",aimEgadsFile);
    remove(aimEgadsFile);
    status = EG_saveModel(model, aimEgadsFile);
    AIM_STATUS(aimInfo, status);

    if (quilting == (int)true) {
      AIM_ERROR  (aimInfo, "Quilting is enabled with 'PW:QuiltName' attribute on faces.");
      AIM_ADDLINE(aimInfo, "       Pointwise input files were generated, but CAPS cannot process the resulting grid.");
      status = CAPS_MISMATCH;
      goto cleanup;
    }

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS)
      printf("Error: pointwiseAIM status %d\n", status);

    // Clean up meshProps
    if (meshProp != NULL) {
        for (i = 0; i < numMeshProp; i++) {
            destroy_meshSizingStruct(&meshProp[i]);
        }
        EG_free(meshProp); meshProp = NULL;
    }

    // delete the model
    if (model != NULL) {
      EG_deleteObject(model);
    } else {
      if (bodyCopy != NULL) {
        for (i = 0;  i < numBody; i++) {
          if (bodyCopy[i] != NULL)  {
            (void) EG_deleteObject(bodyCopy[i]);
          }
        }
      }
    }
    EG_free(bodyCopy);

    return status;
}


// Read back in the resulting grid
int
aimPostAnalysis(void *instStore, void *aimInfo, /*@unused@*/ int restart,
                /*@unused@*/ capsValue *inputs)
{
    int status;

    char attrname[128];

    // Container for volume mesh
    int numVolumeMesh=0;
    meshStruct *volumeMesh=NULL;

    size_t     nline = 0;
    int        i, j, n, ib, ie, it, in;
    int        ivp, npts, nVolPts, iper, iline = 0, cID = 0, nglobal, iglobal;
    int        oclass, mtype, numBody = 0, *senses = NULL, *ivec = NULL;
    int        numFaces = 0, iface, numEdges = 0, ibody, iedge, numNodes = 0, nDegen = 0;
    int        ntri = 0, nquad = 0, elem[4], velem[4], *face_tris, elemIndex;
    int        GMA_MAJOR = 0, GMA_MINOR = 0, numConnector = 0, iCon, numDomain = 0, iDom;
    int        egadsID, edgeID, faceID, bodyID, *bodyIndex=NULL, *faceVertID=NULL;
    const char gmafilename[]   = "caps.GeomToMesh.gma";
    const char ugrid[] = "caps.GeomToMesh.ugrid";
    const char ugrid_lb8[] = "caps.GeomToMesh.lb8.ugrid";
    const char *intents = NULL, *groupName = NULL;
    char       aimFile[PATH_MAX];
    char       *line = NULL, aimEgadsFile[PATH_MAX];
    double     limits[4], trange[4], result[18], uv[2], offuv[2], t, ptol;
    double     *face_uv = NULL, *face_xyz = NULL;
#ifdef CHECKGRID
    double     d;
#endif
    double     *rvec, v1[3], v2[3], faceNormal[3], triNormal[3], ndot;
    const int    *tris = NULL, *tric = NULL, *ptype = NULL, *pindex = NULL;
    const double *pxyz = NULL, *puv = NULL;
    ego        geom, geomBspl, ref, *faces = NULL, *edges = NULL, *nodes = NULL, *bodies = NULL, *objs = NULL;
    ego        context, edge, face, tess=NULL, prev, next, model=NULL;
    bodyData   *bodydata = NULL;
    edgeMapData **edgeMap = NULL;
    meshStruct *surfaceMeshes = NULL;
//    hashElemTable table;
    FILE       *fp = NULL;

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char *string;
    aimStorage *pointwiseInstance;

//    initiate_hashTable(&table);

    if (aim_isFile(aimInfo, ugrid) == CAPS_SUCCESS) {
      status = aim_file(aimInfo, ugrid, aimFile);
      AIM_STATUS(aimInfo, status);
      status = aim_symLink(aimInfo, aimFile, ugrid_lb8);
      AIM_STATUS(aimInfo, status);
    } else if (aim_isFile(aimInfo, ugrid_lb8) == CAPS_SUCCESS) {
      status = aim_file(aimInfo, ugrid_lb8, aimFile);
      AIM_STATUS(aimInfo, status);
      status = aim_symLink(aimInfo, aimFile, ugrid);
      AIM_STATUS(aimInfo, status);
    } else {
      AIM_ERROR(aimInfo, "Pointwise did not generate %s or %s!\n\n", ugrid, ugrid_lb8);
      status = CAPS_IOERR;
      goto cleanup;
    }

    fp = aim_fopen(aimInfo, ugrid_lb8, "rb");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Failed to open %s!\n\n", ugrid_lb8);
      status = CAPS_IOERR;
      goto cleanup;
    }

    pointwiseInstance = (aimStorage *) instStore;

    numVolumeMesh = 1;
    AIM_ALLOC(volumeMesh, numVolumeMesh, meshStruct, aimInfo, status);

    status = initiate_meshStruct(volumeMesh);
    AIM_STATUS(aimInfo, status);

    status = getUGRID(fp, volumeMesh);
    AIM_STATUS(aimInfo, status);
    fclose(fp); fp = NULL;

    nVolPts = volumeMesh->numNode;

#if 0
    // construct the hash table into the surface elements to mark ID's
    allocate_hashTable(nVolPts, volumeMesh->meshQuickRef.numTriangle +
                                volumeMesh->meshQuickRef.numQuadrilateral, &table);
    for (i = 0; i < volumeMesh->meshQuickRef.numTriangle; i++) {
      elemIndex = i+volumeMesh->meshQuickRef.startIndexTriangle;
      status = hash_addElement(3, volumeMesh->element[elemIndex].connectivity,
                               i, &table);
      if (status != CAPS_SUCCESS) goto cleanup;
    }
    for (i = 0; i < volumeMesh->meshQuickRef.numQuadrilateral; i++) {
      elemIndex = i+volumeMesh->meshQuickRef.startIndexQuadrilateral;
      status = hash_addElement(4, volumeMesh->element[elemIndex].connectivity,
                               i+volumeMesh->meshQuickRef.numTriangle, &table);
      if (status != CAPS_SUCCESS) goto cleanup;
    }
#endif

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if ((numBody == 0) || (bodies == NULL)) {
        AIM_ERROR(aimInfo, "numBody = 0!");
        status = CAPS_BADOBJECT;
        goto cleanup;
    }
    AIM_ALLOC(bodyIndex, numBody, int, aimInfo, status);

    // Get context
    status = EG_getContext(bodies[0], &context);
    AIM_STATUS(aimInfo, status);

    status = aim_file(aimInfo, egadsFileName, aimEgadsFile);
    AIM_STATUS(aimInfo, status);

    // read back in the egads file as the bodies might be written in a different order
    status = EG_loadModel(context, 0, aimEgadsFile, &model);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(model, aimInfo, status);

    status = EG_getTopology(model, &ref, &oclass, &mtype, limits,
                            &n, &objs, &senses);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(objs, aimInfo, status);

    for (i = 0; i < numBody; i++) {
      status = EG_attributeRet(objs[i], "_body", &atype, &alen, &ints, &reals,
                               &string);
      if (status != EGADS_SUCCESS || alen != 1 || atype != ATTRINT) {
        AIM_ERROR(aimInfo, "_body attribute is not length 1 or not integer!");
        status = EGADS_ATTRERR;
        goto cleanup;
      }
      bodyID = ints[0];

      for (j = 0; j < numBody; j++) {
        status = EG_attributeRet(bodies[j], "_body", &atype, &alen, &ints, &reals,
                                 &string);
        if (status != EGADS_SUCCESS || alen != 1 || atype != ATTRINT) {
          AIM_ERROR(aimInfo, "_body attribute is not length 1 or not integer!");
          status = EGADS_ATTRERR;
          goto cleanup;
        }

        if (ints[0] == bodyID) {
          bodyIndex[i] = j;
          break;
        }
      }
    }
    EG_deleteObject(model); model=NULL;

    AIM_ALLOC(bodydata, numBody, bodyData, aimInfo, status);
    initiate_bodyData(numBody, bodydata);

    /* get all of the EGADS Objects */
    for (i = 0; i < numBody; i++) {
      bodydata[i].body = bodies[i];

#ifdef DEBUG
      {
        status = EG_attributeRet(bodies[i], "_name", &atype, &alen, &ints,
                                 &reals, &string);
        if (status == EGADS_SUCCESS) {
          printf("Body %d = %s\n", i+1, string);
        }
      }
#endif

      status = EG_getBodyTopos(bodies[i], NULL, NODE, &bodydata[i].nnodes,
                               &bodydata[i].nodes);
      AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(bodies[i], NULL, EDGE, &bodydata[i].nedges,
                               &bodydata[i].edges);
      AIM_STATUS(aimInfo, status);

      numEdges += bodydata[i].nedges; // Accumulate the total number of edges

      status = EG_getBodyTopos(bodies[i], NULL, FACE, &bodydata[i].nfaces,
                               &bodydata[i].faces);
      AIM_STATUS(aimInfo, status);

      numFaces += bodydata[i].nfaces; // Accumulate the total number of faces

      AIM_ALLOC(bodydata[i].surfaces, 2*bodydata[i].nfaces, ego, aimInfo, status);
      for (j = 0; j < 2*bodydata[i].nfaces; j++) bodydata[i].surfaces[j] = NULL;

      AIM_ALLOC(bodydata[i].rvec, bodydata[i].nfaces, double*, aimInfo, status);
      for (j = 0; j < bodydata[i].nfaces; j++) bodydata[i].rvec[j] = NULL;

      AIM_ALLOC(bodydata[i].tedges, bodydata[i].nedges, edgeData, aimInfo, status);
      for (j = 0; j < bodydata[i].nedges; j++) {
          bodydata[i].tedges[j].npts = 0;
          bodydata[i].tedges[j].xyz  = NULL;
          bodydata[i].tedges[j].t    = NULL;
          bodydata[i].tedges[j].ivp  = NULL;
      }

      AIM_ALLOC(bodydata[i].tfaces, bodydata[i].nfaces, faceData, aimInfo, status);
      for (j = 0; j < bodydata[i].nfaces; j++) {
          bodydata[i].tfaces[j].npts  = 0;
          bodydata[i].tfaces[j].xyz   = NULL;
          bodydata[i].tfaces[j].uv    = NULL;
          bodydata[i].tfaces[j].ntri  = 0;
          bodydata[i].tfaces[j].nquad = 0;
          bodydata[i].tfaces[j].tris  = NULL;
          bodydata[i].tfaces[j].ivp   = NULL;
      }

      for (j = 0; j < bodydata[i].nfaces; j++) {
        status = EG_getTopology(bodydata[i].faces[j], &bodydata[i].surfaces[j],
                                &oclass, &mtype, limits, &n, &objs, &senses);
        AIM_STATUS(aimInfo, status);

        geom = bodydata[i].surfaces[j];
        if (geom->mtype == TRIMMED) {
          status = EG_getGeometry(geom, &oclass, &mtype, &ref, &ivec, &rvec);
          AIM_STATUS(aimInfo, status);

          bodydata[i].surfaces[j] = ref;
          EG_free(ivec);
          EG_free(rvec);
          geom = ref;
        }

        status = EG_getGeometry(geom, &oclass, &mtype, &ref, &ivec,
                                &bodydata[i].rvec[j]);
        AIM_STATUS(aimInfo, status);

        AIM_FREE(ivec);

        if (mtype != BSPLINE) {
          status = EG_convertToBSpline(bodydata[i].faces[j],
                                       &bodydata[i].surfaces[j+bodydata[i].nfaces]);
          AIM_STATUS(aimInfo, status, "Face %d ConvertToBSpline", j+1);
        }
      }
    }

    // Count the number of degenerate edges in all bodies
    nDegen = 0;
    for (ib = 0; ib < numBody; ib++) {
      for ( iedge = 0; iedge < bodydata[ib].nedges; iedge++ ) {
        status = EG_getTopology( bodydata[ib].edges[iedge], &ref, &oclass,
                                &mtype, limits, &n, &objs, &senses);
        AIM_STATUS(aimInfo, status);
        if ( mtype == DEGENERATE ) nDegen++;
      }
    }


    /* open and parse the gma file to count surface/edge tessellations points */
    fp = aim_fopen(aimInfo, gmafilename, "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open file: %s!", gmafilename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Read the gma file version
    status = nextHeader(&line, &nline, &iline, fp);
    if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "line %d Could not find next header!", iline);
        goto cleanup;
    }
/*@-nullpass@*/
    status = sscanf(line, "%d %d", &GMA_MAJOR, &GMA_MINOR);
/*@+nullpass@*/
    if (status != 2) {
        AIM_ERROR(aimInfo, "line %d Cannot get gma version number!", iline);
        goto cleanup;
    }
    if (!(GMA_MAJOR == 2 && GMA_MINOR == 0)) {
        AIM_ERROR(aimInfo, "line %d Cannot read gma file version %d.%d",
                  iline, GMA_MAJOR, GMA_MINOR);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // Read the number of connectors and domains
    status = nextHeader(&line, &nline, &iline, fp);
    if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "line %d Could not find next header!", iline);
        goto cleanup;
    }
/*@-nullpass@*/
    status = sscanf(line, "%d %d", &numConnector, &numDomain);
/*@+nullpass@*/
    if (status != 2) {
        AIM_ERROR(aimInfo, "line %d Cannot get Connector and Domain count!", iline);
        goto cleanup;
    }
#if 0
    if (numConnector != numEdges-nDegen) {
        printf(" Error: Number of Connector (%d) must match Edge count (%d)\n",
               numConnector, numEdges-nDegen);
        status = CAPS_IOERR;
        goto cleanup;
    }
#endif
    if (numDomain != numFaces) {
        AIM_ERROR(aimInfo, "Number of Domains (%d) must match Face count (%d)",
                  numDomain, numFaces);
        status = CAPS_IOERR;
        goto cleanup;
    }


    for (iCon = 0; iCon < numConnector; iCon++) {
        // Read the EGADS ID and the number of vertexes
        status = nextHeader(&line, &nline, &iline, fp);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "line %d Could not find next Connector header!", iline);
            goto cleanup;
        }
/*@-nullpass@*/
        status = sscanf(line, "%d %d", &egadsID, &npts);
/*@+nullpass@*/
        if (status != 2) {
            AIM_ERROR(aimInfo, "read line %d return = %d", iline, status);
            goto cleanup;
        }
        decodeEgadsID(egadsID, &it, &ib, &in);
        if (ib < 0 || ib >= numBody) {
          AIM_ERROR(aimInfo, "line %d Body ID (%d) out of bounds (0 - %d)!\n",
                 iline, ib, numBody);
          status = CAPS_MISMATCH;
          goto cleanup;
        }
        ib = bodyIndex[ib];

        // Check the ID type
        if (it != EDGEID) {
            AIM_ERROR(aimInfo, "line %d Expected Edge ID!", iline);
            goto cleanup;
        }
        // Check the Edge index
        if ((in < 0) || (in >= bodydata[ib].nedges)) {
            AIM_ERROR(aimInfo, "line %d Bad Edge index = %d [1-%d]!\n",
                   iline, in+1, bodydata[ib].nedges);
            status = CAPS_MISMATCH;
            goto cleanup;
        }
        // skip the comment line
        status = getline(&line, &nline, fp); iline++;
        if (status == -1) {
          AIM_ERROR(aimInfo, "line %d Failed to read comment!", iline);
          status = CAPS_IOERR;
          goto cleanup;
        }

        edge   = bodydata[ib].edges[in];
        status = EG_getTopology(edge, &geom, &oclass, &mtype, trange, &numNodes,
                                &nodes, &senses);
        AIM_STATUS(aimInfo, status, "line %d Bad Edge status = %d!", iline, status);

        if (geom->mtype == BSPLINE) {
            status = EG_getRange(geom, trange, &iper);
            AIM_STATUS(aimInfo, status, "line %d EG_getRange C = %d!", iline, status);
        }

        // get the tolerance of the edge
        status = EG_getTolerance(bodydata[ib].edges[in], &ptol);
        if (status != EGADS_SUCCESS) goto cleanup;

        bodydata[ib].tedges[in].npts = npts;
        AIM_ALLOC(bodydata[ib].tedges[in].xyz, 3*npts, double, aimInfo, status);
        AIM_ALLOC(bodydata[ib].tedges[in].t  ,   npts, double, aimInfo, status);
        AIM_ALLOC(bodydata[ib].tedges[in].ivp,   npts, int   , aimInfo, status);

        for (i = 0; i < npts; i++) {
            status = fscanf(fp, "%d %d %lf\n", &ivp, &edgeID, &t); iline++;
            if (status != 3) {
                AIM_ERROR(aimInfo, "read line %d return = %d", iline, status);
                status = CAPS_IOERR;
                goto cleanup;
            }
            if (edgeID != egadsID) {
                AIM_ERROR(aimInfo, "line %d Connector vertex ID (%d) does not match Edge ID (%d)\n",
                       iline, edgeID, egadsID);
                status = CAPS_MISMATCH;
                goto cleanup;
            }
            bodydata[ib].tedges[in].ivp[i] = ivp;
            ivp--;

            /* t-values from pointwise are always normalized in the range 0-1 */
            t = trange[0] + t*(trange[1]-trange[0]);

            bodydata[ib].tedges[in].xyz[3*i+0] = volumeMesh->node[ivp].xyz[0];
            bodydata[ib].tedges[in].xyz[3*i+1] = volumeMesh->node[ivp].xyz[1];
            bodydata[ib].tedges[in].xyz[3*i+2] = volumeMesh->node[ivp].xyz[2];
            bodydata[ib].tedges[in].t[i] = t;
#ifdef CHECKGRID
            status = EG_evaluate(edge, &t, result);
            AIM_STATUS(aimInfo, status);

            d = dist_DoubleVal(result, volumeMesh->node[ivp].xyz);
            if (d > MAX(1e-5,ptol)) {
              printf(" line %3d: Body %d Edge %d Edge deviation = %le  t = %lf mtype = %d\n",
                                 iline, ib+1, in+1, d, t, geom->mtype);

              status = EG_invEvaluate(edge, volumeMesh->node[ivp].xyz, &t,
                                      result);
              AIM_STATUS(aimInfo, status);

              d = dist_DoubleVal(result, volumeMesh->node[ivp].xyz);
              printf("           %d  %le  %lf [%lf %lf]\n",
                     status, d, t, limits[0], limits[1]);
            }
#endif
        }

        orientEdgeTess(&bodydata[ib].tedges[in]);
    }

    elemIndex = volumeMesh->numElement;

    for (iDom = 0; iDom < numDomain; iDom++) {
        // Read the EGADS ID and the number of vertexes, triangles, and quads
        status = nextHeader(&line, &nline, &iline, fp);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "line %d Could not find next Domain header!", iline);
            goto cleanup;
        }

/*@-nullpass@*/
        status = sscanf(line, "%d %d %d %d", &egadsID, &npts, &ntri, &nquad);
/*@+nullpass@*/
        if (status != 4) {
            AIM_ERROR(aimInfo, "read line %d return = %d", iline, status);
            status = CAPS_MISMATCH;
            goto cleanup;
        }
        decodeEgadsID(egadsID, &it, &ib, &in);
        ib = bodyIndex[ib];

        // Check the ID type
        if (it != FACEID) {
            AIM_ERROR(aimInfo, "line %d Expected Face ID!", iline);
            status = CAPS_MISMATCH;
            goto cleanup;
        }
        // Check the Edge index
        if ((in < 0) || (in >= bodydata[ib].nfaces)) {
            AIM_ERROR(aimInfo, "line %d Bad Face index = %d [1-%d]!", iline, in+1, bodydata[ib].nfaces);
            status = CAPS_MISMATCH;
            goto cleanup;
        }
        // skip the comment line

        status = getline(&line, &nline, fp); iline++;
        if (status == -1 || nline == 0 || line == NULL || line[0] != '#') {
          AIM_ERROR(aimInfo, "line %d Failed to read comment!", iline);
          status = CAPS_IOERR;
          goto cleanup;
        }
        if (nquad != 0) {
          AIM_ERROR(aimInfo, "line %d Quads are currently not supported!", iline);
          status = CAPS_MISMATCH;
          goto cleanup;
        }

        // allocate new surface elements
        volumeMesh->numElement += ntri+nquad;
        AIM_REALL(volumeMesh->element, volumeMesh->numElement, meshElementStruct, aimInfo, status);

        // initialize the new elements
        for (i = elemIndex; i < volumeMesh->numElement; i++ ) {
          status = initiate_meshElementStruct(&volumeMesh->element[i],
                                              volumeMesh->analysisType);
          AIM_STATUS(aimInfo, status);
        }

        face = bodydata[ib].faces[in];
        geom = bodydata[ib].surfaces[in];

        if (geom->mtype == BSPLINE) {
          geomBspl = geom;
        } else {
          geomBspl = bodydata[ib].surfaces[in+bodydata[ib].nfaces];
        }
        status = EG_getRange(geomBspl, limits, &iper);
        AIM_STATUS(aimInfo, status);

        /* correct UV-values for periodic shapes */
        status = correctUV(face, geom, geomBspl, bodydata[ib].rvec[in], offuv);
        AIM_STATUS(aimInfo, status);

        // get the tolerance of the face
        status = EG_getTolerance(face, &ptol);
        AIM_STATUS(aimInfo, status);

        // Look for component/boundary ID for attribute mapper based on capsGroup
        status = retrieve_CAPSGroupAttr(face, &groupName);
        if (status != CAPS_SUCCESS) {
          AIM_ERROR(aimInfo, "No capsGroup attribute found on Face %d, unable to assign a boundary index value",
                    in+1);
          print_AllAttr( aimInfo, face );
          goto cleanup;
        }

/*@-nullpass@*/
        status = get_mapAttrToIndexIndex(&pointwiseInstance->groupMap,
                                         groupName, &cID);
        AIM_STATUS(aimInfo, status, "Unable to retrieve boundary index from capsGroup: %s",
                   groupName);
/*@+nullpass@*/

        bodydata[ib].tfaces[in].npts  = npts;
        bodydata[ib].tfaces[in].ntri  = ntri;
        bodydata[ib].tfaces[in].nquad = nquad;
        AIM_ALLOC(bodydata[ib].tfaces[in].xyz , 3*npts, double, aimInfo, status);
        AIM_ALLOC(bodydata[ib].tfaces[in].uv  , 2*npts, double, aimInfo, status);
        AIM_ALLOC(bodydata[ib].tfaces[in].tris, 3*ntri, int   , aimInfo, status);
        AIM_ALLOC(bodydata[ib].tfaces[in].ivp ,   npts, int   , aimInfo, status);

        for (i = 0; i < npts; i++) {
          status = getline(&line, &nline, fp); iline++;
          if (status == -1 || nline == 0) {
            AIM_ERROR(aimInfo, "line %d Failed to read line!", iline);
            status = CAPS_IOERR;
            goto cleanup;
          }
/*@-nullpass@*/
          status = sscanf(line, "%d %d %lf %lf", &ivp, &faceID, &uv[0], &uv[1]);
/*@+nullpass@*/
          if (status != 4) {
            AIM_ERROR(aimInfo, "read line %d return = %d (!= 4)", iline, status);
            status = CAPS_IOERR;
            goto cleanup;
          }
          if (faceID != egadsID) {
            AIM_ERROR(aimInfo, "line %d Domain vertex ID (%d) does not match Face ID (%d)\n",
                      iline, faceID, egadsID);
            status = CAPS_MISMATCH;
            goto cleanup;
          }
          bodydata[ib].tfaces[in].ivp[i] = ivp;
          ivp--;

          /* UV-values from pointwise are always normalized in the range 0-1 */
          uv[0] = limits[0] + uv[0]*(limits[1]-limits[0]) + offuv[0];
          uv[1] = limits[2] + uv[1]*(limits[3]-limits[2]) + offuv[1];

          bodydata[ib].tfaces[in].xyz[3*i+0] = volumeMesh->node[ivp].xyz[0];
          bodydata[ib].tfaces[in].xyz[3*i+1] = volumeMesh->node[ivp].xyz[1];
          bodydata[ib].tfaces[in].xyz[3*i+2] = volumeMesh->node[ivp].xyz[2];
          bodydata[ib].tfaces[in].uv[2*i+0] = uv[0];
          bodydata[ib].tfaces[in].uv[2*i+1] = uv[1];

#ifdef CHECKGRID
          // check the face vertexes
          status = EG_evaluate(face, uv, result);
          AIM_STATUS(aimInfo, status);

          d = dist_DoubleVal(result, volumeMesh->node[ivp].xyz);
          if (d > MAX(1e-5,ptol)) {
            printf(" line %d: %d %d Face deviation = %le  %d  %d\n",
                   iline, ib+1, in+1, d, geom->mtype, egadsID);
          }
#endif
        }

        // skip the comment line
        status = getline(&line, &nline, fp); iline++;
        if (status == -1 || line == NULL || line[0] != '#') {
          AIM_ERROR(aimInfo, "line %d Failed to read comment!", iline);
          status = CAPS_IOERR;
          goto cleanup;
        }

        for (i = 0; i < ntri; i++) {
          status = fscanf(fp, "%d %d %d\n", &elem[0], &elem[1], &elem[2]);
          iline++;
          if (status != 3) {
            AIM_ERROR(aimInfo, "read line %d return = %d", iline, status);
            status = CAPS_IOERR;
            goto cleanup;
          }
          for (j = 0; j < 3; j++) {
            if ((elem[j] < 1) || (elem[j] > npts)) {
              AIM_ERROR(aimInfo, "line %d Bad Vertex index = %d [1-%d]!", iline, elem[j], nVolPts);
              status = CAPS_MISMATCH;
              goto cleanup;
            }
          }

          // these triangles map into surfacedata
          bodydata[ib].tfaces[in].tris[3*i  ] = elem[0];
          bodydata[ib].tfaces[in].tris[3*i+1] = elem[1];
          bodydata[ib].tfaces[in].tris[3*i+2] = elem[2];


          // get the volume index of the triangle
          velem[0] = bodydata[ib].tfaces[in].ivp[elem[0]-1];
          velem[1] = bodydata[ib].tfaces[in].ivp[elem[1]-1];
          velem[2] = bodydata[ib].tfaces[in].ivp[elem[2]-1];

          volumeMesh->element[elemIndex].elementType = Triangle;
          volumeMesh->element[elemIndex].elementID   = elemIndex+1;

          status = mesh_allocMeshElementConnectivity(&volumeMesh->element[elemIndex]);
          if (status != CAPS_SUCCESS) goto cleanup;

          volumeMesh->element[elemIndex].connectivity[0] = velem[0];
          volumeMesh->element[elemIndex].connectivity[1] = velem[1];
          volumeMesh->element[elemIndex].connectivity[2] = velem[2];

          volumeMesh->element[elemIndex].markerID = cID;
          volumeMesh->element[elemIndex].topoIndex = in+1;
          elemIndex++;
        }
    }

    // done reading the file
    fclose(fp); fp = NULL;

    // generate the QuickRefLists now that all surface elements have been added
    status = mesh_fillQuickRefList(volumeMesh);
    AIM_STATUS(aimInfo, status);

    // construct a map between coincident Solid/Sheet body Edges
    status = matchSameEdges(bodies, numBody, &edgeMap);
    AIM_STATUS(aimInfo, status);

    // Allocate surfaceMesh from number of bodies
    AIM_ALLOC(volumeMesh->referenceMesh, numBody, meshStruct, aimInfo, status);
    volumeMesh->numReferenceMesh = numBody;
    surfaceMeshes = volumeMesh->referenceMesh;

    // Initiate surface meshes
    for (ib = 0; ib < numBody; ib++){
        status = initiate_meshStruct(&surfaceMeshes[ib]);
        AIM_STATUS(aimInfo, status);
    }

    // Allocate surfaceMesh from number of bodies
    AIM_ALLOC(pointwiseInstance->meshRef.maps, numBody, aimMeshTessMap, aimInfo, status);
    pointwiseInstance->meshRef.nmap = numBody;
    for (ib = 0; ib < numBody; ib++){
        pointwiseInstance->meshRef.maps[ib].tess = NULL;
        pointwiseInstance->meshRef.maps[ib].map = NULL;
    }


    // populate the tess objects
    for (ib = 0; ib < numBody; ib++) {

        // Build up the body tessellation object
        status = EG_initTessBody(bodies[ib], &tess);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(tess, aimInfo, status);

        for ( ie = 0; ie < bodydata[ib].nedges; ie++ ) {

            // Check if the edge is degenerate
            status = EG_getTopology(bodydata[ib].edges[ie], &ref, &oclass,
                                    &mtype, limits, &n, &objs, &senses);
            AIM_STATUS(aimInfo, status);
            if (mtype == DEGENERATE) continue;

            // set the edge tessellation on the tess object
            if (bodydata[ib].tedges[ie].npts > 0) {
              ibody = ib;
              iedge = ie;
            } else {
              AIM_NOTNULL(edgeMap, aimInfo, status);
              ibody = edgeMap[ib][ie].ibody;
              iedge = edgeMap[ib][ie].iedge;
            }

            status = EG_setTessEdge(tess, ie+1, bodydata[ibody].tedges[iedge].npts,
                                    bodydata[ibody].tedges[iedge].xyz,
                                    bodydata[ibody].tedges[iedge].t);
            AIM_STATUS(aimInfo, status, "Failed to set tessellation on Edge %d!", iedge+1);

            // Add the unique indexing of the edge tessellation
            snprintf(attrname, 128, "edgeVertID_%d",ie+1);
            status = EG_attributeAdd(tess, attrname, ATTRINT,
                                     bodydata[ibody].tedges[iedge].npts,
                                     bodydata[ibody].tedges[iedge].ivp, NULL, NULL);
            AIM_STATUS(aimInfo, status);
       }

        for (iface = 0; iface < bodydata[ib].nfaces; iface++) {

            ntri  = bodydata[ib].tfaces[iface].ntri;
            nquad = bodydata[ib].tfaces[iface].nquad;

            face_tris = bodydata[ib].tfaces[iface].tris;
            face_uv   = bodydata[ib].tfaces[iface].uv;
            face_xyz  = bodydata[ib].tfaces[iface].xyz;

            AIM_NOTNULL(face_tris, aimInfo, status);
            AIM_NOTNULL(face_uv  , aimInfo, status);
            AIM_NOTNULL(face_xyz , aimInfo, status);

            // check the normals of the elements match the geometry normals
            // only need to check one element per face to decide for all
            elem[0] = face_tris[0]-1;
            elem[1] = face_tris[1]-1;
            elem[2] = face_tris[2]-1;

            // get the uv centroid
            uv[0] = (face_uv[2*elem[0]+0] + face_uv[2*elem[1]+0] + face_uv[2*elem[2]+0])/3.;
            uv[1] = (face_uv[2*elem[0]+1] + face_uv[2*elem[1]+1] + face_uv[2*elem[2]+1])/3.;

            // get the normal of the face
            status = EG_evaluate(bodydata[ib].faces[iface], uv, result);
            AIM_STATUS(aimInfo, status);

            // use cross dX/du x dX/dv to get geometry normal
            v1[0] = result[3]; v1[1] = result[4]; v1[2] = result[5];
            v2[0] = result[6]; v2[1] = result[7]; v2[2] = result[8];
            CROSS(faceNormal, v1, v2);

            // get mtype=SFORWARD or mtype=SREVERSE for the face to get topology normal
            status = EG_getInfo(bodydata[ib].faces[iface], &oclass, &mtype, &ref, &prev, &next);
            if (status != EGADS_SUCCESS) goto cleanup;
            faceNormal[0] *= mtype;
            faceNormal[1] *= mtype;
            faceNormal[2] *= mtype;

            // get the normal of the mesh triangle
            v1[0] = face_xyz[3*elem[1]+0] - face_xyz[3*elem[0]+0];
            v1[1] = face_xyz[3*elem[1]+1] - face_xyz[3*elem[0]+1];
            v1[2] = face_xyz[3*elem[1]+2] - face_xyz[3*elem[0]+2];

            v2[0] = face_xyz[3*elem[2]+0] - face_xyz[3*elem[0]+0];
            v2[1] = face_xyz[3*elem[2]+1] - face_xyz[3*elem[0]+1];
            v2[2] = face_xyz[3*elem[2]+2] - face_xyz[3*elem[0]+2];

            CROSS(triNormal, v1, v2);

            // get the dot product between the triangle and face normals
            ndot = DOT(faceNormal,triNormal);

            // if the normals are opposite, swap all triangles
            if (ndot < 0) {
              // swap two vertices to reverse the normal
              for (i = 0; i < ntri; i++) {
                swapi(&face_tris[3*i+0], &face_tris[3*i+2]);
              }
            }

            status = EG_setTessFace(tess,
                                    iface+1,
                                    bodydata[ib].tfaces[iface].npts,
                                    face_xyz,
                                    face_uv,
                                    ntri,
                                    face_tris);
            AIM_STATUS(aimInfo, status);


            // The points get reindexed to be consistent with loops in EG_setTessFace
            // This uses the new triangulation to map that index change.
            status = EG_getTessFace(tess, iface+1, &npts, &pxyz, &puv, &ptype,
                                    &pindex, &ntri, &tris, &tric);
            AIM_STATUS(aimInfo, status);
            AIM_NOTNULL(tris, aimInfo, status);

            AIM_ALLOC(faceVertID, bodydata[ib].tfaces[iface].npts, int, aimInfo, status);

            for (i = 0; i < ntri; i++) {
              for (j = 0; j < 3; j++) {
                faceVertID[tris[3*i+j]-1] = bodydata[ib].tfaces[iface].ivp[face_tris[3*i+j]-1];
              }
            }

            // Add the unique indexing of the tessellation
            snprintf(attrname, 128, "faceVertID_%d",iface+1);
            status = EG_attributeAdd(tess, attrname, ATTRINT,
                                     bodydata[ib].tfaces[iface].npts,
                                     faceVertID, NULL, NULL);
            AIM_STATUS(aimInfo, status);

            // replace the shuffled volume ID's
            AIM_FREE(bodydata[ib].tfaces[iface].ivp);
            bodydata[ib].tfaces[iface].ivp = faceVertID;
            faceVertID = NULL;
        }

        // finalize the tessellation
        status = EG_statusTessBody(tess, &bodies[ib], &i, &nglobal);
        AIM_STATUS(aimInfo, status, "Tessellation object was not built correctly!!!");

        status = copy_mapAttrToIndexStruct( &pointwiseInstance->groupMap,
                                            &surfaceMeshes[ib].groupMap );
        AIM_STATUS(aimInfo, status);

        // Create the map from the tessellation global vertex index to the volume mesh vertex index
        AIM_ALLOC(pointwiseInstance->meshRef.maps[ib].map, nglobal, int, aimInfo, status);

        for (iface = 0; iface < bodydata[ib].nfaces; iface++) {
          status = EG_getTessFace(tess, iface+1, &npts, &pxyz, &puv, &ptype,
                                  &pindex, &ntri, &tris, &tric);
          AIM_STATUS(aimInfo, status);

          /* construct global vertex indices */
          for (i = 0; i < npts; i++) {
            status = EG_localToGlobal(tess, iface+1, i+1, &iglobal);
            AIM_STATUS(aimInfo, status);
            pointwiseInstance->meshRef.maps[ib].map[iglobal-1] = bodydata[ib].tfaces[iface].ivp[i];
          }
        }

        // save the tessellation with caps
        status = aim_newTess(aimInfo, tess);
        AIM_STATUS(aimInfo, status);

/*@-kepttrans@*/
        // reference the surface mesh object
        surfaceMeshes[ib].egadsTess = tess;
        pointwiseInstance->meshRef.maps[ib].tess = tess;
        tess = NULL;
/*@+kepttrans@*/

        status = mesh_surfaceMeshEGADSTess(aimInfo, &surfaceMeshes[ib]);
        AIM_STATUS(aimInfo, status);

        printf("Body = %d\n", ib+1);
        printf("\tNumber of nodes    = %d\n", surfaceMeshes[ib].numNode);
        printf("\tNumber of elements = %d\n", surfaceMeshes[ib].numElement);
        if (surfaceMeshes[ib].meshQuickRef.useStartIndex == (int) true ||
            surfaceMeshes[ib].meshQuickRef.useListIndex == (int) true) {
            printf("\tNumber of tris = %d\n", surfaceMeshes[ib].meshQuickRef.numTriangle);
            printf("\tNumber of quad = %d\n", surfaceMeshes[ib].meshQuickRef.numQuadrilateral);
        }
    }

    // write out the mesh if requested
    status = writeMesh(aimInfo, pointwiseInstance,
                       numVolumeMesh, volumeMesh);
    AIM_STATUS(aimInfo, status);

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) {
      AIM_ADDLINE(aimInfo, "Please make sure you are using Pointwise V18.4 or newer.\n");
    }

    // Destroy volume mesh allocated arrays
    for (i = 0; i < numVolumeMesh; i++) {
        for (j = 0; j < volumeMesh[i].numReferenceMesh; j++) {
            (void) destroy_meshStruct(&volumeMesh[i].referenceMesh[j]);
        }
        (void) destroy_meshStruct(&volumeMesh[i]);
        AIM_FREE(volumeMesh);
    }

    AIM_FREE(bodyIndex);
    if (model != NULL) EG_deleteObject(model);
    if (tess  != NULL) EG_deleteObject(tess);
    destroy_bodyData(numBody, bodydata);
    AIM_FREE(bodydata);
    AIM_FREE(faceVertID);
    AIM_FREE(faces);
    AIM_FREE(edges);
    if (edgeMap != NULL) {
      for (i = 0; i < numBody; i++)
        AIM_FREE(edgeMap[i]);
      AIM_FREE(edgeMap);
    }
    if (fp != NULL) fclose(fp);

    if(line != NULL) free(line); // allocated by getline, must use free

    return status;
}


// Available AIM outputs
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsPointwise AIM Outputs
     * The following list outlines the Pointwise AIM outputs available through the AIM interface.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" pointwiseAIM/aimOutputs  index = %d!\n", index);
#endif

    if (index == Volume_Mesh) {

        *aoname           = AIM_NAME(Volume_Mesh);
        form->type        = PointerMesh;
        form->dim         = Scalar;
        form->lfixed      = Fixed;
        form->sfixed      = Fixed;
        form->vals.AIMptr = NULL;
        form->nullVal     = IsNull;

        /*! \page aimOutputsPointwise
         * - <B> Volume_Mesh </B> <br>
         * The volume mesh for a link
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);
    }

    AIM_NOTNULL(*aoname, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*aoname);
    return status;
}


// Get value for a given output variable
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ int index, capsValue *val)
{
    int        status = CAPS_SUCCESS;
    aimStorage *pointwiseInstance;
    aimMesh    mesh;

#ifdef DEBUG
    printf(" pointwiseAIM/aimCalcOutput  index = %d!\n", index);
#endif
    pointwiseInstance = (aimStorage *) instStore;

    if (Volume_Mesh == index) {

        status = aim_queryMeshes( aimInfo, Volume_Mesh, &pointwiseInstance->meshRef );
        if (status > 0) {

/*@-immediatetrans@*/
          mesh.meshData = NULL;
          mesh.meshRef = &pointwiseInstance->meshRef;
/*@+immediatetrans@*/

          status = aim_readBinaryUgrid(aimInfo, &mesh);
          AIM_STATUS(aimInfo, status);

          status = aim_writeMeshes(aimInfo, Volume_Mesh, &mesh);
          AIM_STATUS(aimInfo, status);

          status = aim_freeMeshData(mesh.meshData);
          AIM_STATUS(aimInfo, status);
          AIM_FREE(mesh.meshData);

        }
        else
          AIM_STATUS(aimInfo, status);

/*@-immediatetrans@*/
        // Return the volume mesh reference
        val->nrow        = 1;
        val->vals.AIMptr = &pointwiseInstance->meshRef;
/*@+immediatetrans@*/
    } else {

        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

    }

cleanup:

    return status;
}


void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *pointwiseInstance;

#ifdef DEBUG
    printf(" pointwiseAIM/aimCleanup!\n");
#endif
    pointwiseInstance = (aimStorage *) instStore;
    status = destroy_aimStorage(pointwiseInstance);
    if (status != CAPS_SUCCESS)
      printf("Status = %d -- pointwiseAIM aimStorage cleanup!!!\n", status);

    EG_free(pointwiseInstance);
}
