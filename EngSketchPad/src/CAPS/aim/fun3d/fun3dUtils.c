// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>

#include "egads.h"        // Bring in egads utilss
#include "capsTypes.h"    // Bring in CAPS types
#include "aimUtil.h"      // Bring in AIM utils
#include "aimMesh.h"      // Bring in AIM meshing

#include "miscUtils.h"    // Bring in misc. utility functions
#include "meshUtils.h"    // Bring in meshing utility functions
#include "cfdTypes.h"     // Bring in cfd specific types
#include "cfdUtils.h"
#include "tecplotUtils.h" // Bring in tecplot utility functions
#include "fun3dUtils.h"   // Bring in fun3d utility header
#include "fun3dInputs.h"
#include "ugridWriter.h"


#ifdef WIN32
#define strcasecmp  stricmp
#define strncasecmp _strnicmp
#endif


// Extract the FEPOINT Tecoplot data from a FUN3D Aero-Loads file (connectivity is ignored) - dataMatrix = [numVariable][numDataPoint]
int fun3d_readAeroLoad(void *aimInfo, char *filename, int *numVariable, char **variableName[],
                       int *numDataPoint, double ***dataMatrix)
{

    int status = CAPS_SUCCESS; // Function return
    int i, j; // Indexing

    size_t linecap = 0;

    char *line = NULL; // Temporary line holder
    char *tempStr = NULL, *tempStr2 = NULL; // Temporary strings
    int stringLen = 0; // Length of string holder

    FILE *fp = NULL; // File pointer

    // Open file
    fp = aim_fopen(aimInfo, filename, "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file: %s\n", filename);
        return CAPS_IOERR;
    }

    printf("Reading FUN3D AeroLoad File - %s!!!!!!\n", filename);

    *numVariable = 0;
    *numDataPoint = 0;
    // Loop through file line by line until we have determined how many variables and data points there are
    while (*numVariable == 0 || *numDataPoint == 0) {

        // Get line from file
        status = getline(&line, &linecap, fp);
        if ((status < 0) || (line == NULL)) break;

        // Get variable list if available in file line
        if (strncmp("variables=", line, strlen("variables=")) == 0) {

            // Pull out substring at first occurrence of "
            tempStr = strstr(line, "\"");

            // Create a temperory string of the variables in a the folling format - ["a","ae"]
            stringLen  = strlen(tempStr)-1 + 2;

            tempStr2 = (char *) EG_alloc((stringLen +1)*sizeof(char *));
            if (tempStr2 == NULL) {
                fclose(fp);
                if (line != NULL) EG_free(line);
                return EGADS_MALLOC;
            }

            tempStr2[0] = '[';
            strncpy(tempStr2+1, tempStr, strlen(tempStr)-1);
            tempStr2[stringLen-1] = ']';
            tempStr2[stringLen] = '\0';

            // Sort string into an array of strings
            status =  string_toStringDynamicArray(tempStr2, numVariable, variableName);

            if (tempStr2 != NULL) EG_free(tempStr2);
            tempStr2 = NULL;

            if (status != CAPS_SUCCESS) goto cleanup;

            // Print out list of variables found in load file
            printf("Variables found in file %s:\n",filename);
            for (i = 0; i < *numVariable; i++) printf("Variable %d = %s\n", i, (*variableName)[i]);
        }

        // Get the number of data points in file if available in file line
        if (strncmp("zone t=", line, strlen("zone t=")) == 0) {

            // Pull out substring at first occurrence of i=
            tempStr = strstr(line, "i=");

            // Retrieve the i= value
            sscanf(&tempStr[2], "%d", numDataPoint);

            // Print out the number of data points found in load file
            printf("Number of data points = %d, in file %s\n", *numDataPoint, filename);
        }

    }

    if (*numVariable != 0 && *numDataPoint != 0) {

        // Allocate dataMatrix array
        AIM_FREE(*dataMatrix);

        AIM_ALLOC(*dataMatrix, (*numVariable), double *, aimInfo, status);
        for (i = 0; i < *numVariable; i++) (*dataMatrix)[i] = NULL;

        for (i = 0; i < *numVariable; i++) {
            AIM_ALLOC((*dataMatrix)[i], (*numDataPoint), double, aimInfo, status);
        }

        // Loop through the file and fill up the data matrix
        for (j = 0; j < *numDataPoint; j++) {
            for (i = 0; i < *numVariable; i++) {
                fscanf(fp, "%lf", &(*dataMatrix)[i][j]);
            }
        }

        // Output the first row of the dataMatrix
        //for (i = 0; i < *numVariable; i++) printf("Variable %d - %.6f\n", i, (*dataMatrix)[i][0]);

    } else {

        printf("No data values extracted from file - %s",filename);
        status = CAPS_BADVALUE;
    }

cleanup:
    if (fp != NULL) fclose(fp);

    if (status != CAPS_SUCCESS) {
      if (*dataMatrix != NULL) {
        for (j = 0; j < *numVariable; j++) {
          AIM_FREE((*dataMatrix)[j]);
        }
        AIM_FREE((*dataMatrix));
      }
    }

    EG_free(line);
    return status;
}


static int
fun3d_read2DBinaryUgrid(void *aimInfo, FILE *fp, meshStruct *surfaceMesh)
{
  int    status = CAPS_SUCCESS;

  int    numNode, numLine, numTriangle, numQuadrilateral;
  int    numTetrahedral, numPyramid, numPrism, numHexahedral;
  int    i, elementIndex, numPoint, bcID;
  double *coords = NULL;

  /* we get a binary UGRID file */
  status = fread(&numNode,          sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fread(&numTriangle,      sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fread(&numQuadrilateral, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fread(&numTetrahedral,   sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fread(&numPyramid,       sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fread(&numPrism,         sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fread(&numHexahedral,    sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  if ( numTetrahedral   +
       numPyramid       +
       numPrism         +
       numHexahedral  > 0) {
    AIM_ERROR(aimInfo, "Expecting a 2D ugrid file!!!");
    status = CAPS_IOERR;
    goto cleanup;
  }

  /*
  printf("\n Header from UGRID file: %d  %d %d  %d %d %d %d\n", numNode,
         numTriangle, numQuadrilateral, numTetrahedral, numPyramid, numPrism,
         numHexahedral);
   */

  AIM_ALLOC(coords, 3*numNode, double, aimInfo, status);

  /* read all of the vertices */
  status = fread(coords, sizeof(double), 3*numNode, fp);
  if (status != 3*numNode) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  surfaceMesh->analysisType = UnknownMeshAnalysis;

  // Set that this is a volume mesh
  surfaceMesh->meshType = Surface2DMesh;

  // Numbers
  surfaceMesh->numNode = numNode;
  surfaceMesh->numElement = numTriangle + numQuadrilateral;

  surfaceMesh->meshQuickRef.useStartIndex = (int) true;

  surfaceMesh->meshQuickRef.numTriangle      = numTriangle;
  surfaceMesh->meshQuickRef.numQuadrilateral = numQuadrilateral;

  // Nodes - allocate
  AIM_ALLOC(surfaceMesh->node, surfaceMesh->numNode, meshNodeStruct, aimInfo, status);

  // Initialize
  for (i = 0; i < surfaceMesh->numNode; i++) {
      status = initiate_meshNodeStruct(&surfaceMesh->node[i],
                                        surfaceMesh->analysisType);
      AIM_STATUS(aimInfo, status);
  }

  // Nodes - set
  for (i = 0; i < surfaceMesh->numNode; i++) {

    // Copy node data
    surfaceMesh->node[i].nodeID = i+1;

    surfaceMesh->node[i].xyz[0] = coords[3*i+0];
    surfaceMesh->node[i].xyz[1] = coords[3*i+1];
    surfaceMesh->node[i].xyz[2] = coords[3*i+2];
  }
  AIM_FREE(coords);

  // Elements - allocate
  AIM_ALLOC(surfaceMesh->element, surfaceMesh->numElement, meshElementStruct, aimInfo, status);

  // Initialize
  for (i = 0; i < surfaceMesh->numElement; i++ ) {
      status = initiate_meshElementStruct(&surfaceMesh->element[i],
                                           surfaceMesh->analysisType);
      AIM_STATUS(aimInfo, status);
  }

  // Start of element index
  elementIndex = 0;

  // Elements -Set triangles
  if (numTriangle > 0)
    surfaceMesh->meshQuickRef.startIndexTriangle = elementIndex;

  numPoint = mesh_numMeshConnectivity(Triangle);
  for (i = 0; i < numTriangle; i++) {

      surfaceMesh->element[elementIndex].elementType = Triangle;
      surfaceMesh->element[elementIndex].elementID   = elementIndex+1;

      status = mesh_allocMeshElementConnectivity(&surfaceMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(surfaceMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      elementIndex += 1;
  }

  // Elements -Set quadrilateral
  if (numQuadrilateral > 0)
    surfaceMesh->meshQuickRef.startIndexQuadrilateral = elementIndex;

  numPoint = mesh_numMeshConnectivity(Quadrilateral);
  for (i = 0; i < numQuadrilateral; i++) {

      surfaceMesh->element[elementIndex].elementType = Quadrilateral;
      surfaceMesh->element[elementIndex].elementID   = elementIndex+1;

      status = mesh_allocMeshElementConnectivity(&surfaceMesh->element[elementIndex]);
      if (status != CAPS_SUCCESS) goto cleanup;

      // read the element connectivity
      status = fread(surfaceMesh->element[elementIndex].connectivity,
                     sizeof(int), numPoint, fp);
      if (status != numPoint) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

      elementIndex += 1;
  }

  // skip face ID section of the file
  status = fseek(fp, (numTriangle + numQuadrilateral)*sizeof(int), SEEK_CUR);
  if (status != 0) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  // Get the number of Line elements
  status = fread(&numLine, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  // Elements - re-allocate with Line elements
  surfaceMesh->meshQuickRef.numLine = numLine;
  AIM_REALL(surfaceMesh->element, surfaceMesh->numElement+numLine, meshElementStruct, aimInfo, status);

  surfaceMesh->meshQuickRef.startIndexLine = elementIndex;

  // Initialize
  for (i = surfaceMesh->numElement; i < surfaceMesh->numElement+numLine; i++ ) {
      status = initiate_meshElementStruct(&surfaceMesh->element[i],
                                           surfaceMesh->analysisType);
      AIM_STATUS(aimInfo, status);
  }
  surfaceMesh->numElement += numLine;

  numPoint = mesh_numMeshConnectivity(Line);
  for (i = 0; i < numLine; i++) {

    surfaceMesh->element[elementIndex].elementType = Line;
    surfaceMesh->element[elementIndex].elementID   = elementIndex+1;

    status = mesh_allocMeshElementConnectivity(&surfaceMesh->element[elementIndex]);
    AIM_STATUS(aimInfo, status);

    // read the element connectivity
    status = fread(surfaceMesh->element[elementIndex].connectivity,
                   sizeof(int), numPoint, fp);
    if (status != numPoint) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
    status = fread(&bcID, sizeof(int), 1, fp);
    if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

    surfaceMesh->element[elementIndex].markerID = bcID;

    elementIndex += 1;
  }

  status = CAPS_SUCCESS;

  cleanup:
      if (status != CAPS_SUCCESS)
        printf("Premature exit in getUGRID status = %d\n", status);

      EG_free(coords); coords = NULL;

      return status;
}


// Create a 3D mesh for FUN3D from a 2D mesh
int fun3d_2DMesh(void *aimInfo,
                 aimMeshRef *meshRef,
                 const char *projectName,
                 mapAttrToIndexStruct *groupMap,
                 cfdBoundaryConditionStruct *bcProps)
{

    int status; // Function return status

    int i; // Indexing

    int faceBCIndex = -1, extrusionBCIndex = -1;

    double extrusion = -1.0; // Extrusion length

    // Flip coordinates
    int xMeshConstant = (int) true, yMeshConstant = (int) true, zMeshConstant= (int) true; // 2D mesh checks
    double tempCoord;

    meshStruct surfaceMesh;
    meshStruct volumeMesh;

    char filename[PATH_MAX];
//    int elementIndex;
    FILE *fp = NULL;

    status = initiate_meshStruct(&surfaceMesh);
    AIM_STATUS(aimInfo, status);

    status = initiate_meshStruct(&volumeMesh);
    AIM_STATUS(aimInfo, status);

    sprintf(filename, "%s%s", meshRef->fileName, MESHEXTENSION);

    fp = fopen(filename, "rb");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open file: %s\n", filename);
      status = CAPS_IOERR;
      goto cleanup;
    }

    status = fun3d_read2DBinaryUgrid(aimInfo, fp, &surfaceMesh);
    AIM_STATUS(aimInfo, status);

    // add boundary elements if they are missing
    if (surfaceMesh.meshQuickRef.numLine == 0) {
        status = mesh_addTess2Dbc(aimInfo, &surfaceMesh, groupMap);
        if (status != CAPS_SUCCESS) goto cleanup;
    }


    // Find the faceBCIndex for the symmetry plane
    for (i = 0; i < bcProps->numSurfaceProp-1; i++) {
      if (bcProps->surfaceProp[i].surfaceType == Symmetry) {
        faceBCIndex = bcProps->surfaceProp[i].bcID;
        break;
      }
    }

    if (faceBCIndex == -1) {
      // Add plane boundary condition
      AIM_REALL(bcProps->surfaceProp, bcProps->numSurfaceProp+1, cfdSurfaceStruct, aimInfo, status);
      bcProps->numSurfaceProp += 1;

      status = initiate_cfdSurfaceStruct(&bcProps->surfaceProp[bcProps->numSurfaceProp-1]);
      AIM_STATUS(aimInfo, status);

      bcProps->surfaceProp[bcProps->numSurfaceProp-1].surfaceType = Symmetry;
      bcProps->surfaceProp[bcProps->numSurfaceProp-1].symmetryPlane = 1;

      // Find largest index value for bcID and set it plus 1 to the new surfaceProp
      for (i = 0; i < bcProps->numSurfaceProp-1; i++) {
          if (bcProps->surfaceProp[i].bcID >= faceBCIndex) {
            faceBCIndex = bcProps->surfaceProp[i].bcID + 1;
          }
      }
      bcProps->surfaceProp[bcProps->numSurfaceProp-1].bcID = faceBCIndex;
    }


    // Add extruded plane boundary condition
    AIM_REALL(bcProps->surfaceProp, bcProps->numSurfaceProp+1, cfdSurfaceStruct, aimInfo, status);
    bcProps->numSurfaceProp += 1;

    status = initiate_cfdSurfaceStruct(&bcProps->surfaceProp[bcProps->numSurfaceProp-1]);
    AIM_STATUS(aimInfo, status);

    bcProps->surfaceProp[bcProps->numSurfaceProp-1].surfaceType = Symmetry;
    bcProps->surfaceProp[bcProps->numSurfaceProp-1].symmetryPlane = 2;

    // Find largest index value for bcID and set it plus 1 to the new surfaceProp
    for (i = 0; i < bcProps->numSurfaceProp-1; i++) {
        if (bcProps->surfaceProp[i].bcID >= extrusionBCIndex) {
          extrusionBCIndex = bcProps->surfaceProp[i].bcID + 1;
        }
    }
    bcProps->surfaceProp[bcProps->numSurfaceProp-1].bcID = extrusionBCIndex;


    // Set the symmetry index for all Tri/Quad
    for (i = 0; i < surfaceMesh.numElement; i++) {

        if (surfaceMesh.element[i].elementType != Triangle &&
            surfaceMesh.element[i].elementType != Quadrilateral) {
            continue;
        }

        surfaceMesh.element[i].markerID = faceBCIndex;
    }

#ifdef I_DONT_THINK_WE_NEED_THIS
    // Determine a suitable boundary index of the extruded plane
    *extrusionBCIndex = faceBCIndex;
    for (i = 0; i < surfaceMesh.meshQuickRef.numLine; i++) {

        if (surfaceMesh.meshQuickRef.startIndexLine >= 0) {
            elementIndex = surfaceMesh.meshQuickRef.startIndexLine + i;
        } else {
            elementIndex = surfaceMesh.meshQuickRef.listIndexLine[i];
        }

        marker = surfaceMesh.element[elementIndex].markerID;

        if (marker > *extrusionBCIndex)  {
            *extrusionBCIndex = marker;
        }
    }
    *extrusionBCIndex += 1;
#endif

    // Check to make sure the face is on the y = 0 plane
    for (i = 0; i < surfaceMesh.numNode; i++) {

        if (surfaceMesh.node[i].xyz[1] != 0.0) {
            printf("\nSurface mesh is not on y = 0.0 plane, FUN3D could fail during execution for this 2D mesh!!!\n");
            break;
        }
    }

    // Constant x?
    for (i = 0; i < surfaceMesh.numNode; i++) {
        if ((surfaceMesh.node[i].xyz[0] - surfaceMesh.node[0].xyz[0]) > 1E-7) {
            xMeshConstant = (int) false;
            break;
        }
    }

    // Constant y?
    for (i = 0; i < surfaceMesh.numNode; i++) {
        if ((surfaceMesh.node[i].xyz[1] - surfaceMesh.node[0].xyz[1] ) > 1E-7) {
            yMeshConstant = (int) false;
            break;
        }
    }

    // Constant z?
    for (i = 0; i < surfaceMesh.numNode; i++) {
        if ((surfaceMesh.node[i].xyz[2] - surfaceMesh.node[0].xyz[2]) > 1E-7) {
            zMeshConstant = (int) false;
            break;
        }
    }

    if (yMeshConstant != (int) true) {
        printf("FUN3D expects 2D meshes be in the x-z plane... attempting to rotate mesh!\n");

        if (xMeshConstant == (int) true && zMeshConstant == (int) false) {
            printf("Swapping y and x coordinates!\n");
            for (i = 0; i < surfaceMesh.numNode; i++) {
                tempCoord = surfaceMesh.node[i].xyz[0];
                surfaceMesh.node[i].xyz[0] = surfaceMesh.node[i].xyz[1];
                surfaceMesh.node[i].xyz[1] = tempCoord;
            }

        } else if(xMeshConstant == (int) false && zMeshConstant == (int) true) {

            printf("Swapping y and z coordinates!\n");
            for (i = 0; i < surfaceMesh.numNode; i++) {
                tempCoord = surfaceMesh.node[i].xyz[2];
                surfaceMesh.node[i].xyz[2] = surfaceMesh.node[i].xyz[1];
                surfaceMesh.node[i].xyz[1] = tempCoord;
            }

        } else {
            AIM_ERROR(aimInfo, "Unable to rotate mesh!\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }
    }

    status = extrude_SurfaceMesh(extrusion, extrusionBCIndex, &surfaceMesh, &volumeMesh);
    AIM_STATUS(aimInfo, status);

    strcpy(filename, projectName);

    // Write AFLR3
/*@-nullpass@*/
     status = mesh_writeAFLR3(aimInfo, filename,
                              0, // write binary file
                              &volumeMesh,
                              1.0);
/*@+nullpass@*/
     AIM_STATUS(aimInfo, status);

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) {
        printf("Error: Premature exit in fun3d_2DMesh status = %d\n", status);
    }

/*@-dependenttrans@*/
    if (fp != NULL) fclose(fp);
/*@+dependenttrans@*/

    // Destroy meshes
    (void) destroy_meshStruct(&surfaceMesh);
    (void) destroy_meshStruct(&volumeMesh);

    return status;
}


// Remove unused nodes from dataMatrix and update connectivity matrix
static int fun3d_removeUnused(void *aimInfo, int numVariable, int *numNode,  int used[],
                              double ***data, /*@null@*/ int *numConnect,
                              /*@null@*/ int *dataConnectMatrix)
{
    int i, j, k; //Indexing
    int status;

    int *usedNode=NULL; // Freeable

    double **dataMatrix;

    dataMatrix = *data;

    // Copy used node array
    usedNode = (int *) EG_alloc(*numNode*sizeof(int));
    if (usedNode == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }
    for (i = 0; i < *numNode; i++) usedNode[i] = used[i];

    // Remove unused nodes
    j = 0;
    for (i = 0; i< *numNode; i++ ) {
        if (usedNode[i] == (int) false || usedNode[i] < 0) continue;

        usedNode[i] = j+1; // Set i-th node to essentially the node ID,  1-bias

        j +=1;
    }

    j = 0;
    for (i = 0; i< *numNode; i++ ) {
        if (usedNode[i] == (int) false || usedNode[i] < 0) continue;

        for (k = 0; k < numVariable; k++) {
            dataMatrix[k][j] = dataMatrix[k][i]; //Re-order dataMatrix - bubbling i'th index to j
        }

        j += 1;
    }

    *numNode = j; // New number of nodes

    // Redo connectivity
    if (dataConnectMatrix != NULL) {
        AIM_NOTNULL(numConnect, aimInfo, status);
        j = 0;
        for (i = 0; i < *numConnect; i++) {

            if (usedNode[dataConnectMatrix[4*i+ 0]-1] == (int) false) continue;
            if (usedNode[dataConnectMatrix[4*i+ 1]-1] == (int) false) continue;
            if (usedNode[dataConnectMatrix[4*i+ 2]-1] == (int) false) continue;
            if (usedNode[dataConnectMatrix[4*i+ 3]-1] == (int) false) continue;

            for (k = 0; k < 4; k++) {
               dataConnectMatrix[4*j+ k] = usedNode[dataConnectMatrix[4*i+ k]-1];
            }

            j += 1;
        }

        *numConnect = j; // New number of elements
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in fun3d_removeUnused status = %d\n",
               status);

    if (usedNode != NULL) EG_free(usedNode);
    return status;
}


// Write FUN3D data transfer files
int fun3d_dataTransfer(void *aimInfo,
                       const char *projectName,
                       mapAttrToIndexStruct *groupMap,
                       cfdBoundaryConditionStruct bcProps,
                       aimMeshRef *meshRef,
                       /*@null@*/ cfdModalAeroelasticStruct *eigenVector)
{

    /*! \page dataTransferFUN3D FUN3D Data Transfer
     *
     * \section dataToFUN3D Data transfer to FUN3D (FieldIn)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements (as from a structural solver)
     *   and updates FUN3D's surface mesh; a new [project_name]_body1.dat file is written out which may
     *   be loaded into FUN3D to update the surface mesh/move the volume mesh using the FUN3D command line option
     *   -\-read_surface_from_file
     * </ul>
     *
     * <ul>
     *  <li> <B>"EigenVector_#"</B> </li> <br>
     *   Retrieves modal eigen-vectors from a structural solver, where "#" should be replaced by the
     *   corresponding mode number for the eigen-vector (eg. EigenVector_3 would correspond to the third mode,
     *   while EigenVector_6 would be the sixth mode) . A [project_name]_body1_mode#.dat file is written
     *   out for each mode.
     * </ul>
     *
     */

    int status; // Function return status
    int i, j, ibound, ibody, iface, iglobal, eigenIndex; // Indexing

    int stringLength = 0;

    char *filename = NULL;

    // Discrete data transfer variables
    capsDiscr *discr;
    char **boundName = NULL;
    int numBoundName = 0;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int dataTransferRank;
    double *dataTransferData;
    char   *units;

    int state, nGlobal, *globalOffset=NULL, nFace;
    ego body, *faces=NULL;

    int alen, ntri, atype, itri, ielem;
    const double *face_xyz, *face_uv, *reals;
    const int *face_ptype, *face_pindex, *face_tris, *face_tric, *nquad=NULL;
    const char *string, *groupName = NULL;

    // Variables used in global node mapping
    int ptype, pindex;
    double xyz[3];

    // Data transfer Out variables
    const char *dataOutName[] = {"x","y","z", "id", "dx", "dy", "dz"};
    const int dataOutFormat[] = {Double, Double, Double, Integer, Double, Double, Double};

    double **dataOutMatrix = NULL;
    int *dataConnectMatrix = NULL;

    int numOutVariable = 7;
    int numOutDataPoint = 0;
    int numOutDataConnect = 0;

    const char fileExtBody[] = "_body1";
    const char fileExt[] = ".dat";
    const char fileExtMode[] = "_mode";

    int foundDisplacement = (int) false, foundEigenVector = (int) false;

    int marker;

    int numUsedNode = 0, numUsedConnectivity = 0, usedElems;
    int *usedNode = NULL;

    status = aim_getBounds(aimInfo, &numBoundName, &boundName);
    AIM_STATUS(aimInfo, status);

    foundDisplacement = foundEigenVector = (int) false;
    for (ibound = 0; ibound < numBoundName; ibound++) {
      AIM_NOTNULL(boundName, aimInfo, status);

      status = aim_getDiscr(aimInfo, boundName[ibound], &discr);
      if (status != CAPS_SUCCESS) continue;

      status = aim_getDataSet(discr,
                              "Displacement",
                              &dataTransferMethod,
                              &numDataTransferPoint,
                              &dataTransferRank,
                              &dataTransferData,
                              &units);

      if (status == CAPS_SUCCESS) { // If we do have data ready is the rank correct

        foundDisplacement = (int) true;

        if (dataTransferRank != 3) {
          AIM_ERROR(aimInfo, "Displacement transfer data found however rank is %d not 3!!!!\n", dataTransferRank);
          status = CAPS_BADRANK;
          goto cleanup;
        }
        break;
      }

      if (eigenVector != NULL) {
        for (eigenIndex = 0; eigenIndex < eigenVector->numEigenValue; eigenIndex++) {

          status = aim_getDataSet(discr,
                                  eigenVector->eigenValue[eigenIndex].name,
                                  &dataTransferMethod,
                                  &numDataTransferPoint,
                                  &dataTransferRank,
                                  &dataTransferData,
                                  &units);

          if (status == CAPS_SUCCESS) { // If we do have data ready is the rank correct

            foundEigenVector = (int) true;

            if (dataTransferRank != 3) {
              AIM_ERROR(aimInfo, "EigenVector transfer data found however rank is %d not 3!!!!\n", dataTransferRank);
              status = CAPS_BADRANK;
              goto cleanup;
            }
            break;
          }
        } // Loop through EigenValues

        if (foundEigenVector == (int) true) break;

      } // If eigen-vectors provided
    } // Loop through transfer names

    if (foundDisplacement != (int) true && foundEigenVector != (int) true) {
      printf("No recognized data transfer names found!\n");
      status = CAPS_NOTFOUND;
      goto cleanup;
    }

    // Ok looks like we have displacements/EigenVectors to get so lets continue
    printf("Writing FUN3D data transfer files\n");

    // Allocate data arrays that are going to be output
    AIM_ALLOC(dataOutMatrix, numOutVariable, double*, aimInfo, status);
    for (i = 0; i < numOutVariable; i++) dataOutMatrix[i] = NULL;

    AIM_ALLOC(globalOffset, meshRef->nmap+1, int, aimInfo, status);
    numOutDataPoint = 0;
    ielem = 0;
    globalOffset[0] = 0;
    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &nGlobal);
      AIM_STATUS(aimInfo, status);

      // re-allocate data arrays
      for (j = 0; j < numOutVariable; j++) {
        AIM_REALL(dataOutMatrix[j], numOutDataPoint + nGlobal, double, aimInfo, status);
      }
      AIM_REALL(usedNode, numOutDataPoint + nGlobal, int, aimInfo, status);
      for (j = globalOffset[0]; j < numOutDataPoint + nGlobal; j++) usedNode[j] = (int) false;


      for (iglobal = 0; iglobal < nGlobal; iglobal++) {
        status = EG_getGlobal(meshRef->maps[i].tess,
                              iglobal+1, &ptype, &pindex, xyz);
        AIM_STATUS(aimInfo, status);

        // First just set the Coordinates
        dataOutMatrix[0][globalOffset[i]+iglobal] = xyz[0];
        dataOutMatrix[1][globalOffset[i]+iglobal] = xyz[1];
        dataOutMatrix[2][globalOffset[i]+iglobal] = xyz[2];

        // Volume mesh node ID
        dataOutMatrix[3][globalOffset[i]+iglobal] = meshRef->maps[i].map[iglobal];

        // Delta displacements
        dataOutMatrix[4][i] = 0;
        dataOutMatrix[5][i] = 0;
        dataOutMatrix[6][i] = 0;
      }

      // check if the tessellation has a mixture of quad and tess
      status = EG_attributeRet(meshRef->maps[i].tess, ".mixed",
                               &atype, &alen, &nquad, &reals, &string);
      if (status != EGADS_SUCCESS &&
          status != EGADS_NOTFOUND) AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(body, NULL, FACE, &nFace, &faces);
      AIM_STATUS(aimInfo, status);
      for (iface = 0; iface < nFace; iface++) {
        // get the face tessellation
        status = EG_getTessFace(meshRef->maps[i].tess, iface+1, &alen, &face_xyz, &face_uv,
                                &face_ptype, &face_pindex, &ntri, &face_tris, &face_tric);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(faces, aimInfo, status);

        status = retrieve_CAPSGroupAttr(faces[iface], &groupName);
        if (status == EGADS_SUCCESS) {
          AIM_NOTNULL(groupName, aimInfo, status);
          status = get_mapAttrToIndexIndex(groupMap, groupName, &marker);
          if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
            goto cleanup;
          }
        } else {
          AIM_ERROR(aimInfo, "No capsGroup on face %d", iface+1);
          print_AllAttr(aimInfo, faces[iface]);
          goto cleanup;
        }

        //  To keep with the moving_bodying input we will assume used nodes are all inviscid and viscous surfaces instead
        //                        usedNode[k] = (int) true;
        usedElems = (int) false;
        for (j = 0; j < bcProps.numSurfaceProp; j++) {
          if (marker != bcProps.surfaceProp[j].bcID) continue;

          if (bcProps.surfaceProp[j].surfaceType == Viscous ||
              bcProps.surfaceProp[j].surfaceType == Inviscid) {
              usedElems = (int) true;
          }
          break;
        }

        if (nquad == NULL) { // all triangles

          // re-allocate data arrays
          AIM_REALL(dataConnectMatrix, 4*(numOutDataConnect + ntri), int, aimInfo, status);

          for (itri = 0; itri < ntri; itri++, ielem++) {
            for (j = 0; j < 3; j++) {
              status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+j], &iglobal);
              AIM_STATUS(aimInfo, status);
              dataConnectMatrix[4*ielem+j] = globalOffset[i] + iglobal;
              usedNode[globalOffset[i]+iglobal-1] = usedElems;
            }
            // repeat the last node for triangles
            dataConnectMatrix[4*ielem+3] = dataConnectMatrix[4*ielem+2];
          }

          numOutDataConnect += ntri;

        } else { // mixture of tri and quad elements

          // re-allocate data arrays
          AIM_REALL(dataConnectMatrix, 4*(numOutDataConnect + ntri-nquad[iface]), int, aimInfo, status);

          // process triangles
          for (itri = 0; itri < ntri-2*nquad[iface]; itri++, ielem++) {
            for (j = 0; j < 3; j++) {
              status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+j], &iglobal);
              AIM_STATUS(aimInfo, status);
              dataConnectMatrix[4*ielem+j] = globalOffset[i] + iglobal;
              usedNode[globalOffset[i]+iglobal-1] = usedElems;
            }
            // repeat the last node for triangle
            dataConnectMatrix[4*ielem+3] = dataConnectMatrix[4*ielem+2];
          }
          // process quads
          for (; itri < ntri; itri++, ielem++) {
            for (j = 0; j < 3; j++) {
              status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+j], &iglobal);
              AIM_STATUS(aimInfo, status);
              dataConnectMatrix[4*ielem+j] = globalOffset[i] + iglobal;
              usedNode[globalOffset[i]+iglobal-1] = usedElems;
            }

            // add the last node from the 2nd triangle to make the quad
            itri++;
            status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+2], &iglobal);
            AIM_STATUS(aimInfo, status);
            dataConnectMatrix[4*ielem+3] = globalOffset[i] + iglobal;
            usedNode[globalOffset[i]+iglobal-1] = usedElems;
          }

          numOutDataConnect += ntri-nquad[iface];
        }

      }
      AIM_FREE(faces);

      numOutDataPoint += nGlobal;
      globalOffset[i+1] = globalOffset[i] + nGlobal;
    }

    // Re-loop through transfers - if we are doing displacements
    if (foundDisplacement == (int) true) {

      for (ibound = 0; ibound < numBoundName; ibound++) {
        AIM_NOTNULL(boundName, aimInfo, status);

        status = aim_getDiscr(aimInfo, boundName[ibound], &discr);
        if (status != CAPS_SUCCESS) continue;

        status = aim_getDataSet(discr,
                                "Displacement",
                                &dataTransferMethod,
                                &numDataTransferPoint,
                                &dataTransferRank,
                                &dataTransferData,
                                &units);
        if (status != CAPS_SUCCESS) continue; // If no elements in this object skip to next transfer name

        if (numDataTransferPoint != discr->nPoints &&
            numDataTransferPoint > 1) {
          AIM_ERROR(aimInfo, "Developer error!! %d != %d", numDataTransferPoint, discr->nPoints);
          status = CAPS_MISMATCH;
          goto cleanup;
        }

        for (i = 0; i < discr->nPoints; i++) {

          ibody   = discr->tessGlobal[2*i+0];
          iglobal = discr->tessGlobal[2*i+1];

          status = EG_getGlobal(discr->bodys[ibody-1].tess,
                                iglobal, &ptype, &pindex, xyz);
          AIM_STATUS(aimInfo, status);

          // Find the disc tessellation in the original list of tessellations
          for (j = 0; j < meshRef->nmap; j++) {
            if (discr->bodys[ibody-1].tess == meshRef->maps[j].tess) {
              break;
            }
          }
          if (j == meshRef->nmap) {
            AIM_ERROR(aimInfo, "Could not find matching tessellation!");
            status = CAPS_MISMATCH;
            goto cleanup;
          }

          if (numDataTransferPoint == 1) {
            // A single point means this is an initialization phase

            // Apply delta displacements
            dataOutMatrix[0][globalOffset[j]+iglobal-1] += dataTransferData[0];
            dataOutMatrix[1][globalOffset[j]+iglobal-1] += dataTransferData[1];
            dataOutMatrix[2][globalOffset[j]+iglobal-1] += dataTransferData[2];

            // save delta displacements
            dataOutMatrix[4][globalOffset[j]+iglobal-1] = dataTransferData[0];
            dataOutMatrix[5][globalOffset[j]+iglobal-1] = dataTransferData[1];
            dataOutMatrix[6][globalOffset[j]+iglobal-1] = dataTransferData[2];

          } else {
            // Apply delta displacements
            dataOutMatrix[0][globalOffset[j]+iglobal-1] += dataTransferData[3*i+0];
            dataOutMatrix[1][globalOffset[j]+iglobal-1] += dataTransferData[3*i+1];
            dataOutMatrix[2][globalOffset[j]+iglobal-1] += dataTransferData[3*i+2];

            dataOutMatrix[4][globalOffset[j]+iglobal-1] = dataTransferData[3*i+0];
            dataOutMatrix[5][globalOffset[j]+iglobal-1] = dataTransferData[3*i+1];
            dataOutMatrix[6][globalOffset[j]+iglobal-1] = dataTransferData[3*i+2];
          }
        }
      } // End numBoundName loop

      // Remove unused nodes
      numUsedNode = numOutDataPoint;
      numUsedConnectivity = numOutDataConnect;
      AIM_NOTNULL(usedNode, aimInfo, status);
      status = fun3d_removeUnused(aimInfo, numOutVariable, &numUsedNode, usedNode,
                                  &dataOutMatrix, &numUsedConnectivity,
                                  dataConnectMatrix);
      AIM_STATUS(aimInfo, status);

      stringLength = strlen(projectName) + strlen(fileExtBody) + strlen(fileExt) + 1;
      AIM_ALLOC(filename, stringLength+1, char, aimInfo, status);

      strcpy(filename, projectName);
      strcat(filename, fileExtBody);
      strcat(filename, fileExt);
      filename[stringLength] = '\0';

      // Write out displacement in tecplot file
      /*@-nullpass@*/
      status = tecplot_writeFEPOINT(aimInfo, filename,
                                    "FUN3D AeroLoads",
                                    NULL,
                                    numOutVariable,
                                    (char **)dataOutName,
                                    numUsedNode, // numOutDataPoint,
                                    dataOutMatrix,
                                    dataOutFormat,
                                    numUsedConnectivity, //numOutDataConnect, // numConnectivity
                                    dataConnectMatrix, // connectivity matrix
                                    NULL); // Solution time
      /*@+nullpass@*/
      AIM_STATUS(aimInfo, status);
      AIM_FREE(filename);
    } // End if found displacements

    // Re-loop through transfers - if we are doing eigen-vectors
    if ((foundEigenVector == (int) true) && (eigenVector != NULL)) {

      for (eigenIndex = 0; eigenIndex < eigenVector->numEigenValue; eigenIndex++) {

        // Zero out the eigen-vectors each time we are writing out a new one
        for (i = 0; i < numOutDataPoint; i++ ) {

          // Delta eigen-vectors
          dataOutMatrix[4][i] = 0;
          dataOutMatrix[5][i] = 0;
          dataOutMatrix[6][i] = 0;
        }

        for (ibound = 0; ibound < numBoundName; ibound++) {
          AIM_NOTNULL(boundName, aimInfo, status);

          status = aim_getDiscr(aimInfo, boundName[ibound], &discr);
          if (status != CAPS_SUCCESS) continue;

          status = aim_getDataSet(discr,
                                  eigenVector->eigenValue[eigenIndex].name,
                                  &dataTransferMethod,
                                  &numDataTransferPoint,
                                  &dataTransferRank,
                                  &dataTransferData,
                                  &units);
          if (status != CAPS_SUCCESS) continue; // If no elements in this object skip to next transfer name

          if (numDataTransferPoint != discr->nPoints &&
              numDataTransferPoint > 1) {
            AIM_ERROR(aimInfo, "Developer error!! %d != %d", numDataTransferPoint, discr->nPoints);
            status = CAPS_MISMATCH;
            goto cleanup;
          }

          for (i = 0; i < discr->nPoints; i++) {

            ibody   = discr->tessGlobal[2*i+0];
            iglobal = discr->tessGlobal[2*i+1];

            // Find the disc tessellation in the original list of tessellations
            for (j = 0; j < meshRef->nmap; j++) {
              if (discr->bodys[ibody-1].tess == meshRef->maps[j].tess) {
                break;
              }
            }
            if (j == meshRef->nmap) {
              AIM_ERROR(aimInfo, "Could not find matching tessellation!");
              status = CAPS_MISMATCH;
              goto cleanup;
            }

            if (numDataTransferPoint == 1) {
              // A single point means this is an initialization phase

              // save Eigen-vector
              dataOutMatrix[4][globalOffset[j]+iglobal-1] = dataTransferData[0];
              dataOutMatrix[5][globalOffset[j]+iglobal-1] = dataTransferData[1];
              dataOutMatrix[6][globalOffset[j]+iglobal-1] = dataTransferData[2];

            } else {
              // save Eigen-vector
              dataOutMatrix[4][globalOffset[j]+iglobal-1] = dataTransferData[3*i+0];
              dataOutMatrix[5][globalOffset[j]+iglobal-1] = dataTransferData[3*i+1];
              dataOutMatrix[6][globalOffset[j]+iglobal-1] = dataTransferData[3*i+2];
            }
          }
        } // End dataTransferDiscreteObj loop

        // Remove unused nodes
        numUsedNode = numOutDataPoint;
        AIM_NOTNULL(usedNode, aimInfo, status);

        if (eigenIndex == 0) {
          numUsedConnectivity = numOutDataConnect;
          status = fun3d_removeUnused(aimInfo, numOutVariable, &numUsedNode,
                                      usedNode, &dataOutMatrix,
                                      &numUsedConnectivity, dataConnectMatrix);
          AIM_STATUS(aimInfo, status);
        } else {
          status = fun3d_removeUnused(aimInfo, numOutVariable, &numUsedNode,
                                      usedNode, &dataOutMatrix, NULL, NULL);
          AIM_STATUS(aimInfo, status);
        }

        stringLength = strlen(projectName) +
                       strlen(fileExtBody) +
                       strlen(fileExtMode) +
                       strlen(fileExt) + 5;

        AIM_ALLOC(filename, stringLength+1, char, aimInfo, status);

        sprintf(filename, "%s%s%s%d%s",
                projectName,
                fileExtBody,
                fileExtMode,  // Change modeNumber so it always starts at 1!
                eigenIndex+1, // eigenVector->eigenValue[eigenIndex].modeNumber,
                fileExt);

        // Write out eigen-vector in tecplot file
        /*@-nullpass@*/
        status = tecplot_writeFEPOINT(aimInfo,
                                      filename,
                                      "FUN3D Modal",
                                      NULL,
                                      numOutVariable,
                                      (char **)dataOutName,
                                      numUsedNode, //numOutDataPoint,
                                      dataOutMatrix,
                                      dataOutFormat,
                                      numUsedConnectivity, //numOutDataConnect, // numConnectivity
                                      dataConnectMatrix, // connectivity matrix
                                      NULL); // Solution time
        /*@+nullpass@*/
        AIM_STATUS(aimInfo, status);
        AIM_FREE(filename);

      } // End eigenvector names
    } // End if found eigenvectors

    status = CAPS_SUCCESS;

    // Clean-up
cleanup:

    if (status != CAPS_SUCCESS &&
        status != CAPS_NOTFOUND) printf("Error: Premature exit in fun3d_dataTransfer status = %d\n", status);

    if (dataOutMatrix != NULL) {
      for (i = 0; i < numOutVariable; i++) {
        AIM_FREE(dataOutMatrix[i]);
      }
    }

    AIM_FREE(faces);
    AIM_FREE(dataOutMatrix);
    AIM_FREE(dataConnectMatrix);

    AIM_FREE(filename);
    AIM_FREE(boundName);
    AIM_FREE(usedNode);

    return status;
}


// Write FUN3D fun3d.nml file
int fun3d_writeNML(void *aimInfo, capsValue *aimInputs, cfdBoundaryConditionStruct bcProps)
{

    int status; // Function return status

    int i; // Indexing

    FILE *fnml = NULL;
    char filename[PATH_MAX];
    char fileExt[] ="fun3d.nml";

    printf("Writing fun3d.nml\n");
    if (aimInputs[Design_Functional-1].nullVal == NotNull ||
        aimInputs[Design_SensFile-1].vals.integer == (int)true) {
#ifdef WIN32
        snprintf(filename, PATH_MAX, "Flow\\%s", fileExt);
#else
        snprintf(filename, PATH_MAX, "Flow/%s", fileExt);
#endif
    } else {
        strcpy(filename, fileExt);
    }

    fnml = aim_fopen(aimInfo, filename, "w");
    if (fnml == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file - %s\n", filename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // &project
    fprintf(fnml,"&project\n");
    fprintf(fnml," project_rootname = \"%s\"\n",
            aimInputs[Proj_Name-1].vals.string);
    fprintf(fnml,"/\n\n");

    // &raw_grid
    fprintf(fnml,"&raw_grid\n");
    //fprintf(fnml," grid_format = \"%s\"\n",
    //        aimInputs[Mesh_Format-1].vals.string);

//    if (aimInputs[Mesh_ASCII_Flag-1].vals.integer == (int) true) {
//        fprintf(fnml," data_format = \"ascii\"\n");
//    } else fprintf(fnml," data_format = \"stream\"\n");

    fprintf(fnml," grid_format = \"AFLR3\"\n");
    fprintf(fnml," data_format = \"stream\"\n");

    if (aimInputs[Two_Dimensional-1].vals.integer == (int) true) {
        fprintf(fnml," twod_mode = .true.\n");
      //fprintf(fnml," ignore_euler_number = .true.\n");
    }

    fprintf(fnml,"/\n\n");

    // &reference_physical_properties
    fprintf(fnml,"&reference_physical_properties\n");

    if (aimInputs[Mach-1].nullVal != IsNull) {
        fprintf(fnml," mach_number = %f\n", aimInputs[Mach-1].vals.real);
    }

    if (aimInputs[Re-1].nullVal != IsNull) {
        fprintf(fnml," reynolds_number = %f\n", aimInputs[Re-1].vals.real);
    }

    if (aimInputs[Alpha-1].nullVal != IsNull) {
        fprintf(fnml," angle_of_attack = %f\n", aimInputs[Alpha-1].vals.real);
    }

    if (aimInputs[Beta-1].nullVal != IsNull) {
        fprintf(fnml," angle_of_yaw = %f\n", aimInputs[Beta-1].vals.real);
    }

    if (aimInputs[Reference_Temperature-1].nullVal != IsNull) {
            fprintf(fnml," temperature = %f\n", aimInputs[Reference_Temperature-1].vals.real);

            if (aimInputs[Reference_Temperature-1].units != NULL) {
                fprintf(fnml," temperature_units = \'%s\'\n", aimInputs[Reference_Temperature-1].units);
            }
    }

    fprintf(fnml,"/\n\n");

    // &governing_equations
    fprintf(fnml,"&governing_equations\n");

    if (aimInputs[Viscoux-1].nullVal != IsNull) {
        fprintf(fnml," viscous_terms = \"%s\"\n", aimInputs[Viscoux-1].vals.string);
    }

    if (aimInputs[Equation_Type-1].nullVal != IsNull) {
        fprintf(fnml," eqn_type = \"%s\"\n", aimInputs[Equation_Type-1].vals.string);
    }

    fprintf(fnml,"/\n\n");

    // &nonlinear_solver_parameters
    fprintf(fnml,"&nonlinear_solver_parameters\n");

    if (aimInputs[Time_Accuracy-1].nullVal != IsNull) {
        fprintf(fnml," time_accuracy = \"%s\"\n",
                aimInputs[Time_Accuracy-1].vals.string);
    }

    if (aimInputs[Time_Step-1].nullVal != IsNull) {
        fprintf(fnml," time_step_nondim = %f\n", aimInputs[Time_Step-1].vals.real);
    }

    if (aimInputs[Num_Subiter-1].nullVal != IsNull) {
        fprintf(fnml," subiterations = %d\n",
                aimInputs[Num_Subiter-1].vals.integer);
    }

    if (aimInputs[Temporal_Error-1].nullVal != IsNull) {

        fprintf(fnml," temporal_err_control = .true.\n");
        fprintf(fnml," temporal_err_floor = %f\n",
                aimInputs[Temporal_Error-1].vals.real);

    }

    if (aimInputs[CFL_Schedule-1].nullVal != IsNull) {
        fprintf(fnml," schedule_cfl = %f %f\n",
                aimInputs[CFL_Schedule-1].vals.reals[0],
                aimInputs[CFL_Schedule-1].vals.reals[1]);
    }

    if (aimInputs[CFL_Schedule_Iter-1].nullVal != IsNull) {
        fprintf(fnml," schedule_iteration = %d %d\n",
                aimInputs[CFL_Schedule_Iter-1].vals.integers[0],
                aimInputs[CFL_Schedule_Iter-1].vals.integers[1]);
    }

    fprintf(fnml,"/\n\n");

    // &code_run_control
    fprintf(fnml,"&code_run_control\n");

    if (aimInputs[Num_Iter-1].nullVal != IsNull) {
        fprintf(fnml," steps = %d\n", aimInputs[Num_Iter-1].vals.integer);
    }

    if (aimInputs[Restart_Read-1].nullVal != IsNull) {
        fprintf(fnml," restart_read = '%s'\n",
                aimInputs[Restart_Read-1].vals.string);
    }


    fprintf(fnml,"/\n\n");

    //&force_moment_integ_properties
    fprintf(fnml,"&force_moment_integ_properties\n");

    if (aimInputs[Reference_Area-1].nullVal != IsNull) {
        fprintf(fnml," area_reference = %f\n",
                aimInputs[Reference_Area-1].vals.real);
    }

    if (aimInputs[Moment_Length-1].nullVal != IsNull) {
        fprintf(fnml," x_moment_length = %f\n",
                aimInputs[Moment_Length-1].vals.reals[0]);

        fprintf(fnml," y_moment_length = %f\n",
                aimInputs[Moment_Length-1].vals.reals[1]);
    }

    if (aimInputs[Moment_Center-1].nullVal != IsNull) {
        fprintf(fnml," x_moment_center = %f\n",
                aimInputs[Moment_Center-1].vals.reals[0]);

        fprintf(fnml," y_moment_center = %f\n",
                aimInputs[Moment_Center-1].vals.reals[1]);

        fprintf(fnml," z_moment_center = %f\n",
                aimInputs[Moment_Center-1].vals.reals[2]);
    }

    fprintf(fnml,"/\n\n");

    //&boundary_conditions
    fprintf(fnml,"&boundary_conditions\n");

    // Loop through boundary conditions
    for (i = 0; i < bcProps.numSurfaceProp ; i++) {

        // Temperature
        if (bcProps.surfaceProp[i].wallTemperatureFlag == (int) true) {
            fprintf(fnml," wall_temperature(%d) = %f\n",bcProps.surfaceProp[i].bcID,
                    bcProps.surfaceProp[i].wallTemperature);
            fprintf(fnml," wall_temp_flag(%d) = .true.\n",bcProps.surfaceProp[i].bcID);
        }

        // Total pressure and temperature
        if (bcProps.surfaceProp[i].surfaceType == SubsonicInflow) {

            fprintf(fnml, " total_pressure_ratio(%d) = %f\n", bcProps.surfaceProp[i].bcID,
                    bcProps.surfaceProp[i].totalPressure);

            fprintf(fnml, " total_temperature_ratio(%d) = %f\n", bcProps.surfaceProp[i].bcID,
                    bcProps.surfaceProp[i].totalTemperature);
        }

        // Static pressure
        if (bcProps.surfaceProp[i].surfaceType == BackPressure ||
                bcProps.surfaceProp[i].surfaceType == SubsonicOutflow) {

            fprintf(fnml, " static_pressure_ratio(%d) = %f\n", bcProps.surfaceProp[i].bcID,
                    bcProps.surfaceProp[i].staticPressure);
        }

        // Mach number
        if (bcProps.surfaceProp[i].surfaceType == MachOutflow ||
                bcProps.surfaceProp[i].surfaceType == MassflowOut) {

            fprintf(fnml, " mach_bc(%d) = %f\n", bcProps.surfaceProp[i].bcID,
                    bcProps.surfaceProp[i].machNumber);
        }

        // Massflow
        if (bcProps.surfaceProp[i].surfaceType == MassflowIn ||
                bcProps.surfaceProp[i].surfaceType == MassflowOut) {

            fprintf(fnml, " massflow(%d) = %f\n", bcProps.surfaceProp[i].bcID,
                    bcProps.surfaceProp[i].massflow);
        }

        // Fixed inflow and outflow
        /*if (bcProps.surfaceProp[i].surfaceType == FixedInflow ||
            bcProps.surfaceProp[i].surfaceType == FixedOutflow) {

            fprintf(fnml, " qset(%d,1) = %f\n", bcProps.surfaceProp[i].bcID,
                                                bcProps.surfaceProp[i].staticDensity);

            fprintf(fnml, " qset(%d,2) = %f\n", bcProps.surfaceProp[i].bcID,
                                                bcProps.surfaceProp[i].uVelocity);

            fprintf(fnml, " qset(%d,3) = %f\n", bcProps.surfaceProp[i].bcID,
                                                bcProps.surfaceProp[i].vVelocity);

            fprintf(fnml, " qset(%d,4) = %f\n", bcProps.surfaceProp[i].bcID,
                                                bcProps.surfaceProp[i].wVelocity);

            fprintf(fnml, " qset(%d,5) = %f\n", bcProps.surfaceProp[i].bcID,
                                                bcProps.surfaceProp[i].staticDensity);
        }*/
    }

    fprintf(fnml,"/\n\n");

    // &noninertial_reference_frame
    fprintf(fnml,"&noninertial_reference_frame\n");


    if (aimInputs[NonInertial_Rotation_Rate-1].nullVal   != IsNull ||
        aimInputs[NonInertial_Rotation_Center-1].nullVal != IsNull) {

        fprintf(fnml," noninertial = .true.\n");
    }

    if (aimInputs[NonInertial_Rotation_Center-1].nullVal != IsNull) {
        fprintf(fnml," rotation_center_x = %f\n",
                aimInputs[NonInertial_Rotation_Center-1].vals.reals[0]);
        fprintf(fnml," rotation_center_y = %f\n",
                aimInputs[NonInertial_Rotation_Center-1].vals.reals[1]);
        fprintf(fnml," rotation_center_z = %f\n",
                aimInputs[NonInertial_Rotation_Center-1].vals.reals[2]);
    }

    if (aimInputs[NonInertial_Rotation_Rate-1].nullVal != IsNull) {
        fprintf(fnml," rotation_rate_x = %f\n",
                aimInputs[NonInertial_Rotation_Rate-1].vals.reals[0]);
        fprintf(fnml," rotation_rate_y = %f\n",
                aimInputs[NonInertial_Rotation_Rate-1].vals.reals[1]);
        fprintf(fnml," rotation_rate_z = %f\n",
                aimInputs[NonInertial_Rotation_Rate-1].vals.reals[2]);

    }

    fprintf(fnml,"/\n\n");

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in fun3d_writeNML status = %d\n", status);

    if (fnml != NULL) fclose(fnml);

    return status;
}


// Write FUN3D movingbody.input file
int fun3d_writeMovingBody(void *aimInfo, double fun3dVersion, cfdBoundaryConditionStruct bcProps,
                          cfdModalAeroelasticStruct *modalAeroelastic)
{

    int status; // Function return status

    int i, eigenIndex; // Indexing
    int counter = 0;
    FILE *fp = NULL;
    char *filename = NULL;
    char fileExt[] ="moving_body.input";

    int bodyIndex = 1;
    int stringLength;

    printf("Writing moving_body.input");

    stringLength = strlen(fileExt) + 1;

    AIM_ALLOC(filename,stringLength +1, char, aimInfo, status);

    strcpy(filename, fileExt);
    filename[stringLength] = '\0';

    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file - %s\n", filename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    // &body_definitions
    fprintf(fp,"&body_definitions\n");

    fprintf(fp," n_moving_bodies = %d\n", bodyIndex);

    counter = 0;
    for (i = 0; i < bcProps.numSurfaceProp; i++) {

        if (bcProps.surfaceProp[i].surfaceType == Viscous ||
                bcProps.surfaceProp[i].surfaceType == Inviscid) {

            fprintf(fp," defining_bndry(%d,%d) = %d\n", counter+1,
                    bodyIndex,
                    bcProps.surfaceProp[i].bcID);

            counter += 1;
        }
    }

    fprintf(fp," n_defining_bndry(%d) = %d\n", bodyIndex, counter);

    fprintf(fp," motion_driver(%d) = ", bodyIndex);
    if (modalAeroelastic != NULL) {
        fprintf(fp,"\"aeroelastic\"\n");
    }

    fprintf(fp," mesh_movement(%d) = ", bodyIndex);
    if (modalAeroelastic != NULL) {
        fprintf(fp,"\"deform\"\n");
    }

    fprintf(fp,"/\n\n");

    if (modalAeroelastic != NULL) {

        // &aeroelastic_modal_data
        fprintf(fp,"&aeroelastic_modal_data\n");

        fprintf(fp," nmode(%d) = %d\n", bodyIndex, modalAeroelastic->numEigenValue);

        if (fun3dVersion < 13.1) {
            fprintf(fp," uinf(%d) = %f\n", bodyIndex, modalAeroelastic->freestreamVelocity);
            fprintf(fp," qinf(%d) = %f\n", bodyIndex, modalAeroelastic->freestreamDynamicPressure);
            fprintf(fp," grefl(%d) = %f\n", bodyIndex, modalAeroelastic->lengthScaling);
        } else {
            fprintf(fp," uinf = %f\n", modalAeroelastic->freestreamVelocity);
            fprintf(fp," qinf = %f\n", modalAeroelastic->freestreamDynamicPressure);
            fprintf(fp," grefl = %f\n",modalAeroelastic->lengthScaling);
        }

        fprintf(fp, "\n");


        for (i = 0; i < modalAeroelastic->numEigenValue; i++) {

            eigenIndex = i + 1; // Change mode number so that it always starts at 1
            // modalAeroelastic->eigenValue[i].modeNumber

            fprintf(fp, " ! Mode %d of %d (structural mode %d)\n", eigenIndex,
                    modalAeroelastic->numEigenValue,
                    modalAeroelastic->eigenValue[i].modeNumber);

            fprintf(fp," freq(%d,%d) = %f\n", eigenIndex, bodyIndex, modalAeroelastic->eigenValue[i].frequency);
            fprintf(fp," damp(%d,%d) = %f\n", eigenIndex, bodyIndex, modalAeroelastic->eigenValue[i].damping);

            fprintf(fp," gmass(%d,%d) = %f\n"  , eigenIndex, bodyIndex, modalAeroelastic->eigenValue[i].generalMass);
            fprintf(fp," gdisp0(%d,%d) = %f\n" , eigenIndex, bodyIndex, modalAeroelastic->eigenValue[i].generalDisplacement);
            fprintf(fp," gvel0(%d,%d) = %f\n"  , eigenIndex, bodyIndex, modalAeroelastic->eigenValue[i].generalVelocity);
            fprintf(fp," gforce0(%d,%d) = %f\n", eigenIndex, bodyIndex, modalAeroelastic->eigenValue[i].generalForce);

            fprintf(fp, "\n");
        }
        fprintf(fp,"/\n\n");
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in fun3d_writeMovingBody status = %d\n",
               status);

    if (fp != NULL) fclose(fp);

    if (filename != NULL) EG_free(filename);

    return status;
}


// Write FUN3D parameterization/sensitivity file
int  fun3d_writeParameterization(void *aimInfo,
                                 int numDesignVariable,
                                 cfdDesignVariableStruct designVariable[],
                                 int writeSensitivity,
                                 aimMeshRef *meshRef)
{

    int status; // Function return status

    int i, j, k, m, row, col; // Indexing

    int stringLength = 7;

    // Data transfer Out variables
    char **dataOutName= NULL;

    double ***dataOutMatrix = NULL;
    int *dataOutFormat = NULL;
    int *dataConnectMatrix = NULL;

    int numOutVariable = 4; // x, y, z, id, ... + 3* active GeomIn
    int numOutDataPoint = 0;
    int numOutDataConnect = 0;

    // Variables used in global node mapping
    int ptype, pindex;
    double xyz[3];

    const char *geomInName;
    int numPoint;
    double *dxyz = NULL;

    int index;
    int iface, iglobal;
    int state, nFace;
    ego body, *faces=NULL;

    int alen, ntri, atype, itri, ielem;
    const double *face_xyz, *face_uv, *reals;
    const int *face_ptype, *face_pindex, *face_tris, *face_tric, *nquad=NULL;
    const char *string;

    char message[100];
    char filePre[] = "model.tec.";
    char fileExt[] = ".sd1";
    char *filename = NULL;
    char folder[]  = "Rubberize";
    char zoneTitle[100];

    capsValue *geomInVal;
    int *geomSelect = NULL;

    AIM_ALLOC(geomSelect, numDesignVariable, int, aimInfo, status);

    // Determine number of geometry input variables

    for (i = 0; i < numDesignVariable; i++) {
        geomSelect[i] = (int) false;

        index = aim_getIndex(aimInfo, designVariable[i].name, GEOMETRYIN);
        if (index == CAPS_NOTFOUND) continue;
        if (index < CAPS_SUCCESS ) {
          status = index;
          AIM_STATUS(aimInfo, status);
        }

        if(aim_getGeomInType(aimInfo, index) != 0) {
            AIM_ERROR(aimInfo, "GeometryIn value %s is a configuration parameter and not a valid design parameter - can't get sensitivity\n",
                      designVariable[i].name);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        // Fun3D always requires rubberize files, but don't compute sensitivities if not needed
        geomSelect[i] = writeSensitivity;

        if (writeSensitivity == (int)true) {
          status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
          AIM_STATUS(aimInfo, status);

          numOutVariable += 3*geomInVal->length; // xD1, yD1, zD1, ...
        }
    }

    if (numOutVariable > 99999999) {
        AIM_ERROR(aimInfo, "Array of design variable names will be over-run!");
        status = CAPS_RANGEERR;
        goto cleanup;
    }

    // Allocate our names
    AIM_ALLOC(dataOutName, numOutVariable, char*, aimInfo, status);
    for (i = 0; i < numOutVariable; i++) dataOutName[i] = NULL;

    stringLength = 11;
    j = 1;
    k = 1;
    for (i = 0; i < numOutVariable; i++) {
        AIM_ALLOC(dataOutName[i], stringLength+1, char, aimInfo, status);

        // Set names
        if      (i == 0) sprintf(dataOutName[i], "%s", "x");
        else if (i == 1) sprintf(dataOutName[i], "%s", "y");
        else if (i == 2) sprintf(dataOutName[i], "%s", "z");
        else if (i == 3) sprintf(dataOutName[i], "%s", "id");
        else {

            if      (j == 1) sprintf(dataOutName[i], "%s%d", "xD", k);
            else if (j == 2) sprintf(dataOutName[i], "%s%d", "yD", k);
            else if (j == 3) {
                sprintf(dataOutName[i], "%s%d", "zD", k);
                j = 0;
                k += 1;
            }
            j += 1;
        }
    }

    AIM_ALLOC(dataOutFormat, numOutVariable, int, aimInfo, status);

    // Set data out formatting
    for (i = 0; i < numOutVariable; i++) {
        if (strcasecmp(dataOutName[i], "id") == 0) {
            dataOutFormat[i] = (int) Integer;
        } else {
            dataOutFormat[i] = (int) Double;
        }
    }

    // Allocate data arrays that are going to be output
    AIM_ALLOC(dataOutMatrix, meshRef->nmap, double**, aimInfo, status);
    for (i = 0; i < meshRef->nmap; i++) dataOutMatrix[i] = NULL;

    for (i = 0; i < meshRef->nmap; i++) {
      AIM_ALLOC(dataOutMatrix[i], numOutVariable, double*, aimInfo, status);
      for (j = 0; j < numOutVariable; j++) dataOutMatrix[i][j] = NULL;
    }

    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &numOutDataPoint);
      AIM_STATUS(aimInfo, status);

      // allocate data arrays
      for (j = 0; j < numOutVariable; j++) {
        AIM_ALLOC(dataOutMatrix[i][j], numOutDataPoint, double, aimInfo, status);
      }

      for (iglobal = 0; iglobal < numOutDataPoint; iglobal++) {
        status = EG_getGlobal(meshRef->maps[i].tess,
                              iglobal+1, &ptype, &pindex, xyz);
        AIM_STATUS(aimInfo, status);

        // First just set the Coordinates
        dataOutMatrix[i][0][iglobal] = xyz[0];
        dataOutMatrix[i][1][iglobal] = xyz[1];
        dataOutMatrix[i][2][iglobal] = xyz[2];

        // Volume mesh node ID
        dataOutMatrix[i][3][iglobal] = meshRef->maps[i].map[iglobal];
      }
    }

    // Loop over the geometry in values and compute sensitivities for all bodies
    m = 4;
    for (j = 0; j < numDesignVariable; j++) {

      if (geomSelect[j] == (int) false) continue;

      geomInName = designVariable[j].name;
      index = aim_getIndex(aimInfo, geomInName, GEOMETRYIN);

      status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
      AIM_STATUS(aimInfo, status);

      for (row = 0; row < geomInVal->nrow; row++) {
        for (col = 0; col < geomInVal->ncol; col++) {

          for (i = 0; i < meshRef->nmap; i++) {
            status = aim_tessSensitivity(aimInfo,
                                         geomInName,
                                         row+1, col+1, // row, col
                                         meshRef->maps[i].tess,
                                         &numPoint, &dxyz);
            AIM_STATUS(aimInfo, status, "Sensitivity for: %s\n", geomInName);
            AIM_NOTNULL(dxyz, aimInfo, status);

            for (k = 0; k < numPoint; k++) {
              dataOutMatrix[i][m+0][k] = dxyz[3*k + 0]; // dx/dGeomIn
              dataOutMatrix[i][m+1][k] = dxyz[3*k + 1]; // dy/dGeomIn
              dataOutMatrix[i][m+2][k] = dxyz[3*k + 2]; // dz/dGeomIn
            }
            AIM_FREE(dxyz);
          }
          m += 3;
        }
      }
    }

    // Write sensitivity files for each body tessellation

    for (i = 0; i < meshRef->nmap; i++) {
      status = EG_statusTessBody(meshRef->maps[i].tess, &body, &state, &numOutDataPoint);
      AIM_STATUS(aimInfo, status);

      // check if the tessellation has a mixture of quad and tess
      status = EG_attributeRet(meshRef->maps[i].tess, ".mixed",
                               &atype, &alen, &nquad, &reals, &string);
      if (status != EGADS_SUCCESS &&
          status != EGADS_NOTFOUND) AIM_STATUS(aimInfo, status);

      status = EG_getBodyTopos(body, NULL, FACE, &nFace, &faces);
      AIM_STATUS(aimInfo, status);

      ielem = 0;
      numOutDataConnect = 0;
      for (iface = 0; iface < nFace; iface++) {
        // get the face tessellation
        status = EG_getTessFace(meshRef->maps[i].tess, iface+1, &alen, &face_xyz, &face_uv,
                                &face_ptype, &face_pindex, &ntri, &face_tris, &face_tric);
        AIM_STATUS(aimInfo, status);

        if (nquad == NULL) { // all triangles

          // re-allocate data arrays
          AIM_REALL(dataConnectMatrix, 4*(numOutDataConnect + ntri), int, aimInfo, status);

          for (itri = 0; itri < ntri; itri++, ielem++) {
            for (j = 0; j < 3; j++) {
              status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+j], &iglobal);
              AIM_STATUS(aimInfo, status);
              dataConnectMatrix[4*ielem+j] = iglobal;
            }
            // repeat the last node for triangles
            dataConnectMatrix[4*ielem+3] = dataConnectMatrix[4*ielem+2];
          }

          numOutDataConnect += ntri;

        } else { // mixture of tri and quad elements

          // re-allocate data arrays
          AIM_REALL(dataConnectMatrix, 4*(numOutDataConnect + ntri-nquad[iface]), int, aimInfo, status);

          // process triangles
          for (itri = 0; itri < ntri-2*nquad[iface]; itri++, ielem++) {
            for (j = 0; j < 3; j++) {
              status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+j], &iglobal);
              AIM_STATUS(aimInfo, status);
              dataConnectMatrix[4*ielem+j] = iglobal;
            }
            // repeat the last node for triangle
            dataConnectMatrix[4*ielem+3] = dataConnectMatrix[4*ielem+2];
          }
          // process quads
          for (; itri < ntri; itri++, ielem++) {
            for (j = 0; j < 3; j++) {
              status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+j], &iglobal);
              AIM_STATUS(aimInfo, status);
              dataConnectMatrix[4*ielem+j] = iglobal;
            }

            // add the last node from the 2nd triangle to make the quad
            itri++;
            status = EG_localToGlobal(meshRef->maps[i].tess, iface+1, face_tris[3*itri+2], &iglobal);
            AIM_STATUS(aimInfo, status);
            dataConnectMatrix[4*ielem+3] = iglobal;
          }

          numOutDataConnect += ntri-nquad[iface];
        }

      }
      AIM_FREE(faces);

      if (writeSensitivity == (int)true)
        sprintf(message,"%s %d,", "sensitivity file for body", i+1);
      else
        sprintf(message,"%s %d,", "mesh file for body", i+1);

      stringLength = strlen(folder) + 1 + strlen(filePre) + 7 + strlen(fileExt) + 1 ;

      AIM_REALL(filename, stringLength, char, aimInfo, status);

#ifdef WIN32
      sprintf(filename, "%s\\%s%d%s", folder, filePre, i+1, fileExt);
#else
      sprintf(filename, "%s/%s%d%s",  folder, filePre, i+1, fileExt);
#endif

      sprintf(zoneTitle, "%s_%d", "Body", i+1);
      /*@-nullpass@*/
      status = tecplot_writeFEPOINT(aimInfo,
                                    filename,
                                    message,
                                    zoneTitle,
                                    numOutVariable,
                                    dataOutName,
                                    numOutDataPoint,
                                    dataOutMatrix[i],
                                    dataOutFormat,
                                    numOutDataConnect,
                                    dataConnectMatrix,
                                    NULL);
      /*@+nullpass@*/
      AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;

cleanup:
    (void) string_freeArray(numOutVariable, &dataOutName);
#ifdef S_SPLINT_S
    EG_free(dataOutName);
#endif

    if (dataOutMatrix != NULL) {
        for (i = 0; i < meshRef->nmap; i++) {
          if (dataOutMatrix[i] != NULL) {
              for (j = 0; j < numOutVariable; j++) {
                  AIM_FREE(dataOutMatrix[i][j]);
              }
          }
          AIM_FREE(dataOutMatrix[i]);
        }
    }

    AIM_FREE(dataOutMatrix);
    AIM_FREE(dataOutFormat);
    AIM_FREE(dataConnectMatrix);
    AIM_FREE(faces);

    AIM_FREE(dxyz);
    AIM_FREE(filename);
    AIM_FREE(geomSelect);

    return status;
}


static int _writeFunctinoalComponent(void *aimInfo, FILE *fp, cfdDesignFunctionalCompStruct *comp)
{
    int status = CAPS_SUCCESS;

    const char *names[] = {"cl", "cd",                   // Lift, drag coefficients
                           "clp", "cdp",                 // Lift, drag coefficients: pressure contributions
                           "clv", "cdv",                 // Lift, drag coefficients: shear contributions
                           "cmx", "cmy", "cmz",          // x/y/z-axis moment coefficients
                           "cmxp", "cmyp", "cmzp",       // x/y/z-axis moment coefficients: pressure contributions
                           "cmxv", "cmyv", "cmzv",       // x/y/z-axis moment coefficients: shear contributions
                           "cx", "cy", "cz",             // x/y/z-axis force coefficients
                           "cxp", "cyp", "czp",          // x/y/z-axis force coefficients: pressure contributions
                           "cxv", "cyv", "czv",          // x/y/z-axis force coefficients: shear contributions
                           "powerx", "powery", "powerz", // x/y/z-axis power coefficients
                           "clcd",                       // Lift-to-drag ratio
                           "fom"    ,                    // Rotorcraft figure of merit
                           "propeff",                    // Rotorcraft propulsive efficiency
                           "rtr"    ,                    // thrust Rotorcraft thrust function
                           "pstag"  ,                    // RMS of stagnation pressure in cutting plane disk
                           "distort",                    // Engine inflow distortion
                           "boom"   ,                    // targ Near-field p/p∞ pressure target
                           "sboom"  ,                    // Coupled sBOOM ground-based noise metrics
                           "ae"     ,                    // Supersonic equivalent area target distribution
                           "press"  ,                    // box RMS of pressure in user-defined box, also pointwise dp/dt, dρ/dt
                           "cpstar" ,                    // Target pressure distributions
                           "sgen"                        // Entropy generation
                          };

    int i, found = (int)false;
    char *function=NULL, *c=NULL;

    for (i = 0; i < sizeof(names)/sizeof(char*); i++)
        if ( strcasecmp(names[i], comp->name) == 0 ) {
            found = (int)true;
            break;
        }

    if (found == (int) false) {
        AIM_ERROR(aimInfo, "Unknown function: '%s'", comp->name);
        AIM_ADDLINE(aimInfo, "Available functions:");
        for (i = 0; i < sizeof(names)/sizeof(char*); i++)
            AIM_ADDLINE(aimInfo, "'%s'", names[i]);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // make the name lower case
    AIM_STRDUP(function, comp->name, aimInfo, status);
    for (c = function; *c != '\0'; ++c) *c = tolower(*c);

    // fprintf(fp, "Components of function   1: boundary id (0=all)/name/value/weight/target/power\n");
    fprintf(fp, " %d %s          %f %f %f %f\n", comp->bcID,
                                                 function,
                                                 0.0,
                                                 comp->weight,
                                                 comp->target,
                                                 comp->power);

    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(function);
    return status;
}


// Write FUN3D  rubber.data file
int fun3d_writeRubber(void *aimInfo,
                      cfdDesignStruct design,
                      double fun3dVersion,
                      aimMeshRef *meshRef)
{
    int status; // Function return status

    int i, j, k, m; // Indexing

    int numBody = 0, numShapeVar = 0;

    char file[] = "rubber.data";
    FILE *fp = NULL;

    printf("Writing %s \n", file);
    if (meshRef == NULL) return CAPS_NULLVALUE;

    // Open file
    fp = aim_fopen(aimInfo, file, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file: %s\n", file);
        status = CAPS_IOERR;
        goto cleanup;
    }

    numBody = meshRef->nmap;

    // Determine number of geometry input variables
    for (i = 0; i < design.numDesignVariable; i++) {

        if (design.designVariable[i].type != DesignVariableGeometry) continue;

        numShapeVar += design.designVariable[i].var->length; // xD1, yD1, zD1, ...
    }

    fprintf(fp, "################################################################################\n");
    fprintf(fp, "########################### Design Variable Information ########################\n");
    fprintf(fp, "################################################################################\n");
    fprintf(fp, "Global design variables (Mach number, AOA, Yaw, Noninertial rates)\n");
    fprintf(fp, "Var Active         Value               Lower Bound            Upper Bound\n");

    for (i = 0; i < design.numDesignVariable; i++) {

        if (strcasecmp(design.designVariable[i].name, "Mach") != 0) continue;

        if (design.designVariable[i].var->length < 1) {
            status = CAPS_RANGEERR;
            goto cleanup;
        }

        fprintf(fp, "Mach    1   %.15E  %.15E  %.15E\n", design.designVariable[i].value[0],
                                                         design.designVariable[i].lowerBound[0],
                                                         design.designVariable[i].upperBound[0]);
        break;
    }

    if (i >= design.numDesignVariable) fprintf(fp, "Mach    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");

    for (i = 0; i < design.numDesignVariable; i++) {

        if (strcasecmp(design.designVariable[i].name, "Alpha") != 0) continue;

        if (design.designVariable[i].var->length < 1) {
            status = CAPS_RANGEERR;
            goto cleanup;
        }

        fprintf(fp, "AOA     1   %.15E  %.15E  %.15E\n", design.designVariable[i].value[0],
                                                         design.designVariable[i].lowerBound[0],
                                                         design.designVariable[i].upperBound[0]);
        break;
    }

    if (i >= design.numDesignVariable) fprintf(fp, "AOA     0   0.000000000000000E+00  0.000000000000000E+00  0.00000000000000E+00\n");


    if (fun3dVersion > 12.4) {
        // FUN3D version 13.1 - version 12.4 doesn't have these available

        for (i = 0; i < design.numDesignVariable; i++) {

            if (strcasecmp(design.designVariable[i].name, "Beta") != 0) continue;

            if (design.designVariable[i].var->length < 1) {
                status = CAPS_RANGEERR;
                goto cleanup;
            }

            fprintf(fp, "Yaw     1   %.15E  %.15E  %.15E\n", design.designVariable[i].value[0],
                                                             design.designVariable[i].lowerBound[0],
                                                             design.designVariable[i].upperBound[0]);
            break;
        }

        if (i >= design.numDesignVariable) fprintf(fp, "Yaw     0   0.000000000000000E+00  0.000000000000000E+00  0.00000000000000E+00\n");

        for (i = 0; i < design.numDesignVariable; i++) {

            if (strcasecmp(design.designVariable[i].name, "NonInertial_Rotation_Rate") != 0) continue;

            if (design.designVariable[i].var->length < 3) {
                status = CAPS_RANGEERR;
                goto cleanup;
            }

            fprintf(fp, "xrate   1   %.15E  %.15E  %.15E\n", design.designVariable[i].value[0],
                                                             design.designVariable[i].lowerBound[0],
                                                             design.designVariable[i].upperBound[0]);
            break;
        }

        if (i >= design.numDesignVariable) fprintf(fp, "xrate   0   0.000000000000000E+00  0.000000000000000E+00  0.00000000000000E+00\n");

        for (i = 0; i < design.numDesignVariable; i++) {

            if (strcasecmp(design.designVariable[i].name, "NonInertial_Rotation_Rate") != 0) continue;

            if (design.designVariable[i].var->length < 3) {
                status = CAPS_RANGEERR;
                goto cleanup;
            }

            fprintf(fp, "yrate   1   %.15E  %.15E  %.15E\n", design.designVariable[i].value[1],
                                                             design.designVariable[i].lowerBound[1],
                                                             design.designVariable[i].upperBound[1]);
            break;
        }

        if (i >= design.numDesignVariable) fprintf(fp, "yrate   0   0.000000000000000E+00  0.000000000000000E+00  0.00000000000000E+00\n");

        for (i = 0; i < design.numDesignVariable; i++) {

            if (strcasecmp(design.designVariable[i].name, "NonInertial_Rotation_Rate") != 0) continue;

            if (design.designVariable[i].var->length < 3) {
                status = CAPS_RANGEERR;
                goto cleanup;
            }

            fprintf(fp, "zrate   1   %.15E  %.15E  %.15E\n", design.designVariable[i].value[2],
                                                             design.designVariable[i].lowerBound[2],
                                                             design.designVariable[i].upperBound[2]);
            break;
        }

        if (i >= design.numDesignVariable) fprintf(fp, "zrate   0   0.000000000000000E+00  0.000000000000000E+00  0.00000000000000E+00\n");
    }

    fprintf(fp, "Number of bodies\n");
    fprintf(fp, "%d\n", numBody);
    for (i = 0; i < numBody; i++) {
        fprintf(fp, "Rigid motion design variables for 'Body %d'\n", i+1);
        fprintf(fp, "Var Active         Value               Lower Bound            Upper Bound\n");
        fprintf(fp, "RotRate  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotFreq  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotAmpl  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotOrgx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotOrgy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotOrgz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotVecx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotVecy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "RotVecz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "TrnRate  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "TrnFreq  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "TrnAmpl  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "TrnVecx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "TrnVecy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        fprintf(fp, "TrnVecz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00\n");
        if (fun3dVersion > 12.4) {
            // FUN3D version 13.1 - version 12.4 doesn't have option 5
            fprintf(fp, "Parameterization Scheme (Massoud=1 Bandaids=2 Sculptor=4 User-Defined=5)\n");
            fprintf(fp, "5\n");
        } else {

            fprintf(fp, "Parameterization Scheme (Massoud=1 Bandaids=2 Sculptor=4)\n");
            fprintf(fp, "1\n");
        }

        fprintf(fp, "Number of shape variables for 'Body %d'\n", i+1);
        fprintf(fp, "%d\n", numShapeVar);
        fprintf(fp, "Index Active         Value               Lower Bound            Upper Bound\n");

        m = 1;
        for (j = 0; j < design.numDesignVariable; j++) {

            if (design.designVariable[j].type != DesignVariableGeometry) continue;

            for (k = 0; k < design.designVariable[j].var->length; k++ ) {

                fprintf(fp, "%d    1   %.15E  %.15E  %.15E\n", m, design.designVariable[j].value[k],
                                                                  design.designVariable[j].lowerBound[k],
                                                                  design.designVariable[j].upperBound[k]);

                m += 1;
            }
        }
    }

    fprintf(fp, "################################################################################\n");
    fprintf(fp, "############################### Function Information ###########################\n");
    fprintf(fp, "################################################################################\n");
    fprintf(fp, "Number of composite functions for design problem statement\n");
    fprintf(fp, "%d\n", design.numDesignFunctional);

    for (i = 0; i < design.numDesignFunctional; i++) {
        fprintf(fp, "################################################################################\n");
        fprintf(fp, "Cost function (1) or constraint (2)\n");
        fprintf(fp, "1\n");
        fprintf(fp, "If constraint, lower and upper bounds\n");
        fprintf(fp, "0.0 0.0\n");
        fprintf(fp, "Number of components for function   %d\n", i+1);
        fprintf(fp, "%d\n", design.designFunctional[i].numComponent);
        fprintf(fp, "Physical timestep interval where function is defined\n");
        fprintf(fp, "1 1\n");
        fprintf(fp, "Composite function weight, target, and power\n");
        fprintf(fp, "1.0 0.0 1.0\n");
        fprintf(fp, "Components of function   %d: boundary id (0=all)/name/value/weight/target/power\n", i+1);

        for (j = 0; j < design.designFunctional[i].numComponent; j++) {
            //fprintf(fp, "0 clcd          0.000000000000000    1.000   10.00000 2.000\n");
            status = _writeFunctinoalComponent(aimInfo, fp, &design.designFunctional[i].component[j]);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        fprintf(fp, "Current value of function   %d\n", i+1);
        fprintf(fp, "0.000000000000000\n");
        fprintf(fp, "Current derivatives of function wrt global design variables\n");
        fprintf(fp, "0.000000000000000\n");
        fprintf(fp, "0.000000000000000\n");

        if (fun3dVersion > 12.4) {
            // FUN3D version 13.1 - version 12.4 doesn't have these available
            fprintf(fp, "0.000000000000000\n"); // Yaw
            fprintf(fp, "0.000000000000000\n"); // xrate
            fprintf(fp, "0.000000000000000\n"); // yrate
            fprintf(fp, "0.000000000000000\n"); // zrate
        }

        for (j = 0; j < numBody; j++) {
            fprintf(fp, "Current derivatives of function wrt rigid motion design variables of body %d\n", j+1);
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "0.000000000000000\n");
            fprintf(fp, "Current derivatives of function wrt shape design variables of body %d\n", j+1);

            for (k = 0; k < design.numDesignVariable; k++) {

                if (design.designVariable[k].type != DesignVariableGeometry) continue;

                for (m = 0; m < design.designVariable[k].var->length; m++ ) fprintf(fp, "0.000000000000000\n");
            }
        }
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("Error: Premature exit in fun3d_writeRubber status = %d\n",
               status);

    if (fp != NULL) fclose(fp);

    return status;

    /* Version 13.1 - template
	################################################################################
	########################### Design Variable Information ########################
	################################################################################
	Global design variables (Mach number, AOA, Yaw, Noninertial rates)
	Var Active         Value               Lower Bound            Upper Bound
	 Mach    0   0.800000000000000E+00  0.000000000000000E+00  0.900000000000000E+00
	  AOA    1   1.000000000000000E+00  0.000000000000000E+00  5.000000000000000E+00
	  Yaw    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	xrate    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	yrate    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	zrate    0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	Number of bodies
	2
	Rigid motion design variables for 'wing'
	Var Active         Value               Lower Bound            Upper Bound
	RotRate  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotFreq  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotAmpl  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotOrgx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotOrgy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotOrgz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotVecx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotVecy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotVecz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnRate  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnFreq  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnAmpl  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnVecx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnVecy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnVecz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	Parameterization Scheme (Massoud=1 Bandaids=2 Sculptor=4 User-Defined=5)
	1
	Number of shape variables for 'wing'
	3
	Index Active         Value               Lower Bound            Upper Bound
	1    1   1.000000000000000E+00  0.000000000000000E+00  2.000000000000000E+00
	2    1   1.000000000000000E+00  0.000000000000000E+00  2.000000000000000E+00
	3    1   1.000000000000000E+00  0.000000000000000E+00  2.000000000000000E+00
	Rigid motion design variables for 'tail'
	Var Active         Value               Lower Bound            Upper Bound
	RotRate  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotFreq  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotAmpl  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotOrgx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotOrgy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotOrgz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotVecx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotVecy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	RotVecz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnRate  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnFreq  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnAmpl  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnVecx  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnVecy  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	TrnVecz  0   0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
	Parameterization Scheme (Massoud=1 Bandaids=2 Sculptor=4 User-Defined=5)
	2
	Number of shape variables for 'tail'
	2
	Index Active         Value               Lower Bound            Upper Bound
	1    1   2.000000000000000E+00 -1.000000000000000E+00  5.000000000000000E+00
	2    1   2.000000000000000E+00 -1.000000000000000E+00  5.000000000000000E+00
	################################################################################
	############################### Function Information ###########################
	################################################################################
	Number of composite functions for design problem statement
	2
	################################################################################
	Cost function (1) or constraint (2)
	1
	If constraint, lower and upper bounds
	0.0 0.0
	Number of components for function   1
	1
	Physical timestep interval where function is defined
	1 1
	Composite function weight, target, and power
	1.0 0.0 1.0
	Components of function   1: boundary id (0=all)/name/value/weight/target/power
	0 cl            0.000000000000000    1.000   10.00000 2.000
	Current value of function   1
	0.000000000000000
	Current derivatives of function wrt global design variables
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt rigid motion design variables of body   1
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt shape design variables of body   1
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt rigid motion design variables of body   2
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt shape design variables of body   2
	0.000000000000000
	0.000000000000000
	################################################################################
	Cost function (1) or constraint (2)
	2
	If constraint, lower and upper bounds
	-0.03 -0.01
	Number of components for function   2
	1
	Physical timestep interval where function is defined
	1 1
	Composite function weight, target, and power
	1.0 0.0 1.0
	Components of function   2: boundary id (0=all)/name/value/weight/target/power
	0 cmy           0.000000000000000    1.000    0.00000 1.000
	Current value of function   2
	0.000000000000000
	Current derivatives of function wrt global design variables
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt rigid motion design variables of body   1
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt shape design variables of body   1
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt rigid motion design variables of body   2
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	0.000000000000000
	Current derivatives of function wrt shape design variables of body   2
	0.000000000000000
	0.000000000000000
     */
}

// Finds a line in rubber.data that matches header
static int findHeader(const char *header, char **line, size_t *nline, int *iline, FILE *fp)
{
    int i;
    while (getline(line, nline, fp) != -1) {
        (*iline)++;
        if (nline == 0) continue;
        if ((*line)[0] == '#') continue;
        i = 0;
        while ((*line)[i] == ' ' && (*line)[i] != '\0') i++;
        if (strncasecmp(header, &(*line)[i], strlen(header)) == 0)
          break;
    }
    if (feof(fp)) return CAPS_IOERR;

    return CAPS_SUCCESS;
}


// Read objective value and derivatives from FUN3D rubber.data file
int fun3d_readRubber(void *aimInfo,
                     cfdDesignStruct design,
                     double fun3dVersion)
{
    int status; // Function return status

    int i, j, k, ibody; // Indexing

    size_t nline;
    int iline=0;
    char *line=NULL;
    char file[] = "rubber.data";
    double dfun_dvar = 0;
    FILE *fp = NULL;

    int numBody;
    const char *intents;
    ego *bodies;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    printf("Reading %s \n", file);

    // Open file
    fp = aim_fopen(aimInfo, file, "r");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file: %s\n", file);
        status = CAPS_IOERR;
        goto cleanup;
    }

    for (i = 0; i < design.numDesignFunctional; i++) {
        status = findHeader("Current value of function", &line, &nline, &iline, fp);
        AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);

        // Get the line with the objective value
        status = getline(&line, &nline, fp); iline++;
        if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
        AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
        AIM_NOTNULL(line, aimInfo, status);

        status = sscanf(line, "%lf", &design.designFunctional[i].value);
        if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
        AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);

        // skip: Current derivatives of function wrt rigid motion design variables of
        status = getline(&line, &nline, fp); iline++;
        if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
        AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
        AIM_NOTNULL(line, aimInfo, status);

        // read the dFun/dMach
        status = getline(&line, &nline, fp); iline++;
        if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
        AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
        AIM_NOTNULL(line, aimInfo, status);

        for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {

            if (strcasecmp(design.designFunctional[i].dvar[j].name, "Mach") != 0) continue;

            status = sscanf(line, "%lf", &design.designFunctional[i].dvar[j].value[0]);
            if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
            AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
            break;
        }

        // read the dFun/dAOA
        status = getline(&line, &nline, fp); iline++;
        if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
        AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
        AIM_NOTNULL(line, aimInfo, status);

        for (j = 0; j < design.numDesignVariable; j++) {

            if (strcasecmp(design.designFunctional[i].dvar[j].name, "Alpha") != 0) continue;

            status = sscanf(line, "%lf", &design.designFunctional[i].dvar[j].value[0]);
            if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
            AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
            break;
        }

        if (fun3dVersion > 12.4) {
            // FUN3D version 13.1 - version 12.4 doesn't have these available

            // read the dFun/dBeta
            status = getline(&line, &nline, fp); iline++;
            if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
            AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
            AIM_NOTNULL(line, aimInfo, status);

            for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {

                if (strcasecmp(design.designFunctional[i].dvar[j].name, "Beta") != 0) continue;

                status = sscanf(line, "%lf", &design.designFunctional[i].dvar[j].value[0]);
                if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
                AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
                break;
            }

            // read the dFun/dNonInertial_Rotation_Rate[0]
            status = getline(&line, &nline, fp); iline++;
            if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
            AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
            AIM_NOTNULL(line, aimInfo, status);

            for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {

                if (strcasecmp(design.designFunctional[i].dvar[j].name, "NonInertial_Rotation_Rate") != 0) continue;

                status = sscanf(line, "%lf", &design.designFunctional[i].dvar[j].value[0]);
                if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
                AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
                break;
            }

            // read the dFun/dNonInertial_Rotation_Rate[1]
            status = getline(&line, &nline, fp); iline++;
            if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
            AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
            AIM_NOTNULL(line, aimInfo, status);

            for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {

                if (strcasecmp(design.designFunctional[i].dvar[j].name, "NonInertial_Rotation_Rate") != 0) continue;

                status = sscanf(line, "%lf", &design.designFunctional[i].dvar[j].value[1]);
                if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
                AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
                break;
            }

            // read the dFun/dNonInertial_Rotation_Rate[2]
            status = getline(&line, &nline, fp); iline++;
            if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
            AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
            AIM_NOTNULL(line, aimInfo, status);

            for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {

                if (strcasecmp(design.designFunctional[i].dvar[j].name, "NonInertial_Rotation_Rate") != 0) continue;

                status = sscanf(line, "%lf", &design.designFunctional[i].dvar[j].value[2]);
                if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
                AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
                break;
            }
        }


        // zero out previous derivatives
        for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {
            if (design.designFunctional[i].dvar[j].type != DesignVariableGeometry) continue;
            for (k = 0; k < design.designVariable[j].var->length; k++ ) {
              design.designFunctional[i].dvar[j].value[k] = 0;
            }
        }

        for (ibody = 0; ibody < numBody; ibody++) {

          // Rigid motion design variables
          status = findHeader("Current derivatives of", &line, &nline, &iline, fp);
          AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);

          // Skip reading rigid motion design variables for now
          
          // Read shape design variables
          status = findHeader("Current derivatives of", &line, &nline, &iline, fp);
          AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);

          for (j = 0; j < design.designFunctional[i].numDesignVariable; j++) {

              if (design.designFunctional[i].dvar[j].type != DesignVariableGeometry) continue;

              for (k = 0; k < design.designVariable[j].var->length; k++ ) {

                  // read the dFun/dDesignVar[k]
                  status = getline(&line, &nline, fp); iline++;
                  if (status == -1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
                  AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);
                  AIM_NOTNULL(line, aimInfo, status);

                  status = sscanf(line, "%lf", &dfun_dvar);
                  if (status != 1) status = CAPS_IOERR; else status = CAPS_SUCCESS;
                  AIM_STATUS(aimInfo, status, "rubber.data line %d", iline);

                  // accumulate derivatives across bodies
                  design.designFunctional[i].dvar[j].value[k] += dfun_dvar;
              }
          }
        }
    }

    status = CAPS_SUCCESS;

cleanup:

    if (line != NULL) free(line); // Must use free!
    if (fp   != NULL) fclose(fp);

    return status;
}

// Make FUN3D directory structure/tree
int fun3d_makeDirectory(void *aimInfo)
{
    int status; // Function return status

    char filename[PATH_MAX];

    char flow[] = "Flow";
    char datafile[] = "datafiles";
    char adjoint[] = "Adjoint";
    char rubber[] = "Rubberize";

    printf("Creating FUN3D directory tree\n");

    // Flow
    status = aim_mkDir(aimInfo, flow);
    AIM_STATUS(aimInfo, status);

    // Datafiles
    sprintf(filename, "%s/%s", flow, datafile);
    status = aim_mkDir(aimInfo, filename);
    AIM_STATUS(aimInfo, status);

    // Adjoint
    status = aim_mkDir(aimInfo, adjoint);
    AIM_STATUS(aimInfo, status);

    // Rubber
    status = aim_mkDir(aimInfo, rubber);
    AIM_STATUS(aimInfo, status);

    status = CAPS_SUCCESS;

cleanup:
    return status;
}
