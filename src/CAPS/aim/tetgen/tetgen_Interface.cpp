// Tetgen interface functions - Written by Dr. Ryan Durscher AFRL/RQVC
// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <vector>
#include <set>
#include <iostream>

#include "aimUtil.h"

#include "tetgen.h"
#include "meshTypes.h"
#include "egads.h"
#include "capsTypes.h"

#include "meshUtils.h"

#include "tetgen_Interface.hpp"

#define CROSS(a,b,c)  a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                      a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                      a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)     (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


static int
tetgen_to_MeshStruct(tetgenio *mesh, meshStruct *genUnstrMesh)
{
    int status; // Function return status

    int i, j, elementIndex; // Indexing

    int numPoint;
    int defaultVolID = 1; // Default volume ID

    meshAnalysisTypeEnum analysisType;

    analysisType = genUnstrMesh->analysisType;

    // Cleanup existing node and elements
    (void) destroy_meshNodes(genUnstrMesh);

    (void) destroy_meshElements(genUnstrMesh);

    (void) destroy_meshQuickRefStruct(&genUnstrMesh->meshQuickRef);

    genUnstrMesh->meshType = VolumeMesh;

    // Numbers
    genUnstrMesh->numNode = mesh->numberofpoints;
    genUnstrMesh->numElement = mesh->numberoftrifaces + mesh->numberoftetrahedra;

    genUnstrMesh->meshQuickRef.useStartIndex = (int) true;
    genUnstrMesh->meshQuickRef.numTriangle = mesh->numberoftrifaces;
    genUnstrMesh->meshQuickRef.numTetrahedral = mesh->numberoftetrahedra;
    genUnstrMesh->meshQuickRef.startIndexTriangle = 0;
    genUnstrMesh->meshQuickRef.startIndexTetrahedral = mesh->numberoftrifaces;

    // Nodes - allocate
    genUnstrMesh->node = (meshNodeStruct *) EG_alloc(genUnstrMesh->numNode*sizeof(meshNodeStruct));
    if (genUnstrMesh->node == NULL) return EGADS_MALLOC;

    // Nodes - set
    for (i = 0; i < genUnstrMesh->numNode; i++) {

        // Initiate node
        status = initiate_meshNodeStruct(&genUnstrMesh->node[i], analysisType);
        if (status != CAPS_SUCCESS) return status;

        // Copy node data
        genUnstrMesh->node[i].nodeID = i+1;

        genUnstrMesh->node[i].xyz[0] = mesh->pointlist[3*i+0];
        genUnstrMesh->node[i].xyz[1] = mesh->pointlist[3*i+1];
        genUnstrMesh->node[i].xyz[2] = mesh->pointlist[3*i+2];

        /*
        if (mesh->numberofpointattributes != 0) {
        mesh->pointmarkerlist != NULL
        }
         */
    }

    // Elements - allocate
    genUnstrMesh->element = (meshElementStruct *) EG_alloc(genUnstrMesh->numElement*sizeof(meshElementStruct));
    if (genUnstrMesh->element == NULL) return EGADS_MALLOC;

    elementIndex = 0;
    numPoint = 0;
    // Elements -Set triangles
    for (i = 0; i < mesh->numberoftrifaces; i++) {
        status = initiate_meshElementStruct(&genUnstrMesh->element[elementIndex], analysisType);
        if (status != CAPS_SUCCESS) return status;

        genUnstrMesh->element[elementIndex].elementType = Triangle;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID = mesh->trifacemarkerlist[i];

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) return status;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] = mesh->trifacelist[numPoint*i+j];
        }

        elementIndex += 1;
    }

    //Set tetrahedral
    numPoint = 0;
    for (i = 0; i < mesh->numberoftetrahedra; i++) {
        status = initiate_meshElementStruct(&genUnstrMesh->element[elementIndex], analysisType);
        if (status != CAPS_SUCCESS) return status;

        genUnstrMesh->element[elementIndex].elementType = Tetrahedral;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;

        if (mesh->numberoftetrahedronattributes != 0) {
            genUnstrMesh->element[elementIndex].markerID = (int) mesh->tetrahedronattributelist[mesh->numberoftetrahedronattributes*i + 0];
        } else {
            genUnstrMesh->element[elementIndex].markerID = defaultVolID;
        }

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) return status;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] = mesh->tetrahedronlist[numPoint*i+j];
        }

        elementIndex += 1;
    }

    return CAPS_SUCCESS;
}

static int
writeBinaryUgrid(void *aimInfo, const char *fileName, tetgenio *mesh)
{
  int    status = CAPS_SUCCESS;

  int    nTri, nQuad;
  int    nTet, nPyramid, nPrism, nHex;
  char   aimFile[PATH_MAX];
  FILE *fp = NULL;

  snprintf(aimFile, PATH_MAX, "%s.lb8.ugrid", fileName);

  fp = fopen(aimFile, "wb");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open file: tetgen.lb8.ugrid");
    status = CAPS_IOERR;
    goto cleanup;
  }

  nTri     = mesh->numberoftrifaces;
  nQuad    = 0;
  nTet     = mesh->numberoftetrahedra;
  nPyramid = 0;
  nPrism   = 0;
  nHex     = 0;

  /* read a binary UGRID file */
  status = fwrite(&mesh->numberofpoints, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nTri,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nQuad,    sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nTet,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nPyramid, sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nPrism,   sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }
  status = fwrite(&nHex,     sizeof(int), 1, fp);
  if (status != 1) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  /* write all of the vertices */
  status = fwrite(mesh->pointlist, sizeof(double), 3*mesh->numberofpoints, fp);
  if (status != 3*mesh->numberofpoints) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  /* write the triangle connectivity */
  status = fwrite(mesh->trifacelist, sizeof(int), 3*nTri, fp);
  if (status != 3*nTri) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  /* write the BC ID */
  status = fwrite(mesh->trifacemarkerlist, sizeof(int), nTri, fp);
  if (status != nTri) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  /* write the tetrahedral connectivity */
  status = fwrite(mesh->tetrahedronlist, sizeof(int), 4*nTet, fp);
  if (status != 4*nTet) { status = CAPS_IOERR; AIM_STATUS(aimInfo, status); }

  status = CAPS_SUCCESS;

cleanup:
  if (fp != NULL) fclose(fp);

  return status;
}


extern "C"
int tetgen_VolumeMesh(void *aimInfo,
                      meshInputStruct meshInput,
                      const char *fileName,
                      meshStruct *surfaceMesh,
                      meshStruct *volumeMesh)
{
    // tetgen documentation:
    //
    // -p Generate tetrahedra
    // -Y Preserves the input surface mesh (does not modify it).
    // -V verbose
    // -q mesh quality (maximum radius-edge ratio)/(minimum dihedral angle)
    // -a maximum volume constraint
    // -f provides the interior+boundry triangular faces
    // -nn to get tet neighbors for each triangular face
    // -k dumps to paraview when last argument is NULL
    // -m Applies a mesh sizing function

    // Initialize variables
    int status = 0; //Function status return
    int i, j; // Indexing
    char *inputString = NULL; // Input string to tetgen

    // Create Tetgen input string
    char temp[120] = "p"; // Tetrahedralize a piecewise linear complex flag
    char q[80];


    // TetGen variables
    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *p;

    printf("\nGenerating volume mesh using TetGen.....\n");

    // First index - either 0 or 1
    in.firstnumber = 1;

    // Check inputs
    if (surfaceMesh->numNode    == 0) return CAPS_BADVALUE;
    if (surfaceMesh->numElement == 0) return CAPS_BADVALUE;

    // Set the number of surface nodes
    in.numberofpoints = surfaceMesh->numNode;

    // Create surface point array
    in.pointlist = new REAL[in.numberofpoints *3];

    // Transfer input coordinates to tetgen array
    for(i = 0; i < surfaceMesh->numNode; i++) {
        in.pointlist[3*i  ] = surfaceMesh->node[i].xyz[0];
        in.pointlist[3*i+1] = surfaceMesh->node[i].xyz[1];
        in.pointlist[3*i+2] = surfaceMesh->node[i].xyz[2];
    }

    // Set the number of surface tri
    in.numberoffacets = surfaceMesh->numElement;

    // Create surface tri arrays
    in.facetlist = new tetgenio::facet[in.numberoffacets];

    // Create surface marker/BC arrays
    in.facetmarkerlist = new int[in.numberoffacets];

    // Transfer input surfaceMesh->localTriFaceList array to tetgen array
    for(i = 0; i < surfaceMesh->numElement; i++) {
        f =&in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];

        p->numberofvertices = mesh_numMeshElementConnectivity(&surfaceMesh->element[i]);
        p->vertexlist = new int[p->numberofvertices];

        for (j = 0; j < p->numberofvertices; j++) {
            p->vertexlist[j] = surfaceMesh->element[i].connectivity[j];
        }
    }

    // Transfer input BC array to tetgen array
    std::set<int> uniqueMarker;
    for(i = 0; i < surfaceMesh->numElement; i++) {
        in.facetmarkerlist[i] = surfaceMesh->element[i].markerID;
        uniqueMarker.insert(in.facetmarkerlist[i]);
    }

    // If no input string is provided create a simple one based on exposed parameters
    if (meshInput.tetgenInput.meshInputString == NULL) {

        if (meshInput.preserveSurfMesh !=0) strcat(temp, "Y"); // Preserve surface mesh flag

        if ((meshInput.tetgenInput.meshQuality_rad_edge != 0) ||
            (meshInput.tetgenInput.meshQuality_angle != 0)) {

            strcat(temp, "q"); // Mesh quality flag

            if (meshInput.tetgenInput.meshQuality_rad_edge >= 0) {

                sprintf(q, "%.3f",meshInput.tetgenInput.meshQuality_rad_edge);
                strcat(temp,q);

            } else if (meshInput.tetgenInput.meshQuality_rad_edge < 0)  {

                printf("Not setting meshQuality radius-edge ratio. Value needs to be positive\n");
            }

            if (meshInput.tetgenInput.meshQuality_angle >= 0) {

                sprintf(q, "/%.3f",meshInput.tetgenInput.meshQuality_angle);
                strcat(temp,q);

            } else if (meshInput.tetgenInput.meshQuality_angle < 0) {

                printf("Not setting meshQuality dihedral angle. Value needs to be positive\n");
            }

        }

        if ((meshInput.quiet != 0) && (meshInput.tetgenInput.verbose == 0)) strcat(temp, "Q"); // Quiet: No terminal output except errors
        if (meshInput.tetgenInput.verbose != 0) strcat(temp, "V"); // Verbose: Detailed information, more terminal output.

        if (meshInput.tetgenInput.meshTolerance > 0) {
            sprintf(q, "T%.2e",meshInput.tetgenInput.meshTolerance);
            strcat(temp,q);
        }

        if (meshInput.tetgenInput.regions.size > 0
            || meshInput.tetgenInput.holes.size > 0) {
            strcat(temp, "A");
        }

        // Transfer temp character array to input
        inputString = temp;

    } else {

        inputString = meshInput.tetgenInput.meshInputString;
        // Dont think this is needed - maybe p is the default
        // Check to make sure a 'p' is included in the input string provided by a user
        /*if (strchr(input,'p') == NULL) {
            printf("Warning: No \'p\' was found in the Tetgen input string. Input into Tetgen is currently based on PLC's\n");
        }*/
    }


    /*==============================================================*/

    //Create an "empty mesh" where only surface nodes are connected to create the volume
    //Tet centers of the empty mesh will be used to identify holes
    //
    // -p Generate tetrahedra
    // -Y Preserves the input surface mesh (does not modify it).
    tetgenio emptymesh;
    try {
        tetrahedralize((char*)"pYQ", &in, &emptymesh);
    } catch (...){
        printf("Tetgen failed to generate an empty volume mesh......!!!\n");
        printf("  See Tecplot file tetegenDebugSurface.dat for the surface mesh\n");
        mesh_writeTecplot(aimInfo,"tetegenDebugSurface.dat", 1, surfaceMesh, 1.0);
        return -335;
    }

    std::vector<REAL> holepoints;

    // Only solid bodies can have holes
    for ( std::set<int>::const_iterator marker = uniqueMarker.begin(); marker != uniqueMarker.end(); marker++ )
    {
        int globaltri = 0;
        while ( globaltri < in.numberoffacets && in.facetmarkerlist[globaltri] != *marker ) globaltri++;
        tetgenio::facet *f = in.facetlist + globaltri;
        tetgenio::polygon *p = f->polygonlist;

        // Look for two tets attached to a polygon
        int tet = 0;
        int twotets[2][4] = {{-1, -1, -1, -1}, {-1, -1, -1, -1}};
        for ( int k = 0; k < emptymesh.numberoftetrahedra; k++ )
        {
            bool match[3] = {false, false, false};
            for (int v0 = 0; v0 < 3; v0++)
            {
                for (int n0 = 0; n0 < 4; n0++)
                {
                    if ( p->vertexlist[v0] == emptymesh.tetrahedronlist[4*k+n0] )
                    {
                        match[v0] = true;
                        break;
                    }
                }
                //Found a match, save it of. There should be two tets for a face with holes
                if ( match[0] && match[1] && match[2] )
                {
                    for (int n0 = 0; n0 < 4; n0++)
                        twotets[tet][n0] = emptymesh.tetrahedronlist[4*k+n0];
                    tet++;

                }
            }
            if ( tet == 2 )
                break;
        }


        if ( tet == 2 )
        {
            //Found two tets, the one with a postive normal vector to the cell center is a hole

            //Compute the center of the polygon on the surface
            REAL polycenter[3] = {0, 0, 0};
            REAL polynormal[3] = {0, 0, 0};
            REAL polyedge0[3] = {0, 0, 0};
            REAL polyedge1[3] = {0, 0, 0};
            for (int n = 0; n < 3; n++)
            {
                for (int v = 0; v < 3; v++)
                    polycenter[n] += in.pointlist[(p->vertexlist[v]-1)*3 + n];
                polycenter[n] /= 3;

                polyedge0[n] = in.pointlist[(p->vertexlist[1]-1)*3 + n] - in.pointlist[(p->vertexlist[0]-1)*3 + n];
                polyedge1[n] = in.pointlist[(p->vertexlist[2]-1)*3 + n] - in.pointlist[(p->vertexlist[0]-1)*3 + n];
            }
            CROSS(polynormal, polyedge0, polyedge1);

            //Compute the center of the two tetrahedra
            REAL tetcenters[2][3] = {{0, 0, 0}, {0, 0, 0}};
            for (tet = 0; tet < 2; tet++)
                for (int n = 0; n < 3; n++)
                {
                    for (int v = 0; v < 4; v++)
                        tetcenters[tet][n] += emptymesh.pointlist[(twotets[tet][v]-1)*3+n];
                    tetcenters[tet][n] /= 4;
                }

            REAL diffc[3] = {0, 0, 0};
            //Compute the dot product between the normal vector and the vector to the tet-center
            for (tet = 0; tet < 2; tet++)
            {
                for (int n = 0; n < 3; n++)
                    diffc[n] = tetcenters[tet][n] - polycenter[n];

                //Positive dot product means hole
                if ( DOT(diffc, polynormal) > 0 )
                {
                    for (int n = 0; n < 3; n++)
                        holepoints.push_back(tetcenters[tet][n]);
                }
            }
        }
    }

    const tetgenRegionsStruct* regions = &meshInput.tetgenInput.regions;
    const tetgenHolesStruct* holes = &meshInput.tetgenInput.holes;

    if (regions->size > 0 || holes->size > 0)
    {
      in.numberofregions = regions->size;
      if (in.regionlist != NULL)
      {
        delete [] in.regionlist;
        in.regionlist = NULL;
      }
      if (regions->size > 0)
      {
        in.regionlist = new REAL[5 * regions->size];
        for (int n = 0; n < regions->size; ++n)
        {
          in.regionlist[5 * n + 0] = regions->x[n];
          in.regionlist[5 * n + 1] = regions->y[n];
          in.regionlist[5 * n + 2] = regions->z[n];
          in.regionlist[5 * n + 3] = (REAL) regions->attribute[n];
          in.regionlist[5 * n + 4] = regions->volume_constraint[n];
        }
      }

      in.numberofholes = holes->size;
      if (in.holelist != NULL)
      {
        delete [] in.holelist;
        in.holelist = NULL;
      }
      if (holes->size > 0)
      {
        in.holelist = new REAL[3 * holes->size];
        for (int n = 0; n < holes->size; ++n)
        {
          in.holelist[3 * n + 0] = holes->x[n];
          in.holelist[3 * n + 1] = holes->y[n];
          in.holelist[3 * n + 2] = holes->z[n];
        }
      }
    }
    else
    {
      if ( holepoints.size() > 0 )
      {
          in.numberofholes = holepoints.size()/3;
          in.holelist = new REAL[holepoints.size()];
          for ( std::size_t n = 0; n < holepoints.size(); n++ )
              in.holelist[n] = holepoints[n];
/*
          for ( std::size_t n = 0; n < holepoints.size()/3; n++ )
          {
              std::cout << "hole : " << in.holelist[3*n+0] << ", "
                                     << in.holelist[3*n+1] << ", "
                                     << in.holelist[3*n+2] << std::endl;
          }
*/
      }
    }

    /*==============================================================*/


    printf("\nTetgen input string = %s""\n", inputString);

    // Create volume mesh
    try {
        tetrahedralize(inputString, &in, &out);
    } catch (...){
        AIM_ERROR(aimInfo, "Tetgen failed to generate a volume mesh......!!!");
        AIM_ADDLINE(aimInfo, "  See Tecplot file tetegenDebugSurface.dat for the surface mesh");
        mesh_writeTecplot(aimInfo,"tetegenDebugSurface.dat", 1, surfaceMesh, 1.0);
        return -335;
    }

    // Save data
    //in.save_nodes((char *)"TETGEN_Test");
    //in.save_poly((char *) "TETGEN_Test");
    //out.save_faces((char *) "TETGEN_Test");

    // Transfer tetgen mesh structure to genUnstrMesh format
    status = tetgen_to_MeshStruct(&out, volumeMesh);
    AIM_STATUS(aimInfo, status);

    status = writeBinaryUgrid(aimInfo, fileName, &out);
    AIM_STATUS(aimInfo, status);

    // Populate surface mesh
    /* if (meshInput.preserveSurfMesh == (int) false &&
        meshInput.tetgenInput.meshInputString == NULL &&
        meshInput.tetgenInput.ignoreSurfaceExtract == (int) false) {

        printf("Extracting surface mesh from volume\n");

        int nodeCount = 0;
        int node;
        bool found = false;

        int *unique_Nodelist;

        unique_Nodelist = new int [3*out.numberoftrifaces]; // At most 3*numTriFaces

        for (i = 0; i < 3*out.numberoftrifaces; i++) {

            node = out.trifacelist[i];

            found = false;
            for (j = 0; j < nodeCount; j++){

                if (unique_Nodelist[j] == node) {
                    found = true;
                    break;
                }
            }

            if (found == false) {
                unique_Nodelist[nodeCount] = node;
                nodeCount += 1;
            }
        }

        // Populate local node list
        //printf("Num. nodes = %d\n",nodeCount);
        genUnstrMesh->surfaceMesh.numNodes = nodeCount;

        if (genUnstrMesh->surfaceMesh.nodeLocal2Global != NULL) EG_free(genUnstrMesh->surfaceMesh.nodeLocal2Global);

        genUnstrMesh->surfaceMesh.nodeLocal2Global = (int *) EG_alloc(genUnstrMesh->surfaceMesh.numNodes * sizeof(int));
        if (genUnstrMesh->surfaceMesh.nodeLocal2Global == NULL) return EGADS_MALLOC;

        if (genUnstrMesh->surfaceMesh.xyz != NULL) EG_free(genUnstrMesh->surfaceMesh.xyz);

        genUnstrMesh->surfaceMesh.xyz  = (double *) EG_alloc(3*genUnstrMesh->surfaceMesh.numNodes *sizeof(double));
        if (genUnstrMesh->surfaceMesh.xyz == NULL) return EGADS_MALLOC;

        for (i = 0; i < genUnstrMesh->surfaceMesh.numNodes; i++) {
            genUnstrMesh->surfaceMesh.xyz[3*i + 0] = out.pointlist[(unique_Nodelist[i]-1) + 0]; // -1 -> 1 bias of node numbering
            genUnstrMesh->surfaceMesh.xyz[3*i + 1] = out.pointlist[(unique_Nodelist[i]-1) + 1];
            genUnstrMesh->surfaceMesh.xyz[3*i + 2] = out.pointlist[(unique_Nodelist[i]-1) + 2];

            genUnstrMesh->surfaceMesh.nodeLocal2Global[i] = unique_Nodelist[i];
        }

        // Populate local node connectivity
        //genUnstrMesh->surfaceMesh.numTriFace = out.numberoftrifaces;
        if (genUnstrMesh->surfaceMesh.localTriFaceList != NULL) EG_free(genUnstrMesh->surfaceMesh.localTriFaceList);

        genUnstrMesh->surfaceMesh.localTriFaceList = (int *) EG_alloc(3*genUnstrMesh->surfaceMesh.numTriFace * sizeof(int));
        if (genUnstrMesh->surfaceMesh.localTriFaceList == NULL) return EGADS_MALLOC;

        for (i = 0; i < 3*genUnstrMesh->surfaceMesh.numTriFace; i++) {

            for (j = 0; j < genUnstrMesh->surfaceMesh.numNodes; j++) {

                if (out.trifacelist[i] == genUnstrMesh->surfaceMesh.nodeLocal2Global[j]) {
                    genUnstrMesh->surfaceMesh.localTriFaceList[i] = j+1; // Indexing starts at 1 -> 1 bias of node numbering
                    break;
                }
            }
        }

        // Clean up unique node list
        delete [] unique_Nodelist;
        unique_Nodelist = NULL;

    } else if (meshInput.preserveSurfMesh == (int) true &&
               meshInput.tetgenInput.meshInputString == NULL) {
        printf("Surface mesh as been preserved getting global node connectivity and local-2-global node mapping\n");

        if (genUnstrMesh->surfaceMesh.nodeLocal2Global != NULL) EG_free(genUnstrMesh->surfaceMesh.nodeLocal2Global);

        genUnstrMesh->surfaceMesh.nodeLocal2Global = (int *) EG_alloc(genUnstrMesh->surfaceMesh.numNodes * sizeof(int));
        if (genUnstrMesh->surfaceMesh.nodeLocal2Global == NULL) return EGADS_MALLOC;

        //Get local to global node numbering
        for (i = 0; i < genUnstrMesh->surfaceMesh.numNodes; i++) {
            for(j = 0; j < out.numberofpoints; j++){
                if ( genUnstrMesh->surfaceMesh.xyz[3*i + 0] == out.pointlist[3*j + 0] &&
                     genUnstrMesh->surfaceMesh.xyz[3*i + 1] == out.pointlist[3*j + 1] &&
                     genUnstrMesh->surfaceMesh.xyz[3*i + 2] == out.pointlist[3*j + 2]) {
                    break;
                }
            }

            genUnstrMesh->surfaceMesh.nodeLocal2Global[i] = j+in.firstnumber; // Indexing starts at 1 -> 1 bias of node numbering
        }
    }
     */
    // Release TetGen memory - get released when function goes out of scope
    //in.deinitialize();
    //out.deinitialize();

    status = CAPS_SUCCESS;
    printf("Done meshing using TetGen!\n");
cleanup:
    return status;
}
