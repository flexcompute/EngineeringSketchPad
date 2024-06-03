// This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

#include <string.h>
#include <math.h>

#include "aimUtil.h" // Bring in AIM utility functions

#include "miscUtils.h" // Bring in misc. utility functions
#include "meshUtils.h" // Bring in mesh utility functions
#include "vlmUtils.h"  // Bring in vortex lattice method utility functions

#include "arrayUtils.h" // Bring in array utility functions

#include "feaTypes.h"  // Bring in FEA structures
#include "meshTypes.h" // Bring in mesh structures
#include "capsTypes.h" // Bring in CAPS types

#include "feaUtils.h"
#include "jsonUtils.h"

#ifdef WIN32
#define strcasecmp  stricmp
#endif

//#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))



int fea_createMesh(void *aimInfo,
       /*@null@*/  double paramTess[3],                 // (in)  Tessellation parameters
                   int    edgePointMin,                 // (in)  minimum points along any Edge
                   int    edgePointMax,                 // (in)  maximum points along any Edge
                   int    quadMesh,                     // (in)  only do tris-for faces
                   mapAttrToIndexStruct *groupMap,      // (in)  map from CAPSGroup names to indexes
                   mapAttrToIndexStruct *constraintMap, // (in)  map from CAPSConstraint names to indexes
                   mapAttrToIndexStruct *loadMap,       // (in)  map from CAPSLoad names to indexes
                   mapAttrToIndexStruct *transferMap,   // (in)  map from CAPSTransfer names to indexes
                   mapAttrToIndexStruct *connectMap,    // (in)  map from CAPSConnect names to indexes
                   mapAttrToIndexStruct *responseMap,   // (in)  map from CAPSResponse names to indexes
        /*@null@*/ mapAttrToIndexStruct *referenceMap,  // (in)  map from CAPSReference names to indexes
                   int *numMesh,                        // (out) total number of FEA mesh structures
                   meshStruct **feaMesh,                // (out) FEA mesh structure
                   feaProblemStruct *feaProblem ) {     // (out) FEA problem structure

    int status; // Function return status

    int i, body; // Indexing

    // Bodies
    const char *intents;
    int   numBody; // Number of Bodies
    ego  *bodies;
    const char *discipline;
    int stat, nGlobal;

    // Coordinate system
    mapAttrToIndexStruct coordSystemMap, attrMapTemp1, attrMapTemp2, meshMap;

    int meshInd;
    capsValue *meshVal, *meshMorphVal;

    meshStruct *feaMeshes = NULL, tempMesh;
    int *feaMeshList = NULL; // List to seperate structural meshes for aero

    ego tempBody;
    double capsMeshLength = 0, bbox[6], refLen;

    // Inherited fea/volume mesh related variables
    int numFEAMesh = 0;
    int feaMeshInherited = (int) false;

    // Destroy feaMeshes
    if ((*feaMesh) != NULL) {

        for (i = 0; i < *numMesh; i++) {
            status = destroy_meshStruct(&(*feaMesh)[i]);
            if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        AIM_FREE(*feaMesh);
    }

    (*feaMesh) = NULL;
    *numMesh = 0;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if ((numBody <= 0) || (bodies == NULL)) {
        AIM_ERROR(aimInfo, "No Bodies!\n");
        return CAPS_SOURCEERR;
    }

    initiate_meshStruct(&tempMesh);

    // Initiate our maps
    status = initiate_mapAttrToIndexStruct(&coordSystemMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&attrMapTemp1);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&attrMapTemp2);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&meshMap);
    AIM_STATUS(aimInfo, status);

    // Alloc feaMesh list
    AIM_ALLOC(feaMeshList, numBody, int, aimInfo ,status);

    // Set all indexes to true
    for (i = 0; i < numBody; i++ ) feaMeshList[i] = (int) true;

    // Check for capsDiscipline consistency
    i = 0;
    for (body = 0; body < numBody; body++) {

        status = retrieve_CAPSDisciplineAttr(bodies[body], &discipline);
        if (status != CAPS_SUCCESS) continue; // Need to add an error code

        if (strcasecmp(discipline, "structure") != 0) {
            feaMeshList[body] = (int) false;
            i++;
        }
    }

    if (i == numBody) {
        AIM_ERROR(aimInfo, "No bodies with capsDiscipline Structure!");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    //for (body = 0; body< numBody; body++) printf("Body %i, FeaMeshList = %i\n", body, feaMeshList[body]);

    // Get CoordSystem attribute to index mapping
    status = create_CoordSystemAttrToIndexMap(numBody,
                                              bodies,
                                              3, //>2 - search the body, faces, edges, and all the nodes
                                              &coordSystemMap);
    AIM_STATUS(aimInfo, status);

    status = fea_getCoordSystem(numBody,
                                bodies,
                                coordSystemMap,
                                &feaProblem->numCoordSystem,
                                &feaProblem->feaCoordSystem);
    AIM_STATUS(aimInfo, status);

    // Get capsConstraint name and index mapping
    status = create_CAPSConstraintAttrToIndexMap(numBody,
                                                 bodies,
                                                 3, //>2 - search the body, faces, edges, and all the nodes
                                                 constraintMap);
    AIM_STATUS(aimInfo, status);

    // Get capsLoad name and index mapping
    status = create_CAPSLoadAttrToIndexMap(numBody,
                                           bodies,
                                           3, //>2 - search the body, faces, edges, and all the nodes
                                           loadMap);
    AIM_STATUS(aimInfo, status);

    // Get transfer to index mapping
    status = create_CAPSBoundAttrToIndexMap(numBody,
                                            bodies,
                                            3, //>2 - search the body, faces, edges, and all the nodes
                                            transferMap);
    AIM_STATUS(aimInfo, status);

    if (connectMap != NULL) {
        // Get connect to index mapping
        status = create_CAPSConnectAttrToIndexMap(numBody,
                                                  bodies,
                                                  3, //>2 - search the body, faces, edges, and all the nodes
                                                  connectMap);
        AIM_STATUS(aimInfo, status);
    }

    if (responseMap != NULL) {
        // Get response to index mapping
        status = create_CAPSResponseAttrToIndexMap(numBody,
                                                  bodies,
                                                  3, //>2 - search the body, faces, edges, and all the nodes
                                                  responseMap);
        AIM_STATUS(aimInfo, status);
    }

    if (referenceMap != NULL) {
        // Get response to index mapping
        status = create_CAPSReferenceAttrToIndexMap(numBody,
                                                  bodies,
                                                  3, //>2 - search the body, faces, edges, and all the nodes
                                                  referenceMap);
        AIM_STATUS(aimInfo, status);
    }

    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            3, //>2 - search the body, faces, edges, and all the nodes
                                            groupMap);
    AIM_STATUS(aimInfo, status);

    // Get the Mesh input Value
    meshInd = aim_getIndex(aimInfo, "Mesh", ANALYSISIN);
    if (meshInd < 1)
        meshInd = aim_getIndex(aimInfo, "Surface_Mesh", ANALYSISIN);
    if (meshInd < 1) {
        AIM_ERROR(aimInfo, "Developer error: No 'Mesh' or 'Surface_Mesh' ANALYSISIN Index!");
        status = CAPS_BADINDEX;
        goto cleanup;
    }

    status = aim_getValue(aimInfo, meshInd, ANALYSISIN, &meshVal);
    AIM_STATUS(aimInfo, status);

    if (meshVal->type != PointerMesh) {
        AIM_ERROR(aimInfo, "Developer error! Mesh is not a PointerMesh!");
        status = CAPS_NOTIMPLEMENT;
        goto cleanup;
    }

    // Get the Mesh_Morph input Value
    meshInd = aim_getIndex(aimInfo, "Mesh_Morph", ANALYSISIN);
    if (meshInd < 1) {
        AIM_ERROR(aimInfo, "Developer error: No 'Mesh_Morph' ANALYSISIN Index!");
        status = CAPS_BADINDEX;
        goto cleanup;
    }

    status = aim_getValue(aimInfo, meshInd, ANALYSISIN, &meshMorphVal);
    AIM_STATUS(aimInfo, status);


    feaMeshInherited = (int) false;

    if ((meshVal->nullVal == NotNull) ||
        (meshMorphVal->vals.integer == (int) true &&
         meshVal->nullVal == IsNull)) {

        if ( meshMorphVal->vals.integer == (int) true &&
             meshVal->nullVal == IsNull ) { // If we are mighty morphing

            if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
                feaProblem->meshRefIn == NULL ||
                feaProblem->meshRefIn != &feaProblem->meshRefObj) {

              // Lets "load" the meshRef now since it's not linked
              status = aim_loadMeshRef(aimInfo, &feaProblem->meshRefObj);
              AIM_STATUS(aimInfo, status);

              // Mighty Morph the mesh
              status = aim_morphMeshUpdate(aimInfo, &feaProblem->meshRefObj, numBody, bodies);
              AIM_STATUS(aimInfo, status);
            }
            /*@-immediatetrans@*/
            feaProblem->meshRefIn = &feaProblem->meshRefObj;
            /*@+immediatetrans@*/
        } else {
            feaProblem->meshRefIn = (aimMeshRef *) meshVal->vals.AIMptr;
        }

        numFEAMesh = feaProblem->meshRefIn->nmap;

        // See if a FEA mesh is available from parent
        if (feaProblem->meshRefIn->type == aimAreaMesh ||
            feaProblem->meshRefIn->type == aimSurfaceMesh) {

            if (numFEAMesh != numBody) { // May not be an error if we are doing aero-struct

                // Check for capsDiscipline consistency
                for (body = 0; body < numFEAMesh; body++) {
                                                                                         // Dummy arguments
                    status = EG_statusTessBody(feaProblem->meshRefIn->maps[body].tess, &tempBody, &stat, &nGlobal);
                    if (status != EGADS_SUCCESS) goto cleanup;

                    status = retrieve_CAPSDisciplineAttr(tempBody, &discipline);
                    if (status != EGADS_SUCCESS) {
                        AIM_ERROR  (aimInfo, "Failed to find a capsDiscipline attribute!\n");
                        AIM_ADDLINE(aimInfo, "Number of linked surface meshes does not match the number of bodies, this is only allowed if doing aero-struct analysis\n");
                        status= CAPS_SOURCEERR;
                        goto cleanup;
                    }

                    if (strcasecmp(discipline, "structure") != 0) {
                        AIM_ERROR  (aimInfo, "Failed to find a capsDiscipline attribute - 'structure'!\n");
                        AIM_ADDLINE(aimInfo, "Number of linked surface meshes does not match the number of bodies, this is only allowed if doing aero-struct analysis\n");
                        status= CAPS_SOURCEERR;
                        goto cleanup;
                    }

                    // need a check to make sure all tempMesh->groupMap are identical
                }

                // We need to update our capsGroup attribute map

                // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
                status = create_CAPSGroupAttrToIndexMap(numBody,
                                                        bodies,
                                                        3, //>2 - search the body, faces, edges, and all the nodes
                                                        &attrMapTemp2);
                AIM_STATUS(aimInfo, status);

                // Get attribute to index mapping
                status = create_MeshRefToIndexMap(aimInfo, feaProblem->meshRefIn, &attrMapTemp1);
                AIM_STATUS(aimInfo, status);

                status = merge_mapAttrToIndexStruct(&attrMapTemp1, &attrMapTemp2, groupMap);
                AIM_STATUS(aimInfo, status);
            }

            AIM_ALLOC(feaMeshes, numFEAMesh, meshStruct, aimInfo, status);
            for (body = 0; body < numFEAMesh; body++) (void) initiate_meshStruct(&feaMeshes[body]);

            for (body = 0; body < numFEAMesh; body++) {

                status = copy_mapAttrToIndexStruct( groupMap,
                                                    &tempMesh.groupMap );
                AIM_STATUS(aimInfo, status);

                tempMesh.egadsTess = feaProblem->meshRefIn->maps[body].tess;
                status = mesh_surfaceMeshEGADSTess(aimInfo, &tempMesh, (int)(feaProblem->meshRefIn->type == aimAreaMesh));
                AIM_STATUS(aimInfo, status);

                // Create a new mesh with topology tagged with capsIgnore being removed, if capsIgnore isn't found the mesh is simply copied.
                status = mesh_createIgnoreMesh(aimInfo, &tempMesh, &feaMeshes[body]);
                AIM_STATUS(aimInfo, status);

                destroy_meshStruct(&tempMesh);

                // Change the analysis type of the mesh
                status = mesh_setAnalysisType(MeshStructure, &feaMeshes[body]);
                AIM_STATUS(aimInfo, status);


                // Get FEA Problem from EGADs body
                status = fea_setAnalysisData(aimInfo,
                                             groupMap,
                                             &coordSystemMap,
                                             constraintMap,
                                             loadMap,
                                             transferMap,
                                             connectMap,
                                             responseMap,
                                             referenceMap,
                                             &feaMeshes[body]);
                AIM_STATUS(aimInfo, status);

                status = mesh_fillQuickRefList( aimInfo, &feaMeshes[body] );
                AIM_STATUS(aimInfo, status);

                printf("\tMesh for body = %d\n", body);
                printf("\tNumber of nodal coordinates = %d\n", feaMeshes[body].numNode);
                printf("\tNumber of elements = %d\n", feaMeshes[body].numElement);
                printf("\tElemental Nodes = %d\n", feaMeshes[body].meshQuickRef.numNode);
                printf("\tElemental Rods  = %d\n", feaMeshes[body].meshQuickRef.numLine);
                printf("\tElemental Tria3 = %d\n", feaMeshes[body].meshQuickRef.numTriangle);
                printf("\tElemental Quad4 = %d\n", feaMeshes[body].meshQuickRef.numQuadrilateral);
            }


            if (numFEAMesh > 1) {
                printf("Combining multiple FEA meshes!\n");
            }

            status = mesh_combineMeshStruct(aimInfo,
                                            numFEAMesh,
                                            feaMeshes,
                                            &feaProblem->feaMesh);
            AIM_STATUS(aimInfo, status);

            if (numFEAMesh > 1) {
                printf("\tCombined Number of nodal coordinates = %d\n", feaProblem->feaMesh.numNode);
                printf("\tCombined Number of elements = %d\n", feaProblem->feaMesh.numElement);
                printf("\tCombined Elemental Nodes = %d\n", feaProblem->feaMesh.meshQuickRef.numNode);
                printf("\tCombined Elemental Rods  = %d\n", feaProblem->feaMesh.meshQuickRef.numLine);
                printf("\tCombined Elemental Tria3 = %d\n", feaProblem->feaMesh.meshQuickRef.numTriangle);
                printf("\tCombined Elemental Quad4 = %d\n", feaProblem->feaMesh.meshQuickRef.numQuadrilateral);
            }

            // Set output meshes
            *numMesh = numFEAMesh;
            *feaMesh = feaMeshes;
            feaMeshes = NULL;

            // Set reference meshes
            feaProblem->feaMesh.numReferenceMesh = *numMesh;
            feaProblem->feaMesh.referenceMesh = (meshStruct *) EG_alloc(*numMesh*sizeof(meshStruct));
            if (feaProblem->feaMesh.referenceMesh == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (body = 0; body < *numMesh; body++) {
                feaProblem->feaMesh.referenceMesh[body] = (*feaMesh)[body];
            }

        } else { // Check to see if a general unstructured volume mesh is available

            printf("Found link for a  volume mesh (Volume_Mesh) from parent\n");

            numFEAMesh = 1;
            if (numFEAMesh != 1) {
                AIM_ERROR(aimInfo, "Can not accept multiple volume meshes");
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            if (numFEAMesh != numBody) {
                printf("Number of inherited volume meshes does not match the number of bodies - assuming volume mesh is already combined\n");
            }

            AIM_ERROR(aimInfo, "Volume meshes not yet supported for structural analysis");
            status = CAPS_BADVALUE;
            goto cleanup;

#if 0
            tempMesh = (meshStruct *) meshVal->vals.AIMptr;
            status = mesh_copyMeshStruct( tempMesh, &feaProblem->feaMesh);
            AIM_STATUS(aimInfo, status);

            //feaProblem->feaMesh.egadsTess =
            status = mesh_surfaceMeshEGADSTess(aimInfo, &feaProblem->feaMesh.egadsTess, (int)false);
            AIM_STATUS(aimInfo, status);


            // Set reference meshes
             feaProblem->feaMesh.numReferenceMesh = tempMesh->numReferenceMesh;
             feaProblem->feaMesh.referenceMesh = (meshStruct *) EG_alloc(tempMesh->numReferenceMesh*sizeof(meshStruct));
             if (feaProblem->feaMesh.referenceMesh == NULL) {
                 status = EGADS_MALLOC;
                 goto cleanup;
             }

             for (body = 0; body < tempMesh->numReferenceMesh; body++) {
                 feaProblem->feaMesh.referenceMesh[body] = tempMesh->referenceMesh[body];
             }

            for (i = 0; i < feaProblem->feaMesh.numReferenceMesh; i++) {

                status = aim_newTess(aimInfo,
                                     feaProblem->feaMesh.referenceMesh[i].egadsTess);
                AIM_STATUS(aimInfo, status);
            }

            // Update/change the analysis data in a meshStruct
            status = change_meshAnalysis(&feaProblem->feaMesh, MeshStructure);
            AIM_STATUS(aimInfo, status);

            // Get FEA Problem from EGADs body
            status = fea_setAnalysisData(aimInfo,
                                         groupMap,
                                         &coordSystemMap,
                                         constraintMap,
                                         loadMap,
                                         transferMap,
                                         connectMap,
                                         responseMap,
                                         referenceMap,
                                         &feaProblem->feaMesh);
            AIM_STATUS(aimInfo, status);
#endif
        }

        feaMeshInherited = (int) true;
    }

    // If we didn't inherit a FEA mesh we need to get one ourselves
    if (feaMeshInherited == (int) false) {

        status = check_CAPSMeshLength(numBody, bodies, &capsMeshLength);

        if (status == CAPS_NOTFOUND) capsMeshLength = -1;
        else AIM_STATUS(aimInfo, status);

        /*
        if (capsMeshLength <= 0 || status != CAPS_SUCCESS) {
          printf("**********************************************************\n");
          if (status != CAPS_SUCCESS)
            printf("capsMeshLength is not set on any body.\n");
          else
            printf("capsMeshLength: %f\n", capsMeshLength);
          printf("\n");
          printf("The capsMeshLength attribute must\n"
                 "present on at least one body.\n"
                 "\n"
                 "capsMeshLength should be a a positive value representative\n"
                 "of a characteristic length of the geometry,\n"
                 "e.g. the MAC of a wing or diameter of a fuselage.\n");
          printf("**********************************************************\n");
          status = CAPS_BADVALUE;
          goto cleanup;
        }*/

        if (paramTess == NULL) {
          AIM_ERROR(aimInfo, "Developer error paramTess == NULL");
          status = CAPS_BADVALUE;
          goto cleanup;
        }

        if (edgePointMin < 2) {
            printf("The minimum number of allowable edge points is 2 not %d\n", edgePointMin);
            edgePointMin = 2;
        }

        if (edgePointMax < edgePointMin) {
            printf("The maximum number of edge points must be greater than the current minimum (%d)\n", edgePointMin);
            edgePointMax = edgePointMin+1;
        }

        if (capsMeshLength <= 0) {
          refLen = 0;
          for (body = 0; body < numBody; body++) {
            if (feaMeshList[body] != (int) true) continue;

            status = EG_getBoundingBox(bodies[body], bbox);
            AIM_STATUS(aimInfo, status);

            refLen = MAX(refLen, sqrt( (bbox[3] - bbox[0]) * (bbox[3] - bbox[0])
                                      +(bbox[4] - bbox[1]) * (bbox[4] - bbox[1])
                                      +(bbox[5] - bbox[2]) * (bbox[5] - bbox[2])));
          }
       } else {
         refLen = capsMeshLength;
       }

        // Modify the EGADS body tessellation based on given inputs
    /*@-nullpass@*/
        status =  mesh_modifyBodyTess(0,
                                      NULL,
                                      edgePointMin,
                                      edgePointMax,
                                      quadMesh,
                                      &refLen,
                                      paramTess,
                                      meshMap,
                                      numBody,
                                      bodies);
    /*@+nullpass@*/
        AIM_STATUS(aimInfo, status);

        for (body = 0; body < numBody; body++) {
            if (feaMeshList[body] != (int) true) continue;

            status = copy_mapAttrToIndexStruct( groupMap,
                                                &tempMesh.groupMap );
            AIM_STATUS(aimInfo, status);

            status = mesh_surfaceMeshEGADSBody(aimInfo,
                                               bodies[body],
                                               refLen,
                                               paramTess,
                                               quadMesh,
                                               &tempMesh.egadsTess);
            AIM_STATUS(aimInfo, status, "Problem during surface meshing of body %d", body+1);
            AIM_NOTNULL(tempMesh.egadsTess, aimInfo, status);

            status = mesh_surfaceMeshEGADSTess(aimInfo, &tempMesh, (int)false);
            AIM_STATUS(aimInfo, status);

            *numMesh += 1;
            AIM_REALL(*feaMesh, *numMesh, meshStruct, aimInfo, status);

            status = initiate_meshStruct(&(*feaMesh)[*numMesh-1]);
            AIM_STATUS(aimInfo, status);

            // Create a new mesh with topology tagged with capsIgnore being removed, if capsIgnore isn't found the mesh is simply copied.
            status = mesh_createIgnoreMesh(aimInfo, &tempMesh, &(*feaMesh)[*numMesh-1]);
            AIM_STATUS(aimInfo, status);

            destroy_meshStruct(&tempMesh);

            // Change the analysis type of the mesh
            status = mesh_setAnalysisType(MeshStructure, &(*feaMesh)[*numMesh-1]);
            AIM_STATUS(aimInfo, status);

            // Get FEA Problem from EGADs body
            status = fea_setAnalysisData(aimInfo,
                                         groupMap,
                                         &coordSystemMap,
                                         constraintMap,
                                         loadMap,
                                         transferMap,
                                         connectMap,
                                         responseMap,
                                         referenceMap,
                                         &(*feaMesh)[*numMesh-1]);
            AIM_STATUS(aimInfo, status);

            status = mesh_fillQuickRefList( aimInfo, &(*feaMesh)[*numMesh-1]);
            AIM_STATUS(aimInfo, status);


#if 0
            // Get FEA Problem from EGADs body
            status = fea_bodyToBEM(aimInfo,
                                   bodies[body], // (in)  EGADS Body
                                   paramTess,    // (in)  Tessellation parameters
                                   edgePointMin, // (in)  minimum points along any Edge
                                   edgePointMax, // (in)  maximum points along any Edge
                                   quadMesh,
                                   groupMap,
                                   &coordSystemMap,
                                   constraintMap,
                                   loadMap,
                                   transferMap,
                                   connectMap,
                                   responseMap,
                                   referenceMap,
                                   &(*feaMesh)[*numMesh-1]);
#endif
            AIM_STATUS(aimInfo, status);
            printf("\tMesh for body = %d\n", body);
            printf("\tNumber of nodal coordinates = %d\n", (*feaMesh)[*numMesh-1].numNode);
            printf("\tNumber of elements = %d\n", (*feaMesh)[*numMesh-1].numElement);
            printf("\tElemental Nodes = %d\n", (*feaMesh)[*numMesh-1].meshQuickRef.numNode);
            printf("\tElemental Rods  = %d\n", (*feaMesh)[*numMesh-1].meshQuickRef.numLine);
            printf("\tElemental Tria3 = %d\n", (*feaMesh)[*numMesh-1].meshQuickRef.numTriangle);
            printf("\tElemental Quad4 = %d\n", (*feaMesh)[*numMesh-1].meshQuickRef.numQuadrilateral);

            // set the resulting tessellation
            status = aim_newTess(aimInfo, (*feaMesh)[*numMesh-1].egadsTess);
            AIM_STATUS(aimInfo, status);
        }

        // Only compbine if there are actual meshes
        if (*numMesh > 0) {
            if (*numMesh > 1) printf("Combining multiple FEA meshes!\n");
            AIM_NOTNULL((*feaMesh), aimInfo, status);

            // Combine fea meshes into a single mesh for the problem
            status = mesh_combineMeshStruct(aimInfo,
                                            *numMesh,
                                            (*feaMesh),
                                            &feaProblem->feaMesh);
            AIM_STATUS(aimInfo, status);

            if (*numMesh > 1) {
                printf("\tCombined Number of nodal coordinates = %d\n", feaProblem->feaMesh.numNode);
                printf("\tCombined Number of elements = %d\n", feaProblem->feaMesh.numElement);
                printf("\tCombined Elemental Nodes = %d\n", feaProblem->feaMesh.meshQuickRef.numNode);
                printf("\tCombined Elemental Rods  = %d\n", feaProblem->feaMesh.meshQuickRef.numLine);
                printf("\tCombined Elemental Tria3 = %d\n", feaProblem->feaMesh.meshQuickRef.numTriangle);
                printf("\tCombined Elemental Quad4 = %d\n", feaProblem->feaMesh.meshQuickRef.numQuadrilateral);
            }

            // Set reference meshes
            feaProblem->feaMesh.numReferenceMesh = *numMesh;
            feaProblem->feaMesh.referenceMesh = (meshStruct *) EG_alloc(*numMesh*sizeof(meshStruct));
            if (feaProblem->feaMesh.referenceMesh == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (body = 0; body < *numMesh; body++) {
                feaProblem->feaMesh.referenceMesh[body] = (*feaMesh)[body];
            }
        }
    }
    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(feaMeshes);
    AIM_FREE(feaMeshList);

    (void ) destroy_meshStruct(&tempMesh);

    (void ) destroy_mapAttrToIndexStruct(&coordSystemMap);

    (void ) destroy_mapAttrToIndexStruct(&attrMapTemp1);
    (void ) destroy_mapAttrToIndexStruct(&attrMapTemp2);

    return status;
}

// Convert an EGADS body to a boundary element model, modified by Ryan Durscher (AFRL)
// from code originally written by John Dannenhoffer @ Syracuse University, patterned after code
// written by Bob Haimes  @ MIT
int fea_bodyToBEM(void *aimInfo,                       // (in)  AIM structure
                  ego    ebody,                        // (in)  EGADS Body
                  double paramTess[3],                 // (in)  Tessellation parameters
                  int    edgePointMin,                 // (in)  minimum points along any Edge
                  int    edgePointMax,                 // (in)  maximum points along any Edge
                  int    quadMesh,                     // (in)  0 - only do tris-for faces, 1 - mixed quad/tri, 2 - regularized quads
                  mapAttrToIndexStruct *attrMap,       // (in)  map from CAPSGroup names to indexes
                  mapAttrToIndexStruct *coordSystemMap,// (in)  map from CoordSystem names to indexes
                  mapAttrToIndexStruct *constraintMap, // (in)  map from CAPSConstraint names to indexes
                  mapAttrToIndexStruct *loadMap,       // (in)  map from CAPSLoad names to indexes
                  mapAttrToIndexStruct *transferMap,   // (in)  map from CAPSTransfer names to indexes
                  mapAttrToIndexStruct *connectMap,    // (in)  map from CAPSConnect names to indexes
                  mapAttrToIndexStruct *responseMap,   // (in)  map from CAPSResponse names to indexes
                  mapAttrToIndexStruct *referenceMap,  // (in)  map from CAPSReference names to indexes
                  meshStruct *feaMesh)                 // (out) FEA mesh structure
{
    int status = 0; // Function return status

    int isNodeBody;
    int i, j, k, face, edge, patch; // Indexing
    double scale;

    // Body entities
    int numNode = 0, numEdge = 0, numFace = 0;
    ego *enodes=NULL, *eedges=NULL, *efaces=NULL;

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char *string;

    // Meshing
    int numPoint = 0, numTri = 0;
    const int  *triConn = NULL, *triNeighbor = NULL; // Triangle

    int numPatch = 0; // Patching
    int n1, n2, *qints = NULL;

    int gID; // Global id

    const double *xyz, *uv;

    const int *ptype = NULL, *pindex = NULL, *pvindex = NULL, *pbounds = NULL;

    int       periodic, nchange, oclass, mtype, nchild, nchild2, *senses;

    // Edge point distributions
    int    *points=NULL, *isouth=NULL, *ieast=NULL, *inorth=NULL, *iwest=NULL;
    double  params[3], bbox[6], size, range[2], arclen, data[4], eval[18], eval2[18];
    double *rpos=NULL;
    ego     eref, *echilds, *echilds2, eloop, tempBody, topObj, prev, next;

    int bodySubType = 0; // Body classification

    int pointType, pointTopoIndex;
    double xyzPoint[3], xyzNode[3];

    // Attributues
    const char *attrName;
    int         attrIndex, coordSystemIndex, loadIndex;

    int numElement = 0; // Number of elements
    int coordID = 0; // Default coordinate id for mesh

    feaMeshDataStruct *feaData;

    meshElementStruct *tempElement = NULL;

    // In case our geometry has ignores
    int ignoreFound = (int) false;

    // ---------------------------------------------------------------

    printf("Creating FEA BEM\n");

    // Check for contradiction where quading is requested but dissabled on the body
    if (quadMesh == (int)true) {
        status = EG_attributeRet(ebody, ".qParams", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && (atype != ATTRREAL || (atype == ATTRREAL && reals[0] <= 0 ))) {
            printf("\tQuading on all faces disabled with .qParams attribute on the body\n");
            quadMesh = (int) false;
        }
    }

    // Get number of Nodes, Edges, and Faces in ebody
    status = EG_getBodyTopos(ebody, NULL, NODE, &numNode, &enodes);
    if (status < EGADS_SUCCESS) goto cleanup;
    AIM_NOTNULL(enodes, aimInfo, status);

    status = EG_getBodyTopos(ebody, NULL, EDGE, &numEdge, &eedges);
    if (status != EGADS_SUCCESS) goto cleanup;

    status = EG_getBodyTopos(ebody, NULL, FACE, &numFace, &efaces);
    if (status < EGADS_SUCCESS) goto cleanup;

    // What type of BODY do we have?
    isNodeBody = aim_isNodeBody(ebody, xyzNode);
    if (isNodeBody < EGADS_SUCCESS) goto cleanup;
    if (isNodeBody == EGADS_SUCCESS) {
      // all attributes are on the body rather than the node for a node body
      enodes[0] = ebody;
    }

    // Determine the nominal number of points along each Edge
    points = (int    *) EG_alloc((numEdge+1)     *sizeof(int   ));
    rpos  = (double *) EG_alloc((edgePointMax)*sizeof(double));
    if (points == NULL || rpos == NULL) {
        status = EGADS_MALLOC;
        printf("\tError in fea_bodyToBEM: EG_alloc\n");
        goto cleanup;
    }

    status = EG_getBoundingBox(ebody, bbox);
    if (status < EGADS_SUCCESS) {
        printf("\tError in fea_bodyToBEM: EG_getBoundingBox\n");
        goto cleanup;
    }

    size = sqrt( (bbox[3] - bbox[0]) * (bbox[3] - bbox[0])
                 +(bbox[4] - bbox[1]) * (bbox[4] - bbox[1])
                 +(bbox[5] - bbox[2]) * (bbox[5] - bbox[2]));

    params[0] = paramTess[0] * size;
    params[1] = paramTess[1] * size;
    params[2] = paramTess[2];

    status = EG_attributeAdd(ebody, ".tParam", ATTRREAL, 3, NULL, params, NULL);
    if (status < EGADS_SUCCESS) {
        printf("\tError in fea_bodyToBEM: EG_attributeAdd\n");
        goto cleanup;
    }

    if (isNodeBody == EGADS_SUCCESS)
        params[0] = 1; // does not matter but can't be zero

    if (params[0] <= 0.0) {
        printf("\tError in fea_bodyToBEM: params[0] = %f must be a positive number!\n", params[0]);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    for (i = 1; i <= numEdge; i++) {
        status = EG_getRange(eedges[i-1], range, &periodic);
        if (status < EGADS_SUCCESS) {
            printf("\tError in fea_bodyToBEM: EG_getRange\n");
            goto cleanup;
        }

        status = EG_arcLength(eedges[i-1], range[0], range[1], &arclen);
        if (status < EGADS_SUCCESS) {
            printf("\tError in fea_bodyToBEM: EG_arcLength\n");
            goto cleanup;
        }

        //points[i] = MIN(MAX(MAX(edgePointMin,2), (int)(1+arclen/params[0])), edgePointMax);
        points[i] = (int) min_DoubleVal(
                                        max_DoubleVal(
                                                max_DoubleVal( (double) edgePointMin, 2.0),
                                                (double) (1+arclen/params[0])),
                                        (double) edgePointMax);
    }


    // make arrays for "opposite" sides of four-sided Faces (with only one loop)
    isouth = (int *) EG_alloc((numFace+1)*sizeof(int));
    ieast  = (int *) EG_alloc((numFace+1)*sizeof(int));
    inorth = (int *) EG_alloc((numFace+1)*sizeof(int));
    iwest  = (int *) EG_alloc((numFace+1)*sizeof(int));

    if (isouth == NULL ||
        ieast == NULL  ||
        inorth == NULL ||
        iwest == NULL   ) {

        status = EGADS_MALLOC;
        printf("\tError in fea_bodyToBEM: EG_alloc\n");
        goto cleanup;
    }

    for (i = 1; i <= numFace; i++) {
        isouth[i] = 0;
        ieast [i] = 0;
        inorth[i] = 0;
        iwest [i] = 0;

        // nothing to check if quading isn't requested
        if (quadMesh == (int) false)
            continue;

        // check if quading is disabled with .qParams
        status = EG_attributeRet(efaces[i-1], ".qParams", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && (atype != ATTRREAL || (atype == ATTRREAL && reals[0] <= 0 ))) {
            printf("\tFace %d quading disabled with attribute .qParams\n", i);
            continue;
        }

        // quading only works with one loop
        status = EG_getTopology(efaces[i-1], &eref, &oclass, &mtype, data, &nchild, &echilds, &senses);
        if (status < EGADS_SUCCESS) goto cleanup;

        if (nchild != 1) continue;

        // quading only works if the loop has 4 edges
        eloop = echilds[0];
        status = EG_getTopology(eloop, &eref, &oclass, &mtype, data, &nchild, &echilds, &senses);
        if (status < EGADS_SUCCESS) goto cleanup;

        if (nchild != 4) continue;

        // Check to see if two "straight" edges next to each other are parallel - Don't Quad if so
        for (j = 0; j < 4; j++) {

            status = EG_getTopology(echilds[j], &eref, &oclass, &mtype, data, &nchild2, &echilds2, &senses);
            if (mtype == DEGENERATE) { status = EGADS_DEGEN; break; }
            if (status < EGADS_SUCCESS) goto cleanup;

            if (j < 3) k = j+1;
            else k = 0;

            status = EG_getTopology(echilds[k], &eref, &oclass, &mtype, range, &nchild2, &echilds2, &senses);
            if (mtype == DEGENERATE) { status = EGADS_DEGEN; break; }
            if (status < EGADS_SUCCESS) goto cleanup;

            status = EG_evaluate(echilds[j], data, eval);
            if (status < EGADS_SUCCESS) goto cleanup;

            status = EG_evaluate(echilds[k], range, eval2);
            if (status < EGADS_SUCCESS) goto cleanup;

            scale = dot_DoubleVal(&eval[3], &eval[3]);
            eval[3] /= scale;
            eval[4] /= scale;
            eval[5] /= scale;

            scale = dot_DoubleVal(&eval2[3], &eval2[3]);
            eval2[3] /= scale;
            eval2[4] /= scale;
            eval2[5] /= scale;

            if (fabs(fabs(dot_DoubleVal(&eval[3], &eval2[3])) - 1) < 1E-6) {
                status = EGADS_OUTSIDE;
                break;
            }
        }

        if (status == EGADS_OUTSIDE) {
            if (quadMesh == (int)true) {
                printf("Face %d has parallel edges - not quading\n", i);
            }
            continue;
        }

        if (status == EGADS_DEGEN) {
            if (quadMesh == (int)true) {
                printf("Face %d has a degenerate edge - not quading\n", i);
            }
            continue;
        }

        isouth[i] = status = EG_indexBodyTopo(ebody, echilds[0]);
        if (status < EGADS_SUCCESS) goto cleanup;

        ieast[i]  = status = EG_indexBodyTopo(ebody, echilds[1]);
        if (status < EGADS_SUCCESS) goto cleanup;

        inorth[i] = status = EG_indexBodyTopo(ebody, echilds[2]);
        if (status < EGADS_SUCCESS) goto cleanup;

        iwest[i]  = status = EG_indexBodyTopo(ebody, echilds[3]);
        if (status < EGADS_SUCCESS) goto cleanup;
    }

    // make "opposite" sides of four-sided Faces (with only one loop) match
    nchange = 1;
    for (i = 0; i < MAX(10*numFace,20); i++) {
        nchange = 0;

        for (face = 1; face <= numFace; face++) {
            if (isouth[face] <= 0 || ieast[face] <= 0 ||
                inorth[face] <= 0 || iwest[face] <= 0   ) continue;

            if (points[iwest[face]] < points[ieast[face]]) {
                points[iwest[face]] = points[ieast[face]];
                nchange++;

            } else if (points[ieast[face]] < points[iwest[face]]) {
                       points[ieast[face]] = points[iwest[face]];
                nchange++;
            }

            if (points[isouth[face]] < points[inorth[face]]) {
                points[isouth[face]] = points[inorth[face]];
                nchange++;
            } else if (points[inorth[face]] < points[isouth[face]]) {
                       points[inorth[face]] = points[isouth[face]];
                nchange++;
            }
        }
        if (nchange == 0) break;
    }
    if (nchange > 0) {
        printf("Exceeded number of tries making \"opposite\" sides of four-sided Faces (with only one loop) match\n");
        status = CAPS_MISMATCH;
        goto cleanup;
    }

    // mark the Edges with points[iedge] evenly-spaced points
    for (edge = 1; edge <= numEdge; edge++) {
        for (i = 1; i < points[edge]-1; i++) {
            rpos[i-1] = (double)(i) / (double)(points[edge]-1);
        }

        if (points[edge] == 2) {
            i = 0;
            status = EG_attributeAdd(eedges[edge-1], ".rPos", ATTRINT, 1, &i, NULL, NULL);
            if (status < EGADS_SUCCESS) goto cleanup;

        } else {
            status = EG_attributeAdd(eedges[edge-1], ".rPos", ATTRREAL, points[edge]-2, NULL, rpos, NULL);
            if (status < EGADS_SUCCESS) goto cleanup;
        }
    }

    // Make tessellation
    status = EG_makeTessBody(ebody, params, &feaMesh->egadsTess);
    if (status != EGADS_SUCCESS) {
      printf("\tError in fea_bodyToBEM: EG_makeTessBody\n");
      goto cleanup;
    }

    // Make Quads on each four-sided Face
    params[0] = 0;
    params[1] = 0;
    params[2] = 0;

    // If making quads on faces lets setup an array to keep track of which faces have been quaded.
    if (quadMesh == (int)true) {
        if( numFace > 0) {
            AIM_ALLOC(qints, numFace, int, aimInfo, status);
        }
        // Set default to 0
        for (face = 0; face < numFace; face++) qints[face] = 0;
    }

    if (quadMesh == (int)true) {
        for (face = 1; face <= numFace; face++) {
            if (iwest[face] <= 0) continue;

            status = EG_makeQuads(feaMesh->egadsTess, params, face);
            if (status < EGADS_SUCCESS) {
                printf("Face = %d, failed to make quads\n", face);
                printf("Edges = %d %d %d %d\n", inorth[face], ieast[face], isouth[face], iwest[face]);
                continue;
                //goto cleanup;
            }
        }
    }

    // Set the mesh type information
    feaMesh->meshType = SurfaceMesh;
    feaMesh->analysisType = MeshStructure;

    // Get number of point in the tessellation
    status = EG_statusTessBody(feaMesh->egadsTess, &tempBody, &i, &feaMesh->numNode);
    if (status != EGADS_SUCCESS) goto cleanup;

    feaMesh->node = (meshNodeStruct *) EG_alloc(feaMesh->numNode*sizeof(meshNodeStruct));
    if (feaMesh->node == NULL) {
        feaMesh->numNode = 0;
        status = EGADS_MALLOC;
        goto cleanup;
    }

    for (i = 0; i < feaMesh->numNode; i++) {
        status = initiate_meshNodeStruct(&feaMesh->node[i], feaMesh->analysisType);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Fill up the Attributes for the nodes
    for (i = 0; i < feaMesh->numNode; i++) {

        status = EG_getGlobal(feaMesh->egadsTess, i+1, &pointType, &pointTopoIndex, xyzPoint);
        if (status != EGADS_SUCCESS) goto cleanup;

        feaMesh->node[i].xyz[0] = xyzPoint[0];
        feaMesh->node[i].xyz[1] = xyzPoint[1];
        feaMesh->node[i].xyz[2] = xyzPoint[2];

        feaMesh->node[i].nodeID = i+1;

        feaData = (feaMeshDataStruct *) feaMesh->node[i].analysisData;

        status = fea_setFEADataPoint(efaces, eedges, enodes,
                                     attrMap,
                                     coordSystemMap,
                                     constraintMap,
                                     loadMap,
                                     transferMap,
                                     connectMap,
                                     responseMap,
                                     referenceMap,
                                     pointType, pointTopoIndex,
                                     feaData); // Set the feaData structure
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Fill element information

    // If body is just a single node
    if (isNodeBody == EGADS_SUCCESS) {

        if (numNode != 1) {
            printf("NodeBody found, but more than one node being reported!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        numElement = numNode;
        feaMesh->numElement = numElement;

        feaMesh->element = (meshElementStruct *) EG_alloc(feaMesh->numElement*sizeof(meshElementStruct));
        if (feaMesh->element == NULL) {
            status = EGADS_MALLOC;
            feaMesh->numElement = 0;
            goto cleanup;
        }

        i = 0;
        (void) initiate_meshElementStruct(&feaMesh->element[i], feaMesh->analysisType);

        status = retrieve_CAPSGroupAttr(enodes[i], &attrName);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "No capsGroup attribute found for node - %d!!", i+1);
            print_AllAttr(aimInfo, enodes[i]);
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(attrMap, attrName, &attrIndex);
        if (status != CAPS_SUCCESS) {
            printf("\tError: capsGroup name %s not found in attribute to index map\n", attrName);
            goto cleanup;
        }

        feaMesh->element[i].elementType = Node;

        feaMesh->element[i].elementID = i+1;

        status = mesh_allocMeshElementConnectivity(&feaMesh->element[i]);
        if (status != CAPS_SUCCESS) goto cleanup;

        feaMesh->element[i].markerID = attrIndex;

        feaMesh->element[i].connectivity[0] = i+1;

        feaData = (feaMeshDataStruct *) feaMesh->element[i].analysisData;

        feaData->propertyID = attrIndex;

        feaData->attrIndex = attrIndex;

        status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
        if (status == CAPS_SUCCESS) {
            feaData->coordID = coordSystemIndex;
        } else {
            feaData->coordID = coordID;
        }

        feaMesh->meshQuickRef.numNode += 1; // Add count

        feaMesh->meshQuickRef.startIndexNode = 0;
        feaMesh->meshQuickRef.useStartIndex = (int) true;

        status = CAPS_SUCCESS;
        goto cleanup;

    } //NODEBODY IF

    /* Determine of the body type */
    status = EG_getTopology(ebody, &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;

    // Can only have "free" edges in wire bodies - Don't want to count the edges of the faces
    //   as "free" edges
    if (bodySubType == WIREBODY) {

        for (i = 0; i < numEdge; i++) {

            status = EG_getInfo(eedges[i], &oclass, &mtype, &topObj, &prev, &next);
            if (status != CAPS_SUCCESS) goto cleanup;
            if (mtype == DEGENERATE) continue;

            status = retrieve_CAPSIgnoreAttr(eedges[i], &attrName);
            if (status == CAPS_SUCCESS) {
                printf("\tcapsIgnore attribute found for edge - %d!!\n", i+1);
                continue;
            }

            status = retrieve_CAPSGroupAttr(eedges[i], &attrName);
            if (status != CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "No capsGroup attribute found for edge - %d!!", i+1);
                print_AllAttr(aimInfo, eedges[i] );
                goto cleanup;
            }

            status = get_mapAttrToIndexIndex(attrMap, attrName, &attrIndex);
            if (status != CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "capsGroup name %s not found in attribute to index map\n", attrName);
                goto cleanup;
            }

            status = EG_getTessEdge(feaMesh->egadsTess, i+1, &numPoint, &xyz, &uv);
            AIM_STATUS(aimInfo, status);

            feaMesh->numElement += numPoint-1;

            AIM_REALL(feaMesh->element, feaMesh->numElement, meshElementStruct, aimInfo, status);
            for (i = numElement; i < feaMesh->numElement; i++) {
                (void) initiate_meshElementStruct(&feaMesh->element[i], feaMesh->analysisType);
            }

            for (j = 0; numElement < feaMesh->numElement; numElement++) {

              feaMesh->element[numElement].elementType = Line;

              feaMesh->element[numElement].elementID = numElement;

              status = mesh_allocMeshElementConnectivity(&feaMesh->element[numElement]);
              AIM_STATUS(aimInfo, status);

              feaMesh->element[numElement].markerID = attrIndex;

              status = EG_localToGlobal(feaMesh->egadsTess, -(i+1), j+1, &gID);
              AIM_STATUS(aimInfo, status);

              feaMesh->element[numElement-1].connectivity[0] = gID;

              status = EG_localToGlobal(feaMesh->egadsTess, -(i+1), j+2, &gID);
              AIM_STATUS(aimInfo, status);

              feaMesh->element[numElement-1].connectivity[1] = gID;

              feaData = (feaMeshDataStruct *) feaMesh->element[numElement].analysisData;

              feaData->propertyID = attrIndex;

              feaData->attrIndex = attrIndex;

              status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
              if (status == CAPS_SUCCESS) {
                  feaData->coordID = coordSystemIndex;
              } else {
                  feaData->coordID = coordID;
              }

              feaMesh->meshQuickRef.numLine += 1; // Add count
            }
        }

        feaMesh->meshQuickRef.startIndexLine = 0;
        feaMesh->meshQuickRef.useStartIndex = (int) true;

        status = CAPS_SUCCESS;
        goto cleanup;

    } //WIREBODY IF

    if (quadMesh == (int) true && numFace > 0) {
        printf("\tGetting quads for BEM!\n");

        // Turn off meshQuick guide if you are getting quads
        feaMesh->meshQuickRef.useStartIndex = (int) false;
    } else {
        feaMesh->meshQuickRef.useStartIndex = (int) true;
        feaMesh->meshQuickRef.startIndexTriangle = numElement;
    }

    // Get Tris and Quads from faces
    for (face = 0; face < numFace; face++) {

        status = retrieve_CAPSIgnoreAttr(efaces[face], &attrName);
        if (status == CAPS_SUCCESS) {
            printf("\tcapsIgnore attribute found for face - %d!!\n", face+1);
            ignoreFound = (int) true;
            continue;
        }

        status = retrieve_CAPSGroupAttr(efaces[face], &attrName);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "No capsGroup attribute found for face - %d!!", face+1);
            print_AllAttr(aimInfo, efaces[face]);
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(attrMap, attrName, &attrIndex);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "capsGroup name %s not found in attribute to index map", attrName);
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
        if (status != CAPS_SUCCESS) coordSystemIndex = 0;


        loadIndex = CAPSMAGIC;
        status = retrieve_CAPSLoadAttr(efaces[face], &attrName);
        if (status == CAPS_SUCCESS) {

            status = get_mapAttrToIndexIndex(loadMap, attrName, &loadIndex);

            if (status != CAPS_SUCCESS) {
                printf("Error: capsLoad name %s not found in attribute to index map\n", attrName);
                goto cleanup;
            }
        }

        if (quadMesh == (int) true) {
            status = EG_getQuads(feaMesh->egadsTess, face+1, &numPoint, &xyz, &uv, &ptype, &pindex, &numPatch);
            AIM_STATUS(aimInfo, status);

        } else numPatch = -1;

        if (numPatch > 0) {

            if (numPatch != 1) {
                status = CAPS_NOTIMPLEMENT;
                printf("feaUtils: EG_localToGlobal accidentally only works for a single quad patch! FIXME!\n");
                goto cleanup;
            }

            qints[face] = 0;
            for (patch = 1; patch <= numPatch; patch++) {

                status = EG_getPatch(feaMesh->egadsTess, face+1, patch, &n1, &n2, &pvindex, &pbounds);
                AIM_STATUS(aimInfo, status);

                for (j = 1; j < n2; j++) {
                    for (i = 1; i < n1; i++) {
                        numElement += 1;

                        feaMesh->meshQuickRef.numQuadrilateral += 1;
                        feaMesh->numElement = numElement;

                        AIM_REALL(feaMesh->element, feaMesh->numElement, meshElementStruct, aimInfo, status);

                        status = initiate_meshElementStruct(&feaMesh->element[numElement-1], feaMesh->analysisType);
                        AIM_STATUS(aimInfo, status);

                        qints[face] += 1;

                        feaMesh->element[numElement-1].elementType = Quadrilateral;

                        feaMesh->element[numElement-1].elementID = numElement;

                        status = mesh_allocMeshElementConnectivity(&feaMesh->element[numElement-1]);
                        AIM_STATUS(aimInfo, status);

                        status = EG_localToGlobal(feaMesh->egadsTess, face+1, pvindex[(i-1)+n1*(j-1)], &gID);
                        AIM_STATUS(aimInfo, status);

                        feaMesh->element[numElement-1].connectivity[0] = gID;

                        status = EG_localToGlobal(feaMesh->egadsTess, face+1, pvindex[(i  )+n1*(j-1)], &gID);
                        AIM_STATUS(aimInfo, status);

                        feaMesh->element[numElement-1].connectivity[1] = gID;

                        status = EG_localToGlobal(feaMesh->egadsTess, face+1, pvindex[(i  )+n1*(j  )], &gID);
                        AIM_STATUS(aimInfo, status);

                        feaMesh->element[numElement-1].connectivity[2] = gID;

                        status = EG_localToGlobal(feaMesh->egadsTess, face+1, pvindex[(i-1)+n1*(j  )], &gID);
                        AIM_STATUS(aimInfo, status);

                        feaMesh->element[numElement-1].connectivity[3] = gID;

                        feaMesh->element[numElement-1].markerID = attrIndex;

                        feaData = (feaMeshDataStruct *) feaMesh->element[numElement-1].analysisData;

                        feaData->propertyID = attrIndex;
                        feaData->attrIndex = attrIndex;
                        feaData->coordID = coordSystemIndex;
                        feaData->loadIndex = loadIndex;

                    }
                }

            }
        } else {
            status = EG_getTessFace(feaMesh->egadsTess, face+1,
                                    &numPoint, &xyz, &uv, &ptype, &pindex,
                                    &numTri, &triConn, &triNeighbor);
            if (status < EGADS_SUCCESS) goto cleanup;

            for (i= 0; i < numTri; i++) {

                numElement += 1;

                feaMesh->meshQuickRef.numTriangle += 1;
                feaMesh->numElement = numElement;


                tempElement = (meshElementStruct *) EG_reall(feaMesh->element,
                                                             feaMesh->numElement*sizeof(meshElementStruct));

                if (tempElement == NULL) {
                    status = EGADS_MALLOC;
                    feaMesh->numElement -= 1;
                    goto cleanup;
                }

                feaMesh->element = tempElement;

                status = initiate_meshElementStruct(&feaMesh->element[numElement-1], feaMesh->analysisType);
                if (status != CAPS_SUCCESS) goto cleanup;

                feaMesh->element[numElement-1].elementType = Triangle;

                feaMesh->element[numElement-1].elementID = numElement;

                status = mesh_allocMeshElementConnectivity(&feaMesh->element[numElement-1]);
                AIM_STATUS(aimInfo, status);

                status = EG_localToGlobal(feaMesh->egadsTess, face+1, triConn[3*i + 0], &gID);
                AIM_STATUS(aimInfo, status);

                feaMesh->element[numElement-1].connectivity[0] = gID;

                status = EG_localToGlobal(feaMesh->egadsTess, face+1, triConn[3*i + 1], &gID);
                AIM_STATUS(aimInfo, status);

                feaMesh->element[numElement-1].connectivity[1] = gID;

                status = EG_localToGlobal(feaMesh->egadsTess, face+1, triConn[3*i + 2], &gID);
                AIM_STATUS(aimInfo, status);

                feaMesh->element[numElement-1].connectivity[2] = gID;

                feaMesh->element[numElement-1].markerID = attrIndex;

                feaData = (feaMeshDataStruct *) feaMesh->element[numElement-1].analysisData;

                feaData->propertyID = attrIndex;
                feaData->attrIndex = attrIndex;
                feaData->coordID = coordSystemIndex;
                feaData->loadIndex = loadIndex;

            }
        }

    }

    if (numPatch > 0) {
        status = EG_attributeAdd(feaMesh->egadsTess, ".mixed", ATTRINT, numFace, qints, NULL, NULL);
        AIM_STATUS(aimInfo, status);
    }

    if (ignoreFound == (int) true) {

        // Look at the nodeID for each node and check to see if it is being used in the element connectivity; if not it is removed
        // Note: that the nodeIDs for the nodes and element connectivity isn't changed, as such if using element connectivity to blindly
        // access a given node this could lead to seg-faults!. mesh_nodeID2Array must be used to access the node array index.
        status = mesh_removeUnusedNodes(feaMesh);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("\tPremature exit in fea_bodyToBem, status = %d\n", status);

    AIM_FREE(iwest );
    AIM_FREE(inorth);
    AIM_FREE(ieast );
    AIM_FREE(isouth);
    AIM_FREE(rpos  );
    AIM_FREE(points );

    AIM_FREE(enodes);
    AIM_FREE(eedges);
    AIM_FREE(efaces);

    AIM_FREE(qints);

    enodes = NULL;
    eedges =NULL;
    efaces =NULL;

    return status;
}

// Set the fea analysis meta data in a mesh
int fea_setAnalysisData( void *aimInfo,                       // (in)  AIM structure
                         mapAttrToIndexStruct *attrMap,       // (in)  map from CAPSGroup names to indexes
                         mapAttrToIndexStruct *coordSystemMap,// (in)  map from CoordSystem names to indexes
                         mapAttrToIndexStruct *constraintMap, // (in)  map from CAPSConstraint names to indexes
                         mapAttrToIndexStruct *loadMap,       // (in)  map from CAPSLoad names to indexes
                         mapAttrToIndexStruct *transferMap,   // (in)  map from CAPSTransfer names to indexes
                         mapAttrToIndexStruct *connectMap,    // (in)  map from CAPSConnect names to indexes
                         mapAttrToIndexStruct *responseMap,   // (in)  map from CAPSResponse names to indexes
                         mapAttrToIndexStruct *referenceMap,  // (in)  map from CAPSReference names to indexes
                         meshStruct *feaMesh)                 // (in/out) FEA mesh structure
{
    int status = 0; // Function return status

    int i, face, edge, node, body; // Indexing
    int dummy;

    int nodeOffset=0, elementOffset=0;

    // Body entities
    int numNode = 0, numEdge = 0, numFace = 0;
    ego *enodes=NULL, *eedges=NULL, *efaces=NULL;
    ego ebody;

    // Meshing
    int elem;

    int     oclass, nchild, *senses;

    // Edge point distributions
    int    *points=NULL, *isouth=NULL, *ieast=NULL, *inorth=NULL, *iwest=NULL;
    double  data[4];
    double *rpos=NULL;
    ego     eref, *echilds;//, topObj, prev, next;

    int isNodeBody, bodySubType = 0; // Body classification

    int pointType, pointTopoIndex;
    double xyzPoint[3];

    // Meshing
    int *n2a=NULL;
    int numPoint = 0, numTri = 0;
    const int  *triConn = NULL, *triNeighbor = NULL; // Triangle

    int gID; // Global id

    const double *xyz, *uv;
    double result[18];

    const int *ptype = NULL, *pindex = NULL;

    // Attributues
    const char *attrName;
    int         attrIndex, coordSystemIndex, loadIndex;

    feaMeshDataStruct *feaData = NULL;
    meshGeomDataStruct *geomData;

    // ---------------------------------------------------------------

    if (feaMesh->meshType == SurfaceMesh || feaMesh->meshType == Surface2DMesh ) {
        printf("Setting FEA Data\n");

        // Get body from tessellation
        status = EG_statusTessBody(feaMesh->egadsTess, &ebody, &dummy, &dummy);
        AIM_STATUS(aimInfo, status);

        // Get number of Nodes, Edges, and Faces in ebody
        status = EG_getBodyTopos(ebody, NULL, NODE, &numNode, &enodes);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(enodes, aimInfo, status);

        status = EG_getBodyTopos(ebody, NULL, EDGE, &numEdge, &eedges);
        AIM_STATUS(aimInfo, status);

        status = EG_getBodyTopos(ebody, NULL, FACE, &numFace, &efaces);
        AIM_STATUS(aimInfo, status);

        // What type of BODY do we have?
        isNodeBody = aim_isNodeBody(ebody, xyzPoint);
        if (isNodeBody < EGADS_SUCCESS) goto cleanup;
        if (isNodeBody == EGADS_SUCCESS) {
          // all attributes are on the body rather than the node for a node body
          enodes[0] = ebody;
        }

        // Fill up the Attributes for the nodes
        for (i = 0; i < feaMesh->numNode; i++) {

            status = EG_getGlobal(feaMesh->egadsTess, feaMesh->node[i].nodeID, &pointType, &pointTopoIndex, xyzPoint);
            AIM_STATUS(aimInfo, status);

            feaData = (feaMeshDataStruct *) feaMesh->node[i].analysisData;

            status = fea_setFEADataPoint(efaces, eedges, enodes,
                                         attrMap,
                                         coordSystemMap,
                                         constraintMap,
                                         loadMap,
                                         transferMap,
                                         connectMap,
                                         responseMap,
                                         referenceMap,
                                         pointType, pointTopoIndex,
                                         feaData); // Set the feaData structure
            AIM_STATUS(aimInfo, status);
        }

        // Fill element information

        // If body is just a single node
        if (isNodeBody == EGADS_SUCCESS) {

            if (feaMesh->numNode != 1) {
                AIM_ERROR(aimInfo, "NodeBody found, but more than one node being reported!\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            i = 0;
            status = retrieve_CAPSGroupAttr(enodes[i], &attrName);
            if (status != CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "No capsGroup attribute found for node - %d!!", i+1);
                print_AllAttr(aimInfo, enodes[i] );
                goto cleanup;
            }

            status = get_mapAttrToIndexIndex(attrMap, attrName, &attrIndex);
            if (status != CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "capsGroup name %s not found in attribute to index map", attrName);
                goto cleanup;
            }

            feaData = (feaMeshDataStruct *) feaMesh->element[i].analysisData;

            feaData->propertyID = attrIndex;

            feaData->attrIndex = attrIndex;

            status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
            if (status != CAPS_SUCCESS) coordSystemIndex = 0;

            feaData->coordID = coordSystemIndex;

            status = CAPS_SUCCESS;
            goto cleanup;

        } //NODEBODY IF

#if 0
        /* Determine the body type */
        status = EG_getTopology(ebody, &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses);
        if (status != EGADS_SUCCESS) goto cleanup;

        // Can only have "free" edges in wire bodies - Don't want to count the edges of the faces
        //   as "free" edges
        if (bodySubType == WIREBODY) {

            for (i = 0; i < feaMesh->meshQuickRef.numLine; i++) {

                status = EG_getInfo(eedges[i], &oclass, &mtype, &topObj, &prev, &next);
                AIM_STATUS(aimInfo, status);
                if (mtype == DEGENERATE) continue;

                if (feaMesh->meshQuickRef.useStartIndex == (int)true)
                  elem = i + feaMesh->meshQuickRef.startIndexLine;
                else if (feaMesh->meshQuickRef.useListIndex == (int)true)
                  elem = feaMesh->meshQuickRef.listIndexLine[i];
                else {
                  status = CAPS_BADOBJECT;
                  printf("DEVELOPER ERROR: Both useStartIndex and useListIndex are true!\n");
                  goto cleanup;
                }

                attrIndex = feaMesh->element[elem].markerID;

                // get the capsGroup attribute string value
                status = get_mapAttrToIndexKeyword(attrMap, attrIndex, &attrName);
                if (status != CAPS_SUCCESS) {
                    printf("\tError: capsGroup index '%d' not found in attribute to index map\n", attrIndex);
                    goto cleanup;
                }

                feaData = (feaMeshDataStruct *) feaMesh->element[elem].analysisData;

                status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
                if (status != CAPS_SUCCESS) coordSystemIndex = 0;

                feaData->propertyID = attrIndex;
                feaData->coordID    = coordSystemIndex;
            }

            status = CAPS_SUCCESS;
            goto cleanup;

        }
#endif

        status = mesh_nodeID2Array(feaMesh, &n2a);
        AIM_STATUS(aimInfo, status);

        // Set line, tri and quad analysis data
        for (elem = 0; elem < feaMesh->numElement; elem++) {

            if (feaMesh->element[elem].elementType != Node &&
                feaMesh->element[elem].elementType != Line &&
                feaMesh->element[elem].elementType != Triangle &&
                feaMesh->element[elem].elementType != Triangle_6 &&
                feaMesh->element[elem].elementType != Quadrilateral &&
                feaMesh->element[elem].elementType != Quadrilateral_8) continue;

            attrIndex = feaMesh->element[elem].markerID;

            // get the capsGroup attribute string value
            status = get_mapAttrToIndexKeyword(attrMap, attrIndex, &attrName);
            if (status != CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "capsGroup index '%d' not found in attribute to index map", attrIndex);
                goto cleanup;
            }

            status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
            if (status != CAPS_SUCCESS) coordSystemIndex = 0;

            if (feaMesh->element[elem].elementType == Node) {
                // get the node index
                node = feaMesh->element[elem].topoIndex;
                if (node < 1 || node > numNode) {
                    AIM_ERROR(aimInfo, "Element '%d': Invalid node topological index: %d, [1-%d]", elem, node, numNode);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                eref = enodes[node-1];

            } else if (feaMesh->element[elem].elementType == Line) {
                // get the edge index
                edge = feaMesh->element[elem].topoIndex;
                if (edge < 1 || edge > numEdge) {
                    printf("Error: Element '%d': Invalid edge topological index: %d, [1-%d]\n", elem, edge, numEdge);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                eref = eedges[edge-1];

            } else {

                // get the face index
                face = feaMesh->element[elem].topoIndex;
                if (face < 1 || face > numFace) {
                    printf("Error: Element '%d': Invalid face topological index: %d, [1-%d]\n", elem, face, numFace);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                eref = efaces[face-1];

                // If this is a face in a capsBound, set geomData in nodes
                status = retrieve_CAPSBoundAttr(eref, &attrName);
                if (status == CAPS_SUCCESS) {

                    status = EG_getTessFace(feaMesh->egadsTess, face,
                                            &numPoint, &xyz, &uv, &ptype, &pindex,
                                            &numTri, &triConn, &triNeighbor);
                    if (status < EGADS_SUCCESS) goto cleanup;

                    for (i = 0; i < numPoint; i++) {

                        status = EG_localToGlobal(feaMesh->egadsTess, face, i+1, &gID);
                        if (status != EGADS_SUCCESS) goto cleanup;

                        geomData = feaMesh->node[n2a[gID]].geomData;

                        if (geomData == NULL) {

                            // Get geometry data for node
                            geomData = EG_alloc(sizeof(meshGeomDataStruct));

                            if (geomData == NULL) {
                                status = EGADS_MALLOC;
                                goto cleanup;
                            }

                            status = initiate_meshGeomDataStruct(geomData);
                            if (status != CAPS_SUCCESS) goto cleanup;

                            geomData->type = ptype[i];
                            geomData->topoIndex = pindex[i];

                            // Want the face index to be set for topoIndex
                            if (geomData->topoIndex < 0) geomData->topoIndex = face+1;

                            geomData->uv[0] = uv[2*i + 0];
                            geomData->uv[1] = uv[2*i + 1];

                            status = EG_evaluate(eref, geomData->uv, result);
                            if (status != EGADS_SUCCESS) goto cleanup;

                            // dU
                            geomData->firstDerivative[0] = result[3];
                            geomData->firstDerivative[1] = result[4];
                            geomData->firstDerivative[2] = result[5];

                            // dV
                            geomData->firstDerivative[3] = result[6];
                            geomData->firstDerivative[4] = result[7];
                            geomData->firstDerivative[5] = result[8];

                            feaMesh->node[n2a[gID]].geomData = geomData;
                        }
                    }
                }
            }

            loadIndex = CAPSMAGIC;
            status = retrieve_CAPSLoadAttr(eref, &attrName);
            if (status == CAPS_SUCCESS) {

                status = get_mapAttrToIndexIndex(loadMap, attrName, &loadIndex);

                if (status != CAPS_SUCCESS) {
                    printf("Error: capsLoad name %s not found in attribute to index map\n", attrName);
                    goto cleanup;
                }
            }

            feaData = (feaMeshDataStruct *) feaMesh->element[elem].analysisData;

            feaData->propertyID = attrIndex;
            feaData->attrIndex = attrIndex;
            feaData->coordID    = coordSystemIndex;
            feaData->loadIndex  = loadIndex;
        }

    } else if (feaMesh->meshType == VolumeMesh) {

        // Show warning statement
        printf("Warning - surface nodes are assumed to be packed sequentially in the volume\n");

        // Loop through reference meshes
        for (body = 0; body < feaMesh->numReferenceMesh; body++) {
            printf("Setting FEA Data from reference mesh %d (of %d)\n", body+1, feaMesh->numReferenceMesh);

            // Get body from tessellation
            status = EG_statusTessBody(feaMesh->referenceMesh[body].egadsTess, &ebody, &dummy, &dummy);
            AIM_STATUS(aimInfo, status);

            // Get number of Nodes, Edges, and Faces in ebody
            status = EG_getBodyTopos(ebody, NULL, NODE, &numNode, &enodes);
            if (status < EGADS_SUCCESS) goto cleanup;
            AIM_NOTNULL(enodes, aimInfo, status);

            status = EG_getBodyTopos(ebody, NULL, EDGE, &numEdge, &eedges);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_getBodyTopos(ebody, NULL, FACE, &numFace, &efaces);
            if (status < EGADS_SUCCESS) goto cleanup;

            // What type of BODY do we have?
            isNodeBody = aim_isNodeBody(ebody, xyzPoint);
            if (isNodeBody < EGADS_SUCCESS) goto cleanup;

            // If body is just a single node
              if (isNodeBody == EGADS_SUCCESS) {

                  printf("NodeBody found, not currently supported for VolumeMesh!\n");
                  status = CAPS_BADVALUE;
                  goto cleanup;

              } //NODEBODY IF

            // Fill up the Attributes for the nodes
            for (i = 0; i < feaMesh->referenceMesh[body].numNode; i++) {

                status = EG_getGlobal(feaMesh->referenceMesh[body].egadsTess, feaMesh->referenceMesh[body].node[i].nodeID, &pointType, &pointTopoIndex, xyzPoint);
                if (status != EGADS_SUCCESS) goto cleanup;

                feaData = (feaMeshDataStruct *) feaMesh->node[i+nodeOffset].analysisData;

                status = fea_setFEADataPoint(efaces, eedges, enodes,
                                             attrMap,
                                             coordSystemMap,
                                             constraintMap,
                                             loadMap,
                                             transferMap,
                                             connectMap,
                                             responseMap,
                                             referenceMap,
                                             pointType, pointTopoIndex,
                                             feaData); // Set the feaData structure
                AIM_STATUS(aimInfo, status);
            }

            // Fill element information //

            // Determine the body type
            status = EG_getTopology(ebody, &eref, &oclass, &bodySubType, data, &nchild, &echilds, &senses);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Can only have "free" edges in wire bodies - Don't want to count the edges of the faces
            //   as "free" edges
            if (bodySubType == WIREBODY) {

                printf("WireBody found, not currently supported for VolumeMesh!\n");
                status = CAPS_BADVALUE;
                goto cleanup;

            } //WIREBODY IF

            // Set tri and quad analysis data
            for (elem = 0; elem < feaMesh->referenceMesh[body].numElement; elem++) {

                if (feaMesh->referenceMesh[body].element[elem].elementType != Triangle &&
                    feaMesh->referenceMesh[body].element[elem].elementType != Triangle_6 &&
                    feaMesh->referenceMesh[body].element[elem].elementType != Quadrilateral &&
                    feaMesh->referenceMesh[body].element[elem].elementType != Quadrilateral_8) continue;

//
//                attrIndex = feaMesh->referenceMesh[body].element[elem].markerID;
//
//                // get the capsGroup attribute string value
//                status = get_mapAttrToIndexKeyword(attrMap, attrIndex, &attrName);
//                if (status != CAPS_SUCCESS) {
//                    printf("\tError: capsGroup index '%d' not found in attribute to index map\n", attrIndex);
//                    goto cleanup;
//                }
//
//                status = get_mapAttrToIndexIndex(coordSystemMap, attrName, &coordSystemIndex);
//                if (status != CAPS_SUCCESS) coordSystemIndex = 0;

                // get the face index
                face = feaMesh->referenceMesh[body].element[elem].topoIndex;
                if (face < 1 || face > numFace) {
                    printf("Error: Element '%d': Invalid face topological index: %d, [1-%d]\n", elem, face, numFace);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }

                loadIndex = CAPSMAGIC;
                status = retrieve_CAPSLoadAttr(efaces[face-1], &attrName);
                if (status == CAPS_SUCCESS) {

                    status = get_mapAttrToIndexIndex(loadMap, attrName, &loadIndex);

                    if (status != CAPS_SUCCESS) {
                        printf("Error: capsLoad name %s not found in attribute to index map\n", attrName);
                        goto cleanup;
                    }
                }

                feaData = (feaMeshDataStruct *) feaMesh->element[elem+elementOffset].analysisData;

                //feaData->propertyID = attrIndex; DONT THINK WE WANT TO CHANGE THE PROPERTY TYPE
                // feaData->coordID    = coordSystemIndex; DONT THINK WE WANT TO CHANGE THE Coordinate TYPE
                feaData->loadIndex  = loadIndex;
            }

            nodeOffset += feaMesh->referenceMesh[body].numNode;
            elementOffset += feaMesh->referenceMesh[body].numElement;
        }

    } else {
        printf("Unknown meshType!\n");
        status = CAPS_BADTYPE;
        goto cleanup;
    }

    status = CAPS_SUCCESS;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in fea_setAnalysisData, status = %d\n", status);

        AIM_FREE(n2a);
        AIM_FREE(iwest );
        AIM_FREE(inorth);
        AIM_FREE(ieast );
        AIM_FREE(isouth);
        AIM_FREE(rpos  );
        AIM_FREE(points );

        AIM_FREE(enodes);
        AIM_FREE(eedges);
        AIM_FREE(efaces);

        enodes = NULL;
        eedges = NULL;
        efaces = NULL;

        return status;
}

// Set feaData for a given point index and topology index. Ego faces, edges, and nodes must be provided along with attribute maps
int fea_setFEADataPoint(ego *faces, ego *edges, ego *nodes,
                        mapAttrToIndexStruct *attrMap,
           /*@unused@*/ mapAttrToIndexStruct *coordSystemMap,
                        mapAttrToIndexStruct *constraintMap,
                        mapAttrToIndexStruct *loadMap,
                        mapAttrToIndexStruct *transferMap,
                        mapAttrToIndexStruct *connectMap,
                        mapAttrToIndexStruct *responseMap,
                        mapAttrToIndexStruct *referenceMap,
                        int pointType, int pointTopoIndex,
                        feaMeshDataStruct *feaData) { // Set the feaData structure

    int status;

    int coordID = 0; // Default coordinate id for mesh

    // Attributes
    const char *attrName;
    int         constraintIndex=-1, loadIndex=-1, transferIndex=-1, connectIndex=-1;
    int         connectLinkIndex=-1, attrIndex=-1, responseIndex=-1, referenceIndex=-1; //coordSystemIndex

    ego object;

    feaData->coordID = coordID;

    // Get attribute index on entity
    constraintIndex = CAPSMAGIC;
    loadIndex = CAPSMAGIC;
    transferIndex = CAPSMAGIC;
    connectIndex = CAPSMAGIC;
    connectLinkIndex = CAPSMAGIC;
    responseIndex = CAPSMAGIC;
    referenceIndex = CAPSMAGIC;

    if (pointType == 0) { // Node

        object = nodes[pointTopoIndex-1];

    } else if (pointType > 0 ) { // Edge

        object = edges[pointTopoIndex-1];

    } else { // Face

        object = faces[pointTopoIndex-1];

    }

    status = retrieve_CAPSGroupAttr(object, &attrName);
    if (status == CAPS_SUCCESS) {
        status = get_mapAttrToIndexIndex(attrMap, attrName, &attrIndex);
        if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
    }

    status = retrieve_CAPSConstraintAttr(object, &attrName);
    if (status == CAPS_SUCCESS) {
        status = get_mapAttrToIndexIndex(constraintMap, attrName, &constraintIndex);
        if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
    }

    status = retrieve_CAPSLoadAttr(object, &attrName);
    if (status == CAPS_SUCCESS) {
        status = get_mapAttrToIndexIndex(loadMap, attrName, &loadIndex);
        if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
    }

    status = retrieve_CAPSBoundAttr(object, &attrName);
    if (status == CAPS_SUCCESS) {
        status = get_mapAttrToIndexIndex(transferMap, attrName, &transferIndex);
        if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
    }

    if (connectMap != NULL) {
        status = retrieve_CAPSConnectAttr(object, &attrName);
        if (status == CAPS_SUCCESS) {
            status = get_mapAttrToIndexIndex(connectMap, attrName, &connectIndex);
            if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
        }

        status = retrieve_CAPSConnectLinkAttr(object, &attrName);
        if (status == CAPS_SUCCESS) {
            status = get_mapAttrToIndexIndex(connectMap, attrName, &connectLinkIndex);
            if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
        }
    }

    if (responseMap != NULL) {
        status = retrieve_CAPSResponseAttr(object, &attrName);
        if (status == CAPS_SUCCESS) {
            status = get_mapAttrToIndexIndex(responseMap, attrName, &responseIndex);
            if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
        }
    }

    if (referenceMap != NULL) {
        status = retrieve_CAPSReferenceAttr(object, &attrName);
        if (status == CAPS_SUCCESS) {
            status = get_mapAttrToIndexIndex(referenceMap, attrName, &referenceIndex);
            if (status != CAPS_SUCCESS && status != CAPS_NOTFOUND && status != CAPS_NULLVALUE) goto cleanup;
        }
    }

    feaData->attrIndex = attrIndex;

    feaData->constraintIndex = constraintIndex;
    feaData->loadIndex = loadIndex;
    feaData->transferIndex = transferIndex;
    feaData->connectIndex = connectIndex;
    feaData->connectLinkIndex = connectLinkIndex;
    feaData->responseIndex = responseIndex;
    feaData->referenceIndex = referenceIndex;

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:

        if (status != CAPS_SUCCESS) printf("Error: Premature exit in fea_setFEADataPoint, status %d\n", status);

        return status;
}

// Initiate (0 out all values and NULL all pointers) of feaProperty in the feaPropertyStruct structure format
int initiate_feaPropertyStruct(feaPropertyStruct *feaProperty) {
    int i;

    feaProperty->name = NULL;

    feaProperty->propertyType = UnknownProperty;

    feaProperty->propertyID = 0; // ID number of property

    feaProperty->materialID = 0; // ID number of material
    feaProperty->materialName = NULL; // Name of material associated with material ID

    // Rods
    feaProperty->crossSecArea           = 0.0; // Bar cross-sectional area
    feaProperty->torsionalConst         = 0.0; // Torsional constant
    feaProperty->torsionalStressReCoeff = 0.0; // Torsional stress recovery coefficient
    feaProperty->massPerLength          = 0.0; // Mass per unit length

    // Bar - see rod for additional variables
    feaProperty->zAxisInertia = 0.0; // Section moment of inertia about the z axis
    feaProperty->yAxisInertia = 0.0; // Section moment of inertia about the y axis

    feaProperty->yCoords[0] = 0.0; // Element y, z coordinates, in the bar cross-section, of
    feaProperty->yCoords[1] = 0.0; //    of four points at which to recover stresses
    feaProperty->yCoords[2] = 0.0;
    feaProperty->yCoords[3] = 0.0;

    feaProperty->zCoords[0] = 0.0;
    feaProperty->zCoords[1] = 0.0;
    feaProperty->zCoords[2] = 0.0;
    feaProperty->zCoords[3] = 0.0;

    feaProperty->areaShearFactors[0] = 0.0; // Area factors for shear
    feaProperty->areaShearFactors[1] = 0.0;

    feaProperty->crossProductInertia = 0.0; // Section cross-product of inertia

    feaProperty->crossSecType = NULL; // Type of cross sections
    for (i = 0; i < 10; i++) feaProperty->crossSecDimension[i] = 0.0; // Dimensions

    feaProperty->orientationVec[0] = 0.0; // Orientation vector
    feaProperty->orientationVec[1] = 0.0;
    feaProperty->orientationVec[2] = 0.0;

    // Shear

    // Shell
    feaProperty->membraneThickness = 0.0;     // Membrane thickness
    feaProperty->materialBendingID = 0;     /// ID number of material for bending - if not specified and bendingInertiaRatio > 0 this value defaults to materialID
    feaProperty->bendingInertiaRatio = 1.0; // Ratio of actual bending moment inertia (I) to bending inertia of a solid
                                            //   plate of thickness TM Real  - default 1.0
    feaProperty->materialShearID = 0;       // ID number of material for shear - if not specified and shearMembraneRatio > 0 this value defaults to materialID
    feaProperty->shearMembraneRatio = 5.0/6.0;     // Ratio of shear to membrane thickness  - default 5/6
    feaProperty->massPerArea = 0.0;           // Mass per unit area
    //feaProperty->neutralPlaneDist[0] = 0;   // Distances from the neutral plane of the plate to locations where
    //feaProperty->neutralPlaneDist[1] = 0;   //   stress is calculate

    feaProperty->zOffsetRel = 0.0; // Offset from the surface of grid points to the element reference plane

    feaProperty->compositeShearBondAllowable = 0.0;  // Shear stress limit for bonding between plies
    feaProperty->compositeFailureTheory = NULL;      // HILL, HOFF, TSAI, STRN
    feaProperty->compositeSymmetricLaminate = (int) false;  // 1- SYM only half the plies are specified, for odd number plies 1/2 thickness of center ply is specified
                                                            // the first ply is the bottom ply in the stack, default (0) all plies specified
    feaProperty->numPly = 0;
    feaProperty->compositeMaterialID = NULL;      // Vector of material ID's from bottom to top for all plies
    feaProperty->compositeThickness = NULL;        // Vector of thicknesses from bottom to top for all plies
    feaProperty->compositeOrientation = NULL;      // Vector of orientations from bottom to top for all plies

    // Solid

    // Concentrated Mass
    feaProperty->mass = 0.0; // Mass value
    feaProperty->massOffset[0]  = 0.0; // Offset distance from the grid point to the center of gravity
    feaProperty->massOffset[1]  = 0.0;
    feaProperty->massOffset[2]  = 0.0;
    feaProperty->massInertia[0] = 0.0; // Mass moment of inertia measured at the mass center of gravity
    feaProperty->massInertia[1] = 0.0;
    feaProperty->massInertia[2] = 0.0;
    feaProperty->massInertia[3] = 0.0;
    feaProperty->massInertia[4] = 0.0;
    feaProperty->massInertia[5] = 0.0;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaProperty in the feaPropertyStruct structure format
int destroy_feaPropertyStruct(feaPropertyStruct *feaProperty) {

    int i;

    if (feaProperty->name != NULL) EG_free(feaProperty->name);

    feaProperty->propertyType = UnknownProperty;

    feaProperty->propertyID = 0; // ID number of property

    feaProperty->materialID = 0; // ID number of material

    if (feaProperty->materialName != NULL) EG_free(feaProperty->materialName); // Name of material associated with material ID
    feaProperty->materialName = NULL;

    // Rods
    feaProperty->crossSecArea           = 0; // Bar cross-sectional area
    feaProperty->torsionalConst         = 0; // Torsional constant
    feaProperty->torsionalStressReCoeff = 0; // Torsional stress recovery coefficient
    feaProperty->massPerLength          = 0; // Mass per unit length

    // Bar - see rod for additional variables
    feaProperty->zAxisInertia = 0; // Section moment of inertia about the z axis
    feaProperty->yAxisInertia = 0; // Section moment of inertia about the y axis

    feaProperty->yCoords[0] = 0; // Element y, z coordinates, in the bar cross-section, of
    feaProperty->yCoords[1] = 0; //    of four points at which to recover stresses
    feaProperty->yCoords[2] = 0;
    feaProperty->yCoords[3] = 0;

    feaProperty->zCoords[0] = 0;
    feaProperty->zCoords[1] = 0;
    feaProperty->zCoords[2] = 0;
    feaProperty->zCoords[3] = 0;

    feaProperty->areaShearFactors[0] = 0; // Area factors for shear
    feaProperty->areaShearFactors[1] = 0;

    feaProperty->crossProductInertia = 0; // Section cross-product of inertia

    if (feaProperty->crossSecType != NULL) EG_free(feaProperty->crossSecType);
    feaProperty->crossSecType = NULL; // Type of cross sections
    for (i = 0; i < 10; i++) feaProperty->crossSecDimension[i] = 0.0; // Dimensions

    feaProperty->orientationVec[0] = 0.0; // Orientation vector
    feaProperty->orientationVec[1] = 0.0;
    feaProperty->orientationVec[2] = 0.0;

    // Shear

    // Shell
    feaProperty->membraneThickness = 0; // Membrane thickness
    feaProperty->materialBendingID = 0;     // ID number of material for bending
    feaProperty->bendingInertiaRatio = 1.0; // Ratio of actual bending moment inertia (I) to bending inertia of a solid
    //   plate of thickness TM Real  - default 1.0
    feaProperty->materialShearID = 0;       // ID number of material for shear
    feaProperty->shearMembraneRatio  = 5.0/6.0;  // Ratio of shear to membrane thickness  - default 5/6
    feaProperty->massPerArea = 0; // Mass per unit area
    //feaProperty->neutralPlaneDist[0] = 0;   // Distances from the neutral plane of the plate to locations where
    //feaProperty->neutralPlaneDist[1] = 0;   //   stress is calculate

    feaProperty->zOffsetRel = 0.0; // Offset from the surface of grid points to the element reference plane

    feaProperty->numPly = 0;
    feaProperty->compositeShearBondAllowable = 0.0;  // Shear stress limit for bonding between plies
    if (feaProperty->compositeFailureTheory != NULL) EG_free(feaProperty->compositeFailureTheory); // HILL, HOFF, TSAI, STRN
    feaProperty->compositeFailureTheory = NULL;
    feaProperty->compositeSymmetricLaminate = (int) false;  // 1- SYM only half the plies are specified, for odd number plies 1/2 thickness of center ply is specified
                                                            // the first ply is the bottom ply in the stack, default (0) all plies specified
    if (feaProperty->compositeMaterialID != NULL) EG_free(feaProperty->compositeMaterialID);
    if (feaProperty->compositeThickness != NULL) EG_free(feaProperty->compositeThickness);
    if (feaProperty->compositeOrientation != NULL) EG_free(feaProperty->compositeOrientation);

    feaProperty->compositeMaterialID = NULL;
    feaProperty->compositeThickness = NULL;
    feaProperty->compositeOrientation = NULL;

    // Solid

    // Concentrated Mass
    feaProperty->mass = 0.0; // Mass value
    feaProperty->massOffset[0]  = 0.0; // Offset distance from the grid point to the center of gravity
    feaProperty->massOffset[1]  = 0.0;
    feaProperty->massOffset[2]  = 0.0;
    feaProperty->massInertia[0] = 0.0; // Mass moment of inertia measured at the mass center of gravity
    feaProperty->massInertia[1] = 0.0;
    feaProperty->massInertia[2] = 0.0;
    feaProperty->massInertia[3] = 0.0;
    feaProperty->massInertia[4] = 0.0;
    feaProperty->massInertia[5] = 0.0;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaMaterial in the feaMaterialStruct structure format
int initiate_feaMaterialStruct(feaMaterialStruct *feaMaterial) {

    feaMaterial->name = NULL;  // Material name

    feaMaterial->materialType = UnknownMaterial; // Set

    feaMaterial->materialID = 0; // ID number of material

    feaMaterial->youngModulus = 0.0; // E - Young's Modulus
    feaMaterial->shearModulus = 0.0; // G - Shear Modulus
    feaMaterial->poissonRatio = 0.0; // Poisson's Ratio
    feaMaterial->density      = 0.0;    // Rho - material mass density
    feaMaterial->thermalExpCoeff = 0.0; //Coefficient of thermal expansion
    feaMaterial->temperatureRef  = 0.0; // Reference temperature
    feaMaterial->dampingCoeff    = 0.0; // Damping coefficient
    feaMaterial->tensionAllow    = 0.0; // Tension allowable for the material
    feaMaterial->compressAllow   = 0.0; // Compression allowable for the material
    feaMaterial->shearAllow      = 0.0; // Shear allowable for the material

    feaMaterial->youngModulusLateral = 0.0; // Young's Modulus in the lateral direction
    feaMaterial->shearModulusTrans1Z = 0.0; // Transverse shear modulus in the 1-Z plane
    feaMaterial->shearModulusTrans2Z = 0.0; // Transverse shear modulus in the 2-Z plane
    feaMaterial->tensionAllowLateral    = 0.0; // Tension allowable for the material
    feaMaterial->compressAllowLateral   = 0.0; // Compression allowable for the material
    feaMaterial->thermalExpCoeffLateral = 0.0; //Coefficient of thermal expansion
    feaMaterial->allowType = 0;

    feaMaterial->gmat[0] = 0.0;
    feaMaterial->gmat[1] = 0.0;
    feaMaterial->gmat[2] = 0.0;
    feaMaterial->gmat[3] = 0.0;
    feaMaterial->gmat[4] = 0.0;
    feaMaterial->gmat[5] = 0.0;

    feaMaterial->specificHeat = 0;
    feaMaterial->kappa = 0;

    feaMaterial->K[0] = 0;
    feaMaterial->K[1] = 0;
    feaMaterial->K[2] = 0;
    feaMaterial->K[3] = 0;
    feaMaterial->K[4] = 0;
    feaMaterial->K[5] = 0;

    feaMaterial->honeycombCellSize = -1;     // (CS) Honeycomb sandwich core cell size. Required if material defines the core of a honeycomb sandwich and dimpling stability index is desired (LAM = HCS on the PCOMP entry).
    feaMaterial->honeycombYoungModulus = -1; // (EC) Honeycomb sandwich core Young's modulus used for stability index analysis.
    feaMaterial->honeycombShearModulus = -1; // (GC) Honeycomb sandwich core shear modulus used for stability index analysis.

    feaMaterial->fractureAngle = -1; // (ALPHA0) Fracture angle for uniaxial transverse compression in degrees. Used in the NASA LaRC02 failure theory only (see LARC02 in PCOMP entry).

    feaMaterial->interlaminarShearAllow = -1; // (SB) Allowable inter-laminar shear stress of the composite laminate bonding material

    feaMaterial->fiberYoungModulus = -1; // (EF1) Modulus of elasticity of fiber
    feaMaterial->fiberPoissonRatio = -1; // (NUF12) Poisson's ratio of fiber
    feaMaterial->meanStressFactor = -1;  // (MSMF) Mean stress magnification factor

    feaMaterial->transTensionSlope = -1;     // (PNPT) Failure envelop slope parameter for transverse tension
    feaMaterial->transCompressionSlope = -1; // (PNPC) Failure envelop slope parameter for transverse compression

    feaMaterial->compositeFailureTheory = NULL; // (FT) Composite failure theory

    feaMaterial->interlaminarNormalStressAllow = -1; // (NB) Allowable inter-laminar normal stress of the composite laminate bonding material (allowable interlaminar normal stress).

    feaMaterial->youngModulusThick = -1; // (E3) Modulus of elasticity in thickness direction, also defined as the matrix direction or 3-direction.
    feaMaterial->poissonRatio23 = 0; // (NU23) Poisson's ratio ( for uniaxial loading in 2-direction).
    feaMaterial->poissonRatio31 = 0; // (NU31) Poisson's ratio ( for uniaxial loading in 3-direction).

    feaMaterial->youngModulusFactor = -1;        // (E1RSF) Longitudinal modulus of elasticity reduction scale factor for nonlinear composite Progressive Ply Failure Analysis (PPFA)
    feaMaterial->youngModulusLateralFactor = -1; // (E2RSF) Lateral modulus of elasticity reduction scale factor for nonlinear composite Progressive Ply Failure Analysis (PPFA)

    feaMaterial->shearModulusFactor = -1; // (G12RSF) In-plane shear modulus reduction scale factor for nonlinear composite Progressive Ply Failure Analysis (PPFA)
    feaMaterial->shearModulusTrans1ZFactor = -1; // (G1ZRSF) Transverse shear modulus reduction scale factor in 1-Z plane for nonlinear composite Progressive Ply Failure Analysis (PPFA)
    feaMaterial->shearModulusTrans2ZFactor = -1; // (G1ZRSF) Transverse shear modulus reduction scale factor in 2-Z plane for nonlinear composite Progressive Ply Failure Analysis (PPFA)

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaMaterial in the feaMaterialStruct structure format
int destroy_feaMaterialStruct(feaMaterialStruct *feaMaterial) {

    AIM_FREE(feaMaterial->name);
    AIM_FREE(feaMaterial->compositeFailureTheory);

    initiate_feaMaterialStruct(feaMaterial);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaUnits in the feaUnitsStruct structure format
int initiate_feaUnitsStruct(feaUnitsStruct *feaUnits) {

    feaUnits->densityVol = NULL;
    feaUnits->densityArea = NULL;
    feaUnits->mass = NULL;
    feaUnits->length = NULL;
    feaUnits->pressure = NULL;
    feaUnits->temperature = NULL;
    feaUnits->momentOfInertia = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaUnits in the feaUnitsStruct structure format
int destroy_feaUnitsStruct(feaUnitsStruct *feaUnits) {

    AIM_FREE(feaUnits->densityVol);
    AIM_FREE(feaUnits->densityArea);
    AIM_FREE(feaUnits->mass);
    AIM_FREE(feaUnits->length);
    AIM_FREE(feaUnits->pressure);
    AIM_FREE(feaUnits->temperature);
    AIM_FREE(feaUnits->momentOfInertia);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaConstraint in the feaConstraintStruct structure format
int initiate_feaConstraintStruct(feaConstraintStruct *feaConstraint) {

    feaConstraint->name = NULL;  // Constraint name

    feaConstraint->constraintType = UnknownConstraint; // Constraint type

    feaConstraint->constraintID = 0; // ID number of constraint

    feaConstraint->numGridID = 0; // Component number of grid
    feaConstraint->gridIDSet = NULL; // List of component number of grids to apply constraint to

    feaConstraint->dofConstraint = 0; // Number to indicate DOF constraints
    feaConstraint->gridDisplacement = 0;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaConstraint in the feaConstraintStruct structure format
int destroy_feaConstraintStruct(feaConstraintStruct *feaConstraint) {

    if (feaConstraint->name != NULL) EG_free(feaConstraint->name);
    feaConstraint->name = NULL;  // Constraint name

    feaConstraint->constraintType = UnknownConstraint; // Constraint type

    feaConstraint->constraintID = 0; // ID number of constraint

    feaConstraint->numGridID = 0; // Number of grid IDs in grid ID set

    if (feaConstraint->gridIDSet != NULL) EG_free(feaConstraint->gridIDSet); // List of component number of grids to apply constraint to

    feaConstraint->dofConstraint = 0; // Number to indicate DOF constraints
    feaConstraint->gridDisplacement = 0; // The value for the displacement

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaSupport in the feaSupportStruct structure format
int initiate_feaSupportStruct(feaSupportStruct *feaSupport) {

    feaSupport->name = NULL;  // Support name

    feaSupport->supportID = 0; // ID number of support

    feaSupport->numGridID = 0; // Component number of grid
    feaSupport->gridIDSet = NULL; // List of component number of grids to apply support to

    feaSupport->dofSupport = 0; // Number to indicate DOF supports

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaSupport in the feaSupportStruct structure format
int destroy_feaSupportStruct(feaSupportStruct *feaSupport) {

    if (feaSupport->name != NULL) EG_free(feaSupport->name);
    feaSupport->name = NULL;  // Support name

    feaSupport->supportID = 0; // ID number of support

    feaSupport->numGridID = 0; // Number of grid IDs in grid ID set

    if (feaSupport->gridIDSet != NULL) EG_free(feaSupport->gridIDSet); // List of component number of grids to apply support to

    feaSupport->dofSupport = 0; // Number to indicate DOF supports

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaAnalysis in the feaAnalysisStruct structure format
int initiate_feaAnalysisStruct(feaAnalysisStruct *feaAnalysis) {

    feaAnalysis->name = NULL; // Analysis name

    feaAnalysis->analysisType = UnknownAnalysis; // Type of analysis

    feaAnalysis->analysisID = 0; // ID number of analysis

    // Loads for the analysis
    feaAnalysis->numLoad = 0;     // Number of loads in the analysis
    feaAnalysis->loadSetID = NULL; // List of the load IDSs

    // Constraints for the analysis
    feaAnalysis->numConstraint = 0;   // Number of constraints in the analysis
    feaAnalysis->constraintSetID = NULL; // List of constraint IDs

    // Supports for the analysis
    feaAnalysis->numSupport = 0;   // Number of supports in the analysis
    feaAnalysis->supportSetID = NULL; // List of support IDs

    // Optimization constraints
    feaAnalysis->numDesignConstraint = 0; // Number of design constraints
    feaAnalysis->designConstraintSetID = NULL; // List of design constraint IDs

    // Optimization responses
    feaAnalysis->numDesignResponse = 0; // Number of design responses
    feaAnalysis->designResponseSetID = NULL; // List of design response IDs

    // MASSSET
    feaAnalysis->numMassIncrement = 0; // Number of mass increments
    feaAnalysis->massIncrementSetID = NULL; // List of mass increment IDs

    // Eigenvalue
    feaAnalysis->extractionMethod = NULL;

    feaAnalysis->frequencyRange[0] = 0;
    feaAnalysis->frequencyRange[1] = 0;

    feaAnalysis->numEstEigenvalue = 0;
    feaAnalysis->numDesiredEigenvalue = 0;
    feaAnalysis->eigenNormalization = NULL;

    feaAnalysis->gridNormalization = 0;
    feaAnalysis->componentNormalization = 0;

    feaAnalysis->lanczosMode = 2; //Lanczos mode for calculating eigenvalues
    feaAnalysis->lanczosType = NULL; //Lanczos matrix type (DPB, DGB)

    // Trim
    feaAnalysis->numMachNumber = 0;
    feaAnalysis->machNumber = NULL; // Mach number
    feaAnalysis->dynamicPressure = 0.0; // Dynamic pressure
    //feaAnalysis->numDensity = 0;
    feaAnalysis->density = 0.0; // Density
    feaAnalysis->aeroSymmetryXY = NULL;
    feaAnalysis->aeroSymmetryXZ = NULL;

    feaAnalysis->numRigidVariable = 0; // Number of rigid trim variables
    feaAnalysis->rigidVariable = NULL; // List of character labels identifying rigid trim variables, size=[numRigidVariables]

    feaAnalysis->numRigidConstraint = 0; // Number of rigid trim constrained variables
    feaAnalysis->rigidConstraint = NULL; // List of character labels identifying rigid constrained trim variables, size=[numRigidConstraint]
    feaAnalysis->magRigidConstraint = NULL; // Magnitude of rigid constrained trim variables, size=[numRigidConstraint]

    feaAnalysis->numControlConstraint = 0; // Number of control surface constrained variables
    feaAnalysis->controlConstraint = NULL; // List of character labels identifying control surfaces to be constrained trim variables, size=[numControlConstraint]
    feaAnalysis->magControlConstraint = NULL; // Magnitude of control surface constrained variables, size=[numControlConstraint]

    // Flutter
    feaAnalysis->numReducedFreq = 0;
    feaAnalysis->reducedFreq = NULL;

    feaAnalysis->numFlutterVel = 0;
    feaAnalysis->flutterVel = NULL;

    feaAnalysis->flutterConvergenceParam = 0.0;

    feaAnalysis->visualFlutter = (int) false;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaAnalysis in the feaAnalysisStruct structure format
int destroy_feaAnalysisStruct(feaAnalysisStruct *feaAnalysis) {
    int status; // Function return status

    if (feaAnalysis->name != NULL) EG_free(feaAnalysis->name);
    feaAnalysis->name = NULL; // Analysis name

    feaAnalysis->analysisType = UnknownAnalysis; // Type of analysis

    feaAnalysis->analysisID = 0; // ID number of analysis

    // Loads for the analysis
    feaAnalysis->numLoad = 0;     // Number of loads in the analysis
    if(feaAnalysis->loadSetID != NULL) EG_free(feaAnalysis->loadSetID);
    feaAnalysis->loadSetID = NULL; // List of the load IDSs

    // Constraints for the analysis
    feaAnalysis->numConstraint = 0;   // Number of constraints in the analysis
    if (feaAnalysis->constraintSetID != NULL) EG_free(feaAnalysis->constraintSetID);
    feaAnalysis->constraintSetID = NULL; // List of constraint IDs

    // Supports for the analysis
    feaAnalysis->numSupport = 0;   // Number of supports in the analysis
    if (feaAnalysis->supportSetID != NULL) EG_free(feaAnalysis->supportSetID);
    feaAnalysis->supportSetID = NULL; // List of support IDs

    // Optimization constraints
    feaAnalysis->numDesignConstraint = 0; // Number of design constraints
    if (feaAnalysis->designConstraintSetID != NULL) EG_free(feaAnalysis->designConstraintSetID);
    feaAnalysis->designConstraintSetID = NULL; // List of design constraint IDs

    // Optimization responses
    feaAnalysis->numDesignResponse = 0; // Number of design responses
    if (feaAnalysis->designResponseSetID != NULL) EG_free(feaAnalysis->designResponseSetID);
    feaAnalysis->designResponseSetID = NULL; // List of design response IDs

    // MASSSET
    feaAnalysis->numMassIncrement = 0; // Number of design responses
    AIM_FREE(feaAnalysis->massIncrementSetID);

    // Eigenvalue
    if (feaAnalysis->extractionMethod != NULL) EG_free(feaAnalysis->extractionMethod);
    feaAnalysis->extractionMethod = NULL;

    feaAnalysis->frequencyRange[0] = 0;
    feaAnalysis->frequencyRange[1] = 0;

    feaAnalysis->numEstEigenvalue = 0;
    feaAnalysis->numDesiredEigenvalue = 0;

    AIM_FREE(feaAnalysis->eigenNormalization);

    feaAnalysis->gridNormalization = 0;
    feaAnalysis->componentNormalization = 0;

    feaAnalysis->lanczosMode = 2; //Lanczos mode for calculating eigenvalues

    if (feaAnalysis->lanczosType != NULL) EG_free(feaAnalysis->lanczosType);
    feaAnalysis->lanczosType = NULL; //Lanczos matrix type (DPB, DGB)

    // Trim
    feaAnalysis->numMachNumber = 0;
    AIM_FREE(feaAnalysis->machNumber); // Mach number

    //feaAnalysis->machNumber = 0.0; // Mach number
    feaAnalysis->dynamicPressure = 0.0; // Dynamic pressure
    feaAnalysis->density = 0.0; // Density

    if (feaAnalysis->aeroSymmetryXY != NULL) EG_free(feaAnalysis->aeroSymmetryXY);
    if (feaAnalysis->aeroSymmetryXZ != NULL) EG_free(feaAnalysis->aeroSymmetryXZ);

    if (feaAnalysis->rigidVariable != NULL) {
        status = string_freeArray(feaAnalysis->numRigidVariable, &feaAnalysis->rigidVariable);
        if (status != CAPS_SUCCESS) printf("Status %d during string_freeArray\n", status);
    }

    feaAnalysis->numRigidVariable = 0; // Number of trim rigid trim variables
    feaAnalysis->rigidVariable = NULL; // List of character labels identifying rigid trim variables, size=[numRigidVariables]

    if (feaAnalysis->rigidConstraint != NULL) {
        status = string_freeArray(feaAnalysis->numRigidConstraint, &feaAnalysis->rigidConstraint);
        if (status != CAPS_SUCCESS) printf("Status %d during string_freeArray\n", status);
    }

    feaAnalysis->numRigidConstraint = 0; // Number of rigid trim constrained variables
    feaAnalysis->rigidConstraint = NULL; // List of character labels identifying rigid constrained trim variables, size=[numRigidConstraint]

    if (feaAnalysis->magRigidConstraint != NULL) EG_free(feaAnalysis->magRigidConstraint);
    feaAnalysis->magRigidConstraint = NULL; // Magnitude of rigid constrained trim variables, size=[numRigidConstraint]

    if (feaAnalysis->controlConstraint != NULL) {
        status = string_freeArray(feaAnalysis->numControlConstraint, &feaAnalysis->controlConstraint);
        if (status != CAPS_SUCCESS) printf("Status %d during string_freeArray\n", status);
    }

    feaAnalysis->numControlConstraint = 0; // Number of control surface constrained variables
    feaAnalysis->controlConstraint = NULL; // List of character labels identifying control surfaces to be constrained trim variables, size=[numControlConstraint]

    EG_free(feaAnalysis->magControlConstraint);
    feaAnalysis->magControlConstraint = NULL; // Magnitude of control surface constrained variables, size=[numControlConstraint]

    // Flutter
    feaAnalysis->numReducedFreq = 0;
    EG_free(feaAnalysis->reducedFreq);
    feaAnalysis->reducedFreq = NULL;

    feaAnalysis->numFlutterVel = 0;
    EG_free(feaAnalysis->flutterVel);
    feaAnalysis->flutterVel = NULL;

    feaAnalysis->flutterConvergenceParam = 0.0;

    feaAnalysis->visualFlutter = (int) false;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaLoad in the feaLoadStruct structure format
int initiate_feaLoadStruct(feaLoadStruct *feaLoad) {

    feaLoad->name = NULL;  // Load name

    feaLoad->loadType = UnknownLoad; // Load type

    feaLoad->loadID = 0; // ID number of load

    feaLoad->loadScaleFactor = 1; // Scale factor for when combining loads

    // Concentrated force at a grid point
    feaLoad->numGridID = 0; // Number of grid IDs in grid ID set
    feaLoad->gridIDSet = NULL; // List of grid IDs to apply the constraint to. size = [numGridID]
    feaLoad->coordSystemID= 0; // Component number of coordinate system in which force vector is specified

    feaLoad->forceScaleFactor= 0.0; // Overall scale factor for the force
    feaLoad->directionVector[0] = 0.0;   // [0]-x, [1]-y, [2]-z components of the force vector
    feaLoad->directionVector[1] = 0.0;
    feaLoad->directionVector[2] = 0.0;

    // Concentrated moment at a grid point (also uses coordSystemID and directionVector)
    feaLoad->momentScaleFactor= 0.0; // Overall scale factor for the moment

    // Gravitational load (also uses coordSystemID and directionVector)
    feaLoad->gravityAcceleration= 0.0; // Gravitational acceleration

    // Pressure load
    feaLoad->pressureForce= 0.0; // Pressure value
    feaLoad->pressureDistributeForce[0] = 0.0; // Pressure load at a specified grid location in the element
    feaLoad->pressureDistributeForce[1] = 0.0;
    feaLoad->pressureDistributeForce[2] = 0.0;
    feaLoad->pressureDistributeForce[3] = 0.0;

    feaLoad->pressureMultiDistributeForce = NULL; // Unique pressure load at a specified grid location for
                                                  // each element in elementIDSet size = [numElementID][4] - used in type PressureExternal
                                                  // where the pressure force is being provided by an external source (i.e. data transfer)

    feaLoad->numElementID = 0; // Number of elements IDs in element ID set
    feaLoad->elementIDSet = NULL; // List element IDs in which to apply the load. size = [numElementID]

    // Rotational velocity (also uses coordSystemID and directionVector)
    feaLoad->angularVelScaleFactor = 0.0; // Overall scale factor for the angular velocity
    feaLoad->angularAccScaleFactor = 0.0; // Overall scale factor for the angular acceleration

    // Thermal load - the temperature at a grid point - use gridIDSet
    feaLoad->temperature = 0.0; // Temperature value
    feaLoad->temperatureDefault = 0.0; // Default temperature of grid point explicitly not used

    feaLoad->temperatureMultiDistribute = NULL;

    return CAPS_SUCCESS;
}

// Copy feaLoad in the feaLoadStruct structure format
// assumes that copy has been initialized with initiate_feaLoadStruct
int copy_feaLoadStruct(void *aimInfo, feaLoadStruct *feaLoad, feaLoadStruct *copy) {

    int status = CAPS_SUCCESS;
    int i;

    if (feaLoad->name != NULL) {
      AIM_STRDUP(copy->name, feaLoad->name, aimInfo, status); // Load name
    }

    copy->loadType = feaLoad->loadType; // Load type

    copy->loadID = feaLoad->loadID; // ID number of load

    copy->loadScaleFactor = feaLoad->loadScaleFactor; // Scale factor for when combining loads

    // Concentrated force at a grid point
    copy->numGridID = feaLoad->numGridID; // Number of grid IDs in grid ID set
    // List of grid IDs to apply the constraint to
    AIM_ALLOC(copy->gridIDSet, feaLoad->numGridID, int, aimInfo, status);
    for (i = 0; i < feaLoad->numGridID; i++)
        copy->gridIDSet[i] = feaLoad->gridIDSet[i];

    copy->coordSystemID = feaLoad->coordSystemID; // Component number of coordinate system in which force vector is specified
    copy->forceScaleFactor = feaLoad->forceScaleFactor; // Overall scale factor for the force
    copy->directionVector[0] = feaLoad->directionVector[0];   // [0]-x, [1]-y, [2]-z components of the force vector
    copy->directionVector[1] = feaLoad->directionVector[1];
    copy->directionVector[2] = feaLoad->directionVector[2];

    // Concentrated moment at a grid pofeaLoad->(also uses coordSystemID and directionVector)
    copy->momentScaleFactor = feaLoad->momentScaleFactor; // Overall scale factor for the moment

    // Gravitational load (also uses coordSystemID and directionVector)
    copy->gravityAcceleration = feaLoad->gravityAcceleration; // Gravitational acceleration

    // Pressure load
    copy->pressureForce = feaLoad->pressureForce; // Pressure value

    copy->pressureDistributeForce[0] = feaLoad->pressureDistributeForce[0]; // Pressure load at a specified grid location in the element
    copy->pressureDistributeForce[1] = feaLoad->pressureDistributeForce[1];
    copy->pressureDistributeForce[2] = feaLoad->pressureDistributeForce[2];
    copy->pressureDistributeForce[3] = feaLoad->pressureDistributeForce[3];

    // Unique pressure load at a specified grid location for
    //each element in elementIDSet size = [4*numElementID]- used in type PressureExternal
    if (feaLoad->pressureMultiDistributeForce != NULL) {
        AIM_ALLOC(copy->pressureMultiDistributeForce, feaLoad->numElementID, double, aimInfo, status);
        for (i = 0; i < 4*feaLoad->numElementID; i++)
            copy->pressureMultiDistributeForce[i] = feaLoad->pressureMultiDistributeForce[i];
    } else {
        copy->pressureMultiDistributeForce = NULL;
    }

    copy->numElementID = feaLoad->numElementID;  // Number of elements IDs in element ID set
    AIM_ALLOC(copy->elementIDSet, feaLoad->numElementID, int, aimInfo, status); // List element IDs in which to apply the load
    for (i = 0; i < feaLoad->numElementID; i++)
        copy->elementIDSet[i] = feaLoad->elementIDSet[i];

    // Rotational velocity (also uses coordSystemID and directionVector)
    copy->angularVelScaleFactor = feaLoad->angularVelScaleFactor; // Overall scale factor for the angular velocity
    copy->angularAccScaleFactor = feaLoad->angularAccScaleFactor; // Overall scale factor for the angular acceleration

    // Thermal load - the temperature at a grid point - use gridIDSet
    copy->temperature = feaLoad->temperature; // Temperature value
    copy->temperatureDefault = feaLoad->temperatureDefault; // Default temperature of grid point explicitly not used

    if (feaLoad->temperatureMultiDistribute != NULL) {
        AIM_ALLOC(copy->temperatureMultiDistribute, feaLoad->numGridID, double, aimInfo, status);
        for (i = 0; i < feaLoad->numGridID; i++)
            copy->temperatureMultiDistribute[i] = feaLoad->temperatureMultiDistribute[i];
    } else {
        copy->temperatureMultiDistribute = NULL;
    }

cleanup:
    return status;
}

// Destroy (0 out all values and NULL all pointers) of feaLoad in the feaLoadStruct structure format
int destroy_feaLoadStruct(feaLoadStruct *feaLoad) {

    AIM_FREE(feaLoad->name); // Load name

    feaLoad->loadType = UnknownLoad; // Load type

    feaLoad->loadID = 0; // ID number of load

    feaLoad->loadScaleFactor = 1; // Scale factor for when combining loads

    // Concentrated force at a grid point
    feaLoad->numGridID = 0; // Number of grid IDs in grid ID set
    AIM_FREE(feaLoad->gridIDSet); // List of grid IDs to apply the constraint to

    feaLoad->coordSystemID= 0; // Component number of coordinate system in which force vector is specified
    feaLoad->forceScaleFactor= 0; // Overall scale factor for the force
    feaLoad->directionVector[0] = 0;   // [0]-x, [1]-y, [2]-z components of the force vector
    feaLoad->directionVector[1] = 0;
    feaLoad->directionVector[2] = 0;

    // Concentrated moment at a grid pofeaLoad->(also uses coordSystemID and directionVector)
    feaLoad->momentScaleFactor= 0; // Overall scale factor for the moment

    // Gravitational load (also uses coordSystemID and directionVector)
    feaLoad->gravityAcceleration= 0; // Gravitational acceleration

    // Pressure load
    feaLoad->pressureForce= 0; // Pressure value

    feaLoad->pressureDistributeForce[0] = 0; // Pressure load at a specified grid location in the element
    feaLoad->pressureDistributeForce[1] = 0;
    feaLoad->pressureDistributeForce[2] = 0;
    feaLoad->pressureDistributeForce[3] = 0;

    // Unique pressure load at a specified grid location for
    //each element in elementIDSet size = [4*numElementID]- used in type PressureExternal
    AIM_FREE(feaLoad->pressureMultiDistributeForce);


    feaLoad->numElementID = 0; // Number of elements IDs in element ID set
    AIM_FREE(feaLoad->elementIDSet); // List element IDs in which to apply the load

    // Rotational velocity (also uses coordSystemID and directionVector)
    feaLoad->angularVelScaleFactor = 0.0; // Overall scale factor for the angular velocity
    feaLoad->angularAccScaleFactor = 0.0; // Overall scale factor for the angular acceleration

    // Thermal load - the temperature at a grid point - use gridIDSet
    feaLoad->temperature = 0.0; // Temperature value
    feaLoad->temperatureDefault = 0.0; // Default temperature of grid point explicitly not used

    AIM_FREE(feaLoad->temperatureMultiDistribute);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaDesignVariable in the feaDesignVariableStruct structure format
int initiate_feaDesignVariableStruct(feaDesignVariableStruct *feaDesignVariable) {

    feaDesignVariable->name = NULL;

    feaDesignVariable->designVariableID = 0; //  ID number of design variable


    feaDesignVariable->initialValue = 0.0; // Initial value of design variable
    feaDesignVariable->lowerBound = 0.0;   // Lower bounds of variable
    feaDesignVariable->upperBound = 0.0;   // Upper bounds of variable
    feaDesignVariable->maxDelta = 0.5;     // Change fraction limit for design variable

    feaDesignVariable->numDiscreteValue = 0; // Number of discrete values that a design variable can assume

    feaDesignVariable->discreteValue = NULL; // List of discrete values that a design variable can assume;

    feaDesignVariable->numIndependVariable = 0;  // Number of independent variables this variables depends on
    feaDesignVariable->independVariable = NULL; // List of independent variable names, size[numIndependVariable]
    feaDesignVariable->independVariableID = NULL;// List of independent variable designVariableIDs, size[numIndependVariable]
    feaDesignVariable->independVariableWeight = NULL; // List of independent variable weights, size[numIndependVariable]

    feaDesignVariable->variableWeight[0] = 0.0; // Weight to apply to if variable is dependent
    feaDesignVariable->variableWeight[1] = 0.0;

    feaDesignVariable->numRelation = 0;
    feaDesignVariable->relationSet = NULL;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaDesignVariable in the feaDesignVariableStruct structure format
int destroy_feaDesignVariableStruct(feaDesignVariableStruct *feaDesignVariable) {

    if (feaDesignVariable->name != NULL) EG_free(feaDesignVariable->name);
    feaDesignVariable->name = NULL;

    feaDesignVariable->designVariableID = 0; //  ID number of design variable

    feaDesignVariable->initialValue = 0.0; // Initial value of design variable
    feaDesignVariable->lowerBound = 0.0;   // Lower bounds of variable
    feaDesignVariable->upperBound = 0.0;   // Upper bounds of variable
    feaDesignVariable->maxDelta = 0.5;     // Change fraction limit for design variable

    feaDesignVariable->numDiscreteValue = 0; // Number of discrete values that a design variable can assume;
    if (feaDesignVariable->discreteValue != NULL) EG_free(feaDesignVariable->discreteValue);
    feaDesignVariable->discreteValue = NULL; // List of discrete values that a design variable can assume;

    (void) string_freeArray(feaDesignVariable->numIndependVariable, &feaDesignVariable->independVariable);
    feaDesignVariable->independVariable = NULL; // List of independent variable names, size[numIndependVariable]

    feaDesignVariable->numIndependVariable = 0;  // Number of independent variables this variables depends on

    if (feaDesignVariable->independVariableID != NULL) EG_free(feaDesignVariable->independVariableID);
    feaDesignVariable->independVariableID = NULL;// List of independent variable designVariableIDs

    if (feaDesignVariable->independVariableWeight != NULL) EG_free(feaDesignVariable->independVariableWeight);
    feaDesignVariable->independVariableWeight = NULL;  // List of independent variable weights, size[numIndependVariable]

    feaDesignVariable->variableWeight[0] = 0.0; // Weight to apply to if variable is dependent
    feaDesignVariable->variableWeight[1] = 0.0;

    feaDesignVariable->numRelation = 0;
    if (feaDesignVariable->relationSet != NULL) EG_free(feaDesignVariable->relationSet);

    return CAPS_SUCCESS;

}

// Initiate (0 out all values and NULL all pointers) of feaDesignConstraint in the feaDesignConstraintStruct structure format
int initiate_feaDesignConstraintStruct(feaDesignConstraintStruct *feaDesignConstraint) {

    feaDesignConstraint->name = NULL;

    feaDesignConstraint->designConstraintID = 0; //  ID number of design constraint

    feaDesignConstraint->designConstraintType = UnknownDesignCon;

    feaDesignConstraint->responseType = NULL;  // Response type options for DRESP1 Entry

    feaDesignConstraint->lowerBound = 0.0;   // Lower bounds of design response
    feaDesignConstraint->upperBound = 0.0;   // Upper bounds of design response

    feaDesignConstraint->numPropertyID = 0;   // Number of property ID to apply the design variable to
    feaDesignConstraint->propertySetID = NULL; // List of property IDs
    feaDesignConstraint->propertySetType = NULL; // List of property types corresponding to the propertySetID

    feaDesignConstraint->fieldPosition = 0; //  Position in card to apply design variable to
    feaDesignConstraint->fieldName = NULL; // Name of property/material to apply design variable to

    feaDesignConstraint->velocityType = NULL;
    feaDesignConstraint->scalingFactor = 0.0;

    feaDesignConstraint->numDensity = 0;
    feaDesignConstraint->density = NULL;

    feaDesignConstraint->numMach = 0;
    feaDesignConstraint->Mach = NULL;

    feaDesignConstraint->numModes = 0;
    feaDesignConstraint->modes = NULL;

    feaDesignConstraint->numVelocity = 0;
    feaDesignConstraint->velocity = NULL;

    feaDesignConstraint->numDamping = 0;
    feaDesignConstraint->damping = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaDesignConstraint in the feaDesignConstraintStruct structure format
int destroy_feaDesignConstraintStruct(feaDesignConstraintStruct *feaDesignConstraint) {

    if (feaDesignConstraint->name != NULL) EG_free(feaDesignConstraint->name);
    feaDesignConstraint->name = NULL;

    feaDesignConstraint->designConstraintID = 0; //  ID number of design constraint

    feaDesignConstraint->designConstraintType = UnknownDesignCon;

    if (feaDesignConstraint->responseType != NULL) EG_free(feaDesignConstraint->responseType);
    feaDesignConstraint->responseType = NULL;  // Response type options for DRESP1 Entry

    feaDesignConstraint->lowerBound = 0.0;   // Lower bounds of design response
    feaDesignConstraint->upperBound = 0.0;   // Upper bounds of design response

    feaDesignConstraint->numPropertyID = 0;   // Number of property ID to apply the design variable to

    if (feaDesignConstraint->propertySetID != NULL) EG_free(feaDesignConstraint->propertySetID);
    feaDesignConstraint->propertySetID = NULL; // List of property IDs

    if (feaDesignConstraint->propertySetType != NULL) EG_free(feaDesignConstraint->propertySetType);
    feaDesignConstraint->propertySetType = NULL; // List of property types corresponding to the propertySetID

    feaDesignConstraint->fieldPosition = 0; //  Position in card to apply design variable to
    if (feaDesignConstraint->fieldName != NULL) EG_free(feaDesignConstraint->fieldName);
    feaDesignConstraint->fieldName = NULL; // Name of property/material to apply design variable to

    if (feaDesignConstraint->velocityType != NULL) EG_free(feaDesignConstraint->velocityType);
    feaDesignConstraint->velocityType = NULL;
    feaDesignConstraint->scalingFactor = 0.0;

    if (feaDesignConstraint->density != NULL) EG_free(feaDesignConstraint->density);
    feaDesignConstraint->numDensity = 0;
    feaDesignConstraint->density = NULL;

    if (feaDesignConstraint->Mach != NULL) EG_free(feaDesignConstraint->Mach);
    feaDesignConstraint->numMach = 0;
    feaDesignConstraint->Mach = NULL;

    if (feaDesignConstraint->modes != NULL) EG_free(feaDesignConstraint->modes);
    feaDesignConstraint->numModes = 0;
    feaDesignConstraint->modes = NULL;

    if (feaDesignConstraint->velocity != NULL) EG_free(feaDesignConstraint->velocity);
    feaDesignConstraint->numVelocity = 0;
    feaDesignConstraint->velocity = NULL;

    if (feaDesignConstraint->damping != NULL) EG_free(feaDesignConstraint->damping);
    feaDesignConstraint->numDamping = 0;
    feaDesignConstraint->damping = NULL;

    return CAPS_SUCCESS;

}

// Initiate (0 out all values and NULL all pointers) of the feaMassIncrementStruct structure format
int initiate_feaMassIncrementStruct(feaMassIncrementStruct *feaMassIncrement) {

    feaMassIncrement->name = NULL;

    feaMassIncrement->massIncrementID = 0; //  ID number of design constraint

    feaMassIncrement->numPropertyID = 0;   // Number of property ID to apply the design variable to
    feaMassIncrement->propertySetID = NULL; // List of property IDs
    feaMassIncrement->propertySetType = NULL; // List of property types corresponding to the propertySetID

    feaMassIncrement->increment = 0.0;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of the feaMassIncrementStruct structure format
int destroy_feaMassIncrementStruct(feaMassIncrementStruct *feaMassIncrement) {

    if (feaMassIncrement->name != NULL) EG_free(feaMassIncrement->name);
    feaMassIncrement->name = NULL;

    feaMassIncrement->massIncrementID = 0; //  ID number of design constraint

    feaMassIncrement->numPropertyID = 0;   // Number of property ID to apply the design variable to

    if (feaMassIncrement->propertySetID != NULL) EG_free(feaMassIncrement->propertySetID);
    feaMassIncrement->propertySetID = NULL; // List of property IDs

    if (feaMassIncrement->propertySetType != NULL) EG_free(feaMassIncrement->propertySetType);
    feaMassIncrement->propertySetType = NULL; // List of property types corresponding to the propertySetID

    feaMassIncrement->increment = 0.0;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaOptimizationControl in the feaOptimizationControlStruct structure format
int initiate_feaOptimizationControlStruct(feaOptimizationControlStruct *feaOptimizationControl) {

    feaOptimizationControl->fullyStressedDesign = 0;    // ASTROS specific - fully stressed design number of iterations
    feaOptimizationControl->mathProgramming = 0;        // Math programming number of iteratinons
    feaOptimizationControl->maxIter = 0;                // Maximum number of optimization iterations
    feaOptimizationControl->constraintRetention = 0.0;  // Constraint retention factor
    feaOptimizationControl->eps = 0.0;                  // A different constraint retention factor?
    feaOptimizationControl->moveLimit = 0.0;            // Move limit
    feaOptimizationControl->string = NULL;              // Manual string override

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaOptimizationControl in the feaOptimizationControlStruct structure format
int destroy_feaOptimzationControlStruct(feaOptimizationControlStruct *feaOptimizationControl) {

    feaOptimizationControl->fullyStressedDesign = 0;    // ASTROS specific - fully stressed design number of iterations
    feaOptimizationControl->mathProgramming = 0;        // Math programming number of iteratinons
    feaOptimizationControl->maxIter = 0;                // Maximum number of optimization iterations
    feaOptimizationControl->constraintRetention = 0.0;  // Constraint retention factor
    feaOptimizationControl->eps = 0.0;                  // A different constraint retention factor?
    feaOptimizationControl->moveLimit = 0.0;            // Move limit
    feaOptimizationControl->string = NULL;               // Manual string override

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaCoordSystem in the feaCoordSystemStruct structure format
int initiate_feaCoordSystemStruct(feaCoordSystemStruct *feaCoordSystem) {

    int i; // Indexing

    feaCoordSystem->name = NULL; // Coordinate system name

    feaCoordSystem->coordSystemType = UnknownCoordSystem;  // Coordinate system type

    feaCoordSystem->coordSystemID = 0; // ID number of coordinate system
    feaCoordSystem->refCoordSystemID = 0; // ID of reference coordinate system

    for (i = 0; i < 3; i++) {
        feaCoordSystem->origin[i] = 0; // x, y, and z coordinates for the origin

        feaCoordSystem->normal1[i] = 0; // First normal direction
        feaCoordSystem->normal2[i] = 0; // Second normal direction
        feaCoordSystem->normal3[i] = 0; // Third normal direction - found from normal1 x normal2
    }

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaCoordSystem in the feaCoordSystemStruct structure format
int destroy_feaCoordSystemStruct(feaCoordSystemStruct *feaCoordSystem) {

    int i; // Indexing

    if (feaCoordSystem->name != NULL) EG_free(feaCoordSystem->name);
    feaCoordSystem->name = NULL; // Coordinate system name

    feaCoordSystem->coordSystemType = UnknownCoordSystem;  // Coordinate system type

    feaCoordSystem->coordSystemID = 0; // ID number of coordinate system
    feaCoordSystem->refCoordSystemID = 0; // ID of reference coordinate system

    for (i = 0; i < 3; i++) {
        feaCoordSystem->origin[i] = 0; // x, y, and z coordinates for the origin

        feaCoordSystem->normal1[i] = 0; // First normal direction
        feaCoordSystem->normal2[i] = 0; // Second normal direction
        feaCoordSystem->normal3[i] = 0; // Third normal direction - found from normal1 x normal2
    }

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaAero in the feaAeroStruct structure format
int initiate_feaAeroStruct(feaAeroStruct *feaAero) {

    int status; // Function return status

    feaAero->name = NULL; // Coordinate system name

    feaAero->surfaceID = 0; // Surface ID
    feaAero->coordSystemID = 0; // Coordinate system ID

    feaAero->numGridID = 0; // Number of grid IDs in grid ID set for the spline
    feaAero->gridIDSet = NULL; // List of grid IDs to apply spline to. size = [numGridID]

    status = initiate_vlmSurfaceStruct(&feaAero->vlmSurface);
    if (status != CAPS_SUCCESS) printf("Status %d during initiate_vlmSurfaceStruct\n", status);

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaAero in the feaAeroStruct structure format
int destroy_feaAeroStruct(feaAeroStruct *feaAero) {

    int status; // Function return status

    if (feaAero->name != NULL) EG_free(feaAero->name);
    feaAero->name = NULL; // Coordinate system name

    feaAero->surfaceID = 0; // Surface ID
    feaAero->coordSystemID = 0; // Coordinate system ID

    feaAero->numGridID = 0; // Number of grid IDs in grid ID set for the spline

    if (feaAero->gridIDSet != NULL) EG_free(feaAero->gridIDSet);
    feaAero->gridIDSet = NULL; // List of grid IDs to apply spline to. size = [numGridID]

    status = destroy_vlmSurfaceStruct(&feaAero->vlmSurface);
    if (status != CAPS_SUCCESS) printf("Status %d during destroy_vlmSurfaceStruct\n", status);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaAeroRef in the feaAeroRefStruct structure format
int initiate_feaAeroRefStruct(feaAeroRefStruct *feaAeroRef) {

    feaAeroRef->coordSystemID = 0; // Aerodynamic coordinate sytem id
    feaAeroRef->rigidMotionCoordSystemID = 0; // Reference coordinate system identification for rigid body motions.

    feaAeroRef->refChord = 1.0; // Reference chord length.  Reference span.  (Real > 0.0)
    feaAeroRef->refSpan = 1.0; // Reference span
    feaAeroRef->refArea = 1.0; // Reference area
    feaAeroRef->refGridID = 0;

    feaAeroRef->refVelocity = 1.0; // Reference Velocity
    feaAeroRef->refDensity = 1.0; // Reference Density

    feaAeroRef->symmetryXZ  = 0; // Symmetry key for the aero coordinate x-z plane.  (Integer = +1 for symmetry, 0 for no symmetry,
                                 // and -1 for antisymmetry; Default = 0)
    feaAeroRef->symmetryXY = 0; // The symmetry key for the aero coordinate x-y plane can be used to simulate ground effects.
                                // (Integer = +1 for antisymmetry, 0 for no symmetry, and -1 for symmetry; Default = 0)
    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaAeroRef in the feaAeroRefStruct structure format
int destroy_feaAeroRefStruct(feaAeroRefStruct *feaAeroRef) {

    feaAeroRef->coordSystemID = 0; // Aerodynamic coordinate sytem id
    feaAeroRef->rigidMotionCoordSystemID = 0; // Reference coordinate system identification for rigid body motions.

    feaAeroRef->refChord = 0; // Reference chord length.  Reference span.  (Real > 0.0)
    feaAeroRef->refSpan = 0; // Reference span
    feaAeroRef->refArea = 0; // Reference area
    feaAeroRef->refGridID = 0;

    feaAeroRef->refVelocity = 0; // Reference Velocity
    feaAeroRef->refDensity = 0; // Reference Density

    feaAeroRef->symmetryXZ  = 0; // Symmetry key for the aero coordinate x-z plane.  (Integer = +1 for symmetry, 0 for no symmetry,
                                 // and -1 for antisymmetry; Default = 0)
    feaAeroRef->symmetryXY = 0; // The symmetry key for the aero coordinate x-y plane can be used to simulate ground effects.
                                // (Integer = +1 for antisymmetry, 0 for no symmetry, and -1 for symmetry; Default = 0)
    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaMassProp in the feaMassPropStruct structure format
int initiate_feaMassPropStruct(feaMassPropStruct *feaMassProp) {

    int i;

    feaMassProp->mass = 0;
    for (i = 0; i < 3; i++) feaMassProp->CG[i] = 0;
    for (i = 0; i < 6; i++) feaMassProp->massInertia[i] = 0;

    feaMassProp->massUnit = NULL;
    feaMassProp->lengthUnit = NULL;
    feaMassProp->momentOfInertiaUnit = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaMassProp in the feaMassPropStruct structure format
int destroy_feaMassPropStruct(feaMassPropStruct *feaMassProp) {

    AIM_FREE(feaMassProp->massUnit);
    AIM_FREE(feaMassProp->lengthUnit);
    AIM_FREE(feaMassProp->momentOfInertiaUnit);

    initiate_feaMassPropStruct(feaMassProp);

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaConnect in the feaConnectionStruct structure format
int initiate_feaConnectionStruct(feaConnectionStruct *feaConnect) {

    feaConnect->name = NULL; // Connection name

    feaConnect->connectionID =0; // Connection ID

    feaConnect->connectionType = UnknownConnection; // Connection type

    feaConnect->elementID = 0;

    // RBE2 - dependent degrees of freedom
    feaConnect->connectivity[0] = 0; // Grid IDs - 0 index = Independent grid ID, 1 index = Dependent grid ID
    feaConnect->connectivity[1] = 0;
    feaConnect->dofDependent = 0;

    // Spring
    feaConnect->stiffnessConst  = 0.0;
    feaConnect->componentNumberStart = 0;
    feaConnect->componentNumberEnd = 0;
    feaConnect->dampingConst  = 0.0;
    feaConnect->stressCoeff = 0.0;

    // Damper - see spring for additional entries

    // Mass (scalar) - see spring for additional entries
    feaConnect->mass = 0.0;

    // RBE3 - master/slave
    feaConnect->numMaster = 0;
    feaConnect->masterIDSet = NULL; // Independent
    feaConnect->masterWeighting = NULL;
    feaConnect->masterComponent = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaConnect in the feaConnectionStruct structure format
int destroy_feaConnectionStruct(feaConnectionStruct *feaConnect) {

    if (feaConnect->name != NULL) EG_free(feaConnect->name);
    feaConnect->name = NULL; // Connection name

    feaConnect->connectionID =0; // Connection ID

    feaConnect->connectionType = UnknownConnection; // Connection type

    feaConnect->elementID = 0;

    // RBE2 - dependent degrees of freedom
    feaConnect->connectivity[0] = 0; // Grid IDs - 0 index = Independent grid ID, 1 index = Dependent grid ID
    feaConnect->connectivity[1] = 0;
    feaConnect->dofDependent = 0;

    // Spring
    feaConnect->stiffnessConst = 0.0;
    feaConnect->componentNumberStart = 0;
    feaConnect->componentNumberEnd = 0;
    feaConnect->dampingConst  = 0.0;
    feaConnect->stressCoeff = 0.0;

    // Damper - see spring for additional entries

    // Mass (scalar) - see spring for additional entries
    feaConnect->mass = 0.0;

    // RBE3 - master/slave
    feaConnect->numMaster = 0;
    if (feaConnect->masterIDSet != NULL) EG_free(feaConnect->masterIDSet);
    feaConnect->masterIDSet = NULL; // Independent
    if (feaConnect->masterWeighting != NULL) EG_free(feaConnect->masterWeighting);
    feaConnect->masterWeighting = NULL;
    if (feaConnect->masterComponent != NULL) EG_free(feaConnect->masterComponent);
    feaConnect->masterComponent = NULL;

    return CAPS_SUCCESS;
}

int initiate_feaDesignEquationStruct(feaDesignEquationStruct *equation) {

    if (equation == NULL) return CAPS_NULLVALUE;

    equation->equationID = 0;
    equation->name = NULL;
    equation->equationArraySize = 0;
    equation->equationArray = NULL;

    return CAPS_SUCCESS;
}

int destroy_feaDesignEquationStruct(feaDesignEquationStruct *equation) {

    int i;

    if (equation == NULL) return CAPS_NULLVALUE;

    if (equation->name != NULL) EG_free(equation->name);

    if (equation->equationArray != NULL) {
        for (i = 0; i < equation->equationArraySize; i++) {
            if (equation->equationArray[i] != NULL) {
                EG_free(equation->equationArray[i]);
            }
        }
        EG_free(equation->equationArray);
    }

    return initiate_feaDesignEquationStruct(equation);
}

int initiate_feaDesignResponseStruct(feaDesignResponseStruct *response) {

    if (response == NULL) return CAPS_NULLVALUE;

    response->responseID = 0;

    response->name = NULL;

    response->responseType = NULL;
    response->propertyType = NULL;

    response->region = 0;

    response->component = 0; // Component number
    response->itemCode = 0; // Item code
    response->attb = 0;

    response->lamina = 0; // Lamina number
    response->frequency = 0.0; // Frequency value
    response->time = 0.0; // Time value
    response->restraintFlag = 0; // Restraint flag

    response->gridID = 0; // Grid ID
    response->propertyID = 0; // Property entry ID

    return CAPS_SUCCESS;
}

int destroy_feaDesignResponseStruct(feaDesignResponseStruct *response) {

    if (response == NULL) return CAPS_NULLVALUE;

    AIM_FREE(response->name);

    AIM_FREE(response->responseType);
    AIM_FREE(response->propertyType);

    return initiate_feaDesignResponseStruct(response);
}

int initiate_feaDesignEquationResponseStruct(feaDesignEquationResponseStruct* equationResponse) {

    if (equationResponse == NULL) return CAPS_NULLVALUE;

    equationResponse->equationResponseID = 0;

    equationResponse->name = NULL; // Name of the equation

    equationResponse->equationName = NULL;

    equationResponse->region = 0; // Region identifier for constant screening

    equationResponse->numDesignVariable = 0;
    equationResponse->designVariableNameSet = NULL; // Design variable names, size = [numDesignVariable]

    equationResponse->numConstant = 0;
    equationResponse->constantLabelSet = NULL; // Labels of the table constants, size = [numConstant]

    equationResponse->numResponse = 0;
    equationResponse->responseNameSet = NULL; // Names of design sensitivity response quantities, size = [numResponse]

    equationResponse->numGrid = 0;
    equationResponse->gridIDSet = NULL; // Grid IDs, size = [numGrid]
    equationResponse->dofNumberSet = NULL; // Degree of freedom numbers, size = [numGrid]

    equationResponse->numEquationResponse = 0;
    equationResponse->equationResponseNameSet = NULL; // Names of design sensitivity equation response quantities, size = [numEquationResponse]

    return CAPS_SUCCESS;
}

int destroy_feaDesignEquationResponseStruct(feaDesignEquationResponseStruct* equationResponse) {

    if (equationResponse == NULL) return CAPS_NULLVALUE;

    if (equationResponse->name != NULL) EG_free(equationResponse->name);

    if (equationResponse->equationName != NULL) EG_free(equationResponse->equationName);

    if (equationResponse->designVariableNameSet != NULL) {
        string_freeArray(equationResponse->numDesignVariable, &equationResponse->designVariableNameSet);
    }

    if (equationResponse->constantLabelSet != NULL) {
        string_freeArray(equationResponse->numConstant, &equationResponse->constantLabelSet);
    }

    if (equationResponse->responseNameSet != NULL) {
        string_freeArray(equationResponse->numResponse, &equationResponse->responseNameSet);
    }

    if (equationResponse->gridIDSet != NULL) EG_free(equationResponse->gridIDSet);
    if (equationResponse->dofNumberSet != NULL) EG_free(equationResponse->dofNumberSet);

    if (equationResponse->equationResponseNameSet != NULL) {
        string_freeArray(equationResponse->numEquationResponse, &equationResponse->equationResponseNameSet);
    }

    return initiate_feaDesignEquationResponseStruct(equationResponse);
}

int initiate_feaDesignTableStruct(feaDesignTableStruct *table) {

    if (table == NULL) return CAPS_NULLVALUE;

    table->numConstant = 0;
    table->constantLabel = NULL;
    table->constantValue = NULL;

    return CAPS_SUCCESS;
}

int destroy_feaDesignTableStruct(feaDesignTableStruct *table) {

    if (table == NULL) return CAPS_NULLVALUE;

    if (table->constantLabel != NULL) string_freeArray(table->numConstant, &table->constantLabel);
    if (table->constantValue != NULL) EG_free(table->constantValue);

    return initiate_feaDesignTableStruct(table);
}


int initiate_feaDesignOptParamStruct(feaDesignOptParamStruct *table) {

    if (table == NULL) return CAPS_NULLVALUE;

    table->numParam = 0;
    table->paramLabel = NULL;
    table->paramValue = NULL;
    table->paramType = NULL;

    return CAPS_SUCCESS;
}

int destroy_feaDesignOptParamStruct(feaDesignOptParamStruct *table) {

    int i;

    if (table == NULL) return CAPS_NULLVALUE;

    if (table->paramLabel != NULL) string_freeArray(table->numParam, &table->paramLabel);

    if (table->paramValue != NULL) {
        for (i = 0; i < table->numParam; i++) {
            if (table->paramValue[i] != NULL) {
                EG_free(table->paramValue[i]);
            }
        }
        EG_free(table->paramValue);
    }

    if (table->paramType != NULL) EG_free(table->paramType);

    return initiate_feaDesignOptParamStruct(table);
}

int initiate_feaDesignVariableRelationStruct(feaDesignVariableRelationStruct *relation) {

    if (relation == NULL) return CAPS_NULLVALUE;

    relation->name = NULL;

    relation->componentType = 0;

    relation->relationID = 0;

    relation->numDesignVariable = 0;
    relation->designVariableNameSet = NULL;
    relation->designVariableSet = NULL;

    relation->fieldPosition = 0;
    relation->fieldName = NULL;

    relation->isEquation = 0;
    relation->equationName = NULL;

    relation->constantRelationCoeff = 0.0;
    relation->linearRelationCoeff = NULL;

    relation->numMaterialID = 0; // Number of materials to apply the design variable to
    relation->materialSetID = NULL; // List of materials IDs
    relation->materialSetType = NULL; // List of materials types corresponding to the materialSetID

    relation->numPropertyID = 0;   // Number of property ID to apply the design variable to
    relation->propertySetID = NULL; // List of property IDs
    relation->propertySetType = NULL; // List of property types corresponding to the propertySetID

    relation->numElementID = 0; // Number of element ID to apply the design variable to
    relation->elementSetID = NULL; // List of element IDs
    relation->elementSetType = NULL; // List of element types corresponding to the elementSetID
    relation->elementSetSubType = NULL; // List of element subtypes correspoding to the elementSetID

    return CAPS_SUCCESS;
}

int destroy_feaDesignVariableRelationStruct(feaDesignVariableRelationStruct *relation) {

    if (relation == NULL) return CAPS_NULLVALUE;

    if (relation->designVariableNameSet != NULL) string_freeArray(relation->numDesignVariable, &relation->designVariableNameSet);
    if (relation->name != NULL) EG_free(relation->name);
    if (relation->designVariableSet != NULL) EG_free(relation->designVariableSet);

    if (relation->fieldName != NULL) EG_free(relation->fieldName);

    if (relation->linearRelationCoeff != NULL) EG_free(relation->linearRelationCoeff);

    relation->numMaterialID = 0; // Number of materials to apply the design variable to
    if (relation->materialSetID != NULL) EG_free(relation->materialSetID);
    relation->materialSetID = NULL; // List of materials IDs

    if (relation->materialSetType != NULL) EG_free(relation->materialSetType);
    relation->materialSetType = NULL; // List of materials types corresponding to the materialSetID

    relation->numPropertyID = 0;   // Number of property ID to apply the design variable to
    if (relation->propertySetID != NULL) EG_free(relation->propertySetID);
    relation->propertySetID = NULL; // List of property IDs

    if (relation->propertySetType != NULL) EG_free(relation->propertySetType);
    relation->propertySetType = NULL; // List of property types corresponding to the propertySetID

    relation->numElementID = 0; // Number of element ID to apply the design variable to
    if (relation->elementSetID != NULL) EG_free(relation->elementSetID);
    relation->elementSetID = NULL; // List of element IDs

    if (relation->elementSetType != NULL) EG_free(relation->elementSetType);
    relation->elementSetType = NULL; // List of element types corresponding to the elementSetID

    if (relation->elementSetSubType != NULL) EG_free(relation->elementSetSubType);
    relation->elementSetSubType = NULL; // List of element subtypes correspoding to the elementSetID

    return initiate_feaDesignVariableRelationStruct(relation);
}

// Get the material properties from a capsTuple
int fea_getMaterial(void *aimInfo,
                    int numMaterialTuple,
                    capsTuple materialTuple[],
                    feaUnitsStruct *feaUnits,
                    int *numMaterial,
                    feaMaterialStruct *feaMaterial[]) {

    /*! \page feaMaterial FEA Material
     * Structure for the material tuple  = ("Material Name", "Value").
     * "Material Name" defines the reference name for the material being specified.
     * The "Value" can either be a JSON String dictionary (see Section \ref jsonStringMaterial) or a single string keyword
     * (see Section \ref keyStringMaterial).
     */

    int status; //Function return

    int i; // Indexing

    char *keyValue = NULL;
    char *keyWord = NULL;

    // Destroy our material structures coming in if aren't 0 and NULL already
    if (*feaMaterial != NULL) {
        for (i = 0; i < *numMaterial; i++) {
            status = destroy_feaMaterialStruct(&(*feaMaterial)[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (*feaMaterial != NULL) EG_free(*feaMaterial);
    *feaMaterial = NULL;
    *numMaterial = 0;

    printf("\nGetting FEA materials.......\n");

    *numMaterial = numMaterialTuple;
    printf("\tNumber of materials - %d\n", *numMaterial);

    if (*numMaterial > 0) {
        *feaMaterial = (feaMaterialStruct *) EG_alloc(*numMaterial * sizeof(feaMaterialStruct));
        if (*feaMaterial == NULL) return EGADS_MALLOC;

    } else {
        printf("\tNumber of material values in input tuple is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < *numMaterial; i++) {
        status = initiate_feaMaterialStruct(&(*feaMaterial)[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    for (i = 0; i < *numMaterial; i++) {

        printf("\tMaterial name - %s\n", materialTuple[i].name);

        (*feaMaterial)[i].name = (char *) EG_alloc(((strlen(materialTuple[i].name)) + 1)*sizeof(char));
        if ((*feaMaterial)[i].name == NULL) return EGADS_MALLOC;

        memcpy((*feaMaterial)[i].name, materialTuple[i].name, strlen(materialTuple[i].name)*sizeof(char));
        (*feaMaterial)[i].name[strlen(materialTuple[i].name)] = '\0';

        (*feaMaterial)[i].materialID = i + 1;

        // Do we have a json string?
        if (strncmp(materialTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", materialTuple[i].value);


            /*! \page feaMaterial
             * \section jsonStringMaterial JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS || MASSTRAN)
             *  (e.g. "Value" = {"density": 7850, "youngModulus": 120000.0, "poissonRatio": 0.5, "materialType": "isotropic"})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * <ul>
             * <li> <B>materialType = "Isotropic"</B> </li> <br>
             * \if (MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS)
             *      Material property type. Options: Isotropic, Anisothotropic, Orthotropic, or Anisotropic.
             * \elseif MASSTRAN
             *      Material property type. Options: Isotropic.
             * \endif
             * </ul>
             */

            // Get material Type
            keyWord = "materialType";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                //{UnknownMaterial, Isotropic, Anisothotropic, Orthotropic, Anisotropic}
                if      (strcasecmp(keyValue, "\"Isotropic\"")      == 0) (*feaMaterial)[i].materialType = Isotropic;
                else if (strcasecmp(keyValue, "\"Anisothotropic\"") == 0) (*feaMaterial)[i].materialType = Anisothotropic;
                else if (strcasecmp(keyValue, "\"Orthotropic\"")    == 0) (*feaMaterial)[i].materialType = Orthotropic;
                else if (strcasecmp(keyValue, "\"Anisotropic\"")    == 0) (*feaMaterial)[i].materialType = Anisotropic;
                else {

                    printf("\tUnrecognized \"%s\" specified (%s) for Material tuple %s, defaulting to \"Isotropic\"\n", keyWord,
                            keyValue,
                            materialTuple[i].name);
                    (*feaMaterial)[i].materialType = Isotropic;
                }

            } else {

                printf("\tNo \"%s\" specified for Material tuple %s, defaulting to \"Isotropic\"\n", keyWord,
                        materialTuple[i].name);
                (*feaMaterial)[i].materialType = Isotropic;
            }
            AIM_FREE(keyValue);

            //Fill up material properties

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS)
             *  <ul>
             *	<li> <B>youngModulus = 0.0</B> </li> <br>
             *  Also known as the elastic modulus, defines the relationship between stress and strain.
             *  Default if `shearModulus' and `poissonRatio' != 0, youngModulus = 2*(1+poissonRatio)*shearModulus
             *  </ul>
             * \endif
             *
             */
            keyWord = "youngModulus";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                if ( feaUnits->pressure != NULL ) {
                    status = string_toDoubleUnits(aimInfo, keyValue, feaUnits->pressure, &(*feaMaterial)[i].youngModulus);
                } else {
                    status = string_toDouble(keyValue, &(*feaMaterial)[i].youngModulus);
                }
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS)
             *  <ul>
             *  <li> <B>shearModulus = 0.0</B> </li> <br>
             *  Also known as the modulus of rigidity, is defined as the ratio of shear stress to the shear strain.
             *  Default if `youngModulus' and `poissonRatio' != 0, shearModulus = youngModulus/(2*(1+poissonRatio))
             *  </ul>
             * \endif
             */
            keyWord = "shearModulus";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearModulus);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS)
             *  <ul>
             *  <li> <B>poissonRatio = 0.0</B> </li> <br>
             *  The fraction of expansion divided by the fraction of compression.
             *  Default if `youngModulus' and `shearModulus' != 0, poissonRatio = (2*youngModulus/shearModulus) - 1
             *  </ul>
             * \endif
             */
            keyWord = "poissonRatio";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].poissonRatio);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            // Check Young's modulus, shear modulus, and Poisson's ratio
            if ( ((*feaMaterial)[i].youngModulus == 0 && (*feaMaterial)[i].poissonRatio == 0) ||
                 ((*feaMaterial)[i].shearModulus == 0 && (*feaMaterial)[i].poissonRatio == 0) ||
                 ((*feaMaterial)[i].youngModulus == 0 && (*feaMaterial)[i].shearModulus == 0) ) {
                // Do nothing

            } else if ((*feaMaterial)[i].youngModulus == 0) {
                (*feaMaterial)[i].youngModulus = 2*(1+(*feaMaterial)[i].poissonRatio)*(*feaMaterial)[i].shearModulus;
            } else if ((*feaMaterial)[i].shearModulus == 0) {
                (*feaMaterial)[i].shearModulus = (*feaMaterial)[i].youngModulus/(2*(1+(*feaMaterial)[i].poissonRatio));
            } else if ((*feaMaterial)[i].poissonRatio == 0) {
                (*feaMaterial)[i].poissonRatio = (*feaMaterial)[i].youngModulus/(2*(*feaMaterial)[i].shearModulus) -1;
            }

            /*! \page feaMaterial
             *
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS || MASSTRAN)
             *  <ul>
             *  <li> <B>density = 0.0</B> </li> <br>
             *  Density of the material.
             *  </ul>
             * \endif
             */
            keyWord = "density";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                if ( feaUnits->densityVol != NULL ) {
                    status = string_toDoubleUnits(aimInfo, keyValue, feaUnits->densityVol, &(*feaMaterial)[i].density);
                } else {
                    status = string_toDouble(keyValue, &(*feaMaterial)[i].density);
                }
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  <ul>
             *  <li> <B>thermalExpCoeff = 0.0</B> </li> <br>
             *  Thermal expansion coefficient of the material.
             *  </ul>
             * \endif
             */
            keyWord = "thermalExpCoeff";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].thermalExpCoeff);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  <ul>
             *  <li> <B>thermalExpCoeffLateral = 0.0</B> </li> <br>
             *  Thermal expansion coefficient of the material.
             *  </ul>
             * \endif
             */
            keyWord = "thermalExpCoeffLateral";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].thermalExpCoeffLateral);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  <ul>
             *  <li> <B>temperatureRef = 0.0</B> </li> <br>
             *  Reference temperature for material properties.
             *  </ul>
             * \endif
             */
            keyWord = "temperatureRef";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].temperatureRef);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  <ul>
             *  <li> <B>dampingCoeff = 0.0</B> </li> <br>
             *  Damping coefficient for the material.
             *  </ul>
             * \endif
             */
            keyWord = "dampingCoeff";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].dampingCoeff);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             *
             * \if NASTRAN
             *  <ul>
             *  <li> <B>yieldAllow = 0.0</B> </li> <br>
             *  Yield strength/allowable for the material.
             *  </ul>
             * \endif
             */
            keyWord = "yieldAllow";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                status = string_toDouble(keyValue, &(*feaMaterial)[i].yieldAllow);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             *
             * \if NASTRAN
             *  <ul>
             *  <li> <B>tensionAllow = 0.0</B> </li> <br>
             *  Tension allowable for the material.
             *  </ul>
             * \endif
             */
            keyWord = "tensionAllow";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].tensionAllow);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             *
             * \if NASTRAN
             *  <ul>
             *  <li> <B>tensionAllowLateral = 0.0</B> </li> <br>
             *  Lateral tension allowable for the material.
             *  </ul>
             * \endif
             */
            keyWord = "tensionAllowLateral";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].tensionAllowLateral);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if NASTRAN
             *  <ul>
             *  <li> <B>compressAllow = 0.0</B> </li> <br>
             *  Compression allowable for the material.
             *  </ul>
             * \endif
             */
            keyWord = "compressAllow";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].compressAllow);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if NASTRAN
             *  <ul>
             *  <li> <B>compressAllowLateral = 0.0</B> </li> <br>
             *  Lateral compression allowable for the material.
             *  </ul>
             * \endif
             */
            keyWord = "compressAllowLateral";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].compressAllowLateral);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if NASTRAN
             *  <ul>
             *  <li> <B>shearAllow = 0.0</B> </li> <br>
             *  Shear allowable for the material.
             *  </ul>
             * \endif
             */
            keyWord = "shearAllow";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearAllow);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (NASTRAN)
             * <ul>
             *  <li> <B>allowType = 0 </B> </li> <br>
             *  This flag defines if the above allowables <c>compressAllow</c> etc. are defined in terms of stress (0) or strain (1).  The default is stress (0).
             * </ul>
             * \endif
             */
            keyWord = "allowType";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toInteger(keyValue, &(*feaMaterial)[i].allowType);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>youngModulusLateral = 0.0</B> </li> <br>
             *  Elastic modulus in lateral direction for an orthotropic material
             * </ul>
             * \endif
             */
            keyWord = "youngModulusLateral";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].youngModulusLateral);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>shearModulusTrans1Z = 0.0</B> </li> <br>
             *  Transverse shear modulus in the 1-Z plane for an orthotropic material
             * </ul>
             * \endif
             */
            keyWord = "shearModulusTrans1Z";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearModulusTrans1Z);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>shearModulusTrans2Z = 0.0</B> </li> <br>
             *  Transverse shear modulus in the 2-Z plane for an orthotropic material
             * </ul>
             * \endif
             */
            keyWord = "shearModulusTrans2Z";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearModulusTrans2Z);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>kappa = 0.0</B> </li> <br>
             *  Thermal conductivity for an isotropic solid
             * </ul>
             * \endif
             */
            keyWord = "kappa";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].kappa);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>K = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Thermal conductivity for an anisotropic solid (KXX, KXY, KXZ, KYY, KYZ, KZZ)
             * </ul>
             * \endif
             */
            keyWord = "K";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleArray(keyValue, 6, (*feaMaterial)[i].K);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>specificHeat = 0.0</B> </li> <br>
             *  Specific heat constant pressure (per unit mass) for an isotropic solid
             * </ul>
             * \endif
             */
            keyWord = "specificHeat";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].specificHeat);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B> Gij  = (no default) </B> </li> <br>
             *  List of Gij material properties (e.g. [G11, G12, G13, G22, fG23, G33]). Length must be 6.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             *
             * Something else ....
             *
             * \endif
             */

            keyWord = "gmat";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleArray(keyValue,
                                              6,
                                              (*feaMaterial)[i].gmat);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>honeycombCellSize = NULL</B> </li> <br>
             *  Honeycomb sandwich core cell size. Required if material defines the core of a honeycomb sandwich and dimpling stability index is desired
             * </ul>
             * \endif
             */
            keyWord = "honeycombCellSize";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].honeycombCellSize);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>honeycombYoungModulus = NULL</B> </li> <br>
             *  Honeycomb sandwich core Young's modulus used for stability index analysis
             * </ul>
             * \endif
             */
            keyWord = "honeycombYoungModulus";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].honeycombYoungModulus);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>honeycombShearModulus = NULL</B> </li> <br>
             *  Honeycomb sandwich core shear modulus used for stability index analysis
             * </ul>
             * \endif
             */
            keyWord = "honeycombShearModulus";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].honeycombShearModulus);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>fractureAngle = NULL</B> </li> <br>
             *  Fracture angle for uniaxial transverse compression in degrees. Used in the NASA LaRC02 failure theory only
             * </ul>
             * \endif
             */
            keyWord = "fractureAngle";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].fractureAngle);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>interlaminarShearAllow = NULL</B> </li> <br>
             *  Allowable inter-laminar shear stress of the composite laminate bonding material
             * </ul>
             * \endif
             */
            keyWord = "interlaminarShearAllow";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].interlaminarShearAllow);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>fiberYoungModulus = NULL</B> </li> <br>
             *  Modulus of elasticity of fiber
             * </ul>
             * \endif
             */
            keyWord = "fiberYoungModulus";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].fiberYoungModulus);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>fiberPoissonRatio = NULL</B> </li> <br>
             *  Poisson's ratio of fiber
             * </ul>
             * \endif
             */
            keyWord = "fiberPoissonRatio";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].fiberPoissonRatio);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>meanStressFactor = NULL</B> </li> <br>
             *  Mean stress magnification factor
             * </ul>
             * \endif
             */
            keyWord = "meanStressFactor";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].meanStressFactor);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>transTensionSlope = NULL</B> </li> <br>
             *  Failure envelop slope parameter for transverse tension
             * </ul>
             * \endif
             */
            keyWord = "transTensionSlope";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].transTensionSlope);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>transCompressionSlope = NULL</B> </li> <br>
             *  Failure envelop slope parameter for transverse compression
             * </ul>
             * \endif
             */
            keyWord = "transCompressionSlope";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].transCompressionSlope);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>compositeFailureTheory = NULL</B> </li> <br>
             *  Composite failure theory (string value)
             * </ul>
             * \endif
             */
            keyWord = "compositeFailureTheory";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                AIM_STRDUP((*feaMaterial)[i].compositeFailureTheory, keyValue, aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>interlaminarNormalStressAllow = NULL</B> </li> <br>
             *  Allowable inter-laminar normal stress of the composite laminate bonding material (allowable interlaminar normal stress)
             * </ul>
             * \endif
             */
            keyWord = "interlaminarNormalStressAllow";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].interlaminarNormalStressAllow);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>youngModulusThick = NULL</B> </li> <br>
             *  Modulus of elasticity in thickness direction, also defined as the matrix direction or 3-direction
             * </ul>
             * \endif
             */
            keyWord = "youngModulusThick";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].youngModulusThick);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>poissonRatio23 = 0.0</B> </li> <br>
             *  Poisson's ratio ( for uniaxial loading in 2-direction)
             * </ul>
             * \endif
             */
            keyWord = "poissonRatio23";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].poissonRatio23);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>poissonRatio31 = 0.0</B> </li> <br>
             *  Poisson's ratio ( for uniaxial loading in 3-direction)
             * </ul>
             * \endif
             */
            keyWord = "poissonRatio31";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].poissonRatio31);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>youngModulusFactor = NULL</B> </li> <br>
             *  Longitudinal modulus of elasticity reduction scale factor for nonlinear composite Progressive Ply Failure Analysis (PPFA)
             * </ul>
             * \endif
             */
            keyWord = "youngModulusFactor";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].youngModulusFactor);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>youngModulusLateralFactor = NULL</B> </li> <br>
             *  Lateral modulus of elasticity reduction scale factor for nonlinear composite Progressive Ply Failure Analysis (PPFA)
             * </ul>
             * \endif
             */
            keyWord = "youngModulusLateralFactor";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].youngModulusLateralFactor);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>shearModulusFactor = NULL</B> </li> <br>
             *  In-plane shear modulus reduction scale factor for nonlinear composite Progressive Ply Failure Analysis (PPFA)
             * </ul>
             * \endif
             */
            keyWord = "shearModulusFactor";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearModulusFactor);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>shearModulusTrans1ZFactor = NULL</B> </li> <br>
             *  Transverse shear modulus reduction scale factor in 1-Z plane for nonlinear composite Progressive Ply Failure Analysis (PPFA)
             * </ul>
             * \endif
             */
            keyWord = "shearModulusTrans1ZFactor";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearModulusTrans1ZFactor);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaMaterial
             *
             * \if (TACS || NASTRAN || ABAQUS)
             * <ul>
             *  <li> <B>shearModulusTrans2ZFactor = NULL</B> </li> <br>
             *  Transverse shear modulus reduction scale factor in 2-Z plane for nonlinear composite Progressive Ply Failure Analysis (PPFA)
             * </ul>
             * \endif
             */
            keyWord = "shearModulusTrans2ZFactor";
            status = search_jsonDictionary( materialTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &(*feaMaterial)[i].shearModulusTrans2ZFactor);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

        } else {

            /*! \page feaMaterial
             * \section keyStringMaterial Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined material lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */
            // CALL material look up
            AIM_ERROR(aimInfo, "Material tuple value ('%s') is expected to be a JSON string", materialTuple[i].value);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    printf("\tDone getting FEA materials\n");
    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(keyValue);

    return status;
}

// Get the property properties from a capsTuple
int fea_getProperty(void *aimInfo,
                    int numPropertyTuple,
                    capsTuple propertyTuple[],
                    mapAttrToIndexStruct *groupMap,
                    feaUnitsStruct *feaUnits,
                    feaProblemStruct *feaProblem) {

    /*! \page feaProperty FEA Property
     * Structure for the property tuple  = ("Property Name", "Value").
     * "Property Name" defines the reference <c>capsGroup</c> for the property being specified.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringProperty) or a single string keyword
     *	(see Section \ref keyStringProperty).
     */

    int status; //Function return

    int i, j,  matIndex; // Indexing
    int found; // Bool test function

    int pidIndex; // Property identification index found in attribute map

    int tempInteger;

    char *keyValue = NULL;
    char *keyWord = NULL;
    char *tempString = NULL;
    char **tempStringArray = NULL;

    // double tempDouble;

    // Destroy our property structures coming in if aren't 0 and NULL already
    if (feaProblem->feaProperty != NULL) {
        for (i = 0; i < feaProblem->numProperty; i++) {
            status = destroy_feaPropertyStruct(&feaProblem->feaProperty[i]);
            AIM_STATUS(aimInfo, status);
        }
    }
    AIM_FREE(feaProblem->feaProperty);
    feaProblem->numProperty = 0;

#if 0 // This check does not work for aeroelastic
    // Check that all capsGroups have been given a property
    if (numPropertyTuple != groupMap->numAttribute) {
      AIM_ERROR  (aimInfo, "Number of Property (%d) does not match number of capsGroups (%d)!", numPropertyTuple, groupMap->numAttribute);
      AIM_ADDLINE(aimInfo, "capsGroups:");
      for (i = 0; i < groupMap->numAttribute; i++)
        AIM_ADDLINE(aimInfo, "  %s", groupMap->attributeName[i]);
      AIM_ADDLINE(aimInfo, "Property names:");
      for (i = 0; i < numPropertyTuple; i++)
        AIM_ADDLINE(aimInfo, "  %s", propertyTuple[i].name);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
#endif

    printf("\nGetting FEA properties.......\n");

    feaProblem->numProperty = numPropertyTuple;
    printf("\tNumber of properties - %d\n", feaProblem->numProperty);

    if (feaProblem->numProperty > 0) {

        AIM_ALLOC(feaProblem->feaProperty, feaProblem->numProperty, feaPropertyStruct, aimInfo, status);
        for (i = 0; i < feaProblem->numProperty; i++) {
          status = initiate_feaPropertyStruct(&feaProblem->feaProperty[i]);
          AIM_STATUS(aimInfo, status);
        }
    } else {
        AIM_ERROR(aimInfo, "Number of property values in input tuple is 0\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    for (i = 0; i < feaProblem->numProperty; i++) {
        status = initiate_feaPropertyStruct(&feaProblem->feaProperty[i]);
        AIM_STATUS(aimInfo, status, "Unable to initiate feaProperty structure (number = %d)", i);
    }

    for (i = 0; i < feaProblem->numProperty; i++) {

        printf("\tProperty name - %s\n", propertyTuple[i].name);

        // Set property name from tuple name
        AIM_STRDUP(feaProblem->feaProperty[i].name, propertyTuple[i].name, aimInfo, status);

        // Get to property ID number from the attribute map
        status = get_mapAttrToIndexIndex(groupMap, feaProblem->feaProperty[i].name, &pidIndex);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "Tuple name '%s' not found in attribute map of PIDS!!!!\n", feaProblem->feaProperty[i].name);
            AIM_FREE(keyValue);
            goto cleanup;
        } else {
            feaProblem->feaProperty[i].propertyID = pidIndex;
        }

        // Do we have a json string?
        if (strncmp(propertyTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n",propertyTuple[i].value);

            /*! \page feaProperty
             * \section jsonStringProperty JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN)
             *  (e.g. "Value" = {"shearMembraneRatio": 0.83, "bendingInertiaRatio": 1.0, "membraneThickness": 0.2, "propertyType": "Shell"})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>propertyType = No Default value</B> </li> <br>
             *  Type of property to apply to a given capsGroup <c>Name</c>. Options: ConcentratedMass, Rod,
             *  Bar, Shear, Shell, Composite, and Solid
             * </ul>
             * \elseif (MASSTRAN)
             * <ul>
             *  <li> <B>propertyType = No Default value</B> </li> <br>
             *  Type of property to apply to a given capsGroup <c>Name</c>. Options: ConcentratedMass, Shell
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>propertyType = No Default value</B> </li> <br>
             *  Type of property to apply to a give capsGroup <c>Name</c>. Options: ConcentratedMass, Rod,
             *  Bar, Shear, Shell, Membrane, Composite, and Solid
             * </ul>
             * \endif
             *
             */

            // Get property Type
            keyWord = "propertyType";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                //{UnknownProperty, ConcentratedMass Rod, Bar, Shear, Shell, Composite, Solid}
                if      (strcasecmp(keyValue, "\"ConcentratedMass\"") == 0) feaProblem->feaProperty[i].propertyType = ConcentratedMass;
                else if (strcasecmp(keyValue, "\"Rod\"")              == 0) feaProblem->feaProperty[i].propertyType = Rod;
                else if (strcasecmp(keyValue, "\"Bar\"")              == 0) feaProblem->feaProperty[i].propertyType = Bar;
                else if (strcasecmp(keyValue, "\"Shear\"")            == 0) feaProblem->feaProperty[i].propertyType = Shear;
                else if (strcasecmp(keyValue, "\"Shell\"")            == 0) feaProblem->feaProperty[i].propertyType = Shell;
                else if (strcasecmp(keyValue, "\"Membrane\"")         == 0) feaProblem->feaProperty[i].propertyType = Membrane;
                else if (strcasecmp(keyValue, "\"Composite\"")        == 0) feaProblem->feaProperty[i].propertyType = Composite;
                else if (strcasecmp(keyValue, "\"Solid\"")            == 0) feaProblem->feaProperty[i].propertyType = Solid;
                else {
                    AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for Property tuple %s, current options are "
                                       "\"Rod, Bar, Shear, Shell, Composite, and Solid\"\n", keyWord,
                                       keyValue,
                                       propertyTuple[i].name);
                    AIM_FREE(keyValue);

                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

            } else {
                AIM_ERROR(aimInfo, "\tNo \"%s\" specified for Property tuple %s, this mandatory! Current options are "
                                   "\"ConcentratedMass, Rod, Bar, Shear, Shell, Composite, and Solid\"\n", keyWord,
                                   propertyTuple[i].name);
                AIM_FREE(keyValue);

                status = CAPS_NOTFOUND;
                goto cleanup;
            }
            AIM_FREE(keyValue);

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             *  <ul>
             *  <li> <B>material = "Material Name" (\ref feaMaterial) </B> </li> <br>
             *  "Material Name" from \ref feaMaterial to use for property. If no material is set the first material
             *   created will be used
             *  </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             *  <ul>
             *  <li> <B>material = `Material Name' (\ref feaMaterial) </B> </li> <br>
             *  `Material Name' from \ref feaMaterial to use for property. If no material is set the first material
             *   created will be used
             * </ul>
             * \endif
             */
            keyWord = "material";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                found = (int) false;
                for (matIndex = 0; matIndex < feaProblem->numMaterial; matIndex++ ) {

                    AIM_FREE(tempString);
                    tempString = string_removeQuotation(keyValue);
                    AIM_NOTNULL(tempString, aimInfo, status);

                    if (strcasecmp(feaProblem->feaMaterial[matIndex].name, tempString) == 0) {
                        feaProblem->feaProperty[i].materialID = feaProblem->feaMaterial[matIndex].materialID;

                        AIM_STRDUP(feaProblem->feaProperty[i].materialName, feaProblem->feaMaterial[matIndex].name, aimInfo, status);
                        feaProblem->feaProperty[i].materialName[strlen(feaProblem->feaMaterial[matIndex].name)] = '\0';

                        found = (int) true;
                        break;
                    }
                }

                AIM_FREE(tempString);

                if (found == (int) false) {
                    AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for Property tuple %s. No match in Material tuple\n", keyWord,
                                                                                                                       keyValue,
                                                                                                                       propertyTuple[i].name);
                    AIM_FREE(keyValue);
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

            } else {

                if (feaProblem->feaProperty[i].propertyType != ConcentratedMass &&
                    feaProblem->feaProperty[i].propertyType != Composite) {
                    printf("\tNo \"%s\" specified for Property tuple %s, defaulting to an index of 1\n", keyWord,
                            propertyTuple[i].name);
                }

                feaProblem->feaProperty[i].materialID = 1;
                if (feaProblem->numMaterial > 0) {
                    AIM_STRDUP(feaProblem->feaProperty[i].materialName, feaProblem->feaMaterial[0].name, aimInfo, status);
                }
            }

            AIM_FREE(keyValue);

            //  Fill up properties ///

            // Rods
            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             *  <ul>
             *  <li> <B>crossSecArea = 0.0</B> </li> <br>
             *  Cross sectional area.
             *  </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             *  <ul>
             *  <li> <B>crossSecArea = 0.0</B> </li> <br>
             *  Cross sectional area.
             *  </ul>
             * \endif
             */
            keyWord = "crossSecArea";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].crossSecArea);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             *  <ul>
             *  <li> <B>torsionalConst = 0.0</B> </li> <br>
             *  Torsional constant.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             *  <ul>
             *  <li> <B>torsionalConst = 0.0</B> </li> <br>
             *  Torsional constant.
             * </ul>
             * \endif
             */
            keyWord = "torsionalConst";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].torsionalConst);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             *  <ul>
             *  <li> <B>torsionalStressReCoeff = 0.0</B> </li> <br>
             *  Torsional stress recovery coefficient.
             *  </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>torsionalStressReCoeff = 0.0</B> </li> <br>
             *  Torsional stress recovery coefficient.
             * </ul>
             * \endif
             */
            keyWord = "torsionalStressReCoeff";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].torsionalStressReCoeff);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (HSM)
             *  <ul>
             *  <li> <B>massPerLength = 0.0</B> </li> <br>
             *  Mass per unit length.
             *  </ul>
             *
             * \elseif ( MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             * <ul>
             *  <li> <B>massPerLength = 0.0</B> </li> <br>
             *  Non-structural mass per unit length.
             * </ul>
             *
             * \endif
             */
            keyWord = "massPerLength";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].massPerLength);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            // Bar - see rod for additional variables

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>zAxisInertia = 0.0</B> </li> <br>
             *  Section moment of inertia about the element z-axis.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>zAxisInertia = 0.0</B> </li> <br>
             *  Section moment of inertia about the element z-axis.
             * </ul>
             * \endif
             */
            keyWord = "zAxisInertia";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].zAxisInertia);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             *  <ul>
             *  <li> <B>yAxisInertia = 0.0</B> </li> <br>
             *  Section moment of inertia about the element y-axis.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>yAxisInertia = 0.0</B> </li> <br>
             *  Section moment of inertia about the element y-axis.
             * </ul>
             * \endif
             */
            keyWord = "yAxisInertia";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].yAxisInertia);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>yCoords[4] = [0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Element y-coordinates, in the bar cross-section, of four points at which to recover stresses
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>yCoords[4] = [0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Element y-coordinates, in the bar cross-section, of four points at which to recover stresses
             * </ul>
             * \endif
             */
            keyWord = "yCoords";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaProperty[i].yCoords)/sizeof(double),
                                              feaProblem->feaProperty[i].yCoords);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>zCoords[4] = [0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Element z-coordinates, in the bar cross-section, of four points at which to recover stresses
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>zCoords[4] = [0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Element z-coordinates, in the bar cross-section, of four points at which to recover stresses
             * </ul>
             * \endif
             */
            keyWord = "zCoords";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaProperty[i].zCoords)/sizeof(double),
                                              feaProblem->feaProperty[i].zCoords);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>areaShearFactors[2] = [0.0, 0.0]</B> </li> <br>
             *  Area factors for shear.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>areaShearFactors[2] = [0.0, 0.0]</B> </li> <br>
             *  Area factors for shear.
             * </ul>
             * \endif
             */
            keyWord = "areaShearFactors";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaProperty[i].areaShearFactors)/sizeof(double),
                                              feaProblem->feaProperty[i].areaShearFactors);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>crossProductInertia = 0.0</B> </li> <br>
             *  Section cross-product of inertia.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>crossProductInertia = 0.0</B> </li> <br>
             *  Section cross-product of inertia.
             * </ul>
             * \endif
             */
            keyWord = "crossProductInertia";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].crossProductInertia);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN )
             * <ul>
             *  <li> <B>crossSecType = NULL</B> </li> <br>
             *  Cross-section type. Must be one of following character variables: BAR, BOX, BOX1,
             *  CHAN, CHAN1, CHAN2, CROSS, H, HAT, HEXA, I, I1, ROD, T, T1, T2, TUBE, or Z.
             * </ul>
             * \elseif ASTROS
             * <ul>
             *  <li> <B>crossSecType = NULL</B> </li> <br>
             *  Cross-section type. Must be one of following character variables: I, T, BOX, BAR,
             *  TUBE, ROD, HAT, or GBOX.
             * </ul>
             * \endif
             */
             keyWord = "crossSecType";
             status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
             if (status == CAPS_SUCCESS) {
                 feaProblem->feaProperty[i].crossSecType = string_removeQuotation(keyValue);
                 AIM_FREE(keyValue);
             }

             /*! \page feaProperty
              *
              * \if (MYSTRAN || NASTRAN || ASTROS)
              * <ul>
              *  <li> <B>crossSecDimension = [0,0,0,....]</B> </li> <br>
              *  Cross-sectional dimensions (length of array is dependent on the "crossSecType"). Max
              *  supported length array is 10!
              * </ul>
              * \endif
              */
             keyWord = "crossSecDimension";
             status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
             if (status == CAPS_SUCCESS) {
                 status = string_toDoubleArray(keyValue,
                                               (int) sizeof(feaProblem->feaProperty[i].crossSecDimension)/sizeof(double),
                                               feaProblem->feaProperty[i].crossSecDimension);
                 AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                 AIM_FREE(keyValue);
             }


            /* UNDOCUMENTED orientation vector of bars - MAY be changed in the future without warning
             */
            keyWord = "orientationVec";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaProperty[i].orientationVec)/sizeof(double),
                                              feaProblem->feaProperty[i].orientationVec);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }


            // Shear

            // Shell

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN || MASSTRAN)
             * <ul>
             *  <li> <B>membraneThickness = 0.0</B> </li> <br>
             *  Membrane thickness.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>membraneThickness = 0.0</B> </li> <br>
             *  Membrane thickness.
             * </ul>
             * \endif
             */
            keyWord = "membraneThickness";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
              if ( feaUnits->length != NULL ) {
                  status = string_toDoubleUnits(aimInfo, keyValue, feaUnits->length, &feaProblem->feaProperty[i].membraneThickness);
              } else {
                  status = string_toDouble(keyValue, &feaProblem->feaProperty[i].membraneThickness);
              }
              AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
              AIM_FREE(keyValue);
            }


            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>bendingInertiaRatio = 1.0</B> </li> <br>
             *  Ratio of actual bending moment inertia to the bending inertia of a solid plate of thickness "membraneThickness"
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>bendingInertiaRatio = 1.0</B> </li> <br>
             *  Ratio of actual bending moment inertia to the bending inertia of a solid plate of thickness "membraneThickness"
             * </ul>
             * \endif
             */
            keyWord = "bendingInertiaRatio";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].bendingInertiaRatio);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }


            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>shearMembraneRatio = 5.0/6.0</B> </li> <br>
             *  Ratio shear thickness to membrane thickness.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>shearMembraneRatio = 5.0/6.0</B> </li> <br>
             *  Ratio shear thickness to membrane thickness.
             * </ul>
             * \endif
             */
            keyWord = "shearMembraneRatio";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].shearMembraneRatio);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>materialBending = "Material Name" (\ref feaMaterial)</B> </li> <br>
             *  "Material Name" from \ref feaMaterial to use for property bending. If no material is given and
             *  "bendingInertiaRatio" is greater than 0, the material name provided in "material" is used.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>materialBending = "Material Name" (\ref feaMaterial)</B> </li> <br>
             *  "Material Name" from \ref feaMaterial to use for property bending.
             * </ul>
             * \endif
             */
            keyWord = "materialBending"; // Shell specific materials
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                found = (int) false;
                for (matIndex = 0; matIndex < feaProblem->numMaterial; matIndex++ ) {

                    if (tempString != NULL) EG_free(tempString);
                    tempString = string_removeQuotation(keyValue);
                    AIM_NOTNULL(tempString, aimInfo, status);

                    if (strcasecmp(feaProblem->feaMaterial[matIndex].name, tempString) == 0) {
                        feaProblem->feaProperty[i].materialBendingID = feaProblem->feaMaterial[matIndex].materialID;
                        found = (int) true;

                        AIM_FREE(keyValue);
                        break;
                    }

                }

                AIM_FREE(tempString);

                if (found == (int) false) {
                    AIM_ERROR(aimInfo, "\tUnrecognized \"%s\" specified (%s) for Property tuple %s. No match in Material tuple\n", keyWord,
                              keyValue,
                              propertyTuple[i].name);
                    AIM_FREE(keyValue);
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

            } else { // Don't default to anything - yet

                /*
                printf("No \"%s\" specified for Property tuple %s which is a shell element, "
                       " no bending material will be specified\n", keyWord,
                                                                   propertyTuple[i].name);

                feaProblem->feaProperty[i].materialBendingID = 0;
                 */

                if (feaProblem->feaProperty[i].bendingInertiaRatio > 0) {

                    feaProblem->feaProperty[i].materialBendingID = feaProblem->feaProperty[i].materialID;
                }
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>materialShear = "Material Name" (\ref feaMaterial)</B> </li> <br>
             *  "Material Name" from \ref feaMaterial to use for property shear. If no material is given and
             *  "shearMembraneRatio" is greater than 0, the material name provided in "material" is used.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>materialShear = "Material Name" (\ref feaMaterial)</B> </li> <br>
             *  "Material Name" from \ref feaMaterial to use for property shear.
             * </ul>
             * \endif
             */
            keyWord = "materialShear"; // Shell specific materials
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                found = (int) false;
                for (matIndex = 0; matIndex < feaProblem->numMaterial; matIndex++ ) {

                    AIM_FREE(tempString);
                    tempString = string_removeQuotation(keyValue);
                    AIM_NOTNULL(tempString, aimInfo, status);

                    if (strcasecmp(feaProblem->feaMaterial[matIndex].name, tempString) == 0) {
                        feaProblem->feaProperty[i].materialShearID = feaProblem->feaMaterial[matIndex].materialID;
                        found = (int) true;

                        AIM_FREE(keyValue);
                        break;
                    }

                }

                AIM_FREE(tempString);

                if (found == (int) false) {
                    AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for Property tuple %s. No match in Material tuple\n", keyWord,
                              keyValue,
                              propertyTuple[i].name);
                    AIM_FREE(keyValue);
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

            } else {
                /*
                printf("No \"%s\" specified for Property tuple %s which is a shell element, "
                       " no shear material will be specified\n", keyWord,
                       propertyTuple[i].name);

                feaProblem->feaProperty[i].materialShearID = 0;
                 */

                if (feaProblem->feaProperty[i].shearMembraneRatio > 0) {
                    feaProblem->feaProperty[i].materialShearID = feaProblem->feaProperty[i].materialID;
                }
            }

            /*! \page feaProperty
             *
             * \if (HSM || MASSTRAN)
             * <ul>
             *  <li> <B>massPerArea = 0.0</B> </li> <br>
             *  Mass per unit area.
             * </ul>
             * \elseif (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             *
             * <ul>
             *  <li> <B>massPerArea = 0.0</B> </li> <br>
             *  Non-structural mass per unit area.
             * </ul>
             * \endif
             */
            keyWord = "massPerArea";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                if ( feaUnits->densityArea != NULL ) {
                    status = string_toDoubleUnits(aimInfo, keyValue, feaUnits->densityArea, &feaProblem->feaProperty[i].massPerArea);
                    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                } else {
                    status = string_toDouble(keyValue, &feaProblem->feaProperty[i].massPerArea);
                    AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                }
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>zOffsetRel = 0.0</B> </li> <br>
             *  Relative offset from the surface of grid points to the element reference
             *  plane as a percentage of the thickness. zOffSet = thickness*zOffsetRel/100
             *
             * </ul>
             * \endif
             */
            keyWord = "zOffsetRel";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].zOffsetRel);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            // Composite

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>compositeMaterial = "no default" </B> </li> <br>
             *  List of "Material Name"s, ["Material Name -1", "Material Name -2", ...], from \ref feaMaterial to use for composites.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>compositeMaterial = "no default" </B> </li> <br>
             *  List of "Material Name"s, ["Material Name -1", "Material Name -2", ...], from \ref feaMaterial to use for composites.
             * </ul>
             * \endif
             */
            keyWord = "compositeMaterial";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                status = string_toStringDynamicArray(keyValue, &feaProblem->feaProperty[i].numPly, &tempStringArray);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);

                AIM_ALLOC(feaProblem->feaProperty[i].compositeMaterialID, feaProblem->feaProperty[i].numPly, int, aimInfo, status);

                for (j = 0; j < feaProblem->feaProperty[i].numPly; j++) {
                    found = (int) false;
                    for (matIndex = 0; matIndex < feaProblem->numMaterial; matIndex++ ) {

                        AIM_FREE(tempString);
                        tempString = string_removeQuotation(tempStringArray[j]);
                        AIM_NOTNULL(tempString, aimInfo, status);

                        if (strcasecmp(feaProblem->feaMaterial[matIndex].name, tempString) == 0) {
                            feaProblem->feaProperty[i].compositeMaterialID[j] = feaProblem->feaMaterial[matIndex].materialID;
                            found = (int) true;
                            break;
                        }

                    }

                    if (found == (int) false) {
                        AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for Property tuple %s. No match in Material tuple\n", keyValue,
                                           keyWord,
                                           propertyTuple[i].name);

                        AIM_FREE(keyValue);
                        AIM_FREE(tempString);

                        (void) string_freeArray(feaProblem->feaProperty[i].numPly, &tempStringArray);

                        status = CAPS_NOTFOUND;
                        goto cleanup;
                    }

                    AIM_FREE(tempString);
                }
                AIM_FREE(keyValue);
                status = string_freeArray(feaProblem->feaProperty[i].numPly, &tempStringArray);
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>shearBondAllowable = 0.0 </B> </li> <br>
             *  Allowable interlaminar shear stress.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>shearBondAllowable = 0.0 </B> </li> <br>
             *  Allowable interlaminar shear stress.
             * </ul>
             * \endif
             */
            keyWord = "shearBondAllowable";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].compositeShearBondAllowable);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>symmetricLaminate = False </B> </li> <br>
             *  Symmetric lamination option. True- SYM only half the plies are specified, for odd number plies 1/2 thickness
             *   of center ply is specified with the first ply being the bottom ply in the stack, default (False) all plies specified.
             * </ul>
             * \elseif (ASTROS)
             * <ul>
             *  <li> <B>symmetricLaminate = False </B> </li> <br>
             *  Symmetric lamination option. If "True" only half the plies are specified (the plies will be repeated in
             *  reverse order internally in the PCOMP card). For an odd number of plies, the 1/2 thickness
             *   of the center ply is specified with the first ply being the bottom ply in the stack, default (False) all plies specified.
             * </ul>
             * \endif
             */
            keyWord = "symmetricLaminate";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toBoolean(keyValue, &feaProblem->feaProperty[i].compositeSymmetricLaminate);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>compositeFailureTheory = "(no default)" </B> </li> <br>
             *  Composite failure theory. Options: "HILL", "HOFF", "TSAI", and "STRN"
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>compositeFailureTheory = "(no default)" </B> </li> <br>
             *  Composite failure theory.
             * </ul>
             * \endif
             */
            keyWord = "compositeFailureTheory";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                feaProblem->feaProperty[i].compositeFailureTheory = string_removeQuotation(keyValue);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>compositeThickness = (no default) </B> </li> <br>
             *  List of composite thickness for each layer (e.g. [1.2, 4.0, 3.0]). If the length of this list doesn't match the length
             *  of the "compositeMaterial" list, the list is either truncated [ >length("compositeMaterial")] or expanded [ <length("compositeMaterial")]
             *  in which case the last thickness provided is repeated.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>compositeThickness = (no default) </B> </li> <br>
             *  List of composite thickness for each layer (e.g. [1.2, 4.0, 3.0]). If the length of this list doesn't match the length
             *  of the "compositeMaterial" list, the list is either truncated [ >length("compositeMaterial")] or expanded [ <length("compositeMaterial")]
             *  in which case the last thickness provided is repeated.
             * </ul>
             * \endif
             */
            keyWord = "compositeThickness";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleDynamicArray(keyValue, &tempInteger,
                                                     &feaProblem->feaProperty[i].compositeThickness);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);

                if (tempInteger < feaProblem->feaProperty[i].numPly) {

                    printf("\tThe number of thicknesses provided does not match the number of materials for the composite. "
                            "The last thickness will be repeated %d times\n", feaProblem->feaProperty[i].numPly - tempInteger);

                    AIM_REALL(feaProblem->feaProperty[i].compositeThickness,
                              feaProblem->feaProperty[i].numPly, double, aimInfo, status);

                    for (j = 0; j < feaProblem->feaProperty[i].numPly - tempInteger; j++) {

                        feaProblem->feaProperty[i].compositeThickness[j+tempInteger] = feaProblem->feaProperty[i].compositeThickness[tempInteger-1];
                    }
                }

                if (tempInteger > feaProblem->feaProperty[i].numPly) {

                    printf("\tThe number of thicknesses provided does not match the number of materials for the composite. "
                            "The last %d thicknesses will be not be used\n", tempInteger -feaProblem->feaProperty[i].numPly);

                    AIM_REALL(feaProblem->feaProperty[i].compositeThickness,
                              feaProblem->feaProperty[i].numPly, double, aimInfo, status);
                }
            } else {

                if (feaProblem->feaProperty[i].numPly != 0 &&
                        feaProblem->feaProperty[i].propertyType == Composite) {

                    AIM_ERROR(aimInfo, "\"compositeMaterial\" have been set but no thicknesses (\"compositeThickness\") provided!!!");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>compositeOrientation = (no default) </B> </li> <br>
             *  List of composite orientations (angle relative element material axis) for each layer (eg. [5.0, 10.0, 30.0]).
             *  If the length of this list doesn't match the length of the "compositeMaterial" list, the list is either
             *  truncated [ >length("compositeMaterial")] or expanded [ <length("compositeMaterial")] in which case
             *  the last orientation provided is repeated.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>compositeOrientation = (no default) </B> </li> <br>
             *  List of composite orientations (angle relative element material axis) for each layer (eg. [5.0, 10.0, 30.0]).
             *  If the length of this list doesn't match the length of the "compositeMaterial" list, the list is either
             *  truncated [ >length("compositeMaterial")] or expanded [ <length("compositeMaterial")] in which case
             *  the last orientation provided is repeated.
             * </ul>
             * \endif
             */
            keyWord = "compositeOrientation";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleDynamicArray(keyValue, &tempInteger,
                                                     &feaProblem->feaProperty[i].compositeOrientation);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);

                if (tempInteger < feaProblem->feaProperty[i].numPly) {

                    printf("\tThe number of orientations provided does not match the number of materials for the composite. "
                            "The last orientation will be repeated %d times\n", feaProblem->feaProperty[i].numPly - tempInteger);

                    AIM_REALL(feaProblem->feaProperty[i].compositeOrientation,
                              feaProblem->feaProperty[i].numPly, double, aimInfo, status);

                    for (j = 0; j < feaProblem->feaProperty[i].numPly - tempInteger; j++) {

                        feaProblem->feaProperty[i].compositeOrientation[j+tempInteger] = feaProblem->feaProperty[i].compositeOrientation[tempInteger-1];
                    }
                }

                if (tempInteger > feaProblem->feaProperty[i].numPly) {

                    printf("\tThe number of orientations provided does not match the number of materials for the composite. "
                            "The last %d orientation will be not be used\n", tempInteger -feaProblem->feaProperty[i].numPly);

                    AIM_REALL(feaProblem->feaProperty[i].compositeOrientation,
                              feaProblem->feaProperty[i].numPly, double, aimInfo, status);
                }

            } else {

                if (feaProblem->feaProperty[i].numPly != 0 &&
                        feaProblem->feaProperty[i].propertyType == Composite) {

                    AIM_ERROR(aimInfo, "\"compositeMaterial\" have been set but no Orientation  (\"compositeOrientation\") provided!!!");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            // Mass

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>mass = 0.0</B> </li> <br>
             *  Mass value.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>mass = 0.0</B> </li> <br>
             *  Mass value.
             * </ul>
             * \endif
             */
            keyWord = "mass";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
              if ( feaUnits->mass != NULL ) {
                status = string_toDoubleUnits(aimInfo, keyValue, feaUnits->mass, &feaProblem->feaProperty[i].mass);
              } else {
                status = string_toDouble(keyValue, &feaProblem->feaProperty[i].mass);
              }
              AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
              AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN || MASSTRAN)
             * <ul>
             *  <li> <B>massOffset = [0.0, 0.0, 0.0]</B> </li> <br>
             *  Offset distance from the grid point to the center of gravity for a concentrated mass.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>massOffset = [0.0, 0.0, 0.0]</B> </li> <br>
             *  Offset distance from the grid point to the center of gravity for a concentrated mass.
             * </ul>
             * \endif
             */
            keyWord = "massOffset";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaProperty[i].massOffset)/sizeof(double),
                                              feaProblem->feaProperty[i].massOffset);
                AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
                AIM_FREE(keyValue);
            }

            /*! \page feaProperty
             *
             * \if (MYSTRAN || NASTRAN || MASSTRAN)
             * <ul>
             *  <li> <B>massInertia = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Mass moment of inertia measured at the mass center of gravity.
             * </ul>
             * \elseif ABAQUS
             *v
             * Something else ....
             *
             * \elseif ASTROS
             * <ul>
             *  <li> <B>massInertia = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]</B> </li> <br>
             *  Mass moment of inertia measured at the mass center of gravity.
             * </ul>
             * \endif
             */
            keyWord = "massInertia";
            status = search_jsonDictionary( propertyTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
              if ( feaUnits->momentOfInertia != NULL ) {
                status = string_toDoubleArrayUnits(aimInfo, keyValue, feaUnits->momentOfInertia,
                                                   (int) sizeof(feaProblem->feaProperty[i].massInertia)/sizeof(double),
                                                   feaProblem->feaProperty[i].massInertia);
              } else {
                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaProperty[i].massInertia)/sizeof(double),
                                              feaProblem->feaProperty[i].massInertia);
              }
              AIM_STATUS(aimInfo, status, "While parsing \"%s\":\"%s\"", keyWord, keyValue);
              AIM_FREE(keyValue);
            }

        } else {

            /*! \page feaProperty
             * \section keyStringProperty Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined property lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */

            // CALL property look up
            AIM_ERROR(aimInfo, "Property tuple value ('s') is expected to be a JSON string", propertyTuple[i].value);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    printf("\tDone getting FEA properties\n");
    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(tempString);
    AIM_FREE(keyValue);

    return status;
}

// Get the constraint properties from a capsTuple
int fea_getConstraint(void *aimInfo,
                      int numConstraintTuple,
                      capsTuple constraintTuple[],
                      mapAttrToIndexStruct *attrMap,
                      feaProblemStruct *feaProblem) {

    /*! \page feaConstraint FEA Constraint
     * Structure for the constraint tuple  = ("Constraint Name", "Value").
     * "Constraint Name" defines the reference name for the constraint being specified.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringConstraint) or a single string keyword
     *	(see Section \ref keyStringConstraint).
     */

    int status; //Function return

    int i, groupIndex, attrIndex, nodeIndex; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    int numGroupName = 0;
    char **groupName = NULL;

    feaMeshDataStruct *feaData;

    // Destroy our constraint structures coming in if aren't 0 and NULL already
    if (feaProblem->feaConstraint != NULL) {
        for (i = 0; i < feaProblem->numConstraint; i++) {
            status = destroy_feaConstraintStruct(&feaProblem->feaConstraint[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (feaProblem->feaConstraint != NULL) EG_free(feaProblem->feaConstraint);
    feaProblem->feaConstraint = NULL;
    feaProblem->numConstraint = 0;

    printf("\nGetting FEA constraints.......\n");

    feaProblem->numConstraint = numConstraintTuple;

    printf("\tNumber of constraints - %d\n", feaProblem->numConstraint);

    // Allocate constraints
    if (feaProblem->numConstraint > 0) {
        feaProblem->feaConstraint = (feaConstraintStruct *) EG_alloc(feaProblem->numConstraint * sizeof(feaConstraintStruct));

        if (feaProblem->feaConstraint == NULL ) return EGADS_MALLOC;
    }

    // Initiate constraints to default values
    for (i = 0; i < feaProblem->numConstraint; i++) {
        status = initiate_feaConstraintStruct(&feaProblem->feaConstraint[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    // Loop through tuples and fill out the constraint structures
    for (i = 0; i < feaProblem->numConstraint; i++) {

        printf("\tConstraint name - %s\n", constraintTuple[i].name );

        // Set constraint name to tuple attribute name
        feaProblem->feaConstraint[i].name = (char *) EG_alloc((strlen(constraintTuple[i].name) + 1)*sizeof(char));
        if (feaProblem->feaConstraint[i].name == NULL) return EGADS_MALLOC;

        memcpy(feaProblem->feaConstraint[i].name, constraintTuple[i].name, strlen(constraintTuple[i].name)*sizeof(char));
        feaProblem->feaConstraint[i].name[strlen(constraintTuple[i].name)] = '\0';

        // Set constraint id -> 1 bias
        feaProblem->feaConstraint[i].constraintID = i+1;

        // Do we have a json string?
        if (strncmp(constraintTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", constraintTuple[i].value);

            /*! \page feaConstraint
             * \section jsonStringConstraint JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS)
             *  (eg. "Value" = {"groupName": "plateEdge", "dofConstraint": 123456})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS)
             * <ul>
             *  <li> <B>constraintType = "ZeroDisplacement"</B> </li> <br>
             *  Type of constraint. Options: "Displacement", "ZeroDisplacement".
             * </ul>
             * \endif
             *
             */
            // Get constraint Type
            keyWord = "constraintType";
            status = search_jsonDictionary( constraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                //{UnknownConstraint, Displacement, ZeroDisplacement}
                if      (strcasecmp(keyValue, "\"Displacement\"")      == 0) feaProblem->feaConstraint[i].constraintType = Displacement;
                else if (strcasecmp(keyValue, "\"ZeroDisplacement\"")  == 0) feaProblem->feaConstraint[i].constraintType = ZeroDisplacement;
                else {

                    printf("\tUnrecognized \"%s\" specified (%s) for Constraint tuple %s, defaulting to \"ZeroDisplacement\"\n", keyWord,
                                                                                                                                 keyValue,
                                                                                                                                 constraintTuple[i].name);
                    feaProblem->feaConstraint[i].constraintType = ZeroDisplacement;
                }

            } else {

                printf("\tNo \"%s\" specified for Constraint tuple %s, defaulting to \"ZeroDisplacement\"\n", keyWord,
                        constraintTuple[i].name);
                feaProblem->feaConstraint[i].constraintType = ZeroDisplacement;
            }

            if (keyValue != NULL) {
                EG_free(keyValue);
                keyValue = NULL;
            }

            // Get constraint node set

            /*! \page feaConstraint
             *
             * \if MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsConstraint</c> names on which to apply the constraint
             *  (e.g. "Name1" or ["Name1","Name2",...]. If not provided, the constraint tuple name will be
             *  used.
             * </ul>
             * \endif
             *
             */
            keyWord = "groupName";
            status = search_jsonDictionary( constraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;

                if (status != CAPS_SUCCESS) return status;

            } else {

                printf("\tNo \"%s\" specified for Constraint tuple %s, going to use constraint name\n", keyWord,
                                                                                                        constraintTuple[i].name);

                status = string_toStringDynamicArray(constraintTuple[i].name, &numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;

            }

            // Determine how many point constraints we have
            feaProblem->feaConstraint[i].numGridID = 0;
            for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                status = get_mapAttrToIndexIndex(attrMap, (const char *) groupName[groupIndex], &attrIndex);

                if (status == CAPS_NOTFOUND) {
                    printf("\tName %s not found in attribute map of capsConstraints!!!!\n", groupName[groupIndex]);
                    continue;
                } else if (status != CAPS_SUCCESS) return status;

                // Now lets loop through the grid to see how many grid points have the attrIndex
                for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                    if (feaProblem->feaMesh.node[nodeIndex].analysisType == MeshStructure) {
                        feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;
                    } else {
                        continue;
                    }

                    if (feaData->constraintIndex == attrIndex) {

                        feaProblem->feaConstraint[i].numGridID += 1;

                        // Allocate/Re-allocate grid ID array
                        if (feaProblem->feaConstraint[i].numGridID == 1) {
                            feaProblem->feaConstraint[i].gridIDSet = (int *) EG_alloc(feaProblem->feaConstraint[i].numGridID*sizeof(int));
                        } else {
                            feaProblem->feaConstraint[i].gridIDSet = (int *) EG_reall(feaProblem->feaConstraint[i].gridIDSet,
                                                                                      feaProblem->feaConstraint[i].numGridID*sizeof(int));
                        }

                        if (feaProblem->feaConstraint[i].gridIDSet == NULL) {
                            status = string_freeArray(numGroupName, &groupName);
                            if (status != CAPS_SUCCESS) printf("Status %d during string_freeArray\n", status);

                            return EGADS_MALLOC;
                        }

                        // Set grid ID value -> 1 bias
                        feaProblem->feaConstraint[i].gridIDSet[feaProblem->feaConstraint[i].numGridID-1] = feaProblem->feaMesh.node[nodeIndex].nodeID;
                    }
                }
            }

            status = string_freeArray(numGroupName, &groupName);
            if (status != CAPS_SUCCESS) return status;
            groupName = NULL;


            //Fill up constraint properties
            /*! \page feaConstraint
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS)
             * <ul>
             *  <li> <B>dofConstraint = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom that will be constrained (123 - zero translation in all three
             *  directions).
             * </ul>
             * \endif
             */
            keyWord = "dofConstraint";
            status = search_jsonDictionary( constraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaConstraint[i].dofConstraint);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConstraint
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || HSM || ABAQUS)
             * <ul>
             *  <li> <B>gridDisplacement = 0.0 </B> </li> <br>
             *  Value of displacement for components defined in "dofConstraint".
             * </ul>
             * \endif
             */
            keyWord = "gridDisplacement";
            status = search_jsonDictionary( constraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaConstraint[i].gridDisplacement);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

        } else {

            /*! \page feaConstraint
             * \section keyStringConstraint Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined constraint lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */

            // Call some look up table maybe?
            printf("\tError: Constraint tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;
        }
    }

    if (keyValue != NULL) EG_free(keyValue);

    printf("\tDone getting FEA constraints\n");

    status = CAPS_SUCCESS;
cleanup:
    return status;
}

// Get the support properties from a capsTuple
int fea_getSupport(int numSupportTuple,
                   capsTuple supportTuple[],
                   mapAttrToIndexStruct *attrMap,
                   feaProblemStruct *feaProblem) {

    /*! \page feaSupport FEA Support
     * Structure for the support tuple  = ("Support Name", "Value").
     * "Support Name" defines the reference name for the support being specified.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringSupport) or a single string keyword
     *	(see Section \ref keyStringSupport).
     */

    int status; //Function return

    int i, groupIndex, attrIndex, nodeIndex; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    int numGroupName = 0;
    char **groupName = NULL;

    feaMeshDataStruct *feaData;

    // Destroy our support structures coming in if aren't 0 and NULL already
    if (feaProblem->feaSupport != NULL) {
        for (i = 0; i < feaProblem->numSupport; i++) {
            status = destroy_feaSupportStruct(&feaProblem->feaSupport[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    AIM_FREE(feaProblem->feaSupport);
    feaProblem->numSupport = 0;

    printf("\nGetting FEA supports.......\n");

    feaProblem->numSupport = numSupportTuple;

    printf("\tNumber of supports - %d\n", feaProblem->numSupport);

    // Allocate supports
    if (feaProblem->numSupport > 0) {
        feaProblem->feaSupport = (feaSupportStruct *) EG_alloc(feaProblem->numSupport * sizeof(feaSupportStruct));

        if (feaProblem->feaSupport == NULL ) return EGADS_MALLOC;
    }

    // Initiate supports to default values
    for (i = 0; i < feaProblem->numSupport; i++) {
        status = initiate_feaSupportStruct(&feaProblem->feaSupport[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    // Loop through tuples and fill out the support structures
    for (i = 0; i < feaProblem->numSupport; i++) {

        printf("\tSupport name - %s\n", supportTuple[i].name );

        // Set support name to tuple attribute name
        feaProblem->feaSupport[i].name = (char *) EG_alloc((strlen(supportTuple[i].name) + 1)*sizeof(char));
        if (feaProblem->feaSupport[i].name == NULL) return EGADS_MALLOC;

        memcpy(feaProblem->feaSupport[i].name, supportTuple[i].name, strlen(supportTuple[i].name)*sizeof(char));
        feaProblem->feaSupport[i].name[strlen(supportTuple[i].name)] = '\0';

        // Set support id -> 1 bias
        feaProblem->feaSupport[i].supportID = i+1;

        // Do we have a json string?
        if (strncmp(supportTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", supportTuple[i].value);

            /*! \page feaSupport
             * \section jsonStringSupport JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  (eg. "Value" = {"groupName": "plateEdge", "dofSupport": 123456})
             * \endif
             *  the following keywords ( = default values) may be used:
             */

            // Get support node set

            /*! \page feaSupport
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsConstraint</c> names on which to apply the support
             *  (e.g. "Name1" or ["Name1","Name2",...]. If not provided, the constraint tuple name will be
             *  used.
             * </ul>
             * \endif
             */
            keyWord = "groupName";
            status = search_jsonDictionary( supportTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;

                if (status != CAPS_SUCCESS) return status;

            } else {

                printf("\tNo \"%s\" specified for Support tuple %s, going to use support name\n", keyWord,
                        supportTuple[i].name);

                status = string_toStringDynamicArray(supportTuple[i].name, &numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;

            }

            // Determine how many point supports we have
            feaProblem->feaSupport[i].numGridID = 0;
            for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                status = get_mapAttrToIndexIndex(attrMap, (const char *) groupName[groupIndex], &attrIndex);

                if (status == CAPS_NOTFOUND) {
                    printf("\tName %s not found in attribute map of capsConstraints!!!!\n", groupName[groupIndex]);
                    continue;
                } else if (status != CAPS_SUCCESS) return status;

                // Now lets loop through the grid to see how many grid points have the attrIndex
                for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                    if (feaProblem->feaMesh.node[nodeIndex].analysisType == MeshStructure) {
                        feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;
                    } else {
                        continue;
                    }

                    if (feaData->constraintIndex == attrIndex) {

                        feaProblem->feaSupport[i].numGridID += 1;

                        // Allocate/Re-allocate grid ID array
                        if (feaProblem->feaSupport[i].numGridID == 1) {
                            feaProblem->feaSupport[i].gridIDSet = (int *) EG_alloc(feaProblem->feaSupport[i].numGridID*sizeof(int));
                        } else {
                            feaProblem->feaSupport[i].gridIDSet = (int *) EG_reall(feaProblem->feaSupport[i].gridIDSet,
                                                                                   feaProblem->feaSupport[i].numGridID*sizeof(int));
                        }

                        if (feaProblem->feaSupport[i].gridIDSet == NULL) {
                            status = string_freeArray(numGroupName, &groupName);
                            if (status != CAPS_SUCCESS) printf("Status %d during string_freeArray\n", status);

                            return EGADS_MALLOC;
                        }

                        // Set grid ID value -> 1 bias
                        feaProblem->feaSupport[i].gridIDSet[feaProblem->feaSupport[i].numGridID-1] = feaProblem->feaMesh.node[nodeIndex].nodeID;
                    }
                }
            }

            status = string_freeArray(numGroupName, &groupName);
            if (status != CAPS_SUCCESS) return status;
            groupName = NULL;

            //Fill up support properties
            /*! \page feaSupport
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>dofSupport = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom that will be supported (123 - zero translation in all three
             *  directions).
             * </ul>
             * \endif
             */
            keyWord = "dofSupport";
            status = search_jsonDictionary( supportTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaSupport[i].dofSupport);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

        } else {

            /*! \page feaSupport
             * \section keyStringSupport Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined support lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */

            // Call some look up table maybe?
            printf("\tError: Support tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;

        }
    }

    if (keyValue != NULL) {
        EG_free(keyValue);
        keyValue = NULL;
    }

    printf("\tDone getting FEA supports\n");
    return CAPS_SUCCESS;
}


static int fea_setConnection(char *connectionName,
                             feaConnectionTypeEnum connectionType,
                             int connectionID,
                             int elementOffSet,
                             int dofDependent,
                             double stiffnessConst,
                             double dampingConst,
                             double stressCoeff,
                             int componentNumberStart,
                             int componentNumberEnd,
                             int srcNodeID,
                             double masterWeight,
                             int masterComponent,
                             int numNode, int *node,
                             int *numConnect, feaConnectionStruct *feaConnect[]) { // Out

    int status;
    int i;

    if (numNode == 0) {
        status = CAPS_BADVALUE; // Nothing to do
        goto cleanup;
    }

    if (connectionType == RigidBodyInterpolate) {
        *numConnect += 1;

        (*feaConnect) = (feaConnectionStruct *) EG_reall((*feaConnect),*numConnect*sizeof(feaConnectionStruct));
        if ((*feaConnect) == NULL) {
            *numConnect = 0;
            status = EGADS_MALLOC;
            goto cleanup;
        }

        status = initiate_feaConnectionStruct(&(*feaConnect)[*numConnect-1]);
        if (status != CAPS_SUCCESS) return status;

        (*feaConnect)[*numConnect-1].connectionID = connectionID; // ConnectionTuple index
        (*feaConnect)[*numConnect-1].connectionType = connectionType;

        (*feaConnect)[*numConnect-1].elementID = *numConnect + elementOffSet;

        (*feaConnect)[*numConnect-1].dofDependent = dofDependent;

        (*feaConnect)[*numConnect-1].connectivity[1] = srcNodeID; // Dependent
        (*feaConnect)[*numConnect-1].numMaster = numNode; // Independent

        (*feaConnect)[*numConnect-1].masterIDSet = (int *) EG_alloc(numNode*sizeof(int)); // [numMaster]
        if ((*feaConnect)[*numConnect-1].masterIDSet == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }
        (*feaConnect)[*numConnect-1].masterWeighting = (double *) EG_alloc(numNode*sizeof(double));; // [numMaster]
        if ((*feaConnect)[*numConnect-1].masterWeighting == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        (*feaConnect)[*numConnect-1].masterComponent =(int *) EG_alloc(numNode*sizeof(int));; // [numMaster]
        if ((*feaConnect)[*numConnect-1].masterComponent == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        //printf("\tMasters (%d)", numNode);
        // Master values;
        for (i = 0; i < numNode; i++) {
          //  printf(" %d", node[i]);

            (*feaConnect)[*numConnect-1].masterIDSet[i] = node[i];
            (*feaConnect)[*numConnect-1].masterWeighting[i] = masterWeight;
            (*feaConnect)[*numConnect-1].masterComponent[i] = masterComponent;
        }
        //printf("\n");

    } else { // For all other types of connections create an individual connection - single pair of nodes

        (*feaConnect) = (feaConnectionStruct *) EG_reall((*feaConnect),(*numConnect+numNode)*sizeof(feaConnectionStruct));
        if ((*feaConnect) == NULL) {
            *numConnect = 0;
            status = EGADS_MALLOC;
            goto cleanup;
        }

        for (i = 0; i < numNode; i++) {
            status = initiate_feaConnectionStruct(&(*feaConnect)[*numConnect+i]);
            if (status != CAPS_SUCCESS) return status;
        }

        for (i = 0; i < numNode; i++) {

            *numConnect += 1;

            (*feaConnect)[*numConnect-1].name = (char *) EG_alloc((strlen(connectionName) + 1)*sizeof(char));
            if ((*feaConnect)[*numConnect-1].name == NULL) return EGADS_MALLOC;

            memcpy((*feaConnect)[*numConnect-1].name,
                   connectionName,
                   strlen(connectionName)*sizeof(char));
            (*feaConnect)[*numConnect-1].name[strlen(connectionName)] = '\0';

            (*feaConnect)[*numConnect-1].connectionID = connectionID; // ConnectionTuple index

            (*feaConnect)[*numConnect-1].connectionType = connectionType;
            (*feaConnect)[*numConnect-1].elementID = *numConnect + elementOffSet;

            (*feaConnect)[*numConnect-1].connectivity[0] = srcNodeID;
            (*feaConnect)[*numConnect-1].connectivity[1] = node[i];

            (*feaConnect)[*numConnect-1].dofDependent = dofDependent;
            (*feaConnect)[*numConnect-1].stiffnessConst = stiffnessConst;
            (*feaConnect)[*numConnect-1].dampingConst = dampingConst;
            (*feaConnect)[*numConnect-1].stressCoeff = stressCoeff;
            (*feaConnect)[*numConnect-1].componentNumberStart = componentNumberStart;
            (*feaConnect)[*numConnect-1].componentNumberEnd = componentNumberEnd;
        }
    }

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in fea_setConnection, status = %d\n", status);
        return status;
}
// Get the Connections properties from a capsTuple and create connections based on the mesh
int fea_getConnection(void *aimInfo,
                      int numConnectionTuple,
                      capsTuple connectionTuple[],
                      mapAttrToIndexStruct *attrMap,
                      feaProblemStruct *feaProblem) {

    /*! \page feaConnection FEA Connection
     * Structure for the connection tuple  = ("Connection Name", "Value").
     * "Connection Name" defines the reference name to the capsConnect being specified and denotes the "source" node
     * for the connection.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringConnection) or a single string keyword
     *	(see Section \ref keyStringConnection).
     */

    int status; //Function return

    int i, groupIndex, nodeIndex, nodeIndexDest; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    int numGroupName = 0;
    char **groupName = NULL;
    int attrIndex, attrIndexDest;

    // Values to set
    feaConnectionTypeEnum connectionType;
    int dofDependent = 0, componentNumberStart = 0, componentNumberEnd = 0;
    double stiffnessConst = 0.0, dampingConst  = 0.0, stressCoeff = 0.0, mass = 0.0;

    double weighting=1;

    int glue = (int) false, glueNumMaster = 5;
    double glueSearchRadius=0;

    feaMeshDataStruct *feaData, *feaDataDest;

    int numDestNode=0, *destNode=NULL;

    // Destroy our support structures coming in if aren't 0 and NULL already
    if (feaProblem->feaConnect != NULL) {
        for (i = 0; i < feaProblem->numConnect; i++) {
            status = destroy_feaConnectionStruct(&feaProblem->feaConnect[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (feaProblem->feaConnect != NULL) EG_free(feaProblem->feaConnect);
    feaProblem->feaConnect = NULL;
    feaProblem->numConnect = 0;

    printf("\nGetting FEA connections.......\n");

    printf("\tNumber of connection tuples - %d\n", numConnectionTuple);

    // Loop through tuples and fill out the support structures
    for (i = 0; i < numConnectionTuple; i++) {

        // Reset defaults
        dofDependent  = componentNumberStart = componentNumberEnd = 0;
        stiffnessConst = dampingConst = stressCoeff = mass = 0.0;

        weighting=1;

        glue = (int) false;
        glueNumMaster = 5;
        glueSearchRadius = 0;

        printf("\tConnection name - %s\n", connectionTuple[i].name );

        // Look for connection name in connection map
        status = get_mapAttrToIndexIndex(attrMap, (const char *) connectionTuple[i].name, &attrIndex);
        if (status == CAPS_NOTFOUND) {
            printf("\tName %s not found in attribute map of capsConnect!!!!\n", connectionTuple[i].name);
            continue;
        } else if (status != CAPS_SUCCESS) return status;

        // Do we have a json string?
        if (strncmp(connectionTuple[i].value, "{", 1) == 0) {

            /*! \page feaConnection
             * \section jsonStringConnection JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  (e.g. "Value" = {"dofDependent": 1, "propertyType": "RigidBody"})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>connectionType = RigidBody</B> </li> <br>
             *  Type of connection to apply to a given capsConnect pair defined by "Connection Name" and the "groupName".
             *  Options: Mass (scalar), Spring (scalar), Damper (scalar), RigidBody, RigidBodyInterpolate.
             * </ul>
             * \elseif (ASTROS)
             * <ul>
             *  <li> <B>connectionType = RigidBody</B> </li> <br>
             *  Type of connection to apply to a given capsConnect pair defined by "Connection Name" and the "groupName".
             *  Options: Mass (scalar), Spring (scalar), RigidBody, RigidBodyInterpolate.
             * </ul>
             * \endif
             *
             */
            // Get connection Type
            keyWord = "connectionType";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                //{UnknownConnection, Mass, Spring, Damper, RigidBody}
                if      (strcasecmp(keyValue, "\"Mass\"")      == 0) connectionType = Mass;
                else if (strcasecmp(keyValue, "\"Spring\"")    == 0) connectionType = Spring;
                else if (strcasecmp(keyValue, "\"Damper\"")    == 0) connectionType = Damper;
                else if (strcasecmp(keyValue, "\"RigidBody\"") == 0) connectionType = RigidBody;
                else if (strcasecmp(keyValue, "\"RigidBodyInterpolate\"") == 0) connectionType = RigidBodyInterpolate;
                else {
                    printf("\tUnrecognized \"%s\" specified (%s) for Connection tuple %s, current options are "
                            "\"Mass, Spring, Damper, RigidBody, and RigidBodyInterpolate\"\n", keyWord,
                            keyValue,
                            connectionTuple[i].name);
                    AIM_FREE(keyValue);

                    return CAPS_NOTFOUND;
                }

            } else {

                printf("\tNo \"%s\" specified for Connection tuple %s, defaulting to RigidBody\n", keyWord,
                                                                                                   connectionTuple[i].name);
                connectionType = RigidBody;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>dofDependent = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the dependent end of rigid body connections (ex. 123 - translation in all three
             *  directions).
             * </ul>
             * \endif
             *
             */
            keyWord = "dofDependent";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &dofDependent);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }


            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>componentNumberStart = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the starting point of the connection for mass,
             *  spring, and damper elements (scalar) ( 0 <= Integer <= 6).
             * </ul>
             * \endif
             */
            keyWord = "componentNumberStart";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &componentNumberStart);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>componentNumberEnd= 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the ending point of the connection for mass,
             *  spring, damper elements (scalar), and rigid body interpolative connection ( 0 <= Integer <= 6).
             * </ul>
             * \endif
             *
             */
            keyWord = "componentNumberEnd";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &componentNumberEnd);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>stiffnessConst = 0.0 </B> </li> <br>
             *  Stiffness constant of a spring element (scalar).
             * </ul>
             * \endif
             */
            keyWord = "stiffnessConst";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &stiffnessConst);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>dampingConst = 0.0 </B> </li> <br>
             *  Damping coefficient/constant of a spring or damping element (scalar).
             * </ul>
             * \endif
             */
            keyWord = "dampingConst";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &dampingConst);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>stressCoeff = 0.0 </B> </li> <br>
             *  Stress coefficient of a spring element (scalar).
             * </ul>
             * \endif
             */
            keyWord = "stressCoeff";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &stressCoeff);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>mass = 0.0 </B> </li> <br>
             *  Mass of a mass element (scalar).
             * </ul>
             * \endif
             */
            keyWord = "mass";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &mass);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>glue = False </B> </li> <br>
             *  Turn on gluing for the connection.
             * </ul>
             * \endif
             */
            keyWord = "glue";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toBoolean(keyValue, &glue);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>glueNumMaster = 5 </B> </li> <br>
             *  Maximum number of the masters for a glue connections.
             * </ul>
             * \endif
             */
            keyWord = "glueNumMaster";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &glueNumMaster);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>glueSearchRadius = 0 </B> </li> <br>
             *  Search radius when looking for masters for a glue connections.
             * </ul>
             * \endif
             */
            keyWord = "glueSearchRadius";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &glueSearchRadius);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>weighting = 1 </B> </li> <br>
             *  Weighting factor for a rigid body interpolative connections.
             * </ul>
             * \endif
             */
            keyWord = "weighting";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &weighting);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsConnect</c> names on which to connect the nodes found with the
             *  tuple name ("Connection Name") to. (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             */
            numGroupName = 0;
            keyWord = "groupName";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                AIM_FREE(keyValue);
                AIM_NOTNULL(groupName, aimInfo, status);

                if (status != CAPS_SUCCESS) return status;
                if (glue == (int) true && connectionType != RigidBodyInterpolate) {
                    printf("\tInvalid connectionType while glue = True, setting glue to False!\n");
                    glue = (int) false;
                }

                if (glue == (int) true && connectionType == RigidBodyInterpolate) {
                    status = fea_glueMesh(aimInfo,
                                          &feaProblem->feaMesh,
                                          i+1, //connectionID
                                          connectionType,
                                          dofDependent,
                                          connectionTuple[i].name,
                                          numGroupName,
                                          groupName,
                                          attrMap,
                                          glueNumMaster,
                                          glueSearchRadius,
                                          &feaProblem->numConnect,
                                          &feaProblem->feaConnect);
                    if (status != CAPS_SUCCESS) return status;

                } else {


                    destNode = (int *) EG_alloc(feaProblem->feaMesh.numNode*sizeof(int));
                    if (destNode == NULL) {
                        (void) string_freeArray(numGroupName, &groupName);
                        return EGADS_MALLOC;
                    }

                    // Find the "source" node in the mesh
                    for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                        feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;

                        // If "source" doesn't match - continue
                        if (feaData->connectIndex != attrIndex) continue;

                        numDestNode = 0;

                        // Lets loop through the group names and create the connections
                        for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                            status = get_mapAttrToIndexIndex(attrMap, (const char *) groupName[groupIndex], &attrIndexDest);
                           if (status == CAPS_NOTFOUND) {
                                printf("\tName %s not found in attribute map of capsConnects!!!!\n", groupName[groupIndex]);
                                continue;
                            } else if (status != CAPS_SUCCESS) {
                                (void) string_freeArray(numGroupName, &groupName);
                                AIM_FREE(destNode);
                                return status;
                            }

                            // Find the "destination" node in the mesh
                            for (nodeIndexDest = 0; nodeIndexDest < feaProblem->feaMesh.numNode; nodeIndexDest++) {

                                feaDataDest = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndexDest].analysisData;

                                // If the "destination" doesn't match - continue
                                if (feaDataDest->connectIndex != attrIndexDest) continue;

                                destNode[numDestNode] = feaProblem->feaMesh.node[nodeIndexDest].nodeID;
                                numDestNode += 1;
                            } // End destination loop

                        } // End group loop

                        if (numDestNode <= 0) {
                            printf("\tNo destination nodes found for connection %s\n", connectionTuple[i].name);
                        } else {

                            status = fea_setConnection(connectionTuple[i].name,
                                                       connectionType,
                                                       i+1,
                                                       feaProblem->feaMesh.numElement,
                                                       dofDependent,
                                                       stiffnessConst,
                                                       dampingConst,
                                                       stressCoeff,
                                                       componentNumberStart,
                                                       componentNumberEnd,
                                                       feaProblem->feaMesh.node[nodeIndex].nodeID,
                                                       weighting,
                                                       componentNumberEnd,
                                                       numDestNode, destNode,
                                                       &feaProblem->numConnect,
                                                       &feaProblem->feaConnect);
                            if (status != CAPS_SUCCESS) {
                                (void) string_freeArray(numGroupName, &groupName);
                                if (destNode !=NULL) EG_free(destNode);
                                return status;
                            }
                        }
                    } // End source loop
                } // Glue ifelse

                if (destNode !=NULL) EG_free(destNode);
                destNode = NULL;

                status =  string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;

            } else {
                printf("\tNo \"%s\" specified for Connection tuple %s!\n", keyWord, connectionTuple[i].name);
            }

            // Create automatic connections from the "capsConnectLink" tag
            printf("\tLooking for automatic connections from the use of capsConnectLink for %s\n", connectionTuple[i].name);

            destNode = (int *) EG_alloc(feaProblem->feaMesh.numNode*sizeof(int));
            if (destNode == NULL) return EGADS_MALLOC;

            // Find the "source" node in the mesh
            for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;

                // If "source" doesn't match - continue
                if (feaData->connectIndex != attrIndex) continue;

                numDestNode = 0;
                // Find the "destination" node in the mesh
                for (nodeIndexDest = 0; nodeIndexDest < feaProblem->feaMesh.numNode; nodeIndexDest++) {

                    feaDataDest = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndexDest].analysisData;

                    // If the "destination" doesn't match - continue
                    if (feaDataDest->connectLinkIndex != attrIndex) continue;

                    destNode[numDestNode] = feaProblem->feaMesh.node[nodeIndexDest].nodeID;
                    numDestNode += 1;
                } // End destination loop

                if (numDestNode > 0) {
                    status = fea_setConnection(connectionTuple[i].name,
                                           connectionType,
                                           i+1,
                                           feaProblem->feaMesh.numElement,
                                           dofDependent,
                                           stiffnessConst,
                                           dampingConst,
                                           stressCoeff,
                                           componentNumberStart,
                                           componentNumberEnd,
                                           feaProblem->feaMesh.node[nodeIndex].nodeID,
                                           weighting,
                                           componentNumberEnd,
                                           numDestNode, destNode,
                                           &feaProblem->numConnect,
                                           &feaProblem->feaConnect);
                    if (status != CAPS_SUCCESS) {
                        if (destNode !=NULL) EG_free(destNode);
                        return status;
                    }
                }

                if (numDestNode > 0) {
                    printf("\t%d automatic connections were made for capsConnect %s (node id %d)\n", numDestNode,
                                                                                                     connectionTuple[i].name,
                                                                                                     feaProblem->feaMesh.node[nodeIndex].nodeID);
                }
            } // End source loop

            if (destNode !=NULL) EG_free(destNode);
            destNode = NULL;

        } else {

            /*! \page feaConnection
             * \section keyStringConnection Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined connection lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */

            // Call some look up table maybe?
            printf("\tError: Connection tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;

        }
    }

    if (keyValue != NULL) {
        EG_free(keyValue);
        keyValue = NULL;
    }

    printf("\tDone getting FEA connections\n");
    status = CAPS_SUCCESS;
cleanup:
    return status;
}

#ifdef FEA_GETCONNECTIONORIG_DEPREICATE
// Get the Connections properties from a capsTuple and create connections based on the mesh
static int fea_getConnectionOrig(int numConnectionTuple,
                      capsTuple connectionTuple[],
                      mapAttrToIndexStruct *attrMap,
                      feaProblemStruct *feaProblem) {

    /* \page feaConnection FEA Connection
     * Structure for the connection tuple  = ("Connection Name", "Value").
     * "Connection Name" defines the reference name to the capsConnect being specified and denotes the "source" node
     * for the connection.
     *  The "Value" can either be a JSON String dictionary (see Section \ref jsonStringConnection) or a single string keyword
     *  (see Section \ref keyStringConnection).
     */

    int status; //Function return

    int i, groupIndex, nodeIndex, nodeIndexDest, counter = 0; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    int numGroupName = 0;
    char **groupName = NULL;
    int attrIndex, attrIndexDest;

    // Values to set
    feaConnectionTypeEnum connectionType;
    int dofDependent = 0, componentNumberStart = 0, componentNumberEnd = 0;
    double stiffnessConst = 0.0, dampingConst  = 0.0, stressCoeff = 0.0, mass = 0.0;

    feaMeshDataStruct *feaData, *feaDataDest;

    // Destroy our support structures coming in if aren't 0 and NULL already
    if (feaProblem->feaConnect != NULL) {
        for (i = 0; i < feaProblem->numConnect; i++) {
            status = destroy_feaConnectionStruct(&feaProblem->feaConnect[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (feaProblem->feaConnect != NULL) EG_free(feaProblem->feaConnect);
    feaProblem->feaConnect = NULL;
    feaProblem->numConnect = 0;

    printf("\nGetting FEA connections.......\n");

    printf("\tNumber of connection tuples - %d\n", numConnectionTuple);

    // Loop through tuples and fill out the support structures
    for (i = 0; i < numConnectionTuple; i++) {

        // Reset defaults
        dofDependent  = componentNumberStart = componentNumberEnd = 0;
        stiffnessConst = dampingConst = stressCoeff = mass = 0.0;

        printf("\tConnection name - %s\n", connectionTuple[i].name );

        // Look for connection name in connection map
        status = get_mapAttrToIndexIndex(attrMap, (const char *) connectionTuple[i].name, &attrIndex);
        if (status == CAPS_NOTFOUND) {
            printf("\tName %s not found in attribute map of capsConnect!!!!\n", connectionTuple[i].name);
            continue;
        } else if (status != CAPS_SUCCESS) return status;

        // Do we have a json string?
        if (strncmp(connectionTuple[i].value, "{", 1) == 0) {

            /* \page feaConnection
             * \section jsonStringConnection JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || ASTROS)
             *  (e.g. "Value" = {"dofDependent": 1, "propertyType": "RigidBody"})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>connectionType = RigidBody</B> </li> <br>
             *  Type of connection to apply to a given capsConnect pair defined by "Connection Name" and the "groupName".
             *  Options: Mass (scalar), Spring (scalar), Damper (scalar), RigidBody.
             * </ul>
             * \elseif (ASTROS)
             * <ul>
             *  <li> <B>connectionType = RigidBody</B> </li> <br>
             *  Type of connection to apply to a given capsConnect pair defined by "Connection Name" and the "groupName".
             *  Options: Mass (scalar), Spring (scalar), RigidBody.
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             * <li> <B>connectionType = RigidBody</B> </li> <br>
             *  Type of connection to apply to a given capsConnect pair defined by "Connection Name" and the "groupName".
             *  Options: Mass (scalar), Spring (scalar), Damper (scalar), RigidBody.
             *  </ul>
             * \endif
             *
             */

            // Get connection Type
            keyWord = "connectionType";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                //{UnknownConnection, Mass, Spring, Damper, RigidBody}
                if      (strcasecmp(keyValue, "\"Mass\"")      == 0) connectionType = Mass;
                else if (strcasecmp(keyValue, "\"Spring\"")    == 0) connectionType = Spring;
                else if (strcasecmp(keyValue, "\"Damper\"")    == 0) connectionType = Damper;
                else if (strcasecmp(keyValue, "\"RigidBody\"") == 0) connectionType = RigidBody;
                else {
                    printf("\tUnrecognized \"%s\" specified (%s) for Connection tuple %s, current options are "
                            "\"Mass, Spring, Damper, and RigidBody\"\n", keyWord,
                            keyValue,
                            connectionTuple[i].name);
                    if (keyValue != NULL) {
                        EG_free(keyValue);
                        keyValue = NULL;
                    }

                    return CAPS_NOTFOUND;
                }

            } else {

                printf("\tNo \"%s\" specified for Connection tuple %s, defaulting to RigidBody\n", keyWord,
                                                                                                   connectionTuple[i].name);
                connectionType = RigidBody;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>dofDependent = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the dependent end of rigid body connections (ex. 123 - translation in all three
             *  directions).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>dofDependent = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the dependent end of rigid body connections (ex. 123 - translation in all three
             *  directions).
             * </ul>
             * \endif
             *
             */
            keyWord = "dofDependent";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &dofDependent);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }


            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>componentNumberStart = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the starting point of the connection for mass,
             *  spring, and damper elements (scalar) ( 0 <= Integer <= 6).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>componentNumberStart = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the starting point of the connection for mass,
             *  spring, and damper elements (scalar) ( 0 <= Integer <= 6).
             * </ul>
             * \endif
             *
             */
            keyWord = "componentNumberStart";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &componentNumberStart);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>componentNumberEnd= 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the ending point of the connection for mass,
             *  spring, and damper elements (scalar) ( 0 <= Integer <= 6).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>componentNumberEnd = 0 </B> </li> <br>
             *  Component numbers / degrees of freedom of the ending point of the connection for mass,
             *  spring, and damper elements (scalar) ( 0 <= Integer <= 6).
             * </ul>
             * \endif
             *
             */
            keyWord = "componentNumberEnd";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &componentNumberEnd);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>stiffnessConst = 0.0 </B> </li> <br>
             *  Stiffness constant of a spring element (scalar).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>stiffnessConst = 0.00 </B> </li> <br>
             *  Stiffness constant of a spring element (scalar).
             * </ul>
             * \endif
             *
             */
            keyWord = "stiffnessConst";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &stiffnessConst);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>dampingConst = 0.0 </B> </li> <br>
             *  Damping coefficient/constant of a spring or damping element (scalar).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>dampingConst = 0.0 </B> </li> <br>
             *  Damping constant of a spring or damping element (scalar).
             * </ul>
             * \endif
             *
             */
            keyWord = "dampingConst";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &dampingConst);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>stressCoeff = 0.0 </B> </li> <br>
             *  Stress coefficient of a spring element (scalar).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>stressCoeff = 0.0 </B> </li> <br>
             * Stress coefficient of a spring element (scalar).
             * </ul>
             * \endif
             *
             */
            keyWord = "stressCoeff";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &stressCoeff);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>mass = 0.0 </B> </li> <br>
             *  Mass of a mass element (scalar).
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             *  <li> <B>mass = 0.0 </B> </li> <br>
             * Mass of a mass element (scalar).
             * </ul>
             * \endif
             *
             */
            keyWord = "mass";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &mass);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /* \page feaConnection
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsConnect</c> names on which to connect the nodes found with the
             *  tuple name ("Connection Name") to. (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \elseif ABAQUS
             *
             * Something else ....
             *
             * \else
             * <ul>
             * <li> <B>groupName = "(no default)"</B> </li> <br>
             * Single or list of <c>capsConnect</c> names on which to connect the nodes found with the
             *  tuple name ("Connection Name") to. (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             *
             */
            numGroupName = 0;
            keyWord = "groupName";
            status = search_jsonDictionary( connectionTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;

                if (status != CAPS_SUCCESS) return status;

                // Lets loop through the group names and create the connections
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    status = get_mapAttrToIndexIndex(attrMap, (const char *) groupName[groupIndex], &attrIndexDest);
                    if (status == CAPS_NOTFOUND) {
                        printf("\tName %s not found in attribute map of capsConnects!!!!\n", groupName[groupIndex]);
                        continue;
                    } else if (status != CAPS_SUCCESS) {
                        (void) string_freeArray(numGroupName, &groupName);
                        return status;
                    }

                    // Find the "source" node in the mesh
                    for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                        feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;

                        // If "source" doesn't match - continue
                        if (feaData->connectIndex != attrIndex) continue;

                        // Find the "destination" node in the mesh
                        for (nodeIndexDest = 0; nodeIndexDest < feaProblem->feaMesh.numNode; nodeIndexDest++) {

                            feaDataDest = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndexDest].analysisData;

                            // If the "destination" doesn't match - continue
                            if (feaDataDest->connectIndex != attrIndexDest) continue;

                            feaProblem->numConnect += 1;

                            if (feaProblem->numConnect == 1) {
                                feaProblem->feaConnect = (feaConnectionStruct *) EG_alloc(sizeof(feaConnectionStruct));
                            } else {
                                feaProblem->feaConnect = (feaConnectionStruct *) EG_reall(feaProblem->feaConnect ,
                                                                                          feaProblem->numConnect*sizeof(feaConnectionStruct));
                            }

                            if (feaProblem->feaConnect == NULL) {
                                feaProblem->numConnect = 0;
                                (void) string_freeArray(numGroupName, &groupName);
                                return EGADS_MALLOC;
                            }

                            status = initiate_feaConnectionStruct(&feaProblem->feaConnect[feaProblem->numConnect-1]);
                            if (status != CAPS_SUCCESS) return status;

                            feaProblem->feaConnect[feaProblem->numConnect-1].name = (char *) EG_alloc((strlen(connectionTuple[i].name) + 1)*sizeof(char));
                            if (feaProblem->feaConnect[feaProblem->numConnect-1].name == NULL) {
                                (void) string_freeArray(numGroupName, &groupName);
                                return EGADS_MALLOC;
                            }

                            memcpy(feaProblem->feaConnect[feaProblem->numConnect-1].name,
                                   connectionTuple[i].name,
                                   strlen(connectionTuple[i].name)*sizeof(char));
                            feaProblem->feaConnect[feaProblem->numConnect-1].name[strlen(connectionTuple[i].name)] = '\0';

                            feaProblem->feaConnect[feaProblem->numConnect-1].connectionID = i+1; // ConnectionTuple index

                            feaProblem->feaConnect[feaProblem->numConnect-1].connectionType = connectionType;
                            feaProblem->feaConnect[feaProblem->numConnect-1].elementID = feaProblem->numConnect + feaProblem->feaMesh.numElement;
                            feaProblem->feaConnect[feaProblem->numConnect-1].connectivity[0] = feaProblem->feaMesh.node[nodeIndex].nodeID;
                            feaProblem->feaConnect[feaProblem->numConnect-1].connectivity[1] = feaProblem->feaMesh.node[nodeIndexDest].nodeID;

                            feaProblem->feaConnect[feaProblem->numConnect-1].dofDependent = dofDependent;
                            feaProblem->feaConnect[feaProblem->numConnect-1].stiffnessConst = stiffnessConst;
                            feaProblem->feaConnect[feaProblem->numConnect-1].dampingConst = dampingConst;
                            feaProblem->feaConnect[feaProblem->numConnect-1].stressCoeff = stressCoeff;
                            feaProblem->feaConnect[feaProblem->numConnect-1].mass = mass;
                            feaProblem->feaConnect[feaProblem->numConnect-1].componentNumberStart = componentNumberStart;
                            feaProblem->feaConnect[feaProblem->numConnect-1].componentNumberEnd = componentNumberEnd;
                        } // End destination loop
                    } // End source loop
                } // End group loop

                status =  string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;

            } else {
                printf("\tNo \"%s\" specified for Connection tuple %s!\n", keyWord,
                        connectionTuple[i].name);
            }


            // Create automatic connections from the "capsConnectLink" tag
            printf("\tLooking for automatic connections from the use of capsConnectLink for %s\n", connectionTuple[i].name);

            // Find the "source" node in the mesh
            for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;

                // If "source" doesn't match - continue
                if (feaData->connectIndex != attrIndex) continue;

                counter = 0;
                // Find the "destination" node in the mesh
                for (nodeIndexDest = 0; nodeIndexDest < feaProblem->feaMesh.numNode; nodeIndexDest++) {

                    feaDataDest = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndexDest].analysisData;

                    // If the "destination" doesn't match - continue
                    if (feaDataDest->connectLinkIndex != attrIndex) continue;

                    counter +=1;

                    feaProblem->numConnect += 1;

                    if (feaProblem->numConnect == 1) {
                        feaProblem->feaConnect = (feaConnectionStruct *) EG_alloc(sizeof(feaConnectionStruct));
                    } else {
                        feaProblem->feaConnect = (feaConnectionStruct *) EG_reall(feaProblem->feaConnect ,
                                                                                  feaProblem->numConnect*sizeof(feaConnectionStruct));
                    }

                    if (feaProblem->feaConnect == NULL) {
                        feaProblem->numConnect = 0;
                        return EGADS_MALLOC;
                    }

                    feaProblem->feaConnect[feaProblem->numConnect-1].name = (char *) EG_alloc((strlen(connectionTuple[i].name) + 1)*sizeof(char));
                    if (feaProblem->feaConnect[feaProblem->numConnect-1].name == NULL) return EGADS_MALLOC;

                    memcpy(feaProblem->feaConnect[feaProblem->numConnect-1].name,
                           connectionTuple[i].name,
                           strlen(connectionTuple[i].name)*sizeof(char));
                    feaProblem->feaConnect[feaProblem->numConnect-1].name[strlen(connectionTuple[i].name)] = '\0';

                    feaProblem->feaConnect[feaProblem->numConnect-1].connectionID = i+1; // ConnectionTuple index

                    feaProblem->feaConnect[feaProblem->numConnect-1].connectionType = connectionType;
                    feaProblem->feaConnect[feaProblem->numConnect-1].elementID = feaProblem->numConnect + feaProblem->feaMesh.numElement;
                    feaProblem->feaConnect[feaProblem->numConnect-1].connectivity[0] = feaProblem->feaMesh.node[nodeIndex].nodeID;
                    feaProblem->feaConnect[feaProblem->numConnect-1].connectivity[1] = feaProblem->feaMesh.node[nodeIndexDest].nodeID;

                    feaProblem->feaConnect[feaProblem->numConnect-1].dofDependent = dofDependent;
                    feaProblem->feaConnect[feaProblem->numConnect-1].stiffnessConst = stiffnessConst;
                    feaProblem->feaConnect[feaProblem->numConnect-1].dampingConst = dampingConst;
                    feaProblem->feaConnect[feaProblem->numConnect-1].stressCoeff = stressCoeff;
                    feaProblem->feaConnect[feaProblem->numConnect-1].componentNumberStart = componentNumberStart;
                    feaProblem->feaConnect[feaProblem->numConnect-1].componentNumberEnd = componentNumberEnd;
                } // End destination loop

                printf("\t%d automatic connections were made for capsConnect %s\n", counter,
                                                                                    connectionTuple[i].name);
            } // End source loop
        } else {

            /* \page feaConnection
             * \section keyStringConnection Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined connection lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */

            // Call some look up table maybe?
            printf("\tError: Connection tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;

        }
    }

    if (keyValue != NULL) {
        EG_free(keyValue);
        keyValue = NULL;
    }

    printf("\tDone getting FEA connections\n");
    return CAPS_SUCCESS;
}
#endif

// Get the load properties from a capsTuple
int fea_getLoad(void *aimInfo,
                int numLoadTuple,
                capsTuple loadTuple[],
                mapAttrToIndexStruct *attrMap,
                feaProblemStruct *feaProblem) {

    /*! \page feaLoad FEA Load
     * Structure for the load tuple  = ("Load Name", "Value").
     * "Load Name" defines the reference name for the load being specified.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringLoad) or a single string keyword
     *	(see Section \ref keyStringLoad).
     */

    int status; //Function return

    int i, groupIndex, attrIndex, nodeIndex, elementIndex; // Indexing

    char *keyValue = NULL; // Key values from tuple searches
    char *keyWord = NULL; // Key words to find in the tuples

    int numGroupName = 0;
    char **groupName = NULL;

    char *tempString = NULL; // Temporary string holder

    feaMeshDataStruct *feaData;

    // Destroy our load structures coming in if aren't 0 and NULL already
    if (feaProblem->feaLoad != NULL) {
        for (i = 0; i < feaProblem->numLoad; i++) {
            status = destroy_feaLoadStruct(&feaProblem->feaLoad[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    AIM_FREE(feaProblem->feaLoad);
    feaProblem->numLoad = 0;

    printf("\nGetting FEA loads.......\n");

    feaProblem->numLoad = numLoadTuple;

    printf("\tNumber of loads - %d\n", feaProblem->numLoad);

    // Allocate loads
    if (feaProblem->numLoad > 0) {
        feaProblem->feaLoad = (feaLoadStruct *) EG_alloc(feaProblem->numLoad * sizeof(feaLoadStruct));

        if (feaProblem->feaLoad == NULL ) return EGADS_MALLOC;
    }

    // Initiate loads to default values
    for (i = 0; i < feaProblem->numLoad; i++) {
        status = initiate_feaLoadStruct(&feaProblem->feaLoad[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    // Loop through tuples and fill out the load structures
    for (i = 0; i < feaProblem->numLoad; i++) {

        printf("\tLoad name - %s\n", loadTuple[i].name );

        // Set load name to tuple attribute name
        feaProblem->feaLoad[i].name = (char *) EG_alloc((strlen(loadTuple[i].name) + 1)*sizeof(char));
        if (feaProblem->feaLoad[i].name == NULL) return EGADS_MALLOC;

        memcpy(feaProblem->feaLoad[i].name, loadTuple[i].name, strlen(loadTuple[i].name)*sizeof(char));
        feaProblem->feaLoad[i].name[strlen(loadTuple[i].name)] = '\0';

        // Set load id -> 1 bias
        feaProblem->feaLoad[i].loadID = i+1;

        // Do we have a json string?
        if (strncmp(loadTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", constraintTuple[i].value);

            /*! \page feaLoad
             * \section jsonStringLoad JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS || TACS)
             *  (e.g. "Value" = {"groupName": "plate", "loadType": "Pressure", "pressureForce": 2000000.0})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>loadType = "(no default)"</B> </li> <br>
             *  Type of load. Options: "GridForce", "GridMoment", "Rotational", "Thermal",
             *  "Pressure", "PressureDistribute", "PressureExternal", "TermalExternal", "Gravity".
             * </ul>
             * \elseif ABAQUS
              * <ul>
             *   <li> <B>loadType = "(no default)"</B> </li> <br>
             *   Type of load. Options: "Pressure", "Gravity".
             * </ul>
             * \elseif HSM
             *  <ul>
             *   <li> <B>loadType = "(no default)"</B> </li> <br>
             *   Type of load. Options: "GridForce", "GridMoment", "LineForce", "LineMoment", "Rotational",
             *   "Pressure", "PressureDistribute", "PressureExternal", "Gravity".
             * </ul>
             * \else
             * \endif
             *
             */
            // Get load Type
            keyWord = "loadType";
            status = search_jsonDictionary(loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                //{UnknownLoad, GridForce, GridMoment, Gravity, Pressure, Rotational, Thermal}
                if      (strcasecmp(keyValue, "\"GridForce\"")          == 0) feaProblem->feaLoad[i].loadType = GridForce;
                else if (strcasecmp(keyValue, "\"GridMoment\"")         == 0) feaProblem->feaLoad[i].loadType = GridMoment;
                else if (strcasecmp(keyValue, "\"LineForce\"")          == 0) feaProblem->feaLoad[i].loadType = LineForce;
                else if (strcasecmp(keyValue, "\"LineMoment\"")         == 0) feaProblem->feaLoad[i].loadType = LineMoment;
                else if (strcasecmp(keyValue, "\"Rotational\"")         == 0) feaProblem->feaLoad[i].loadType = Rotational;
                else if (strcasecmp(keyValue, "\"Thermal\"")            == 0) feaProblem->feaLoad[i].loadType = Thermal;
                else if (strcasecmp(keyValue, "\"Pressure\"")           == 0) feaProblem->feaLoad[i].loadType = Pressure;
                else if (strcasecmp(keyValue, "\"PressureDistribute\"") == 0) feaProblem->feaLoad[i].loadType = PressureDistribute;
                else if (strcasecmp(keyValue, "\"PressureExternal\"")   == 0) feaProblem->feaLoad[i].loadType = PressureExternal;
                else if (strcasecmp(keyValue, "\"ThermalExternal\"")    == 0) feaProblem->feaLoad[i].loadType = ThermalExternal;
                else if (strcasecmp(keyValue, "\"Gravity\"")            == 0) feaProblem->feaLoad[i].loadType = Gravity;
                else {
                    printf("\tUnrecognized \"%s\" specified (%s) for Load tuple %s\n", keyWord,
                            keyValue,
                            loadTuple[i].name);
                    if (keyValue != NULL) EG_free(keyValue);
                    return CAPS_NOTFOUND;
                }

            } else {
                printf("\t\"loadType\" variable not found in tuple %s, this is required input!!\n", loadTuple[i].name);

                if (keyValue != NULL) EG_free(keyValue);
                return status;
            }

            AIM_FREE(keyValue);

            // Get load node/element set
            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS || TACS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsLoad</c> names on which to apply the load
             *  (e.g. "Name1" or ["Name1","Name2",...]. If not provided, the load tuple name will be
             *  used.
             * </ul>
             * \endif
             *
             */
            keyWord = "groupName";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) EG_free(keyValue);
                keyValue = NULL;

                if (status != CAPS_SUCCESS) return status;

            } else {

                printf("\tNo \"%s\" specified for Load tuple %s, going to use load name\n", keyWord,
                        loadTuple[i].name);

                status = string_toStringDynamicArray(loadTuple[i].name, &numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;

            }

            // Determine how many element/point loads we have
            feaProblem->feaLoad[i].numGridID = 0;
            feaProblem->feaLoad[i].numElementID = 0;
            for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                // Do nothing for PressureExternal, ThermalExternal, and Gravity loads
                if (feaProblem->feaLoad[i].loadType == PressureExternal) continue;
                if (feaProblem->feaLoad[i].loadType == ThermalExternal) continue;
                if (feaProblem->feaLoad[i].loadType == Gravity) continue;

                status = get_mapAttrToIndexIndex(attrMap, (const char *) groupName[groupIndex], &attrIndex);

                if (status == CAPS_NOTFOUND) {
                    printf("\tName %s not found in attribute map of capsLoads!!!!\n", groupName[groupIndex]);
                    continue;

                } else if (status != CAPS_SUCCESS) {

                    (void)  string_freeArray(numGroupName, &groupName);
                    return status;
                }

                //{UnknownLoad, GridForce, GridMoment, Gravity, Pressure, Rotational, Thermal}
                if (feaProblem->feaLoad[i].loadType == GridForce  ||
                    feaProblem->feaLoad[i].loadType == GridMoment ||
                    feaProblem->feaLoad[i].loadType == Rotational ||
                    feaProblem->feaLoad[i].loadType == Thermal) {

                    // Now lets loop through the grid to see how many grid points have the attrIndex
                    for (nodeIndex = 0; nodeIndex < feaProblem->feaMesh.numNode; nodeIndex++ ) {

                        if (feaProblem->feaMesh.node[nodeIndex].analysisType == MeshStructure) {
                            feaData = (feaMeshDataStruct *) feaProblem->feaMesh.node[nodeIndex].analysisData;
                        } else {
                            continue;
                        }

                        if (feaData->loadIndex == attrIndex) {

                            feaProblem->feaLoad[i].numGridID += 1;

                            // Allocate/Re-allocate grid ID array
                            if (feaProblem->feaLoad[i].numGridID == 1) {

                                feaProblem->feaLoad[i].gridIDSet = (int *) EG_alloc(feaProblem->feaLoad[i].numGridID*sizeof(int));

                            } else {

                                feaProblem->feaLoad[i].gridIDSet = (int *) EG_reall(feaProblem->feaLoad[i].gridIDSet,
                                                                                    feaProblem->feaLoad[i].numGridID*sizeof(int));
                            }

                            if (feaProblem->feaLoad[i].gridIDSet == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status == CAPS_SUCCESS) return status;

                                return EGADS_MALLOC;
                            }

                            // Set grid ID value -> 1 bias
                            feaProblem->feaLoad[i].gridIDSet[feaProblem->feaLoad[i].numGridID-1] = feaProblem->feaMesh.node[nodeIndex].nodeID;

                            //printf("GroupName = %s %d\n", groupName[groupIndex], feaProblem->feaMesh.node[nodeIndex].nodeID);

                        }
                    }

                } else if (feaProblem->feaLoad[i].loadType == LineForce ||
                           feaProblem->feaLoad[i].loadType == LineMoment ) {

                  for (elementIndex = 0; elementIndex < feaProblem->feaMesh.numElement; elementIndex++) {

                      if (feaProblem->feaMesh.element[elementIndex].elementType != Line) {
                          continue;
                      }

                      if (feaProblem->feaMesh.element[elementIndex].analysisType == MeshStructure) {
                          feaData = (feaMeshDataStruct *) feaProblem->feaMesh.element[elementIndex].analysisData;
                      } else {
                          continue;
                      }

                      if (feaData->loadIndex == attrIndex) {

                          feaProblem->feaLoad[i].numElementID += 1;

                          feaProblem->feaLoad[i].elementIDSet = (int *) EG_reall(feaProblem->feaLoad[i].elementIDSet,
                                                                                 feaProblem->feaLoad[i].numElementID*sizeof(int));

                          if (feaProblem->feaLoad[i].elementIDSet == NULL) {
                              string_freeArray(numGroupName, &groupName);
                              return EGADS_MALLOC;
                          }

                          // Set element ID value -> 1 bias
                          feaProblem->feaLoad[i].elementIDSet[feaProblem->feaLoad[i].numElementID-1] = feaProblem->feaMesh.element[elementIndex].elementID;
                      }
                  }

                } else if (feaProblem->feaLoad[i].loadType == LineForce ||
                           feaProblem->feaLoad[i].loadType == LineMoment ||
                           feaProblem->feaLoad[i].loadType == Pressure ||
                           feaProblem->feaLoad[i].loadType == PressureDistribute) {

                    // Now lets loop through the elements to see how many elements have the attrIndex

                    // Element types - CTRIA3, CTRIA3K, CQUAD4, CQUAD4K
                    for (elementIndex = 0; elementIndex < feaProblem->feaMesh.numElement; elementIndex++) {

                        if (feaProblem->feaMesh.element[elementIndex].elementType != Triangle &&
                            feaProblem->feaMesh.element[elementIndex].elementType != Triangle_6 &&
                            feaProblem->feaMesh.element[elementIndex].elementType != Quadrilateral &&
                            feaProblem->feaMesh.element[elementIndex].elementType != Quadrilateral_8) {
                            continue;
                        }

                        if (feaProblem->feaMesh.element[elementIndex].analysisType == MeshStructure) {
                            feaData = (feaMeshDataStruct *) feaProblem->feaMesh.element[elementIndex].analysisData;
                        } else {
                            continue;
                        }

                        if (feaData->loadIndex == attrIndex) {

                            feaProblem->feaLoad[i].numElementID += 1;

                            // Allocate/Re-allocate element ID array
                            if (feaProblem->feaLoad[i].numElementID == 1) {

                                feaProblem->feaLoad[i].elementIDSet = (int *) EG_alloc(feaProblem->feaLoad[i].numElementID*sizeof(int));

                            } else {

                                feaProblem->feaLoad[i].elementIDSet = (int *) EG_reall(feaProblem->feaLoad[i].elementIDSet,
                                                                                       feaProblem->feaLoad[i].numElementID*sizeof(int));
                            }

                            if (feaProblem->feaLoad[i].elementIDSet == NULL) {
                                string_freeArray(numGroupName, &groupName);
                                return EGADS_MALLOC;
                            }

                            // Set element ID value -> 1 bias
                            feaProblem->feaLoad[i].elementIDSet[feaProblem->feaLoad[i].numElementID-1] = feaProblem->feaMesh.element[elementIndex].elementID;
                        }
                    }

                }
            } // End attr search loop

            status = string_freeArray(numGroupName, &groupName);
            if (status != CAPS_SUCCESS) return status;
            groupName = NULL;

            // Free keyValue (just in case)
            AIM_FREE(keyValue);

            //Fill up load properties

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>loadScaleFactor = 1.0 </B> </li> <br>
             *  Scale factor to use when combining loads.
             * </ul>
             * \endif
             */
            keyWord = "loadScaleFactor";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].loadScaleFactor);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || TACS)
             * <ul>
             *  <li> <B>forceScaleFactor = 0.0 </B> </li> <br>
             *  Overall scale factor for the force for a "GridForce" load.
             * </ul>
             * \endif
             */
            keyWord = "forceScaleFactor";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].forceScaleFactor);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || TACS)
             * <ul>
             *  <li> <B>directionVector = [0.0, 0.0, 0.0] </B> </li> <br>
             *  X-, y-, and z- components of the force vector for a "GridForce", "GridMoment", or "Gravity" load.
             * </ul>
             * \elseif ABAQUS
             * <ul>
             *  <li> <B>directionVector = [0.0, 0.0, 0.0] </B> </li> <br>
             *  X-, y-, and z- components of the force vector for a "Gravity" load.
             * </ul>
             * \endif
             *
             */
            keyWord = "directionVector";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleArray(keyValue,
                        (int) sizeof(feaProblem->feaLoad[i].directionVector)/sizeof(double),
                        feaProblem->feaLoad[i].directionVector);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || TACS)
             * <ul>
             *  <li> <B>momentScaleFactor = 0.0 </B> </li> <br>
             *  Overall scale factor for the moment for a "GridMoment" load.
             * </ul>
             * \endif
             *
             */
            keyWord = "momentScaleFactor";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].momentScaleFactor);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || HSM || ASTROS || ABAQUS || TACS)
             * <ul>
             *  <li> <B>gravityAcceleration = 0.0 </B> </li> <br>
             *  Acceleration value for a "Gravity" load.
             * </ul>
             * \endif
             */
            keyWord = "gravityAcceleration";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].gravityAcceleration);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if MYSTRAN
             * <ul>
             *  <li> <B>pressureForce = 0.0 </B> </li> <br>
             *  Uniform pressure force for a "Pressure" load (only applicable to 2D elements).
             * </ul>
             * \elseif (NASTRAN || HSM || ASTROS || ABAQUS || TACS)
             *  <ul>
             *  <li> <B>pressureForce = 0.0 </B> </li> <br>
             *  Uniform pressure force for a "Pressure" load.
             * </ul>
             * \endif
             */
            keyWord = "pressureForce";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].pressureForce);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if MYSTRAN
             * <ul>
             *  <li> <B>pressureDistributeForce = [0.0, 0.0, 0.0, 0.0] </B> </li> <br>
             *  Distributed pressure force for a "PressureDistribute" load (only applicable to 2D elements). The four values
             *  correspond to the 4 (quadrilateral elements) or 3 (triangle elements) node locations.
             * </ul>
             * \elseif (NASTRAN || HSM || ASTROS || TACS)
             * <ul>
             *  <li> <B>pressureDistributeForce = [0.0, 0.0, 0.0, 0.0] </B> </li> <br>
             *  Distributed pressure force for a "PressureDistribute" load. The four values
             *  correspond to the 4 (quadrilateral elements) or 3 (triangle elements) node locations.
             * </ul>
             * \endif
             */
            keyWord = "pressureDistributeForce";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleArray(keyValue,
                        (int) sizeof(feaProblem->feaLoad[i].pressureDistributeForce)/sizeof(double),
                        feaProblem->feaLoad[i].pressureDistributeForce);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>angularVelScaleFactor = 0.0 </B> </li> <br>
             *  An overall scale factor for the angular velocity in revolutions per unit time for a "Rotational" load.
             * </ul>
             * \elseif HSM
             * <ul>
             *  <li> <B>angularVelScaleFactor = 0.0 </B> </li> <br>
             *  An overall scale factor for the angular velocity in revolutions per unit time for a "Rotational" load - applied in a global sense.
             * </ul>
             * \endif
             */
            keyWord = "angularVelScaleFactor";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].angularVelScaleFactor);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>angularAccScaleFactor = 0.0 </B> </li> <br>
             *  An overall scale factor for the angular acceleration in revolutions per unit time squared for a "Rotational" load.
             * </ul>
             * \elseif HSM
             * <ul>
             *  <li> <B>angularAccScaleFactor = 0.0 </B> </li> <br>
             *  An overall scale factor for the angular acceleration in revolutions per unit time squared for a "Rotational" load - applied in a global sense.
             * </ul>
             * \endif
             */
            keyWord = "angularAccScaleFactor";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].angularAccScaleFactor);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>coordinateSystem = "(no default)" </B> </li> <br>
             *  Name of coordinate system in which defined force components are in reference to. If no value
             *  is provided the global system is assumed.
             * </ul>
             * \endif
             *
             */
            keyWord = "coordinateSystem";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                for (attrIndex = 0; attrIndex < feaProblem->numCoordSystem; attrIndex++) {

                    if (strcasecmp(tempString, feaProblem->feaCoordSystem[attrIndex].name) == 0) {
                        feaProblem->feaLoad[i].coordSystemID  = feaProblem->feaCoordSystem[attrIndex].coordSystemID;
                        break;
                    }
                }

                if (feaProblem->feaLoad[i].coordSystemID == 0) {
                    printf("\tCoordinate system %s not found, defaulting to global system!!", keyValue);
                }

                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                AIM_FREE(tempString);
            }


            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>temperature = 0.0 </B> </li> <br>
             *  Temperature at a given node for a "Temperature" load.
             * </ul>
             * \endif
             */
            keyWord = "temperature";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].temperature);
                AIM_FREE(keyValue);
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaLoad
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || TACS)
             * <ul>
             *  <li> <B>temperatureDefault = 0.0 </B> </li> <br>
             *  Default temperature at a node not explicitly being used for a "Temperature" load.
             * </ul>
             * \endif
             */
            keyWord = "temperatureDefault";
            status = search_jsonDictionary( loadTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaLoad[i].temperatureDefault);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }


        } else { // if not a JSON string search

            /*! \page feaLoad
             * \section keyStringLoad Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined load lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             */

            // Call some look up table maybe?
            printf("\tError: Load tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;
        }
    }

    /// Free keyValue and tempString (just in case)
    if (keyValue != NULL) EG_free(keyValue);

    if (tempString != NULL) EG_free(tempString);

    printf("\tDone getting FEA loads\n");
    status = CAPS_SUCCESS;
cleanup:
    return status;
}

// Get the analysis properties from a capsTuple
int fea_getAnalysis(void *aimInfo,
                    int numAnalysisTuple,
                    capsTuple analysisTuple[],
                    feaProblemStruct *feaProblem) {

    /*! \page feaAnalysis FEA Analysis
     * Structure for the analysis tuple  = (`Analysis Name', `Value').
     * 'Analysis Name' defines the reference name for the analysis being specified.
     *	The "Value" can either be a JSON String dictionary (see Section \ref jsonStringAnalysis) or a single string keyword
     *	(see Section \ref keyStringAnalysis).
     */

    int status; //Function return

    int i, j, groupIndex, attrIndex; // Indexing


    char *keyValue = NULL;
    char *keyWord = NULL;

    char *tempString = NULL;

    char **groupName = NULL;
    int  numGroupName = 0;

    int tempInt=0; // Temporary integer

    // Destroy our analysis structures coming in if aren't 0 and NULL already
    if (feaProblem->feaAnalysis != NULL) {
        for (i = 0; i < feaProblem->numAnalysis; i++) {
            status = destroy_feaAnalysisStruct(&feaProblem->feaAnalysis[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }

    if (feaProblem->feaAnalysis != NULL) EG_free(feaProblem->feaAnalysis);
    feaProblem->feaAnalysis = NULL;
    feaProblem->numAnalysis = 0;

    printf("\nGetting FEA analyses.......\n");

    feaProblem->numAnalysis = numAnalysisTuple;

    printf("\tNumber of analyses - %d\n", feaProblem->numAnalysis);

    if (feaProblem->numAnalysis > 0) {
        feaProblem->feaAnalysis = (feaAnalysisStruct *) EG_alloc(feaProblem->numAnalysis * sizeof(feaAnalysisStruct));
    } else {
        printf("\tNumber of analysis values in input tuple is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < feaProblem->numAnalysis; i++) {
        status = initiate_feaAnalysisStruct(&feaProblem->feaAnalysis[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    for (i = 0; i < feaProblem->numAnalysis; i++) {

        printf("\tAnalysis name - %s\n", analysisTuple[i].name);

        feaProblem->feaAnalysis[i].name = (char *) EG_alloc(((strlen(analysisTuple[i].name)) + 1)*sizeof(char));
        if (feaProblem->feaAnalysis[i].name == NULL) return EGADS_MALLOC;

        memcpy(feaProblem->feaAnalysis[i].name, analysisTuple[i].name, strlen(analysisTuple[i].name)*sizeof(char));
        feaProblem->feaAnalysis[i].name[strlen(analysisTuple[i].name)] = '\0';

        feaProblem->feaAnalysis[i].analysisID = i + 1;

        // Do we have a json string?
        if (strncmp(analysisTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", analysisTuple[i].value);

            /*! \page feaAnalysis
             * \section jsonStringAnalysis JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             *  (e.g. "Value" = {"numDesiredEigenvalue": 10, "eigenNormalization": "MASS", "numEstEigenvalue": 1,
             * "extractionMethod": "GIV", "frequencyRange": [0, 10000]})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * \if MYSTRAN
             * <ul>
             *  <li> <B>analysisType = "Modal"</B> </li> <br>
             *  Type of load. Options: "Modal", "Static".
             * </ul>
             * \elseif (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>analysisType = "Modal"</B> </li> <br>
             *  Type of load. Options: "Modal", "Static", "AeroelasticTrim", "AeroelasticFlutter"
             *  Note: "AeroelasticStatic" is still supported but refers to "AeroelasticTrim"
             *  Note: "Optimization" and "StaticOpt" are not valid - Optimization is initialized by the Analysis_Type AIM Input
             * </ul>
             * \elseif ABAQUS
             * <ul>
             *  <li> <B>analysisType = "Modal"</B> </li> <br>
             *  Type of load. Options: "Modal", "Static"
             * </ul>
             * \endif
             *
             */
            // Get analysis Type
            keyWord = "analysisType";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                //{UnknownAnalysis, Modal, Static}
                if      (strcasecmp(keyValue, "\"Modal\"")  == 0) feaProblem->feaAnalysis[i].analysisType = Modal;
                else if (strcasecmp(keyValue, "\"Static\"") == 0) feaProblem->feaAnalysis[i].analysisType = Static;
                else if (strcasecmp(keyValue, "\"StaticOpt\"") == 0) feaProblem->feaAnalysis[i].analysisType = Optimization;
                else if (strcasecmp(keyValue, "\"Optimization\"") == 0) feaProblem->feaAnalysis[i].analysisType = Optimization;
                else if (strcasecmp(keyValue, "\"AeroelasticTrim\"") == 0) feaProblem->feaAnalysis[i].analysisType = AeroelasticTrim;
                else if (strcasecmp(keyValue, "\"AeroelasticStatic\"") == 0) feaProblem->feaAnalysis[i].analysisType = AeroelasticTrim;
                else if (strcasecmp(keyValue, "\"AeroelasticFlutter\"") == 0) feaProblem->feaAnalysis[i].analysisType = AeroelasticFlutter;
                else if (strcasecmp(keyValue, "\"Gust\"") == 0) feaProblem->feaAnalysis[i].analysisType = Gust;
                else {

                    printf("\tUnrecognized \"%s\" specified (%s) for Analysis tuple %s, defaulting to \"Modal\"\n", keyWord,
                                                                                                                    keyValue,
                                                                                                                    analysisTuple[i].name);
                    feaProblem->feaAnalysis[i].analysisType = Modal;
                }

            } else {

                printf("\tNo \"%s\" specified for Analysis tuple %s, defaulting to \"Modal\"\n", keyWord, analysisTuple[i].name);

                feaProblem->feaAnalysis[i].analysisType = Modal;
            }
            AIM_FREE(keyValue);

            // Get loads to be applied for a given analysis

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             * <ul>
             *  <li> <B>analysisLoad = "(no default)"</B> </li> <br>
             *  Single or list of "Load Name"s defined in \ref feaLoad in which to use for the analysis (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             */
            keyWord = "analysisLoad";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaAnalysis[i].numLoad = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numLoad; attrIndex++) {

                        if (strcasecmp(feaProblem->feaLoad[attrIndex].name, groupName[groupIndex]) == 0) {

//                            if (feaProblem->feaLoad[attrIndex].loadType == Thermal) {
//                                printf("Combining Thermal loads in a subcase isn't supported yet!\n");
//                                status = string_freeArray(numGroupName, &groupName);
//                                if (status != CAPS_SUCCESS) return status;
//                                groupName = NULL;
//                                return CAPS_BADVALUE;
//                            }

                            feaProblem->feaAnalysis[i].numLoad += 1;

                            if (feaProblem->feaAnalysis[i].numLoad == 1)  {

                                feaProblem->feaAnalysis[i].loadSetID = (int *) EG_alloc(feaProblem->feaAnalysis[i].numLoad *sizeof(int));

                            } else {

                                feaProblem->feaAnalysis[i].loadSetID = (int *) EG_reall(feaProblem->feaAnalysis[i].loadSetID,
                                                                               feaProblem->feaAnalysis[i].numLoad *sizeof(int));
                            }

                            if (feaProblem->feaAnalysis[i].loadSetID == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaAnalysis[i].loadSetID[feaProblem->feaAnalysis[i].numLoad-1] = feaProblem->feaLoad[attrIndex].loadID;
                            break;
                        }
                    }

                    if (feaProblem->feaAnalysis[i].numLoad != groupIndex+1) {

                        printf("\tWarning: Analysis load name, %s, not found in feaLoad structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;
                groupName = NULL;

            }

            // Get constraints to be applied for a given analysis
            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             * <ul>
             *  <li> <B>analysisConstraint = "(no default)"</B> </li> <br>
             *  Single or list of "Constraint Name"s defined in \ref feaConstraint in which to use for the analysis (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             */
            keyWord = "analysisConstraint";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaAnalysis[i].numConstraint = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numConstraint; attrIndex++) {

                        if (strcasecmp(feaProblem->feaConstraint[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaAnalysis[i].numConstraint += 1;

                            if (feaProblem->feaAnalysis[i].numConstraint == 1)  {

                                feaProblem->feaAnalysis[i].constraintSetID = (int *) EG_alloc(feaProblem->feaAnalysis[i].numConstraint *sizeof(int));

                            } else {

                                feaProblem->feaAnalysis[i].constraintSetID = (int *) EG_reall(feaProblem->feaAnalysis[i].constraintSetID,
                                                                                              feaProblem->feaAnalysis[i].numConstraint *sizeof(int));
                            }

                            if (feaProblem->feaAnalysis[i].constraintSetID == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaAnalysis[i].constraintSetID[feaProblem->feaAnalysis[i].numConstraint-1] = feaProblem->feaConstraint[attrIndex].constraintID;
                            break;
                        }
                    }

                    if (feaProblem->feaAnalysis[i].numConstraint != groupIndex+1) {

                        printf("\tWarning: Analysis constraint name, %s, not found in feaConstraint structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;
                groupName = NULL;

            }

            // Get supports to be applied for a given analysis
            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>analysisSupport = "(no default)"</B> </li> <br>
             *  Single or list of "Support Name"s defined in \ref feaSupport in which to use for the analysis (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             *
             */
            keyWord = "analysisSupport";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaAnalysis[i].numSupport = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numSupport; attrIndex++) {

                        if (strcasecmp(feaProblem->feaSupport[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaAnalysis[i].numSupport += 1;

                            if (feaProblem->feaAnalysis[i].numSupport == 1)  {

                                feaProblem->feaAnalysis[i].supportSetID = (int *) EG_alloc(feaProblem->feaAnalysis[i].numSupport *sizeof(int));

                            } else {

                                feaProblem->feaAnalysis[i].supportSetID = (int *) EG_reall(feaProblem->feaAnalysis[i].supportSetID,
                                                                                           feaProblem->feaAnalysis[i].numSupport *sizeof(int));
                            }

                            if (feaProblem->feaAnalysis[i].supportSetID == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaAnalysis[i].supportSetID[feaProblem->feaAnalysis[i].numSupport-1] = feaProblem->feaSupport[attrIndex].supportID;
                            break;
                        }
                    }

                    if (feaProblem->feaAnalysis[i].numSupport != groupIndex+1) {

                        printf("\tWarning: Analysis support name, %s, not found in feaSupport structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;
                groupName = NULL;

            }

            // Get design constraints to be applied for a given analysis

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>analysisDesignConstraint = "(no default)"</B> </li> <br>
             *  Single or list of "Design Constraint Name"s defined in \ref feaDesignConstraint in which to use for the analysis (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             */
            keyWord = "analysisDesignConstraint";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaAnalysis[i].numDesignConstraint = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numDesignConstraint; attrIndex++) {

                        if (strcasecmp(feaProblem->feaDesignConstraint[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaAnalysis[i].numDesignConstraint += 1;

                            if (feaProblem->feaAnalysis[i].numDesignConstraint == 1)  {

                                feaProblem->feaAnalysis[i].designConstraintSetID = (int *) EG_alloc(feaProblem->feaAnalysis[i].numDesignConstraint *sizeof(int));

                            } else {

                                feaProblem->feaAnalysis[i].designConstraintSetID = (int *) EG_reall(feaProblem->feaAnalysis[i].designConstraintSetID,
                                                                                                    feaProblem->feaAnalysis[i].numDesignConstraint *sizeof(int));
                            }

                            if (feaProblem->feaAnalysis[i].designConstraintSetID == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaAnalysis[i].designConstraintSetID[feaProblem->feaAnalysis[i].numDesignConstraint-1] = feaProblem->feaDesignConstraint[attrIndex].designConstraintID;
                            break;
                        }
                    }

                    if (feaProblem->feaAnalysis[i].numDesignConstraint != groupIndex+1) {

                        printf("\tWarning: Analysis design constraint name, %s, not found in feaDesignConstraint structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) {
                    if (tempString != NULL) EG_free(tempString);

                    return status;
                }
                groupName = NULL;

            }

            // Get mass increment for a concentrated mass to be applied for a given analysis subcase

            /*! \page feaAnalysis
             *
             * \if (NASTRAN)
             * <ul>
             *  <li> <B>analysisCaseMassIncrement = "(no default)"</B> </li> <br>
             *  Single or list of "Concentrated Mass Property Name"s defined in in which to use for the analysis (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             */
            keyWord = "analysisMassIncrement";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaAnalysis[i].numMassIncrement = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numMassIncrement; attrIndex++) {

//                        printf("%s : %s\n", feaProblem->feaMassIncrement[attrIndex].name, groupName[groupIndex]);

                        if (strcasecmp(feaProblem->feaMassIncrement[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaAnalysis[i].numMassIncrement += 1;

                            if (feaProblem->feaAnalysis[i].numMassIncrement == 1)  {

                                feaProblem->feaAnalysis[i].massIncrementSetID = (int *) EG_alloc(feaProblem->feaAnalysis[i].numMassIncrement *sizeof(int));

                            } else {

                                feaProblem->feaAnalysis[i].massIncrementSetID = (int *) EG_reall(feaProblem->feaAnalysis[i].massIncrementSetID,
                                                                                                 feaProblem->feaAnalysis[i].numMassIncrement *sizeof(int));
                            }

                            if (feaProblem->feaAnalysis[i].massIncrementSetID == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaAnalysis[i].massIncrementSetID[feaProblem->feaAnalysis[i].numMassIncrement-1] = feaProblem->feaMassIncrement[attrIndex].massIncrementID;
                            break;
                        }
                    }

                    if (feaProblem->feaAnalysis[i].numMassIncrement != groupIndex+1) {

                        printf("\tWarning: Analysis mass increment name, %s, not found in feaProperty structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) {
                    if (tempString != NULL) EG_free(tempString);

                    return status;
                }
                groupName = NULL;

            }

            //Fill up analysis properties
            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             * <ul>
             *  <li> <B>extractionMethod = "(no default)"</B> </li> <br>
             *  Extraction method for modal analysis.
             * </ul>
             * \endif
             */
            keyWord = "extractionMethod";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);
                AIM_FREE(keyValue);

                feaProblem->feaAnalysis[i].extractionMethod = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaAnalysis[i].extractionMethod == NULL) {
                    AIM_FREE(tempString);
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                memcpy(feaProblem->feaAnalysis[i].extractionMethod, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaAnalysis[i].extractionMethod[strlen(tempString)] = '\0';

                AIM_FREE(tempString);
            }


            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             * <ul>
             *  <li> <B>frequencyRange = [0.0, 0.0] </B> </li> <br>
             *  Frequency range of interest for modal analysis.
             * </ul>
             * \endif
             */
            keyWord = "frequencyRange";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleArray(keyValue,
                                              (int) sizeof(feaProblem->feaAnalysis[i].frequencyRange)/sizeof(double),
                                              feaProblem->feaAnalysis[i].frequencyRange);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>numEstEigenvalue = 0 </B> </li> <br>
             *  Number of estimated eigenvalues for modal analysis.
             * </ul>
             * \endif
             */
            keyWord = "numEstEigenvalue";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaAnalysis[i].numEstEigenvalue);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS || ABAQUS)
             * <ul>
             *  <li> <B>numDesiredEigenvalue = 0</B> </li> <br>
             *  Number of desired eigenvalues for modal analysis.
             * </ul>
             * \endif
             *
             */
            keyWord = "numDesiredEigenvalue";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaAnalysis[i].numDesiredEigenvalue);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>eigenNormalization = "(no default)"</B> </li> <br>
             *  Method of eigenvector renormalization. Options: "POINT", "MAX", "MASS"
             * </ul>
             * \elseif ABAQUS
             * <ul>
             *  <li> <B>eigenNormalization = "(no default)"</B> </li> <br>
             *  Method of eigenvector renormalization. Options: "DISPLACEMENT", "MASS"
             * </ul>
             * \endif
             */
            keyWord = "eigenNormalization";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                AIM_STRDUP(feaProblem->feaAnalysis[i].eigenNormalization, tempString, aimInfo, status);

                AIM_FREE(keyValue);
                AIM_FREE(tempString);
            }
            // Look for old misspelling and ask user to fix it.
            keyWord = "eigenNormaliztion";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
              AIM_FREE(keyValue);
              AIM_ERROR(aimInfo, "The incorrectly spelled 'eigenNormaliztion' (missing 'a') has been been fixed to 'eigenNormalization'. Please correct your input!");
              status = CAPS_BADVALUE;
              goto cleanup;
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>gridNormalization = 0 </B> </li> <br>
             *  Grid point to be used in normalizing eigenvector to 1.0 when using eigenNormalization = "POINT"
             * </ul>
             * \endif
             */
            keyWord = "gridNormalization";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaAnalysis[i].gridNormalization);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>componentNormalization = 0</B> </li> <br>
             *  Degree of freedom about "gridNormalization" to be used in normalizing eigenvector to 1.0
             *  when using eigenNormalization = "POINT"
             * </ul>
             * \endif
             *
             */
            keyWord = "componentNormalization";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaAnalysis[i].componentNormalization);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>lanczosMode = 2</B> </li> <br>
             *  Mode refers to the Lanczos mode type to be used in the solution. In mode 3 the mass matrix, Maa,must
             *  be nonsingular whereas in mode 2 the matrix K aa - sigma*Maa must be nonsingular
             *  </ul>
             * \endif
             */
            keyWord = "lanczosMode";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaAnalysis[i].lanczosMode);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN)
             * <ul>
             *  <li> <B>lanczosType = "(no default)"</B> </li> <br>
             *  Lanczos matrix type. Options: DPB, DGB.
             *  </ul>
             * \endif
             */
            keyWord = "lanczosType";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                feaProblem->feaAnalysis[i].lanczosType = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaAnalysis[i].lanczosType == NULL) {
                    if (tempString != NULL) EG_free(tempString);

                    return EGADS_MALLOC;
                }

                memcpy(feaProblem->feaAnalysis[i].lanczosType, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaAnalysis[i].lanczosType[strlen(tempString)] = '\0';

                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                AIM_FREE(tempString);
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN)
             * <ul>
             *  <li> <B>machNumber = 0.0 or [0.0, ..., 0.0]</B> </li> <br>
             *  Mach number used in trim analysis OR Mach numbers used in flutter analysis..
             *  </ul>
             * \elseif (ASTROS)
             * <ul>
             *  <li> <B>machNumber = 0.0 or [0.0, ..., 0.0]</B> </li> <br>
             *  Mach number used in trim analysis OR Mach up to 6 values used in flutter analysis..
             *  </ul>
             * \endif
             */
            keyWord = "machNumber";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleDynamicArray(keyValue, &feaProblem->feaAnalysis[i].numMachNumber, &feaProblem->feaAnalysis[i].machNumber);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>dynamicPressure = 0.0</B> </li> <br>
             *  Dynamic pressure used in trim analysis.
             *  </ul>
             * \endif
             */
            keyWord = "dynamicPressure";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaAnalysis[i].dynamicPressure);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>density = 0.0</B> </li> <br>
             *  Density used in trim analysis to determine true velocity, or flutter analysis.
             *  </ul>
             * \endif
             */
            keyWord = "density";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaAnalysis[i].density);
                AIM_FREE(keyValue);
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>aeroSymmetryXY = "(no default)"</B> </li> <br>
             *  Aerodynamic symmetry about the XY Plane. Options: SYM, ANTISYM, ASYM.
             *  Aerodynamic symmetry about the XY Plane. Options: SYM, ANTISYM, ASYM. SYMMETRIC Indicates that a half span aerodynamic model
             *  is moving in a symmetric manner with respect to the XY plane.
             *  ANTISYMMETRIC Indicates that a half span aerodynamic model is moving in an antisymmetric manner with respect to the XY plane.
             *  ASYMMETRIC Indicates that a full aerodynamic model is provided.
             *  </ul>
             * \endif
             */
            keyWord = "aeroSymmetryXY";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                feaProblem->feaAnalysis[i].aeroSymmetryXY = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaAnalysis[i].aeroSymmetryXY == NULL) {
                    if (tempString != NULL) EG_free(tempString);

                    return EGADS_MALLOC;
                }

                memcpy(feaProblem->feaAnalysis[i].aeroSymmetryXY, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaAnalysis[i].aeroSymmetryXY[strlen(tempString)] = '\0';

                AIM_FREE(keyValue);
                AIM_FREE(tempString);
            }

            // check for the old option trimSymmetry
            keyWord = "trimSymmetry";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                feaProblem->feaAnalysis[i].aeroSymmetryXY = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaAnalysis[i].aeroSymmetryXY == NULL) {
                    if (tempString != NULL) EG_free(tempString);

                    return EGADS_MALLOC;
                }

                memcpy(feaProblem->feaAnalysis[i].aeroSymmetryXY, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaAnalysis[i].aeroSymmetryXY[strlen(tempString)] = '\0';

                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                AIM_FREE(tempString);
            }

            /*! \page feaAnalysis
             *
             * \if (MYSTRAN || NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>aeroSymmetryXZ = "(no default)"</B> </li> <br>
             *  Aerodynamic symmetry about the XZ Plane. Options: SYM, ANTISYM, ASYM. SYMMETRIC Indicates that a half span aerodynamic model
             *  is moving in a symmetric manner with respect to the XZ plane.
             *  ANTISYMMETRIC Indicates that a half span aerodynamic model is moving in an antisymmetric manner with respect to the XZ plane.
             *  ASYMMETRIC Indicates that a full aerodynamic model is provided.
             *  </ul>
             * \endif
             */
            keyWord = "aeroSymmetryXZ";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                feaProblem->feaAnalysis[i].aeroSymmetryXZ = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaAnalysis[i].aeroSymmetryXZ == NULL) {
                    if (tempString != NULL) EG_free(tempString);

                    return EGADS_MALLOC;
                }

                memcpy(feaProblem->feaAnalysis[i].aeroSymmetryXZ, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaAnalysis[i].aeroSymmetryXZ[strlen(tempString)] = '\0';

                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                AIM_FREE(tempString);
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN )
             * <ul>
             *  <li> <B>rigidVariable = ["no default"]</B> </li> <br>
             *  List of rigid body motions to be used as trim variables during a trim analysis. Nastran
             *  valid labels are: ANGLEA, SIDES, ROLL, PITCH, YAW, URDD1, URDD2, URDD3, URDD4, URDD5, URDD6
             *  </ul>
             * \elseif ASTROS
             * <ul>
             *  <li> <B>rigidVariable = ["no default"]</B> </li> <br>
             *  List of rigid body motions to be used as trim variables during a trim analysis. Nastran format
             *  labels are used and will be converted by the AIM automatically.
             *  Expected inputs: ANGLEA, SIDES, ROLL, PITCH, YAW, URDD1, URDD2, URDD3, URDD4, URDD5, URDD6
             *  </ul>
             * \endif
             */
            keyWord = "rigidVariable";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue,
                        &feaProblem->feaAnalysis[i].numRigidVariable,
                        &feaProblem->feaAnalysis[i].rigidVariable);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN)
             * <ul>
             *  <li> <B>rigidConstraint = ["no default"]</B> </li> <br>
             *  List of rigid body motions to be used as trim constraint variables during a trim analysis. Nastran
             *  valid labels are: ANGLEA, SIDES, ROLL, PITCH, YAW, URDD1, URDD2, URDD3, URDD4, URDD5, URDD6
             *  </ul>
             * \elseif ASTROS
             * <ul>
             *  <li> <B>rigidConstraint = ["no default"]</B> </li> <br>
             *  List of rigid body motions to be used as trim constraint variables during a trim analysis. Nastran format
             *  labels are used and will be converted by the AIM automatically.
             *  Expected inputs: ANGLEA, SIDES, ROLL, PITCH, YAW, URDD1, URDD2, URDD3, URDD4, URDD5, URDD6
             * </ul>
             * \endif
             */
            keyWord = "rigidConstraint";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue,
                        &feaProblem->feaAnalysis[i].numRigidConstraint,
                        &feaProblem->feaAnalysis[i].rigidConstraint);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            if (feaProblem->feaAnalysis[i].numRigidConstraint != 0) {
                /*! \page feaAnalysis
                 *
                 * \if (NASTRAN || ASTROS)
                 * <ul>
                 *  <li> <B>magRigidConstraint = [0.0 , 0.0, ...]</B> </li> <br>
                 *  List of magnitudes of trim constraint variables. If none and 'rigidConstraint'(s) are specified
                 *  then 0.0 is assumed for each rigid constraint.
                 *  </ul>
                 * \endif
                 */
                keyWord = "magRigidConstraint";
                status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
                if (status == CAPS_SUCCESS) {

                    status = string_toDoubleDynamicArray(keyValue,
                            &tempInt,
                            &feaProblem->feaAnalysis[i].magRigidConstraint);
                    if (keyValue != NULL) {
                        EG_free(keyValue);
                        keyValue = NULL;
                    }
                    if (status != CAPS_SUCCESS) return status;
                } else {

                    tempInt = feaProblem->feaAnalysis[i].numRigidConstraint;

                    feaProblem->feaAnalysis[i].magRigidConstraint = (double *) EG_alloc(tempInt*sizeof(double));
                    if (feaProblem->feaAnalysis[i].magRigidConstraint == NULL) return EGADS_MALLOC;

                    for (j = 0; j < feaProblem->feaAnalysis[i].numRigidConstraint; j++) {

                        feaProblem->feaAnalysis[i].magRigidConstraint[j] = 0.0;
                    }
                }

                if (tempInt != feaProblem->feaAnalysis[i].numRigidConstraint) {
                    printf("\tDimensional mismatch between 'magRigidConstraint' and 'rigidConstraint'.\n");
                    printf("\t 'magRigidConstraint' will be resized.\n");

                    feaProblem->feaAnalysis[i].magRigidConstraint = (double *) EG_reall(feaProblem->feaAnalysis[i].magRigidConstraint,
                            feaProblem->feaAnalysis[i].numRigidConstraint*sizeof(double));

                    if (feaProblem->feaAnalysis[i].magRigidConstraint == NULL) return EGADS_MALLOC;

                    for (j = tempInt; j < feaProblem->feaAnalysis[i].numRigidConstraint; j++) {
                        feaProblem->feaAnalysis[i].magRigidConstraint[j] = 0.0;
                    }
                }
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>controlConstraint = ["no default"]</B> </li> <br>
             *  List of controls surfaces to be used as trim constraint variables during a trim analysis.
             *  </ul>
             * \endif
             */
            keyWord = "controlConstraint";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue,
                        &feaProblem->feaAnalysis[i].numControlConstraint,
                        &feaProblem->feaAnalysis[i].controlConstraint);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            if (feaProblem->feaAnalysis[i].numControlConstraint != 0) {
                /*! \page feaAnalysis
                 *
                 * \if (NASTRAN || ASTROS)
                 * <ul>
                 *  <li> <B>magControlConstraint = [0.0 , 0.0, ...]</B> </li> <br>
                 *  List of magnitudes of trim control surface constraint variables. If none and 'controlConstraint'(s) are specified
                 *   then 0.0 is assumed for each control surface constraint.
                 *  </ul>
                 * \endif
                 */
                keyWord = "magControlConstraint";
                status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
                if (status == CAPS_SUCCESS) {

                    status = string_toDoubleDynamicArray(keyValue,
                            &tempInt,
                            &feaProblem->feaAnalysis[i].magControlConstraint);
                    if (keyValue != NULL) {
                        EG_free(keyValue);
                        keyValue = NULL;
                    }
                    if (status != CAPS_SUCCESS) return status;
                } else {

                    tempInt = feaProblem->feaAnalysis[i].numControlConstraint;

                    feaProblem->feaAnalysis[i].magControlConstraint = (double *) EG_alloc(tempInt*sizeof(double));
                    if (feaProblem->feaAnalysis[i].magControlConstraint == NULL) return EGADS_MALLOC;

                    for (j = 0; j < feaProblem->feaAnalysis[i].numControlConstraint; j++) {

                        feaProblem->feaAnalysis[i].magControlConstraint[j] = 0.0;
                    }
                }

                if (tempInt != feaProblem->feaAnalysis[i].numControlConstraint) {
                    printf("\tDimensional mismatch between 'magControlConstraint' and 'controlConstraint'.\n");
                    printf("\t 'magControlConstraint' will be resized.\n");

                    feaProblem->feaAnalysis[i].magControlConstraint = (double *) EG_reall(feaProblem->feaAnalysis[i].magControlConstraint,
                            feaProblem->feaAnalysis[i].numControlConstraint*sizeof(double));

                    if (feaProblem->feaAnalysis[i].magControlConstraint == NULL) return EGADS_MALLOC;

                    for (j = tempInt; j < feaProblem->feaAnalysis[i].numControlConstraint; j++) {
                        feaProblem->feaAnalysis[i].magControlConstraint[j] = 0.0;
                    }
                }
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>reducedFreq = [0.1, ..., 20.0], No Default Values are defined.</B> </li> <br>
             *  Reduced Frequencies to be used in Flutter Analysis.  Up to 8 values can be defined.
             *  </ul>
             * \endif
             */
            keyWord = "reducedFreq";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleDynamicArray(keyValue, &feaProblem->feaAnalysis[i].numReducedFreq, &feaProblem->feaAnalysis[i].reducedFreq);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (feaProblem->feaAnalysis[i].numReducedFreq > 8) {
                    printf("\tError: The number of reduced frequencies (reducedFreq) entered in an Analysis AIM Input must be eight or less\n");
                    return CAPS_BADVALUE;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>flutterVel = [0.1, ..., 20.0]</B> </li> <br>
             *  Velocities to be used in Flutter Analysis. If no values are provided the
             *  following relation is used
             *
             *   v = sqrt(2*dynamicPressure/density)
             *   dv = (v*2 - v/2) / 20;
             *
             *   flutterVel[0]  = v/10
             *   flutterVel[i] = v/2 + i*dv; where i = 1....21
             *   flutterVel[22] = v*10;
             *  </ul>
             * \endif
             */
            keyWord = "flutterVel";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleDynamicArray(keyValue, &feaProblem->feaAnalysis[i].numFlutterVel, &feaProblem->feaAnalysis[i].flutterVel);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaAnalysis
             *
             * \if (ASTROS)
             * <ul>
             *  <li> <B>flutterConvergenceParam = 1e-5 </B> </li> <br>
             *  Convergence parameter for flutter eigenvalue.
             *  </ul>
             * \endif
             */
            keyWord = "flutterConvergenceParam";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaAnalysis[i].flutterConvergenceParam);
                AIM_FREE(keyValue);
                if (status != CAPS_SUCCESS) return status;
            }
            else {
                feaProblem->feaAnalysis[i].flutterConvergenceParam = 1.e-5;
            }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>visualFlutter = False </B> </li> <br>
             *  Turn on flutter visualization f06 output.
             *  </ul>
             * \endif
             */
            // keyWord = "visualFlutter";
            // status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            // if (status == CAPS_SUCCESS) {
            //
            //     status = string_toBoolean(keyValue, &feaProblem->feaAnalysis[i].visualFlutter);
            //     if (keyValue != NULL) {
            //         EG_free(keyValue);
            //         keyValue = NULL;
            //     }
            //     if (status != CAPS_SUCCESS) return status;
            // }

            /*! \page feaAnalysis
             *
             * \if (NASTRAN)
             * <ul>
             *  <li> <B>analysisResponse = "(no default)"</B> </li> <br>
             *  Single or list of "DesignResponse Name"s defined in \ref feaDesignResponse to use for the analysis response spanning sets (e.g. "Name1" or ["Name1","Name2",...].
             * </ul>
             * \endif
             */
            keyWord = "analysisResponse";
            status = search_jsonDictionary( analysisTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaAnalysis[i].numDesignResponse = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numDesignResponse; attrIndex++) {

                        if (strcasecmp(feaProblem->feaDesignResponse[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaAnalysis[i].numDesignResponse += 1;

                            if (feaProblem->feaAnalysis[i].numDesignResponse == 1)  {

                                feaProblem->feaAnalysis[i].designResponseSetID = (int *) EG_alloc(feaProblem->feaAnalysis[i].numDesignResponse *sizeof(int));

                            } else {

                                feaProblem->feaAnalysis[i].designResponseSetID = (int *) EG_reall(feaProblem->feaAnalysis[i].designResponseSetID,
                                                                                                    feaProblem->feaAnalysis[i].numDesignResponse *sizeof(int));
                            }

                            if (feaProblem->feaAnalysis[i].designResponseSetID == NULL) {
                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaAnalysis[i].designResponseSetID[feaProblem->feaAnalysis[i].numDesignResponse-1] = feaProblem->feaDesignResponse[attrIndex].responseID;
                            break;
                        }
                    }

                    if (feaProblem->feaAnalysis[i].numDesignResponse != groupIndex+1) {

                        printf("\tWarning: Analysis design response name, %s, not found in feaDesignResponse structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) {
                    if (tempString != NULL) EG_free(tempString);

                    return status;
                }
                groupName = NULL;
            }

            if (feaProblem->feaAnalysis[i].analysisType == AeroelasticFlutter) {
                // Setup the default flutter velocities  if not specified
                status = fea_defaultFlutterVelocity(&feaProblem->feaAnalysis[i]);
                if (status != CAPS_SUCCESS) return status;
            }

        } else { // Not JSONstring

            /*! \page feaAnalysis
             * \section keyStringAnalysis Single Value String
             *
             * If "Value" is a string, the string value may correspond to an entry in a predefined analysis lookup
             * table. NOT YET IMPLEMENTED!!!!
             *
             *
             */

            // CALL analysis look up
            printf("\tError: Analysis tuple value is expected to be a JSON string\n");
            //printf("\t%s", analysisTuple[i].value);
            return CAPS_BADVALUE;
        }
    }


    printf("\tDone getting FEA analyses\n");
    status = CAPS_SUCCESS;
cleanup:
    AIM_FREE(keyValue);
    AIM_FREE(tempString);
    return status;
}

// Get the design variables from a capsTuple
int fea_getDesignVariable(void *aimInfo,
             /*@unused@*/ int requireGroup,
                          int numDesignVariableTuple,
                          capsTuple designVariableTuple[],
                          int numDesignVariableRelationTuple,
                          capsTuple designVariableRelationTuple[],
                          mapAttrToIndexStruct *attrMap,
                          feaProblemStruct *feaProblem) {

    /*! \page feaDesignVariable FEA Design Variables
     * Structure for the design variable tuple  = ("DesignVariable Name", "Value").
     * "DesignVariable Name" defines the reference name for the design variable being specified.
     *  This string will be used in the FEA input directly. The "Value" must be a JSON String dictionary
     *  (see Section \ref jsonStringDesignVariable).
     *  \if NASTRAN
     *  In Nastran the DesignVariable Name will be the LABEL used in the DESVAR input.
     *  For this reason the user should keep the length of this input to a minimum number of characters, ideally 7 or less.
     *
     *  - <c>DESVAR  ID      LABEL   XINIT   XLB     XUB     DELXV   DDVAL</c>
     *  \endif
     *
     */
    int status; //Function return

    int i, j, k; // Indexing
    int found;

    char *keyValue = NULL;
    char *keyWord = NULL;

    int tempInteger = 0;
    int designVariableID = 0;

    char **groupName = NULL;
    int  numGroupName = 0;

    feaDesignVariableStruct *designVariable;
    feaDesignVariableRelationStruct *designVariableRelation;

    // Destroy our design variable structures coming in if aren't 0 and NULL already
    if (feaProblem->feaDesignVariable != NULL) {
        for (i = 0; i < feaProblem->numDesignVariable; i++) {
            status = destroy_feaDesignVariableStruct(&feaProblem->feaDesignVariable[i]);
            AIM_STATUS(aimInfo, status);
        }
    }
    AIM_FREE(feaProblem->feaDesignVariable);
    feaProblem->feaDesignVariable = NULL;
    feaProblem->numDesignVariable = 0;

    // Destroy our design variable relation structures coming in if aren't 0 and NULL already
    if (feaProblem->feaDesignVariableRelation != NULL) {
        for (i = 0; i < feaProblem->numDesignVariableRelation; i++) {
            status = destroy_feaDesignVariableRelationStruct(&feaProblem->feaDesignVariableRelation[i]);
            AIM_STATUS(aimInfo, status);
        }
    }
    AIM_FREE(feaProblem->feaDesignVariableRelation);
    feaProblem->feaDesignVariableRelation = NULL;
    feaProblem->numDesignVariableRelation = 0;


    printf("\nGetting FEA design variables.......\n");

    feaProblem->numDesignVariable = numDesignVariableTuple;

    printf("\tNumber of design variables          - %d\n", feaProblem->numDesignVariable);

    if (feaProblem->numDesignVariable > 0) {
        AIM_ALLOC(feaProblem->feaDesignVariable,
                  feaProblem->numDesignVariable, feaDesignVariableStruct, aimInfo, status);

    } else {
        AIM_ERROR(aimInfo, "Number of design variable values in input tuple is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < feaProblem->numDesignVariable; i++) {
        status = initiate_feaDesignVariableStruct(&feaProblem->feaDesignVariable[i]);
        AIM_STATUS(aimInfo, status);
    }

    for (i = 0; i < feaProblem->numDesignVariable; i++) {

        designVariable = &feaProblem->feaDesignVariable[i];

        printf("\tDesign_Variable name - %s\n", designVariableTuple[i].name);

        AIM_STRDUP(designVariable->name, designVariableTuple[i].name, aimInfo, status);

        // Do we have a json string?
        if (json_isDict(designVariableTuple[i].value)) {

            /*! \page feaDesignVariable
             * \section jsonStringDesignVariable JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (NASTRAN || ASTROS || MASSTRAN)
             *  (eg. "Value" = {"initialValue": 5.0, "upperBound": 10.0})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             */

            /*! \page feaDesignVariable
             *
             * \if (NASTRAN || ASTROS || MASSTRAN)
             * <ul>
             *  <li> <B>initialValue = 0.0</B> </li> <br>
             *  Initial value for the design variable.
             * </ul>
             * \endif
             *
             */
            keyWord = "initialValue";
            status = json_getDouble(
                designVariableTuple[i].value, keyWord,
                &designVariable->initialValue);

            if (status != CAPS_SUCCESS) {
                designVariable->initialValue = 0.0;
            }

            /*! \page feaDesignVariable
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>lowerBound = 0.0</B> </li> <br>
             *  Lower bound for the design variable.
             * </ul>
             * \endif
             */
            keyWord = "lowerBound";
            status = json_getDouble(
                designVariableTuple[i].value, keyWord,
                &designVariable->lowerBound);

            if (status != CAPS_SUCCESS) {
                designVariable->lowerBound = 0.0;
            }

            /*! \page feaDesignVariable
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>upperBound = 0.0</B> </li> <br>
             *  Upper bound for the design variable.
             * </ul>
             * \endif
             */
            keyWord = "upperBound";
            status = json_getDouble(
                designVariableTuple[i].value, keyWord,
                &designVariable->upperBound);

            if (status != CAPS_SUCCESS) {
                designVariable->upperBound = 0.0;
            }

            /*! \page feaDesignVariable
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>maxDelta = 0.5</B> </li> <br>
             *  Change fraction limit for the design variable.
             * </ul>
             * \endif
             */
            keyWord = "maxDelta";
            status = json_getDouble(
                designVariableTuple[i].value, keyWord,
                &designVariable->maxDelta);

            if (status != CAPS_SUCCESS) {
                designVariable->maxDelta = 0.5;
            }

            /*! \page feaDesignVariable
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>discreteValue = 0.0</B> </li> <br>
             *  List of discrete values to use for the design variable (e.g. [0.0,1.0,1.5,3.0].
             * </ul>
             * \endif
             */
            keyWord = "discreteValue";
            status = json_getDoubleDynamicArray(
                designVariableTuple[i].value, keyWord,
                &designVariable->numDiscreteValue, &designVariable->discreteValue);

            if (status != CAPS_SUCCESS) {
                designVariable->numDiscreteValue = 0;
                designVariable->discreteValue = NULL;
            }

            /*! \page feaDesignVariable
             *
             * \if NASTRAN
             * <ul>
             *  <li> <B>independentVariable = "(no default)"</B> </li> <br>
             *  Single or list of "DesignVariable Name"s  (that is the Tuple name) used to create/designate a
             *  dependent design variable.
             *  - independentValue = variableWeight[1] + variableWeight[2] * SUM{independentVariableWeight[i] * independentVariable[i]}
             * </ul>
             * \endif
             */

            keyWord = "independentVariable";
            status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue,
                                                     &designVariable->numIndependVariable,
                                                     &designVariable->independVariable);
                AIM_FREE(keyValue);
                AIM_STATUS(aimInfo, status);

            }

            /*! \page feaDesignVariable
             *
             * \if NASTRAN
             * <ul>
             *  <li> <B>independentVariableWeight = 1.0 or [1.0, 1.0, ...]</B> </li> <br>
             *  Single or list of weighting constants with respect to the variables set for "independentVariable".
             *  If the length of this list doesn't match the length
             *  of the "independentVariable" list, the list is either truncated [ >length("independentVariable")] or expanded [ <length("independentVariable")]
             *  in which case the <b>last weight is repeated</b>.
             * </ul>
             * \endif
             */

            keyWord = "independentVariableWeight";
            status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleDynamicArray(keyValue,
                                                     &tempInteger,
                                                     &designVariable->independVariableWeight);
                AIM_FREE(keyValue);
                AIM_STATUS(aimInfo, status);

                // We have weights, but no variables
                if (designVariable->numIndependVariable == 0) {

                    printf("\tWeighting constants have been provided, but no independent design variables were set!\n");

                    // Less weights than variables
                } else if( tempInteger < designVariable->numIndependVariable) {

                    printf("\tThe number of weighting constants provided does not match the number of independent design variables. "
                            "The last weight will be repeated %d times\n", designVariable->numIndependVariable - tempInteger);

                    AIM_REALL(designVariable->independVariableWeight,
                              designVariable->numIndependVariable, double, aimInfo, status);

                    for (j = 0; j < designVariable->numIndependVariable - tempInteger; j++) {

                        designVariable->independVariableWeight[j+tempInteger] = designVariable->independVariableWeight[tempInteger-1];
                    }

                    // More weights than variables
                } else if (tempInteger > designVariable->numIndependVariable) {

                    printf("\tThe number of weighting constants provided does not match the number of independent design variables. "
                            "The last %d weights will be not be used\n", tempInteger -designVariable->numIndependVariable);

                    AIM_REALL(designVariable->independVariableWeight,
                              designVariable->numIndependVariable, double, aimInfo, status);
                }

            } else { // No weights provided - set default value of 1.0

                if (designVariable->numIndependVariable != 0) {
                    AIM_ALLOC(designVariable->independVariableWeight, designVariable->numIndependVariable, double, aimInfo, status);

                    for (j = 0; j < designVariable->numIndependVariable; j++) {
                        designVariable->independVariableWeight[j] = 1.0;
                    }
                }
            }

            /*! \page feaDesignVariable
             *
             * \if NASTRAN
             * <ul>
             *  <li> <B>variableWeight = [1.0, 1.0]</B> </li> <br>
             *  Weighting constants for a dependent variable - used if "independentVariable"(s) have been provided.
             * </ul>
             * \endif
             */
            keyWord = "variableWeight";
            status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDoubleArray(keyValue,
                                              2,
                                              designVariable->variableWeight);
                AIM_FREE(keyValue);

                if (status != CAPS_SUCCESS) {
                    AIM_ERROR(aimInfo, "Retrieving variableWeight - status %d\n", status);
                    goto cleanup;
                }
            } else {
                designVariable->variableWeight[0] = 1.0;
                designVariable->variableWeight[1] = 1.0;
            }

            // check if design variable relation info also included in input
            keyWord = "fieldName";
            status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);

            AIM_FREE(keyValue);

            if (status != CAPS_SUCCESS) {
                keyWord = "fieldPosition";
                status = search_jsonDictionary( designVariableTuple[i].value, keyWord, &keyValue);
                AIM_FREE(keyValue);
            }

            if (status == CAPS_SUCCESS) {

//                // The main relation associated with this design variable, if defined in input
//                AIM_REALL(feaProblem->feaDesignVariableRelation, relationIndex+1, feaDesignVariableRelationStruct, aimInfo, status);
//
//                designVariableRelation = &feaProblem->feaDesignVariableRelation[relationIndex++];
//
//                feaProblem->numDesignVariableRelation = relationIndex;
//
//                status = initiate_feaDesignVariableRelationStruct(designVariableRelation);
//                AIM_STATUS(aimInfo, status);
//
//                printf("\tWarning: the ability to provide design variable relation data "
//                       "within Design_Variable input is deprecated and "
//                       "will be removed in the future. Please use provide relation data "
//                       "via \"Design_Variable_Relation\" instead.\n");
//
//                status = fea_getDesignVariableRelationEntry( &designVariableTuple[i], designVariableRelation,
//                                                             attrMap, feaProblem, designVariable->name);
//                AIM_STATUS(aimInfo, status);
//
//                designVariableRelation->relationID = relationIndex;

                AIM_ERROR(aimInfo, "The ability to provide design variable relation data "
                                   "within Design_Variable input is deprecated. Please provide relation data "
                                   "via \"Design_Variable_Relation\" instead.\n");
                status = CAPS_BADVALUE;
                goto cleanup;
            }


        } else {

            AIM_ERROR(aimInfo, "Design_Variable tuple value is expected to be a JSON string\n");
            status = CAPS_BADVALUE;
            goto cleanup;

        }
    }

    AIM_FREE(keyValue);

    // Order the design variable ID's such that Analysis and Geometry ID's are continuous.
    // This is required by the TACS AIM

    designVariableID = 1;
    for( i = 0; i < feaProblem->numDesignVariable; i++) {

        // only analysis design variables
        if (aim_getIndex(aimInfo, feaProblem->feaDesignVariable[i].name, GEOMETRYIN) <= 0)
          feaProblem->feaDesignVariable[i].designVariableID = designVariableID++;
    }
    for( i = 0; i < feaProblem->numDesignVariable; i++) {

        // only geometry design variables
        if (aim_getIndex(aimInfo, feaProblem->feaDesignVariable[i].name, GEOMETRYIN) > 0)
          feaProblem->feaDesignVariable[i].designVariableID = designVariableID++;
    }

    // Now that we are done going through all the tuples we need to populate/create the independVaraiableID array
    // if independentVariable was set for any of them.
    for (i = 0; i < feaProblem->numDesignVariable; i++) {

        designVariable = &feaProblem->feaDesignVariable[i];

        if (designVariable->numIndependVariable != 0) {

            AIM_ALLOC(designVariable->independVariableID, designVariable->numIndependVariable, int, aimInfo, status);

            // Loop through the independent variable names
            for (j = 0; j < designVariable->numIndependVariable; j++) {

                // Compare the independent variable names with design variable name
                found = (int) false;
                for (k = 0; k < feaProblem->numDesignVariable; k++) {
                    if (strcasecmp(designVariable->independVariable[j], feaProblem->feaDesignVariable[k].name) == 0) {
                        found = (int) true;
                        break;
                    }
                }

                // If NOT found
                if (found != (int) true) {
                    AIM_ERROR(aimInfo, "\tDesign variable name, \"%s\", not found when searching for independent design variables for "
                                       "variable %s!!!\n", designVariable->independVariable[j],
                                                           designVariable->name);
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

                designVariable->independVariableID[j] = feaProblem->feaDesignVariable[k].designVariableID;
            }
        }
    }

    if (designVariableRelationTuple != NULL) {

        AIM_ALLOC(feaProblem->feaDesignVariableRelation, numDesignVariableRelationTuple, feaDesignVariableRelationStruct, aimInfo, status);
        feaProblem->numDesignVariableRelation = numDesignVariableRelationTuple;

        for (i = 0; i < numDesignVariableRelationTuple; i++) {
          status = initiate_feaDesignVariableRelationStruct(&feaProblem->feaDesignVariableRelation[i]);
          AIM_STATUS(aimInfo, status);
        }

        // Go through design variable relations defined via 'Design_Variable_Relation'
        for (i = 0; i < numDesignVariableRelationTuple; i++) {

            /*! \page feaDesignVariableRelation FEA DesignVariableRelation
             * Structure for the design variable tuple  = ("DesignVariableRelation Name", "Value").
             * "DesignVariableRelation Name" defines the reference name for the design variable being specified.
             *  This string will be used in the FEA input directly. The "Value" must be a JSON String dictionary
             *  (see Section \ref jsonStringDesignVariableRelation).
             *
             */

            designVariableRelation = &feaProblem->feaDesignVariableRelation[i];

            designVariableRelation->relationID = i+1;

            status = fea_getDesignVariableRelationEntry(
                aimInfo, &designVariableRelationTuple[i], designVariableRelation,
                attrMap, feaProblem, NULL
            );
            AIM_STATUS(aimInfo, status);
        }

    }
    status = fea_linkDesignVariablesAndRelations(feaProblem);
    AIM_STATUS(aimInfo, status);

    printf("\tNumber of design variable relations - %d\n", feaProblem->numDesignVariableRelation);
    printf("\tDone getting FEA design variables\n");
    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(keyValue);

    if (groupName != NULL)
        string_freeArray(numGroupName, &groupName);

    return status;
}

int fea_getDesignVariableRelationEntry(void *aimInfo,
                                       capsTuple *designVariableInput,
                          /*@unused@*/ feaDesignVariableRelationStruct *designVariableRelation,
                          /*@unused@*/ mapAttrToIndexStruct *attrMap,
                          /*@unused@*/ feaProblemStruct *feaProblem,
                          /*@unused@*/ char *forceGroupName) {


    int status; //Function return

    int i, j; // Indexing

    int numLinearCoeff;

    char *keyValue = NULL;
    char *keyWord = NULL;

    char **groupName = NULL;
    int  numGroupName = 0;

    feaMeshDataStruct *feaData = NULL;

    int numMaterial;
    feaMaterialStruct **materialSet = NULL, *material;

    int numProperty;
    feaPropertyStruct **propertySet = NULL, *property;

    int numElement;
    meshElementStruct **elementSet = NULL, *element;

    printf("\tDesign_Variable_Relation name - %s\n", designVariableInput->name);

    designVariableRelation->name = EG_strdup(designVariableInput->name);

    /*! \page feaDesignVariableRelation FEA Design Variable Relations
     * \section jsonStringDesignVariableRelation JSON String Dictionary
     *
     * If "Value" is JSON string dictionary
     * \if (NASTRAN || TACS || MASSTRAN)
     *  (eg. "Value" = {"componentType": "Property", "componentName": "plate", "fieldName": "TM", "variableName": "MyDesVar"})
     * \endif
     * \if (ASTROS)
     *  (eg. "Value" = {"componentType": "Property", "componentName": "plate", "variableName": "MyDesVar"})
     * \endif
     *  the following keywords ( = default values) may be used:
     *
     */
    if (!json_isDict(designVariableInput->value)) {
        AIM_ERROR(aimInfo,
            "'Design_Variable_Relation' tuple value must be a JSON dictionary");
        return CAPS_BADVALUE;
    }

    /*! \page feaDesignVariableRelation
     * \if (NASTRAN || TACS || MASSTRAN)
     * <ul>
     * <li> <B>componentType = "Property"</B> </li> <br>
     *  The type of component for this design variable relation. Options: "Material", "Property", "Element".
     * </ul>
     * \endif
     * \if (ASTROS)
     * <ul>
     * <li> <B>componentType = "Property"</B> </li> <br>
     *  The type of component for this design variable relation. Options: "Property".
     * </ul>
     * \endif
     *
     */
    // Get component type for design variable relation
    keyWord = "componentType";
    status = json_getString(
        designVariableInput->value, keyWord, &keyValue);

    // If "componentType" not found, check "variableType" for legacy purposes
    // Warn that "variableType" is deprecated and will be removed in the future!
    if (status != CAPS_SUCCESS) {

        status = json_getString(
            designVariableInput->value, "variableType", &keyValue);

        if (status == CAPS_SUCCESS) {
            printf("\tWarning: \"variableType\" is deprecated and "
                    "will be removed in the future. Please use \"componentType\" "
                    "instead.\n");
        }
    }

    if (status == CAPS_SUCCESS) {
        AIM_NOTNULL(keyValue, aimInfo, status);

        if      (strcasecmp(keyValue, "Material")  == 0) designVariableRelation->componentType = MaterialDesignVar;
        else if (strcasecmp(keyValue, "Property")  == 0) designVariableRelation->componentType = PropertyDesignVar;
        else if (strcasecmp(keyValue, "Element")  == 0) designVariableRelation->componentType = ElementDesignVar;
        else {

            printf("\tUnrecognized \"%s\" specified (%s) for Design_Variable_Relation tuple %s, defaulting to \"Property\"\n", keyWord,
                                                                                                                        keyValue,
                                                                                                                        designVariableInput->name);
            designVariableRelation->componentType = PropertyDesignVar;
        }

    }
    else {
        printf("\tNo \"%s\" specified for Design_Variable_Relation tuple %s, defaulting to \"Property\"\n", keyWord,
                                                                                                    designVariableInput->name);

        designVariableRelation->componentType = PropertyDesignVar;
    }

    AIM_FREE(keyValue);

    /*! \page feaDesignVariableRelation
     *
     * \if (NASTRAN || TACS || MASSTRAN)
     * <ul>
     *  <li> <B>componentName = "(no default)"</B> </li> <br>
     *  Single or list of FEA Property(ies), or FEA Material name(s) linked
     *  to the design variable relation (e.g. "Name1" or ["Name1","Name2",...].
     *  - For <c>componentType</c> Property a \ref feaProperty name (or names) is given.
     *  - For <c>componentType</c> Material a \ref feaMaterial name (or names) is given.
     *  \if (NASTRAN || TACS)
     *  - For <c>componentType</c> Element a <c>capsGroup</c> Name (or names) is given.
     *  \endif
     * </ul>
     *  \endif
     * \if (ASTROS)
     * <ul>
     *  <li> <B>componentName = "(no default)"</B> </li> <br>
     *  Single FEA Property linked to the design variable (e.g. "Name1").
     *  - For <c>componentType</c> Property a \ref feaProperty name (or names) is given.
     * </ul>
     *  \endif
     *
     */
    keyWord = "componentName";
    status = json_getStringDynamicArray(
        designVariableInput->value, keyWord,
        &numGroupName, &groupName);

    if (status != CAPS_SUCCESS) {
        // required
        AIM_ERROR(aimInfo, "No \"%s\" specified for Design_Variable tuple %s",
                    keyWord, designVariableInput->name);
        goto cleanup;
    }
    AIM_NOTNULL(groupName, aimInfo, status);

    // collect materials associated with groupName
    status = fea_findMaterialsByNames(feaProblem, numGroupName, groupName, &numMaterial, &materialSet);

    if (status == CAPS_SUCCESS) {

        AIM_ALLOC(designVariableRelation->materialSetID, numMaterial, int, aimInfo, status);
        AIM_ALLOC(designVariableRelation->materialSetType, numMaterial, int, aimInfo, status);
        designVariableRelation->numMaterialID = numMaterial;

        for (j = 0; j < numMaterial; j++) {

            material = materialSet[j];
            designVariableRelation->materialSetID[j] = material->materialID;
            designVariableRelation->materialSetType[j] = material->materialType;

        }
    }

    // collect properties associated with groupName
    status = fea_findPropertiesByNames(feaProblem, numGroupName, groupName,
                                        &numProperty, &propertySet);

    if (status == CAPS_SUCCESS) {

        AIM_ALLOC(designVariableRelation->propertySetID, numProperty, int, aimInfo, status);
        AIM_ALLOC(designVariableRelation->propertySetType, numProperty, int, aimInfo, status);
        designVariableRelation->numPropertyID = numProperty;

        for (j = 0; j < numProperty; j++) {

            property = propertySet[j];
            designVariableRelation->propertySetID[j] = property->propertyID;
            designVariableRelation->propertySetType[j] = property->propertyType;

        }
    }

    // collect elements associated with groupName
    status = mesh_findGroupElements(
        &feaProblem->feaMesh, attrMap, numGroupName, groupName, &numElement, &elementSet);

    if (status == CAPS_SUCCESS) {

        AIM_ALLOC(designVariableRelation->elementSetID, numElement, int, aimInfo, status);
        AIM_ALLOC(designVariableRelation->elementSetType, numElement, int, aimInfo, status);
        AIM_ALLOC(designVariableRelation->elementSetSubType, numElement, int, aimInfo, status);
        designVariableRelation->numElementID = numElement;

        for (j = 0; j < numElement; j++) {

            element = elementSet[j];
            designVariableRelation->elementSetID[j] = element->elementID;
            designVariableRelation->elementSetType[j] = element->elementType;
            feaData = (feaMeshDataStruct *) element->analysisData;
            designVariableRelation->elementSetSubType[j] = feaData->elementSubType;

        }
    }

    if (groupName != NULL)
        string_freeArray(numGroupName, &groupName);
    numGroupName = 0;
    groupName = NULL;

    AIM_FREE(materialSet);
    AIM_FREE(propertySet);
    AIM_FREE(elementSet);

    /*! \page feaDesignVariableRelation
     *
     * \if (NASTRAN || TACS || MASSTRAN)
     * <ul>
     *  <li> <B>variableName = "(no default)"</B> </li> <br>
     *  Single or list of names of design variables linked to this relation
     * </ul>
     * \endif
     * \if (ASTROS)
     * <ul>
     *  <li> <B>variableName = "(no default)"</B> </li> <br>
     *  Name of design variable linked to this relation
     * </ul>
     * \endif
     */
    keyWord = "variableName";
    status = json_getStringDynamicArray(
        designVariableInput->value, keyWord,
        &designVariableRelation->numDesignVariable, &designVariableRelation->designVariableNameSet);
    if (status != CAPS_SUCCESS) {
        // required
        AIM_ERROR(aimInfo, "No \"%s\" specified for Design_Variable_Relation tuple %s",
                    keyWord, designVariableInput->name);
        goto cleanup;
    }

    /*! \page feaDesignVariableRelation
     *
     * \if (NASTRAN || TACS || MASSTRAN)
     * <ul>
     *  <li> <B>fieldName = "(no default)"</B> </li> <br>
     *  Fieldname of variable relation (e.g. "E" for Young's Modulus). Design Variable Relations can be defined as three types based on the <c>variableType</c> value.
     *  These are Material, Property, or Element.  This means that an aspect of a material, property, or element input can change in the optimization problem.  This input
     *  specifies what aspect of the Material, Property, or Element is changing.
     *  \if (MASSTRAN)
     *  -# <b> Material Types</b> Selected based on the material type (see \ref feaMaterial, materialType) referenced  in the <c>componentName</c> above.
     *      - <c> materialType = "Isotropic" </c>
     *          - "density"
     *  -# <b> Property Types</b> (see \ref feaProperty)
     *      - <c>propertyType = "ConcentratedMass"</c>
     *          - "mass",
     *          - "massOffset1", "massOffset2", "massOffset2"
     *          - "Ixx", "Iyy", "Izz", "Ixy", "Ixz", "Iyz"
     *      - <c>propertyType = "Shell"</c>
     *          - "membraneThickness", "massPerArea"
     *  \else
     *  -# <b> Material Types</b> Selected based on the material type (see \ref feaMaterial, materialType) referenced  in the <c>componentName</c> above.
     *      - <c><b> MAT1</b>,     materialType = "Isotropic" </c>
     *          - "E", "G", "NU", "RHO", "A"
     *      - <c><b>  MAT2</b>, materialType = "Anisothotropic"  </c>
     *          - "G11", "G12", "G13", "G22", "G23", "G33", "RHO", "A1", "A2", "A3"
     *      - <c><b>  MAT8</b>, materialType = "Orthotropic"  </c>
     *          - "E1", "E2", "NU12", "G12", "G1Z", "G2Z", "RHO", "A1", "A2"
     *      - <c><b>  MAT9</b>, materialType = "Anisotropic"  </c>
     *          - "G11", "G12", "G13", "G14", "G15", "G16"
     *          - "G22", "G23", "G24", "G25", "G26"
     *          - "G33", "G34", "G35", "G36"
     *          - "G44", "G45", "G46"
     *          - "G55", "G56", "G66"
     *          - "RHO", "A1", "A2", "A3", "A4", "A5", "A6"
     *  -# <b> Property Types</b> (see \ref feaProperty)
     *      - <c><b> PROD</b> </c> <c>propertyType = "Rod"</c>
     *          - "A", "J"
     *      - <c><b> PBAR</b> </c> <c>propertyType = "Bar"</c>
     *          - "A", "I1", "I2", "J"
     *      - <c><b> PSHELL</b> </c> <c>propertyType = "Shell"</c>
     *          - "T"
     *      - <c><b> PCOMP</b> </c> <c>propertyType = "Composite"</c>
     *          - "T1", "THETA1", "T2", "THETA2", ... "Ti", "THETAi"
     *      - <c><b> PSOLID</b> </c> <c>propertyType = "Solid"</c>
     *          - not supported
     *  -# <b> Element Types</b>
     *      - <c><b> CTRIA3, CQUAD4</b> </c> <c>propertyType = "Shell"</c>
     *          - "ZOFFS"
     *  \endif
     * </ul>
     * \endif
     */
    keyWord = "fieldName";
    status = json_getString(
        designVariableInput->value, keyWord,
        &designVariableRelation->fieldName);
    if (status != CAPS_SUCCESS) {
        // optional
    }

    /*! \page feaDesignVariableRelation
     *
     * \if (NASTRAN || TACS)
     * <ul>
     *  <li> <B>fieldPosition = 0</B> </li> <br>
     *  This input is ignored if not defined.  The user may use this field instead of the <c>fieldName</c> input defined above to
     *  relate design variables and property, material, or elements.  This requires knowledge of Nastran bulk data input format for material,
     *   property, and element input cards.
     * </ul>
     * \endif
     */

    keyWord = "fieldPosition";
    status = json_getInteger(
        designVariableInput->value, keyWord,
        &designVariableRelation->fieldPosition);
    if (status != CAPS_SUCCESS) {
        // optional
    }

    /*! \page feaDesignVariableRelation
     *
     * \if (NASTRAN || TACS || ASTROS || MASSTRAN)
     * <ul>
     *  <li> <B>constantCoeff = 0.0</B> </li> <br>
     *  Constant term of relation.
     * </ul>
     * \endif
     */
    keyWord = "constantCoeff";
    status = json_getDouble(
        designVariableInput->value, keyWord,
        &designVariableRelation->constantRelationCoeff);
    if (status != CAPS_SUCCESS) {
        // default
        designVariableRelation->constantRelationCoeff = 0.0;
    }

    /*! \page feaDesignVariableRelation
     *
     * \if (NASTRAN || TACS || ASTROS || MASSTRAN)
     * <ul>
     *  <li> <B>linearCoeff = 1.0</B> </li> <br>
     *  Single or list of coefficients of linear relation. Must be same length as <c>variableName</c>.
     * </ul>
     * \endif
     */
    keyWord = "linearCoeff";
    status = json_getDoubleDynamicArray(
        designVariableInput->value, keyWord,
        &numLinearCoeff, &designVariableRelation->linearRelationCoeff);

    if (status != CAPS_SUCCESS) {
        // default
        numLinearCoeff = designVariableRelation->numDesignVariable;
        designVariableRelation->linearRelationCoeff = EG_alloc(numLinearCoeff * sizeof(double));
        if (designVariableRelation->linearRelationCoeff == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        for (i = 0; i < numLinearCoeff; i++) {
            designVariableRelation->linearRelationCoeff[i] = 1.0;
        }
    }

    if (numLinearCoeff != designVariableRelation->numDesignVariable) {
        AIM_ERROR(aimInfo, "Number of \"linearCoeff\" values (%d) does not match"
                    " number of \"variableName\" values (%d)",
                    numLinearCoeff, designVariableRelation->numDesignVariable);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(keyValue);

    return status;
}

// Get the design constraints from a capsTuple
int fea_getDesignConstraint(void *aimInfo,
                            int numDesignConstraintTuple,
                            capsTuple designConstraintTuple[],
                            feaProblemStruct *feaProblem) {

    /*! \page feaDesignConstraint FEA Design Constraints
     * Structure for the design constraint tuple  = (`DesignConstraint Name', `Value').
     * 'DesignConstraint Name' defines the reference name for the design constraint being specified.
     * The "Value" must be a JSON String dictionary (see Section \ref jsonStringDesignConstraint).
     */

    int status; //Function return

    int i, groupIndex, attrIndex; // Indexing

    char *keyValue = NULL;
    char *keyWord = NULL;

    char *tempString = NULL;

    char **groupName = NULL;
    int  numGroupName = 0;

    // Destroy our design constraints structures coming in if aren't 0 and NULL already
    if (feaProblem->feaDesignConstraint != NULL) {
        for (i = 0; i < feaProblem->numDesignConstraint; i++) {
            status = destroy_feaDesignConstraintStruct(&feaProblem->feaDesignConstraint[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (feaProblem->feaDesignConstraint != NULL) EG_free(feaProblem->feaDesignConstraint);
    feaProblem->feaDesignConstraint = NULL;
    feaProblem->numDesignConstraint = 0;

    printf("\nGetting FEA design constraints.......\n");

    feaProblem->numDesignConstraint = numDesignConstraintTuple;

    printf("\tNumber of design constraints - %d\n", feaProblem->numDesignConstraint);

    if (feaProblem->numDesignConstraint > 0) {
        feaProblem->feaDesignConstraint = (feaDesignConstraintStruct *) EG_alloc(feaProblem->numDesignConstraint * sizeof(feaDesignConstraintStruct));

        if (feaProblem->feaDesignConstraint == NULL) return EGADS_MALLOC;

    } else {
        printf("\tNumber of design constraint values in input tuple is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < feaProblem->numDesignConstraint; i++) {
        status = initiate_feaDesignConstraintStruct(&feaProblem->feaDesignConstraint[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    for (i = 0; i < feaProblem->numDesignConstraint; i++) {

        printf("\tDesign_Constraint name - %s\n", designConstraintTuple[i].name);

        feaProblem->feaDesignConstraint[i].name = (char *) EG_alloc(((strlen(designConstraintTuple[i].name)) + 1)*sizeof(char));
        if (feaProblem->feaDesignConstraint[i].name == NULL) return EGADS_MALLOC;

        memcpy(feaProblem->feaDesignConstraint[i].name, designConstraintTuple[i].name, strlen(designConstraintTuple[i].name)*sizeof(char));
        feaProblem->feaDesignConstraint[i].name[strlen(designConstraintTuple[i].name)] = '\0';

        feaProblem->feaDesignConstraint[i].designConstraintID = i + 1;

        // Do we have a json string?
        if (strncmp(designConstraintTuple[i].value, "{", 1) == 0) {
            //printf("JSON String - %s\n", designConstraintTuple[i].value);

            /*! \page feaDesignConstraint
             * \section jsonStringDesignConstraint JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (NASTRAN || ASTROS)
             *  (eg. "Value" = {"groupName": "plate", "upperBound": 10.0})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsGroup</c> name(s)
             *  to the design variable (e.g. "Name1" or ["Name1","Name2",...].The property (see \ref feaProperty) also
             *  assigned to the same <c>capsGroup</c> will be automatically related to this constraint entry.
             * </ul>
             * \endif
             *
             */

            // Get material/properties that the design constraint should be applied to
            keyWord = "groupName";
            status = search_jsonDictionary( designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaDesignConstraint[i].numPropertyID = 0;
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numProperty; attrIndex++) {

                        if (strcasecmp(feaProblem->feaProperty[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaDesignConstraint[i].numPropertyID += 1;

                            if (feaProblem->feaDesignConstraint[i].numPropertyID == 1)  {

                                feaProblem->feaDesignConstraint[i].propertySetID   = (int *) EG_alloc(feaProblem->feaDesignConstraint[i].numPropertyID *sizeof(int));
                                feaProblem->feaDesignConstraint[i].propertySetType = (int *) EG_alloc(feaProblem->feaDesignConstraint[i].numPropertyID *sizeof(int));

                            } else {

                                feaProblem->feaDesignConstraint[i].propertySetID   = (int *) EG_reall(feaProblem->feaDesignConstraint[i].propertySetID,
                                                                                                      feaProblem->feaDesignConstraint[i].numPropertyID *sizeof(int));
                                feaProblem->feaDesignConstraint[i].propertySetType = (int *) EG_reall(feaProblem->feaDesignConstraint[i].propertySetType,
                                                                                                      feaProblem->feaDesignConstraint[i].numPropertyID *sizeof(int));

                            }

                            if (feaProblem->feaDesignConstraint[i].propertySetID   == NULL ||
                                feaProblem->feaDesignConstraint[i].propertySetType == NULL    ) {

                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaDesignConstraint[i].propertySetID[feaProblem->feaDesignConstraint[i].numPropertyID-1]   = feaProblem->feaProperty[attrIndex].propertyID;
                            feaProblem->feaDesignConstraint[i].propertySetType[feaProblem->feaDesignConstraint[i].numPropertyID-1] = feaProblem->feaProperty[attrIndex].propertyType;

                            break;
                        }
                    }

                    if (feaProblem->feaDesignConstraint[i].numPropertyID != groupIndex+1) {

                        printf("\tWarning: DesignConstraint property name, %s, not found in feaProperty structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;
                groupName = NULL;
            }

            //Fill up designConstraint properties
            /*! \page feaDesignConstraint
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B> constraintType = "Property" </B> </li> <br>
             *  The type of design constraint. Options: "Property", "Flutter"
             * </ul>
             * \endif
             */
            keyWord = "constraintType";
            status = json_getString(
                designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);
                if (strcasecmp(keyValue, "Property") == 0) {
                    feaProblem->feaDesignConstraint[i].designConstraintType = PropertyDesignCon;
                }
                else if (strcasecmp(keyValue, "Flutter") == 0) {
                    feaProblem->feaDesignConstraint[i].designConstraintType = FlutterDesignCon;
                }
                else if (strcasecmp(keyValue, "EquationResponse") == 0) {
                    feaProblem->feaDesignConstraint[i].designConstraintType = EquationResponseDesignCon;
                }
                else {
                    PRINT_WARNING("Unknown constraintType %s... defaulting to Property.", keyValue);
                    feaProblem->feaDesignConstraint[i].designConstraintType = PropertyDesignCon;
                }
            }
            else {
                feaProblem->feaDesignConstraint[i].designConstraintType = PropertyDesignCon;
            }

            if (keyValue != NULL) EG_free(keyValue);
            keyValue = NULL;

            /*! \page feaDesignConstraint
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>lowerBound = 0.0</B> </li> <br>
             *  Lower bound for the design constraint.
             * </ul>
             * \endif
             *
             */
            keyWord = "lowerBound";
            status = search_jsonDictionary( designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaDesignConstraint[i].lowerBound);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaDesignConstraint
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>upperBound = 0.0</B> </li> <br>
             *  Upper bound for the design constraint.
             * </ul>
             * \endif
             *
             */
            keyWord = "upperBound";
            status = search_jsonDictionary( designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaDesignConstraint[i].upperBound);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

            /*! \page feaDesignConstraint
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>responseType = "(no default)"</B> </li> <br>
             *  Response type options for DRESP1 Entry (see Nastran manual).
             *      - Implemented Options
             *          -# <c>STRESS</c>, for <c>propertyType = "Rod" or "Shell"</c> (see \ref feaProperty)
             *          -# <c>CFAILURE</c>, for <c>propertyType = "Composite"</c> (see \ref feaProperty)
             * </ul>
             * \endif
             *
             */
            keyWord = "responseType";
            status = search_jsonDictionary( designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                feaProblem->feaDesignConstraint[i].responseType = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaDesignConstraint[i].responseType == NULL) {
                    if (tempString != NULL) EG_free(tempString);

                    return EGADS_MALLOC;
                }

                memcpy(feaProblem->feaDesignConstraint[i].responseType, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaDesignConstraint[i].responseType[strlen(tempString)] = '\0';

                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                AIM_FREE(tempString);
            }

            /*! \page feaDesignConstraint
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>fieldName = "(no default)"</B> </li> <br>
             *  For constraints, this field is only used currently when applying constraints to composites.  This field is used to identify
             *  the specific lamina in a stacking sequence that a constraint is being applied too.  Note if the user has design variables
             *  for both THEATA1 and T1 it is likely that only a single constraint on the first lamina is required.  For this reason, the user
             *  can simply enter LAMINA1 in addition to the possible entries defined in the \ref feaDesignVariable section.
             *  Additionally, the <c>fieldPosition</c> integer entry below can be used.  In this case <c>"LAMINA1" = 1</c>.
             *
             *  *  -# <b> Property Types</b> (see \ref feaProperty)
             *      - <c><b> PCOMP</b> </c> <c>propertyType = "Composite"</c>
             *          - "T1", "THETA1", "T2", "THETA2", ... "Ti", "THETAi"
             *          - "LAMINA1", "LAMINA2", ... "LAMINAi"
             *
             * </ul>
             * \endif
             */

            keyWord = "fieldName";
            status = search_jsonDictionary( designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                tempString = string_removeQuotation(keyValue);
                AIM_NOTNULL(tempString, aimInfo, status);

                feaProblem->feaDesignConstraint[i].fieldName = EG_alloc((strlen(tempString) + 1)*sizeof(char));
                if(feaProblem->feaDesignConstraint[i].fieldName == NULL) {
                    if (tempString != NULL) EG_free(tempString);

                    return EGADS_MALLOC;
                }

                memcpy(feaProblem->feaDesignConstraint[i].fieldName, tempString, strlen(tempString)*sizeof(char));
                feaProblem->feaDesignConstraint[i].fieldName[strlen(tempString)] = '\0';

                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (tempString != NULL) {
                    EG_free(tempString);
                    tempString = NULL;
                }
            }

            /*! \page feaDesignConstraint
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>fieldPosition = 0</B> </li> <br>
             *  This input is ignored if not defined.  The user may use this field instead of the <c>fieldName</c> input defined above to
             *  identify a specific lamina in a composite stacking sequence where a constraint is applied.  Please read the <c>fieldName</c>
             *  information above for more information.
             * </ul>
             * \endif
             */

            keyWord = "fieldPosition";
            status = search_jsonDictionary( designConstraintTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toInteger(keyValue, &feaProblem->feaDesignConstraint[i].fieldPosition);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                AIM_STATUS(aimInfo, status);
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>velocityType = "TRUE" </B> </li> <br>
             *  The nature of the velocity values defined with "velocity" keyword. Can be either "TRUE" for true velocity
             *  or "EQUIV" for equivalent air speed.
             * </ul>
             * \endif
             */
            keyWord = "velocityType";
            status = json_getString(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].velocityType);
            if (status != CAPS_SUCCESS) {
                feaProblem->feaDesignConstraint[i].velocityType = EG_strdup("TRUE");
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>scalingFactor = 0.10 </B> </li> <br>
             *  The constraint scaling factor.
             * </ul>
             * \endif
             */
            keyWord = "scalingFactor";
            status = json_getDouble(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].scalingFactor);
            if (status != CAPS_SUCCESS) {
                feaProblem->feaDesignConstraint[i].scalingFactor = 0.10;
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>density = (No Default) </B> </li> <br>
             *  The density values for flutter constraint.
             *  Must be an array with equal length to defined "damping" values.
             * </ul>
             * \endif
             */
            keyWord = "density";
            status = json_getDoubleDynamicArray(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].numDensity,
                &feaProblem->feaDesignConstraint[i].density);
            if (status != CAPS_SUCCESS) {
                // only required if constraintType == FlutterCon
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>Mach = (No Default) </B> </li> <br>
             *  The Mach values for flutter constraint.
             *  Must be an array with equal length to defined "damping" values.
             * </ul>
             * \endif
             */
            keyWord = "Mach";
            status = json_getDoubleDynamicArray(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].numMach,
                &feaProblem->feaDesignConstraint[i].Mach);
            if (status != CAPS_SUCCESS) {
                // only required if constraintType == FlutterCon
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>velocity = (No Default) </B> </li> <br>
             *  The velocity values for flutter constraint.
             *  Must be an array with equal length to defined "damping" values.
             * </ul>
             * \endif
             */
            keyWord = "velocity";
            status = json_getDoubleDynamicArray(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].numVelocity,
                &feaProblem->feaDesignConstraint[i].velocity);
            if (status != CAPS_SUCCESS) {
                // only required if constraintType == FlutterCon
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>modes = (No Default) </B> </li> <br>
             *  The modes values for flutter constraint.
             *  Must be an array with equal length to defined "damping" values.
             * </ul>
             * \endif
             */
            keyWord = "modes";
            status = json_getIntegerDynamicArray(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].numModes,
                &feaProblem->feaDesignConstraint[i].modes);
            if (status != CAPS_SUCCESS) {
                // only required if constraintType == FlutterCon
            }

            /*! \page feaDesignConstraint
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>damping = (No Default) </B> </li> <br>
             *  The damping values for flutter constraint.
             *  Must be an array with equal length to defined "velocity" values.
             * </ul>
             * \endif
             */
            keyWord = "damping";
            status = json_getDoubleDynamicArray(
                designConstraintTuple[i].value, keyWord,
                &feaProblem->feaDesignConstraint[i].numDamping,
                &feaProblem->feaDesignConstraint[i].damping);
            if (status != CAPS_SUCCESS) {
                // only required if constraintType == FlutterCon
            }

        } else {

            printf("\tError: Design_Constraint tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;
            // CALL designConstraint look up

        }
    }

    printf("Done getting FEA design constraints\n");
    status = CAPS_SUCCESS;
cleanup:
    AIM_FREE(keyValue);
    AIM_FREE(tempString);
    return status;
}

// Get the design constraints from a capsTuple
int fea_getOptimizationControl(char *optimizationControlInput,
                               feaProblemStruct *feaProblem) {

    /*! \page feaOptimizationControl FEA Optimization Control
     * Structure for the optimization control dictionary  = `Value'.
     * The "Value" must be a JSON String dictionary (see Section \ref jsonStringDesignConstraint).
     */

    int status; //Function return

    char *keyValue = NULL;
    char *keyWord = NULL;

    // Destroy our optimization control structures coming in if aren't 0 and NULL already
    status = destroy_feaOptimzationControlStruct(&feaProblem->feaOptimizationControl);
    if (status != CAPS_SUCCESS) return status;

    printf("\nGetting FEA optimization control.......\n");

    status = initiate_feaOptimizationControlStruct(&feaProblem->feaOptimizationControl);
    if (status != CAPS_SUCCESS) return status;

    if (optimizationControlInput == NULL) {
        optimizationControlInput = "{}";
    }

    // Do we have a json string?
    if (strncmp(optimizationControlInput, "{", 1) == 0) {

        /*! \page feaOptimizationControl
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>fullyStressedDesign = 0</B> </li> <br>
         *  Number of iterations with fully stressed design.
         * </ul>
         * \endif
         *
         */
        keyWord = "fullyStressedDesign";
        status = json_getInteger(
            optimizationControlInput, keyWord,
            &feaProblem->feaOptimizationControl.fullyStressedDesign);

        // Set default value if value is not in JSON
        if (status != CAPS_SUCCESS) {
            feaProblem->feaOptimizationControl.fullyStressedDesign = 0;
        }

        /*! \page feaOptimizationControl
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>mathProgramming = 30</B> </li> <br>
         *  Number of iterations for math programming methods.
         * </ul>
         * \endif
         *
         */
        keyWord = "mathProgramming";
        status = json_getInteger(
            optimizationControlInput, keyWord,
            &feaProblem->feaOptimizationControl.mathProgramming);

        // Set default value if value is not in JSON
        if (status != CAPS_SUCCESS) {
            feaProblem->feaOptimizationControl.mathProgramming = 30;
        }

        /*! \page feaOptimizationControl
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>maxIter = 30</B> </li> <br>
         *  Maximum number of optimization iterations.
         * </ul>
         * \endif
         *
         */
        keyWord = "maxIter";
        status = json_getInteger(
            optimizationControlInput, keyWord,
            &feaProblem->feaOptimizationControl.maxIter);

        // Set default value if value is not in JSON
        if (status != CAPS_SUCCESS) {
            feaProblem->feaOptimizationControl.maxIter = 30;
        }

        /*! \page feaOptimizationControl
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>constraintRetention = 1.5</B> </li> <br>
         *  Constraint retention factor.
         *  Will be at least 1.5 times the number of design variables
         * </ul>
         * \endif
         *
         */
        keyWord = "constraintRetention";
        status = json_getDouble(
            optimizationControlInput, keyWord,
            &feaProblem->feaOptimizationControl.constraintRetention);

        // Set default value if value is not in JSON
        if (status != CAPS_SUCCESS) {
            feaProblem->feaOptimizationControl.constraintRetention = 1.5;
        }

        /*! \page feaOptimizationControl
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>eps = 1.0</B> </li> <br>
         *  Constraint retention parameter in which all constraints having a
         *  value greater than "eps" will be considered active.
         *
         * </ul>
         * \endif
         *
         */
        keyWord = "eps";
        status = json_getDouble(
            optimizationControlInput, keyWord,
            &feaProblem->feaOptimizationControl.eps);

        // Set default value if value is not in JSON
        if (status != CAPS_SUCCESS) {
            feaProblem->feaOptimizationControl.eps = 1.0;
        }

        /*! \page feaOptimizationControl
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>moveLimit = 1.0</B> </li> <br>
         *  Move limit bound.
         * </ul>
         * \endif
         *
         */
        keyWord = "moveLimit";
        status = json_getDouble(
            optimizationControlInput, keyWord,
            &feaProblem->feaOptimizationControl.moveLimit);

        // Set default value if value is not in JSON
        if (status != CAPS_SUCCESS) {
            feaProblem->feaOptimizationControl.moveLimit = 1.5;
        }

    } else {

        printf("\tError: OptimizationControl tuple value is expected to be a JSON string\n");
        return CAPS_BADVALUE;

    }

    if (keyValue != NULL) {
        EG_free(keyValue);
        keyValue = NULL;
    }

    printf("Done getting FEA Optimization Control\n");
    return CAPS_SUCCESS;
}

// Get the design constraints from a capsTuple
int fea_getMassIncrement(int numMassIncrementTuple,
                         capsTuple massIncrementTuple[],
                         feaProblemStruct *feaProblem) {

    /*! \page feaMassIncrement FEA Mass Increments
     * Structure for the mass increment tuple  = (`MassIncrement Name', `Value').
     * 'MassIncrement Name' defines the reference name for the mass increment being specified.
     * The "Value" must be a JSON String dictionary (see Section \ref jsonStringDesignConstraint).
     */

    int status; //Function return

    int i, groupIndex, attrIndex; // Indexing

    char *keyValue = NULL;
    char *keyWord = NULL;

    char **groupName = NULL;
    int  numGroupName = 0;

    // Destroy our design constraints structures coming in if aren't 0 and NULL already
    if (feaProblem->feaMassIncrement != NULL) {
        for (i = 0; i < feaProblem->numMassIncrement; i++) {
            status = destroy_feaMassIncrementStruct(&feaProblem->feaMassIncrement[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (feaProblem->feaMassIncrement != NULL) EG_free(feaProblem->feaMassIncrement);
    feaProblem->feaMassIncrement = NULL;
    feaProblem->numMassIncrement = 0;

    printf("\nGetting FEA Mass Increment.......\n");

    feaProblem->numMassIncrement = numMassIncrementTuple;

    printf("\tNumber of analysis case mass increments - %d\n", feaProblem->numMassIncrement);

    if (feaProblem->numMassIncrement > 0) {
        feaProblem->feaMassIncrement = (feaMassIncrementStruct *) EG_alloc(feaProblem->numMassIncrement * sizeof(feaMassIncrementStruct));

        if (feaProblem->feaMassIncrement == NULL) return EGADS_MALLOC;

    } else {
        printf("\tNumber of mass increment values in input tuple is 0\n");
        return CAPS_NOTFOUND;
    }

    for (i = 0; i < feaProblem->numMassIncrement; i++) {
        status = initiate_feaMassIncrementStruct(&feaProblem->feaMassIncrement[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    for (i = 0; i < feaProblem->numMassIncrement; i++) {

        printf("\tMass_Increment name - %s\n", massIncrementTuple[i].name);

        feaProblem->feaMassIncrement[i].name = (char *) EG_alloc(((strlen(massIncrementTuple[i].name)) + 1)*sizeof(char));
        if (feaProblem->feaMassIncrement[i].name == NULL) return EGADS_MALLOC;

        memcpy(feaProblem->feaMassIncrement[i].name, massIncrementTuple[i].name, strlen(massIncrementTuple[i].name)*sizeof(char));
        feaProblem->feaMassIncrement[i].name[strlen(massIncrementTuple[i].name)] = '\0';

        feaProblem->feaMassIncrement[i].massIncrementID = i + 1;

        // Do we have a json string?
        if (strncmp(massIncrementTuple[i].value, "{", 1) == 0) {
//            printf("JSON String - %s\n", massIncrementTuple[i].value);

            /*! \page feaDesignConstraint
             * \section jsonStringDesignConstraint JSON String Dictionary
             *
             * If "Value" is JSON string dictionary
             * \if (NASTRAN || ASTROS)
             *  (eg. "Value" = {"groupName": "plate", "upperBound": 10.0})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * \if (NASTRAN || ASTROS)
             * <ul>
             *  <li> <B>groupName = "(no default)"</B> </li> <br>
             *  Single or list of <c>capsGroup</c> name(s)
             *  to the design variable (e.g. "Name1" or ["Name1","Name2",...].The property (see \ref feaProperty) also
             *  assigned to the same <c>capsGroup</c> will be automatically related to this constraint entry.
             * </ul>
             * \endif
             *
             */

            // Get material/properties that the design constraint should be applied to
            keyWord = "groupName";
            status = search_jsonDictionary( massIncrementTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toStringDynamicArray(keyValue, &numGroupName, &groupName);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }

                if (status != CAPS_SUCCESS) return status;

                feaProblem->feaMassIncrement[i].numPropertyID = 0;
//                printf("numGroup: %d, numProp: %d\n", numGroupName, feaProblem->numProperty);
                for (groupIndex = 0; groupIndex < numGroupName; groupIndex++) {

                    for (attrIndex = 0; attrIndex < feaProblem->numProperty; attrIndex++) {

//                        printf("%s : %s\n",groupName[groupIndex], feaProblem->feaProperty[attrIndex].name);
                        if (strcasecmp(feaProblem->feaProperty[attrIndex].name, groupName[groupIndex]) == 0) {

                            feaProblem->feaMassIncrement[i].numPropertyID += 1;

                            if (feaProblem->feaMassIncrement[i].numPropertyID == 1)  {

                                feaProblem->feaMassIncrement[i].propertySetID   = (int *) EG_alloc(feaProblem->feaMassIncrement[i].numPropertyID *sizeof(int));
                                feaProblem->feaMassIncrement[i].propertySetType = (int *) EG_alloc(feaProblem->feaMassIncrement[i].numPropertyID *sizeof(int));

                            } else {

                                feaProblem->feaMassIncrement[i].propertySetID   = (int *) EG_reall(feaProblem->feaMassIncrement[i].propertySetID,
                                                                                                   feaProblem->feaMassIncrement[i].numPropertyID *sizeof(int));
                                feaProblem->feaMassIncrement[i].propertySetType = (int *) EG_reall(feaProblem->feaMassIncrement[i].propertySetType,
                                                                                                   feaProblem->feaMassIncrement[i].numPropertyID *sizeof(int));

                            }

                            if (feaProblem->feaMassIncrement[i].propertySetID   == NULL ||
                                feaProblem->feaMassIncrement[i].propertySetType == NULL    ) {

                                status = string_freeArray(numGroupName, &groupName);
                                if (status != CAPS_SUCCESS) return status;
                                groupName = NULL;
                                return EGADS_MALLOC;
                            }

                            feaProblem->feaMassIncrement[i].propertySetID[feaProblem->feaMassIncrement[i].numPropertyID-1]   = feaProblem->feaProperty[attrIndex].propertyID;
                            feaProblem->feaMassIncrement[i].propertySetType[feaProblem->feaMassIncrement[i].numPropertyID-1] = feaProblem->feaProperty[attrIndex].propertyType;

                            break;
                        }
                    }

                    if (feaProblem->feaMassIncrement[i].numPropertyID != groupIndex+1) {

                        printf("\tWarning: MassIncrement property name, %s, not found in feaProperty structure\n", groupName[groupIndex]);

                    }
                }

                status = string_freeArray(numGroupName, &groupName);
                if (status != CAPS_SUCCESS) return status;
                groupName = NULL;
            }

            /*! \page feaMassIncrement
             *
             * \if ASTROS
             * <ul>
             *  <li> <B>increment = (No Default) </B> </li> <br>
             *  The increment value for analysis mass increment.
             * </ul>
             * \endif
             */

            keyWord = "increment";
            status = search_jsonDictionary(massIncrementTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &feaProblem->feaMassIncrement[i].increment);
                if (keyValue != NULL) {
                    EG_free(keyValue);
                    keyValue = NULL;
                }
                if (status != CAPS_SUCCESS) return status;
            }

        } else {

            printf("\tError: Mass_Increment tuple value is expected to be a JSON string\n");
            return CAPS_BADVALUE;
            // CALL designConstraint look up

        }
    }

    if (keyValue != NULL) {
        EG_free(keyValue);
        keyValue = NULL;
    }

    printf("Done getting FEA mass increment\n");
    return CAPS_SUCCESS;
}

// Get the coordinate system information from the bodies and an attribute map (of CoordSystem)
int fea_getCoordSystem(int numBody,
                       ego bodies[],
                       mapAttrToIndexStruct coordSystemMap,
                       int *numCoordSystem,
                       feaCoordSystemStruct *feaCoordSystem[]) {

    int status; // Function return

    int i, body, face, edge, node; // Indexing

    int numFace = 0, numEdge = 0, numNode = 0; // Number of egos
    ego *faces = NULL, *edges = NULL, *nodes = NULL; // Geometry

    int atype, alen; // EGADS return variables
    const int    *ints;
    const double *reals;
    const char *string;

    int found = (int) false;

    // Destroy our CoordSystem structures coming in if aren't 0 and NULL already
    if (*feaCoordSystem != NULL) {
        for (i = 0; i < *numCoordSystem; i++) {
            status = destroy_feaCoordSystemStruct(&(*feaCoordSystem)[i]);
            if (status != CAPS_SUCCESS) return status;
        }
    }
    if (*feaCoordSystem != NULL) EG_free(*feaCoordSystem);
    *feaCoordSystem = NULL;
    *numCoordSystem = 0;

    printf("\nGetting FEA coordinate systems.......\n");

    *numCoordSystem = coordSystemMap.numAttribute;
    printf("\tNumber of coordinate systems - %d\n", *numCoordSystem);

    if (*numCoordSystem > 0) {
        *feaCoordSystem = (feaCoordSystemStruct *) EG_alloc((*numCoordSystem)*sizeof(feaCoordSystemStruct));
        if (*feaCoordSystem == NULL) return EGADS_MALLOC;

    } else {
        printf("\tNo coordinate systems found - defaulting to global\n");
        return CAPS_SUCCESS;
    }

    for (i = 0; i < *numCoordSystem; i++) {
        status = initiate_feaCoordSystemStruct(&(*feaCoordSystem)[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    for (i = 0; i < *numCoordSystem; i++) {

        if (faces != NULL) EG_free(faces);
        if (edges != NULL) EG_free(edges);
        if (nodes != NULL) EG_free(nodes);

        faces = NULL;
        edges = NULL;
        nodes = NULL;

        printf("\tCoordinate system name - %s\n", coordSystemMap.attributeName[i]);

        (*feaCoordSystem)[i].name = EG_strdup(coordSystemMap.attributeName[i]);
        if ((*feaCoordSystem)[i].name == NULL) { status = EGADS_MALLOC; goto cleanup; }

        (*feaCoordSystem)[i].coordSystemID = coordSystemMap.attributeIndex[i];

        (*feaCoordSystem)[i].refCoordSystemID = 0;

        (*feaCoordSystem)[i].coordSystemType = RectangularCoordSystem;

        found = (int) false;

        // Search through bodies
        for (body = 0; body < numBody; body++) {

            // Look at the body level
            status = EG_attributeRet(bodies[body], (*feaCoordSystem)[i].name, &atype, &alen, &ints, &reals, &string);
            if (status != EGADS_SUCCESS &&
                status != EGADS_NOTFOUND) goto cleanup;

            if (atype == ATTRCSYS) {

                // Save the origin
                (*feaCoordSystem)[i].origin[0] = reals[alen+0];
                (*feaCoordSystem)[i].origin[1] = reals[alen+1];
                (*feaCoordSystem)[i].origin[2] = reals[alen+2];

                (*feaCoordSystem)[i].normal1[0] = reals[alen+3];
                (*feaCoordSystem)[i].normal1[1] = reals[alen+4];
                (*feaCoordSystem)[i].normal1[2] = reals[alen+5];

                (*feaCoordSystem)[i].normal2[0] = reals[alen+6];
                (*feaCoordSystem)[i].normal2[1] = reals[alen+7];
                (*feaCoordSystem)[i].normal2[2] = reals[alen+8];

                (*feaCoordSystem)[i].normal3[0] = reals[alen+9];
                (*feaCoordSystem)[i].normal3[1] = reals[alen+10];
                (*feaCoordSystem)[i].normal3[2] = reals[alen+11];

                found = (int) true;
            }

            if (found == (int) true)  break;

            // Determine the number of faces, edges, and nodes
            status = EG_getBodyTopos(bodies[body], NULL, FACE, &numFace, &faces);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_getBodyTopos(bodies[body], NULL, EDGE, &numEdge, &edges);
            if (status != EGADS_SUCCESS) goto cleanup;

            status = EG_getBodyTopos(bodies[body], NULL, NODE, &numNode, &nodes);
            if (status != EGADS_SUCCESS) goto cleanup;

            // Loop through faces
            for (face = 0; face < numFace; face++) {

                status = EG_attributeRet(faces[face], (const char *) (*feaCoordSystem)[i].name, &atype, &alen, &ints, &reals, &string);

                if (status == EGADS_NOTFOUND) continue;
                if (status != EGADS_SUCCESS) goto cleanup;

                if (atype == ATTRCSYS) {

                    // Save the origin
                    (*feaCoordSystem)[i].origin[0] = reals[alen+0];
                    (*feaCoordSystem)[i].origin[1] = reals[alen+1];
                    (*feaCoordSystem)[i].origin[2] = reals[alen+2];

                    (*feaCoordSystem)[i].normal1[0] = reals[alen+3];
                    (*feaCoordSystem)[i].normal1[1] = reals[alen+4];
                    (*feaCoordSystem)[i].normal1[2] = reals[alen+5];

                    (*feaCoordSystem)[i].normal2[0] = reals[alen+6];
                    (*feaCoordSystem)[i].normal2[1] = reals[alen+7];
                    (*feaCoordSystem)[i].normal2[2] = reals[alen+8];

                    (*feaCoordSystem)[i].normal3[0] = reals[alen+9];
                    (*feaCoordSystem)[i].normal3[1] = reals[alen+10];
                    (*feaCoordSystem)[i].normal3[2] = reals[alen+11];

                    found = (int) true;

                    break;
                }
            } // End face loop

            if (found == (int) true)  break;

            // Loop through edges
            for (edge = 0; edge < numEdge; edge++) {

                status = EG_attributeRet(edges[edge], (*feaCoordSystem)[i].name, &atype, &alen, &ints, &reals, &string);

                if (status == EGADS_NOTFOUND) continue;
                if (status != EGADS_SUCCESS) goto cleanup;

                if (atype == ATTRCSYS) {

                    // Save the origin
                    (*feaCoordSystem)[i].origin[0] = reals[alen+0];
                    (*feaCoordSystem)[i].origin[1] = reals[alen+1];
                    (*feaCoordSystem)[i].origin[2] = reals[alen+2];

                    (*feaCoordSystem)[i].normal1[0] = reals[alen+3];
                    (*feaCoordSystem)[i].normal1[1] = reals[alen+4];
                    (*feaCoordSystem)[i].normal1[2] = reals[alen+5];

                    (*feaCoordSystem)[i].normal2[0] = reals[alen+6];
                    (*feaCoordSystem)[i].normal2[1] = reals[alen+7];
                    (*feaCoordSystem)[i].normal2[2] = reals[alen+8];

                    (*feaCoordSystem)[i].normal3[0] = reals[alen+9];
                    (*feaCoordSystem)[i].normal3[1] = reals[alen+10];
                    (*feaCoordSystem)[i].normal3[2] = reals[alen+11];

                    found = (int) true;

                    break;
                }

            } // End edge loop

            if (found == (int) true)  break;

            // Loop through nodes
            for (node = 0; node < numNode; node++) {

                status = EG_attributeRet(nodes[node], (*feaCoordSystem)[i].name, &atype, &alen, &ints, &reals, &string);

                if (status == EGADS_NOTFOUND) continue;
                if (status != EGADS_SUCCESS) goto cleanup;

                if (atype == ATTRCSYS) {

                    // Save the origin
                    (*feaCoordSystem)[i].origin[0] = reals[alen+0];
                    (*feaCoordSystem)[i].origin[1] = reals[alen+1];
                    (*feaCoordSystem)[i].origin[2] = reals[alen+2];

                    (*feaCoordSystem)[i].normal1[0] = reals[alen+3];
                    (*feaCoordSystem)[i].normal1[1] = reals[alen+4];
                    (*feaCoordSystem)[i].normal1[2] = reals[alen+5];

                    (*feaCoordSystem)[i].normal2[0] = reals[alen+6];
                    (*feaCoordSystem)[i].normal2[1] = reals[alen+7];
                    (*feaCoordSystem)[i].normal2[2] = reals[alen+8];

                    (*feaCoordSystem)[i].normal3[0] = reals[alen+9];
                    (*feaCoordSystem)[i].normal3[1] = reals[alen+10];
                    (*feaCoordSystem)[i].normal3[2] = reals[alen+11];

                    found = (int) true;

                    break;
                }
            } // End node loop

            EG_free(faces); faces = NULL;
            EG_free(edges); edges = NULL;
            EG_free(nodes); nodes = NULL;

            if (found == (int) true)  break;

        } // End body loop
    }

    status = CAPS_SUCCESS;

    cleanup:

        if (status != CAPS_SUCCESS) printf("\tError in fea_getCoordSystem = %d\n", status);

        EG_free(faces);
        EG_free(edges);
        EG_free(nodes);

        return status;
}

// Get the design equations from a capsTuple
int fea_getDesignEquation(int numEquationTuple,
                          capsTuple equationTuple[],
                          feaProblemStruct *feaProblem) {

    /*! \page feaDesignEquation FEA Design Equations
     * Structure for the design equation tuple  = ("DesignEquation Name", ["Value1", ... , "ValueN"]).
     * "DesignEquation Name" defines the reference name for the design equation being specified.
     *  This string will be used in the FEA input directly. The values "Value1", ... , "ValueN" are a
     *  list of strings containing the equation defintions.
     *  (see Section \ref tupleValueDesignEquation).
     *
     */
    int i, status;

    feaDesignEquationStruct *equation;

    // Ensure we are starting with no equations
    if (feaProblem->feaEquation != NULL) {
        for (i = 0; i < feaProblem->numEquation; i++) {
            status = destroy_feaDesignEquationStruct(&feaProblem->feaEquation[i]);
            if (status != CAPS_SUCCESS) return status;
        }
        EG_free(feaProblem->feaEquation);
    }
    feaProblem->numEquation = 0;
    feaProblem->feaEquation = NULL;

    printf("\nGetting Equations.......\n");

    feaProblem->numEquation = numEquationTuple;
    printf("\tNumber of Equations - %d\n", feaProblem->numEquation);

    if (feaProblem->numEquation > 0) {

        feaProblem->feaEquation = EG_alloc(
            feaProblem->numEquation * sizeof(feaDesignEquationStruct));
        if (feaProblem->feaEquation == NULL) {
            return EGADS_MALLOC;
        }

    } else {
        printf("\tNumber of equations in Analysis tuple is %d\n",
               feaProblem->numEquation);
        return CAPS_NOTFOUND;
    }

    // for each analysis equation tuple
    for (i = 0; i < feaProblem->numEquation; i++) {

        printf("\tDesign_Equation name - %s\n", equationTuple[i].name);

        equation = &feaProblem->feaEquation[i];

        // initiate equation structure
        status = initiate_feaDesignEquationStruct(equation);
        if (status != CAPS_SUCCESS) return status;

        // set name
        equation->name = EG_strdup(equationTuple[i].name);

        // set equation ID
        equation->equationID = i+1;

        /*! \page feaDesignEquation
        * \section tupleValueDesignEquation List of equation strings
        *
        * Each design equation tuple value is a list of strings containing the equation definitions
        * \if NASTRAN
        *  (eg. ["dispsum3(s1,s2,s3)=sum(s1,s2,s3)"]
        * \endif
        *
        */
        // set the equation array
        status = string_toStringDynamicArray(
            equationTuple[i].value,
            &equation->equationArraySize,
            &equation->equationArray);

        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Get the design table constants from a capsTuple
int fea_getDesignTable(int numConstantTuple,
                         capsTuple constantTuple[],
                         feaProblemStruct *feaProblem) {


    /*! \page feaDesignTable FEA Table Constants
     * Structure for the table constant tuple  = ("TableConstant Name", "Value").
     * "TableConstant Name" defines the reference name for the table constant being specified.
     *  This string will be used in the FEA input directly. The "Value" is the value of the
     *  table constant.
     *  \if NASTRAN
     *  In Nastran the TableConstant Name will be the LABLi used in the DTABLE input.
     *  For this reason the user should keep the length of this input to a minimum number of characters, ideally 7 or less.
     *
     *  - <c>DTABLE  LABL1   VALU1   LABL2   VALU2   LABL3   VALU3   -etc-   </c>
     *  \endif
     *
     */
    int i, status;

    feaDesignTableStruct *table = &feaProblem->feaDesignTable;

    // Ensure we are starting with no constants
    status = destroy_feaDesignTableStruct(table);
    if (status != CAPS_SUCCESS) return status;

    printf("\nGetting Design Table Constants.......\n");

    table->numConstant = numConstantTuple;
    printf("\tNumber of Design Table Constants - %d\n", table->numConstant);

    if (table->numConstant > 0) {

        table->constantLabel = EG_alloc(table->numConstant * sizeof(char *));
        if (table->constantLabel == NULL) {
            return EGADS_MALLOC;
        }

        table->constantValue = EG_alloc(table->numConstant * sizeof(double));
        if (table->constantValue == NULL) {
            return EGADS_MALLOC;
        }

    } else {
        printf("\tNumber of design table constants in Analysis tuple is %d\n",
               table->numConstant);
        return CAPS_NOTFOUND;
    }

    // for each analysis table constant tuple
    for (i = 0; i < table->numConstant; i++) {

        printf("\tDesign_Table - %s: %s\n",
               constantTuple[i].name, constantTuple[i].value);

        // set constant label
        // TODO: ensure label <= 8 chars
        table->constantLabel[i] = EG_strdup(constantTuple[i].name);

        // set constant value
        status = string_toDouble(constantTuple[i].value, &table->constantValue[i]);
        if (status != CAPS_SUCCESS) return status;
    }

    return CAPS_SUCCESS;
}

// Function used by fea_getDesignResponse to determine which nodes are in response group
static int _matchResponseNode(meshNodeStruct *node, void *responseIndex) {

    feaMeshDataStruct *feaData;

    if (node->analysisType == MeshStructure) {

        feaData = (feaMeshDataStruct *) node->analysisData;

        if (feaData->responseIndex == *((int *) responseIndex)) {
            return (int) true;
        }
    }
    return (int) false;
}

int fea_getDesignResponse(/*@unused@*/ void *aimInfo,
                          int numDesignResponseTuple,
                          capsTuple designResponseTuple[],
                          mapAttrToIndexStruct *responseMap,
                          feaProblemStruct *feaProblem) {

    /*! \page feaDesignResponse FEA Design Responses
     * Structure for the design response tuple  = ("DesignResponse Name", "Value").
     * "DesignResponse Name" defines the reference name for the design response being specified.
     *  This string will be used in the FEA input directly. The "Value" must be a JSON String dictionary
     *  (see Section \ref jsonStringDesignResponse).
     *  \if NASTRAN
     *  In Nastran the DesignResponse Name will be the LABEL used in the DRESP1 input.
     *  For this reason the user should keep the length of this input to a minimum number of characters, ideally 7 or less.
     *
     *  - <c>DRESP1  ID      LABEL   RTYPE   PTYPE   REGION   ATTA    ATTB    ATT1\n        ATT2    -etc-</c>
     *  \endif
     *
     */
    int i, status;

    int attrIndex;

    char *keyword;
    char *groupName = NULL;

    int numNode;
    meshNodeStruct **nodeSet = NULL;

    feaDesignResponseStruct *response;

    // Ensure we are starting with no design responses
    if (feaProblem->feaDesignResponse != NULL) {
        for (i = 0; i < feaProblem->numDesignResponse; i++) {
            status = destroy_feaDesignResponseStruct(&feaProblem->feaDesignResponse[i]);
            if (status != CAPS_SUCCESS) return status;
        }
        EG_free(feaProblem->feaDesignResponse);
    }
    feaProblem->numDesignResponse = 0;
    feaProblem->feaDesignResponse = NULL;

    printf("\nGetting Design Responses.......\n");

    feaProblem->numDesignResponse = numDesignResponseTuple;
    printf("\tNumber of Design Responses - %d\n", feaProblem->numDesignResponse);

    if (feaProblem->numDesignResponse > 0) {

        feaProblem->feaDesignResponse = EG_alloc(
            feaProblem->numDesignResponse * sizeof(feaDesignResponseStruct));
        if (feaProblem->feaDesignResponse == NULL) {
            return EGADS_MALLOC;
        }

    } else {
        printf("\tNumber of design responses in Analysis tuple is %d\n",
               feaProblem->numDesignResponse);
        return CAPS_NOTFOUND;
    }

    // for each analysis design response tuple
    for (i = 0; i < feaProblem->numDesignResponse; i++) {

        printf("\tDesign_Response name - %s\n", designResponseTuple[i].name);

        response = &feaProblem->feaDesignResponse[i];

        // initiate design response structure
        status = initiate_feaDesignResponseStruct(response);
        if (status != CAPS_SUCCESS) return status;

        // set name
        response->name = EG_strdup(designResponseTuple[i].name);

        // set response ID
        response->responseID = i+1;

        /*! \page feaDesignResponse
         * \section jsonStringDesignResponse JSON String Dictionary
         *
         * If "Value" is JSON string dictionary
         * \if NASTRAN
         *  (eg. "Value" = {"responseType": "DISP", groupName": "plate", "component": 3})
         * \endif
         *  the following keywords ( = default values) may be used:
         *
         */
        // make sure design response tuple value is json string
        if (!json_isDict(designResponseTuple[i].value)) {
            AIM_ERROR(aimInfo,
                "'Design_Response' tuple value must be a JSON dictionary");
            return CAPS_BADVALUE;
        }

        /*! \page feaDesignResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>responseType</B> </li> <br>
         *  Type of design sensitivity response. For options,
         *  see NASTRAN User Guide DRESP1 Design Sensitivity Response Attributes table.
         * </ul>
         * \endif
         */
        keyword = "responseType";
        status = json_getString(
            designResponseTuple[i].value, keyword, &response->responseType);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "Missing required entry \"responseType\" "
                        "in 'Design_Response' tuple value");
            return status;
        }

        // keyword = "propertyType";
        // status = json_getString(
        //     designResponseTuple[i].value, keyword, &response->propertyType);
        // if (status != CAPS_SUCCESS) {
        //     // optional
        // }

        /*! \page feaDesignResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>component = "(no default)"</B> </li> <br>
         *  Component flag.
         * </ul>
         * \endif
         */
        keyword = "component";
        status = json_getInteger(
            designResponseTuple[i].value, keyword, &response->component);
        if (status != CAPS_SUCCESS) {
            // optional
        }

         /*! \page feaDesignResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>attb = "(no default)"</B> </li> <br>
         *  ATTB Inputs.
         * </ul>
         * \endif
         */
        keyword = "attb";
        status = json_getInteger(
            designResponseTuple[i].value, keyword, &response->attb);
        if (status != CAPS_SUCCESS) {
            // optional
        }

        // TODO: temporary input
        keyword = "grid";
        status = json_getInteger(
            designResponseTuple[i].value, keyword, &response->gridID);
        if (status != CAPS_SUCCESS) {
            // optional
        }

        /*! \page feaDesignResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>groupName = "(no default)"</B> </li> <br>
         *  Defines the reference <c>capsGroup</c> for the node being specified for the response.
         * </ul>
         * \endif
         */
        keyword = "groupName";
        status = json_getString(
            designResponseTuple[i].value, keyword, &groupName);

        if (status == CAPS_SUCCESS) {
            AIM_NOTNULL(groupName, aimInfo, status);

            // Get the corresponding response index
            status = get_mapAttrToIndexIndex(responseMap, (const char *) groupName, &attrIndex);

            if (status == CAPS_NOTFOUND) {
                printf("\tName %s not found in attribute map!!!!\n", groupName);
                continue;
            } else if (status != CAPS_SUCCESS) return status;

            status = mesh_findNodes(&feaProblem->feaMesh,
                                    _matchResponseNode, (void *) &attrIndex,
                                    &numNode, &nodeSet);
            if (status != CAPS_SUCCESS) return status;

            if (numNode == 0) {
                AIM_ERROR(aimInfo, "No node found for capsGroup %s", groupName);
                return CAPS_NOTFOUND;
            }
            else if (numNode > 1) { // TODO: would there ever be more than 1 node expected?
                PRINT_WARNING("More than 1 node found for capsGroup %s"
                            "... using first matching node.",
                            groupName);
            }

            response->gridID = nodeSet[0]->nodeID;

            EG_free(nodeSet);
            nodeSet = NULL;

        }
        else {
            // required
            // AIM_ERROR(aimInfo, "Missing required entry \"groupName\" "
            //             "in 'Design_Response' tuple value");
            // return status;
        }

        if (groupName != NULL) {
            EG_free(groupName);
            groupName = NULL;
        }

    }

    status = CAPS_SUCCESS;
cleanup:
    return status;
}

int fea_getDesignEquationResponse(int numDesignEquationResponseTuple,
                            capsTuple designEquationResponseTuple[],
                            feaProblemStruct *feaProblem) {

    /*! \page feaDesignEquationResponse FEA Design Equation Responses
     * Structure for the design equation response tuple  = ("DesignEquationResponse Name", "Value").
     * "DesignEquationResponse Name" defines the reference name for the design equation response being specified.
     *  This string will be used in the FEA input directly. The "Value" must be a JSON String dictionary
     *  (see Section \ref jsonStringDesignEquationResponse).
     *  \if NASTRAN
     *  In Nastran the DesignEquationResponse Name will be the LABEL used in the DRESP2 input.
     *  For this reason the user should keep the length of this input to a minimum number of characters, ideally 7 or less.
     *
     *  - <c>DRESP2  ID      LABEL   EQID    REGION ...</c>
     *  \endif
     *
     */
    int i, status;

    char *keyword;

    feaDesignEquationResponseStruct *equationResponse;

    // Ensure we are starting with no design equation responses
    if (feaProblem->feaEquationResponse != NULL) {
        for (i = 0; i < feaProblem->numEquationResponse; i++) {
            status = destroy_feaDesignEquationResponseStruct(&feaProblem->feaEquationResponse[i]);
            if (status != CAPS_SUCCESS) return status;
        }
        EG_free(feaProblem->feaEquationResponse);
    }
    feaProblem->numEquationResponse = 0;
    feaProblem->feaEquationResponse = NULL;

    printf("\nGetting Design Equation Responses.......\n");

    feaProblem->numEquationResponse = numDesignEquationResponseTuple;
    printf("\tNumber of Design Equation Responses - %d\n", feaProblem->numEquationResponse);

    if (feaProblem->numEquationResponse > 0) {

        feaProblem->feaEquationResponse = EG_alloc(
            feaProblem->numEquationResponse * sizeof(feaDesignEquationResponseStruct));
        if (feaProblem->feaEquationResponse == NULL) {
            return EGADS_MALLOC;
        }

    } else {
        printf("\tNumber of design equation responses in Analysis tuple is %d\n",
               feaProblem->numEquationResponse);
        return CAPS_NOTFOUND;
    }

    // for each analysis design equation response tuple
    for (i = 0; i < feaProblem->numEquationResponse; i++) {

        printf("\tDesign_Equation_Response name - %s\n", designEquationResponseTuple[i].name);

        equationResponse = &feaProblem->feaEquationResponse[i];

        // initiate design response structure
        status = initiate_feaDesignEquationResponseStruct(equationResponse);
        if (status != CAPS_SUCCESS) return status;

        // set name
        equationResponse->name = EG_strdup(designEquationResponseTuple[i].name);

        // set equation response ID
        equationResponse->equationResponseID = i+1;

        /*! \page feaDesignEquationResponse
         * \section jsonStringDesignEquationResponse JSON String Dictionary
         *
         * If "Value" is JSON string dictionary
         * \if NASTRAN
         *  (eg. "Value" = {"equation": "EQ1", "constant": ["PI", "YM", "L"]})
         * \endif
         *  the following keywords ( = default values) may be used:
         *
         */
        // make sure design response tuple value is json string
        if (!json_isDict(designEquationResponseTuple[i].value)) {
            PRINT_ERROR(
                "'Design_Equation_Response' tuple value must be a JSON dictionary");
            return CAPS_BADVALUE;
        }

        /*! \page feaDesignEquationResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>equation</B> </li> <br>
         *  The name of the equation referenced by this equation response.
         * </ul>
         * \endif
         */
        keyword = "equation";
        status = json_getString(
            designEquationResponseTuple[i].value, keyword, &equationResponse->equationName);
        if (status != CAPS_SUCCESS) {
            PRINT_ERROR("Missing required entry \"equation\" "
                        "in 'Design_Equation_Response' tuple value");
            return status;
        }

        /*! \page feaDesignEquationResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>variable = "(no default)"</B> </li> <br>
         *  Single or list of names of design variable equation parameters.
         * </ul>
         * \endif
         */
        keyword = "variable";
        status = json_getStringDynamicArray(
            designEquationResponseTuple[i].value, keyword,
            &equationResponse->numDesignVariable, &equationResponse->designVariableNameSet);
        if (status != CAPS_SUCCESS) {
            // optional
        }

        /*! \page feaDesignEquationResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>constant = "(no default)"</B> </li> <br>
         *  Single or list of names of table constant equation parameters.
         * </ul>
         * \endif
         */
        keyword = "constant";
        status = json_getStringDynamicArray(
            designEquationResponseTuple[i].value, keyword,
            &equationResponse->numConstant, &equationResponse->constantLabelSet);
        if (status != CAPS_SUCCESS) {
            // optional
        }

        /*! \page feaDesignEquationResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>response = "(no default)"</B> </li> <br>
         *  Single or list of names of design response equation parameters.
         * </ul>
         * \endif
         */
        keyword = "response";
        status = json_getStringDynamicArray(
            designEquationResponseTuple[i].value, keyword,
            &equationResponse->numResponse, &equationResponse->responseNameSet);
        if (status != CAPS_SUCCESS) {
            // optional
        }

        /*! \page feaDesignEquationResponse
         *
         * \if NASTRAN
         * <ul>
         *  <li> <B>equationResponse = "(no default)"</B> </li> <br>
         *  Single or list of names of design equation response equation parameters.
         * </ul>
         * \endif
         */
        keyword = "equationResponse";
        status = json_getStringDynamicArray(
            designEquationResponseTuple[i].value, keyword,
            &equationResponse->numEquationResponse, &equationResponse->equationResponseNameSet);
        if (status != CAPS_SUCCESS) {
            // optional
        }
    }

    return CAPS_SUCCESS;
}

// Get the design optimization parameters from a capsTuple
int fea_getDesignOptParam(int numParamTuple,
                          capsTuple paramTuple[],
                          feaProblemStruct *feaProblem) {

    /*! \page feaDesignOptParam FEA Design Optimization Parameters
     * Structure for the design optimization parameter tuple  = ("DesignOptParam Name", "Value").
     * "DesignOptParam Name" defines the reference name for the design optimization parameter being specified.
     *  This string will be used in the FEA input directly. The "Value" is the value of the
     *  design optimization parameter.
     *  \if NASTRAN
     *  In Nastran the DesignOptParam Name will be the PARAMi used in the DOPTPRM input.
     *  For this reason the user should keep the length of this input to a minimum number of characters, ideally 7 or less.
     *
     *  - <c>DOPTPRM PARAM1  VAL1    PARAM2  VAL2    PARAM3  VAL3       -etc-   </c>
     *  \endif
     *
     */
    int i, status;

    int *paramInt = NULL;
    double *paramDouble = NULL;

    feaDesignOptParamStruct *table = &feaProblem->feaDesignOptParam;

    // Ensure we are starting with no params
    status = destroy_feaDesignOptParamStruct(table);
    if (status != CAPS_SUCCESS) return status;

    printf("\nGetting Design Optimization Parameters.......\n");

    table->numParam = numParamTuple;
    printf("\tNumber of Design Optimization Parameters - %d\n", table->numParam);

    if (table->numParam > 0) {

        table->paramLabel = EG_alloc(table->numParam * sizeof(char *));
        if (table->paramLabel == NULL) {
            return EGADS_MALLOC;
        }

        table->paramValue = EG_alloc(table->numParam * sizeof(void *));
        if (table->paramValue == NULL) {
            return EGADS_MALLOC;
        }

        table->paramType = EG_alloc(table->numParam * sizeof(int));
        if (table->paramType == NULL) {
            return EGADS_MALLOC;
        }

    } else {
        printf("\tNumber of design optimization parameters in Analysis tuple is %d\n",
               table->numParam);
        return CAPS_NOTFOUND;
    }

    // for each analysis design optimization parameter tuple
    for (i = 0; i < table->numParam; i++) {

        printf("\tDesign_Opt_Param - %s: %s\n",
               paramTuple[i].name, paramTuple[i].value);

        // set param label
        // TODO: ensure label <= 8 chars
        table->paramLabel[i] = EG_strdup(paramTuple[i].name);

        // set param value

        // if param value is real
        if (strchr(paramTuple[i].value, '.')
            || strchr(paramTuple[i].value, 'e')
            || strchr(paramTuple[i].value, 'E')) {

            paramDouble = EG_alloc(sizeof(double));
            if (paramDouble == NULL) return EGADS_MALLOC;

            status = string_toDouble(paramTuple[i].value, paramDouble);
            if (status != CAPS_SUCCESS) {
                EG_free(paramDouble);
                return status;
            }

            table->paramType[i] = Double;
            table->paramValue[i] = (void*) paramDouble;
        }
        // else param is integer
        else {

            paramInt = EG_alloc(sizeof(int));
            if (paramInt == NULL) return EGADS_MALLOC;

            status = string_toInteger(paramTuple[i].value, paramInt);
            if (status != CAPS_SUCCESS) {
                EG_free(paramInt);
                return status;
            }

            table->paramType[i] = Integer;
            table->paramValue[i] = (void*) paramInt;
        }
    }

    return CAPS_SUCCESS;
}

// Function used by fea_getAeroReference to determine which nodes are in reference group
static int _matchReferenceNode(meshNodeStruct *node, void *referenceIndex) {

    feaMeshDataStruct *feaData;

    if (node->analysisType == MeshStructure) {

        feaData = (feaMeshDataStruct *) node->analysisData;

        if (feaData->referenceIndex == *((int *) referenceIndex)) {
            return (int) true;
        }
    }
    return (int) false;
}

// Get aero reference information
int fea_getAeroReference(void *aimInfo,
                         int numAeroRefTuple,
                         capsTuple aeroRefTuple[],
                         mapAttrToIndexStruct *attrMap,
                         feaProblemStruct *feaProblem) {

    /*! \page feaAeroReference FEA Aerodynamic References
     * Tuple of the aerodynamic reference input
     *  (see Section \ref jsonStringAeroReference).
     *
     */
    int status, attrIndex = 0;
    int i;

    char *refNodeName = NULL;
    char *keyValue = NULL, *keyWord = NULL;

    int numNodes;
    meshNodeStruct **nodes = NULL;

    printf("\nGetting Aero Reference Parameters.......\n");

    /*! \page feaAeroReference
     * \section jsonStringAeroReference JSON String Dictionary
     *
     * The following keywords ( = default values) may be used:
     *
     */
    for (i = 0; i < numAeroRefTuple; i++) {

        /*! \page feaAeroReference
         *
         * \if ASTROS
         * <ul>
         *  <li> <B>referenceNode</B> </li> <br>
         *  Defines the reference <c>capsReference</c> for the node to be used for stability derivative calculations.
         * </ul>
         * \endif
         */
        keyWord = "referenceNode";
        status = search_jsonDictionary(aeroRefTuple[i].value, keyWord, &refNodeName);
        if (status == CAPS_SUCCESS) {
            AIM_NOTNULL(refNodeName, aimInfo, status);

            // find node with capsReference == refNodeName
            status = get_mapAttrToIndexIndex(attrMap, refNodeName, &attrIndex);
            if (status != CAPS_SUCCESS) {
                PRINT_ERROR("capsGroup name %s not found in attribute to index map", refNodeName);
                goto cleanup;
            }

            status = mesh_findNodes(&feaProblem->feaMesh, _matchReferenceNode, &attrIndex, &numNodes, &nodes);
            if (status == CAPS_NOTFOUND) {
                PRINT_ERROR("No node found with capsReference name %s", refNodeName);
            }
            if (status != CAPS_SUCCESS) {
                goto cleanup;
            }

            // set refGridID to the nodeID
            feaProblem->feaAeroRef.refGridID = nodes[0]->nodeID;
        }

        /*! \page feaAeroReference
         *
         * \if ASTROS
         * <ul>
         *  <li> <B>referenceVelocity</B> </li> <br>
         *  Defines the reference <c>capsReference</c> for the node to be used for stability derivative calculations.
         * </ul>
         * \endif
         */
        keyWord = "referenceVelocity";
        if (strcmp(aeroRefTuple[i].name, keyWord) == 0) {
            status = string_toDouble(aeroRefTuple[i].value, &feaProblem->feaAeroRef.refVelocity);
            if (keyValue != NULL) {
                EG_free(keyValue);
                keyValue = NULL;
            }
            if (status != CAPS_SUCCESS) return status;
        }

        /*! \page feaAeroReference
         *
         * \if ASTROS
         * <ul>
         *  <li> <B>referenceDensity</B> </li> <br>
         *  Defines the reference <c>capsReference</c> for the node to be used for stability derivative calculations.
         * </ul>
         * \endif
         */
        keyWord = "referenceDensity";
        if (strcmp(aeroRefTuple[i].name, keyWord) == 0) {
            status = string_toDouble(aeroRefTuple[i].value, &feaProblem->feaAeroRef.refDensity);
            if (keyValue != NULL) {
                EG_free(keyValue);
                keyValue = NULL;
            }
            if (status != CAPS_SUCCESS) return status;
        }
    }

    status = CAPS_SUCCESS;

    cleanup:

        AIM_FREE(refNodeName);
        AIM_FREE(nodes);

        return status;
}

// Find feaPropertyStructs with given names in feaProblem
int fea_findPropertiesByNames(feaProblemStruct *feaProblem,
                             int numPropertyNames,
                             char **propertyNames,
                             int *numProperties,
                             feaPropertyStruct ***properties) {
    int i, status;

    int numFound;
    feaPropertyStruct **propertiesFound, *property;

    propertiesFound = EG_alloc(sizeof(feaPropertyStruct *) * numPropertyNames);

    numFound = 0;
    for (i = 0; i < feaProblem->numProperty; i++) {

        property = &feaProblem->feaProperty[i];

        if (string_isInArray(property->name, numPropertyNames, propertyNames)) {

            propertiesFound[numFound++] = property;

            if (numFound == numPropertyNames) {
                // found all
                break;
            }
        }
    }

    if (numFound < numPropertyNames) {
        propertiesFound = EG_reall(propertiesFound, sizeof(feaPropertyStruct *) * numFound);
        status = CAPS_NOTFOUND;
    }
    else {
        status = CAPS_SUCCESS;
    }

    *numProperties = numFound;
    *properties = propertiesFound;

    return status;
}

// Find feaMaterialStructs with given names in feaProblem
int fea_findMaterialsByNames(feaProblemStruct *feaProblem,
                             int numMaterialNames,
                             char **materialNames,
                             int *numMaterials,
                             feaMaterialStruct ***materials) {
    int i, status;

    int numFound;
    feaMaterialStruct **materialsFound, *material;

    materialsFound = EG_alloc(sizeof(feaMaterialStruct *) * numMaterialNames);

    numFound = 0;
    for (i = 0; i < feaProblem->numMaterial; i++) {

        material = &feaProblem->feaMaterial[i];

        if (string_isInArray(material->name, numMaterialNames, materialNames)) {

            materialsFound[numFound++] = material;

            if (numFound == numMaterialNames) {
                // found all
                break;
            }
        }
    }

    if (numFound < numMaterialNames) {
        materialsFound = EG_reall(materialsFound, sizeof(feaMaterialStruct *) * numFound);
        status = CAPS_NOTFOUND;
    }
    else {
        status = CAPS_SUCCESS;
    }

    *numMaterials = numFound;
    *materials = materialsFound;

    return status;
}

// Find feaDesignVariableStructs with given names in feaProblem
int fea_findDesignVariablesByNames(const feaProblemStruct *feaProblem,
                                   int numDesignVariableNames,
                                   char **designVariableNames,
                                   int *numDesignVariables,
                                   feaDesignVariableStruct ***designVariables) {
    int i, status;

    int numFound;
    feaDesignVariableStruct **designVariablesFound=NULL, *designVariable=NULL;

    designVariablesFound = (feaDesignVariableStruct **) EG_alloc(numDesignVariableNames * sizeof(feaDesignVariableStruct *));

    numFound = 0;
    for (i = 0; i < feaProblem->numDesignVariable; i++) {

        designVariable = &feaProblem->feaDesignVariable[i];

        if (string_isInArray(designVariable->name, numDesignVariableNames, designVariableNames)) {

            designVariablesFound[numFound++] = designVariable;

            if (numFound == numDesignVariableNames) {
                // found all
                break;
            }
        }
    }

    if (numFound < numDesignVariableNames) {
        designVariablesFound = EG_reall(designVariablesFound, numFound * sizeof(feaDesignVariableStruct *));
        status = CAPS_NOTFOUND;
    } else {
        status = CAPS_SUCCESS;
    }

    *numDesignVariables = numFound;
    *designVariables = designVariablesFound;

    return status;
}

// Find feaDesignResponseStructs with given names in feaProblem
int fea_findDesignResponsesByNames(const feaProblemStruct *feaProblem,
                                   int numDesignResponseNames,
                                   char **designResponseNames,
                                   int *numDesignResponses,
                                   feaDesignResponseStruct ***designResponses) {
    int i, status;

    int numFound;
    feaDesignResponseStruct **designResponsesFound, *designResponse;

    designResponsesFound = EG_alloc(sizeof(feaDesignResponseStruct *) * numDesignResponseNames);

    numFound = 0;
    for (i = 0; i < feaProblem->numDesignResponse; i++) {

        designResponse = &feaProblem->feaDesignResponse[i];

        if (string_isInArray(designResponse->name, numDesignResponseNames, designResponseNames)) {

            designResponsesFound[numFound++] = designResponse;

            if (numFound == numDesignResponseNames) {
                // found all
                break;
            }
        }
    }

    if (numFound < numDesignResponseNames) {
        designResponsesFound = EG_reall(designResponsesFound, sizeof(feaDesignResponseStruct *) * numFound);
        status = CAPS_NOTFOUND;
    }
    else {
        status = CAPS_SUCCESS;
    }

    *numDesignResponses = numFound;
    *designResponses = designResponsesFound;

    return status;
}

// Find feaDesignEquationResponseStructs with given names in feaProblem
int fea_findEquationResponsesByNames(const feaProblemStruct *feaProblem,
                                     int numEquationResponseNames,
                                     char *const* equationResponseNames,
                                     int *numEquationResponses,
                                     feaDesignEquationResponseStruct ***equationResponses) {
    int i, status;

    int numFound;
    feaDesignEquationResponseStruct **equationResponsesFound, *equationResponse;

    equationResponsesFound = EG_alloc(sizeof(feaDesignEquationResponseStruct *) * numEquationResponseNames);

    numFound = 0;
    for (i = 0; i < feaProblem->numEquationResponse; i++) {

        equationResponse = &feaProblem->feaEquationResponse[i];

        if (string_isInArray(equationResponse->name, numEquationResponseNames, equationResponseNames)) {

            equationResponsesFound[numFound++] = equationResponse;

            if (numFound == numEquationResponseNames) {
                // found all
                break;
            }
        }
    }

    if (numFound < numEquationResponseNames) {
        equationResponsesFound = EG_reall(equationResponsesFound, sizeof(feaDesignEquationResponseStruct *) * numFound);
        status = CAPS_NOTFOUND;
    }
    else {
        status = CAPS_SUCCESS;
    }

    *numEquationResponses = numFound;
    *equationResponses = equationResponsesFound;

    return status;
}

// Find feaDesignEquationStruct with given equationName in feaProblem
int fea_findEquationByName(const feaProblemStruct *feaProblem, char *equationName, feaDesignEquationStruct **equation) {

    int i;

    for (i = 0; i < feaProblem->numEquation; i++) {

        if (strcmp(feaProblem->feaEquation[i].name, equationName) == 0) {
            *equation = &feaProblem->feaEquation[i];
            return CAPS_SUCCESS;
        }
    }

    return CAPS_NOTFOUND;
}

// Find vlmControlStruct with given controlSurfName in feaProble
int fea_findControlSurfaceByName(const feaProblemStruct *feaProblem, char *controlSurfName, vlmControlStruct **controlSurf) {

    int isection, icontrol;

    vlmSurfaceStruct *surface;
    vlmSectionStruct *section;

    surface = &feaProblem->feaAero->vlmSurface;

    for (isection = 0; isection < surface->numSection; isection++) {

        section = &surface->vlmSection[isection];

        for (icontrol = 0; icontrol < section->numControl; icontrol++) {

            if (strcmp(section->vlmControl[icontrol].name, controlSurfName) == 0) {
                *controlSurf = &section->vlmControl[icontrol];
                return CAPS_SUCCESS;
            }
        }
    }

    return CAPS_NOTFOUND;
}

// Find feaDesignConstraintStruct with given constraintType in feaProblem
int fea_findDesignConstraintByType(const feaProblemStruct *feaProblem, feaDesignConstraintTypeEnum constraintType, feaDesignConstraintStruct **designConstraint) {

    int i;

    for (i = 0; i < feaProblem->numDesignConstraint; i++) {

        if (feaProblem->feaDesignConstraint[i].designConstraintType == constraintType) {
            *designConstraint = &feaProblem->feaDesignConstraint[i];
            return CAPS_SUCCESS;
        }
    }

    return CAPS_NOTFOUND;
}

static inline int _feaDesignVariableRelation_linkDesignVariable(feaDesignVariableRelationStruct *relation,
                                                  feaDesignVariableStruct *desvar,
                                                  int desvarIndex) {
    if (relation->designVariableSet == NULL) {
        relation->designVariableSet = (
            (feaDesignVariableStruct **)
            EG_alloc(sizeof(feaDesignVariableStruct *) * relation->numDesignVariable)
        );
    }
    if (relation->designVariableSet == NULL) {
        return EGADS_MALLOC;
    }

    relation->designVariableSet[desvarIndex] = desvar;

    return CAPS_SUCCESS;
}

static inline int _feaDesignVariable_linkDesignVariableRelation(feaDesignVariableStruct *desvar,
                                                  feaDesignVariableRelationStruct *relation) {

    if (desvar->numRelation == 0) {
        desvar->relationSet = (
            (feaDesignVariableRelationStruct **)
            EG_alloc(sizeof(feaDesignVariableRelationStruct *))
        );
    }
    else {
        desvar->relationSet = (
            (feaDesignVariableRelationStruct **)
            EG_reall(desvar->relationSet, sizeof(feaDesignVariableRelationStruct *) * (desvar->numRelation + 1))
        );
    }
    if (desvar->relationSet == NULL) {
        return EGADS_MALLOC;
    }

    desvar->relationSet[desvar->numRelation] = relation;
    desvar->numRelation++;

    return CAPS_SUCCESS;
}

// Populate the feaDesignVariable.relationSet and feaDesignVariableRelation.variableSet members
// in all feaDesignVariables and feaDesignVariableRelations
int fea_linkDesignVariablesAndRelations(const feaProblemStruct *feaProblem) {

    int status;
    int idvGlobal, idv, irel;

    char *desvarName = NULL, **desvarNameSet = NULL;

    feaDesignVariableStruct *desvar = NULL;
    feaDesignVariableRelationStruct *desvarRelation = NULL;

    for (irel = 0; irel < feaProblem->numDesignVariableRelation; irel++) {

        desvarRelation = &feaProblem->feaDesignVariableRelation[irel];
        desvarNameSet = desvarRelation->designVariableNameSet;

        for (idv = 0; idv < desvarRelation->numDesignVariable; idv++) {

            desvarName = desvarNameSet[idv];

            for (idvGlobal = 0; idvGlobal < feaProblem->numDesignVariable; idvGlobal++) {

                desvar = &feaProblem->feaDesignVariable[idvGlobal];

                if (strcmp(desvar->name, desvarName) == 0) {

                    status = _feaDesignVariableRelation_linkDesignVariable(desvarRelation, desvar, idv);
                    if (status != CAPS_SUCCESS) return status;

                    status = _feaDesignVariable_linkDesignVariableRelation(desvar, desvarRelation);
                    if (status != CAPS_SUCCESS) return status;
                    break;
                }
            }
        }
    }

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaProblem in the feaProblemStruct structure format
int initiate_feaProblemStruct(feaProblemStruct *feaProblem) {

    int status = 0;

    // Problem analysis
    feaProblem->numAnalysis = 0;
    feaProblem->feaAnalysis = NULL;

    // Materials
    feaProblem->numMaterial = 0;
    feaProblem->feaMaterial = NULL;

    // Properties
    feaProblem->numProperty = 0;
    feaProblem->feaProperty = NULL;

    // Constraints
    feaProblem->numConstraint = 0;
    feaProblem->feaConstraint = NULL;

    // Supports
    feaProblem->numSupport = 0;
    feaProblem->feaSupport = NULL;

    // Loads
    feaProblem->numLoad = 0;
    feaProblem->feaLoad = NULL;

    // Connections
    feaProblem->numConnect = 0;
    feaProblem->feaConnect = NULL;

    // Mesh
    status = initiate_meshStruct(&feaProblem->feaMesh);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Output formatting
    status = initiate_feaFileFormatStruct(&feaProblem->feaFileFormat);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Optimization - Design Variables
    feaProblem->numDesignVariable = 0;
    feaProblem->feaDesignVariable = NULL;

    // Optimization - Design Variable Relations
    feaProblem->numDesignVariableRelation = 0;
    feaProblem->feaDesignVariableRelation = NULL;

    // Optimization - Design Constraints
    feaProblem->numDesignConstraint = 0;
    feaProblem->feaDesignConstraint = NULL;

    // Optimization - Design Constraints
    feaProblem->numMassIncrement = 0;
    feaProblem->feaMassIncrement = NULL;

    // Optimization - Equations
    feaProblem->numEquation = 0;
    feaProblem->feaEquation = NULL;

    // Optimization - Table Constants
    status = initiate_feaDesignTableStruct(&feaProblem->feaDesignTable);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Optimization - Design Optimization Parameters
    status = initiate_feaDesignOptParamStruct(&feaProblem->feaDesignOptParam);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Optimization - Design Sensitivity Response Quantities
    feaProblem->numDesignResponse = 0;
    feaProblem->feaDesignResponse = NULL;

    // Optimization - Design Sensitivity Equation Response Quantities
    feaProblem->numEquationResponse = 0;
    feaProblem->feaEquationResponse = NULL;

    // Coordinate Systems
    feaProblem->numCoordSystem = 0;
    feaProblem->feaCoordSystem = NULL;

    // Aerodynamics
    feaProblem->numAero = 0;
    feaProblem->feaAero = NULL;
    status = initiate_feaAeroRefStruct(&feaProblem->feaAeroRef);
    if (status != CAPS_SUCCESS) goto cleanup;

    aim_initMeshRef(&feaProblem->meshRefObj, aimUnknownMeshType);
    feaProblem->meshRefIn = NULL;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS) printf("Error: Status %d during initiate_feaProblemStruct!\n", status);

    return status;
}

// Destroy (0 out all values and NULL all pointers) of feaProblem in the feaProblemStruct structure format
int destroy_feaProblemStruct(feaProblemStruct *feaProblem) {

    int i; // Indexing

    int status = 0; // Status return

    // Analysis
    if (feaProblem->feaAnalysis != NULL) {
        for (i = 0; i < feaProblem->numAnalysis; i++) {
            status = destroy_feaAnalysisStruct(&feaProblem->feaAnalysis[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaAnalysisStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaAnalysis);
    feaProblem->numAnalysis = 0;

    // Materials
    if (feaProblem->feaMaterial != NULL) {
        for (i = 0; i < feaProblem->numMaterial; i++) {
            status = destroy_feaMaterialStruct(&feaProblem->feaMaterial[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaMaterialStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaMaterial);
    feaProblem->numMaterial = 0;

    // Properties
    if (feaProblem->feaProperty != NULL) {
        for (i = 0; i < feaProblem->numProperty; i++) {
            status = destroy_feaPropertyStruct(&feaProblem->feaProperty[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaPropertyStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaProperty);
    feaProblem->numProperty = 0;

    // Constraints
    if (feaProblem->feaConstraint != NULL) {
        for (i = 0; i < feaProblem->numConstraint; i++) {
            status = destroy_feaConstraintStruct(&feaProblem->feaConstraint[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaConstraintStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaConstraint);
    feaProblem->numConstraint = 0;

    // Supports
    if (feaProblem->feaSupport != NULL) {
        for (i = 0; i < feaProblem->numSupport; i++) {
            status = destroy_feaSupportStruct(&feaProblem->feaSupport[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaSupportStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaSupport);
    feaProblem->numSupport = 0;

    // Loads
    if (feaProblem->feaLoad != NULL) {
        for (i = 0; i < feaProblem->numLoad; i++) {
            status = destroy_feaLoadStruct(&feaProblem->feaLoad[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaLoadStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaLoad);
    feaProblem->numLoad = 0;

    // Connections
    if (feaProblem->feaConnect != NULL) {
        for (i = 0; i < feaProblem->numConnect; i++) {
            status = destroy_feaConnectionStruct(&feaProblem->feaConnect[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaConnectStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaConnect);
    feaProblem->numConnect = 0;

    // Mesh
    status = destroy_meshStruct(&feaProblem->feaMesh);
    if (status != CAPS_SUCCESS) printf("Status %d during destroy_meshStruct\n", status);

    // Output formatting
    status = destroy_feaFileFormatStruct(&feaProblem->feaFileFormat);
    if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaFileFormatStruct\n", status);


    if (feaProblem->feaDesignVariable != NULL) {
        for (i = 0; i < feaProblem->numDesignVariable; i++) {
            status = destroy_feaDesignVariableStruct(&feaProblem->feaDesignVariable[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignVariableStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaDesignVariable);
    feaProblem->numDesignVariable = 0;

    // Optimization - design variable relations
    if (feaProblem->feaDesignVariableRelation != NULL) {
        for (i = 0; i < feaProblem->numDesignVariableRelation; i++) {
            status = destroy_feaDesignVariableRelationStruct(&feaProblem->feaDesignVariableRelation[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignVariableRelationStruct\n", status);
        }
        AIM_FREE(feaProblem->feaDesignVariableRelation);
    }
    feaProblem->numDesignVariableRelation = 0;

    // Optimization - design constraint
    if (feaProblem->feaDesignConstraint != NULL) {
        for (i = 0; i < feaProblem->numDesignConstraint; i++) {
            status = destroy_feaDesignConstraintStruct(&feaProblem->feaDesignConstraint[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignConstraintStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaDesignConstraint);
    feaProblem->numDesignConstraint = 0;

    // Mass Increment
    if (feaProblem->feaMassIncrement != NULL) {
        for (i = 0; i < feaProblem->numMassIncrement; i++) {
            status = destroy_feaMassIncrementStruct(&feaProblem->feaMassIncrement[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaMassIncrementStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaMassIncrement);
    feaProblem->numMassIncrement = 0;

    // Optimization - Equations
    if (feaProblem->feaEquation != NULL) {
        for (i = 0; i < feaProblem->numEquation; i++) {
            status = destroy_feaDesignEquationStruct(&feaProblem->feaEquation[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignEquationStruct\n", status);
        }
        AIM_FREE(feaProblem->feaEquation);
    }
    feaProblem->numEquation = 0;

    // Optimization - Table Constants
    status = destroy_feaDesignTableStruct(&feaProblem->feaDesignTable);
    if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignTableStruct\n", status);

    // Optimization - Design Optimization Parameters
    status = destroy_feaDesignOptParamStruct(&feaProblem->feaDesignOptParam);
    if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignOptParamStruct\n", status);

    // Optimization - Design Sensitivity Response Quantities
    if (feaProblem->feaDesignResponse != NULL) {
        for (i = 0; i < feaProblem->numDesignResponse; i++) {
            status = destroy_feaDesignResponseStruct(&feaProblem->feaDesignResponse[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaDesignResponseStruct\n", status);
        }
        AIM_FREE(feaProblem->feaDesignResponse);
    }
    feaProblem->numDesignResponse = 0;

    // Optimization - Design Sensitivity Equation Response Quantities
    if (feaProblem->feaEquationResponse != NULL) {
        for (i = 0; i < feaProblem->numEquationResponse; i++) {
            status = destroy_feaDesignEquationResponseStruct(&feaProblem->feaEquationResponse[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaEquationResponseStruct\n", status);
        }
        AIM_FREE(feaProblem->feaEquationResponse);
    }
    feaProblem->numEquationResponse = 0;

    // Coordinate Systems
    if (feaProblem->feaCoordSystem != NULL) {

        for (i = 0; i < feaProblem->numCoordSystem; i++) {
            status = destroy_feaCoordSystemStruct(&feaProblem->feaCoordSystem[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaCoordSystemStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaCoordSystem);
    feaProblem->numCoordSystem = 0;

    // Aerodynamics
    if (feaProblem->feaAero != NULL) {

        for (i = 0; i < feaProblem->numAero; i++) {
            status = destroy_feaAeroStruct(&feaProblem->feaAero[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaAeroStruct\n", status);
        }
    }
    AIM_FREE(feaProblem->feaAero);
    feaProblem->numAero = 0;

    (void) destroy_feaAeroRefStruct(&feaProblem->feaAeroRef);

    aim_freeMeshRef(&feaProblem->meshRefObj);
    feaProblem->meshRefIn = NULL;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaFileFormat in the feaFileFormatStruct structure format
int initiate_feaFileFormatStruct(feaFileFormatStruct *feaFileFormat) {

    feaFileFormat->fileType = SmallField;

    feaFileFormat->gridFileType = LargeField;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaFileFormat in the feaFileFormatStruct structure format
int destroy_feaFileFormatStruct(feaFileFormatStruct *feaFileFormat) {

    feaFileFormat->fileType = SmallField;

    feaFileFormat->gridFileType = LargeField;

    return CAPS_SUCCESS;
}

// Initiate (0 out all values and NULL all pointers) of feaFileFormat in the feaFileFormatStruct structure format
int initiate_feaSolFileStruct(feaSolFileStruct *feaSolFile) {

    feaSolFile->fileForm = UnknownSolFileForm;

    feaSolFile->filename = NULL;

    return CAPS_SUCCESS;
}

// Destroy (0 out all values and NULL all pointers) of feaFileFormat in the feaFileFormatStruct structure format
int destroy_feaSolFileStruct(feaSolFileStruct *feaSolFile) {

    feaSolFile->fileForm = UnknownSolFileForm;

    AIM_FREE(feaSolFile->filename);

    return CAPS_SUCCESS;
}

// Transfer external pressure from the discrObj into the feaLoad structure
int fea_transferExternalPressure(void *aimInfo, const meshStruct *feaMesh, feaLoadStruct *feaLoad) {

    // [in/out] feaLoad
    // [in] feaMesh
    // [in] aimInfo

    int status; // Function status return

    int i, j, bIndex; // Indexing

    // Variables used in global node mapping
    int globalNodeID;
    int nodeIndex[4], transferIndex[4], elementID, elementIndex, elementCount;

    // Discrete data transfer variables
    capsDiscr *dataTransferDiscreteObj;
    char **transferName = NULL;
    int numTransferName, transferNameIndex;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int numDataTransferElement = 0, discElements;
    int dataTransferRank;
    double *dataTransferData;
    char   *units;

    feaLoad->numElementID = 0;
    AIM_FREE(feaLoad->pressureMultiDistributeForce);
    AIM_FREE(feaLoad->elementIDSet);

    //See if we have data transfer information
    status = aim_getBounds(aimInfo, &numTransferName, &transferName);
    AIM_STATUS(aimInfo, status);

    numDataTransferElement = 0;
    elementIndex = 0;
    elementCount = 0;

    for (transferNameIndex = 0; transferNameIndex < numTransferName; transferNameIndex++) {

        status = aim_getDiscr(aimInfo, transferName[transferNameIndex], &dataTransferDiscreteObj);
        if (status == CAPS_NOTFOUND) continue;
        AIM_STATUS(aimInfo, status);

        status = aim_getDataSet(dataTransferDiscreteObj,
                                "Pressure",
                                &dataTransferMethod,
                                &numDataTransferPoint,
                                &dataTransferRank,
                                &dataTransferData,
                                &units);
        if (status == CAPS_NOTFOUND) continue; // If no elements in this object skip to next transfer name
        AIM_STATUS(aimInfo, status);

        printf("Extracting external pressure loads from data transfer....\n");

        // If we do have data ready, how many elements there are?
        if (numDataTransferPoint == 1) {
            AIM_ERROR(aimInfo, "Pressures not initialized!");
            status = CAPS_BADINIT;
            goto cleanup;
        }

        if (dataTransferRank != 1) {
            AIM_ERROR(aimInfo, "Pressure transfer data found however rank is %d not 1!!!!", dataTransferRank);
            status = CAPS_BADRANK;
            goto cleanup;
        }

        discElements = 0;
        for (bIndex = 0; bIndex < dataTransferDiscreteObj->nBodys; bIndex++)
          discElements += dataTransferDiscreteObj->bodys[bIndex].nElems;

        numDataTransferElement += discElements;
        printf("\tTransferName = %s\n", transferName[transferNameIndex]);
        printf("\tNumber of Elements = %d (total = %d)\n", discElements, numDataTransferElement);

        // allocate arrays
        AIM_REALL(feaLoad->elementIDSet, numDataTransferElement, int, aimInfo, status);
        AIM_REALL(feaLoad->pressureMultiDistributeForce, 4*numDataTransferElement, double, aimInfo, status);

        //Now lets loop through our ctria3 mesh and get the node indexes for each element
        for (i = 0; i < feaMesh->numElement; i++) {

            if (feaMesh->element[i].elementType != Triangle) continue;

            elementID = feaMesh->element[i].elementID;
            //printf("Element Id = %d\n", elementID);

            // elementID is 1 bias
            nodeIndex[0] = feaMesh->element[i].connectivity[0];
            nodeIndex[1] = feaMesh->element[i].connectivity[1];
            nodeIndex[2] = feaMesh->element[i].connectivity[2];

            //printf("Node Index = %d %d %d\n", nodeIndex[0], nodeIndex[1], nodeIndex[2]);

            transferIndex[0] = -1;
            transferIndex[1] = -1;
            transferIndex[2] = -1;
            transferIndex[3] = -1;

            // Loop through the nodeMap of the data set getting nodeIDs trying
            //   to match the nodes in the element
            for (j = 0; j < numDataTransferPoint; j++) {

                bIndex       = dataTransferDiscreteObj->tessGlobal[2*j  ];
                globalNodeID = dataTransferDiscreteObj->tessGlobal[2*j+1] +
                               dataTransferDiscreteObj->bodys[bIndex-1].globalOffset;

                if (nodeIndex[0] == globalNodeID) transferIndex[0] = j;
                if (nodeIndex[1] == globalNodeID) transferIndex[1] = j;
                if (nodeIndex[2] == globalNodeID) transferIndex[2] = j;

                // If the nodes completely match the nodes on the element - break
                if (transferIndex[0] >= 0 && transferIndex[1] >= 0 && transferIndex[2] >= 0) {
                    break;
                }
            }

            // If all the nodeIndexes match the transferIndex the element is in the data set
            //  so transfer the pressure forces
            if (transferIndex[0] >= 0 && transferIndex[1] >= 0 && transferIndex[2] >= 0) {

                feaLoad->elementIDSet[elementIndex] = elementID;

                feaLoad->pressureMultiDistributeForce[4*elementIndex+0] = dataTransferData[transferIndex[0]];
                feaLoad->pressureMultiDistributeForce[4*elementIndex+1] = dataTransferData[transferIndex[1]];
                feaLoad->pressureMultiDistributeForce[4*elementIndex+2] = dataTransferData[transferIndex[2]];
                feaLoad->pressureMultiDistributeForce[4*elementIndex+3] = 0.0;

                elementIndex += 1;
                elementCount += 1;
                feaLoad->numElementID  += 1;

            }
        }

        //Now lets loop through our cquad4 mesh and get the node indexes for each element
        for (i = 0; i < feaMesh->numElement; i++) {

            if (feaMesh->element[i].elementType != Quadrilateral) continue;

            elementID = feaMesh->element[i].elementID;

            // elementID is 1 bias
            nodeIndex[0] = feaMesh->element[i].connectivity[0];
            nodeIndex[1] = feaMesh->element[i].connectivity[1];
            nodeIndex[2] = feaMesh->element[i].connectivity[2];
            nodeIndex[3] = feaMesh->element[i].connectivity[3];

            transferIndex[0] = -1;
            transferIndex[1] = -1;
            transferIndex[2] = -1;
            transferIndex[3] = -1;

            // Loop through the nodeMap of the data set getting nodeIDs trying
            //   to match the nodes in the element
            for (j = 0; j < numDataTransferPoint; j++) {

                bIndex       = dataTransferDiscreteObj->tessGlobal[2*j  ];
                globalNodeID = dataTransferDiscreteObj->tessGlobal[2*j+1] +
                               dataTransferDiscreteObj->bodys[bIndex-1].globalOffset;

                if (nodeIndex[0] == globalNodeID) transferIndex[0] = j;
                if (nodeIndex[1] == globalNodeID) transferIndex[1] = j;
                if (nodeIndex[2] == globalNodeID) transferIndex[2] = j;
                if (nodeIndex[3] == globalNodeID) transferIndex[3] = j;

                // If the nodes completely match the nodes on the element - break
                if (transferIndex[0] >= 0 && transferIndex[1] >= 0 &&
                    transferIndex[2] >= 0 && transferIndex[3] >= 0) {
                    break;
                }
            }

            // If all the nodeIndexes match the transferIndex the element is in the data set
            //  so transfer the pressure forces
            if (transferIndex[0] >= 0 && transferIndex[1] >= 0 &&
                transferIndex[2] >= 0 && transferIndex[3] >= 0) {

                feaLoad->elementIDSet[elementIndex] = elementID;

                feaLoad->pressureMultiDistributeForce[4*elementIndex+0] = dataTransferData[transferIndex[0]];
                feaLoad->pressureMultiDistributeForce[4*elementIndex+1] = dataTransferData[transferIndex[1]];
                feaLoad->pressureMultiDistributeForce[4*elementIndex+2] = dataTransferData[transferIndex[2]];
                feaLoad->pressureMultiDistributeForce[4*elementIndex+3] = dataTransferData[transferIndex[3]];

                elementIndex += 1;
                elementCount += 1;
                feaLoad->numElementID  += 1;
            }
        }

        if (elementCount != numDataTransferElement) {
            AIM_ERROR(aimInfo, "Element transfer mismatch: number of elements found = %d, number"
                               " of elements in transfer data set %d", elementCount, numDataTransferElement);
            AIM_FREE(transferName);

            status = CAPS_MISMATCH;
            goto cleanup;
        }

        // Resize
        if (feaLoad->numElementID != numDataTransferElement) {
            AIM_REALL(feaLoad->pressureMultiDistributeForce, 4*feaLoad->numElementID, double, aimInfo, status);
        }

    } // End data transfer name loop


    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(transferName);

    return status;
}

// Transfer external nodal pressure from the discrObj into the feaLoad structure
int fea_transferExternalPressureNode(void *aimInfo, feaLoadStruct *feaLoad) {

    // [in/out] feaLoad
    // [in] aimInfo

    int status; // Function status return

    int i, j, bIndex; // Indexing

    // Variables used in global node mapping
    int globalNodeID;

    // Discrete data transfer variables
    capsDiscr *dataTransferDiscreteObj;
    char **transferName = NULL;
    int numTransferName, transferNameIndex, *navg=NULL;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int dataTransferRank;
    double *dataTransferData;
    char   *units;

    feaLoad->numGridID = 0;
    AIM_FREE(feaLoad->pressureMultiDistributeForce);
    AIM_FREE(feaLoad->gridIDSet);

    //See if we have data transfer information
    status = aim_getBounds(aimInfo, &numTransferName, &transferName);
    AIM_STATUS(aimInfo, status);

    for (transferNameIndex = 0; transferNameIndex < numTransferName; transferNameIndex++) {

        status = aim_getDiscr(aimInfo, transferName[transferNameIndex], &dataTransferDiscreteObj);
        if (status == CAPS_NOTFOUND) continue;
        AIM_STATUS(aimInfo, status);

        status = aim_getDataSet(dataTransferDiscreteObj,
                                "Pressure",
                                &dataTransferMethod,
                                &numDataTransferPoint,
                                &dataTransferRank,
                                &dataTransferData,
                                &units);
        if (status == CAPS_NOTFOUND) continue; // If no elements in this object skip to next transfer name
        AIM_STATUS(aimInfo, status);

        printf("Extracting external pressure loads from data transfer....\n");

        // If we do have data ready, how many elements there are?
        if (numDataTransferPoint == 1) {
            AIM_ERROR(aimInfo, "Pressure not initialized!");
            status = CAPS_BADINIT;
            goto cleanup;
        }

        if (dataTransferRank != 1) {
            AIM_ERROR(aimInfo, "Pressure transfer data found however rank is %d not 1!!!!", dataTransferRank);
            status = CAPS_BADRANK;
            goto cleanup;
        }

        printf("\tTransferName = %s\n", transferName[transferNameIndex]);
        printf("\tNumber of Nodes = %d (total = %d)\n", numDataTransferPoint, feaLoad->numGridID+numDataTransferPoint);

        // allocate arrays
        AIM_REALL(feaLoad->gridIDSet, feaLoad->numGridID+numDataTransferPoint, int, aimInfo, status);
        AIM_REALL(feaLoad->pressureMultiDistributeForce, feaLoad->numGridID+numDataTransferPoint, double, aimInfo, status);
        AIM_REALL(navg, feaLoad->numGridID+numDataTransferPoint, int, aimInfo, status);

        // Loop through the nodeMap of the data set getting nodeIDs trying
        //   to match the nodes in the element
        for (i = 0; i < numDataTransferPoint; i++) {

            bIndex       = dataTransferDiscreteObj->tessGlobal[2*i  ];
            globalNodeID = dataTransferDiscreteObj->tessGlobal[2*i+1] +
                           dataTransferDiscreteObj->bodys[bIndex-1].globalOffset;

            for (j = 0; j < feaLoad->numGridID; j++) {
              if (feaLoad->gridIDSet[j] == globalNodeID) {
                feaLoad->pressureMultiDistributeForce[j] += dataTransferData[i];
                navg[j]++;
                break;
              }
            }

            if (j < feaLoad->numGridID) continue;

            feaLoad->gridIDSet[feaLoad->numGridID] = globalNodeID;
            feaLoad->pressureMultiDistributeForce[feaLoad->numGridID] = dataTransferData[i];
            navg[feaLoad->numGridID] = 1;
            feaLoad->numGridID += 1;
        }
    } // End data transfer name loop

    for (j = 0; j < feaLoad->numGridID; j++) {
        feaLoad->pressureMultiDistributeForce[j] /= navg[j];
    }

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(transferName);
    AIM_FREE(navg);

    return status;
}

// Transfer external temperature from the discrObj into the feaLoad structure
int fea_transferExternalTemperature(void *aimInfo, feaLoadStruct *feaLoad) {

    // [in/out] feaLoad
    // [in] aimInfo

    int status; // Function status return

    int i, j, bIndex; // Indexing

    // Variables used in global node mapping
    int globalNodeID;

    // Discrete data transfer variables
    capsDiscr *dataTransferDiscreteObj;
    char **transferName = NULL;
    int numTransferName, transferNameIndex, *navg=NULL;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int dataTransferRank;
    double *dataTransferData;
    char   *units;

    feaLoad->numGridID = 0;
    AIM_FREE(feaLoad->temperatureMultiDistribute);
    AIM_FREE(feaLoad->gridIDSet);

    //See if we have data transfer information
    status = aim_getBounds(aimInfo, &numTransferName, &transferName);
    AIM_STATUS(aimInfo, status);

    for (transferNameIndex = 0; transferNameIndex < numTransferName; transferNameIndex++) {

        status = aim_getDiscr(aimInfo, transferName[transferNameIndex], &dataTransferDiscreteObj);
        if (status == CAPS_NOTFOUND) continue;
        AIM_STATUS(aimInfo, status);

        status = aim_getDataSet(dataTransferDiscreteObj,
                                "Temperature",
                                &dataTransferMethod,
                                &numDataTransferPoint,
                                &dataTransferRank,
                                &dataTransferData,
                                &units);
        if (status == CAPS_NOTFOUND) continue; // If no elements in this object skip to next transfer name
        AIM_STATUS(aimInfo, status);

        printf("Extracting external thermal loads from data transfer....\n");

        // If we do have data ready, how many elements there are?
        if (numDataTransferPoint == 1) {
            AIM_ERROR(aimInfo, "Temperature not initialized!");
            status = CAPS_BADINIT;
            goto cleanup;
        }

        if (dataTransferRank != 1) {
            AIM_ERROR(aimInfo, "Temperature transfer data found however rank is %d not 1!!!!", dataTransferRank);
            status = CAPS_BADRANK;
            goto cleanup;
        }

        printf("\tTransferName = %s\n", transferName[transferNameIndex]);
        printf("\tNumber of Nodes = %d (total = %d)\n", numDataTransferPoint, feaLoad->numGridID+numDataTransferPoint);

        // allocate arrays
        AIM_REALL(feaLoad->gridIDSet, feaLoad->numGridID+numDataTransferPoint, int, aimInfo, status);
        AIM_REALL(feaLoad->temperatureMultiDistribute, feaLoad->numGridID+numDataTransferPoint, double, aimInfo, status);
        AIM_REALL(navg, feaLoad->numGridID+numDataTransferPoint, int, aimInfo, status);

        // Loop through the nodeMap of the data set getting nodeIDs trying
        //   to match the nodes in the element
        for (i = 0; i < numDataTransferPoint; i++) {

            bIndex       = dataTransferDiscreteObj->tessGlobal[2*i  ];
            globalNodeID = dataTransferDiscreteObj->tessGlobal[2*i+1] +
                           dataTransferDiscreteObj->bodys[bIndex-1].globalOffset;

            for (j = 0; j < feaLoad->numGridID; j++) {
              if (feaLoad->gridIDSet[j] == globalNodeID) {
                feaLoad->temperatureMultiDistribute[j] += dataTransferData[i];
                navg[j]++;
                break;
              }
            }

            if (j < feaLoad->numGridID) continue;

            feaLoad->gridIDSet[feaLoad->numGridID] = globalNodeID;
            feaLoad->temperatureMultiDistribute[feaLoad->numGridID] = dataTransferData[i];
            navg[feaLoad->numGridID] = 1;
            feaLoad->numGridID += 1;
        }
    } // End data transfer name loop

    for (j = 0; j < feaLoad->numGridID; j++) {
        feaLoad->temperatureMultiDistribute[j] /= navg[j];
    }

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(transferName);
    AIM_FREE(navg);

    return status;
}

// Retrieve aerodynamic reference quantities from bodies
int fea_retrieveAeroRef(int numBody, ego *bodies, feaAeroRefStruct *feaAeroRef) {

    int status; // Function return status
    int body;

    // EGADS return values
    int          atype, alen;
    const int    *ints;
    const char   *string;
    const double *reals;

    // Get reference quantities from the bodies
    for (body = 0; body < numBody; body++) {

        status = EG_attributeRet(bodies[body], "capsReferenceArea", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && atype == ATTRREAL) {

            feaAeroRef->refArea = (double) reals[0];
        }

        status = EG_attributeRet(bodies[body], "capsReferenceChord", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && atype == ATTRREAL){

            feaAeroRef->refChord = (double) reals[0];
        }

        status = EG_attributeRet(bodies[body], "capsReferenceSpan", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && atype == ATTRREAL) {

            feaAeroRef->refSpan = (double) reals[0];
        }

        /*
        status = EG_attributeRet(bodies[body], "capsReferenceX", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && atype == Double) {
            XX = (double) reals[0];
        }

        status = EG_attributeRet(bodies[body], "capsReferenceY", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && atype == Double) {
            XX = (double) reals[0];
        }

        status = EG_attributeRet(bodies[body], "capsReferenceZ", &atype, &alen, &ints, &reals, &string);
        if (status == EGADS_SUCCESS && atype == Double) {
            XX = (double) reals[0];
        }
         */
    }

    return CAPS_SUCCESS;
}


// Assign element "subtypes" based on properties set
int fea_assignElementSubType(int numProperty, feaPropertyStruct *feaProperty, meshStruct *feaMesh)
{
    int propertyIndex, i ;

    feaMeshDataStruct *feaData;

    if (numProperty > 0 && feaProperty == NULL) return CAPS_NULLVALUE;
    if (feaMesh == NULL) return CAPS_NULLVALUE;

    printf("Updating mesh element types based on properties input\n");

    for (propertyIndex = 0; propertyIndex < numProperty; propertyIndex++ ) {

        // Types that don't need subtypes - setting this correctly requires knowledge of what mesh_writeNastran
        //  writes by default for each mesh element type
        if (feaProperty[propertyIndex].propertyType == Rod   ||
            //feaProperty[propertyIndex].propertyType == Shell ||
            //feaProperty[propertyIndex].propertyType == Composite ||
            feaProperty[propertyIndex].propertyType == Solid) continue;

        for (i = 0; i < feaMesh->numElement; i++) {

            // What if this is a volume mesh we inherited ?
            if (feaMesh->element[i].analysisType != MeshStructure) {
                printf("Developer error: Analysis type not set to MeshStructure for element %d\n", feaMesh->element[i].elementID);
                return CAPS_BADVALUE;
            }

            // Matching proptery id
            if (feaMesh->element[i].analysisData != NULL) {
                feaData = (feaMeshDataStruct *) feaMesh->element[i].analysisData;

                if (feaData->propertyID != feaProperty[propertyIndex].propertyID) continue;

            } else {
                printf("Developer error: No 'feaData' set on element %d\n", feaMesh->element[i].elementID);
                return CAPS_BADVALUE;

            }

//            if (feaMesh->element[i].markerID != feaProperty[propertyIndex].propertyID) continue;
//
//            feaData = (feaMeshDataStruct *) feaMesh->element[i].analysisData;
//            if (feaData->propertyID != feaProperty[propertyIndex].propertyID) {
//                printf("Developer error: Property ID mismatch between element \"markerID\" (%d) and feaData \"propertyID\" (%d) for element %d\n",
//                        feaData->propertyID, feaProperty[propertyIndex].propertyID, feaMesh->element[i].elementID);
//                return CAPS_BADVALUE;
//            }


            if (feaProperty[propertyIndex].propertyType == ConcentratedMass &&
                feaMesh->element[i].elementType == Node) {

                feaData->elementSubType = ConcentratedMassElement;
            }

            if (feaProperty[propertyIndex].propertyType == Bar &&
                feaMesh->element[i].elementType == Line) {

                feaData->elementSubType = BarElement;
            }

            if (feaProperty[propertyIndex].propertyType == Beam &&
                feaMesh->element[i].elementType == Line) {

                feaData->elementSubType = BeamElement;
            }

            if (feaProperty[propertyIndex].propertyType == Shear &&
                feaMesh->element[i].elementType == Quadrilateral){

                feaData->elementSubType = ShearElement;
            }

            if (feaProperty[propertyIndex].propertyType == Membrane &&
                feaMesh->element[i].elementType == Quadrilateral){

                feaData->elementSubType = MembraneElement;
            }

            // Only need to set these if the zOffset is needed based on mesh_writeNASTRAN
            if ( ( feaProperty[propertyIndex].propertyType == Shell || feaProperty[propertyIndex].propertyType == Composite )  &&
                   feaProperty[propertyIndex].zOffsetRel != 0.0 &&
                 ( feaMesh->element[i].elementType == Quadrilateral ||
                  feaMesh->element[i].elementType == Triangle ||
                  feaMesh->element[i].elementType == Triangle_6 ||
                  feaMesh->element[i].elementType == Quadrilateral_8)) {

                feaData->elementSubType = ShellElement;
            }
        }
    }

    return CAPS_SUCCESS;
}

// Create connections for gluing - Connections are appended
int fea_glueMesh(void *aimInfo,
                 meshStruct *mesh,
                 int connectionID,
    /*@unused@*/ int connectionType,
                 int dofDependent,
                 char *slaveName,
                 int numMasterName,
                 char *masterName[],
                 mapAttrToIndexStruct *attrMap,
                 int maxNumMaster,
                 double searchRadius,
                 int *numConnect,
                 feaConnectionStruct *feaConnect[]) {

    // Input parameters
    //double exactTol=1E-9;
    //int weigthing = 1; // 1=linear

    int masterWeight = 1;
    int masterComponent = 123;

    int status;
    int i, masterIndex, slaveIndex, distIndex, attrIndex;

    int numMaster;

    int *glueConn=NULL;
    double *glueDist=NULL;

    double dist, maxDist;

    feaMeshDataStruct *feaData, *feaDataMaster;

//    int numConnect;
//    feaConnectionStruct *feaConnect;


// maps from nodeID to mesh->node index
//      mesh_nodeID2Array(nasMesh, &n2a);

    printf("\tCreating glue connections\n");

    if (searchRadius <=0) {
        printf("\tSearch radius must be greater than 0 when gluing, current value = %g\n", searchRadius);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    status = array_allocIntegerVector(maxNumMaster, -1, &glueConn);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = array_allocDoubleVector(maxNumMaster, 1E9, &glueDist);
    if (status != CAPS_SUCCESS) goto cleanup;

    for (slaveIndex = 0; slaveIndex < mesh->numNode; slaveIndex++) {

        feaData = (feaMeshDataStruct *) mesh->node[slaveIndex].analysisData;

        status = get_mapAttrToIndexIndex(attrMap, (const char *) slaveName, &attrIndex);
        if (status == CAPS_NOTFOUND) {
            printf("\tName %s not found in attribute map of capsConnect!!!!\n", slaveName);
            continue;
        } else if (status != CAPS_SUCCESS) return status;

        if (feaData->connectIndex != attrIndex) continue;

        // We have a slave at this point

       // printf("Slave - node ID %d\n", mesh->node[slaveIndex].nodeID);

        status = array_setDoubleVectorValue(maxNumMaster, 1E9, &glueDist);
        if (status != CAPS_SUCCESS) goto cleanup;

        status = array_setIntegerVectorValue(maxNumMaster, -1, &glueConn);
        if (status != CAPS_SUCCESS) goto cleanup;

        for (masterIndex = 0; masterIndex < mesh->numNode; masterIndex++) {

            if (slaveIndex == masterIndex) continue;

            feaDataMaster = (feaMeshDataStruct *) mesh->node[masterIndex].analysisData;

            for (i = 0; i < numMasterName; i++) {
                status = get_mapAttrToIndexIndex(attrMap, (const char *) masterName[i], &attrIndex);
                if (status == CAPS_NOTFOUND) {
                    printf("\tName %s not found in attribute map of capsConnect!!!!\n", masterName[i]);
                    continue;
                } else if (status != CAPS_SUCCESS) return status;

                if (feaDataMaster->connectIndex != attrIndex) continue;

                break;
            }

            if (i >= numMasterName) continue;

            // We have a potential master at this point

            dist = dist_DoubleVal(mesh->node[masterIndex].xyz, mesh->node[slaveIndex].xyz);

            //printf("dist = %f (radius = %f), masterIndex %d\n", dist, searchRadius, masterIndex);

            if ( dist > searchRadius) continue;
            AIM_NOTNULL(glueDist, aimInfo, status);

            status = array_maxDoubleValue(maxNumMaster, glueDist, &distIndex, &maxDist);
            if (status != CAPS_SUCCESS) goto cleanup;

            //printf("dist = %f, maxdist = %e (index = %d)\n", dist, maxDist, distIndex);

            if (dist < maxDist) {
                glueDist[distIndex] = dist;
                glueConn[distIndex] = masterIndex;
            }
        } // End of master loop

        // How many masters were found
        numMaster = 0;
        for (i = 0; i < maxNumMaster; i++) {
            if (glueDist[i] > searchRadius) continue; // Nothing was every set
            numMaster += 1;
        }

        if (numMaster <= 0) {
            printf("\tWarning: no masters were found for slave node (id = %d, slave name = %s)!\n", mesh->node[slaveIndex].nodeID, slaveName);
        } else {
            AIM_NOTNULL(glueConn, aimInfo, status);
            // Create and set connections
            status = fea_setConnection(slaveName,
                                       RigidBodyInterpolate,
                                       connectionID,
                                       mesh->numElement,
                                       dofDependent,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0,
                                       0,
                                       mesh->node[slaveIndex].nodeID,
                                       masterWeight,
                                       masterComponent,
                                       numMaster, glueConn,
                                       numConnect,feaConnect);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
//
//        *numConnect += 1;
//
//        (*feaConnect) = (feaConnectionStruct *) EG_reall((*feaConnect),*numConnect*sizeof(feaConnectionStruct));
//        if ((*feaConnect) == NULL) {
//            *numConnect = 0;
//            status = EGADS_MALLOC;
//            goto cleanup;
//        }
//
//        status = initiate_feaConnectionStruct(&(*feaConnect)[*numConnect-1]);
//        if (status != CAPS_SUCCESS) goto cleanup;
//
//        (*feaConnect)[*numConnect-1].name = (char *) EG_alloc((strlen(slaveName) + 1)*sizeof(char));
//        if ((*feaConnect)[*numConnect-1].name == NULL) {
//            status = EGADS_MALLOC;
//            goto cleanup;
//        }
//
//        memcpy((*feaConnect)[*numConnect-1].name, slaveName, strlen(slaveName)*sizeof(char));
//        (*feaConnect)[*numConnect-1].name[strlen(slaveName)] = '\0';
//
//        (*feaConnect)[*numConnect-1].connectionID = connectionID; // ConnectionTuple index
//        (*feaConnect)[*numConnect-1].connectionType = connectionType;
//
//        (*feaConnect)[*numConnect-1].elementID = *numConnect + mesh->numElement;
//
//        (*feaConnect)[*numConnect-1].dofDependent = dofDependent;
//
//        (*feaConnect)[*numConnect-1].connectivity[1] = mesh->node[slaveIndex].nodeID; // Dependent
//        (*feaConnect)[*numConnect-1].numMaster = numMaster; // Independent
//
//        (*feaConnect)[*numConnect-1].masterIDSet = (int *) EG_alloc(numMaster*sizeof(int)); // [numMaster]
//        if ((*feaConnect)[*numConnect-1].masterIDSet == NULL) {
//            status = EGADS_MALLOC;
//            goto cleanup;
//        }
//        (*feaConnect)[*numConnect-1].masterWeighting = (double *) EG_alloc(numMaster*sizeof(double));; // [numMaster]
//        if ((*feaConnect)[*numConnect-1].masterWeighting == NULL) {
//            status = EGADS_MALLOC;
//            goto cleanup;
//        }
//
//        (*feaConnect)[*numConnect-1].masterComponent =(int *) EG_alloc(numMaster*sizeof(int));; // [numMaster]
//        if ((*feaConnect)[*numConnect-1].masterComponent == NULL) {
//            status = EGADS_MALLOC;
//            goto cleanup;
//        }
//
//        printf("\tMasters (%d)", numMaster);
//
//        // Master values;
//        numMaster = 0;
//        for (i = 0; i < maxNumMaster; i++) {
//            if (glueDist[i] > searchRadius) continue; // Nothing was every set
//            printf(" %d(%d)", mesh->node[glueConn[i]].nodeID, glueConn[i]);
//            (*feaConnect)[*numConnect-1].masterIDSet[numMaster] = mesh->node[glueConn[i]].nodeID;
//            (*feaConnect)[*numConnect-1].masterWeighting[numMaster] = masterWeight;
//            (*feaConnect)[*numConnect-1].masterComponent[numMaster] = masterComponent;
//            numMaster += 1;
//        }
//        printf("\n");

    } // End of slave loop

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in fea_glueMesh, status = %d\n", status);

        if (glueDist != NULL) EG_free(glueDist);
        if (glueConn != NULL) EG_free(glueConn);

        return status;
}

//int fea_glueMesh2(meshStruct *meshIn,
//
//                 mapAttrToIndexStruct *attrMap) {
//
//    // Input parameters
//    int numGlue=5;
//    double exactTol=1E-9;
//    double radiusOfInfluence = 0.1;//1E-3;
//    int weigthing = 1; // 1=linear
//
//    int status;
//    int i, j, k, m, masterIndex, slaveIndex, distIndex;
//    int found;
//
//    int numMesh;
//    meshStruct *mesh, *master, *slave, *search;
//
//    int **glueConn=NULL;
//    double **glueDist=NULL;
//
//    double *tempDist=NULL;
//    int *tempGlue=NULL;
//    int numNode = 0;
//    int *meshIndexStart=NULL;
//
//    double dist, maxDist;
//
//    feaMeshDataStruct *feaData;
//
//// maps from nodeID to mesh->node index
////      mesh_nodeID2Array(nasMesh, &n2a);
//
//    printf("Welcome to glueing!\n");
//    numMesh = meshIn->numReferenceMesh;
//    mesh = meshIn->referenceMesh;
//
//    printf("Start with allocation!\n");
//    status = array_allocIntegerVector(numMesh, 0, &meshIndexStart);
//    if (status != CAPS_SUCCESS) goto cleanup;
//
//    // Determine number of nodes in reference meshes - Remember capsIgnore may have already been applied
//    for (i=0; i < numMesh; i++) {
//        meshIndexStart[i] = numNode;
//        numNode += mesh[i].numNode;
//    }
//
//    if (numNode != meshIn->numNode) {
//        printf("Inconsistent number of nodes determined!");
//        status = CAPS_MISMATCH;
//        goto cleanup;
//    }
//
//    // glueConn structure:
//    // glueConn[row][col] where [row] is node index and
//    // [row][0] : 1=master, 0=slave, -1=unused
//    // [row][1:numGlue+1] : nearest node indices, -1=unused
//    status = array_allocIntegerMatrix(numNode, numGlue+1, -1, &glueConn);
//    if (status != CAPS_SUCCESS) goto cleanup;
//
//    // glueDist structure:
//    // glueDist[row][col] where [row] is node index and
//    // [row][0:numGlue] : distance to nearest node indices found in glueConn[row][1:numGlue+1]
//    status = array_allocDoubleMatrix(numNode, numGlue, -1, &glueDist);
//    if (status != CAPS_SUCCESS) goto cleanup;
//
//    status = array_allocIntegerVector(numGlue, 0, &tempGlue);
//    if (status != CAPS_SUCCESS) goto cleanup;
//
//    status = array_allocDoubleVector(numGlue, 0, &tempDist);
//    if (status != CAPS_SUCCESS) goto cleanup;
//
//    printf("Done with allocation!\n");
//
//    // Lets mark all slaves
//    //
//
//    int attrIndex;
//    int numSlave=1;
//    char slaveName[] = "slave";
//    int numMasterName;
//    char **masterName;
//
//    for (i = 0; i < numNode; i++) {
//
//        feaData = (feaMeshDataStruct *) meshIn->node[i].analysisData;
//
//        for (j = 0; j < numSlave; j++) {
//            status = get_mapAttrToIndexIndex(attrMap, (const char *) slaveName, &attrIndex);
//            if (status == CAPS_NOTFOUND) {
//                printf("\tName %s not found in attribute map of capsConnect!!!!\n", slaveName);
//                continue;
//            } else if (status != CAPS_SUCCESS) return status;
//
//
//            if (feaData->connectIndex != attrIndex) continue;
//
//            glueConn[i][0] = 0;
//            break;
//        }
//
//        if (glueConn[i][0] != 0) continue;
//
//        status = string_toStringDynamicArray('["master"]', &numMasterName, &masterName);
//        if (status != CAPS_SUCCESS); goto cleanup;
//
//        for (k = 0; k < numNode; k++) {
//
//            if (i == k) continue;
//            for (m = 0; m < numMasterName; m++) {
//
//            }
//        }
//
//
//        (void) string_freeArray(numMasterName, &masterName);
//        masterName = NULL;
//        numMasterName = 0;
//    }
//
//
//    //    for (i = 0; i < numMesh; i++) { // THIS IS MOSTLY correct - error in master slave relation of RBE3 elements
////        master = &mesh[i];
////
//////        printf("Master Mesh %d (of %d)\n", i, numMesh);
////        //Exhaustive search
////        for (masterIndex = 0; masterIndex < master->numNode; masterIndex++) {
////
//////            printf("master index %d (of %d)\n", masterIndex, master->numNode);
////            // Slaves may be promoted to masters, but cannot be slaves to multiple masters
////            //if (glueConn[masterIndex+meshIndexStart[masterIndex]][0] != 0) continue;
////
////            status = array_setDoubleVectorValue(numGlue, 1E9, &tempDist);
////            if (status != CAPS_SUCCESS) goto cleanup;
////
////            status = array_setIntegerVectorValue(numGlue, -1, &tempGlue);
////            if (status != CAPS_SUCCESS) goto cleanup;
////
////            for (j = 0; j < numMesh; j++) {
////
////                if (i == j) continue;
////
////                slave = &mesh[j];
////
////                for (slaveIndex = 0; slaveIndex < slave->numNode; slaveIndex++) {
////
////                    // Slaves cannot have multiple masters - FALSE
////                    if (glueConn[slaveIndex+meshIndexStart[j]][0] == 0) continue;
////
////                    // Make sure master slave relation doesn't already exist
////                    found = (int) false;
////                    for (k =1; k < numGlue+1; k++) {
////                        if (masterIndex+meshIndexStart[i] == glueConn[slaveIndex+meshIndexStart[j]][k]){
////                            found = (int) true;
////                            break;
////                        }
////                    }
////
////                    if (found == (int) true) continue;
////
//////                    printf("glueConn %d\n", glueConn[slaveIndex+meshIndexStart[j]][0]);
////
////                    dist = dist_DoubleVal(master->node[masterIndex].xyz, slave->node[slaveIndex].xyz);
////
//////                    printf("dist = %f (radius = %f)\n", dist, radiusOfInfluence);
////
////                    if ( dist > radiusOfInfluence) continue;
////
////                    status = array_maxDoubleValue(numGlue, tempDist, &distIndex, &maxDist);
////                    if (status != CAPS_SUCCESS) goto cleanup;
////
////                    //printf("dist = %f, maxdist = %e (index = %d)\n", dist, maxDist, distIndex);
////
////                    if (dist < maxDist) {
////                        tempDist[distIndex] = dist;
////                        tempGlue[distIndex] = slaveIndex + meshIndexStart[j] ;
////                    }
////                }
////            }
////
////            if (tempDist[0] > radiusOfInfluence) continue; // Nothing was every set
////
////            // glueConn structure:
////            // glueConn[row][col] where [row] is node index and
////            // [row][0] : 1=master, 0=slave, -1=unused
////            // [row][1:numGlue+1] : nearest node indices, -1=unused
////
////            // glueDist[row][col] where [row] is node index and
////            // [row][0:numGlue] : distance to nearest node indices found in glueConn[row][1:numGlue+1]
////
////            glueConn[masterIndex+meshIndexStart[i]][0] = 1; //Mark as Master
////            distIndex = 0;
////            for (j = 0; j < numGlue; j++) {
////                if (tempDist[j] > radiusOfInfluence) continue; // Nothing was every set
////
////                glueDist[masterIndex+meshIndexStart[i]][distIndex] = tempDist[j]; // Distances for a master
////                glueConn[masterIndex+meshIndexStart[i]][distIndex+1] = tempGlue[j]; // Glue connections/slaves for a master
////
////                // Mark the masters slave's as slaves - A slave can become a master to another slave,
////                // but it can not be a slave to another master
////                if (glueConn[tempGlue[j]][0] < 0) glueConn[tempGlue[j]][0] = 0;
////
////                distIndex +=1;
////            }
////
////            printf("Set %d : %d: ",masterIndex+meshIndexStart[i], glueConn[masterIndex+meshIndexStart[i]][0]);
////            for(j = 0; j < numGlue; j++) {
////                printf("%d ", glueConn[masterIndex+meshIndexStart[i]][j+1]);
////            }
////            printf("\n");
////        }
////    }
////
////    for (i = 0; i < numNode;i++) {
////        if (glueConn[i][0] < 0) continue;
////        printf("Node %d : %d : ", i, glueConn[i][0]);
////        for(j = 0; j < numGlue; j++) {
//////            printf("%d (%.2f)", glueConn[i][j+1], glueDist[i][j]);
////            printf("%d ", glueConn[i][j+1]);
////        }
////        printf("\n");
////    }
////    // Finished getting all my glue
//
//    cleanup:
//        if (status != CAPS_SUCCESS) printf("\tPremature exit in fea_glueMesh, status = %d\n", status);
//
//        if (tempDist != NULL) EG_free(tempDist);
//        if (tempGlue != NULL) EG_free(tempGlue);
//        if (meshIndexStart != NULL) EG_free(meshIndexStart);
//
//        if (glueConn != NULL) (void) array_freeIntegerMatrix(numNode, numGlue, &glueConn);
//
//        if (glueDist != NULL) (void) array_freeDoubleMatrix(numNode, numGlue, &glueDist);
//
//        (void) string_freeArray(numMasterName, &masterName);
//        masterName = NULL;
//
//        return status;
//}

// Create a default analysis structure based on previous inputs
int fea_createDefaultAnalysis(void *aimInfo, feaProblemStruct *feaProblem, const char *analysisType) {

    int status;
    int i;

    capsTuple *tupleVal;

    char *json=NULL;
    int maxSize=2048;

    if (feaProblem == NULL) return CAPS_NULLVALUE;

    if (feaProblem->numAnalysis != 0 || feaProblem->feaAnalysis != NULL) {
        // Destroy our analysis structures coming in if aren't 0 and NULL already
        for (i = 0; i < feaProblem->numAnalysis; i++) {
            status = destroy_feaAnalysisStruct(&feaProblem->feaAnalysis[i]);
            if (status != CAPS_SUCCESS) return status;
        }

        if (feaProblem->feaAnalysis != NULL) EG_free(feaProblem->feaAnalysis);
        feaProblem->feaAnalysis = NULL;
        feaProblem->numAnalysis = 0;
    }

    tupleVal = (capsTuple *) EG_alloc(1*sizeof(capsTuple));
    if (tupleVal == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }
    tupleVal[0].name = NULL;
    tupleVal[0].value = NULL;

    json = (char *) EG_alloc(maxSize*sizeof(char));
    if (json == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    strcpy(json,"{\"analysisType\":\"");
    if (strcasecmp(analysisType, "Optimization") == 0) strcat(json,"Static");
    else strcat(json, analysisType);
    strcat(json,"\"");

    if (feaProblem->numLoad != 0) {
        strcat(json,",\"analysisLoad\":[");
        for (i = 0; i < feaProblem->numLoad; i++) {
            if (i != 0) strcat(json,",");
            strcat(json,"\"");
            strcat(json,feaProblem->feaLoad[i].name);
            strcat(json,"\"");
        }
        strcat(json, "]");
    }

    if (feaProblem->numConstraint != 0) {
        strcat(json,",\"analysisConstraint\":[");
        for (i = 0; i < feaProblem->numConstraint; i++) {
            if (i != 0) strcat(json,",");
            strcat(json,"\"");
            strcat(json,feaProblem->feaConstraint[i].name);
            strcat(json,"\"");
        }
        strcat(json, "]");
    }

    if (feaProblem->numSupport != 0) {
        strcat(json,",\"analysisSupport\":[");
        for (i = 0; i < feaProblem->numSupport; i++) {
            if (i != 0) strcat(json,",");
            strcat(json,"\"");
            strcat(json,feaProblem->feaSupport[i].name);
            strcat(json,"\"");
        }
        strcat(json, "]");
    }

    if (feaProblem->numDesignConstraint != 0) {
        strcat(json,",\"analysisDesignConstraint\":[");
        for (i = 0; i < feaProblem->numDesignConstraint; i++) {
            if (i != 0) strcat(json,",");
            strcat(json,"\"");
            strcat(json,feaProblem->feaDesignConstraint[i].name);
            strcat(json,"\"");
        }
        strcat(json, "]");
    }

    strcat(json,"}");

    tupleVal[0].name = EG_strdup("Default");
    tupleVal[0].value = json;

    //printf("Default analysis tuple - %s\n", tupleVal[0].value);

    status = fea_getAnalysis(aimInfo,
                             1,
                             tupleVal,
                             feaProblem);
    if (status != CAPS_SUCCESS) goto cleanup;

//    feaProblem->numAnalysis = 1;
//    feaProblem->feaAnalysis = (feaAnalysisStruct *) EG_alloc(feaProblem->numAnalysis*sizeof(feaAnalysisStruct));
//    if (feaProblem->feaAnalysis == NULL) return EGADS_MALLOC;
//
//    feaAnalysis=&feaProblem->feaAnalysis[0];
//    status = initiate_feaAnalysisStruct(feaAnalysis);
//    if (status != CAPS_SUCCESS) return status;
//
//
//    feaAnalysis->name = (char *) EG_alloc(((strlen(name)) + 1)*sizeof(char));
//    if (feaProblem->feaAnalysis->name == NULL) return EGADS_MALLOC;
//
//    memcpy(feaAnalysis->name, name, strlen(name)*sizeof(char));
//    feaAnalysis->name[strlen(name)] = '\0';
//
//    feaAnalysis->analysisID = 1;


    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in fea_createDefaultAnalysis, status = %d\n", status);

        if (tupleVal != NULL) {
            if (tupleVal[0].name != NULL) EG_free(tupleVal[0].name);
            if (tupleVal[0].value != NULL) EG_free(tupleVal[0].value);
            EG_free(tupleVal);
        }

        return status;
}

// Setup the default flutter velocities if not specified
int fea_defaultFlutterVelocity(feaAnalysisStruct *feaAnalysis) {
    int i;

    int defaultNum = 23;
    double velocity=0, dv=1, vmin=0, vmax=0;

    if (feaAnalysis->numFlutterVel == 0) {
        AIM_FREE(feaAnalysis->flutterVel);

        feaAnalysis->numFlutterVel = defaultNum;
        feaAnalysis->flutterVel = EG_alloc(feaAnalysis->numFlutterVel*sizeof(double));
        if (feaAnalysis->flutterVel == NULL) {
            return EGADS_MALLOC;
        }

        velocity = sqrt(2*feaAnalysis->dynamicPressure/feaAnalysis->density);
        vmin = velocity / 2.0;
        vmax = 2 * velocity;
        dv = (vmax - vmin) / (double) (feaAnalysis->numFlutterVel-3);

        for (i = 0; i < feaAnalysis->numFlutterVel-2; i++) {
            feaAnalysis->flutterVel[i+1] = vmin + (double) i * dv;
        }

        feaAnalysis->flutterVel[0] = velocity/10;
        feaAnalysis->flutterVel[feaAnalysis->numFlutterVel-1] = velocity*10;
    }

    return CAPS_SUCCESS;
}
