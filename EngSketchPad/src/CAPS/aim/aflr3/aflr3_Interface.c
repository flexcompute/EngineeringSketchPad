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

#include <AFLR3_LIB.h> // AFLR3_API Library include
#include <aflr4/AFLR4_LIB.h> // Bring in AFLR4 API library
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


int aflr3_to_MeshStruct( INT_ Number_of_Nodes,
                         INT_ Number_of_Surf_Trias,
                         INT_ Number_of_Surf_Quads,
                         INT_ Number_of_Vol_Tets,
                         INT_ Number_of_Vol_Pents_5,
                         INT_ Number_of_Vol_Pents_6,
                         INT_ Number_of_Vol_Hexs,
                         INT_1D *Surf_ID_Flag,
                         INT_3D *Surf_Tria_Connectivity,
                         INT_4D *Surf_Quad_Connectivity,
                         INT_4D *Vol_Tet_Connectivity,
                         INT_5D *Vol_Pent_5_Connectivity,
                         INT_6D *Vol_Pent_6_Connectivity,
                         INT_8D *Vol_Hex_Connectivity,
                         DOUBLE_3D *Coordinates,
                         meshStruct *genUnstrMesh) {

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
    genUnstrMesh->numNode = Number_of_Nodes;
    genUnstrMesh->numElement = Number_of_Surf_Trias  +
                               Number_of_Surf_Quads  +
                               Number_of_Vol_Tets    +
                               Number_of_Vol_Pents_5 +
                               Number_of_Vol_Pents_6 +
                               Number_of_Vol_Hexs;

    genUnstrMesh->meshQuickRef.useStartIndex = (int) true;

    genUnstrMesh->meshQuickRef.numTriangle      = Number_of_Surf_Trias;
    genUnstrMesh->meshQuickRef.numQuadrilateral = Number_of_Surf_Quads;

    genUnstrMesh->meshQuickRef.numTetrahedral = Number_of_Vol_Tets;
    genUnstrMesh->meshQuickRef.numPyramid     = Number_of_Vol_Pents_5;
    genUnstrMesh->meshQuickRef.numPrism       = Number_of_Vol_Pents_6;
    genUnstrMesh->meshQuickRef.numHexahedral  = Number_of_Vol_Hexs;

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

        genUnstrMesh->node[i].xyz[0] = Coordinates[i+1][0];
        genUnstrMesh->node[i].xyz[1] = Coordinates[i+1][1];
        genUnstrMesh->node[i].xyz[2] = Coordinates[i+1][2];
    }


    // Start of element index
    elementIndex = 0;

    // Elements-Set triangles
    if (Number_of_Surf_Trias > 0)
        genUnstrMesh->meshQuickRef.startIndexTriangle = elementIndex;

    numPoint = 0;
    for (i = 0; i < Number_of_Surf_Trias; i++) {

        genUnstrMesh->element[elementIndex].elementType = Triangle;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;

        genUnstrMesh->element[elementIndex].markerID = Surf_ID_Flag[i+1];

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                Surf_Tria_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set quadrilateral
    if (Number_of_Surf_Quads > 0)
        genUnstrMesh->meshQuickRef.startIndexQuadrilateral = elementIndex;

    for (i = 0; i < Number_of_Surf_Quads; i++) {

        genUnstrMesh->element[elementIndex].elementType = Quadrilateral;
        genUnstrMesh->element[elementIndex].elementID   = elementIndex+1;
        genUnstrMesh->element[elementIndex].markerID    =
            Surf_ID_Flag[Number_of_Surf_Trias+i+1];

        status = mesh_allocMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (i == 0) { // Only need this once
            numPoint = mesh_numMeshElementConnectivity(&genUnstrMesh->element[elementIndex]);
        }

        for (j = 0; j < numPoint; j++ ) {
            genUnstrMesh->element[elementIndex].connectivity[j] =
                Surf_Quad_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Tetrahedral
    if (Number_of_Vol_Tets > 0)
        genUnstrMesh->meshQuickRef.startIndexTetrahedral = elementIndex;

    for (i = 0; i < Number_of_Vol_Tets; i++) {

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
                Vol_Tet_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Pyramid
    if (Number_of_Vol_Pents_5 > 0)
        genUnstrMesh->meshQuickRef.startIndexPyramid = elementIndex;

    for (i = 0; i < Number_of_Vol_Pents_5; i++) {

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
                Vol_Pent_5_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Prism
    if (Number_of_Vol_Pents_6 > 0)
        genUnstrMesh->meshQuickRef.startIndexPrism = elementIndex;

    for (i = 0; i < Number_of_Vol_Pents_6; i++) {

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
                Vol_Pent_6_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    // Elements -Set Hexa
    if (Number_of_Vol_Hexs > 0)
        genUnstrMesh->meshQuickRef.startIndexHexahedral = elementIndex;

    for (i = 0; i < Number_of_Vol_Hexs; i++) {

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
                Vol_Hex_Connectivity[i+1][j];
        }

        elementIndex += 1;
    }

    status = CAPS_SUCCESS;


cleanup:
    if (status != CAPS_SUCCESS)
        printf("Premature exit in aflr3_to_MeshStruct status = %d\n", status);

    return status;
}


int aflr3_Volume_Mesh (void *aimInfo,
                       capsValue *aimInputs,
                       int ibodyOffset,
                       meshInputStruct meshInput,
                       const char *fileName,
                       int boundingBoxIndex,
                       int createBL,
                       double globalBLSpacing,
                       double globalBLThickness,
                       double capsMeshLength,
                       const mapAttrToIndexStruct *groupMap,
                       const mapAttrToIndexStruct *meshMap,
                       int numMeshProp,
                       meshSizingStruct *meshProp,
                       meshStruct *volumeMesh)
{
    int i, d; // Indexing

    int propIndex;
    int bodyIndex, state, np, ibody, nbody;

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

    int transp_intrnl = (int)false;
    INT_ Input_Surf_Trias = 0;
    int face_ntri, isurf, itri;
    const double *face_xyz, *face_uv;
    const int *face_ptype, *face_pindex, *face_tris, *face_tric;

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

    INT_ create_tess_mode = 2;
    INT_ set_node_map = 1;

    int aflr4_argc = 1;
    char **aflr4_argv = NULL;

    UG_Param_Struct *AFLR4_Param_Struct_Ptr = NULL;

    const char* bcType = NULL;
    INT_ nbl, nbldiff;
    INT_ index=0;

    // Declare AFLR3 grid generation variables.
    INT_1D *Edge_ID_Flag = NULL;
    INT_1D *Surf_Error_Flag= NULL;
    INT_1D *Surf_Grid_BC_Flag = NULL;
    INT_1D *Surf_ID_Flag = NULL;
    INT_1D *Surf_Reconnection_Flag = NULL;
    INT_2D *Surf_Edge_Connectivity = NULL;
    INT_3D *Surf_Tria_Connectivity = NULL;
    INT_4D *Surf_Quad_Connectivity = NULL;
    INT_1D *Vol_ID_Flag = NULL;
    INT_4D *Vol_Tet_Connectivity = NULL;
    INT_5D *Vol_Pent_5_Connectivity = NULL;
    INT_6D *Vol_Pent_6_Connectivity = NULL;
    INT_8D *Vol_Hex_Connectivity = NULL;

    DOUBLE_3D *Coordinates = NULL;
    DOUBLE_1D *BL_Normal_Spacing = NULL;
    DOUBLE_1D *BL_Thickness = NULL;

    INT_4D *BG_Vol_Tet_Neigbors = NULL;
    INT_4D *BG_Vol_Tet_Connectivity = NULL;

    DOUBLE_3D *BG_Coordinates = NULL;
    DOUBLE_1D *BG_Spacing = NULL;
    DOUBLE_6D *BG_Metric = NULL;

    DOUBLE_3D *Source_Coordinates = NULL;
    DOUBLE_1D *Source_Spacing = NULL;
    DOUBLE_6D *Source_Metric = NULL;

    INT_ status= 0;
    INT_ Message_Flag = 1;
    INT_ Number_of_BL_Vol_Tets = 0;
    INT_ Number_of_Nodes = 0;
    INT_ Number_of_Surf_Edges = 0;
    INT_ Number_of_Surf_Quads = 0;
    INT_ Number_of_Surf_Trias = 0;
    INT_ Number_of_Vol_Hexs = 0;
    INT_ Number_of_Vol_Pents_5 = 0;
    INT_ Number_of_Vol_Pents_6 = 0;
    INT_ Number_of_Vol_Tets = 0;

    INT_ Number_of_BG_Nodes = 0;
    INT_ Number_of_BG_Vol_Tets = 0;

    INT_ Number_of_Source_Nodes = 0;

    // Declare main program variables.

    INT_ idef, noquad = 0;
    INT_ mclosed = 1;
    INT_ glue_trnsp = 1, nfree;
    DOUBLE_2D *u = NULL;

    // INT_ M_BL_Thickness, M_Initial_Normal_Spacing;

    CHAR_UG_MAX Output_Case_Name;

    DOUBLE_1D *BG_U_Scalars = NULL;
    DOUBLE_6D *BG_U_Metrics = NULL;

    INT_ * bc_ids_vector = NULL;
    DOUBLE_1D *bl_ds_vector = NULL;
    DOUBLE_1D *bl_del_vector = NULL;

    void *ext_cad_data = NULL;
    egads_struct *ptr = NULL;

    // Set and register program parameter functions.

    ug_set_prog_param_code (3);

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
    aflr4_register_cad_geom_data_cleanup (egads_cad_geom_data_cleanup);
    aflr4_register_cad_geom_file_read (egads_cad_geom_file_read);
    aflr4_register_cad_geom_file_write (egads_cad_geom_file_write);
    aflr4_register_cad_geom_create_tess (egads_aflr4_create_tess);
    aflr4_register_cad_geom_reset_attr (egads_cad_geom_reset_attr);
    aflr4_register_cad_geom_setup (egads_cad_geom_setup);
    aflr4_register_cad_tess_to_dgeom (egads_aflr4_tess_to_dgeom);
    aflr4_register_set_ext_cad_data (egads_set_ext_cad_data);

    dgeom_register_cad_eval_curv_at_uv (egads_eval_curv_at_uv);
    dgeom_register_cad_eval_edge_arclen (egads_eval_edge_arclen);
    dgeom_register_cad_eval_uv_bounds (egads_eval_uv_bounds);
    dgeom_register_cad_eval_uv_at_t (egads_eval_uv_at_t);
    dgeom_register_cad_eval_uv_at_xyz (egads_eval_uv_at_xyz);
    dgeom_register_cad_eval_xyz_at_t (egads_eval_xyz_at_u);
    dgeom_register_cad_eval_xyz_at_uv (egads_eval_xyz_at_uv);
    dgeom_register_discrete_eval_xyz_at_t (surfgen_discrete_eval_xyz_at_t);


    status = ug_add_new_arg (&aflr3_argv, (char*)"allocate_and_initialize_argv");
    AIM_STATUS(aimInfo, status);

    // Set other command options
    if (createBL == (int) true) {
        status = ug_add_flag_arg ((char*)"mbl=1", &aflr3_argc, &aflr3_argv);
        AIM_STATUS(aimInfo, status);

        if (aimInputs[BL_Max_Layers-1].nullVal == NotNull) {
          nbl = aimInputs[BL_Max_Layers-1].vals.integer;
          status = ug_add_flag_arg ("nbl", &aflr3_argc, &aflr3_argv); AIM_STATUS(aimInfo, status);
          status = ug_add_int_arg (  nbl , &aflr3_argc, &aflr3_argv); AIM_STATUS(aimInfo, status);
        }

        if (aimInputs[BL_Max_Layer_Diff-1].nullVal == NotNull) {
          nbldiff = aimInputs[BL_Max_Layer_Diff-1].vals.integer;
          status = ug_add_flag_arg ("nbldiff", &aflr3_argc, &aflr3_argv); AIM_STATUS(aimInfo, status);
          status = ug_add_int_arg (  nbldiff , &aflr3_argc, &aflr3_argv); AIM_STATUS(aimInfo, status);
        }

        status = ug_add_flag_arg ((char*)"mblelc=1", &aflr3_argc, &aflr3_argv);
        AIM_STATUS(aimInfo, status);

    } else {
        status = ug_add_flag_arg ((char*)"mbl=0", &aflr3_argc, &aflr3_argv);
        AIM_STATUS(aimInfo, status);
    }

    status = ug_add_flag_arg ((char*)"mrecm=3" , &aflr3_argc, &aflr3_argv);
    AIM_STATUS(aimInfo, status);
    status = ug_add_flag_arg ((char*)"mrecqm=3", &aflr3_argc, &aflr3_argv);
    AIM_STATUS(aimInfo, status);

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

    // Set meshInputs
    if (meshInput.quiet == 1) Message_Flag = 0;
    else Message_Flag = 1;

    // find all bodies with all AFLR_GBC TRANSP SRC/INTRNL

    AIM_ALLOC(transpBody, volumeMesh->numReferenceMesh, int, aimInfo, status);
    for (i = 0; i < volumeMesh->numReferenceMesh; i++) transpBody[i] = 0;

    for (bodyIndex = 0; bodyIndex < volumeMesh->numReferenceMesh; bodyIndex++) {

      status = EG_statusTessBody(volumeMesh->referenceMesh[bodyIndex].egadsTess, &body, &state, &np);
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

        // HACK to release ESP 1.21
        if (bcType != NULL &&
            strncasecmp(bcType, "TRANSP_INTRNL_UG3_GBC", 20) == 0) {
          transp_intrnl = (int)true;
        }

        // check to see if all faces on a body are TRANSP SRC/INTRNL
        if (bcType != NULL &&
            (strncasecmp(bcType, "TRANSP_SRC_UG3_GBC", 18) == 0 ||
             strncasecmp(bcType, "TRANSP_INTRNL_UG3_GBC", 20) == 0)) {
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

    for (bodyIndex = 0, nbody = 0; bodyIndex < volumeMesh->numReferenceMesh; bodyIndex++) {
        if (transpBody[bodyIndex] != 1) nbody++;
    }

    AIM_ALLOC(copy_body_tess, 2*volumeMesh->numReferenceMesh, ego, aimInfo, status);

    // Copy bodies making sure TRANSP bodies are last
    ibody = 0;
    for (itransp = 1; itransp >= -1; itransp -= 2) {
      for (bodyIndex = 0; bodyIndex < volumeMesh->numReferenceMesh; bodyIndex++) {
        if (transpBody[bodyIndex] == itransp) continue;

        status = EG_statusTessBody(volumeMesh->referenceMesh[bodyIndex].egadsTess, &body, &state, &np);
        AIM_STATUS(aimInfo, status);

        status = EG_copyObject(body, NULL, &copy_body_tess[ibody]);
        AIM_STATUS(aimInfo, status);

        status = EG_copyObject(volumeMesh->referenceMesh[bodyIndex].egadsTess, copy_body_tess[ibody],
                               &copy_body_tess[volumeMesh->numReferenceMesh+ibody]);
        AIM_STATUS(aimInfo, status);

        status = EG_getBodyTopos(copy_body_tess[ibody], NULL, FACE, &numFace, &faces);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(faces, aimInfo, status);

        if (transpBody[bodyIndex] == -1) {
            for (iface = 0; iface < numFace; iface++) {
                status = EG_getTessFace(volumeMesh->referenceMesh[bodyIndex].egadsTess,
                                        iface+1, &nnode_face, &face_xyz, &face_uv,
                                        &face_ptype, &face_pindex, &face_ntri, &face_tris, &face_tric);
                AIM_STATUS(aimInfo, status);

                Input_Surf_Trias += face_ntri;
            }
        }

        // Build up an array of all faces in the model
        AIM_REALL(modelFaces, numModelFace+numFace, ego, aimInfo, status);
        for (i = 0; i < numFace; i++) modelFaces[numModelFace+i] = faces[i];
        AIM_REALL(faceBodyIndex, numModelFace+numFace, int, aimInfo, status);
        for (i = 0; i < numFace; i++) faceBodyIndex[numModelFace+i] = bodyIndex;

        AIM_REALL(faceGroupIndex, numModelFace+numFace, int, aimInfo, status);
        for (i = 0; i < numFace; i++) {
          status = retrieve_CAPSGroupAttr(faces[i], &groupName);
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
      bl_ds_vector[iface]  = globalBLSpacing;
      bl_del_vector[iface] = globalBLThickness;

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

    status = EG_makeTopology(context, NULL, MODEL, 2*volumeMesh->numReferenceMesh, NULL, volumeMesh->numReferenceMesh,
                             copy_body_tess, NULL, &model);
    AIM_STATUS(aimInfo, status);


    // allocate and initialize a new argument vector

    status = ug_add_new_arg (&aflr4_argv, "allocate_and_initialize_argv");
    AIM_STATUS(aimInfo, status);

    // setup input parameter structure

    status = aflr4_setup_param (Message_Flag, 0, aflr4_argc, aflr4_argv,
                                &AFLR4_Param_Struct_Ptr);
    AIM_STATUS(aimInfo, status);

    (void)ug_set_int_param ("geom_type", 1, AFLR4_Param_Struct_Ptr);
    (void)ug_set_int_param ("mmsg", Message_Flag, AFLR4_Param_Struct_Ptr);

    // set CAD geometry data structure

    status = aflr4_set_ext_cad_data (&model);
    AIM_STATUS(aimInfo, status);

    // setup geometry data

    status = aflr4_setup_and_grid_gen (0, AFLR4_Param_Struct_Ptr);
    AIM_STATUS(aimInfo, status);

    // set dgeom from input data

    status = aflr4_cad_tess_to_dgeom ();
    AIM_STATUS(aimInfo, status);

    dgeom_def_get_idef (index, &idef);


    // add and glue all individual surfaces within composite surface definition

    status = dgeom_add_and_glue_comp (glue_trnsp, idef, mclosed, Message_Flag, &nfree);
    AIM_STATUS(aimInfo, status);

    //------------------------------------------------ GET COPY OF SURFACE MESH DATA

    status = aflr4_get_def (idef, noquad,
                            &Number_of_Surf_Edges,
                            &Number_of_Surf_Trias,
                            &Number_of_Nodes,
                            &Number_of_Surf_Quads,
                            &Surf_Grid_BC_Flag,
                            &Edge_ID_Flag,
                            &Surf_ID_Flag,
                            &Surf_Edge_Connectivity,
                            &Surf_Tria_Connectivity,
                            &Surf_Quad_Connectivity,
                            &u, &Coordinates);
    AIM_STATUS(aimInfo, status);

    aflr4_free_all (0);

    // check that all the inputs

    status = ug_check_prog_param(aflr3_argv, aflr3_argc, Message_Flag);
    AIM_STATUS(aimInfo, status);

    if (createBL == (int)true) {

      status = ug_add_flag_arg ("BC_IDs", &aflr3_argc, &aflr3_argv);                             AIM_STATUS(aimInfo, status);
      status = ug_add_int_vector_arg (numModelFace, bc_ids_vector, &aflr3_argc, &aflr3_argv);    AIM_STATUS(aimInfo, status);
      status = ug_add_flag_arg ("BL_DS", &aflr3_argc, &aflr3_argv);                              AIM_STATUS(aimInfo, status);
      status = ug_add_double_vector_arg (numModelFace, bl_ds_vector, &aflr3_argc, &aflr3_argv);  AIM_STATUS(aimInfo, status);
      status = ug_add_flag_arg ("BL_DEL", &aflr3_argc, &aflr3_argv);                             AIM_STATUS(aimInfo, status);
      status = ug_add_double_vector_arg (numModelFace, bl_del_vector, &aflr3_argc, &aflr3_argv); AIM_STATUS(aimInfo, status);

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
    AIM_STATUS(aimInfo, status);

    status = ug_add_flag_arg ((char *) "mmsg", &aflr3_argc, &aflr3_argv); AIM_STATUS(aimInfo, status);
    status = ug_add_int_arg(Message_Flag, &aflr3_argc, &aflr3_argv);      AIM_STATUS(aimInfo, status);

    // note that if mpfrmt is not set to 0 and BL mesh generation is on then only
    // the surface mesh is returned

    if (transp_intrnl == (int)true) {
        // HACK for TRANSP_INTRNL

        // This is needed aflr3 will segfault without it
        Surf_Reconnection_Flag = (INT_1D*) ug_malloc (&status, (Number_of_Surf_Trias+Number_of_Surf_Quads+1) * sizeof (INT_1D));
        if (status != 0) { status = EGADS_MALLOC; goto cleanup; }
        for (i = 0; i < Number_of_Surf_Trias+Number_of_Surf_Quads; i++) {
            Surf_Reconnection_Flag[i+1] = 7; // Do not allow swapping on the surface triangles
        }

        status = aflr3_grid_generator (aflr3_argc, aflr3_argv,
                                       &Number_of_BL_Vol_Tets,
                                       &Number_of_BG_Nodes,
                                       &Number_of_BG_Vol_Tets,
                                       &Number_of_Nodes,
                                       &Number_of_Source_Nodes,
                                       &Number_of_Surf_Quads,
                                       &Number_of_Surf_Trias,
                                       &Number_of_Vol_Hexs,
                                       &Number_of_Vol_Pents_5,
                                       &Number_of_Vol_Pents_6,
                                       &Number_of_Vol_Tets,
                                       &Surf_Error_Flag,
                                       &Surf_Grid_BC_Flag,
                                       &Surf_ID_Flag,
                                       &Surf_Reconnection_Flag,
                                       &Surf_Quad_Connectivity,
                                       &Surf_Tria_Connectivity,
                                       &Vol_Hex_Connectivity,
                                       &Vol_Pent_5_Connectivity,
                                       &Vol_Pent_6_Connectivity,
                                       &Vol_Tet_Connectivity,
                                       &BG_Vol_Tet_Neigbors,
                                       &BG_Vol_Tet_Connectivity,
                                       &Coordinates,
                                       &BL_Normal_Spacing,
                                       &BL_Thickness,
                                       &BG_Coordinates,
                                       &BG_Spacing,
                                       &BG_Metric,
                                       &Source_Coordinates,
                                       &Source_Spacing,
                                       &Source_Metric);
    } else {
      status = aflr3_vol_gen (aflr3_argc, aflr3_argv, Message_Flag,
                              &Number_of_Surf_Edges,
                              &Number_of_Surf_Trias,
                              &Number_of_Surf_Quads,
                              &Number_of_BL_Vol_Tets,
                              &Number_of_Vol_Tets,
                              &Number_of_Vol_Pents_5,
                              &Number_of_Vol_Pents_6,
                              &Number_of_Vol_Hexs,
                              &Number_of_Nodes,
                              &Number_of_BG_Vol_Tets,
                              &Number_of_BG_Nodes,
                              &Number_of_Source_Nodes,
                              &Edge_ID_Flag,
                              &Surf_Edge_Connectivity,
                              &Surf_Grid_BC_Flag,
                              &Surf_ID_Flag,
                              &Surf_Error_Flag,
                              &Surf_Reconnection_Flag,
                              &Surf_Tria_Connectivity,
                              &Surf_Quad_Connectivity,
                              &Vol_ID_Flag,
                              &Vol_Tet_Connectivity,
                              &Vol_Pent_5_Connectivity,
                              &Vol_Pent_6_Connectivity,
                              &Vol_Hex_Connectivity,
                              &BG_Vol_Tet_Neigbors,
                              &BG_Vol_Tet_Connectivity,
                              &Coordinates,
                              &BL_Normal_Spacing,
                              &BL_Thickness,
                              &BG_Coordinates,
                              &BG_Spacing,
                              &BG_Metric,
                              &Source_Coordinates,
                              &Source_Spacing,
                              &Source_Metric);
    }

    if (status != 0) {
        FILE *fp;

        strcpy (Output_Case_Name, "debug");
        ug3_write_surf_grid_error_file (Output_Case_Name,
                                        status,
                                        Number_of_Nodes,
                                        Number_of_Surf_Trias,
                                        Surf_Error_Flag,
                                        Surf_Grid_BC_Flag,
                                        Surf_ID_Flag,
                                        Surf_Tria_Connectivity,
                                        Coordinates);

        strcpy(aimFile, "aflr3_surf_debug.tec");

        fp = fopen(aimFile, "w");
        if (fp == NULL) goto cleanup;
        fprintf(fp, "VARIABLES = X, Y, Z, BC, ID\n");

        fprintf(fp, "ZONE N=%d, E=%d, F=FEBLOCK, ET=Triangle\n",
                Number_of_Nodes, Number_of_Surf_Trias);
        fprintf(fp, ", VARLOCATION=([1,2,3]=NODAL,[4,5]=CELLCENTERED)\n");
        // write nodal coordinates
        if (Coordinates != NULL)
            for (d = 0; d < 3; d++) {
                for (i = 0; i < Number_of_Nodes; i++) {
                    if (i % 5 == 0) fprintf( fp, "\n");
                    fprintf(fp, "%22.15e ", Coordinates[i+1][d]);
                }
                fprintf(fp, "\n");
            }

        if (Surf_Grid_BC_Flag != NULL)
            for (i = 0; i < Number_of_Surf_Trias; i++) {
                if (i % 5 == 0) fprintf( fp, "\n");
                fprintf(fp, "%d ", Surf_Grid_BC_Flag[i+1]);
            }

        if (Surf_ID_Flag != NULL)
            for (i = 0; i < Number_of_Surf_Trias; i++) {
                if (i % 5 == 0) fprintf( fp, "\n");
                if ((Surf_Error_Flag != NULL) && (Surf_Error_Flag[i+1] < 0))
                    fprintf(fp, "-1 ");
                else
                    fprintf(fp, "%d ", Surf_ID_Flag[i+1]);
            }

        // cell connectivity
        if (Surf_Tria_Connectivity != NULL)
            for (i = 0; i < Number_of_Surf_Trias; i++)
                fprintf(fp, "%d %d %d\n", Surf_Tria_Connectivity[i+1][0],
                                          Surf_Tria_Connectivity[i+1][1],
                                          Surf_Tria_Connectivity[i+1][2]);
/*@-dependenttrans@*/
        fclose(fp);
/*@+dependenttrans@*/
        AIM_ERROR(aimInfo, "AFLR3 Grid generation error. The input surfaces mesh has been written to: %s", aimFile);
        goto cleanup;
    }


    if (transp_intrnl == (int)true &&
        Input_Surf_Trias != Number_of_Surf_Trias) {
        printf("\nInfo: Use of TRANSP_INTRNL_UG3_GBC when the surface mesh is modified precludes mesh sensitivities and data transfer.\n");
        printf("      Surface Mesh Number of Triangles: %d\n", Input_Surf_Trias);
        printf("      Volume  Mesh Number of Triangles: %d\n", Number_of_Surf_Trias);
        printf("\n");

    } else if (transp_intrnl == (int)true) {

      // TRANSP_INTRNL_UG3_GBC is used, so aflr4_cad_geom_create_tess does not work,
      // but the surface mesh was preserved

      AIM_NOTNULL(Surf_ID_Flag, aimInfo, status);
      AIM_NOTNULL(Surf_Tria_Connectivity, aimInfo, status);

      isurf = 1;
      itri = 1;

      for (bodyIndex = 0; bodyIndex < volumeMesh->numReferenceMesh; bodyIndex++) {
        if (transpBody[bodyIndex] == 1) continue;

        // set the file name to write the egads file
        snprintf(bodyNumber, 42, AFLR3TESSFILE, bodyIndex+ibodyOffset);
        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        status = EG_statusTessBody(volumeMesh->referenceMesh[bodyIndex].egadsTess, &body, &state, &np);
        AIM_STATUS(aimInfo, status);

        status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
        AIM_STATUS(aimInfo, status);

        for (iface = 0; iface < nface; iface++, isurf++) {

          status = EG_getTessFace(volumeMesh->referenceMesh[bodyIndex].egadsTess,
                                  iface+1, &nnode_face, &face_xyz, &face_uv,
                                  &face_ptype, &face_pindex, &face_ntri, &face_tris, &face_tric);
          AIM_STATUS(aimInfo, status);

          AIM_ALLOC(face_node_map, nnode_face, int, aimInfo, status);

          for (i = 0; i < face_ntri; i++, itri++) {
            if (Surf_ID_Flag[itri] != isurf) {
              AIM_ERROR(aimInfo, "Developer error Surf_ID_Flag missmatch!");
              status = CAPS_BADTYPE;
              goto cleanup;
            }

            for (d = 0; d < 3; d++)
              face_node_map[face_tris[3*i+d]-1] = Surf_Tria_Connectivity[itri][d];
          }

          // Add the unique indexing of the tessellation
          snprintf(attrname, 128, "face_node_map_%d",iface+1);
          status = EG_attributeAdd(volumeMesh->referenceMesh[bodyIndex].egadsTess, attrname, ATTRINT,
                                   nnode_face, face_node_map, NULL, NULL);
          AIM_FREE (face_node_map);
          AIM_STATUS(aimInfo, status);
        }

        remove(aimFile);
        status = EG_saveTess(volumeMesh->referenceMesh[bodyIndex].egadsTess, aimFile);
        AIM_STATUS(aimInfo, status);
      }

    } else {

      // extract surface mesh data from an existing CAD Model with previously
      // generated Tess Objects and save the data in the DGEOM data structure

      // re-cerate the model without any TRANSP SRC/INTRNL bodies if necessary
      if (nbody < volumeMesh->numReferenceMesh) {
          status = EG_deleteObject(model);
          AIM_STATUS(aimInfo, status);

          for (bodyIndex = 0, ibody = 0; bodyIndex < volumeMesh->numReferenceMesh; bodyIndex++) {
              if (transpBody[bodyIndex] == 1) continue;

              status = EG_statusTessBody(volumeMesh->referenceMesh[bodyIndex].egadsTess, &body, &state, &np);
              AIM_STATUS(aimInfo, status);

              status = EG_copyObject(body, NULL, &copy_body_tess[ibody]);
              AIM_STATUS(aimInfo, status);

              status = EG_copyObject(volumeMesh->referenceMesh[bodyIndex].egadsTess, copy_body_tess[ibody], &copy_body_tess[nbody+ibody]);
              AIM_STATUS(aimInfo, status);

              ibody++;
          }

          status = EG_makeTopology(context, NULL, MODEL, 2*nbody, NULL, nbody,
                                   copy_body_tess, NULL, &model);
          AIM_STATUS(aimInfo, status);
      }

      status = aflr4_set_ext_cad_data (&model);
      AIM_STATUS(aimInfo, status);

      status = aflr4_setup_and_grid_gen (0, AFLR4_Param_Struct_Ptr);
      AIM_STATUS(aimInfo, status);

      // set dgeom from input data

      status = aflr4_cad_tess_to_dgeom ();
      AIM_STATUS(aimInfo, status);

      // save the surface mesh from AFLR3

      status = aflr4_input_to_dgeom (Number_of_Surf_Edges, Number_of_Surf_Trias, Number_of_Surf_Quads,
                                     Edge_ID_Flag, Surf_Edge_Connectivity, Surf_ID_Flag, Surf_Tria_Connectivity, Surf_Quad_Connectivity, Coordinates);
      AIM_STATUS(aimInfo, status);

      // create tess object from input surface mesh

      status = aflr4_cad_geom_create_tess (Message_Flag, create_tess_mode, set_node_map);
      AIM_STATUS(aimInfo, status);

      ext_cad_data = dgeom_get_ext_cad_data ();
      ptr = (egads_struct *) ext_cad_data;

      iface = 1;

      for (bodyIndex = 0, ibody = 0; bodyIndex < volumeMesh->numReferenceMesh; bodyIndex++) {
          if (transpBody[bodyIndex] == 1) continue;

          // set the file name to write the egads file
          snprintf(bodyNumber, 42, AFLR3TESSFILE, bodyIndex+ibodyOffset);
          status = aim_file(aimInfo, bodyNumber, aimFile);
          AIM_STATUS(aimInfo, status);

          status = EG_getBodyTopos(ptr->bodies[ibody], NULL, FACE, &nface, NULL);
          AIM_STATUS(aimInfo, status);

          for (i = 0; i < nface; i++, iface++) {

              status = egads_face_node_map_get (iface, &nnode_face, &face_node_map);
              AIM_STATUS(aimInfo, status);

              // Add the unique indexing of the tessellation
              snprintf(attrname, 128, "face_node_map_%d",i+1);
              status = EG_attributeAdd(ptr->tess[ibody], attrname, ATTRINT,
                                       nnode_face, face_node_map+1, NULL, NULL); // face_node_map is index on [i+1]
              AIM_STATUS(aimInfo, status);

              ug_free (face_node_map);
              face_node_map = NULL;
          }

          remove(aimFile);
          status = EG_saveTess(ptr->tess[ibody], aimFile);
          AIM_STATUS(aimInfo, status);
          ibody++;
      }
    }

    // map the face index to the capsGroup index
    AIM_NOTNULL(Surf_ID_Flag, aimInfo, status);
    AIM_NOTNULL(faceGroupIndex, aimInfo, status);

    for (i = 0; i < Number_of_Surf_Trias + Number_of_Surf_Quads; i++) {
      Surf_ID_Flag[i+1] = faceGroupIndex[Surf_ID_Flag[i+1]-1];
    }

    // Write the mesh to disk

    snprintf(aimFile, PATH_MAX, "%s.lb8.ugrid", fileName);

    status = ug_io_write_grid_file(aimFile,
                                   Message_Flag,
                                   Number_of_BL_Vol_Tets,
                                   Number_of_Nodes,
                                   Number_of_Surf_Quads,
                                   Number_of_Surf_Trias,
                                   Number_of_Vol_Hexs,
                                   Number_of_Vol_Pents_5,
                                   Number_of_Vol_Pents_6,
                                   Number_of_Vol_Tets,
                                   Surf_Grid_BC_Flag,
                                   Surf_ID_Flag,
                                   Surf_Reconnection_Flag,
                                   Surf_Quad_Connectivity,
                                   Surf_Tria_Connectivity,
                                   Vol_Hex_Connectivity,
                                   Vol_ID_Flag,
                                   Vol_Pent_5_Connectivity,
                                   Vol_Pent_6_Connectivity,
                                   Vol_Tet_Connectivity,
                                   Coordinates,
                                   BL_Normal_Spacing,
                                   BL_Thickness);
    AIM_STATUS(aimInfo, status);

    status = aflr3_to_MeshStruct(Number_of_Nodes,
                                 Number_of_Surf_Trias,
                                 Number_of_Surf_Quads,
                                 Number_of_Vol_Tets,
                                 Number_of_Vol_Pents_5,
                                 Number_of_Vol_Pents_6,
                                 Number_of_Vol_Hexs,
                                 Surf_ID_Flag,
                                 Surf_Tria_Connectivity,
                                 Surf_Quad_Connectivity,
                                 Vol_Tet_Connectivity,
                                 Vol_Pent_5_Connectivity,
                                 Vol_Pent_6_Connectivity,
                                 Vol_Hex_Connectivity,
                                 Coordinates,
                                 volumeMesh);
    AIM_STATUS(aimInfo, status);

    // Remove the temporary grid created by AFLR
    remove(".tmp.b8.ugrid");

    status = CAPS_SUCCESS;
cleanup:

    ug_free (face_node_map);
    face_node_map = NULL;

    // Free program arguements
    ug_free_argv(aflr3_argv); aflr3_argv = NULL;
    ug_free_argv(aflr4_argv); aflr4_argv = NULL;
    ug_free_param (AFLR4_Param_Struct_Ptr); AFLR4_Param_Struct_Ptr = NULL;

    // Free grid generation and parameter array space in structures.
    // Note that aflr3_grid_generator frees all background data.
    ug_io_free_grid(Surf_Grid_BC_Flag,
                    Surf_ID_Flag,
                    Surf_Reconnection_Flag,
                    Surf_Quad_Connectivity,
                    Surf_Tria_Connectivity,
                    Vol_Hex_Connectivity,
                    Vol_ID_Flag,
                    Vol_Pent_5_Connectivity,
                    Vol_Pent_6_Connectivity,
                    Vol_Tet_Connectivity,
                    Coordinates,
                    BL_Normal_Spacing,
                    BL_Thickness);
    Surf_Grid_BC_Flag = NULL;
    Surf_ID_Flag = NULL;
    Surf_Reconnection_Flag = NULL;
    Surf_Quad_Connectivity = NULL;
    Surf_Tria_Connectivity = NULL;
    Vol_Hex_Connectivity = NULL;
    Vol_ID_Flag = NULL;
    Vol_Pent_5_Connectivity = NULL;
    Vol_Pent_6_Connectivity = NULL;
    Vol_Tet_Connectivity = NULL;
    Coordinates = NULL;
    BL_Normal_Spacing = NULL;
    BL_Thickness = NULL;

    ug_free(Surf_Error_Flag);
    Surf_Error_Flag= NULL;

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

    ug_free (Edge_ID_Flag);
    ug_free (Surf_ID_Flag);
    ug_free (Surf_Edge_Connectivity);
    ug_free (u);

    Edge_ID_Flag = NULL;
    Surf_ID_Flag = NULL;
    Surf_Edge_Connectivity = NULL;
    u = NULL;

    ug_io_free_node(Source_Coordinates, Source_Spacing, Source_Metric);

    Source_Coordinates = NULL;
    Source_Spacing = NULL;
    Source_Metric = NULL;

    ug_free (bc_ids_vector);
    ug_free (bl_ds_vector);
    ug_free (bl_del_vector);

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

    return status;
}
