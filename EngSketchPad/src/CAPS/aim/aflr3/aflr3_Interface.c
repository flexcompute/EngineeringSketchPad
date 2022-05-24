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

#include "AFLR3_LIB.h" // AFLR3_API Library include

// UG_IO Library (file I/O) include
// UG_CPP Library (c++ code) include

#include "ug_io/UG_IO_LIB_INC.h"
#include "ug_cpp/UG_CPP_LIB_INC.h"

// UG_GQ Library (grid quality) include
// This is optional and not required for implementation.

#include "ug_gq/UG_GQ_LIB_INC.h"

// Includes for version functions
// This is optional and not required for implementation.

#include "otb/OTB_LIB_INC.h"
#include "rec3/REC3_LIB_INC.h"
#include "ug3/UG3_LIB_INC.h"
#include "dftr3/DFTR3_LIB_INC.h"
#include "ice3/ICE3_LIB_INC.h"

#ifdef _ENABLE_BL_
#include "aflr2c/AFLR2_LIB.h"
#include "anbl3/ANBL3_LIB.h"

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
                       meshInputStruct meshInput,
                       const char *fileName,
                       int createBL,
                       int blFlag[],
                       double blSpacing[],
                       double blThickness[],
                       double capsMeshLength,
                       int numMeshProp,
                       meshSizingStruct *meshProp,
                       meshStruct *surfaceMesh,
                       meshStruct *volumeMesh)
{
    int i, d, triIndex = 1, quadIndex = 1; // Indexing

    int propIndex, pointIndex[4] = {-1,-1,-1,-1};

    // Command line variables
    int  prog_argc   = 1;    // Number of arguments
    char **prog_argv = NULL; // String arrays
    char *meshInputString = NULL;
    char *rest = NULL, *token = NULL;
    char aimFile[PATH_MAX];

    INT_ bcType;
    INT_ nbl, nbldiff;

    // Declare AFLR3 grid generation variables.
    INT_1D *Surf_Error_Flag= NULL;
    INT_1D *Surf_Grid_BC_Flag = NULL;
    INT_1D *Surf_ID_Flag = NULL;
    INT_1D *Surf_Reconnection_Flag = NULL;
    INT_3D *Surf_Tria_Connectivity = NULL;
    INT_4D *Surf_Quad_Connectivity = NULL;
    INT_1D *Vol_ID_Flag = NULL;
    INT_4D *Vol_Tet_Connectivity = NULL;
    INT_5D *Vol_Pent_5_Connectivity = NULL;
    INT_6D *Vol_Pent_6_Connectivity = NULL;
    INT_8D *Vol_Hex_Connectivity = NULL;

    DOUBLE_3D *Coordinates = NULL;
    DOUBLE_1D *Initial_Normal_Spacing = NULL;
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
    // This is optional and not required for implementation.

    //CHAR_UG_MAX Text;

    //INT_ mbl, mblchki, mmet, mtr;
    INT_ M_BL_Thickness, M_Initial_Normal_Spacing;
    //INT_ M_Vol_ID_Flag;
    //INT_ Number_of_Vol_Elems;
    //INT_ PSDATA_File_Flag, Set_Vol_ID_Flag, Tags_Data_File_Flag;
    //INT_ Write_BL_Only;

    //INT_ BG_Flag = 1;
    //INT_ Check_Input_Flag = 0;
    //INT_ Help_Flag = 0;
    //INT_ Help_UG_IO_Flag = 0;
    INT_ Interior_Flag = 1;
    //INT_ mpfrmt = 0;
    //INT_ Number_of_Calls = 1;
    //INT_ Output_File_Flag = 0;
    INT_ Program_Flag = 1;
    INT_ Reconnection_Flag = 1;

    // Declare grid data I/O variables.
    // This is optional and not required for implementation.

    //CHAR_UG_MAX Arg_File_Name;
    //CHAR_UG_MAX Input_Case_Name;
    //CHAR_UG_MAX Input_Grid_File_Name;
    CHAR_UG_MAX Output_Case_Name;
    //CHAR_UG_MAX Output_Grid_File_Name;
    //CHAR_UG_MAX Output_Grid_File_Type_Suffix;
    //CHAR_UG_MAX BG_Case_Name;
    //CHAR_UG_MAX BG_Grid_File_Name;
    //CHAR_UG_MAX BG_Func_File_Name;
    //CHAR_UG_MAX Source_Data_File_Name;

    //INT_1D *BG_Surf_Grid_BC_Flag = NULL;
    //INT_1D *BG_Surf_ID_Flag = NULL;
    //INT_1D *BG_Surf_Reconnection_Flag = NULL;
    //INT_1D *BG_Vol_ID_Flag = NULL;
    //INT_3D *BG_Surf_Tria_Connectivity = NULL;
    //INT_4D *BG_Surf_Quad_Connectivity = NULL;
    //INT_5D *BG_Vol_Pent_5_Connectivity = NULL;
    //INT_6D *BG_Vol_Pent_6_Connectivity = NULL;
    //INT_8D *BG_Vol_Hex_Connectivity = NULL;

    //DOUBLE_1D *BG_Initial_Normal_Spacing = NULL;
    //DOUBLE_1D *BG_BL_Thickness = NULL;

    //CHAR_21 *BG_U_Scalar_Labels = NULL;
    //CHAR_21 *BG_U_Vector_Labels = NULL;
    //CHAR_21 *BG_U_Matrix_Labels = NULL;
    //CHAR_21 *BG_U_Metric_Labels = NULL;
    //INT_1D *BG_U_Scalar_Flags = NULL;
    //INT_1D *BG_U_Vector_Flags = NULL;
    //INT_1D *BG_U_Matrix_Flags = NULL;
    //INT_1D *BG_U_Metric_Flags = NULL;
    DOUBLE_1D *BG_U_Scalars = NULL;
    //DOUBLE_3D *BG_U_Vectors = NULL;
    //DOUBLE_9D *BG_U_Matrixes = NULL;
    DOUBLE_6D *BG_U_Metrics = NULL;

    //INT_ Number_of_BG_BL_Vol_Tets = 0;
    //INT_ Number_of_BG_Surf_Quads = 0;
    //INT_ Number_of_BG_Surf_Trias = 0;
    //INT_ Number_of_BG_Vol_Hexs = 0;
    //INT_ Number_of_BG_Vol_Pents_5 = 0;
    //INT_ Number_of_BG_Vol_Pents_6 = 0;

    //INT_ Number_of_BG_U_Nodes = 0;
    //INT_ Number_of_BG_U_Scalars = 0;
    //INT_ Number_of_BG_U_Vectors = 0;
    //INT_ Number_of_BG_U_Matrixes = 0;
    //INT_ Number_of_BG_U_Metrics = 0;

    //INT_ File_Data_Type, Output_Grid_File_Format;

    // Declare grid quality measure variables.
    // This is optional and not required for implementation.

    //INT_ GQ_Max_Dist_Increments, GQ_Surf_Measure_Flag, GQ_Vol_Measure_Flag;

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
    aflr3_register_anbl3(anbl3_grid_generator, anbl3_initialize_param);
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

    // Register CGNS I/O library grid file read and write routines.
    // This is optional and not required for implementation.

#ifdef _ENABLE_CGNS_LIB_
    ug_io_cgns_register_read_grid (cgns_ug_io_read_grid);
    ug_io_cgns_register_write_grid (cgns_ug_io_write_grid);
#endif

    // Register MESH I/O library grid and function file read and write routines.
    // This is optional and not required for implementation.

#ifdef _ENABLE_MESH_LIB_
    ug_io_mesh_register_read_grid (mesh_ug_io_read_grid);
    ug_io_mesh_register_read_func (mesh_ug_io_read_func);
    ug_io_mesh_register_write_grid (mesh_ug_io_write_grid);
#endif

    status = ug_add_new_arg (&prog_argv, (char*)"allocate_and_initialize_argv");
    if (status != CAPS_SUCCESS) goto cleanup;


    // Set other command options
    if (createBL == (int) true) {
        status = ug_add_flag_arg ((char*)"mbl=1", &prog_argc, &prog_argv);
        if (status != 0) goto cleanup;

        if (aimInputs[BL_Max_Layers-1].nullVal == NotNull) {
          nbl = aimInputs[BL_Max_Layers-1].vals.integer;
          status = ug_add_flag_arg ("nbl", &prog_argc, &prog_argv); if (status != 0) goto cleanup;
          status = ug_add_int_arg (  nbl , &prog_argc, &prog_argv); if (status != 0) goto cleanup;
        }

        if (aimInputs[BL_Max_Layer_Diff-1].nullVal == NotNull) {
          nbldiff = aimInputs[BL_Max_Layer_Diff-1].vals.integer;
          status = ug_add_flag_arg ("nbldiff", &prog_argc, &prog_argv); if (status != 0) goto cleanup;
          status = ug_add_int_arg (  nbldiff , &prog_argc, &prog_argv); if (status != 0) goto cleanup;
        }

        status = ug_add_flag_arg ((char*)"mblelc=1", &prog_argc, &prog_argv);
        if (status != 0) goto cleanup;

    } else {
        status = ug_add_flag_arg ((char*)"mbl=0", &prog_argc, &prog_argv);
        if (status != 0) goto cleanup;
    }

    status = ug_add_flag_arg ((char*)"mrecm=3" , &prog_argc, &prog_argv);
    if (status != 0) goto cleanup;
    status = ug_add_flag_arg ((char*)"mrecqm=3", &prog_argc, &prog_argv);
    if (status != 0) goto cleanup;

    // Parse input string
    if (meshInput.aflr3Input.meshInputString != NULL) {

        rest = meshInputString = EG_strdup(meshInput.aflr3Input.meshInputString);
        while ((token = strtok_r(rest, " ", &rest))) {
            status = ug_add_flag_arg (token, &prog_argc, &prog_argv);
            if (status != 0) {
                printf("Error: Failed to parse input string: %s\n", token);
                if (meshInputString != NULL)
                    printf("Complete input string: %s\n", meshInputString);
                goto cleanup;
            }
        }
    }

    //printf("Number of args = %d\n",prog_argc);
    //for (i = 0; i <prog_argc ; i++) printf("Arg %d = %s\n", i, prog_argv[i]);

    // Set meshInputs
    if (meshInput.quiet == 1) Message_Flag = 0;
    else Message_Flag = 1;

    // check that all the inputs

    status = ug_check_prog_param(prog_argv, prog_argc, Message_Flag);
    if (status != 0) goto cleanup; // Add error code

    // Start memory and file status monitors and create TMP file directory.
    // This is optional and not required for implementation.

    if (status == 0) (void) ug_startup(prog_argc, prog_argv);

    // Set program parameter functions.
    // This is optional and not required for implementation.
    /*
    if (status== 0) ug_set_prog_param_functions (3, ug_initialize_aflr_param,
                                                    ug_gq_initialize_param, NULL, NULL, NULL,
                                                    aflr3_initialize_param,
                                                    aflr3_anbl3_initialize_param,
                                                    ice3_initialize_param, NULL, NULL);
     */

    // Get program parameters.
    // This is optional and not required for implementation.

    /*    if (status== 0) status= aflr3_get_prog_param (prog_argv, prog_argc,
                                                  Input_Grid_File_Name,
                                                  Output_Grid_File_Name,
                                                  &BG_Flag, &Check_Input_Flag,
                                                  &Help_Flag, &Help_UG_IO_Flag,
                                                  &Message_Flag, &Number_of_Calls,
                                                  &Output_File_Flag, &Program_Flag,
                                                  &PSDATA_File_Flag, &Set_Vol_ID_Flag,
                                                  &Tags_Data_File_Flag, &Write_BL_Only,
                                                  &mbl, &mblchki, &mmet, &mtr,
                                                  &GQ_Max_Dist_Increments,
                                                  &GQ_Surf_Measure_Flag,
                                                  &GQ_Vol_Measure_Flag);

    printf("MBL = %d\n", mbl);
    printf("mblchk = %d\n", mblchki);
    printf("mmet = %d\n", mmet);
    printf("mtr = %d\n", mtr);
    printf("Input_Grid_File_Name = %s\n", Input_Grid_File_Name);
    printf("Output_Grid_File_Name = %s\n", Output_Grid_File_Name);*/


    // Write input parameter documentation.
    // This is optional and not required for implementation.

    //    if (status== 0) status= ug_io_write_param_info (0, 1, 1, 0, 1, 1, Help_UG_IO_Flag);

    /*    if (Help_Flag || Help_UG_IO_Flag)
    {
        printf("HELP !!!!\n");

        ug_free_argv (prog_argv);
        ug_free_argv (argv_copy);


        ug_shutdown ();

        exit (status);
    }*/

    // Check input and output grid file names and set case name.
    // This is optional and not required for implementation.

    /*if (Program_Flag)
        File_Data_Type = UG_IO_SURFACE_GRID;
    else
        File_Data_Type = UG_IO_VOLUME_GRID;

    strcpy (Input_Case_Name, "");

    if (status== 0)
        status= ug_io_check_file_name (Input_Case_Name,
                                      Input_Grid_File_Name,
                                      2, 3, File_Data_Type, -1);

    if (status== 0 && strcmp (Output_Grid_File_Name, "_null_"))
    {
        strcpy (Output_Case_Name, "");

        status= ug_io_check_file_name (Output_Case_Name,
                                       Output_Grid_File_Name,
                                    1, 3, UG_IO_SURFACE_GRID, 1);

        if (status<= 0 && strcmp (Output_Case_Name, "") == 0)
        {
            strcpy (Output_Case_Name, Input_Case_Name);

            status= ug_io_check_file_name (Output_Case_Name,
                                      Output_Grid_File_Name,
                                      3, 3, UG_IO_SURFACE_GRID, 1);
        }
    }
    else if (status== 0)
        strcpy (Output_Case_Name, Input_Case_Name);

    // Set background grid and function file names.
    // This is optional and not required for implementation.

    strcpy (BG_Case_Name, Input_Case_Name);
    strcat (BG_Case_Name, ".back");
    strcpy (BG_Grid_File_Name, "");

    if (status== 0 && Program_Flag && BG_Flag)
    {
        status= ug_io_find_file (BG_Case_Name, BG_Grid_File_Name,
                                UG_IO_VOLUME_GRID);

        if (status== 0)
            status= ug_io_find_file (BG_Case_Name, BG_Func_File_Name,
                                    UG_IO_FUNCTION_DATA);

        if (status== -1)
        {
            status= 0;

            strcpy (BG_Grid_File_Name, "");
            strcpy (BG_Func_File_Name, "");
        }
    }

    // Set input source data file name.
    // This is optional and not required for implementation.

    strcpy (Source_Data_File_Name, "");

    if (status== 0 && Program_Flag)
    {
        status= ug_io_find_file (Input_Case_Name, Source_Data_File_Name,
                                     UG_IO_NODE_DATA);

        if (status== -1)
        {
            status= 0;

            strcpy (Source_Data_File_Name, "");
        }
    }

    // Check input parameters in parameter structure from argument list.
    // This is optional and not required for implementation.

    if (status== 0 && Check_Input_Flag)
    {
        status= ug_check_prog_param (prog_argv, prog_argc, 1);

        ug_free_argv (prog_argv);
    //    ug_free_argv (argv_copy);

        ug_shutdown ();

        if (Message_Flag >= 1)
            ug_message (" ");

        exit (status);
    }

    // Open output file and set flag.
    // This is optional and not required for implementation.

    if (status== 0) ug_open_output_file (Output_Case_Name, "aflr3", Output_File_Flag);

     */

    // Generate program output message
    // This is optional and not required for implementation.

    /*    if (status== 0)
    {
        ug_message (" ");
        ug_message ("AFLR3    : ---------------------------------------");
        ug_message ("AFLR3    : AFLR3");
        ug_message ("AFLR3    : ADVANCING-FRONT/LOCAL-RECONNECTION");
        ug_message ("AFLR3    : UNSTRUCTURED GRID GENERATOR");
        sprintf (Text, "AFLR3    : Version Number %s", Version_Number);
        ug_message (Text);
        sprintf (Text, "AFLR3    : Version Date   %s", Version_Date);
        ug_message (Text);
        ug_message ("AFLR3    : ---------------------------------------");
        ug_message (" ");
    }*/

    // Set input parameters in parameter structure from argument list.
    // This is optional and not required for implementation.

    //    if (status== 0) status= ug_check_prog_param (prog_argv, prog_argc, Message_Flag);

    // Read grid data.
    // This is optional and not required for implementation.

    /*
      M_BL_Thickness = (mbl) ? 1: 0;
        M_Initial_Normal_Spacing = (mbl) ? 1: 0;
        M_Vol_ID_Flag = (Program_Flag == 0 && Set_Vol_ID_Flag == 1) ? 1: 0;
     */

    /*        status= ug_io_read_grid_file (Input_Grid_File_Name,
                                   Message_Flag,
                                   M_BL_Thickness,
                                   M_Initial_Normal_Spacing, 1, 1, 1,
                                   M_Vol_ID_Flag, 0, 0, 0, 0,
                                   &Number_of_BL_Vol_Tets,
                                   &Number_of_Nodes,
                                   &Number_of_Surf_Quads,
                                   &Number_of_Surf_Trias,
                                   &Number_of_Vol_Hexs,
                                   &Number_of_Vol_Pents_5,
                                   &Number_of_Vol_Pents_6,
                                   &Number_of_Vol_Tets,
                                   &Surf_Grid_BC_Flag,
                                   &Surf_ID_Flag,
                                   &Surf_Reconnection_Flag,
                                   &Surf_Quad_Connectivity,
                                   &Surf_Tria_Connectivity,
                                   &Vol_Hex_Connectivity,
                                   &Vol_ID_Flag,
                                   &Vol_Pent_5_Connectivity,
                                   &Vol_Pent_6_Connectivity,
                                   &Vol_Tet_Connectivity,
                                   &Coordinates,
                                   &Initial_Normal_Spacing,
                                   &BL_Thickness);*/



    // Set flags for BL parameters allocation
    if (createBL == (int) true) M_Initial_Normal_Spacing = 1;
    else                        M_Initial_Normal_Spacing = 0;

    if (createBL == (int) true) M_BL_Thickness = 1;
    else                        M_BL_Thickness = 0;

    //printf("BL thickness flag = %d\n", M_BL_Thickness);
    //printf("BL spacing flag = %d\n", M_Initial_Normal_Spacing);

    if (surfaceMesh->meshQuickRef.useStartIndex == (int) false &&
        surfaceMesh->meshQuickRef.useListIndex  == (int) false) {

        status = mesh_fillQuickRefList( surfaceMesh);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    Number_of_Nodes = surfaceMesh->numNode;
    Number_of_Surf_Quads = surfaceMesh->meshQuickRef.numQuadrilateral;
    Number_of_Surf_Trias = surfaceMesh->meshQuickRef.numTriangle;

    status = ug_io_malloc_grid (M_BL_Thickness,
                                M_Initial_Normal_Spacing,
                                1, //M_Surf_Grid_BC_Flag,
                                1, //M_Surf_ID_Flag,
                                1, //M_Surf_Reconnection_Flag, // This is needed aflr3 will seg-fault without it
                                0, //M_Vol_ID_Flag,
                                Number_of_Nodes,
                                Number_of_Surf_Quads,
                                Number_of_Surf_Trias,
                                Number_of_Vol_Hexs,
                                Number_of_Vol_Pents_5,
                                Number_of_Vol_Pents_6,
                                Number_of_Vol_Tets,
                                &Surf_Grid_BC_Flag,
                                &Surf_ID_Flag,
                                &Surf_Reconnection_Flag,
                                &Surf_Quad_Connectivity,
                                &Surf_Tria_Connectivity,
                                &Vol_Hex_Connectivity,
                                &Vol_ID_Flag,
                                &Vol_Pent_5_Connectivity,
                                &Vol_Pent_6_Connectivity,
                                &Vol_Tet_Connectivity,
                                &Coordinates,
                                &Initial_Normal_Spacing,
                                &BL_Thickness);
    if (status != EGADS_SUCCESS) goto cleanup;

    if (Coordinates != NULL)
        for (i = 1; i <= Number_of_Nodes; i++){
            Coordinates[i][0] = surfaceMesh->node[i-1].xyz[0] ;
            Coordinates[i][1] = surfaceMesh->node[i-1].xyz[1] ;
            Coordinates[i][2] = surfaceMesh->node[i-1].xyz[2] ;
        }

    for (i = 0; i < surfaceMesh->numElement; i++) {

        // Tri connectivity
        if (Surf_Tria_Connectivity != NULL)
            if (surfaceMesh->element[i].elementType == Triangle) {
                Surf_Tria_Connectivity[triIndex][0] =
                                        surfaceMesh->element[i].connectivity[0];
                Surf_Tria_Connectivity[triIndex][1] =
                                        surfaceMesh->element[i].connectivity[1];
                Surf_Tria_Connectivity[triIndex][2] =
                                        surfaceMesh->element[i].connectivity[2];

                triIndex += 1;
            }

        // Quad connectivity
        if (Surf_Quad_Connectivity != NULL)
            if (surfaceMesh->element[i].elementType == Quadrilateral) {
                Surf_Quad_Connectivity[quadIndex][0] =
                                        surfaceMesh->element[i].connectivity[0];
                Surf_Quad_Connectivity[quadIndex][1] =
                                        surfaceMesh->element[i].connectivity[1];
                Surf_Quad_Connectivity[quadIndex][2] =
                                        surfaceMesh->element[i].connectivity[2];
                Surf_Quad_Connectivity[quadIndex][3] =
                                        surfaceMesh->element[i].connectivity[3];

                quadIndex += 1;
            }
        }

    // BC flag index for Triangles
    triIndex = 1;
    for (i = 0; i < surfaceMesh->numElement; i++) {

        if (surfaceMesh->element[i].elementType != Triangle) continue;
        if (Surf_ID_Flag != NULL)
            Surf_ID_Flag[triIndex] = surfaceMesh->element[i].markerID;

        if (Surf_Grid_BC_Flag != NULL)
            Surf_Grid_BC_Flag[triIndex] = blFlag[i];

        if (Surf_Reconnection_Flag != NULL)
            Surf_Reconnection_Flag[triIndex] = 7; // Do not allow swapping on the surface triangles

        triIndex += 1;
    }

    // BC flag index for Quadrilaterals
    quadIndex = triIndex;
    for (i = 0; i < surfaceMesh->numElement; i++) {

        if (surfaceMesh->element[i].elementType != Quadrilateral) continue;
        if (Surf_ID_Flag != NULL)
            Surf_ID_Flag[quadIndex] = surfaceMesh->element[i].markerID;

        if (Surf_Grid_BC_Flag != NULL)
            Surf_Grid_BC_Flag[quadIndex] = blFlag[i];

        if (Surf_Reconnection_Flag != NULL)
            Surf_Reconnection_Flag[quadIndex] = 7; // Do not allow swapping on the surface triangles

        quadIndex += 1;
    }

    // BL desired thickness at each node
    if ((M_BL_Thickness != 0) && (BL_Thickness != NULL)) {
        for (i = 1; i <= Number_of_Nodes; i++) BL_Thickness[i] = blThickness[i-1];
    }

    // BL desired initial wall spacing at each node
    if ((M_Initial_Normal_Spacing != 0) && (Initial_Normal_Spacing != NULL)) {
        for (i = 1; i <= Number_of_Nodes; i++) Initial_Normal_Spacing[i] = blSpacing[i-1];
    }

    // Loop through Elements and see if marker match
    for (i = 0; i < surfaceMesh->numElement; i++) {

        // Loop through meshing properties and see if boundaryLayerThickness and boundaryLayerSpacing have been specified
        for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

            //If they don't match continue
            if (surfaceMesh->element[i].markerID != meshProp[propIndex].attrIndex) continue;

            // If bcType specified in meshProp
            if ((meshProp[propIndex].bcType != NULL) &&
                (Surf_Grid_BC_Flag != NULL)) {

                bcType = (createBL == (int)true) ? -STD_UG3_GBC : STD_UG3_GBC;
                if      (strncasecmp(meshProp[propIndex].bcType, "Farfield"             ,  8) == 0 ||
                         strncasecmp(meshProp[propIndex].bcType, "Freestream"           , 10) == 0 ||
                         strncasecmp(meshProp[propIndex].bcType, "FARFIELD_UG3_GBC"     , 16) == 0)
                    bcType = FARFIELD_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "Viscous"              ,  7) == 0 ||
                         strncasecmp(meshProp[propIndex].bcType, "-STD_UG3_GBC"         , 12) == 0)
                    bcType = -STD_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "Inviscid"             ,  8) == 0 ||
                         strncasecmp(meshProp[propIndex].bcType, "STD_UG3_GBC"          , 11) == 0)
                    bcType = STD_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "Symmetry"             ,  8) == 0 ||
                         strncasecmp(meshProp[propIndex].bcType, "BL_INT_UG3_GBC"       , 14) == 0)
                    bcType = BL_INT_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_SRC_UG3_GBC"   , 18) == 0)
                    bcType = TRANSP_SRC_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_BL_INT_UG3_GBC", 21) == 0)
                    bcType = TRANSP_BL_INT_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_UG3_GBC"       , 14) == 0)
                    bcType = TRANSP_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "-TRANSP_UG3_GBC"      , 15) == 0)
                    bcType = -TRANSP_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "TRANSP_INTRNL_UG3_GBC", 20) == 0)
                    bcType = TRANSP_INTRNL_UG3_GBC;
                else if (strncasecmp(meshProp[propIndex].bcType, "FIXED_BL_INT_UG3_GBC" , 19) == 0)
                    bcType = FIXED_BL_INT_UG3_GBC;

                // Set face BC flag
                Surf_Grid_BC_Flag[i+1] = bcType;
            }

            if ((surfaceMesh->element[i].elementType != Triangle &&
                 surfaceMesh->element[i].elementType != Quadrilateral) ||
                createBL == (int)false) continue;

            // Get face indexing for Triangles - 1 bias
            if (surfaceMesh->element[i].elementType == Triangle) {
                pointIndex[0] = surfaceMesh->element[i].connectivity[0];
                pointIndex[1] = surfaceMesh->element[i].connectivity[1];
                pointIndex[2] = surfaceMesh->element[i].connectivity[2];
                pointIndex[3] = surfaceMesh->element[i].connectivity[2]; //Repeat last point for simplicity
            }

            if (surfaceMesh->element[i].elementType == Quadrilateral) {
                pointIndex[0] = surfaceMesh->element[i].connectivity[0];
                pointIndex[1] = surfaceMesh->element[i].connectivity[1];
                pointIndex[2] = surfaceMesh->element[i].connectivity[2];
                pointIndex[3] = surfaceMesh->element[i].connectivity[3];
            }

            // Set boundary layer spacing
            if (Initial_Normal_Spacing != NULL) {
                Initial_Normal_Spacing[pointIndex[0]] =
                        meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
                Initial_Normal_Spacing[pointIndex[1]] =
                        meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
                Initial_Normal_Spacing[pointIndex[2]] =
                        meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
                Initial_Normal_Spacing[pointIndex[3]] =
                        meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
            }

            // Set boundary layer thickness
            if (BL_Thickness != NULL) {
                BL_Thickness[pointIndex[0]] =
                      meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
                BL_Thickness[pointIndex[1]] =
                      meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
                BL_Thickness[pointIndex[2]] =
                      meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
                BL_Thickness[pointIndex[3]] =
                      meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
            }

            // Set face BC flag if not already set
            if (Surf_Grid_BC_Flag != NULL)
                if (Surf_Grid_BC_Flag[i+1] == 0 &&
                    meshProp[propIndex].boundaryLayerSpacing > 0 &&
                    meshProp[propIndex].boundaryLayerThickness > 0)
                  Surf_Grid_BC_Flag[i+1] = -STD_UG3_GBC;

            break;
        }
    }

    // Read background grid and function data.
    // This is optional and not required for implementation.
    /*
    if (status== 0 && Program_Flag && strcmp (BG_Grid_File_Name, ""))
    {
        status= ug_io_read_grid_file (BG_Grid_File_Name,
                                   Message_Flag,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   &Number_of_BG_BL_Vol_Tets,
                                   &Number_of_BG_Nodes,
                                   &Number_of_BG_Surf_Quads,
                                   &Number_of_BG_Surf_Trias,
                                   &Number_of_BG_Vol_Hexs,
                                   &Number_of_BG_Vol_Pents_5,
                                   &Number_of_BG_Vol_Pents_6,
                                   &Number_of_BG_Vol_Tets,
                                   &BG_Surf_Grid_BC_Flag,
                                   &BG_Surf_ID_Flag,
                                   &BG_Surf_Reconnection_Flag,
                                   &BG_Surf_Quad_Connectivity,
                                   &BG_Surf_Tria_Connectivity,
                                   &BG_Vol_Hex_Connectivity,
                                   &BG_Vol_ID_Flag,
                                   &BG_Vol_Pent_5_Connectivity,
                                   &BG_Vol_Pent_6_Connectivity,
                                   &BG_Vol_Tet_Connectivity,
                                   &BG_Coordinates,
                                   &BG_Initial_Normal_Spacing,
                                   &BG_BL_Thickness);

        ug_free (BG_Surf_Quad_Connectivity);
        ug_free (BG_Surf_Tria_Connectivity);
        ug_free (BG_Vol_Hex_Connectivity);
        ug_free (BG_Vol_Pent_5_Connectivity);
        ug_free (BG_Vol_Pent_6_Connectivity);

        BG_Surf_Quad_Connectivity = NULL;
        BG_Surf_Tria_Connectivity = NULL;
        BG_Vol_Hex_Connectivity = NULL;
        BG_Vol_Pent_5_Connectivity = NULL;
        BG_Vol_Pent_6_Connectivity = NULL;

        if (Number_of_BG_Vol_Hexs || Number_of_BG_Vol_Pents_5 || Number_of_BG_Vol_Pents_6)
        {
            ug_error_message ("*** ERROR : background grid file may not contain hex, prism or pyramid elements ***");
            status= 199;
        }

        if (status== 0)
        {
            status= ug_io_read_func_file (BG_Func_File_Name,
                                     Message_Flag, 0, 3,
                                     &Number_of_BG_U_Nodes,
                                     &Number_of_BG_U_Scalars,
                                     &Number_of_BG_U_Vectors,
                                     &Number_of_BG_U_Matrixes,
                                     &Number_of_BG_U_Metrics,
                                     &BG_U_Scalar_Labels,
                                     &BG_U_Vector_Labels,
                                     &BG_U_Matrix_Labels,
                                     &BG_U_Metric_Labels,
                                     &BG_U_Scalar_Flags,
                                     &BG_U_Vector_Flags,
                                     &BG_U_Matrix_Flags,
                                     &BG_U_Metric_Flags,
                                     &BG_U_Scalars,
                                     (DOUBLE_1D **) &BG_U_Vectors,
                                     (DOUBLE_1D **) &BG_U_Matrixes,
                                     (DOUBLE_1D **) &BG_U_Metrics);

            ug_io_free_func_flag (BG_U_Scalar_Labels, BG_U_Vector_Labels,
                        BG_U_Matrix_Labels, BG_U_Metric_Labels,
                        BG_U_Scalar_Flags, BG_U_Vector_Flags,
                        BG_U_Matrix_Flags, BG_U_Metric_Flags);

            BG_U_Scalar_Labels = NULL;
            BG_U_Vector_Labels = NULL;
            BG_U_Matrix_Labels = NULL;
            BG_U_Metric_Labels = NULL;
            BG_U_Scalar_Flags = NULL;
            BG_U_Vector_Flags = NULL;
            BG_U_Matrix_Flags = NULL;
            BG_U_Metric_Flags = NULL;

            ug_free (BG_U_Vectors);
            ug_free (BG_U_Matrixes);

            BG_U_Vectors = NULL;
            BG_U_Matrixes = NULL;

            if (status== 0)
            {
                if (Number_of_BG_U_Nodes != Number_of_BG_Nodes)
                {
                    ug_error_message ("*** ERROR : background function must contain the same number of nodes as the background grid ***");
                    status= 199;
                }

                else if (Number_of_BG_U_Scalars > 1 || Number_of_BG_U_Vectors ||
                        Number_of_BG_U_Matrixes || Number_of_BG_U_Metrics > 1)
                {
                    ug_error_message ("*** ERROR : background function file may only contain either isotropic spacing or metric ***");
                    status= 199;
                }

                else if (Number_of_BG_U_Scalars == 0 && Number_of_BG_U_Metrics == 0)
                {
                    ug_error_message ("*** ERROR : background function file contains neither isotropic spacing or metric ***");
                    status= 199;
                }

                else if (status== 0 && Number_of_BG_U_Scalars == 0 && Number_of_BG_U_Metrics == 1)
                {
                    Number_of_BG_U_Scalars = 1;

                    BG_U_Scalars = (DOUBLE_1D *) ug_malloc (&status, (Number_of_BG_Nodes+1) * sizeof (DOUBLE_1D));

                    if (status== 0)
                        dftr3_met_met2df_n (Number_of_BG_Nodes, BG_U_Metrics, BG_U_Scalars);

                    else
                    {
                        ug_error_message ("*** ERROR : unable to allocate scalr function for background function file ***");
                        status= 100199;
                    }
                }

                if (status== 0 && Number_of_BG_U_Metrics == 0 && mmet)
                {
                    ug_error_message ("*** ERROR : if the metric option is on (mmet != 0) then the ***");
                    ug_error_message ("*** ERROR : background function file must include a metric ***");
                    status= 199;
                }

                if (status== 0 && mtr)
                {
                    ug_error_message ("*** ERROR : transformation vector option (mtr != 0) can not allowed ***");
                    ug_error_message ("*** ERROR : with a background grid ***");
                    status= 199;
                }

                BG_Spacing = BG_U_Scalars;
                BG_Metric = BG_U_Metrics;
            }
        }
    }
     */

    // Read source node data.
    // This is optional and not required for implementation.

    /*
    if (status== 0 && Program_Flag && strcmp (Source_Data_File_Name, ""))
        status= ug_io_read_node_file (Source_Data_File_Name,
                                   Message_Flag, &Number_of_Source_Nodes,
                                   &Source_Coordinates, &Source_Spacing,
                                   &Source_Metric);
     */

    // Set case name for debug files.
    // This is required for implementation in parallel mode.

    //if (status== 0 && Program_Flag > 0) ug_case_name (Output_Case_Name);


    // Read tags data file (if found) and set grid BC arguments.
    // This is optional and not required for implementation.
    /*

    if (status== 0 && Program_Flag && Tags_Data_File_Flag == 1)
    {
        status= ug_io_read_tags_file (&prog_argc, &prog_argv,
                                   Input_Case_Name, Message_Flag);
        status= MAX (status, 0);
    }
     */

    // Read periodic surface data file (if found) and set related arguments.
    // This is optional and not required for implementation.
    /*

    if (status== 0 && Program_Flag == 1 && PSDATA_File_Flag && mbl)
    {
        status= ug_io_read_psdata_file (&prog_argc, &prog_argv,
                                     Input_Case_Name, Message_Flag);
        status= MAX (status, 0);
    }
     */

    // Read specified normal spacing data file (if found) and set related
    // arguments.
    // This is optional and not required for implementation.
    /*

    if (status== 0 && Program_Flag == 1 && mbl == -1)
    {
        status= ug_io_read_snsdata_file (&prog_argc, &prog_argv,
                                      Input_Case_Name, Message_Flag);
        status= MAX (status, 0);
    }

     */
    // Set memory reduction output file format flag parameter, mpfrmt.
    //
    // If mpfrmt=1,2,3,4,5 then the output grid file is created within AFLR3 and
    // only the surface grid is passed back from AFLR3.
    //
    // If mpfrmt=0 then the output grid file is not created within AFLR3 and
    // the complete volume grid is passed back from AFLR3.
    //
    // This is optional and not required for implementation.

    /*    if (status== 0 && Program_Flag == 1 && meshInput.outputFileName != NULL)// && mbl)
    {

        strcpy(Output_Grid_File_Name, meshInput.outputFileName);

        if (strcasecmp(meshInput.outputFormat,"UGRID")  == 0) {

            if (meshInput.outputASCIIFlag != 0)  strcpy(Output_Grid_File_Type_Suffix, ".ugrid");
            else                                  strcpy(Output_Grid_File_Type_Suffix, ".lb8.ugrid");


        } else if (strcasecmp(meshInput.outputFormat,"SU2")    == 0) {
            strcpy(Output_Grid_File_Name, "_null_");

        } else if (strcasecmp(meshInput.outputFormat,"CGNS")    == 0) strcpy(Output_Grid_File_Type_Suffix, ".cgns");
          else if (strcasecmp(meshInput.outputFormat,"NASTRAN")== 0)  strcpy(Output_Grid_File_Type_Suffix, ".nas");


        // Concat filename and suffix
        strcat(Output_Grid_File_Name, Output_Grid_File_Type_Suffix);

        // Pull just the suffix back off
        if (status == 0) status = ug_io_file_type_suffix (Output_Grid_File_Name, 3,
                                                          Output_Grid_File_Type_Suffix);

        // Get file format
        if (status== 0) status = ug_io_file_format (Output_Grid_File_Name,
                                                   &Output_Grid_File_Format);


        if (status== 0)
        {
            if (strcmp (Output_Grid_File_Type_Suffix, ".ugrid") == 0)
                if (Output_Grid_File_Format == UG_IO_FORMATTED) mpfrmt =  1;
                else if (Output_Grid_File_Format == UG_IO_BINARY_DOUBLE) mpfrmt =  2;
                else if (Output_Grid_File_Format == -UG_IO_BINARY_DOUBLE) mpfrmt = 3;
                else if (Output_Grid_File_Format == UG_IO_UNFORMATTED_DOUBLE) mpfrmt =  4;
                else if (Output_Grid_File_Format == -UG_IO_UNFORMATTED_DOUBLE) mpfrmt =  5;
                else mpfrmt = 0;
            else
                mpfrmt = 0;
        }


    }*/

    // Set write mesh flag - do not write out the mesh internally
    status = ug_add_flag_arg ((char *) "mpfrmt=0", &prog_argc, &prog_argv);
    if (status != 0) goto cleanup;

    // Generate unstructured volume grid.
    if (status== 0 && Program_Flag == 1) {
        status = aflr3_grid_generator (prog_argc,
                                       prog_argv,
                                       Message_Flag,
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
                                       &Initial_Normal_Spacing,
                                       &BL_Thickness,
                                       &BG_Coordinates,
                                       &BG_Spacing,
                                       &BG_Metric,
                                       &Source_Coordinates,
                                       &Source_Spacing,
                                       &Source_Metric);
        if (status != EGADS_SUCCESS) {
            FILE *fp;

            printf("Volume grid generation failed!\n");

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

            if ((Surf_ID_Flag != NULL) && (Surf_Error_Flag != NULL))
                for (i = 0; i < Number_of_Surf_Trias; i++) {
                    if (i % 5 == 0) fprintf( fp, "\n");
                    if (Surf_Error_Flag[i+1] < 0)
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

    } else if (status== 0 && Program_Flag == 2) {

        // Generate unstructured volume grid using simplified interface.
        // This is optional and not required for implementation.

        status= aflr3_tet_grid (Interior_Flag,
                                Message_Flag,
                                Reconnection_Flag,
                                &Number_of_Nodes,
                                &Number_of_Surf_Trias,
                                &Number_of_Vol_Tets,
                                &Surf_Tria_Connectivity,
                                &Vol_Tet_Connectivity,
                                &Coordinates);

        if (status != EGADS_SUCCESS) {
            printf("Volume grid generation failed!\n");

            ug3_write_surf_grid_error_file (Output_Case_Name,
                                            status,
                                            Number_of_Nodes,
                                            Number_of_Surf_Trias,
                                            Surf_Error_Flag,
                                            Surf_Grid_BC_Flag,
                                            Surf_ID_Flag,
                                            Surf_Tria_Connectivity,
                                            Coordinates);
            goto cleanup;
        }
    }

    // Skip remaining tasks except for freeing memory
    // if BL parameter checking is on.
    // This is optional and not required for implementation.

    //if (mblchki == 1) Program_Flag = -1;

    /*    // Save AFLR3 argument file as a local file.
    // This is optional and not required for implementation.

    if (status== 0 && Program_Flag > 0 && Call == Number_of_Calls)
    {
        strcpy (Arg_File_Name, Output_Case_Name);
        strcat (Arg_File_Name, ".aflr3.arg");

        ug_backup_file (Arg_File_Name);

        ug_write_arg_file (argc_copy, argv_copy, Arg_File_Name);
    }*/

    // If memory reduction output grid file is already of the correct type and
    // format then skip remaining tasks except for freeing memory.
    // This is optional and not required for implementation.

    //if (mpfrmt) Program_Flag = -1;

    // Ignore non-BL elements.
    // This is optional and not required for implementation.

    /*    if (status== 0 && Program_Flag >= 0 && Write_BL_Only == 1)
        Number_of_Vol_Tets = Number_of_BL_Vol_Tets;*/

    // Allocate and initialize volume element ID flag.
    // This is optional and not required for implementation.

    /*    if (status== 0 && Program_Flag > 0 && Set_Vol_ID_Flag == 1)
    {
        Number_of_Vol_Elems = Number_of_Vol_Hexs + Number_of_Vol_Pents_5
                    + Number_of_Vol_Pents_6 + Number_of_Vol_Tets;

        status= 0;

        Vol_ID_Flag = (INT_1D *) ug_malloc (&status,
                                    (Number_of_Vol_Elems+1)
     * sizeof (INT_1D));

        if (status)
            status= 100199;

        ug_set_int (1, Number_of_Vol_Elems, -123456, Vol_ID_Flag);
    }*/

    // Write the mesh to disk
    if (status == 0) {

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
                                       Initial_Normal_Spacing,
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
        if (status != CAPS_SUCCESS) {
            printf("Error occurred while transferring volume grid, status = %d\n",
                   status);
            goto cleanup;
        }
    }
    /*if (status== 0 &&
        (Program_Flag > 0 ||
                (Program_Flag == 0 && strcmp (Output_Grid_File_Name, "_null_") == 0)))
        status= ug_gq_write (Output_Case_Name,
                          GQ_Max_Dist_Increments, GQ_Surf_Measure_Flag,
                          GQ_Vol_Measure_Flag, Message_Flag,
                          Number_of_Surf_Quads, Number_of_Surf_Trias,
                          Number_of_Vol_Hexs, Number_of_Vol_Pents_5,
                          Number_of_Vol_Pents_6, Number_of_Vol_Tets,
                          Surf_Quad_Connectivity, Surf_Tria_Connectivity,
                          Vol_Hex_Connectivity,
                          Vol_Pent_5_Connectivity, Vol_Pent_6_Connectivity,
                          Vol_Tet_Connectivity, Coordinates);*/

    // Remove the temporary grid created by AFLR
    remove(".tmp.b8.ugrid");

    status = CAPS_SUCCESS;
cleanup:

    if (status != CAPS_SUCCESS)
      printf("Error: Premature exit in aflr3_Volume_Mesh status = %d\n", status);

    EG_free(meshInputString); meshInputString = NULL;

    // Free program arguements
    ug_free_argv(prog_argv); prog_argv = NULL;

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
                    Initial_Normal_Spacing,
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
        Initial_Normal_Spacing = NULL;
        BL_Thickness = NULL;

        ug_free(Surf_Error_Flag);
        Surf_Error_Flag= NULL;

        ug_free (BG_Vol_Tet_Neigbors);
        ug_free (BG_Vol_Tet_Connectivity); // only needed if there is an error before aflr3
        ug_free (BG_Coordinates); // only needed if there is an error before aflr3
        ug_free (BG_Spacing);
        ug_free (BG_Metric);

        BG_Vol_Tet_Neigbors = NULL;
        BG_Vol_Tet_Connectivity = NULL;
        BG_Coordinates = NULL;
        BG_Spacing = NULL;
        BG_Metric = NULL;

        ug_free(BG_U_Scalars); // only needed if there is an error before aflr3
        ug_free(BG_U_Metrics); // only needed if there is an error before aflr3

        BG_U_Scalars = NULL;
        BG_U_Metrics = NULL;

        ug_io_free_node(Source_Coordinates, Source_Spacing, Source_Metric);

        Source_Coordinates = NULL;
        Source_Spacing = NULL;
        Source_Metric = NULL;

        // Shut off memory and file status monitors and close output file.
        // This is required for implementation!
        ug_shutdown ();

        return status;
}
