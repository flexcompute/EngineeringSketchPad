/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Abaqus AIM
 *
 *     Written by Dr. Ryan Durscher AFRL/RQVC
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

/*! \mainpage Introduction
 *
 * \section overviewAbaqus Abaqus AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with the finite element structural solver Abaqus \cite Abaqus.
 *
 * Current issues include:
 *  - A thorough bug testing needs to be undertaken.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsAbaqus and
 * \ref aimOutputsAbaqus and \ref attributeAbaqus, respectively.
 *
 * Details of the AIM's shareable data structures are outlined in \ref sharableDataAbaqus if
 * connecting this AIM to other AIMs in a parent-child like manner.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferAbaqus
 *
 */

 /* \section abaqusExamples Examples
 * Example problems using the Abaqus AIM may be found at \ref examplesAbaqus .
 *  - \ref abaqusSingleLoadEx
 *  - \ref abaqusMultipleLoadEx
 *  - \ref abaqusModalEx
 *  - \ref abaqusOptEx
 *  - \ref abaqusCompositeEx
 *  - \ref abaqusCompOptimizationEx
 */

/*! \page attributeAbaqus Abaqus AIM attributes
 * The following list of attributes are required for the Abaqus AIM inside the geometry input.
 *
 * - <b> capsGroup</b> This is a name assigned to any geometric body.  This body could be a solid, surface, face, wire, edge or node.
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsGroup $Wing</c>.
 *
 * - <b> capsLoad</b> This is a name assigned to any geometric body where a load is applied.  This attribute was separated from the <c>capsGroup</c>
 * attribute to allow the user to define a local area to apply a load on without adding multiple <c>capsGroup</c> attributes.
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsLoad $force</c>.
 *
 * - <b> capsConstraint</b> This is a name assigned to any geometric body where a constraint/boundary condition is applied.
 * This attribute was separated from the <c>capsGroup</c> attribute to allow the user to define a local area to apply a boundary condition
 * without adding multiple <c>capsGroup</c> attributes. Recall that a string in ESP starts with a $.  For example, attribute <c>capsConstraint $fixed</c>.
 *
 * - <b> capsIgnore</b> It is possible that there is a geometric body (or entity) that you do not want the Astros AIM to pay attention to when creating
 * a finite element model. The capsIgnore attribute allows a body (or entity) to be in the geometry and ignored by the AIM.  For example,
 * because of limitations in OpenCASCADE a situation where two edges are overlapping may occur; capsIgnore allows the user to only
 *  pay attention to one of the overlapping edges.
 *
 */

 /* - <b> capsConnect</b> This is a name assigned to any geometric body where the user wishes to create
 * "fictitious" connections such as springs, dampers, and/or rigid body connections to. The user must manually specify
 * the connection between two <c>capsConnect</c> entities using the "Connect" tuple (see \ref aimInputsAbaqus).
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsConnect $springStart</c>.
 *
 * - <b> capsConnectLink</b> Similar to <c>capsConnect</c>, this is a name assigned to any geometric body where
 * the user wishes to create "fictitious" connections to. A connection is automatically made if a <c>capsConnectLink</c>
 * matches a <c>capsConnect</c> group. Again further specifics of the connection are input using the "Connect"
 * tuple (see \ref aimInputsAbaqus). Recall that a string in ESP starts with a $.
 * For example, attribute <c>capsConnectLink $springEnd</c>.
 *
 * - <b> capsBound </b> This is used to mark surfaces on the structural grid in which data transfer with an external
 * solver will take place. See \ref dataTransferAbaqus for additional details.
 *
 *
 */

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"    // Meshing utilities
#include "miscUtils.h"    // Miscellaneous utilities
#include "feaUtils.h"     // FEA utilities
#include "abaqusUtils.h" // Abaqus utilities

#ifdef WIN32
#define getcwd     _getcwd
#define snprintf   _snprintf
#define strcasecmp stricmp
#define PATH_MAX   _MAX_PATH
#else
#include <unistd.h>
#include <limits.h>
#endif

//#define DEBUG
//#define WRITE_VONMISES_VTK

enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Property,
  Material,
  Constraint,
  Load,
  Analysix,
  Analysis_Type,
  Support,
  Connect,
  Mesh_Morph,
  Mesh,
  NUMINPUT = Mesh              /* Total number of inputs */
};

 enum aimOutputs
 {
//   EigenValue = 1,           /* index is 1-based */
//   EigenRadian,
//   EigenFrequency,
//   EigenGeneralMass,
   Tmax = 1,
   T1max,
   T2max,
   T3max,
   NUMOUTPUT = T3max            /* Total number of outputs */
 };

typedef struct {

    // Project name
    char *projectName; // Project name

    feaUnitsStruct units; // units system

    feaProblemStruct feaProblem;

    // Attribute to index map
    mapAttrToIndexStruct attrMap;

    // Attribute to constraint index map
    mapAttrToIndexStruct constraintMap;

    // Attribute to load index map
    mapAttrToIndexStruct loadMap;

    // Attribute to transfer map
    mapAttrToIndexStruct transferMap;

    // Attribute to connect map
    mapAttrToIndexStruct connectMap;

    // Attribute to response map
    mapAttrToIndexStruct responseMap;

    // Mesh holders
    int numMesh;
    meshStruct *feaMesh;

    double Tmax, T1max, T2max, T3max;

} aimStorage;


static int initiate_aimStorage(void *aimInfo, aimStorage *abaqusInstance) {

    int status;

    // Set initial values for abaqusInstance
    abaqusInstance->projectName = NULL;

    status = initiate_feaUnitsStruct(&abaqusInstance->units);
    AIM_STATUS(aimInfo, status);

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&abaqusInstance->attrMap);
    AIM_STATUS(aimInfo, status);

    // Container for attribute to constraint index map
    status = initiate_mapAttrToIndexStruct(&abaqusInstance->constraintMap);
    AIM_STATUS(aimInfo, status);

    // Container for attribute to load index map
    status = initiate_mapAttrToIndexStruct(&abaqusInstance->loadMap);
    AIM_STATUS(aimInfo, status);

    // Container for transfer to index map
    status = initiate_mapAttrToIndexStruct(&abaqusInstance->transferMap);
    AIM_STATUS(aimInfo, status);

    // Container for connect to index map
    status = initiate_mapAttrToIndexStruct(&abaqusInstance->connectMap);
    AIM_STATUS(aimInfo, status);

    // Container for response to index map
    status = initiate_mapAttrToIndexStruct(&abaqusInstance->responseMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_feaProblemStruct(&abaqusInstance->feaProblem);
    AIM_STATUS(aimInfo, status);

    // Mesh holders
    abaqusInstance->numMesh = 0;
    abaqusInstance->feaMesh = NULL;

    abaqusInstance->Tmax = 0;
    abaqusInstance->T1max = 0;
    abaqusInstance->T2max = 0;
    abaqusInstance->T3max = 0;

cleanup:
    return status;
}

static int destroy_aimStorage(aimStorage *abaqusInstance) {

    int status;
    int i;

    status = destroy_feaUnitsStruct(&abaqusInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&abaqusInstance->attrMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to constraint index map
    status = destroy_mapAttrToIndexStruct(&abaqusInstance->constraintMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to load index map
    status = destroy_mapAttrToIndexStruct(&abaqusInstance->loadMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Transfer to index map
    status = destroy_mapAttrToIndexStruct(&abaqusInstance->transferMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Connect to index map
    status = destroy_mapAttrToIndexStruct(&abaqusInstance->connectMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Response to index map
    status = destroy_mapAttrToIndexStruct(&abaqusInstance->responseMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (abaqusInstance->feaMesh != NULL) {

        for (i = 0; i < abaqusInstance->numMesh; i++) {
            status = destroy_meshStruct(&abaqusInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        EG_free(abaqusInstance->feaMesh);
    }

    abaqusInstance->feaMesh = NULL;
    abaqusInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&abaqusInstance->feaProblem);
    if (status != CAPS_SUCCESS)  printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    // NULL projetName
    abaqusInstance->projectName = NULL;

    abaqusInstance->Tmax = 0;
    abaqusInstance->T1max = 0;
    abaqusInstance->T2max = 0;
    abaqusInstance->T3max = 0;

    return CAPS_SUCCESS;
}

static int checkAndCreateMesh(void* aimInfo, aimStorage *abaqusInstance)
{
  // Function return flag
  int status = CAPS_SUCCESS;

  status = fea_createMesh(aimInfo,
                          NULL,
                          0,
                          0,
                          (int)false,
                          &abaqusInstance->attrMap,
                          &abaqusInstance->constraintMap,
                          &abaqusInstance->loadMap,
                          &abaqusInstance->transferMap,
                          &abaqusInstance->connectMap,
                          &abaqusInstance->responseMap,
                          NULL,
                          &abaqusInstance->numMesh,
                          &abaqusInstance->feaMesh,
                          &abaqusInstance->feaProblem );
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}

/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS;
    int i;
    int *ints=NULL;
    char **strs=NULL;

    aimStorage *abaqusInstance=NULL;

    #ifdef DEBUG
        printf("\n abaqusAIM/aimInitialize   ngIn = %d!\n", ngIn);
    #endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == 1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate */
    *nFields = 4;

    /* specify the name of each field variable */
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);

    strs[0]  = EG_strdup("Displacement");
    strs[1]  = EG_strdup("EigenVector");
    strs[2]  = EG_strdup("EigenVector_*");
    strs[3]  = EG_strdup("Pressure");
    for (i = 0; i < *nFields; i++)
      if (strs[i] == NULL) { status = EGADS_MALLOC; goto cleanup; }
    *fnames  = strs;

    /* specify the dimension of each field variable */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);
    ints[0]  = 3;
    ints[1]  = 3;
    ints[2]  = 3;
    ints[3]  = 1;
    *franks  = ints;
    ints = NULL;

    /* specify if a field is an input field or output field */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = FieldOut;
    ints[1]  = FieldOut;
    ints[2]  = FieldOut;
    ints[3]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;

    // Allocate abaqusInstance
    AIM_ALLOC(abaqusInstance, 1, aimStorage, aimInfo, status);
    *instStore = abaqusInstance;

    // Initialize instance storage
    (void) initiate_aimStorage(aimInfo, abaqusInstance);

cleanup:
    if (status != CAPS_SUCCESS) {
      /* release all possibly allocated memory on error */
      if (*fnames != NULL)
        for (i = 0; i < *nFields; i++) AIM_FREE((*fnames)[i]);
      AIM_FREE(*franks);
      AIM_FREE(*fInOut);
      AIM_FREE(*fnames);
      AIM_FREE(*instStore);
      *nFields = 0;
    }

    return status;
}

int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{

    /*! \page aimInputsAbaqus AIM Inputs
     * The following list outlines the Abaqus inputs along with their default value available
     * through the AIM interface. Unless noted these values will be not be linked to
     * any parent AIMs with variables of the same name.
     */

    int status = CAPS_SUCCESS;

    #ifdef DEBUG
        printf(" abaqusAIM/aimInputs instance = %d  index = %d!\n", inst, index);
    #endif

    *ainame = NULL;

    // Abaqus Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("abaqus_CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsAbaqus
         * - <B> Proj_Name = "Abaqus_CAPS"</B> <br>
         * This corresponds to the project name used for file naming.
         */

    } else if (index == Property) {
        *ainame              = EG_strdup("Property");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAbaqus
         * - <B> Property = NULL</B> <br>
         * Property tuple used to input property information for the model, see \ref feaProperty for additional details.
         */
    } else if (index == Material) {
        *ainame              = EG_strdup("Material");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAbaqus
         * - <B> Material = NULL</B> <br>
         * Material tuple used to input material information for the model, see \ref feaMaterial for additional details.
         */
    } else if (index == Constraint) {
        *ainame              = EG_strdup("Constraint");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAbaqus
         * - <B> Constraint = NULL</B> <br>
         * Constraint tuple used to input constraint information for the model, see \ref feaConstraint for additional details.
         */
    } else if (index == Load) {
        *ainame              = EG_strdup("Load");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAbaqus
         * - <B> Load = NULL</B> <br>
         * Load tuple used to input load information for the model, see \ref feaLoad for additional details.
         */
    } else if (index == Analysix) {
        *ainame              = EG_strdup("Analysis");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAbaqus
         * - <B> Analysis = NULL</B> <br>
         * Analysis tuple used to input analysis/case information for the model, see \ref feaAnalysis for additional details.
         */
    } else if (index == Analysis_Type) {
        *ainame              = EG_strdup("Analysis_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Modal");
        defval->lfixed       = Change;

        /*! \page aimInputsAbaqus
         * - <B> Analysis_Type = "Modal"</B> <br>
         * Type of analysis to generate files for, options include "Modal", "Static".
         */
    } else if (index == Support) {
        *ainame              = EG_strdup("Support");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /* \page aimInputsAbaqus
         * - <B> Support = NULL</B> <br>
         * Support tuple used to input support information for the model, see \ref feaSupport for additional details.
         */
    } else if (index == Connect) {
        *ainame              = EG_strdup("Connect");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /* \page aimInputsAbaqus
         * - <B> Connect = NULL</B> <br>
         * Connect tuple used to define connection to be made in the, see \ref feaConnection for additional details.
         */

    } else if (index == Mesh_Morph) {
        *ainame              = EG_strdup("Mesh_Morph");
        defval->type         = Boolean;
        defval->lfixed       = Fixed;
        defval->vals.integer = (int) false;
        defval->dim          = Scalar;
        defval->nullVal      = NotNull;

        /*! \page aimInputsAbaqus
         * - <B> Mesh_Morph = False</B> <br>
         * Project previous surface mesh onto new geometry.
         */

    } else if (index == Mesh) {
        *ainame             = AIM_NAME(Mesh);
        defval->type        = PointerMesh;
        defval->dim         = Vector;
        defval->lfixed      = Change;
        defval->sfixed      = Change;
        defval->vals.AIMptr = NULL;
        defval->nullVal     = IsNull;

        /*! \page aimInputsAbaqus
         * - <B>Mesh = NULL</B> <br>
         * A Mesh link.
         */

    } else {
      AIM_ERROR(aimInfo, "Unknown input index %d", index);
      status = CAPS_BADINDEX;
      goto cleanup;
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}


// ********************** AIM Function Break *****************************
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
    int status; // Function return status

    // Analysis information
    const char *analysisType = NULL;

    aimStorage *abaqusInstance;

    abaqusInstance = (aimStorage *) instStore;
    AIM_NOTNULL(aimInputs, aimInfo, status);

    if (aimInputs[Mesh-1].nullVal == IsNull &&
        aimInputs[Mesh_Morph-1].vals.integer == (int) false) {
        AIM_ANALYSISIN_ERROR(aimInfo, Mesh, "'Mesh' input must be linked to an output 'Surface_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get project name
    abaqusInstance->projectName = aimInputs[Proj_Name-1].vals.string;

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    // Get FEA mesh if we don't already have one
    status = checkAndCreateMesh(aimInfo, abaqusInstance);
    AIM_STATUS(aimInfo, status);

    // Note: Setting order is important here.
    // 1. Materials should be set before properties.
    // 2. Coordinate system should be set before mesh and loads
    // 3. Mesh should be set before loads, constraints, supports, and connections
    // 4. Constraints and loads should be set before analysis
    // 5. Optimization should be set after properties, but before analysis

    // Set material properties
    if (aimInputs[Material-1].nullVal == NotNull) {
        status = fea_getMaterial(aimInfo,
                                 aimInputs[Material-1].length,
                                 aimInputs[Material-1].vals.tuple,
                                 &abaqusInstance->units,
                                 &abaqusInstance->feaProblem.numMaterial,
                                 &abaqusInstance->feaProblem.feaMaterial);
        AIM_STATUS(aimInfo, status);
    } else printf("Material tuple is NULL - No materials set\n");

    // Set property properties
    if (aimInputs[Property-1].nullVal == NotNull) {
        status = fea_getProperty(aimInfo,
                                 aimInputs[Property-1].length,
                                 aimInputs[Property-1].vals.tuple,
                                 &abaqusInstance->attrMap,
                                 &abaqusInstance->units,
                                 &abaqusInstance->feaProblem);
        AIM_STATUS(aimInfo, status);


        // Assign element "subtypes" based on properties set
        status = fea_assignElementSubType(abaqusInstance->feaProblem.numProperty,
                                          abaqusInstance->feaProblem.feaProperty,
                                          &abaqusInstance->feaProblem.feaMesh);
        AIM_STATUS(aimInfo, status);

    } else printf("Property tuple is NULL - No properties set\n");

    // Set constraint properties
    if (aimInputs[Constraint-1].nullVal == NotNull) {
        status = fea_getConstraint(aimInfo,
                                   aimInputs[Constraint-1].length,
                                   aimInputs[Constraint-1].vals.tuple,
                                   &abaqusInstance->constraintMap,
                                   &abaqusInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Constraint tuple is NULL - No constraints applied\n");

    // Set support properties
    if (aimInputs[Support-1].nullVal == NotNull) {
        status = fea_getSupport(aimInputs[Support-1].length,
                                aimInputs[Support-1].vals.tuple,
                                &abaqusInstance->constraintMap,
                                &abaqusInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Support tuple is NULL - No supports applied\n");

    // Set connection properties
    if (aimInputs[Connect-1].nullVal == NotNull) {
        status = fea_getConnection(aimInfo,
                                   aimInputs[aim_getIndex(aimInfo, "Connect", ANALYSISIN)-1].length,
                                   aimInputs[aim_getIndex(aimInfo, "Connect", ANALYSISIN)-1].vals.tuple,
                                   &abaqusInstance->connectMap,
                                   &abaqusInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Connect tuple is NULL - Using defaults\n");

    // Set load properties
    if (aimInputs[Load-1].nullVal == NotNull) {
        status = fea_getLoad(aimInfo,
                             aimInputs[Load-1].length,
                             aimInputs[Load-1].vals.tuple,
                             &abaqusInstance->loadMap,
                             &abaqusInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Load tuple is NULL - No loads applied\n");

    // Set analysis settings
    if (aimInputs[Analysix-1].nullVal == NotNull) {
        status = fea_getAnalysis(aimInfo,
                                 aimInputs[Analysix-1].length,
                                 aimInputs[Analysix-1].vals.tuple,
                                 &abaqusInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else {
       printf("Analysis tuple is NULL\n"); // Its ok to not have an analysis tuple we will just create one

       status = fea_createDefaultAnalysis(aimInfo, &abaqusInstance->feaProblem, analysisType);
       AIM_STATUS(aimInfo, status);
   }

cleanup:
    return status;
}


// ********************** AIM Function Break *****************************
int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{

    int i; // Indexing

    int status; // Status return

    int *tempIntegerArray = NULL; // Temporary array to store a list of integers

    // Optimization Information
    //char *objectiveMinMax = NULL, *objectiveResp = NULL;

    // File format information
    //char *tempString = NULL, *delimiter = NULL;

    // File IO
    char filename[PATH_MAX]; // Output file name
    FILE *fp = NULL; // Output file pointer

    const char *ext[] = {".com", ".dat", ".inp", ".msg", ".odb", ".prt", ".sta"};

    // Load information
    feaLoadStruct *feaLoad = NULL;  // size = [numLoad]

    const aimStorage *abaqusInstance;

    abaqusInstance = (const aimStorage *) instStore;
    AIM_NOTNULL(aimInputs, aimInfo, status);

    /* remove previous files */
    for (i = 0; i < sizeof(ext)/sizeof(char*); i++) {
      snprintf(filename, PATH_MAX, "%s%s", abaqusInstance->projectName, ext[i]);
      aim_rmFile(aimInfo, filename);
    }

    if (abaqusInstance->feaProblem.numLoad > 0) {
      AIM_ALLOC(feaLoad, abaqusInstance->feaProblem.numLoad, feaLoadStruct, aimInfo, status);
      for (i = 0; i < abaqusInstance->feaProblem.numLoad; i++) initiate_feaLoadStruct(&feaLoad[i]);
      for (i = 0; i < abaqusInstance->feaProblem.numLoad; i++) {
        status = copy_feaLoadStruct(aimInfo, &abaqusInstance->feaProblem.feaLoad[i], &feaLoad[i]);
        AIM_STATUS(aimInfo, status);

        if (feaLoad[i].loadType == PressureExternal) {

          // Transfer external pressures from the AIM discrObj
          status = fea_transferExternalPressure(aimInfo,
                                                &abaqusInstance->feaProblem.feaMesh,
                                                &feaLoad[i]);
          AIM_STATUS(aimInfo, status);
        }
      }
    }

    // Write Abaqus Mesh
    status = mesh_writeAbaqus(aimInfo,
                              abaqusInstance->projectName,
                              1,
                              &abaqusInstance->feaProblem.feaMesh,
                              1.0); // Scale factor for coordinates
    AIM_STATUS(aimInfo, status);

    // Write Abaqus input file
    snprintf(filename, PATH_MAX, "%s.inp", abaqusInstance->projectName);

    printf("\nWriting Abaqus instruction file....\n");
    fp = aim_fopen(aimInfo, filename, "w");

    if (fp == NULL) {
      printf("Unable to open file: %s\n", filename);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // Heading
    fprintf(fp, "*HEADING\n");
    fprintf(fp, "CAPS generated problem for Abaqus\n");
    //fprintf(fp, "*FILE FORMAT, ASCII\n");


/*
    // Analysis type
    analysisType = aimInputs[aim_getIndex(aimInfo, "Analysis_Type", ANALYSISIN)-1].vals.string;

    if     (strcasecmp(analysisType, "Modal")         == 0) fprintf(fp, "SOL 3\n");
    else if(strcasecmp(analysisType, "Static")        == 0) fprintf(fp, "SOL 1\n");
    else if(strcasecmp(analysisType, "Craig-Bampton") == 0) fprintf(fp, "SOL 31\n");
    else {
        printf("Unrecognized \"Analysis_Type\", %s, defaulting to \"Modal\" analysis\n", analysisType);
        analysisType = "Modal";
        fprintf(fp, "SOL 3\n");
    }
*/

    // Model definition
    fprintf(fp, "**\n**Model Definition\n**\n");

    // Include mesh file
    //fprintf(fp, "*PART, name=TEST\n");
    fprintf(fp,"*INCLUDE, INPUT=%s_Mesh.inp\n", abaqusInstance->projectName);

    // Write all sets
    status = abaqus_writeAllSets(aimInfo, fp, &abaqusInstance->feaProblem);
    AIM_STATUS(aimInfo, status);

    // Properties
    fprintf(fp, "**\n**Properties\n**\n");
    printf("Writing properties\n");
    for (i = 0; i < abaqusInstance->feaProblem.numProperty; i++) {
        status = abaqus_writePropertyCard(fp, &abaqusInstance->feaProblem.feaProperty[i]);
        AIM_STATUS(aimInfo, status);
    }
    //fprintf(fp, "*END PART\n");

    // Materials
    fprintf(fp, "**\n**Materials\n**\n");
    printf("Writing materials\n");
    for (i = 0; i < abaqusInstance->feaProblem.numMaterial; i++) {
        status = abaqus_writeMaterialCard(fp, &abaqusInstance->feaProblem.feaMaterial[i]);
        AIM_STATUS(aimInfo, status);
    }
    // Constraints
    fprintf(fp, "**\n**Constraints\n**\n");
    printf("Writing constraints\n");
    for (i = 0; i < abaqusInstance->feaProblem.numConstraint; i++) {
        status = abaqus_writeConstraintCard(fp, &abaqusInstance->feaProblem.feaConstraint[i]);
        AIM_STATUS(aimInfo, status);
    }

    // Analysis Cards - Eigenvalue included
    fprintf(fp, "**\n**Steps\n**\n");
    printf("Writing analysis (steps)\n");
    for (i = 0; i < abaqusInstance->feaProblem.numAnalysis; i++) {
        status = abaqus_writeAnalysisCard(aimInfo, fp,
                                          abaqusInstance->feaProblem.numLoad,
                                          feaLoad,
                                          &abaqusInstance->feaProblem.feaAnalysis[i],
                                          &abaqusInstance->feaProblem.feaMesh);
        AIM_STATUS(aimInfo, status);

        if (i < abaqusInstance->feaProblem.numAnalysis-1) fprintf(fp, "**\n");
    }

    //fprintf(fp, "*RESTART, WRITE, FREQUENCY=1\n");
    //fprintf(fp, "*OUTPUT, FIELD\n");


//    // Testing
//    int arr[18] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
//    char **string;
//    string = (char **) EG_alloc(3*sizeof(char));
//    string[0] = EG_strdup("a");
//    string[1] = EG_strdup("b");
//    string[2] = EG_strdup("c");
//
//    status = abaqus_writeNodeSet(fp, "test1", 18, arr, 0, NULL);
//    AIM_STATUS(aimInfo, status);
//    status = abaqus_writeNodeSet(fp, "test2", 15, arr, 3, string);
//    AIM_STATUS(aimInfo, status);
//
//    status = abaqus_writeElementSet(fp, "test3", 18, arr, 0, NULL);
//    AIM_STATUS(aimInfo, status);
//    status = abaqus_writeElementSet(fp, "test4", 15, arr, 3, string);
//    AIM_STATUS(aimInfo, status);
//

    status = CAPS_SUCCESS;

cleanup:
    if (feaLoad != NULL) {
        for (i = 0; i < abaqusInstance->feaProblem.numLoad; i++) {
            destroy_feaLoadStruct(&feaLoad[i]);
        }
        AIM_FREE(feaLoad);
    }

    AIM_FREE(tempIntegerArray);

    if (fp != NULL) fclose(fp);

    return status;
}


// ********************** AIM Function Break *****************************
/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int  status;
  int i;
  capsValue val, *analysisType;
  char filename[PATH_MAX]; // File to open
  //char extFIL[] = ".fil";
  char extDAT[] = ".dat";

  int *elementIDs = NULL;
  double *elemData = NULL;
  int *n2a=NULL;

  int *nodeIDs = NULL;
  DOUBLE_6 *dispMatrix = NULL;
  double TT=0;

#ifdef WRITE_VONMISES_VTK
  capsValue vonMises;
  DOUBLE_6 *U = NULL;
#endif

  FILE *fp = NULL; // File pointer
  aimStorage *abaqusInstance;

#ifdef DEBUG
  printf(" abaqusAIM/aimPostAnalysis!\n");
#endif
  abaqusInstance = (aimStorage *) instStore;

  // Open abaqus *.dat file
  snprintf(filename, PATH_MAX, "%s%s", abaqusInstance->projectName, extDAT);

  fp = aim_fopen(aimInfo, filename, "r");
  if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open Output file: %s!", filename);
      status = CAPS_IOERR;
      goto cleanup;
  }
  fclose(fp); fp = NULL;

  status = aim_getValue(aimInfo, Analysis_Type, ANALYSISIN, &analysisType);
  AIM_STATUS(aimInfo, status);

  // Static analysis - only
  if (strcasecmp(analysisType->vals.string, "Static") == 0) {

    status = aim_initValue(&val);
    AIM_STATUS(aimInfo, status);

    /* read stress */
    status = abaqus_readDATMises(aimInfo, filename,
                                 abaqusInstance->feaProblem.feaMesh.numElement,
                                 &elementIDs, &elemData);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(elementIDs, aimInfo, status);
    AIM_NOTNULL(elemData, aimInfo, status);

    status = mesh_gridAvg(aimInfo, &abaqusInstance->feaProblem.feaMesh,
                          abaqusInstance->feaProblem.feaMesh.numElement,
                          elementIDs, 1, elemData,
                          &val.vals.reals);
    AIM_STATUS(aimInfo, status);

    val.dim  = Vector;
    val.type = Double;
    val.nrow = abaqusInstance->feaProblem.feaMesh.numNode;
    val.ncol = 1;
    val.length = val.nrow*val.ncol;

#ifdef WRITE_VONMISES_VTK
    vonMises = val;
#endif

    /* create the dynamic output */
    status = aim_makeDynamicOutput(aimInfo, "vonMises_Grid", &val);
    AIM_STATUS(aimInfo, status);

    snprintf(filename, PATH_MAX, "%s%s", abaqusInstance->projectName, extDAT);

    status = abaqus_readDATDisplacement(aimInfo, filename, abaqusInstance->feaProblem.feaMesh.numNode, &nodeIDs, &dispMatrix);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(nodeIDs, aimInfo, status);
    AIM_NOTNULL(dispMatrix, aimInfo, status);

    val.dim  = Vector;
    val.type = Double;
    val.nrow = abaqusInstance->feaProblem.feaMesh.numNode;
    val.ncol = 1;
    val.length = val.nrow*val.ncol;
    AIM_ALLOC(val.vals.reals, abaqusInstance->feaProblem.feaMesh.numNode, double, aimInfo, status);
    for (i = 0; i < abaqusInstance->feaProblem.feaMesh.numNode; i++) val.vals.reals[i] = 0;

    status = mesh_nodeID2Array(&abaqusInstance->feaProblem.feaMesh, &n2a);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(n2a, aimInfo, status);

    abaqusInstance->Tmax  = 0;
    abaqusInstance->T1max = 0;
    abaqusInstance->T2max = 0;
    abaqusInstance->T3max = 0;

    // extract maximum displacements
    for (i = 0; i < abaqusInstance->feaProblem.feaMesh.numNode; i++) {
      if (nodeIDs[i] == -1) continue;
      TT = sqrt(pow(dispMatrix[i][0], 2)
              + pow(dispMatrix[i][1], 2)
              + pow(dispMatrix[i][2], 2));

      val.vals.reals[n2a[nodeIDs[i]]] = TT;

      if (fabs(dispMatrix[i][0]) > abaqusInstance->T1max) abaqusInstance->T1max = fabs(dispMatrix[i][0]);
      if (fabs(dispMatrix[i][1]) > abaqusInstance->T2max) abaqusInstance->T2max = fabs(dispMatrix[i][1]);
      if (fabs(dispMatrix[i][2]) > abaqusInstance->T3max) abaqusInstance->T3max = fabs(dispMatrix[i][2]);
      if (TT                     > abaqusInstance->Tmax ) abaqusInstance->Tmax  = TT;
    }

    /* create the dynamic output */
    status = aim_makeDynamicOutput(aimInfo, "Displacement", &val);
    AIM_STATUS(aimInfo, status);


#ifdef WRITE_VONMISES_VTK
    status = mesh_writeVTK(aimInfo, "vonMises",
                           1,
                           &abaqusInstance->feaProblem.feaMesh,
                           1.0);
    AIM_STATUS(aimInfo, status);

    fp = aim_fopen(aimInfo, "vonMises.vtk","a");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open Output file: vonMises.vtk!");
        status = CAPS_IOERR;
        goto cleanup;
    }

    fprintf(fp,"POINT_DATA %d\n", vonMises.length);
    fprintf(fp,"SCALARS vonMises float 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (i = 0; i < vonMises.length; i++)
      fprintf(fp,"%le\n", vonMises.vals.reals[i]);

    fprintf(fp,"SCALARS disp float 3\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    AIM_ALLOC(U, abaqusInstance->feaProblem.feaMesh.numNode, DOUBLE_6, aimInfo, status);
    for (i = 0; i < abaqusInstance->feaProblem.feaMesh.numNode; i++)
         for (int j = 0; j < 6; j++) U[i][j] = 0;
    for (i = 0; i < abaqusInstance->feaProblem.feaMesh.numNode; i++) {
      if (nodeIDs[i] == -1) continue;
      for (int j = 0; j < 6; j++)
        U[n2a[nodeIDs[i]]][j] = dispMatrix[i][j];
    }
    for (i = 0; i < abaqusInstance->feaProblem.feaMesh.numNode; i++) {
      fprintf(fp,"%le %le %le\n", U[i][0], U[i][1], U[i][2]);
    }
#endif
  }

  status = CAPS_SUCCESS;

cleanup:
    if (fp != NULL) fclose(fp);

    AIM_FREE(n2a);
    AIM_FREE(nodeIDs);
    AIM_FREE(dispMatrix);
    AIM_FREE(elementIDs);
    AIM_FREE(elemData);
#ifdef WRITE_VONMISES_VTK
    AIM_FREE(U);
#endif

    return status;
}


// ********************** AIM Function Break *****************************
// Set Abaqus output variables
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{

    /*! \page aimOutputsAbaqus AIM Outputs
     * The following list outlines the Abaqus outputs available through the AIM interface.
     *
     * None to date.
     */

#ifdef DEBUG
    printf(" abaqusAIM/aimOutputs index = %d!\n", index);
#endif
    int status = CAPS_SUCCESS;

    /*
     * - <B>EigenValue</B> = List of Eigen-Values (\f$ \lambda\f$) after a modal solve.
     * - <B>EigenRadian</B> = List of Eigen-Values in terms of radians (\f$ \omega = \sqrt{\lambda}\f$ ) after a modal solve.
     * - <B>EigenFrequency</B> = List of Eigen-Values in terms of frequencies (\f$ f = \frac{\omega}{2\pi}\f$) after a modal solve.
     * - <B>EigenGeneralMass</B> = List of generalized masses for the Eigen-Values.
*/

    /*! \page aimOutputsAbaqus AIM Outputs
     * - <B>Tmax</B> = Maximum displacement.
     * - <B>T1max</B> = Maximum x-coordinate displacement.
     * - <B>T2max</B> = Maximum y-coordinate displacement.
     * - <B>T3max</B> = Maximum z-coordinate displacement.
     * - <B>vonMises_Grid</B> = Grid coordinate von Mises stress (dynamic output, only static analysis).
     * - <B>Displacement</B> = Grid coordinate displacement (dynamic output, only static analysis).
     */

//    if (index == EigenValue) {
//        *aoname = EG_strdup("EigenValue");
//
//    } else if (index == EigenRadian) {
//        *aoname = EG_strdup("EigenRadian");
//
//    } else if (index == EigenFrequency) {
//        *aoname = EG_strdup("EigenFrequency");
//
//    } else if (index == EigenGeneralMass) {
//        *aoname = EG_strdup("EigenGeneralMass");
//
//    } else
      if (index == Tmax) {
        *aoname = EG_strdup("Tmax");

    } else if (index == T1max) {
        *aoname = EG_strdup("T1max");

    } else if (index == T2max) {
        *aoname = EG_strdup("T2max");

    } else if (index == T3max) {
        *aoname = EG_strdup("T3max");

    } else {
      AIM_ERROR(aimInfo, "Unknown output index %d", index);
      status = CAPS_NOTIMPLEMENT;
    }

//    if (index >= EigenValue && index <= EigenGeneralMass) {
//      form->type   = Double;
//      form->dim    = Vector;
//      form->lfixed = Change;
//      form->sfixed = Change;
//    } else
      if (index >= Tmax && index <= T3max) {
      form->type   = Double;
      form->dim    = Scalar;
      form->nrow   = 1;
      form->ncol   = 1;
      form->lfixed = Fixed;
      form->sfixed = Fixed;
    } else {
      AIM_ERROR(aimInfo, "Unknown output index %d", index);
      status = CAPS_NOTIMPLEMENT;
    }

    return status;
}


// ********************** AIM Function Break *****************************
// Calculate Abaqus output
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, /*@unused@*/ int index,
                  /*@unused@*/ capsValue *val)
{
  int status = CAPS_SUCCESS;

  aimStorage *abaqusInstance;

  abaqusInstance = (aimStorage *) instStore;

  if (index >= Tmax && index <= T3max) {

    val->dim  = Scalar;
    val->nrow = 1;
    val->ncol = 1;
    val->length = val->nrow * val->ncol;

    if        (index == Tmax) {
      val->vals.real = abaqusInstance->Tmax;
    } else if (index == T1max) {
      val->vals.real = abaqusInstance->T1max;
    } else if (index == T2max) {
      val->vals.real = abaqusInstance->T2max;
    } else if (index == T3max) {
      val->vals.real = abaqusInstance->T3max;
    }
  }

  status = CAPS_SUCCESS;

//cleanup:

  return status;
}

void aimCleanup(void *instStore)
{
    int status; // Returning status

    aimStorage *abaqusInstance;

    #ifdef DEBUG
        printf(" abaqusAIM/aimCleanup!\n");
    #endif
    abaqusInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(abaqusInstance);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during clean up of instance\n", status);

    AIM_FREE(abaqusInstance);

}

/************************************************************************/
// CAPS transferring functions


int aimTransfer(capsDiscr *discr, const char *dataName, int numPoint,
                int dataRank, double *dataVal, /*@unused@*/ char **units)
{

    /*! \page dataTransferAbaqus Abaqus Data Transfer
     *
     * The Nastran AIM has the ability to transfer displacements from the AIM and pressure
     * distributions to the AIM using the conservative and interpolative data transfer schemes in CAPS.
     *
     * \section dataTransferAbaqus Data transfer from Abaqus
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements from the *.fil file.
     * </ul>
     *
     * \section dataTransferAbaqus Data transfer to Abaqus
     * <ul>
     *  <li> <B>"Pressure"</B> </li> <br>
     *  Writes appropriate load cards using the provided pressure distribution.
     * </ul>
     *
     */

    int status; // Function return status
    int i, j, dataPoint; // Indexing
    aimStorage *abaqusInstance;

    char *extFIL= ".fil";

    // FIL data variables
    int numGridPoint = 0;
    int numEigenVector = 0;

    double **dataMatrix = NULL;

    // Specific EigenVector to use
    int eigenVectorIndex = 0;

    // Variables used in global node mapping
    int bIndex, globalNodeID;

    // Filename stuff
    char filename[PATH_MAX];
//    FILE *fp; // File pointer

#ifdef DEBUG
    printf(" abaqusAIM/aimTransfer name = %s  instance = %d  npts = %d/%d!\n", dataName, discr->instance, numPoint, dataRank);
#endif
    abaqusInstance = (aimStorage *) discr->instStore;

    if (strcasecmp(dataName, "Displacement") != 0){// &&
        //strncmp(dataName, "EigenVector", 11) != 0) {

        printf("Unrecognized data transfer variable - %s\n", dataName);
        return CAPS_NOTFOUND;
    }

    status = aim_file(discr->aInfo, abaqusInstance->projectName, filename);
    AIM_STATUS(discr->aInfo, status);
    strcat(filename, extFIL);

    if (strcasecmp(dataName, "Displacement") == 0) {

        if (dataRank != 3) {

            AIM_ERROR(discr->aInfo, "Invalid rank for dataName \"%s\" - excepted a rank of 3!!!", dataName);
            status = CAPS_BADRANK;
            goto cleanup;

        } else {

            status = abaqus_readFIL(discr->aInfo, filename, 0, &numGridPoint,  &dataMatrix);
            AIM_STATUS(discr->aInfo, status);
        }

    } else if (strncmp(dataName, "EigenVector", 11) == 0) { // NOT SUPPOSED TO BE IMPLEMENTED YET

        // Which EigenVector do we want ?
        for (i = 0; i < strlen(dataName); i++) {
            if (dataName[i] == '_' ) break;
        }

        if (i == strlen(dataName)) {
            eigenVectorIndex = 1;
        } else {

            status = sscanf(dataName, "EigenVector_%d", &eigenVectorIndex);
            if (status != 1) {
                printf("Unable to determine which EigenVector to use - Defaulting the first EigenVector!!!\n");
                eigenVectorIndex = 1;
            }
        }

        if (dataRank != 3) {

            AIM_ERROR(discr->aInfo, "Invalid rank for dataName \"%s\" - excepted a rank of 3!!!", dataName);
            status = CAPS_BADRANK;
            goto cleanup;

        } else {

            // NEED A FUNCTION TO RETRIEVE EIGEN-VECTORS
        }



    } else {

        AIM_ERROR(discr->aInfo, "Unknown dataName \"%s\"!", dataName);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Check EigenVector range
    if (strncmp(dataName, "EigenVector", 11) == 0)  {
        if (eigenVectorIndex > numEigenVector) {
            AIM_ERROR(discr->aInfo, "Only %d EigenVectors found but index %d requested!\n", numEigenVector, eigenVectorIndex);
            status = CAPS_RANGEERR;
            goto cleanup;
        }

        if (eigenVectorIndex < 1) {
            AIM_ERROR(discr->aInfo, "For EigenVector_X notation, X must be >= 1, currently X = %d\n", eigenVectorIndex);
            status = CAPS_RANGEERR;
            goto cleanup;
        }

    }
    AIM_NOTNULL(dataMatrix, discr->aInfo, status);

    for (i = 0; i < numPoint; i++) {

        bIndex       = discr->tessGlobal[2*i  ];
        globalNodeID = discr->tessGlobal[2*i+1] +
                       discr->bodys[bIndex-1].globalOffset;

        if (strcasecmp(dataName, "Displacement") == 0) {

            for (dataPoint = 0; dataPoint < numGridPoint; dataPoint++) {
                if ((int) dataMatrix[dataPoint][0] == globalNodeID) break;
            }

            if (dataPoint == numGridPoint) {
                printf("Unable to locate global ID = %d in the data matrix\n", globalNodeID);
                status = CAPS_NOTFOUND;
                goto cleanup;
            }

            dataVal[dataRank*i+0] = dataMatrix[dataPoint][1]; // T1
            dataVal[dataRank*i+1] = dataMatrix[dataPoint][2]; // T2
            dataVal[dataRank*i+2] = dataMatrix[dataPoint][3]; // T3

        } else if (strncmp(dataName, "EigenVector", 11) == 0) { // NEED A FUNCTION TO RETRIEVE EIGEN-VECTORS


            for (dataPoint = 0; dataPoint < numGridPoint; dataPoint++) {
                if ((int) dataMatrix[eigenVectorIndex - 1][8*dataPoint + 0] ==  globalNodeID) break;
            }

            if (dataPoint == numGridPoint) {
                printf("Unable to locate global ID = %d in the data matrix\n", globalNodeID);
                status = CAPS_NOTFOUND;
                goto cleanup;
            }

            dataVal[dataRank*i+0] = dataMatrix[eigenVectorIndex- 1][8*dataPoint + 2]; // T1
            dataVal[dataRank*i+1] = dataMatrix[eigenVectorIndex- 1][8*dataPoint + 3]; // T2
            dataVal[dataRank*i+2] = dataMatrix[eigenVectorIndex- 1][8*dataPoint + 4]; // T3
            //dataVal[dataRank*i+3] = dataMatrix[eigenVectorIndex- 1][8*dataPoint + 5]; // R1 - Don't use rotations
            //dataVal[dataRank*i+4] = dataMatrix[eigenVectorIndex- 1][8*dataPoint + 6]; // R2
            //dataVal[dataRank*i+5] = dataMatrix[eigenVectorIndex- 1][8*dataPoint + 7]; // R3

        }
    }

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        // Free data matrix
        if (dataMatrix != NULL) {
            j = 0;
            if      (strcasecmp(dataName, "Displacement") == 0) j = numGridPoint;
            else if (strncmp(dataName, "EigenVector", 11) == 0) j = numEigenVector;

            for (i = 0; i < j; i++) {
                if (dataMatrix[i] != NULL) AIM_FREE(dataMatrix[i]);
            }
            AIM_FREE(dataMatrix);
        }

        return status;
}

int aimDiscr(char *tname, capsDiscr *discr) {

    int status; // Function return status
    int        i;
    ego        *tess = NULL;
    capsValue  *valMesh;
    aimStorage *abaqusInstance;

    abaqusInstance = (aimStorage *) discr->instStore;

    #ifdef DEBUG
        printf(" abaqusAIM/aimDiscr: tname = %s, instance = %d!\n", tname, iIndex);
    #endif

    if (tname == NULL) return CAPS_NOTFOUND;

    status = aim_getValue(discr->aInfo, Mesh, ANALYSISIN, &valMesh);
    AIM_STATUS(discr->aInfo, status);

    if (valMesh->nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(discr->aInfo, Mesh, "'Mesh' input must be linked to an output 'Surface_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Check and generate/retrieve the mesh
    status = checkAndCreateMesh(discr->aInfo, abaqusInstance);
    AIM_STATUS(discr->aInfo, status);

    AIM_ALLOC(tess, abaqusInstance->numMesh, ego, discr->aInfo, status);
    for (i = 0; i < abaqusInstance->numMesh; i++) {
      tess[i] = abaqusInstance->feaMesh[i].egadsTess;
    }

    status = mesh_fillDiscr(tname, &abaqusInstance->attrMap, abaqusInstance->numMesh, tess, discr);
    AIM_STATUS(discr->aInfo, status);

    #ifdef DEBUG
        printf(" abaqusAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
    #endif

    status = CAPS_SUCCESS;

cleanup:
    AIM_FREE(tess);

    return status;
}

void aimFreeDiscrPtr(void *ptr)
{
  // Extra information to store into the discr void pointer - just a int array
  EG_free(ptr);
}

int aimLocateElement(capsDiscr *discr,  double *params, double *param,
                     int *bIndex, int *eIndex, double *bary)
{
#ifdef DEBUG
    printf(" abaqusAIM/aimLocateElement instance = %d!\n", discr->instance);
#endif

    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}

int aimInterpolation(capsDiscr *discr, const char *name,
                     int bIndex, int eIndex, double *bary,
                     int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" abaqusAIM/aimInterpolation  %s  instance = %d!\n", name, discr->instance);
#endif

    return aim_interpolation(discr, name, bIndex, eIndex, bary, rank, data,
                             result);
}


int aimInterpolateBar(capsDiscr *discr, const char *name,
                      int bIndex, int eIndex, double *bary,
                      int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" abaqusAIM/aimInterpolateBar  %s  instance = %d!\n",  name, discr->instance);
#endif

    return aim_interpolateBar(discr, name, bIndex, eIndex, bary, rank, r_bar,
                              d_bar);
}


int aimIntegration(capsDiscr *discr, const char *name,int bIndex, int eIndex,
                   int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" abaqusAIM/aimIntegration  %s  instance = %d!\n", name, discr->instance);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                    int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" abaqusAIM/aimIntegrateBar  %s  instance = %d!\n", name, discr->instance);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank, r_bar, d_bar);
}
