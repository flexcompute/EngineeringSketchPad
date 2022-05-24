/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             TACS AIM
 *
 *     Written by Dr. Marshall C. Galbraith MIT, and Dr. Ryan Durscher and Dr. Ed Alyanak
 *
 */

/*! \mainpage Introduction
 *
 * \section overviewTACS TACS AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with the finite element structural solver TACS \cite TACS.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsTACS and
 * \ref aimOutputsTACS and \ref attributeTACS, respectively.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferTACS
 */

/*! \page attributeTACS AIM attributes
 * The following list of attributes are required for the TACS AIM inside the geometry input.
 *
 * - <b> capsDiscipline</b> This attribute is a requirement if doing aeroelastic analysis within TACS. capsDiscipline allows
 * the AIM to determine which bodies are meant for structural analysis and which are used for aerodynamics. Options
 * are:  Structure and Aerodynamic (case insensitive).
 *
 * - <b> capsGroup</b> This is a name assigned to any geometric body to denote a property.  This body could be a solid, surface, face, wire, edge or node.
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
 * - <b> capsIgnore</b> It is possible that there is a geometric body (or entity) that you do not want the TACS AIM to pay attention to when creating
 * a finite element model. The capsIgnore attribute allows a body (or entity) to be in the geometry and ignored by the AIM.  For example,
 * because of limitations in OpenCASCADE a situation where two edges are overlapping may occur; capsIgnore allows the user to only
 *  pay attention to one of the overlapping edges.
 *
 * - <b> capsConnect</b> This is a name assigned to any geometric body where the user wishes to create
 * "fictitious" connections such as springs, dampers, and/or rigid body connections to. The user must manually specify
 * the connection between two <c>capsConnect</c> entities using the "Connect" tuple (see \ref aimInputsTACS).
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsConnect $springStart</c>.
 *
 * - <b> capsConnectLink</b> Similar to <c>capsConnect</c>, this is a name assigned to any geometric body where
 * the user wishes to create "fictitious" connections to. A connection is automatically made if a <c>capsConnectLink</c>
 * matches a <c>capsConnect</c> group. Again further specifics of the connection are input using the "Connect"
 * tuple (see \ref aimInputsTACS). Recall that a string in ESP starts with a $.
 * For example, attribute <c>capsConnectLink $springEnd</c>.
 *
 * - <b> capsResponse</b> This is a name assigned to any geometric body that will be used to define design sensitivity
 * responses for optimization. Specific information for the responses are input using the "Design_Response" tuple (see
 * \ref aimInputsTACS). Recall that a string in ESP starts with a $. For examples,
 * attribute <c>capsResponse $displacementNode</c>.
 *
 * - <b> capsBound </b> This is used to mark surfaces on the structural grid in which data transfer with an external
 * solver will take place. See \ref dataTransferTACS for additional details.
 *
 * Internal Aeroelastic Analysis
 *
 * - <b> capsBound </b> This is used to mark surfaces on the structural grid in which a spline will be created between
 * the structural and aero-loads.
 *
 * - <b> capsReferenceArea</b>  [Optional: Default 1.0] Reference area to use when doing aeroelastic analysis.
 * This attribute may exist on any aerodynamic cross-section.
 *
 * - <b> capsReferenceChord</b>  [Optional: Default 1.0] Reference chord to use when doing aeroelastic analysis.
 * This attribute may exist on any aerodynamic cross-section.
 *
 * - <b> capsReferenceSpan</b>  [Optional: Default 1.0] Reference span to use when doing aeroelastic analysis.
 * This attribute may exist on any aerodynamic cross-section.
 *
 */

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"    // Meshing utilities
#include "miscUtils.h"    // Miscellaneous utilities
#include "vlmUtils.h"     // Vortex lattice method utilities
#include "vlmSpanSpace.h"
#include "feaUtils.h"     // FEA utilities
#include "nastranUtils.h" // Nastran utilities

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

#define MXCHAR     255

//#define DEBUG

enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Property,
  Material,
  Constraint,
  Load,
  Analysix,
  Analysis_Type,
  File_Format,
  Mesh_File_Format,
  inDesign_Variable,
  Design_Variable_Relation,
  Design_Constraint,
  Design_Equation,
  Design_Table,
  Design_Response,
  Design_Equation_Response,
  Design_Opt_Param,
  Support,
  Connect,
  Parameter,
  Mesh,
  NUMINPUT = Mesh              /* Total number of inputs */
};

#define NUMOUTPUT  0


typedef struct {

    // Project name
    char *projectName; // Project name

    feaUnitsStruct units; // units system

    feaProblemStruct feaProblem;

    // Attribute to capsGroup index map
    mapAttrToIndexStruct groupMap;

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

} aimStorage;



static int initiate_aimStorage(aimStorage *tacsInstance)
{

    int status;

    // Set initial values for tacsInstance
    tacsInstance->projectName = NULL;

/*
    // Check to make sure data transfer is ok
    tacsInstance->dataTransferCheck = (int) true;
*/

    status = initiate_feaUnitsStruct(&tacsInstance->units);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&tacsInstance->groupMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to constraint index map
    status = initiate_mapAttrToIndexStruct(&tacsInstance->constraintMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to load index map
    status = initiate_mapAttrToIndexStruct(&tacsInstance->loadMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for transfer to index map
    status = initiate_mapAttrToIndexStruct(&tacsInstance->transferMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for connect to index map
    status = initiate_mapAttrToIndexStruct(&tacsInstance->connectMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for response to index map
    status = initiate_mapAttrToIndexStruct(&tacsInstance->responseMap);
    if (status != CAPS_SUCCESS) return status;

    status = initiate_feaProblemStruct(&tacsInstance->feaProblem);
    if (status != CAPS_SUCCESS) return status;

    // Mesh holders
    tacsInstance->numMesh = 0;
    tacsInstance->feaMesh = NULL;

    return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *tacsInstance)
{
    int i, status;

    status = destroy_feaUnitsStruct(&tacsInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&tacsInstance->groupMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to constraint index map
    status = destroy_mapAttrToIndexStruct(&tacsInstance->constraintMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to load index map
    status = destroy_mapAttrToIndexStruct(&tacsInstance->loadMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Transfer to index map
    status = destroy_mapAttrToIndexStruct(&tacsInstance->transferMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Connect to index map
    status = destroy_mapAttrToIndexStruct(&tacsInstance->connectMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Response to index map
    status = destroy_mapAttrToIndexStruct(&tacsInstance->responseMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (tacsInstance->feaMesh != NULL) {

        for (i = 0; i < tacsInstance->numMesh; i++) {
            status = destroy_meshStruct(&tacsInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        EG_free(tacsInstance->feaMesh);
    }

    tacsInstance->feaMesh = NULL;
    tacsInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&tacsInstance->feaProblem);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    // NULL projetName
    tacsInstance->projectName = NULL;

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(void *aimInfo, aimStorage *tacsInstance)
{
  // Function return flag
  int status = CAPS_SUCCESS;

  status = fea_createMesh(aimInfo,
                          NULL,
                          0,
                          0,
                          (int)false,
                          &tacsInstance->groupMap,
                          &tacsInstance->constraintMap,
                          &tacsInstance->loadMap,
                          &tacsInstance->transferMap,
                          &tacsInstance->connectMap,
                          &tacsInstance->responseMap,
                          &tacsInstance->numMesh,
                          &tacsInstance->feaMesh,
                          &tacsInstance->feaProblem );
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
    int  *ints=NULL, i, status = CAPS_SUCCESS;
    char **strs=NULL;

    aimStorage *tacsInstance=NULL;

#ifdef DEBUG
    printf("nastranAIM/aimInitialize   instance = %d!\n", inst);
#endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;

    /* specify the name of each field variable */
//    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);

//    strs[0]  = EG_strdup("Displacement");
//    strs[1]  = EG_strdup("EigenVector");
//    strs[2]  = EG_strdup("EigenVector_#");
//    strs[3]  = EG_strdup("Pressure");
//    for (i = 0; i < *nFields; i++)
//      if (strs[i] == NULL) { status = EGADS_MALLOC; goto cleanup; }
    *fnames  = strs;

    /* specify the dimension of each field variable */
//    AIM_ALLOC(ints, *nFields, int, aimInfo, status);
//    ints[0]  = 3;
//    ints[1]  = 3;
//    ints[2]  = 3;
//    ints[3]  = 1;
    *franks  = ints;
    ints = NULL;

    /* specify if a field is an input field or output field */
//    AIM_ALLOC(ints, *nFields, int, aimInfo, status);
//
//    ints[0]  = FieldOut;
//    ints[1]  = FieldOut;
//    ints[2]  = FieldOut;
//    ints[3]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;

    // Allocate tacsInstance
    AIM_ALLOC(tacsInstance, 1, aimStorage, aimInfo, status);
    *instStore = tacsInstance;

    // Initialize instance storage
    (void) initiate_aimStorage(tacsInstance);

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

    /*! \page aimInputsTACS AIM Inputs
     * The following list outlines the TACS inputs along with their default value available
     * through the AIM interface. Unless noted these values will be not be linked to
     * any parent AIMs with variables of the same name.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" nastranAIM/aimInputs index = %d!\n", index);
#endif

    *ainame = NULL;

    // TACS Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("nastran_CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsTACS
         * - <B> Proj_Name = "nastran_CAPS"</B> <br>
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

        /*! \page aimInputsTACS
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

        /*! \page aimInputsTACS
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

        /*! \page aimInputsTACS
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

        /*! \page aimInputsTACS
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

        /*! \page aimInputsTACS
         * - <B> Analysis = NULL</B> <br>
         * Analysis tuple used to input analysis/case information for the model, see \ref feaAnalysis for additional details.
         */
    } else if (index == Analysis_Type) {
        *ainame              = EG_strdup("Analysis_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Modal");
        defval->lfixed       = Change;

        /*! \page aimInputsTACS
         * - <B> Analysis_Type = "Modal"</B> <br>
         * Type of analysis to generate files for, options include "Modal", "Static", "AeroelasticTrim", "AeroelasticFlutter", and "Optimization".
         * Note: "Aeroelastic" and "StaticOpt" are still supported and refer to "AeroelasticTrim" and "Optimization".
         */
    } else if (index == File_Format) {
        *ainame              = EG_strdup("File_Format");
        defval->type         = String;
        defval->vals.string  = EG_strdup("Small"); // Small, Large, Free
        defval->lfixed       = Change;

        /*! \page aimInputsTACS
         * - <B> File_Format = "Small"</B> <br>
         * Formatting type for the bulk file. Options: "Small", "Large", "Free".
         */

    } else if (index == Mesh_File_Format) {
        *ainame              = EG_strdup("Mesh_File_Format");
        defval->type         = String;
        defval->vals.string  = EG_strdup("Free"); // Small, Large, Free
        defval->lfixed       = Change;

        /*! \page aimInputsTACS
         * - <B> Mesh_File_Format = "Small"</B> <br>
         * Formatting type for the mesh file. Options: "Small", "Large", "Free".
         */

    } else if (index == inDesign_Variable) {
        *ainame              = EG_strdup("Design_Variable");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Variable = NULL</B> <br>
         * The design variable tuple is used to input design variable information for the model optimization, see \ref feaDesignVariable for additional details.
         */

    } else if (index == Design_Variable_Relation) {
        *ainame              = EG_strdup("Design_Variable_Relation");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Variable_Relation = NULL</B> <br>
         * The design variable relation tuple is used to input design variable relation information for the model optimization, see \ref feaDesignVariableRelation for additional details.
         */

    } else if (index == Design_Constraint) {
        *ainame              = EG_strdup("Design_Constraint");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Constraint = NULL</B> <br>
         * The design constraint tuple is used to input design constraint information for the model optimization, see \ref feaDesignConstraint for additional details.
         */

    } else if (index == Design_Equation) {
        *ainame              = EG_strdup("Design_Equation");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Equation = NULL</B> <br>
         * The design equation tuple used to input information defining equations for use in design sensitivity, see \ref feaDesignEquation for additional details.
         */

    } else if (index == Design_Table) {
        *ainame              = EG_strdup("Design_Table");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Table = NULL</B> <br>
         * The design table tuple used to input table of real constants used in equations, see \ref feaDesignTable for additional details.
         */

    } else if (index == Design_Response) {
        *ainame              = EG_strdup("Design_Response");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Response = NULL</B> <br>
         * The design response tuple used to input design sensitivity response information, see \ref feaDesignResponse for additional details.
         */

    } else if (index == Design_Equation_Response) {
        *ainame              = EG_strdup("Design_Equation_Response");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Equation_Response = NULL</B> <br>
         * The design equation response tuple used to input design sensitivity equation response information, see \ref feaDesignEquationResponse for additional details.
         */

    } else if (index == Design_Opt_Param) {
        *ainame              = EG_strdup("Design_Opt_Param");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Design_Opt_Param = NULL</B> <br>
         * The design optimization parameter tuple used to input parameters used in design optimization.
         */

    } else if (index == Support) {
        *ainame              = EG_strdup("Support");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
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

        /*! \page aimInputsTACS
         * - <B> Connect = NULL</B> <br>
         * Connect tuple used to define connection to be made in the, see \ref feaConnection for additional details.
         */
    } else if (index == Parameter) {
        *ainame              = EG_strdup("Parameter");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsTACS
         * - <B> Parameter = NULL</B> <br>
         * Parameter tuple used to define PARAM entries. Note, entries are output exactly as inputed, that is, if the PARAM entry
         * requires an integer entry the user must input an integer!
         */

    } else if (index == Mesh) {
        *ainame             = AIM_NAME(Mesh);
        defval->type        = Pointer;
        defval->dim         = Vector;
        defval->lfixed      = Change;
        defval->sfixed      = Change;
        defval->vals.AIMptr = NULL;
        defval->nullVal     = IsNull;
        AIM_STRDUP(defval->units, "meshStruct", aimInfo, status);

        /*! \page aimInputsTACS
         * - <B>Mesh = NULL</B> <br>
         * A Mesh link.
         */
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}


int aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{

    int i, j, k, l; // Indexing

    int status; // Status return

    int found; // Boolean operator

    int *tempIntegerArray = NULL; // Temporary array to store a list of integers

    // Analysis information
    char *analysisType = NULL;
    int haveSubAeroelasticTrim = (int) false;
    int haveSubAeroelasticFlutter = (int) false;

    // Aeroelastic Information
    int numAEStatSurf = 0;
    //char **aeStatSurf = NULL;

    // File format information
    char *tempString = NULL, *delimiter = NULL;

    // File IO
    char *filename = NULL; // Output file name
    FILE *fp = NULL; // Output file pointer

    feaLoadStruct *feaLoad;
    int numThermalLoad=0;

    int numSetID;
    int tempID, *setID = NULL;

    aimStorage *tacsInstance;

    tacsInstance = (aimStorage *) instStore;
    if (aimInputs == NULL) return CAPS_NULLVALUE;

    // Get project name
    tacsInstance->projectName = aimInputs[Proj_Name-1].vals.string;

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    if (aimInputs[Mesh-1].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Mesh, "'Mesh' input must be linked to an output 'Surface_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get FEA mesh if we don't already have one
    status = checkAndCreateMesh(aimInfo, tacsInstance);
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
                                 &tacsInstance->units,
                                 &tacsInstance->feaProblem.numMaterial,
                                 &tacsInstance->feaProblem.feaMaterial);
        AIM_STATUS(aimInfo, status);
    } else printf("Material tuple is NULL - No materials set\n");

    // Set property properties
    if (aimInputs[Property-1].nullVal == NotNull) {
        status = fea_getProperty(aimInfo,
                                 aimInputs[Property-1].length,
                                 aimInputs[Property-1].vals.tuple,
                                 &tacsInstance->groupMap,
                                 &tacsInstance->units,
                                 &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);


        // Assign element "subtypes" based on properties set
        status = fea_assignElementSubType(tacsInstance->feaProblem.numProperty,
                                          tacsInstance->feaProblem.feaProperty,
                                          &tacsInstance->feaProblem.feaMesh);
        AIM_STATUS(aimInfo, status);

    } else printf("Property tuple is NULL - No properties set\n");

    // Set constraint properties
    if (aimInputs[Constraint-1].nullVal == NotNull) {
        status = fea_getConstraint(aimInputs[Constraint-1].length,
                                   aimInputs[Constraint-1].vals.tuple,
                                   &tacsInstance->constraintMap,
                                   &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Constraint tuple is NULL - No constraints applied\n");

    // Set support properties
    if (aimInputs[Support-1].nullVal == NotNull) {
        status = fea_getSupport(aimInputs[Support-1].length,
                                aimInputs[Support-1].vals.tuple,
                                &tacsInstance->constraintMap,
                                &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Support tuple is NULL - No supports applied\n");

    // Set connection properties
    if (aimInputs[Connect-1].nullVal == NotNull) {
        status = fea_getConnection(aimInputs[Connect-1].length,
                                   aimInputs[Connect-1].vals.tuple,
                                   &tacsInstance->connectMap,
                                   &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Connect tuple is NULL - Using defaults\n");

    // Set load properties
    if (aimInputs[Load-1].nullVal == NotNull) {
        status = fea_getLoad(aimInputs[Load-1].length,
                             aimInputs[Load-1].vals.tuple,
                             &tacsInstance->loadMap,
                             &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);

        // Loop through loads to see if any of them are supposed to be from an external source
        for (i = 0; i < tacsInstance->feaProblem.numLoad; i++) {

            if (tacsInstance->feaProblem.feaLoad[i].loadType == PressureExternal) {

                // Transfer external pressures from the AIM discrObj
                status = fea_transferExternalPressure(aimInfo,
                                                      &tacsInstance->feaProblem.feaMesh,
                                                      &tacsInstance->feaProblem.feaLoad[i]);
                AIM_STATUS(aimInfo, status);

            } // End PressureExternal if
        } // End load for loop
    } else printf("Load tuple is NULL - No loads applied\n");

    // Set design variables
    if (aimInputs[inDesign_Variable-1].nullVal == NotNull) {
        status = fea_getDesignVariable(aimInfo,
                                       (int)false,
                                       aimInputs[inDesign_Variable-1].length,
                                       aimInputs[inDesign_Variable-1].vals.tuple,
                                       aimInputs[Design_Variable_Relation-1].length,
                                       aimInputs[Design_Variable_Relation-1].vals.tuple,
                                       &tacsInstance->groupMap,
                                       &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);

    } else printf("Design_Variable tuple is NULL - No design variables applied\n");

    // Set design constraints
    if (aimInputs[Design_Constraint-1].nullVal == NotNull) {
        status = fea_getDesignConstraint(aimInputs[Design_Constraint-1].length,
                                         aimInputs[Design_Constraint-1].vals.tuple,
                                         &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Constraint tuple is NULL - No design constraints applied\n");

    // Set design equations
    if (aimInputs[Design_Equation-1].nullVal == NotNull) {
        status = fea_getDesignEquation(aimInputs[Design_Equation-1].length,
                                       aimInputs[Design_Equation-1].vals.tuple,
                                       &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Equation tuple is NULL - No design equations applied\n");

    // Set design table constants
    if (aimInputs[Design_Table-1].nullVal == NotNull) {
        status = fea_getDesignTable(aimInputs[Design_Table-1].length,
                                    aimInputs[Design_Table-1].vals.tuple,
                                    &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Table tuple is NULL - No design table constants applied\n");

    // Set design optimization parameters
    if (aimInputs[Design_Opt_Param-1].nullVal == NotNull) {
        status = fea_getDesignOptParam(aimInputs[Design_Opt_Param-1].length,
                                       aimInputs[Design_Opt_Param-1].vals.tuple,
                                       &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Opt_Param tuple is NULL - No design optimization parameters applied\n");

    // Set design responses
    if (aimInputs[Design_Response-1].nullVal == NotNull) {
        status = fea_getDesignResponse(aimInputs[Design_Response-1].length,
                                       aimInputs[Design_Response-1].vals.tuple,
                                       &tacsInstance->responseMap,
                                       &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Response tuple is NULL - No design responses applied\n");

    // Set design equation responses
    if (aimInputs[Design_Equation_Response-1].nullVal == NotNull) {
        status = fea_getDesignEquationResponse(aimInputs[Design_Equation_Response-1].length,
                                               aimInputs[Design_Equation_Response-1].vals.tuple,
                                               &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Equation_Response tuple is NULL - No design equation responses applied\n");

    // Set analysis settings
    if (aimInputs[Analysix-1].nullVal == NotNull) {
        status = fea_getAnalysis(aimInputs[Analysix-1].length,
                                 aimInputs[Analysix-1].vals.tuple,
                                 &tacsInstance->feaProblem);
        AIM_STATUS(aimInfo, status); // It ok to not have an analysis tuple
    } else {
        printf("Analysis tuple is NULL\n");
        status = fea_createDefaultAnalysis(&tacsInstance->feaProblem, analysisType);
        AIM_STATUS(aimInfo, status);
    }


    // Set file format type
    if        (strcasecmp(aimInputs[File_Format-1].vals.string, "Small") == 0) {
    	tacsInstance->feaProblem.feaFileFormat.fileType = SmallField;
    } else if (strcasecmp(aimInputs[File_Format-1].vals.string, "Large") == 0) {
    	tacsInstance->feaProblem.feaFileFormat.fileType = LargeField;
    } else if (strcasecmp(aimInputs[File_Format-1].vals.string, "Free") == 0)  {
    	tacsInstance->feaProblem.feaFileFormat.fileType = FreeField;
    } else {
    	printf("Unrecognized \"File_Format\", valid choices are [Small, Large, or Free]. Reverting to default\n");
    }

    // Set grid file format type
    if        (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Small") == 0) {
        tacsInstance->feaProblem.feaFileFormat.gridFileType = SmallField;
    } else if (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Large") == 0) {
        tacsInstance->feaProblem.feaFileFormat.gridFileType = LargeField;
    } else if (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Free") == 0)  {
        tacsInstance->feaProblem.feaFileFormat.gridFileType = FreeField;
    } else {
        printf("Unrecognized \"Mesh_File_Format\", valid choices are [Small, Large, or Free]. Reverting to default\n");
    }

    // Write TACS Mesh
    filename = EG_alloc(MXCHAR +1);
    if (filename == NULL) return EGADS_MALLOC;

    strcpy(filename, tacsInstance->projectName);

    status = mesh_writeNASTRAN(aimInfo,
                               filename,
                               1,
                               &tacsInstance->feaProblem.feaMesh,
                               tacsInstance->feaProblem.feaFileFormat.gridFileType,
                               1.0);
    if (status != CAPS_SUCCESS) {
        EG_free(filename);
        return status;
    }

    // Write TACS subElement types not supported by mesh_writeNASTRAN
    strcat(filename, ".bdf");
    fp = aim_fopen(aimInfo, filename, "a");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        EG_free(filename);
        return CAPS_IOERR;
    }
    EG_free(filename);

    printf("Writing subElement types (if any) - appending mesh file\n");
    status = nastran_writeSubElementCard(fp,
                                         &tacsInstance->feaProblem.feaMesh,
                                         tacsInstance->feaProblem.numProperty,
                                         tacsInstance->feaProblem.feaProperty,
                                         &tacsInstance->feaProblem.feaFileFormat);
    AIM_STATUS(aimInfo, status);

    // Connections
    for (i = 0; i < tacsInstance->feaProblem.numConnect; i++) {

        if (i == 0) {
            printf("Writing connection cards - appending mesh file\n");
        }

        status = nastran_writeConnectionCard(fp,
                                             &tacsInstance->feaProblem.feaConnect[i],
                                             &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }
    if (fp != NULL) fclose(fp);
    fp = NULL;

    // Write nastran input file
    filename = EG_alloc(MXCHAR +1);
    if (filename == NULL) return EGADS_MALLOC;
    strcpy(filename, tacsInstance->projectName);
    strcat(filename, ".dat");

    printf("\nWriting TACS instruction file....\n");
    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file: %s\n", filename);
        EG_free(filename);
        return CAPS_IOERR;
    }
    EG_free(filename);

    // define file format delimiter type
    if (tacsInstance->feaProblem.feaFileFormat.fileType == FreeField) {
        delimiter = ",";
    } else {
        delimiter = " ";
    }

    //////////////// Executive control ////////////////
    fprintf(fp, "ID CAPS generated Problem FOR TACS\n");

    // Analysis type
    if     (strcasecmp(analysisType, "Modal")         		== 0) fprintf(fp, "SOL 3\n");
    else if(strcasecmp(analysisType, "Static")        		== 0) fprintf(fp, "SOL 1\n");
    else if(strcasecmp(analysisType, "Craig-Bampton") 		== 0) fprintf(fp, "SOL 31\n");
    else if(strcasecmp(analysisType, "StaticOpt")     		== 0) fprintf(fp, "SOL 200\n");
    else if(strcasecmp(analysisType, "Optimization")  		== 0) fprintf(fp, "SOL 200\n");
    else if(strcasecmp(analysisType, "Aeroelastic")   		== 0) fprintf(fp, "SOL 144\n");
    else if(strcasecmp(analysisType, "AeroelasticTrim") 	== 0) fprintf(fp, "SOL 144\n");
    else if(strcasecmp(analysisType, "AeroelasticFlutter") 	== 0) fprintf(fp, "SOL 145\n");
    else {
        AIM_ERROR(aimInfo, "Unrecognized \"Analysis_Type\", %s", analysisType);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    fprintf(fp, "CEND\n\n");

    if (tacsInstance->feaProblem.feaMesh.numNode> 10000) fprintf(fp, "LINE=%d\n", tacsInstance->feaProblem.feaMesh.numNode*10);
    else                                                           fprintf(fp, "LINE=10000\n");


    // Set up the case information
     if (tacsInstance->feaProblem.numAnalysis == 0) {
         printf("Error: No analyses in the feaProblem! (this shouldn't be possible)\n");
         status = CAPS_BADVALUE;
         goto cleanup;
     }

    //////////////// Case control ////////////////

    // Write output request information
    fprintf(fp, "DISP (PRINT,PUNCH) = ALL\n"); // Output all displacements

    fprintf(fp, "STRE (PRINT,PUNCH) = ALL\n"); // Output all stress

    fprintf(fp, "STRA (PRINT,PUNCH) = ALL\n"); // Output all strain

    // Write sub-case information if multiple analysis tuples were provide - will always have at least 1
    for (i = 0; i < tacsInstance->feaProblem.numAnalysis; i++) {
        //printf("SUBCASE = %d\n", i);

        fprintf(fp, "SUBCASE %d\n", i+1);
        fprintf(fp, "\tLABEL = %s\n", tacsInstance->feaProblem.feaAnalysis[i].name);

        if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == Static) {
            fprintf(fp,"\tANALYSIS = STATICS\n");
        } else if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == Modal) {
            fprintf(fp,"\tANALYSIS = MODES\n");
        } else if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim) {
            fprintf(fp,"\tANALYSIS = SAERO\n");
            haveSubAeroelasticTrim = (int) true;
        } else if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {
            fprintf(fp,"\tANALYSIS = FLUTTER\n");
            haveSubAeroelasticFlutter = (int) true;
        } else if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == Optimization) {
            printf("\t *** WARNING :: INPUT TO ANALYSIS CASE INPUT analysisType should NOT be Optimization or StaticOpt - Defaulting to Static\n");
            fprintf(fp,"\tANALYSIS = STATICS\n");
            }

        // Write support for sub-case
        if (tacsInstance->feaProblem.feaAnalysis[i].numSupport != 0) {
            if (tacsInstance->feaProblem.feaAnalysis[i].numSupport > 1) {
                printf("\tWARNING: More than 1 support is not supported at this time for sub-cases!\n");

            } else {
                fprintf(fp, "\tSUPORT1 = %d\n", tacsInstance->feaProblem.feaAnalysis[i].supportSetID[0]);
            }
        }

        // Write constraint for sub-case
        if (tacsInstance->feaProblem.numConstraint != 0) {
            fprintf(fp, "\tSPC = %d\n", tacsInstance->feaProblem.numConstraint+i+1); //TODO - change to i+1 to just i
        }

        // Issue some warnings regarding constraints if necessary
        if (tacsInstance->feaProblem.feaAnalysis[i].numConstraint == 0 && tacsInstance->feaProblem.numConstraint !=0) {

            printf("\tWarning: No constraints specified for case %s, assuming "
                    "all constraints are applied!!!!\n", tacsInstance->feaProblem.feaAnalysis[i].name);

        } else if (tacsInstance->feaProblem.numConstraint == 0) {

            printf("\tWarning: No constraints specified for case %s!!!!\n", tacsInstance->feaProblem.feaAnalysis[i].name);
        }

  //        // Write MPC for sub-case - currently only supported when we have RBE2 elements - see above for unification - TODO - investigate
  //        for (j = 0; j < astrosInstance[iIndex].feaProblem.numConnect; j++) {
  //
  //            if (astrosInstance[iIndex].feaProblem.feaConnect[j].connectionType == RigidBody) {
  //
  //                if (addComma == (int) true) fprintf(fp,",");
  //
  //                fprintf(fp, " MPC = %d ", astrosInstance[iIndex].feaProblem.feaConnect[j].connectionID);
  //                addComma = (int) true;
  //                break;
  //            }
  //        }

        if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == Modal) {
            fprintf(fp,"\tMETHOD = %d\n", tacsInstance->feaProblem.feaAnalysis[i].analysisID);
        }

        if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {
            fprintf(fp,"\tMETHOD = %d\n", tacsInstance->feaProblem.feaAnalysis[i].analysisID);
            fprintf(fp,"\tFMETHOD = %d\n", 100+tacsInstance->feaProblem.feaAnalysis[i].analysisID);
        }

        if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim) {
            fprintf(fp,"\tTRIM = %d\n", tacsInstance->feaProblem.feaAnalysis[i].analysisID);
        }

        if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim ||
            tacsInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {

            if (tacsInstance->feaProblem.feaAnalysis[i].aeroSymmetryXY != NULL) {

                if(strcmp("SYM",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"\tAESYMXY = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYM",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"\tAESYMXY = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYM",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"\tAESYMXY = %s\n","ASYMMETRIC");
                } else if(strcmp("SYMMETRIC",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"\tAESYMXY = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYMMETRIC",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"\tAESYMXY = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYMMETRIC",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"\tAESYMXY = %s\n","ASYMMETRIC");
                } else {
                    printf("\t*** Warning *** aeroSymmetryXY Input %s to nastranAIM not understood!\n",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXY );
                }
            }

            if (tacsInstance->feaProblem.feaAnalysis[i].aeroSymmetryXZ != NULL) {

                if(strcmp("SYM",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"\tAESYMXZ = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYM",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"\tAESYMXZ = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYM",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"\tAESYMXZ = %s\n","ASYMMETRIC");
                } else if(strcmp("SYMMETRIC",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"\tAESYMXZ = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYMMETRIC",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"\tAESYMXZ = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYMMETRIC",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"\tAESYMXZ = %s\n","ASYMMETRIC");
                } else {
                    printf("\t*** Warning *** aeroSymmetryXZ Input %s to nastranAIM not understood!\n",tacsInstance->feaProblem.feaAnalysis->aeroSymmetryXZ );
                }

            }
        }

        // Issue some warnings regarding loads if necessary
        if (tacsInstance->feaProblem.feaAnalysis[i].numLoad == 0 && tacsInstance->feaProblem.numLoad !=0) {
            printf("\tWarning: No loads specified for case %s, assuming "
                    "all loads are applied!!!!\n", tacsInstance->feaProblem.feaAnalysis[i].name);
        } else if (tacsInstance->feaProblem.numLoad == 0) {
            printf("\tWarning: No loads specified for case %s!!!!\n", tacsInstance->feaProblem.feaAnalysis[i].name);
        }

        // Write loads for sub-case
        if (tacsInstance->feaProblem.numLoad != 0) {

            found = (int) false;

            for (k = 0; k < tacsInstance->feaProblem.numLoad; k++) {

                feaLoad = &tacsInstance->feaProblem.feaLoad[k];

                if (tacsInstance->feaProblem.feaAnalysis[i].numLoad != 0) { // if loads specified in analysis

                    for (j = 0; j < tacsInstance->feaProblem.feaAnalysis[i].numLoad; j++) { // See if the load is in the loadSet

                        if (feaLoad->loadID == tacsInstance->feaProblem.feaAnalysis[i].loadSetID[j]) break;
                    }

                    if (j >= tacsInstance->feaProblem.feaAnalysis[i].numLoad) continue; // If it isn't in the loadSet move on
                } else {
                    //pass
                }

                if (feaLoad->loadType == Thermal && numThermalLoad == 0) {

                    fprintf(fp, "\tTemperature = %d\n", feaLoad->loadID);
                    numThermalLoad += 1;
                    if (numThermalLoad > 1) {
                        printf("More than 1 Thermal load found - nastranAIM does NOT currently doesn't support multiple thermal loads in a given case!\n");
                    }

                    continue;
                }

                found = (int) true;
            }

            if (found == (int) true) {
                fprintf(fp, "\tLOAD = %d\n", tacsInstance->feaProblem.numLoad+i+1);
            }
        }

        //if (tacsInstance->feaProblem.feaAnalysis[i].analysisType == Optimization) {
            // Write objective function
            //fprintf(fp, "\tDESOBJ(%s) = %d\n", tacsInstance->feaProblem.feaAnalysis[i].objectiveMinMax,
            //                                   tacsInstance->feaProblem.feaAnalysis[i].analysisID);
            // Write optimization constraints
        if (tacsInstance->feaProblem.feaAnalysis[i].numDesignConstraint != 0) {
            fprintf(fp, "\tDESSUB = %d\n", tacsInstance->feaProblem.numDesignConstraint+i+1);
        }
        //}

        // Write response spanning for sub-case
        if (tacsInstance->feaProblem.feaAnalysis[i].numDesignResponse != 0) {

            numSetID = tacsInstance->feaProblem.feaAnalysis[i].numDesignResponse;
            // setID = tacsInstance->feaProblem.feaAnalysis[i].designResponseSetID;
            setID = EG_alloc(numSetID * sizeof(int));
            if (setID == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (j = 0; j < numSetID; j++) {
                tempID = tacsInstance->feaProblem.feaAnalysis[i].designResponseSetID[j];
                setID[j] = tempID + 100000;
            }

            tempID = i+1;
            status = nastran_writeSetCard(fp, tempID, numSetID, setID);

            EG_free(setID);

            AIM_STATUS(aimInfo, status);
            fprintf(fp, "\tDRSPAN = %d\n", tempID);
        }

    }

    //////////////// Bulk data ////////////////
    fprintf(fp, "\nBEGIN BULK\n");
    fprintf(fp, "$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|\n");
    //PRINT PARAM ENTRIES IN BULK DATA

    if (aimInputs[Parameter-1].nullVal == NotNull) {
        for (i = 0; i < aimInputs[Parameter-1].length; i++) {
            fprintf(fp, "PARAM, %s, %s\n", aimInputs[Parameter-1].vals.tuple[i].name,
                                           aimInputs[Parameter-1].vals.tuple[i].value);
        }
    }

    fprintf(fp, "PARAM, %s\n", "POST, -1\n"); // Output OP2 file

    // Write AEROS, AESTAT and AESURF cards
    if (strcasecmp(analysisType, "AeroelasticFlutter") == 0 ||
        haveSubAeroelasticFlutter == (int) true) {

        printf("\tWriting aero card\n");
        status = nastran_writeAEROCard(fp,
                                       &tacsInstance->feaProblem.feaAeroRef,
                                       &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Write AEROS, AESTAT and AESURF cards
    if (strcasecmp(analysisType, "Aeroelastic") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
        haveSubAeroelasticTrim == (int) true) {

        printf("\tWriting aeros card\n");
        status = nastran_writeAEROSCard(fp,
                                        &tacsInstance->feaProblem.feaAeroRef,
                                        &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);

        numAEStatSurf = 0;
        for (i = 0; i < tacsInstance->feaProblem.numAnalysis; i++) {

            if (tacsInstance->feaProblem.feaAnalysis[i].analysisType != AeroelasticTrim) continue;

            if (i == 0) printf("\tWriting aestat cards\n");

            // Loop over rigid variables
            for (j = 0; j < tacsInstance->feaProblem.feaAnalysis[i].numRigidVariable; j++) {

                found = (int) false;

                // Loop over previous rigid variables
                for (k = 0; k < i; k++) {
                    for (l = 0; l < tacsInstance->feaProblem.feaAnalysis[k].numRigidVariable; l++) {

                        // If current rigid variable was previous written - mark as found
                        if (strcmp(tacsInstance->feaProblem.feaAnalysis[i].rigidVariable[j],
                                   tacsInstance->feaProblem.feaAnalysis[k].rigidVariable[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // If not write out an aestat card
                numAEStatSurf += 1;

                fprintf(fp,"%-8s", "AESTAT");

                tempString = convert_integerToString(numAEStatSurf, 7, 1);
                AIM_NOTNULL(tempString, aimInfo, status);
                fprintf(fp, "%s%s", delimiter, tempString);
                AIM_FREE(tempString);

                fprintf(fp, "%s%7s\n", delimiter, tacsInstance->feaProblem.feaAnalysis[i].rigidVariable[j]);
            }

            // Loop over rigid Constraints
            for (j = 0; j < tacsInstance->feaProblem.feaAnalysis[i].numRigidConstraint; j++) {

                found = (int) false;

                // Loop over previous rigid constraints
                for (k = 0; k < i; k++) {
                    for (l = 0; l < tacsInstance->feaProblem.feaAnalysis[k].numRigidConstraint; l++) {

                        // If current rigid constraint was previous written - mark as found
                        if (strcmp(tacsInstance->feaProblem.feaAnalysis[i].rigidConstraint[j],
                                   tacsInstance->feaProblem.feaAnalysis[k].rigidConstraint[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // Make sure constraint isn't already in rigid variables too!
                for (k = 0; k < i; k++) {
                    for (l = 0; l < tacsInstance->feaProblem.feaAnalysis[k].numRigidVariable; l++) {

                        // If current rigid constraint was previous written - mark as found
                        if (strcmp(tacsInstance->feaProblem.feaAnalysis[i].rigidConstraint[j],
                                   tacsInstance->feaProblem.feaAnalysis[k].rigidVariable[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // If not write out an aestat card
                numAEStatSurf += 1;

                fprintf(fp,"%-8s", "AESTAT");

                tempString = convert_integerToString(numAEStatSurf, 7, 1);
                AIM_NOTNULL(tempString, aimInfo, status);
                fprintf(fp, "%s%s", delimiter, tempString);
                AIM_FREE(tempString);

                fprintf(fp, "%s%7s\n", delimiter, tacsInstance->feaProblem.feaAnalysis[i].rigidConstraint[j]);
            }
        }

        fprintf(fp,"\n");
    }

    // Analysis Cards - Eigenvalue and design objective included, as well as combined load, constraint, and design constraints
    for (i = 0; i < tacsInstance->feaProblem.numAnalysis; i++) {

        if (i == 0) printf("\tWriting analysis cards\n");

        status = nastran_writeAnalysisCard(fp,
                                           &tacsInstance->feaProblem.feaAnalysis[i],
                                           &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);

        if (tacsInstance->feaProblem.feaAnalysis[i].numLoad != 0) {

            // Create a temporary list of load IDs
            tempIntegerArray = (int *) EG_alloc(tacsInstance->feaProblem.feaAnalysis[i].numLoad*sizeof(int));
            if (tempIntegerArray == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            k = 0;
            for (j = 0; j <  tacsInstance->feaProblem.feaAnalysis[i].numLoad; j++) {
                for (l = 0; l < tacsInstance->feaProblem.numLoad; l++) {
                    if (tacsInstance->feaProblem.feaAnalysis[i].loadSetID[j] == tacsInstance->feaProblem.feaLoad[l].loadID) break;
                }
                if (tacsInstance->feaProblem.feaLoad[l].loadType == Thermal) continue;
                tempIntegerArray[j] = tacsInstance->feaProblem.feaLoad[l].loadID;
                k += 1;
            }

            tempIntegerArray = (int *) EG_reall(tempIntegerArray, k*sizeof(int));
            if (tempIntegerArray == NULL)  {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            // Write combined load card
            printf("\tWriting load ADD cards\n");
            status = nastran_writeLoadADDCard(fp,
                                               tacsInstance->feaProblem.numLoad+i+1,
                                               k,
                                               tempIntegerArray,
                                               tacsInstance->feaProblem.feaLoad,
                                               &tacsInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

            // Free temporary load ID list
            EG_free(tempIntegerArray);
            tempIntegerArray = NULL;

        } else { // If no loads for an individual analysis are specified assume that all loads should be applied

            if (tacsInstance->feaProblem.numLoad != 0) {

                // Create a temporary list of load IDs
                tempIntegerArray = (int *) EG_alloc(tacsInstance->feaProblem.numLoad*sizeof(int));
                if (tempIntegerArray == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                k = 0;
                for (j = 0; j < tacsInstance->feaProblem.numLoad; j++) {
                    if (tacsInstance->feaProblem.feaLoad[j].loadType == Gravity) continue;
                    tempIntegerArray[j] = tacsInstance->feaProblem.feaLoad[j].loadID;
                    k += 1;
                }

                AIM_ERROR(aimInfo, "Writing load ADD cards is not properly implemented!");
                status = CAPS_NOTIMPLEMENT;
                goto cleanup;

#ifdef FIX_tempIntegerArray_INIT
                //      tempIntegerArray needs to be initialized!!!

                tempIntegerArray = (int *) EG_reall(tempIntegerArray, k*sizeof(int));
                if (tempIntegerArray == NULL)  {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                //TOOO: eliminate load add card?
                // Write combined load card
                printf("\tWriting load ADD cards\n");
                status = nastran_writeLoadADDCard(fp,
                                                  tacsInstance->feaProblem.numLoad+i+1,
                                                  k,
                                                  tempIntegerArray,
                                                  tacsInstance->feaProblem.feaLoad,
                                                  &tacsInstance->feaProblem.feaFileFormat);
                AIM_STATUS(aimInfo, status);

                // Free temporary load ID list
                EG_free(tempIntegerArray);
                tempIntegerArray = NULL;
#endif
            }

        }

        if (tacsInstance->feaProblem.feaAnalysis[i].numConstraint != 0) {
            // Write combined constraint card
            printf("\tWriting constraint ADD cards\n");
            status = nastran_writeConstraintADDCard(fp,
                                                    tacsInstance->feaProblem.numConstraint+i+1,
                                                    tacsInstance->feaProblem.feaAnalysis[i].numConstraint,
                                                    tacsInstance->feaProblem.feaAnalysis[i].constraintSetID,
                                                    &tacsInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

        } else { // If no constraints for an individual analysis are specified assume that all constraints should be applied

            if (tacsInstance->feaProblem.numConstraint != 0) {

                printf("\tWriting combined constraint cards\n");

                // Create a temporary list of constraint IDs
                AIM_ALLOC(tempIntegerArray, tacsInstance->feaProblem.numConstraint, int, aimInfo, status);

                for (j = 0; j < tacsInstance->feaProblem.numConstraint; j++) {
                    tempIntegerArray[j] = tacsInstance->feaProblem.feaConstraint[j].constraintID;
                }

                // Write combined constraint card
                status = nastran_writeConstraintADDCard(fp,
                                                        tacsInstance->feaProblem.numConstraint+i+1,
                                                        tacsInstance->feaProblem.numConstraint,
                                                        tempIntegerArray,
                                                        &tacsInstance->feaProblem.feaFileFormat);
                AIM_STATUS(aimInfo, status);

                // Free temporary constraint ID list
                AIM_FREE(tempIntegerArray);
            }
        }

        if (tacsInstance->feaProblem.feaAnalysis[i].numDesignConstraint != 0) {

            // Write combined design constraint card
            printf("\tWriting design constraint ADD cards\n");
            status = nastran_writeDesignConstraintADDCard(fp,
                                                          tacsInstance->feaProblem.numDesignConstraint+i+1,
                                                          tacsInstance->feaProblem.feaAnalysis[i].numDesignConstraint,
                                                          tacsInstance->feaProblem.feaAnalysis[i].designConstraintSetID,
                                                          &tacsInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

        } else { // If no design constraints for an individual analysis are specified assume that all design constraints should be applied

            if (tacsInstance->feaProblem.numDesignConstraint != 0) {

                // Create a temporary list of design constraint IDs
                AIM_ALLOC(tempIntegerArray, tacsInstance->feaProblem.numDesignConstraint, int, aimInfo, status);

                for (j = 0; j < tacsInstance->feaProblem.numDesignConstraint; j++) {
                    tempIntegerArray[j] = tacsInstance->feaProblem.feaDesignConstraint[j].designConstraintID;
                }

                // Write combined design constraint card
                printf("\tWriting design constraint ADD cards\n");
                status = nastran_writeDesignConstraintADDCard(fp,
                                                              tacsInstance->feaProblem.numDesignConstraint+i+1,
                                                              tacsInstance->feaProblem.numDesignConstraint,
                                                              tempIntegerArray,
                                                              &tacsInstance->feaProblem.feaFileFormat);
                // Free temporary design constraint ID list
                if (tempIntegerArray != NULL) EG_free(tempIntegerArray);
                tempIntegerArray = NULL;
                AIM_STATUS(aimInfo, status);
            }

        }
    }


    // Loads
    for (i = 0; i < tacsInstance->feaProblem.numLoad; i++) {

        if (i == 0) printf("\tWriting load cards\n");

        status = nastran_writeLoadCard(fp,
                                       &tacsInstance->feaProblem.feaLoad[i],
                                       &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Constraints
    for (i = 0; i < tacsInstance->feaProblem.numConstraint; i++) {

        if (i == 0) printf("\tWriting constraint cards\n");

        status = nastran_writeConstraintCard(fp,
                                             &tacsInstance->feaProblem.feaConstraint[i],
                                             &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Supports
    for (i = 0; i < tacsInstance->feaProblem.numSupport; i++) {

        if (i == 0) printf("\tWriting support cards\n");
        j = (int) true;
        status = nastran_writeSupportCard(fp,
                                          &tacsInstance->feaProblem.feaSupport[i],
                                          &tacsInstance->feaProblem.feaFileFormat,
                                          &j);
        AIM_STATUS(aimInfo, status);
    }


    // Materials
    for (i = 0; i < tacsInstance->feaProblem.numMaterial; i++) {

        if (i == 0) printf("\tWriting material cards\n");

        status = nastran_writeMaterialCard(fp,
                                           &tacsInstance->feaProblem.feaMaterial[i],
                                           &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Properties
    for (i = 0; i < tacsInstance->feaProblem.numProperty; i++) {

        if (i == 0) printf("\tWriting property cards\n");

        status = nastran_writePropertyCard(fp,
                                           &tacsInstance->feaProblem.feaProperty[i],
                                           &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Coordinate systems
    // printf("DEBUG: Number of coord systems: %d\n", tacsInstance->feaProblem.numCoordSystem);
    for (i = 0; i < tacsInstance->feaProblem.numCoordSystem; i++) {

        if (i == 0) printf("\tWriting coordinate system cards\n");

        status = nastran_writeCoordinateSystemCard(fp,
                                                   &tacsInstance->feaProblem.feaCoordSystem[i],
                                                   &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design variables
    for( i = 0; i < tacsInstance->feaProblem.numDesignVariable; i++) {

        if (i == 0) printf("\tWriting design variable cards\n");

        status = nastran_writeDesignVariableCard(fp,
                                                 &tacsInstance->feaProblem.feaDesignVariable[i],
                                                 &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design variable relations
    for( i = 0; i < tacsInstance->feaProblem.numDesignVariableRelation; i++) {

        if (i == 0) printf("\tWriting design variable relation cards\n");

        status = nastran_writeDesignVariableRelationCard(aimInfo,
                                                         fp,
                                                         &tacsInstance->feaProblem.feaDesignVariableRelation[i],
                                                         &tacsInstance->feaProblem,
                                                         &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design constraints
    for( i = 0; i < tacsInstance->feaProblem.numDesignConstraint; i++) {

        if (i == 0) printf("\tWriting design constraints and responses cards\n");

        status = nastran_writeDesignConstraintCard(fp,
                                                   &tacsInstance->feaProblem.feaDesignConstraint[i],
                                                   &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design equations
    for( i = 0; i < tacsInstance->feaProblem.numEquation; i++) {

        if (i == 0) printf("\tWriting design equation cards\n");

        status = nastran_writeDesignEquationCard(fp,
                                                 &tacsInstance->feaProblem.feaEquation[i],
                                                 &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design table constants
    if (tacsInstance->feaProblem.feaDesignTable.numConstant > 0)
        printf("\tWriting design table card\n");
    status = nastran_writeDesignTableCard(fp,
                                          &tacsInstance->feaProblem.feaDesignTable,
                                          &tacsInstance->feaProblem.feaFileFormat);
    AIM_STATUS(aimInfo, status);

    // Optimization - design optimization parameters
    if (tacsInstance->feaProblem.feaDesignOptParam.numParam > 0)
        printf("\tWriting design optimization parameters card\n");
    status = nastran_writeDesignOptParamCard(fp,
                                             &tacsInstance->feaProblem.feaDesignOptParam,
                                             &tacsInstance->feaProblem.feaFileFormat);
    AIM_STATUS(aimInfo, status);

    // Optimization - design responses
    for( i = 0; i < tacsInstance->feaProblem.numDesignResponse; i++) {

        if (i == 0) printf("\tWriting design response cards\n");

        status = nastran_writeDesignResponseCard(fp,
                                                 &tacsInstance->feaProblem.feaDesignResponse[i],
                                                 &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design equation responses
    for( i = 0; i < tacsInstance->feaProblem.numEquationResponse; i++) {

        if (i == 0) printf("\tWriting design equation response cards\n");

        status = nastran_writeDesignEquationResponseCard(fp,
                                                 &tacsInstance->feaProblem.feaEquationResponse[i],
                                                 &tacsInstance->feaProblem,
                                                 &tacsInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Include mesh file
    fprintf(fp,"\nINCLUDE \'%s.bdf\'\n\n", tacsInstance->projectName);

    // End bulk data
    fprintf(fp,"ENDDATA\n");

    fclose(fp);
    fp = NULL;

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit in tacsAIM preAnalysis, status = %d\n",
               status);

    if (tempIntegerArray != NULL) EG_free(tempIntegerArray);

    if (fp != NULL) fclose(fp);

    return status;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(void *instStore, void *aimInfo,
       /*@unused@*/ int restart, capsValue *inputs)
{
  int status = CAPS_SUCCESS; // Function return status

  int i, j, k, idv, irow, icol, ibody; // Indexing
  int index, nGeomIn=0;

  char tmp[128];
  int numFunctional=0;
  int **functional_map=NULL;
  double **functional_xyz=NULL;
  double functional_dvar;

  int numNode = 0, *numPoint=NULL, *surf2tess=NULL;

  const char *name;
  char **names=NULL;
  double **dxyz = NULL;

  FILE *fp=NULL;
  capsValue *values=NULL, *geomInVal;
  aimStorage *tacsInstance;
  tacsInstance = (aimStorage*)instStore;
  AIM_NOTNULL(inputs, aimInfo, status);

  if (inputs[inDesign_Variable-1].nullVal == NotNull) {

    /* check for GeometryIn variables*/
    nGeomIn = 0;
    for (i = 0; i < tacsInstance->feaProblem.numDesignVariable; i++) {

      name = tacsInstance->feaProblem.feaDesignVariable[i].name;

      // Loop over the geometry in values and compute sensitivities for all bodies
      index = aim_getIndex(aimInfo, name, GEOMETRYIN);
      if (index == CAPS_NOTFOUND){
       continue;
      }
      if (index < CAPS_SUCCESS ) {
        status = index;
        AIM_STATUS(aimInfo, status);
      }

      if(aim_getGeomInType(aimInfo, index) != 0) {
          AIM_ERROR(aimInfo, "GeometryIn value %s is a configuration parameter and not a valid design parameter - can't get sensitivity\n",
                    name);
          status = CAPS_BADVALUE;
          goto cleanup;
      }

      nGeomIn++;
    }

    // No geometry sensitivties, so don't look for .sens file
    if (nGeomIn == 0) {
      status = CAPS_SUCCESS;
      goto cleanup;
    }

    numNode = 0;
    for (ibody = 0; ibody < tacsInstance->numMesh; ibody++) {
      numNode += tacsInstance->feaMesh[ibody].numNode;
    }

    AIM_ALLOC(surf2tess, 2*numNode, int, aimInfo, status);
    for (i = 0; i < 2*numNode; i++) surf2tess[i] = 0;

    numNode = 0;
    for (ibody = 0; ibody < tacsInstance->numMesh; ibody++) {
      for (i = 0; i < tacsInstance->feaMesh[ibody].numNode; i++) {
        surf2tess[2*(numNode+i)+0] = ibody;
        surf2tess[2*(numNode+i)+1] = i;
      }
      numNode += tacsInstance->feaMesh[ibody].numNode;
    }

    // Read <Proj_Name>.sens
    snprintf(tmp, 128, "%s%s", tacsInstance->projectName, ".sens");
    fp = aim_fopen(aimInfo, tmp, "r");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Unable to open: %s", tmp);
      status = CAPS_IOERR;
      goto cleanup;
    }

    // Number of nodes and functinoals in the file
    status = fscanf(fp, "%d", &numFunctional);
    if (status == EOF || status != 1) {
      AIM_ERROR(aimInfo, "Failed to read sens file number of functionals");
      status = CAPS_IOERR; goto cleanup;
    }

    AIM_ALLOC(numPoint, numFunctional, int, aimInfo, status);
    for (i = 0; i < numFunctional; i++) numPoint[i] = 0;

    AIM_ALLOC(functional_map, numFunctional, int*, aimInfo, status);
    for (i = 0; i < numFunctional; i++) functional_map[i] = NULL;

    AIM_ALLOC(functional_xyz, numFunctional, double*, aimInfo, status);
    for (i = 0; i < numFunctional; i++) functional_xyz[i] = NULL;

    AIM_ALLOC(names, numFunctional, char*, aimInfo, status);
    for (i = 0; i < numFunctional; i++) names[i] = NULL;

    AIM_ALLOC(values, numFunctional, capsValue, aimInfo, status);
    for (i = 0; i < numFunctional; i++) aim_initValue(&values[i]);

    for (i = 0; i < numFunctional; i++) {
      values[i].type = DoubleDeriv;

      /* allocate derivatives */
      AIM_ALLOC(values[i].derivs, tacsInstance->feaProblem.numDesignVariable, capsDeriv, aimInfo, status);
      for (idv = 0; idv < tacsInstance->feaProblem.numDesignVariable; idv++) {
        values[i].derivs[idv].name  = NULL;
        values[i].derivs[idv].deriv = NULL;
        values[i].derivs[idv].len_wrt = 0;
      }
      values[i].nderiv = tacsInstance->feaProblem.numDesignVariable;
    }

    // Read in Functional name, value, numPoint, and dFunctinoal/dxyz
    for (i = 0; i < numFunctional; i++) {

      status = fscanf(fp, "%s", tmp);
      if (status == EOF) {
        AIM_ERROR(aimInfo, "Failed to read sens file functional name");
        status = CAPS_IOERR; goto cleanup;
      }

      AIM_STRDUP(names[i], tmp, aimInfo, status);

      status = fscanf(fp, "%lf", &values[i].vals.real);
      if (status == EOF || status != 1) {
        AIM_ERROR(aimInfo, "Failed to read sens file functional value");
        status = CAPS_IOERR; goto cleanup;
      }

      status = fscanf(fp, "%d", &numPoint[i]);
      if (status == EOF || status != 1) {
        AIM_ERROR(aimInfo, "Failed to read sens file number of points");
        status = CAPS_IOERR; goto cleanup;
      }

      AIM_ALLOC(functional_map[i],   numPoint[i], int   , aimInfo, status);
      AIM_ALLOC(functional_xyz[i], 3*numPoint[i], double, aimInfo, status);

      for (j = 0; j < numPoint[i]; j++) {
        status = fscanf(fp, "%d %lf %lf %lf", &functional_map[i][j],
                                              &functional_xyz[i][3*j+0],
                                              &functional_xyz[i][3*j+1],
                                              &functional_xyz[i][3*j+2]);
        if (status == EOF || status != 4) {
          AIM_ERROR(aimInfo, "Failed to read sens file data");
          status = CAPS_IOERR; goto cleanup;
        }

        if (functional_map[i][j] < 1 || functional_map[i][j] > numNode) {
          AIM_ERROR(aimInfo, "sens file volume mesh vertex index: %d out-of-range [1-%d]", functional_map[i][j], numNode);
          status = CAPS_IOERR; goto cleanup;
        }
      }
    }

    AIM_ALLOC(dxyz, tacsInstance->numMesh, double*, aimInfo, status);
    for (ibody = 0; ibody < tacsInstance->numMesh; ibody++) dxyz[ibody] = NULL;

    /* set derivatives */
    for (idv = 0; idv < tacsInstance->feaProblem.numDesignVariable; idv++) {

      name = tacsInstance->feaProblem.feaDesignVariable[idv].name;

      // Loop over the geometry in values and compute sensitivities for all bodies
      index = aim_getIndex(aimInfo, name, GEOMETRYIN);
      status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
      if (status == CAPS_BADINDEX) continue;
      AIM_STATUS(aimInfo, status);

      for (i = 0; i < numFunctional; i++) {
        AIM_STRDUP(values[i].derivs[idv].name, name, aimInfo, status);

        AIM_ALLOC(values[i].derivs[idv].deriv, geomInVal->length, double, aimInfo, status);
        values[i].derivs[idv].len_wrt  = geomInVal->length;
        for (j = 0; j < geomInVal->length; j++)
          values[i].derivs[idv].deriv[j] = 0;
      }

      for (irow = 0; irow < geomInVal->nrow; irow++) {
        for (icol = 0; icol < geomInVal->ncol; icol++) {

          for (ibody = 0; ibody < tacsInstance->numMesh; ibody++) {
            status = aim_tessSensitivity(aimInfo,
                                         name,
                                         irow+1, icol+1, // row, col
                                         tacsInstance->feaMesh[ibody].egadsTess,
                                         &numNode, &dxyz[ibody]);
            AIM_STATUS(aimInfo, status, "Sensitivity for: %s\n", name);
            AIM_NOTNULL(dxyz[ibody], aimInfo, status);
          }

          for (i = 0; i < numFunctional; i++) {
            functional_dvar = values[i].derivs[idv].deriv[geomInVal->ncol*irow + icol];

            for (j = 0; j < numPoint[i]; j++) {
              k = functional_map[i][j]-1; // 1-based indexing
              ibody = surf2tess[2*k+0];
              k     = surf2tess[2*k+1];

              functional_dvar += functional_xyz[i][3*j+0]*dxyz[ibody][3*k + 0]  // dx/dGeomIn
                               + functional_xyz[i][3*j+1]*dxyz[ibody][3*k + 1]  // dy/dGeomIn
                               + functional_xyz[i][3*j+2]*dxyz[ibody][3*k + 2]; // dz/dGeomIn
            }
            values[i].derivs[idv].deriv[geomInVal->ncol*irow + icol] = functional_dvar;
          }

          for (ibody = 0; ibody < tacsInstance->numMesh; ibody++)
            AIM_FREE(dxyz[ibody]);
        }
      }
    }

    /* create the dynamic output */
    for (i = 0; i < numFunctional; i++) {
      status = aim_makeDynamicOutput(aimInfo, names[i], &values[i]);
      AIM_STATUS(aimInfo, status);
    }
  }

cleanup:

  if (fp != NULL) fclose(fp);

  if (functional_xyz != NULL)
    for (i = 0; i < numFunctional; i++)
      AIM_FREE(functional_xyz[i]);

  if (functional_map != NULL)
    for (i = 0; i < numFunctional; i++)
      AIM_FREE(functional_map[i]);

  if (names != NULL)
    for (i = 0; i < numFunctional; i++)
      AIM_FREE(names[i]);

  if (dxyz != NULL)
    for (i = 0; i < tacsInstance->numMesh; i++)
      AIM_FREE(dxyz[i]);

  AIM_FREE(functional_xyz);
  AIM_FREE(functional_map);
  AIM_FREE(names);
  AIM_FREE(dxyz);
  AIM_FREE(values);
  AIM_FREE(surf2tess);
  AIM_FREE(numPoint);

  return status;
}


// Set TACS output variables
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
    /*@unused@*/ int index, /*@unused@*/ char **aoname, /*@unused@*/ capsValue *form)
{
    /*! \page aimOutputsTACS AIM Outputs
     * The following list outlines the TACS outputs available through the AIM interface.
     */

    #ifdef DEBUG
        printf(" nastranAIM/aimOutputs index = %d!\n", index);
    #endif

    /*<--! \page aimOutputsTACS AIM Outputs-->
     * - <B>EigenValue</B> = List of Eigen-Values (\f$ \lambda\f$) after a modal solve.
     * - <B>EigenRadian</B> = List of Eigen-Values in terms of radians (\f$ \omega = \sqrt{\lambda}\f$ ) after a modal solve.
     * - <B>EigenFrequency</B> = List of Eigen-Values in terms of frequencies (\f$ f = \frac{\omega}{2\pi}\f$) after a modal solve.
     * - <B>EigenGeneralMass</B> = List of generalized masses for the Eigen-Values.
     * - <B>EigenGeneralStiffness</B> = List of generalized stiffness for the Eigen-Values.
     * - <B>Objective</B> = Final objective value for a design optimization case.
     * - <B>ObjectiveHistory</B> = List of objective value for the history of a design optimization case.
     */

#if 0
    if (index == 1) {
        *aoname = EG_strdup("EigenValue");

    } else if (index == 2) {
        *aoname = EG_strdup("EigenRadian");

    } else if (index == 3) {
        *aoname = EG_strdup("EigenFrequency");

    } else if (index == 4) {
        *aoname = EG_strdup("EigenGeneralMass");

    } else if (index == 5) {
        *aoname = EG_strdup("EigenGeneralStiffness");

    } else if (index == 6) {
        *aoname = EG_strdup("Objective");

    } else if (index == 7) {
        *aoname = EG_strdup("ObjectiveHistory");
    }

    form->type       = Double;
    form->units      = NULL;
    form->lfixed     = Change;
    form->sfixed     = Change;
    form->vals.reals = NULL;
    form->vals.real  = 0;
#endif

    return CAPS_SUCCESS;
}


// Calculate TACS output
int aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, /*@unused@*/ int index,
                  /*@unused@*/ capsValue *val)
{
    int status = CAPS_SUCCESS; // Function return status

#if 0
    int i; //Indexing

    int numData = 0;
    double *dataVector = NULL;
    double **dataMatrix = NULL;
    aimStorage *tacsInstance;

    char *filename = NULL; // File to open
    char extF06[] = ".f06", extOP2[] = ".op2";
    FILE *fp = NULL; // File pointer

    tacsInstance = (aimStorage *) instStore;

    if (index <= 5) {
        filename = (char *) EG_alloc((strlen(tacsInstance->projectName) +
                                      strlen(extF06) +1)*sizeof(char));
        if (filename == NULL) return EGADS_MALLOC;

        sprintf(filename, "%s%s", tacsInstance->projectName, extF06);

        fp = aim_fopen(aimInfo, filename, "r");

        EG_free(filename); // Free filename allocation

        if (fp == NULL) {
#ifdef DEBUG
            printf(" nastranAIM/aimCalcOutput Cannot open Output file!\n");
#endif
            return CAPS_IOERR;
        }

        status = nastran_readF06EigenValue(fp, &numData, &dataMatrix);
        if ((status == CAPS_SUCCESS) && (dataMatrix != NULL)) {

            val->nrow = numData;
            val->ncol = 1;
            val->length = val->nrow*val->ncol;
            if (val->length == 1) val->dim = Scalar;
            else val->dim = Vector;

            if (val->length == 1) {
                val->vals.real = dataMatrix[0][index-1];
            } else {

                val->vals.reals = (double *) EG_alloc(val->length*sizeof(double));
                if (val->vals.reals != NULL) {

                    for (i = 0; i < val->length; i++) {

                        val->vals.reals[i] = dataMatrix[i][index-1];
                    }

                } else status = EGADS_MALLOC;
            }
        }
    } else if (index == 6 || index == 7) {
        filename = (char *) EG_alloc((strlen(tacsInstance->projectName) +
                                      strlen(extOP2) +1)*sizeof(char));
        if (filename == NULL) return EGADS_MALLOC;

        sprintf(filename, "%s%s", tacsInstance->projectName, extOP2);

       status = nastran_readOP2Objective(filename, &numData, &dataVector);

        EG_free(filename); // Free filename allocation

        if (status == CAPS_SUCCESS && dataVector != NULL && numData > 0) {

            if (index == 6) val->nrow = 1;
            else val->nrow = numData;

            val->ncol = 1;
            val->length = val->nrow*val->ncol;

            if (val->length == 1) val->dim = Scalar;
            else val->dim = Vector;

            if (val->length == 1) {
                val->vals.real = dataVector[numData-1];
            } else {

                val->vals.reals = (double *) EG_alloc(val->length*sizeof(double));
                if (val->vals.reals != NULL) {

                    for (i = 0; i < val->length; i++) {

                        val->vals.reals[i] = dataVector[i];
                    }

                } else status = EGADS_MALLOC;
            }
        }
    }

    // Restore the path we came in with
    if (fp != NULL) fclose(fp);

    if (dataMatrix != NULL) {
        for (i = 0; i < numData; i++) {
            if (dataMatrix[i] != NULL) EG_free(dataMatrix[i]);
        }
        EG_free(dataMatrix);
    }

    if (dataVector != NULL) EG_free(dataVector);
#endif
    return status;
}


void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *tacsInstance;

#ifdef DEBUG
    printf(" nastranAIM/aimCleanup!\n");
#endif
    tacsInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(tacsInstance);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during clean up of instance\n", status);

    EG_free(tacsInstance);
}


int aimDiscr(char *tname, capsDiscr *discr)
{
    int        status = CAPS_SUCCESS; // Function return status
    int        i;
    ego        *tess = NULL;
    capsValue  *valMesh;
    aimStorage *tacsInstance;

#ifdef DEBUG
    printf(" nastranAIM/aimDiscr: tname = %s!\n", tname);
#endif
    if (tname == NULL) return CAPS_NOTFOUND;

    tacsInstance = (aimStorage *) discr->instStore;

    status = aim_getValue(discr->aInfo, Mesh, ANALYSISIN, &valMesh);
    AIM_STATUS(discr->aInfo, status);

    if (valMesh->nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(discr->aInfo, Mesh, "'Mesh' input must be linked to an output 'Surface_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get mesh
    tacsInstance->numMesh = valMesh->length;
    tacsInstance->feaMesh = (meshStruct *)valMesh->vals.AIMptr;
    AIM_NOTNULL(tacsInstance->feaMesh, discr->aInfo, status);
    AIM_STATUS(discr->aInfo, status);

    AIM_ALLOC(tess, tacsInstance->numMesh, ego, discr->aInfo, status);
    for (i = 0; i < tacsInstance->numMesh; i++) {
      tess[i] = tacsInstance->feaMesh[i].egadsTess;
    }

    status = mesh_fillDiscr(tname, &tacsInstance->groupMap, tacsInstance->numMesh, tess, discr);
    AIM_STATUS(discr->aInfo, status);

#ifdef DEBUG
    printf(" nastranAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
#endif

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit: function aimDiscr nastranAIM status = %d",
               status);

    AIM_FREE(tess);
    return status;
}

int aimTransfer(capsDiscr *discr, const char *dataName, int numPoint,
                int dataRank, double *dataVal, /*@unused@*/ char **units)
{

    /*! \page dataTransferTACS TACS Data Transfer
     *
     * The TACS AIM has the ability to transfer displacements and eigenvectors from the AIM and pressure
     * distributions to the AIM using the conservative and interpolative data transfer schemes in CAPS.
     */

  /*
     * \section dataFromTACS Data transfer from TACS (FieldOut)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements from the *.f06 file.
     * </ul>
     *
     * <ul>
     *  <li> <B>"EigenVector_#"</B> </li> <br>
     *   Retrieves modal eigen-vectors from the *.f06 file, where "#" should be replaced by the
     *   corresponding mode number for the eigen-vector (eg. EigenVector_3 would correspond to the third mode,
     *   while EigenVector_6 would be the sixth mode).
     * </ul>
     *
     * \section dataToTACS Data transfer to TACS (FieldIn)
     * <ul>
     *  <li> <B>"Pressure"</B> </li> <br>
     *  Writes appropriate load cards using the provided pressure distribution.
     * </ul>
     *
     */

    int status; // Function return status
    int i, j, dataPoint, bIndex; // Indexing
    aimStorage *tacsInstance;

    char *extF06 = ".f06";

    // FO6 data variables
    int numGridPoint = 0;
    int numEigenVector = 0;

    double **dataMatrix = NULL;

    // Specific EigenVector to use
    int eigenVectorIndex = 0;

    // Variables used in global node mapping
    //int *storage;
    int globalNodeID;

    // Filename stuff
    char *filename = NULL;
    FILE *fp; // File pointer

#ifdef DEBUG
    printf(" nastranAIM/aimTransfer name = %s  npts = %d/%d!\n",
           dataName, numPoint, dataRank);
#endif
    tacsInstance = (aimStorage *) discr->instStore;

    //Get the appropriate parts of the tessellation to data
    //storage = (int *) discr->ptrm;
    //capsGroupList = &storage[0]; // List of boundary ID (groupMap) in the transfer

    if (strcasecmp(dataName, "Displacement") != 0 &&
        strncmp(dataName, "EigenVector", 11) != 0) {

        printf("Unrecognized data transfer variable - %s\n", dataName);
        return CAPS_NOTFOUND;
    }

    filename = (char *) EG_alloc((strlen(tacsInstance->projectName) +
                                  strlen(extF06) + 1)*sizeof(char));
    if (filename == NULL) return EGADS_MALLOC;
    sprintf(filename,"%s%s", tacsInstance->projectName, extF06);

    // Open file
    fp = aim_fopen(discr->aInfo, filename, "r");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        if (filename != NULL) EG_free(filename);
        return CAPS_IOERR;
    }

    if (filename != NULL) EG_free(filename);
    filename = NULL;

    if (strcasecmp(dataName, "Displacement") == 0) {

        if (dataRank != 3) {

            printf("Invalid rank for dataName \"%s\" - excepted a rank of 3!!!\n",
                   dataName);
            status = CAPS_BADRANK;

        } else {

            status = nastran_readF06Displacement(fp,
                                                 -1,
                                                 &numGridPoint,
                                                 &dataMatrix);
            fclose(fp);
            fp = NULL;
        }

    } else if (strncmp(dataName, "EigenVector", 11) == 0) {

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

            printf("Invalid rank for dataName \"%s\" - excepted a rank of 3!!!\n",
                   dataName);
            status = CAPS_BADRANK;

        } else {

            status = nastran_readF06EigenVector(fp,
                                                &numEigenVector,
                                                &numGridPoint,
                                                &dataMatrix);
        }

        fclose(fp);
        fp = NULL;

    } else {

        status = CAPS_NOTFOUND;
    }

    AIM_STATUS(discr->aInfo, status);


    // Check EigenVector range
    if (strncmp(dataName, "EigenVector", 11) == 0)  {
        if (eigenVectorIndex > numEigenVector) {
            printf("Only %d EigenVectors found but index %d requested!\n",
                   numEigenVector, eigenVectorIndex);
            status = CAPS_RANGEERR;
            goto cleanup;
        }

        if (eigenVectorIndex < 1) {
            printf("For EigenVector_X notation, X must be >= 1, currently X = %d\n",
                   eigenVectorIndex);
            status = CAPS_RANGEERR;
            goto cleanup;
        }

    }
    if (dataMatrix == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    for (i = 0; i < numPoint; i++) {

        bIndex       = discr->tessGlobal[2*i  ];
        globalNodeID = discr->tessGlobal[2*i+1] +
                       discr->bodys[bIndex-1].globalOffset;

        if (strcasecmp(dataName, "Displacement") == 0) {

            for (dataPoint = 0; dataPoint < numGridPoint; dataPoint++) {
                if ((int) dataMatrix[dataPoint][0] ==  globalNodeID) break;
            }

            if (dataPoint == numGridPoint) {
                printf("Unable to locate global ID = %d in the data matrix\n",
                       globalNodeID);
                status = CAPS_NOTFOUND;
                goto cleanup;
            }

            dataVal[dataRank*i+0] = dataMatrix[dataPoint][2]; // T1
            dataVal[dataRank*i+1] = dataMatrix[dataPoint][3]; // T2
            dataVal[dataRank*i+2] = dataMatrix[dataPoint][4]; // T3

        } else if (strncmp(dataName, "EigenVector", 11) == 0) {


            for (dataPoint = 0; dataPoint < numGridPoint; dataPoint++) {
                if ((int) dataMatrix[eigenVectorIndex - 1][8*dataPoint + 0] ==
                    globalNodeID) break;
            }

            if (dataPoint == numGridPoint) {
                printf("Unable to locate global ID = %d in the data matrix\n",
                       globalNodeID);
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

cleanup:

    if (fp != NULL) fclose(fp);
    // Free data matrix
    if (dataMatrix != NULL) {
        j = 0;
        if      (strcasecmp(dataName, "Displacement") == 0) j = numGridPoint;
        else if (strncmp(dataName, "EigenVector", 11) == 0) j = numEigenVector;

        for (i = 0; i < j; i++) {
            if (dataMatrix[i] != NULL) EG_free(dataMatrix[i]);
        }
        EG_free(dataMatrix);
    }

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
    printf(" nastranAIM/aimLocateElement !\n");
#endif

    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}


int aimInterpolation(capsDiscr *discr, const char *name,
                     int bIndex, int eIndex, double *bary,
                     int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" nastranAIM/aimInterpolation  %s!\n", name);
#endif

    return aim_interpolation(discr, name, bIndex, eIndex, bary, rank, data,
                             result);
}


int aimInterpolateBar(capsDiscr *discr, const char *name,
                      int bIndex, int eIndex, double *bary,
                      int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" nastranAIM/aimInterpolateBar  %s!\n", name);
#endif

    return aim_interpolateBar(discr, name, bIndex, eIndex, bary, rank, r_bar,
                              d_bar);
}


int aimIntegration(capsDiscr *discr, const char *name,int bIndex, int eIndex,
                   int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" nastranAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                    int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" nastranAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank, r_bar, d_bar);
}
