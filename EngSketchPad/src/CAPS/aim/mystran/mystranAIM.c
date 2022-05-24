/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             mystran AIM
 *
 *     Written by Dr. Ryan Durscher AFRL/RQVC
 *
 *     This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */


/*! \mainpage Introduction
 * \tableofcontents
 * \section overviewMYSTRAN MYSTRAN AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with the finite element structural solver MYSTRAN \cite MYSTRAN. MYSTRAN is an open source,
 * general purpose, linear finite element analysis computer program written by Dr. Bill Case. Available at,
 * http://www.mystran.com/, MYSTRAN currently supports Linux and Windows operating systems.
 *
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsMYSTRAN and
 * \ref aimOutputsMYSTRAN and \ref attributeMYSTRAN, respectively.
 *
 * The MYSTRAN AIM can automatically execute MYSTRAN, with details provided in \ref aimExecuteMYSTRAN.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferMYSTRAN
 *
 *\section mystranExamples Examples
 * An example problem using the MYSTRAN AIM may be found at \ref mystranExample.
 *
 * \section clearanceMYSTRAN Clearance Statement
 *  This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */


/*! \page attributeMYSTRAN MYSTRAN AIM attributes
 * The following list of attributes are required for the MYSTRAN AIM inside the geometry input.
 *
 * - <b> capsAIM</b> This attribute is a CAPS requirement to indicate the analysis the geometry
 * representation supports.
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
 * - <b> capsIgnore</b> It is possible that there is a geometric body (or entity) that you do not want the MYSTRAN AIM to pay attention to when creating
 * a finite element model. The capsIgnore attribute allows a body (or entity) to be in the geometry and ignored by the AIM.  For example,
 * because of limitations in OpenCASCADE a situation where two edges are overlapping may occur; capsIgnore allows the user to only
 *  pay attention to one of the overlapping edges.
 *
 * - <b> capsBound </b> This is used to mark surfaces on the structural grid in which data transfer with an external
 * solver will take place. See \ref dataTransferMYSTRAN for additional details.
 *
 *
 */

// No connections yet
/* - <b> capsConnect</b> This is a name assigned to any geometric body where the user wishes to create
 * "fictitious" connections such as springs, dampers, and/or rigid body connections to. The user must manually specify
 * the connection between two <c>capsConnect</c> entities using the "Connect" tuple (see \ref aimInputsMYSTRAN).
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsConnect $springStart</c>.
 *
 * - <b> capsConnectLink</b> Similar to <c>capsConnect</c>, this is a name assigned to any geometric body where
 * the user wishes to create "fictitious" connections to. A connection is automatically made if a <c>capsConnectLink</c>
 * matches a <c>capsConnect</c> group. Again further specifics of the connection are input using the "Connect"
 * tuple (see \ref aimInputsMYSTRAN). Recall that a string in ESP starts with a $.
 * For example, attribute <c>capsConnectLink $springEnd</c>.
 */

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"
#include "miscUtils.h"
#include "feaUtils.h"
#include "nastranUtils.h"
#include "mystranUtils.h"

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

//#define DEBUG

enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Tess_Params,
  Edge_Point_Min,
  Edge_Point_Max,
  Quad_Mesh,
  Property,
  Material,
  Constraint,
  Load,
  Analysix,
  Analysis_Type,
  Support,
  Mesh,
  NUMINPUT = Mesh           /* Total number of inputs */
};

#define NUMOUTPUT  4


typedef struct {
    // Project name
    char *projectName;

    feaUnitsStruct units; // units system

    feaProblemStruct feaProblem;

    // Attribute to index map
    mapAttrToIndexStruct attrMap;

    // Attribute to constraint index map
    mapAttrToIndexStruct constraintMap;

    // Attribute to load index map
    mapAttrToIndexStruct loadMap;

    /*
    // Check to make sure data transfer is ok
    int dataTransferCheck;
     */

    // Mesh holders
    int numMesh;
    meshStruct *feaMesh;

} aimStorage;


static int initiate_aimStorage(aimStorage *mystranInstance)
{

    int status;

    // Set initial values for mystranInstance
    mystranInstance->projectName = NULL;
  
    // Mesh holders
    mystranInstance->numMesh = 0;
    mystranInstance->feaMesh = NULL;
  
    /*
    // Check to make sure data transfer is ok
    mystranInstance->dataTransferCheck = (int) true;
     */

    status = initiate_feaUnitsStruct(&mystranInstance->units);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&mystranInstance->attrMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to constraint index map
    status = initiate_mapAttrToIndexStruct(&mystranInstance->constraintMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to load index map
    status = initiate_mapAttrToIndexStruct(&mystranInstance->loadMap);
    if (status != CAPS_SUCCESS) return status;

    status = initiate_feaProblemStruct(&mystranInstance->feaProblem);
    if (status != CAPS_SUCCESS) return status;

    return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *mystranInstance)
{

    int status;
    int i;
  
    status = destroy_feaUnitsStruct(&mystranInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&mystranInstance->attrMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to constraint index map
    status = destroy_mapAttrToIndexStruct(&mystranInstance->constraintMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to load index map
    status = destroy_mapAttrToIndexStruct(&mystranInstance->loadMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (mystranInstance->feaMesh != NULL) {

        for (i = 0; i < mystranInstance->numMesh; i++) {
            status = destroy_meshStruct(&mystranInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        AIM_FREE(mystranInstance->feaMesh);
    }

    mystranInstance->feaMesh = NULL;
    mystranInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&mystranInstance->feaProblem);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    // NULL projetName
    mystranInstance->projectName = NULL;

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(void *aimInfo, aimStorage *mystranInstance)
{
  // Function return flag
  int status = CAPS_SUCCESS;
  int i, remesh = (int)true;

  // Meshing related variables
  double tessParam[3] = {0.025, 0.001, 15};
  int edgePointMin = 2;
  int edgePointMax = 50;
  int quadMesh = (int) false;

  // Attribute to transfer map
  mapAttrToIndexStruct transferMap;

  // Attribute to connect map
  mapAttrToIndexStruct connectMap;

  // analysis input values
  capsValue *TessParams;
  capsValue *EdgePoint_Min;
  capsValue *EdgePoint_Max;
  capsValue *QuadMesh;

  for (i = 0; i < mystranInstance->numMesh; i++) {
      remesh = remesh && (mystranInstance->feaMesh[i].egadsTess->oclass == EMPTY);
  }
  if (remesh == (int) false) return CAPS_SUCCESS;

  // retrieve or create the mesh from fea_createMesh
  status = aim_getValue(aimInfo, Tess_Params,    ANALYSISIN, &TessParams);
  if (status != CAPS_SUCCESS) return status;

  status = aim_getValue(aimInfo, Edge_Point_Min, ANALYSISIN, &EdgePoint_Min);
  if (status != CAPS_SUCCESS) return status;

  status = aim_getValue(aimInfo, Edge_Point_Max, ANALYSISIN, &EdgePoint_Max);
  if (status != CAPS_SUCCESS) return status;

  status = aim_getValue(aimInfo, Quad_Mesh,      ANALYSISIN, &QuadMesh);
  if (status != CAPS_SUCCESS) return status;

  if (TessParams != NULL) {
      tessParam[0] = TessParams->vals.reals[0];
      tessParam[1] = TessParams->vals.reals[1];
      tessParam[2] = TessParams->vals.reals[2];
  }

  // Max and min number of points
  if (EdgePoint_Min != NULL && EdgePoint_Min->nullVal != IsNull) {
      edgePointMin = EdgePoint_Min->vals.integer;
      if (edgePointMin < 2) {
        AIM_ANALYSISIN_ERROR(aimInfo, Edge_Point_Min, "Edge_Point_Min = %d must be greater or equal to 2\n", edgePointMin);
        return CAPS_BADVALUE;
      }
  }

  if (EdgePoint_Max != NULL && EdgePoint_Max->nullVal != IsNull) {
      edgePointMax = EdgePoint_Max->vals.integer;
      if (edgePointMax < 2) {
        AIM_ANALYSISIN_ERROR(aimInfo, Edge_Point_Max, "Edge_Point_Max = %d must be greater or equal to 2\n", edgePointMax);
        return CAPS_BADVALUE;
      }
  }

  if (edgePointMin >= 2 && edgePointMax >= 2 && edgePointMin > edgePointMax) {
    AIM_ERROR  (aimInfo, "Edge_Point_Max must be greater or equal Edge_Point_Min");
    AIM_ADDLINE(aimInfo, "Edge_Point_Max = %d, Edge_Point_Min = %d\n",edgePointMax,edgePointMin);
    return CAPS_BADVALUE;
  }

  if (QuadMesh != NULL) quadMesh = QuadMesh->vals.integer;

  status = initiate_mapAttrToIndexStruct(&transferMap);
  if (status != CAPS_SUCCESS) return status;

  status = initiate_mapAttrToIndexStruct(&connectMap);
  if (status != CAPS_SUCCESS) return status;
/*@-nullpass@*/
  status = fea_createMesh(aimInfo,
                          tessParam,
                          edgePointMin,
                          edgePointMax,
                          quadMesh,
                          &mystranInstance->attrMap,
                          &mystranInstance->constraintMap,
                          &mystranInstance->loadMap,
                          &transferMap,
                          &connectMap,
                          NULL,
                          &mystranInstance->numMesh,
                          &mystranInstance->feaMesh,
                          &mystranInstance->feaProblem);
/*@+nullpass@*/
  if (status != CAPS_SUCCESS) return status;

  status = destroy_mapAttrToIndexStruct(&transferMap);
  if (status != CAPS_SUCCESS) return status;

  status = destroy_mapAttrToIndexStruct(&connectMap);
  if (status != CAPS_SUCCESS) return status;

  return CAPS_SUCCESS;
}


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int  status = CAPS_SUCCESS, *ints=NULL, i;
    char **strs=NULL;

    aimStorage *mystranInstance=NULL;

#ifdef DEBUG
    printf("mystranAIM/aimInitialize  inst = %d!\n", inst);
#endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 4;

    /* specify the name of each field variable */
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);

    strs[0]  = EG_strdup("Displacement");
    strs[1]  = EG_strdup("EigenVector");
    strs[2]  = EG_strdup("EigenVector_#");
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

    // Allocate mystranInstance
    AIM_ALLOC(mystranInstance, 1, aimStorage, aimInfo, status);
    *instStore = mystranInstance;

    // Initialize instance storage
    status = initiate_aimStorage(mystranInstance);
    AIM_STATUS(aimInfo, status);

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
    /*! \page aimInputsMYSTRAN AIM Inputs
     * The following list outlines the MYSTRAN inputs along with their default value available
     * through the AIM interface. Unless noted these values will be not be linked to
     * any parent AIMs with variables of the same name.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" mystranAIM/aimInputs index = %d!\n", index);
#endif

    *ainame = NULL;

    // MYSTRAN Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("mystran_CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsMYSTRAN
         * - <B> Proj_Name = "mystran_CAPS"</B> <br>
         * This corresponds to the project name used for file naming.
         */

    } else if (index == Tess_Params) {
        *ainame               = EG_strdup("Tess_Params");
        defval->type          = Double;
        defval->dim           = Vector;
        defval->nrow          = 3;
        defval->ncol          = 1;
        defval->units         = NULL;
        defval->lfixed        = Fixed;
        defval->vals.reals    = (double *) EG_alloc(defval->nrow*sizeof(double));
        if (defval->vals.reals != NULL) {
            defval->vals.reals[0] = 0.025;
            defval->vals.reals[1] = 0.001;
            defval->vals.reals[2] = 15.00;
        } else return EGADS_MALLOC;

        /*! \page aimInputsMYSTRAN
         * - <B> Tess_Params = [0.025, 0.001, 15.0]</B> <br>
         * Body tessellation parameters used when creating a boundary element model.
         * Tess_Params[0] and Tess_Params[1] get scaled by the bounding
         * box of the body. (From the EGADS manual) A set of 3 parameters that drive the EDGE discretization
         * and the FACE triangulation. The first is the maximum length of an EDGE segment or triangle side
         * (in physical space). A zero is flag that allows for any length. The second is a curvature-based
         * value that looks locally at the deviation between the centroid of the discrete object and the
         * underlying geometry. Any deviation larger than the input value will cause the tessellation to
         * be enhanced in those regions. The third is the maximum interior dihedral angle (in degrees)
         * between triangle facets (or Edge segment tangents for a WIREBODY tessellation), note that a
         * zero ignores this phase
         */

    } else if (index == Edge_Point_Min) {
        *ainame               = EG_strdup("Edge_Point_Min");
        defval->type          = Integer;
        defval->vals.integer  = 2;
        defval->lfixed        = Fixed;
        defval->nrow          = 1;
        defval->ncol          = 1;
        defval->nullVal       = NotNull;

        /*! \page aimInputsMYSTRAN
         * - <B> Edge_Point_Min = 2</B> <br>
         * Minimum number of points on an edge including end points to use when creating a surface mesh (min 2).
         */

    } else if (index == Edge_Point_Max) {
        *ainame               = EG_strdup("Edge_Point_Max");
        defval->type          = Integer;
        defval->vals.integer  = 50;
        defval->lfixed        = Fixed;
        defval->nrow          = 1;
        defval->ncol          = 1;
        defval->nullVal       = NotNull;

        /*! \page aimInputsMYSTRAN
         * - <B> Edge_Point_Max = 50</B> <br>
         * Maximum number of points on an edge including end points to use when creating a surface mesh (min 2).
         */

    } else if (index == Quad_Mesh) {
        *ainame               = EG_strdup("Quad_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = (int) false;

        /*! \page aimInputsMYSTRAN
         * - <B> Quad_Mesh = False</B> <br>
         * Create a quadratic mesh on four edge faces when creating the boundary element model.
         */

    } else if (index == Property) {
        *ainame              = EG_strdup("Property");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsMYSTRAN
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

        /*! \page aimInputsMYSTRAN
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

        /*! \page aimInputsMYSTRAN
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

        /*! \page aimInputsMYSTRAN
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

        /*! \page aimInputsMYSTRAN
         * - <B> Analysis = NULL</B> <br>
         * Analysis tuple used to input analysis/case information for the model, see \ref feaAnalysis for additional details.
         */
    } else if (index == Analysis_Type) {
        *ainame              = EG_strdup("Analysis_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Modal");
        defval->lfixed       = Change;

        /*! \page aimInputsMYSTRAN
         * - <B> Analysis_Type = "Modal"</B> <br>
         * Type of analysis to generate files for, options include "Modal", "Static", and "Craig-Bampton".
         */
    } else if (index == Support) {
        *ainame              = EG_strdup("Support");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsMYSTRAN
         * - <B> Support = NULL</B> <br>
         * Support tuple used to input support information for the model, see \ref feaSupport for additional details.
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

        /*! \page aimInputsMYSTRAN
         * - <B>Mesh = NULL</B> <br>
         * A Mesh link.
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}



int aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{

    int i, j, k; // Indexing

    int status; // Status return

    int found; // Boolean operator

    int *tempIntegerArray = NULL; // Temporary array to store a list of integers

    // Analyis information
    char *analysisType = NULL;

    // File IO
    char filename[PATH_MAX]; // Output file name
    FILE *fp = NULL; // Output file pointer
  
    aimStorage *mystranInstance;
  
    mystranInstance = (aimStorage *) instStore;
    if (aimInputs == NULL) return CAPS_NULLVALUE;

    // Get project name
    mystranInstance->projectName = aimInputs[Proj_Name-1].vals.string;

    // Check and generate/retrieve the mesh
    status = checkAndCreateMesh(aimInfo, mystranInstance);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Note: Setting order is important here.
    // 1. Materials should be set before properties.
    // 2. Coordinate system should be set before mesh and loads
    // 3. Mesh should be set before loads and constraints
    // 4. Constraints and loads should be set before analysis

    // Set material properties
    if (aimInputs[Material-1].nullVal == NotNull) {
        status = fea_getMaterial(aimInfo,
                                 aimInputs[Material-1].length,
                                 aimInputs[Material-1].vals.tuple,
                                 &mystranInstance->units,
                                 &mystranInstance->feaProblem.numMaterial,
                                 &mystranInstance->feaProblem.feaMaterial);
        if (status != CAPS_SUCCESS) goto cleanup;
    } else printf("\nLoad tuple is NULL - No materials set\n");

    // Set property properties
    if (aimInputs[Property-1].nullVal == NotNull) {
        status = fea_getProperty(aimInfo,
                                 aimInputs[Property-1].length,
                                 aimInputs[Property-1].vals.tuple,
                                 &mystranInstance->attrMap,
                                 &mystranInstance->units,
                                 &mystranInstance->feaProblem);
        if (status != CAPS_SUCCESS) goto cleanup;
    } else printf("\nProperty tuple is NULL - No properties set\n");

    // Set constraint properties
    if (aimInputs[Constraint-1].nullVal == NotNull) {
        status = fea_getConstraint(aimInputs[Constraint-1].length,
                                   aimInputs[Constraint-1].vals.tuple,
                                   &mystranInstance->constraintMap,
                                   &mystranInstance->feaProblem);
        if (status != CAPS_SUCCESS) goto cleanup;
    } else printf("\nConstraint tuple is NULL - No constraints applied\n");

    // Set support properties
    if (aimInputs[Support-1].nullVal == NotNull) {
        status = fea_getSupport(aimInputs[Support-1].length,
                                aimInputs[Support-1].vals.tuple,
                                &mystranInstance->constraintMap,
                                &mystranInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Support tuple is NULL - No supports applied\n");

    // Set load properties
    if (aimInputs[Load-1].nullVal == NotNull) {
        status = fea_getLoad(aimInputs[Load-1].length,
                             aimInputs[Load-1].vals.tuple,
                             &mystranInstance->loadMap,
                             &mystranInstance->feaProblem);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Loop through loads to see if any of them are supposed to be from an external source
        for (i = 0; i < mystranInstance->feaProblem.numLoad; i++) {

            if (mystranInstance->feaProblem.feaLoad[i].loadType == PressureExternal) {

                // Transfer external pressures from the AIM discrObj
                status = fea_transferExternalPressure(aimInfo,
                                                      &mystranInstance->feaProblem.feaMesh,
                                                      &mystranInstance->feaProblem.feaLoad[i]);
                if (status != CAPS_SUCCESS) goto cleanup;
            }
        }
    } else printf("\nLoad tuple is NULL - No loads applied\n");

    // Set analysis settings
    if (aimInputs[Analysix-1].nullVal == NotNull) {
        status = fea_getAnalysis(aimInputs[Analysix-1].length,
                                 aimInputs[Analysix-1].vals.tuple,
                                 &mystranInstance->feaProblem);
        if (status != CAPS_SUCCESS) goto cleanup; // It ok to not have an analysis tuple
    } else printf("\nAnalysis tuple is NULL\n");

    // Write Nastran Mesh
    status = mesh_writeNASTRAN(aimInfo,
                               mystranInstance->projectName,
                               1,
                               &mystranInstance->feaProblem.feaMesh,
                               mystranInstance->feaProblem.feaFileFormat.gridFileType,
                               1.0);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Write mystran input file
    strcpy(filename, mystranInstance->projectName);
    strcat(filename, ".dat");

    printf("\nWriting MYSTRAN instruction file....\n");
    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        status = CAPS_IOERR;
        goto cleanup;
    }

    //////////////// Executive control ////////////////
    fprintf(fp, "ID CAPS generated Problem FOR MYSTRAN\n");

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    if     (strcasecmp(analysisType, "Modal")         == 0) fprintf(fp, "SOL 3\n");
    else if(strcasecmp(analysisType, "Static")        == 0) fprintf(fp, "SOL 1\n");
    else if(strcasecmp(analysisType, "Craig-Bampton") == 0) fprintf(fp, "SOL 31\n");
    else {
        AIM_ERROR(aimInfo, "Unrecognized \"Analysis_Type\", %s", analysisType);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (strcasecmp(analysisType, "Modal") == 0) fprintf(fp, "OUTPUT4 EIGEN_VAL, EIGEN_VEC, GEN_MASS, , // -1/21 $\n"); // Binary output of eigenvalues and vectors

    fprintf(fp, "CEND\n\n");

    //////////////// Case control ////////////////

    // Write output request information
    fprintf(fp, "DISP = ALL\n"); // Output all displacements

    if(strcasecmp(analysisType, "Static") == 0) fprintf(fp, "STRE = ALL\n"); // Output all stress

    if(strcasecmp(analysisType, "Static") == 0) fprintf(fp, "STRA = ALL\n"); // Output all strain

    //fprintf(fp, "ELDATA(6,PRINT) = ALL\n");

    // Check thermal load - currently only a single thermal load is supported - also it only works for the global level - no subcase
    found = (int) false;
    for (i = 0; i < mystranInstance->feaProblem.numLoad; i++) {

        if (mystranInstance->feaProblem.feaLoad[i].loadType != Thermal) continue;

        if (found == (int) true) {
            printf("More than 1 Thermal load found - mystranAIM does NOT currently doesn't support multiple thermal loads!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        found = (int) true;

        fprintf(fp, "TEMPERATURE = %d\n",
                mystranInstance->feaProblem.feaLoad[i].loadID);
    }

    // Write constraint information
    if (mystranInstance->feaProblem.numConstraint != 0) {
        fprintf(fp, "SPC = %d\n", mystranInstance->feaProblem.numConstraint+1);
    } else {
        printf("Warning: No constraints specified for job!!!!\n");
    }

    // Modal analysis - only
    // If modal - we are only going to use the first analysis structure we come across that has its type as modal
    if (strcasecmp(analysisType, "Modal") == 0) {

        // Look through analysis structures for a modal one
        found = (int) false;
        for (i = 0; i < mystranInstance->feaProblem.numAnalysis; i++) {
            if (mystranInstance->feaProblem.feaAnalysis[i].analysisType == Modal) {
                found = (int) true;
                break;
            }
        }

        // Write out analysis ID if a modal analysis structure was found
        if (found == (int) true)  {
            fprintf(fp, "METHOD = %d\n",
                    mystranInstance->feaProblem.feaAnalysis[i].analysisID);
        } else {
            printf("Warning: No eigenvalue analysis information specified in \"Analysis\" tuple, though "
                    "AIM input \"Analysis_Type\" is set to \"Modal\"!!!!\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

    }
    // Static analysis - only
    if (strcasecmp(analysisType, "Static") == 0) {

        // If we have multiple analysis structures
        if (mystranInstance->feaProblem.numAnalysis != 0) {

            // Write subcase information if multiple analysis tuples were provide
            for (i = 0; i < mystranInstance->feaProblem.numAnalysis; i++) {

                if (mystranInstance->feaProblem.feaAnalysis[i].analysisType == Static) {
                    fprintf(fp, "SUBCASE %d\n", i);
                    fprintf(fp, "\tLABEL %s\n",
                            mystranInstance->feaProblem.feaAnalysis[i].name);

                    // Write loads for sub-case
                    if (mystranInstance->feaProblem.numLoad != 0) {
                        fprintf(fp, "\tLOAD = %d\n",
                                mystranInstance->feaProblem.numLoad+i+1);
                    }

                    // Issue some warnings regarding loads if necessary
                    if (mystranInstance->feaProblem.feaAnalysis[i].numLoad == 0 &&
                        mystranInstance->feaProblem.numLoad !=0) {
                        printf("Warning: No loads specified for static case %s, assuming "
                                "all loads are applied!!!!\n",
                               mystranInstance->feaProblem.feaAnalysis[i].name);
                    } else if (mystranInstance->feaProblem.numLoad == 0) {
                        printf("Warning: No loads specified for static case %s!!!!\n",
                               mystranInstance->feaProblem.feaAnalysis[i].name);
                    }
                }
            }

        } else { // // If no sub-cases

            if (mystranInstance->feaProblem.numLoad != 0) {
                fprintf(fp, "LOAD = %d\n", mystranInstance->feaProblem.numLoad+1);
            } else {
                printf("Warning: No loads specified for static a job!!!!\n");
            }
        }
    }

    //////////////// Bulk data ////////////////
    fprintf(fp, "\nBEGIN BULK\n");

    // Turn off auto SPC
    //fprintf(fp, "%-8s %7s %7s\n", "PARAM", "AUTOSPC", "N");

    //fprintf(fp, "%-8s %7s %7d\n", "PARAM", "PRTOU4", 1);
    //fprintf(fp, "%-8s %7s %7d\n", "PARAM", "PRTDLR", 1);

    //fprintf(fp, "%-8s %7s %7d\n", "PARAM", "PRTDOF", 1);

    //fprintf(fp, "%-8s %7d %7d\n", "DEBUG", 183, 1);
    // Turn on FEMAP
    //fprintf(fp, "%-8s %7s %7d\n", "PARAM", "POST", -1);

    // Analysis Cards - Eigenvalue included, as well as combined load
    if (mystranInstance->feaProblem.numAnalysis != 0) {

        for (i = 0; i < mystranInstance->feaProblem.numAnalysis; i++) {

            status = nastran_writeAnalysisCard(fp,
                                               &mystranInstance->feaProblem.feaAnalysis[i],
                                               &mystranInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            if (mystranInstance->feaProblem.feaAnalysis[i].numLoad != 0) {

                status = nastran_writeLoadADDCard(fp,
                                                  mystranInstance->feaProblem.numLoad+i+1,
                                                  mystranInstance->feaProblem.feaAnalysis[i].numLoad,
                                                  mystranInstance->feaProblem.feaAnalysis[i].loadSetID,
                                                  mystranInstance->feaProblem.feaLoad,
                                                  &mystranInstance->feaProblem.feaFileFormat);
                if (status != CAPS_SUCCESS) goto cleanup;

            } else { // If no loads for an individual analysis are specified assume that all loads should be applied

                if (mystranInstance->feaProblem.numLoad != 0) {

                    // Ignore thermal loads
                    k = 0;
                    for (j = 0; j < mystranInstance->feaProblem.numLoad; j++) {
                        if (mystranInstance->feaProblem.feaLoad[j].loadType == Thermal) continue;
                        k += 1;
                    }

                    if (k != 0) {
                        // Create a temporary list of load IDs
                        tempIntegerArray = (int *) EG_alloc(k*sizeof(int));
                        if (tempIntegerArray == NULL) {
                            status = EGADS_MALLOC;
                            goto cleanup;
                        }

                        k = 0;
                        for (j = 0; j < mystranInstance->feaProblem.numLoad; j++) {

                            if (mystranInstance->feaProblem.feaLoad[j].loadType == Thermal) continue;
                            tempIntegerArray[k] = mystranInstance->feaProblem.feaLoad[j].loadID;

                            k += 1;
                        }

                        // Write combined load card
                        status = nastran_writeLoadADDCard(fp,
                                                          mystranInstance->feaProblem.numLoad+i+1,
                                                          k,
                                                          tempIntegerArray,
                                                          mystranInstance->feaProblem.feaLoad,
                                                          &mystranInstance->feaProblem.feaFileFormat);
                        if (status != CAPS_SUCCESS) goto cleanup;

                        // Free temporary load ID list
                        if (tempIntegerArray != NULL) EG_free(tempIntegerArray);
                        tempIntegerArray = NULL;
                    }
                }
            }
        }

    } else { // If there aren't any analysis structures just write a single combined load card

        // Combined loads
        if (mystranInstance->feaProblem.numLoad != 0) {

            // Ignore thermal loads
            k = 0;
            for (j = 0; j < mystranInstance->feaProblem.numLoad; j++) {
                if (mystranInstance->feaProblem.feaLoad[j].loadType == Thermal) continue;
                k += 1;
            }

            if (k != 0) {
                // Create a temporary list of load IDs
                tempIntegerArray = (int *) EG_alloc(k*sizeof(int));
                if (tempIntegerArray == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                k = 0;
                for (j = 0; j < mystranInstance->feaProblem.numLoad; j++) {

                    if (mystranInstance->feaProblem.feaLoad[j].loadType == Thermal) continue;
                    tempIntegerArray[k] = mystranInstance->feaProblem.feaLoad[j].loadID;

                    k += 1;
                }

                // Write combined load card
                status = nastran_writeLoadADDCard(fp,
                                                  mystranInstance->feaProblem.numLoad+1,
                                                  k,
                                                  tempIntegerArray,
                                                  mystranInstance->feaProblem.feaLoad,
                                                  &mystranInstance->feaProblem.feaFileFormat);

                if (status != CAPS_SUCCESS) goto cleanup;

                // Free temporary load ID list
                if (tempIntegerArray != NULL) EG_free(tempIntegerArray);
                tempIntegerArray = NULL;
            }
        }
    }

    // Combined constraints
    if (mystranInstance->feaProblem.numConstraint != 0) {

        // Create a temporary list of constraint IDs
        tempIntegerArray = (int *) EG_alloc(mystranInstance->feaProblem.numConstraint*sizeof(int));
        if (tempIntegerArray == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        for (i = 0; i < mystranInstance->feaProblem.numConstraint; i++) {
            tempIntegerArray[i] = mystranInstance->feaProblem.feaConstraint[i].constraintID;
        }

        // Write combined constraint card
        status = nastran_writeConstraintADDCard(fp,
                                                mystranInstance->feaProblem.numConstraint+1,
                                                mystranInstance->feaProblem.numConstraint,
                                                tempIntegerArray,
                                                &mystranInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Free temporary constraint ID list
        if (tempIntegerArray != NULL) EG_free(tempIntegerArray);
        tempIntegerArray = NULL;
    }

    // Loads
    for (i = 0; i < mystranInstance->feaProblem.numLoad; i++) {
        status = nastran_writeLoadCard(fp,
                                       &mystranInstance->feaProblem.feaLoad[i],
                                       &mystranInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Constraints
    for (i = 0; i < mystranInstance->feaProblem.numConstraint; i++) {
        status = nastran_writeConstraintCard(fp,
                                             &mystranInstance->feaProblem.feaConstraint[i],
                                             &mystranInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Supports
    for (i = 0; i < mystranInstance->feaProblem.numSupport; i++) {
/*@-nullpass@*/
        status = nastran_writeSupportCard(fp,
                                          &mystranInstance->feaProblem.feaSupport[i],
                                          &mystranInstance->feaProblem.feaFileFormat,
                                          NULL);
/*@+nullpass@*/
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Materials
    for (i = 0; i < mystranInstance->feaProblem.numMaterial; i++) {
        status = nastran_writeMaterialCard(fp,
                                           &mystranInstance->feaProblem.feaMaterial[i],
                                           &mystranInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Properties
    for (i = 0; i < mystranInstance->feaProblem.numProperty; i++) {
        status = nastran_writePropertyCard(fp,
                                           &mystranInstance->feaProblem.feaProperty[i],
                                           &mystranInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Coordinate systems
    for (i = 0; i < mystranInstance->feaProblem.numCoordSystem; i++) {
        status = nastran_writeCoordinateSystemCard(fp,
                                                   &mystranInstance->feaProblem.feaCoordSystem[i],
                                                   &mystranInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Include mesh file
    fprintf(fp,"INCLUDE \'%s.bdf\'\n", mystranInstance->projectName);

    // End bulk data
    fprintf(fp,"ENDDATA\n");

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
      printf("Premature exit in mystranAIM preAnalysis status = %d\n", status);

    if (fp != NULL) fclose(fp);

    AIM_FREE(tempIntegerArray);

    return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteMYSTRAN AIM Execution
   *
   * If auto execution is enabled when creating an MYSTRAN AIM,
   * the AIM will execute MYSTRAN just-in-time with the command line:
   *
   * \code{.sh}
   * mystran $Proj_Name.dat > Info.out
   * \endcode
   *
   * where preAnalysis generated the file Proj_Name + ".dat" which contains the input information.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an MYSTRAN AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or MYSTRAN can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * mystran.preAnalysis()
   *
   * print ("\n\nRunning......")
   * mystran.system("mystran " + mystran.input.Proj_Name + ".dat > Info.out"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * mystran.postAnalysis()
   * \endcode
   */

  char command[PATH_MAX];
  aimStorage *mystranInstance;
  *state = 0;

  mystranInstance = (aimStorage *) instStore;
  if (mystranInstance == NULL) return CAPS_NULLVALUE;

  snprintf(command, PATH_MAX, "mystran %s.dat > Info.out",
           mystranInstance->projectName);

  return aim_system(aimInfo, NULL, command);
}


// ********************** AIM Function Break *****************************
// Check that MYSTRAN ran without errors
int aimPostAnalysis(void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int  status;
  char *filename = NULL; // File to open
  char extF06[] = ".F06";
  FILE *fp = NULL; // File pointer
  aimStorage *mystranInstance;

#ifdef DEBUG
  printf(" astrosAIM/aimPostAnalysis!\n");
#endif
  mystranInstance = (aimStorage *) instStore;

  // check that the mystran *.F06 file was created
  filename = (char *) EG_alloc((strlen(mystranInstance->projectName) +
                                strlen(extF06)+1)*sizeof(char));
  if (filename == NULL) return EGADS_MALLOC;

  sprintf(filename, "%s%s", mystranInstance->projectName, extF06);

  fp = aim_fopen(aimInfo, filename, "r");
  if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open Output file: %s!", filename);
      status = CAPS_IOERR;
      goto cleanup;
  }
  status = CAPS_SUCCESS;

cleanup:
    EG_free(filename); // Free filename allocation
    if (fp != NULL) fclose(fp);

    return status;
}


// Set MYSTRAN output variables
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsMYSTRAN AIM Outputs
     * The following list outlines the MYSTRAN outputs available through the AIM interface.
     */

#ifdef DEBUG
    printf(" mystranAIM/aimOutputs index = %d!\n", index);
#endif

    /*! \page aimOutputsMYSTRAN AIM Outputs
     * - <B>EigenValue</B> = List of Eigen-Values (\f$ \lambda\f$) after a modal solve.
     * - <B>EigenRadian</B> = List of Eigen-Values in terms of radians (\f$ \omega = \sqrt{\lambda}\f$ ) after a modal solve.
     * - <B>EigenFrequency</B> = List of Eigen-Values in terms of frequencies (\f$ f = \frac{\omega}{2\pi}\f$) after a modal solve.
     * - <B>EigenGeneralMass</B> = List of generalized masses for the Eigen-Values.
     * .
     */

    if (index == 1) {
        *aoname = EG_strdup("EigenValue");

    } else if (index == 2) {
        *aoname = EG_strdup("EigenRadian");

    } else if (index == 3) {
        *aoname = EG_strdup("EigenFrequency");

    } else if (index == 4) {
        *aoname = EG_strdup("EigenGeneralMass");
    }

    form->type       = Double;
    form->units      = NULL;
    form->lfixed     = Change;
    form->sfixed     = Change;
    form->vals.reals = NULL;
    form->vals.real  = 0;

    /*else if (index == 3) {
     *aoname = EG_strdup("EigenVector");
    	form->type          = Double;
    	form->units         = NULL;
    	form->lfixed        = Change;
    	form->sfixed        = Change;
    }*/

    return CAPS_SUCCESS;
}


// Calculate MYSTRAN output
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int status = CAPS_SUCCESS; // Function return status

    int i; // Indexing

    char *filename = NULL; // File to open
    char extOP1[] = ".OP1";
    char extF06[] = ".F06";
    FILE *fp = NULL; // File pointer
    aimStorage *mystranInstance;

    const double pi = 4.0 * atan(1.0);
  
    mystranInstance = (aimStorage *) instStore;

    if (index == 1 || index == 2 || index == 3) {
        // Open mystran *.OP1 file - OUTPUT4
        filename = (char *) EG_alloc((strlen(mystranInstance->projectName) +
                                      strlen(extOP1)+1)*sizeof(char));
        if (filename == NULL) return EGADS_MALLOC;

        sprintf(filename, "%s%s", mystranInstance->projectName, extOP1);

        fp = aim_fopen(aimInfo, filename, "rb");

        EG_free(filename); // Free filename allocation

        if (fp == NULL) {
#ifdef DEBUG
            printf(" mystranAIM/aimCalcOutput Cannot open Output file!\n");
#endif
            return CAPS_IOERR;
        }

        status = mystran_readOutput4Data(fp, "EIGEN_VA", val);
        if (status == CAPS_SUCCESS) {
            if (index == 2) {
                for (i = 0; i < val->length; i++) {
                    val->vals.reals[i] = sqrt(val->vals.reals[i]);
                }
                //val->units = EG_strdup("rads/s");
            }

            if (index == 3) {
                for (i = 0; i < val->length; i++) {
                    val->vals.reals[i] = sqrt(val->vals.reals[i])/(2*pi);
                }
                //val->units = EG_strdup("Hz");
            }
        }

    } else if (index == 4) {
        // Open mystran *.OP1 file - OUTPUT4
        filename = (char *) EG_alloc((strlen(mystranInstance->projectName) +
                                      strlen(extOP1) +1)*sizeof(char));
        if (filename == NULL) return EGADS_MALLOC;

        sprintf(filename, "%s%s", mystranInstance->projectName, extOP1);

        fp = aim_fopen(aimInfo, filename, "rb");

        EG_free(filename); // Free filename allocation

        if (fp == NULL) {
#ifdef DEBUG
            printf(" mystranAIM/aimCalcOutput Cannot open Output file!\n");
#endif
            return CAPS_IOERR;
        }

        status = mystran_readOutput4Data(fp, "GEN_MASS", val);

    } else if (index == 5) {

        // Open mystran *.F06 file
        filename = (char *) EG_alloc((strlen(mystranInstance->projectName) +
                                      strlen(extF06) +1)*sizeof(char));
        if (filename == NULL) return EGADS_MALLOC;

        sprintf(filename, "%s%s", mystranInstance->projectName, extF06);

        fp = aim_fopen(aimInfo, filename, "r");

        EG_free(filename); // Free filename allocation
        if (fp == NULL) {
#ifdef DEBUG
            printf(" mystranAIM/aimCalcOutput Cannot open Output file!\n");
#endif
            return CAPS_IOERR;
        }
/*
        double **dataMatrix = NULL;
        int numEigenVector = 0;
        int numGridPoint = 0;
        status = mystran_readF06Eigen(fp, &numEigenVector, &numGridPoint,
                                      &dataMatrix);
*/
    }

    if (fp != NULL) fclose(fp);

    return status;
}


void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *mystranInstance;

#ifdef DEBUG
    printf(" mystranAIM/aimCleanup!\n");
#endif
    mystranInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(mystranInstance);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during clean up of instance\n", status);

    EG_free(mystranInstance);
}


int aimDiscr(char *tname, capsDiscr *discr)
{
    int status; // Function return status

    aimStorage *mystranInstance;

    int i, numBody;

    // EGADS objects
    ego *bodies = NULL, *tess = NULL;

    const char   *intents;

#ifdef DEBUG
    printf(" mystranAIM/aimDiscr: tname = %s!\n", tname);
#endif

    if (tname == NULL) return CAPS_NOTFOUND;
    mystranInstance = (aimStorage *) discr->instStore;

/*  if (mystranInstance->dataTransferCheck == (int) false) {
        printf("The volume is not suitable for data transfer - possibly the volume mesher "
               "added unaccounted for points\n");
        return CAPS_BADVALUE;
    }*/

    // Currently this ONLY works if the capsTranfer lives on single body!
    status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
    if ((status != CAPS_SUCCESS) || (bodies == NULL)) {
        if (status == CAPS_SUCCESS) status = CAPS_NULLOBJ;
        printf(" mystranAIM/aimDiscr: aim_getBodies = %d!\n", status);
        return status;
    }

    // Check and generate/retrieve the mesh
    status = checkAndCreateMesh(discr->aInfo, mystranInstance);
    AIM_STATUS(discr->aInfo, status);

    AIM_ALLOC(tess, mystranInstance->numMesh, ego, discr->aInfo, status);
    for (i = 0; i < mystranInstance->numMesh; i++) {
        tess[i] = mystranInstance->feaMesh[i].egadsTess;
    }

    status = mesh_fillDiscr(tname, &mystranInstance->attrMap, mystranInstance->numMesh, tess, discr);
    AIM_STATUS(discr->aInfo, status);

#ifdef DEBUG
    printf(" mystranAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
#endif

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit: function aimDiscr mystranAIM status = %d\n",
               status);

    AIM_FREE(tess);
    return status;
}


int aimTransfer(capsDiscr *discr, const char *dataName, int numPoint,
                int dataRank, double *dataVal, /*@unused@*/ char **units)
{

    /*! \page dataTransferMYSTRAN MYSTRAN Data Transfer
     *
     * The MYSTRAN AIM has the ability to transfer displacements and eigenvectors from the AIM and pressure
     * distributions to the AIM using the conservative and interpolative data transfer schemes in CAPS.
     *
     * \section dataFromMYSTRAN Data transfer from MYSTRAN  (FieldOut)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements from the *.F06 file.
     * </ul>
     *
     * <ul>
     *  <li> <B>"EigenVector_#"</B> </li> <br>
     *   Retrieves modal eigen-vectors from the *.F06 file, where "#" should be replaced by the
     *   corresponding mode number for the eigen-vector (eg. EigenVector_3 would correspond to the third mode,
     *   while EigenVector_6 would be the sixth mode).
     * </ul>
     *
     * \section dataToMYSTRAN Data transfer to MYSTRAN  (FieldIn)
     * <ul>
     *  <li> <B>"Pressure"</B> </li> <br>
     *  Writes appropriate load cards using the provided pressure distribution.
     * </ul>
     *
     */

    int status; // Function return status
    int i, j, dataPoint, bIndex; // Indexing
    aimStorage *mystranInstance;

    char *extF06 = ".F06";

    // FO6 data variables
    int numGridPoint = 0;
    int numEigenVector = 0;

    double **dataMatrix = NULL;

    // Specific EigenVector to use
    int eigenVectorIndex = 0;

    // Variables used in global node mapping
    int globalNodeID;

    // Filename stuff
    char *filename = NULL;
    FILE *fp; // File pointer

#ifdef DEBUG
    printf(" mystranAIM/aimTransfer name = %s  npts = %d/%d!\n",
           dataName, numPoint, dataRank);
#endif
    mystranInstance = (aimStorage *) discr->instStore;

    //capsGroupList = &storage[0]; // List of boundary ID (attrMap) in the transfer

    if (strcasecmp(dataName, "Displacement")    != 0 &&
           strncmp(dataName, "EigenVector", 11) != 0) {

        printf("Unrecognized data transfer variable - %s\n", dataName);
        return CAPS_NOTFOUND;
    }

    filename = (char *) EG_alloc((strlen(mystranInstance->projectName) +
                                  strlen(extF06) + 1)*sizeof(char));
    if (filename == NULL) return EGADS_MALLOC;

    sprintf(filename,"%s%s", mystranInstance->projectName, extF06);

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

            status = mystran_readF06Displacement(fp,
                                                 1,
                                                 &numGridPoint,
                                                 &dataMatrix);
            fclose(fp);
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

            status = mystran_readF06EigenVector(fp,
                                                &numEigenVector,
                                                &numGridPoint,
                                                &dataMatrix);
        }

        fclose(fp);

    } else {

        status = CAPS_NOTFOUND;
    }

    if (status != CAPS_SUCCESS) return status;

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
    EG_free(ptr); // Extra information to store into the discr void pointer - just a int array
}


int aimLocateElement(capsDiscr *discr,  double *params, double *param,
                     int *bIndex, int *eIndex, double *bary)
{
#ifdef DEBUG
    printf(" mystranAIM/aimLocateElement!\n");
#endif

    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}


int aimInterpolation(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                     double *bary, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" mystranAIM/aimInterpolation  %s!\n", name);
#endif

    return aim_interpolation(discr, name, bIndex, eIndex, bary, rank, data,
                             result);

}


int aimInterpolateBar(capsDiscr *discr,  const char *name, int bIndex,
                      int eIndex,  double *bary, int rank,
                      double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" mystranAIM/aimInterpolateBar  %s!\n", name);
#endif

    return aim_interpolateBar(discr, name, bIndex, eIndex, bary, rank, r_bar,
                              d_bar);
}


int aimIntegration(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                   int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" mystranAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar(capsDiscr *discr, const char *name, int bIndex, int eIndex,
                    int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" mystranAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank, r_bar, d_bar);
}
