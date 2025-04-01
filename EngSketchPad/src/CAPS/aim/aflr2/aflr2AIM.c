/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AFLR2 AIM
 *
 *      Modified from code written by Dr. Ryan Durscher AFRL/RQVC
 *
 */

/*!\mainpage Introduction
 *
 * \section overviewAFLR2 AFLR2 AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact with the
 * unstructured, surface grid generator AFLR2 \cite Marcum1995 \cite Marcum1998.
 *
 * The AFLR2 AIM provides the CAPS users with the ability to generate "unstructured, 2D grids" using an
 * "Advancing-Front/Local-Reconnection (AFLR) procedure." Both triangular and quadrilateral elements may be generated.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsAFLR2 and
 * \ref aimOutputsAFLR2 and \ref attributeAFLR2, respectively.
 *
 * The complete AFLR documentation is available at the <a href="https://www.simcenter.msstate.edu/software/documentation/system/index.html">SimCenter</a>.
 *
 */

/*! \page attributeAFLR2 AIM attributes
 * The following list of attributes are required for the AFLR2 AIM inside the geometry input.
 *
 * - <b> capsGroup</b> This is a name assigned to any geometric entity to denote a "boundary" for further analysis. Recall that a
 *  string in ESP starts with a $.  For example, attribute <c>capsGroup $Wing</c>.
 *
 * - <b> capsMesh</b> This is a name assigned to any geometric entity in order to control meshing related parameters. Recall that a
 *  string in ESP starts with a $.  For example, attribute <c>capsMesh $Wing</c>.
 *
 */

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "aimMesh.h"

#include "meshUtils.h"       // Collection of helper functions for meshing
#include "miscUtils.h"
#include "deprecateUtils.h"

#include "aflr2_Interface.h" // Bring in AFLR2 'interface' functions

#ifdef WIN32
#define strcasecmp stricmp
#endif

#define MXCHAR  255

//#define DEBUG


enum aimInputs
{
  Proj_Name = 1,               /* index is 1-based */
  Tess_Params,
  Mesh_Quiet_Flag,
  Mesh_Format,
  Mesh_Gen_Input_String,
  Edge_Point_Min,
  Edge_Point_Max,
  Mesh_Sizing,
  NUMINPUT = Mesh_Sizing       /* Total number of inputs */
};

enum aimOutputs
{
  NumberOfElement = 1,      /* index is 1-based */
  NumberOfNode,
  Area_Mesh,
  NUMOUT = Area_Mesh        /* Total number of outputs */
};


typedef struct {

    // Container for mesh input
    meshInputStruct meshInput;

    // Attribute to index map for capsGroup
    mapAttrToIndexStruct groupMap;

    // Attribute to index map for capsMesh
    mapAttrToIndexStruct meshMap;

    aimMeshRef meshRef;

} aimStorage;


static int destroy_aimStorage(aimStorage *aflr2Instance, int inUpdate)
{
    int status; // Function return status

    // Free the meshRef
    aim_freeMeshRef(&aflr2Instance->meshRef);
    aflr2Instance->meshRef.type = aimAreaMesh;

    // Destroy meshInput
    status = destroy_meshInputStruct(&aflr2Instance->meshInput);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, aflr2AIM meshInput cleanup!!!\n", status);

    if (inUpdate == (int)true) return CAPS_SUCCESS;

    // Destroy attribute to index map
    status = destroy_mapAttrToIndexStruct(&aflr2Instance->groupMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, aflr2AIM attributeMap cleanup!!!\n", status);

    status = destroy_mapAttrToIndexStruct(&aflr2Instance->meshMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, aflr2AIM attributeMap cleanup!!!\n", status);

    return CAPS_SUCCESS;
}


#ifdef WIN32
    #if  _MSC_VER >= 1900
    // FILE _iob[] = {*stdin, *stdout, *stderr};
    FILE _iob[3];

    extern FILE * __cdecl __iob_func(void)
    {
        return _iob;
    }
    #endif
#endif


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    aimStorage *aflr2Instance=NULL;

    #ifdef WIN32
        #if  _MSC_VER >= 1900
            _iob[0] = *stdin;
            _iob[1] = *stdout;
            _iob[2] = *stderr;
        #endif
    #endif

    int status; // Function return

#ifdef DEBUG
    printf("\n aflr2AIM/aimInitialize   inst = %d!\n", inst);
#endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames  = NULL;
    *franks  = NULL;
    *fInOut  = NULL;

    // Allocate aflrInstance
    AIM_ALLOC(aflr2Instance, 1, aimStorage, aimInfo, status);
    *instStore = aflr2Instance;

    // Set initial values for aflrInstance //

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&aflr2Instance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&aflr2Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    // Container for mesh input
    status = initiate_meshInputStruct(&aflr2Instance->meshInput);
    AIM_STATUS(aimInfo, status);

    // Mesh reference passed to solver
    status = aim_initMeshRef(&aflr2Instance->meshRef, aimAreaMesh);
    AIM_STATUS(aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);
    return CAPS_SUCCESS;
}


int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{

    /*! \page aimInputsAFLR2 AIM Inputs
     * The following list outlines the AFLR2 meshing options along with their default value available
     * through the AIM interface.
     */

#ifdef DEBUG
    printf(" aflr2AIM/aimInputs instance = %d  index = %d!\n", index);
#endif

    // Meshing Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name"); // If NULL a volume grid won't be written by the AIM
        defval->type         = String;
        defval->nullVal      = NotAllowed;
        defval->vals.string  = EG_strdup("aflr2_CAPS");

        /*! \page aimInputsAFLR2
         * - <B> Proj_Name = "aflr2_CAPS"</B> <br>
         * Output name prefix for meshes to be written in formats specified by Mesh_Format.
         * These meshes are not linked to any analysis, but may be useful exploring meshing parameters.
         */

    } else if (index == Mesh_Format) {
        *ainame               = EG_strdup("Mesh_Format");
        defval->type          = String;
        defval->vals.string   = NULL;
        defval->nullVal       = IsNull;
        defval->dim           = Vector;
        defval->lfixed        = Change;

        /*! \page aimInputsAFLR2
         * \include{doc} Mesh_Format.dox
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

        /*! \page aimInputsAFLR2
         * - <B> Tess_Params = [0.025, 0.001, 15.0]</B> <br>
         * Body tessellation parameters. Tess_Params[0] and Tess_Params[1] get scaled by the capsMeshLength
         * attribute. (From the EGADS manual) A set of 3 parameters that drive the EDGE discretization
         * and the FACE triangulation. The first is the maximum length of an EDGE segment or triangle side
         * (in physical space). A zero is flag that allows for any length. The second is a curvature-based
         * value that looks locally at the deviation between the centroid of the discrete object and the
         * underlying geometry. Any deviation larger than the input value will cause the tessellation to
         * be enhanced in those regions. The third is the maximum interior dihedral angle (in degrees)
         * between triangle facets (or Edge segment tangents for a WIREBODY tessellation), note that a
         * zero ignores this phase.
         */

    } else if (index == Mesh_Quiet_Flag) {
        *ainame               = EG_strdup("Mesh_Quiet_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsAFLR2
         * - <B> Mesh_Quiet_Flag = False</B> <br>
         * Complete suppression of mesh generator (not including errors)
         */

    } else if (index == Mesh_Gen_Input_String) {
        *ainame               = EG_strdup("Mesh_Gen_Input_String");
        defval->type          = String;
        defval->nullVal       = IsNull;
        defval->vals.string   = NULL;

        /*! \page aimInputsAFLR2
         * - <B> Mesh_Gen_Input_String = NULL</B> <br>
         * Meshing program command line string (as if called in bash mode). Use this to specify more
         * complicated options/use features of the mesher not currently exposed through other AIM input
         * variables. Note that this is the exact string that will be provided to the volume mesher; no
         * modifications will be made. If left NULL an input string will be created based on default values
         * of the relevant AIM input variables.
         */

    } else if (index == Edge_Point_Min) {
        *ainame               = EG_strdup("Edge_Point_Min");
        defval->type          = Integer;
        defval->vals.integer  = 0;
        defval->lfixed        = Fixed;
        defval->nrow          = 1;
        defval->ncol          = 1;
        defval->nullVal       = IsNull;

        /*! \page aimInputsAFLR2
         * - <B> Edge_Point_Min = NULL</B> <br>
         * Minimum number of points on an edge including end points to use when creating a surface mesh (min 2).
         */

    } else if (index == Edge_Point_Max) {
        *ainame               = EG_strdup("Edge_Point_Max");
        defval->type          = Integer;
        defval->vals.integer  = 0;
        defval->lfixed        = Fixed;
        defval->nrow          = 1;
        defval->ncol          = 1;
        defval->nullVal       = IsNull;

        /*! \page aimInputsAFLR2
         * - <B> Edge_Point_Max = NULL</B> <br>
         * Maximum number of points on an edge including end points to use when creating a surface mesh (min 2).
         */

    } else if (index == Mesh_Sizing) {
        *ainame              = EG_strdup("Mesh_Sizing");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAFLR2
         * - <B>Mesh_Sizing = NULL </B> <br>
         * See \ref meshSizingProp for additional details.
         */
    }

    return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
    int status; // Function return status

    int i;
    char aimFile[PATH_MAX];

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;

    // Global settings
    int minEdgePoint = -1, maxEdgePoint = -1;
    double refLen = -1.0, capsMeshLength;

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    aimStorage *aflr2Instance;

    aflr2Instance = (aimStorage *) instStore;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if (numBody <= 0 || bodies == NULL) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }
    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Cleanup previous aimStorage for the instance in case this is the
    //         second time through aimUpdateState for the same instance
    status = destroy_aimStorage(aflr2Instance, (int)true);
    AIM_STATUS(aimInfo, status);

    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ) {
      // Get capsMesh name and index mapping
      status = create_CAPSMeshAttrToIndexMap(numBody,
                                             bodies,
                                             3, // Node level
                                             &aflr2Instance->meshMap);
      AIM_STATUS(aimInfo, status);
    }

    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ) {
      // Get capsGroup name and index mapping to make sure all edges have a capsGroup value
      status = create_CAPSGroupAttrToIndexMap(numBody,
                                              bodies,
                                              -2, // Edges only
                                              &aflr2Instance->groupMap);
      AIM_STATUS(aimInfo, status);

      // offset by 1 to account for 'volume' group
      for (i = 0; i < aflr2Instance->groupMap.numAttribute; i++) {
        aflr2Instance->groupMap.attributeIndex[i]++;
      }
    }

    // Setup meshing input structure

    status = check_CAPSMeshLength(numBody, bodies, &capsMeshLength);
    if (capsMeshLength <= 0 || status != CAPS_SUCCESS) {
      if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "capsMeshLength is not set on any body.");
      } else {
        AIM_ERROR(aimInfo, "capsMeshLength: %f", capsMeshLength);
      }
      AIM_ADDLINE(aimInfo, "The capsMeshLength attribute must\n"
                           "present on at least one body.\n"
                           "\n"
                           "capsMeshLength should be a a positive value representative\n"
                           "of a characteristic length of the geometry,\n"
                           "e.g. the MAC of a wing or diameter of a fuselage.\n");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Get Tessellation parameters -Tess_Params
    aflr2Instance->meshInput.paramTess[0] = aimInputs[Tess_Params-1].vals.reals[0]*capsMeshLength;
    aflr2Instance->meshInput.paramTess[1] = aimInputs[Tess_Params-1].vals.reals[1]*capsMeshLength;
    aflr2Instance->meshInput.paramTess[2] = aimInputs[Tess_Params-1].vals.reals[2];

    aflr2Instance->meshInput.quiet           = aimInputs[Mesh_Quiet_Flag-1].vals.integer;

    // Set aflr2 specific mesh inputs
    if (aimInputs[Mesh_Gen_Input_String-1].nullVal != IsNull) {
        AIM_STRDUP(aflr2Instance->meshInput.aflr4Input.meshInputString,
                   aimInputs[Mesh_Gen_Input_String-1].vals.string, aimInfo, status);
    }

    // Max and min number of points
    if (aimInputs[Edge_Point_Min-1].nullVal != IsNull) {
        minEdgePoint = aimInputs[Edge_Point_Min-1].vals.integer;
        if (minEdgePoint < 2) {
          AIM_ERROR(aimInfo, "Edge_Point_Min = %d must be greater or equal to 2\n", minEdgePoint);
          status = CAPS_BADVALUE;
          goto cleanup;
        }
    }

    if (aimInputs[Edge_Point_Max-1].nullVal != IsNull) {
        maxEdgePoint = aimInputs[Edge_Point_Max-1].vals.integer;
        if (maxEdgePoint < 2) {
          AIM_ERROR(aimInfo, "Edge_Point_Max = %d must be greater or equal to 2\n", maxEdgePoint);
          status = CAPS_BADVALUE;
          goto cleanup;
        }
    }

    if (maxEdgePoint >= 2 && minEdgePoint >= 2 && minEdgePoint > maxEdgePoint) {
      AIM_ERROR  (aimInfo, "Edge_Point_Max must be greater or equal Edge_Point_Min\n");
      AIM_ADDLINE(aimInfo, "Edge_Point_Max = %d, Edge_Point_Min = %d\n",maxEdgePoint,minEdgePoint);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Get mesh sizing parameters
    if (aimInputs[Mesh_Sizing-1].nullVal != IsNull &&
        aim_newAnalysisIn(aimInfo, Mesh_Sizing) == CAPS_SUCCESS) {

        status = deprecate_SizingAttr(aimInfo,
                                      aimInputs[Mesh_Sizing-1].length,
                                      aimInputs[Mesh_Sizing-1].vals.tuple,
                                      &aflr2Instance->meshMap,
                                      &aflr2Instance->groupMap);
        AIM_STATUS(aimInfo, status);

        status = mesh_getSizingProp(aimInfo,
                                    aimInputs[Mesh_Sizing-1].length,
                                    aimInputs[Mesh_Sizing-1].vals.tuple,
                                    &aflr2Instance->meshMap,
                                    &numMeshProp,
                                    &meshProp);
        AIM_STATUS(aimInfo, status);
    }

    // Modify the EGADS body tessellation based on given inputs
    status =  mesh_modifyBodyTess(aimInfo,
                                  numMeshProp,
                                  meshProp,
                                  minEdgePoint,
                                  maxEdgePoint,
                                  (int) false, // quadMesh
                                  &refLen,
                                  aflr2Instance->meshInput.paramTess,
                                  aflr2Instance->meshMap,
                                  numBody,
                                  bodies);
    AIM_STATUS(aimInfo, status);


    // set the filename without extensions where the grid is written for solvers
    if (aimInputs[Proj_Name-1].nullVal != IsNull) {
      status = aim_file(aimInfo, aimInputs[Proj_Name-1].vals.string, aimFile);
      AIM_STATUS(aimInfo, status);
    } else {
      status = aim_file(aimInfo, AFLR2FILE, aimFile);
      AIM_STATUS(aimInfo, status);
    }
    AIM_STRDUP(aflr2Instance->meshRef.fileName, aimFile, aimInfo, status);

    AIM_ALLOC(aflr2Instance->meshRef.bnds, aflr2Instance->groupMap.numAttribute, aimMeshBnd, aimInfo, status);
    aflr2Instance->meshRef.nbnd = aflr2Instance->groupMap.numAttribute;
    for (i = 0; i < aflr2Instance->meshRef.nbnd; i++) {
      status = aim_initMeshBnd(aflr2Instance->meshRef.bnds + i);
      AIM_STATUS(aimInfo, status);
    }

    for (i = 0; i < aflr2Instance->meshRef.nbnd; i++) {
      AIM_STRDUP(aflr2Instance->meshRef.bnds[i].groupName, aflr2Instance->groupMap.attributeName[i], aimInfo, status);
      aflr2Instance->meshRef.bnds[i].ID = aflr2Instance->groupMap.attributeIndex[i];
    }

    status = CAPS_SUCCESS;
cleanup:

    if (meshProp != NULL) {
        for (i = 0; i < numMeshProp; i++) {
            (void) destroy_meshSizingStruct(&meshProp[i]);
        }
        AIM_FREE(meshProp);
    }

    return status;
}

// ********************** AIM Function Break *****************************
int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{
    int status; // Status return

    int i, bodyIndex; // Indexing

    int Message_Flag;

    const char *cProj_Name=NULL;

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    // Container for surface mesh
    int numSurface=0;
    meshStruct *surfaceMesh=NULL;

    // File output
    char *filename = NULL;
    //char bodyNumber[11];

    const aimStorage *aflr2Instance;

    aflr2Instance = (const aimStorage *) instStore;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if (numBody <= 0 || bodies == NULL) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }
    AIM_NOTNULL(aimInputs, aimInfo, status);

    // remove previous meshes
    status = aim_deleteMeshes(aimInfo, &aflr2Instance->meshRef);
    AIM_STATUS(aimInfo, status);

    // Allocate surfaceMesh from number of bodies
    numSurface = numBody;
    AIM_ALLOC(surfaceMesh, numSurface, meshStruct, aimInfo, status);

    // Initiate surface meshes
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++){
        status = initiate_meshStruct(&surfaceMesh[bodyIndex]);
        AIM_STATUS(aimInfo, status);

        status = copy_mapAttrToIndexStruct(&aflr2Instance->groupMap, &surfaceMesh[bodyIndex].groupMap);
        AIM_STATUS(aimInfo, status);
    }

    if (aimInputs[Proj_Name-1].nullVal != IsNull) {
      cProj_Name = aimInputs[Proj_Name-1].vals.string;
    } else {
      cProj_Name = AFLR2FILE;
    }

    // Run AFLR2 for each body
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {

        printf("Getting 2D mesh for body %d (of %d)\n", bodyIndex+1, numBody);

        Message_Flag = aflr2Instance->meshInput.quiet == (int)true ? 0 : 1;
        status = aflr2_Surface_Mesh(aimInfo, Message_Flag, bodies[bodyIndex],
                                    &aflr2Instance->meshInput,
                                    &aflr2Instance->groupMap,
                                    &aflr2Instance->meshMap,
                                    cProj_Name,
                                    &surfaceMesh[bodyIndex],
                                    &aflr2Instance->meshRef);
        if (status != CAPS_SUCCESS) {
            AIM_ERROR(aimInfo, "Problem during meshing of body %d\n", bodyIndex+1);
            goto cleanup;
        }

        printf("Number of nodes = %d\n", surfaceMesh[bodyIndex].numNode);
        printf("Number of elements = %d\n", surfaceMesh[bodyIndex].numElement);
        if (surfaceMesh[bodyIndex].meshQuickRef.useStartIndex == (int) true ||
            surfaceMesh[bodyIndex].meshQuickRef.useListIndex == (int) true) {
            printf("Number of tris = %d\n", surfaceMesh[bodyIndex].meshQuickRef.numTriangle);
            printf("Number of quad = %d\n", surfaceMesh[bodyIndex].meshQuickRef.numQuadrilateral);
        }

    }

    status = CAPS_SUCCESS;

cleanup:

    // Destroy surface mesh allocated arrays
    for (i = 0; i < numSurface; i++) {
        (void) destroy_meshStruct(&surfaceMesh[i]);
    }
    AIM_FREE(surfaceMesh);

    AIM_FREE(filename);
    return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ const void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS;
  int i;

  // Body parameters
  const char *intents;
  int numBody = 0; // Number of bodies
  ego *bodies = NULL; // EGADS body objects

  ego body, tess=NULL;
  int state, nglobal;
  int iglobal;
  int localIndex, topoIndex;

  aimStorage *aflr2Instance;
  aimMesh    mesh;
  char aimFile[PATH_MAX];
  double params[3];

  aflr2Instance = (aimStorage *) instStore;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(bodies, aimInfo, status);

  AIM_ALLOC(aflr2Instance->meshRef.maps, 1, aimMeshTessMap, aimInfo, status);
  aflr2Instance->meshRef.nmap = 1;
  aflr2Instance->meshRef.maps[0].map = NULL;
  aflr2Instance->meshRef.maps[0].tess = NULL;

  status = aim_file(aimInfo, AFLR2TESSFILE, aimFile);
  AIM_STATUS(aimInfo, status);

  // read the eto file
  status = EG_loadTess(bodies[0], aimFile, &aflr2Instance->meshRef.maps[0].tess);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(aflr2Instance->meshRef.maps[0].tess, aimInfo, status);

  // make a copy of the tessellation with all the triangulation
  status = EG_copyObject(aflr2Instance->meshRef.maps[0].tess, NULL, &tess);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(tess, aimInfo, status);

  status = EG_openTessBody(aflr2Instance->meshRef.maps[0].tess);
  AIM_STATUS(aimInfo, status);

  // Negating the first parameter triggers EGADS to only put vertexes on edges
  // This removes the face tessellation
  params[0] = -fabs(aflr2Instance->meshInput.paramTess[0]);
  params[1] =       aflr2Instance->meshInput.paramTess[1];
  params[2] =       aflr2Instance->meshInput.paramTess[2];

  status = EG_finishTess( aflr2Instance->meshRef.maps[0].tess, params );
  AIM_STATUS(aimInfo, status);

  // get the body from the input tessellation
  status = EG_statusTessBody(aflr2Instance->meshRef.maps[0].tess, &body, &state, &nglobal);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(body, aimInfo, status);

  // allocate the boundary to area mesh mapping
  AIM_ALLOC(aflr2Instance->meshRef.maps[0].map, nglobal, int, aimInfo, status);

  // Find the boundary mesh in the global tessellation
  for (i = 0; i < nglobal; i++) {

    // Get the local indexes from the boundary mesh
    status = EG_getGlobal(aflr2Instance->meshRef.maps[0].tess, i+1,
                          &localIndex, &topoIndex, NULL);
    AIM_STATUS(aimInfo, status);

    // Get the global index in the full 2D mesh
    if (localIndex == 0) {
      status = EG_localToGlobal(tess, 0, topoIndex, &iglobal);
      AIM_STATUS(aimInfo, status);
    } else if (topoIndex > 0) {
      status = EG_localToGlobal(tess, -topoIndex, localIndex, &iglobal);
      AIM_STATUS(aimInfo, status);
    } else {
      AIM_ERROR(aimInfo, "Developer exception! Should not find Face index!");
      status = CAPS_NOTIMPLEMENT;
      goto cleanup;
    }


    aflr2Instance->meshRef.maps[0].map[i] = iglobal;
  }

  status = aim_queryMeshes( aimInfo, Mesh_Format, ANALYSISIN, &aflr2Instance->meshRef );
  if (status > 0) {

/*@-immediatetrans@*/
      mesh.meshData = NULL;
      mesh.meshRef = &aflr2Instance->meshRef;
/*@+immediatetrans@*/

      status = aim_readBinaryUgrid(aimInfo, &mesh);
      AIM_STATUS(aimInfo, status);

      status = aim_writeMeshes(aimInfo, Mesh_Format, ANALYSISIN, &mesh);
      AIM_STATUS(aimInfo, status);

      status = aim_freeMeshData(mesh.meshData);
      AIM_STATUS(aimInfo, status);
      AIM_FREE(mesh.meshData);

  }
  else
      AIM_STATUS(aimInfo, status);


  status = CAPS_SUCCESS;
cleanup:

  EG_deleteObject(tess);

  return status;
}


// ********************** AIM Function Break *****************************
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               /*@unused@*/ int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsAFLR2 AIM Outputs
     * The following list outlines the AFLR2 AIM outputs available through the AIM interface.
     */

    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" aflr2AIM/aimOutputs index = %d!\n", index);
#endif

    if (index == NumberOfElement) {
        *aoname = EG_strdup("NumberOfElement");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsAFLR2
         * - <B> NumberOfElement </B> <br>
         * Number of elements in the surface mesh
         */

    } else if (index == NumberOfNode) {
        *aoname = EG_strdup("NumberOfNode");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsAFLR2
         * - <B> NumberOfNode </B> <br>
         * Number of vertices in the surface mesh
         */

    } else if (index == Area_Mesh) {
        *aoname           = AIM_NAME(Area_Mesh);
        form->type        = PointerMesh;
        form->dim         = Scalar;
        form->lfixed      = Fixed;
        form->sfixed      = Fixed;
        form->vals.AIMptr = NULL;
        form->nullVal     = IsNull;

        /*! \page aimOutputsAFLR2
         * - <B> Area_Mesh </B> <br>
         * The resulting mesh that can be linked to an anlaysis input.
         */
    }

    AIM_NOTNULL(*aoname, aimStruc, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*aoname);
    return status;
}


// ********************** AIM Function Break *****************************
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ int index, capsValue *val)
{
    int        status = CAPS_SUCCESS;
    int        nVertex, nTri, nQuad, nTet, nPyramid, nPrism, nHex;
    aimStorage *aflr2Instance;
    aimMesh    mesh;

#ifdef DEBUG
    printf(" aflr2AIM/aimCalcOutput index = %d!\n",  index);
#endif
    aflr2Instance = (aimStorage *) instStore;

    if (NumberOfElement == index) {

        // Count the total number of surface elements
        status = aim_readBinaryUgridHeader(aimInfo, &aflr2Instance->meshRef,
                                           &nVertex, &nTri, &nQuad,
                                           &nTet, &nPyramid, &nPrism, &nHex);
        AIM_STATUS(aimInfo, status);

        val->vals.integer = nTri + nQuad;

    } else if (NumberOfNode == index) {

        // Count the total number of surface vertices
        status = aim_readBinaryUgridHeader(aimInfo, &aflr2Instance->meshRef,
                                           &nVertex, &nTri, &nQuad,
                                           &nTet, &nPyramid, &nPrism, &nHex);
        AIM_STATUS(aimInfo, status);

        val->vals.integer = nVertex;

    } else if (Area_Mesh == index) {

        status = aim_queryMeshes( aimInfo, Area_Mesh, ANALYSISOUT, &aflr2Instance->meshRef );
        if (status > 0) {

/*@-immediatetrans@*/
            mesh.meshData = NULL;
            mesh.meshRef = &aflr2Instance->meshRef;
/*@+immediatetrans@*/

            status = aim_readBinaryUgrid(aimInfo, &mesh);
            AIM_STATUS(aimInfo, status);

            status = aim_writeMeshes(aimInfo, Area_Mesh, ANALYSISOUT, &mesh);
            AIM_STATUS(aimInfo, status);

            status = aim_freeMeshData(mesh.meshData);
            AIM_STATUS(aimInfo, status);
            AIM_FREE(mesh.meshData);

        }
        else
            AIM_STATUS(aimInfo, status);

/*@-immediatetrans@*/
        // Return the area mesh reference
        val->nrow        = 1;
        val->vals.AIMptr = &aflr2Instance->meshRef;
/*@+immediatetrans@*/

    } else {

        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

    }

cleanup:

    return status;
}

/************************************************************************/

/*
 * since this AIM does not support field variables or CAPS bounds, the
 * following functions are not required to be filled in except for aimDiscr
 * which just stores away our bodies and aimFreeDiscr that does some cleanup
 * when CAPS terminates
 */


void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *aflr2Instance;

#ifdef DEBUG
    printf(" aflr2AIM/aimCleanup!\n");
#endif
    aflr2Instance = (aimStorage *) instStore;

    // Clean up aflr2Instance data
    status = destroy_aimStorage(aflr2Instance, (int)false);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr2AIM aimStorage cleanup!!!\n", status);

    AIM_FREE(aflr2Instance);
}
