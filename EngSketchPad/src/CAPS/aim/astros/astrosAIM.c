/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             ASTROS AIM
 *
 *     Written by Dr. Ryan Durscher and Dr. Ed Alyanak AFRL/RQVC
 *
 *     This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 *
 */

/*! \mainpage Introduction
 *
 * \section overviewAstros Astros AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with the finite element structural solver ASTROS.
 *
 * Current issues include:
 *  - A thorough bug testing needs to be undertaken.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsAstros and
 * \ref aimOutputsAstros and \ref attributeAstros, respectively.
 *
 * The Astros AIM can automatically execute Astros, with details provided in \ref aimExecuteAstros.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferAstros
 *
 * \section clearanceAstros Clearance Statement
 *  This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 *
 */


/*! \page attributeAstros Astros AIM attributes
 * The following list of attributes are required for the Astros AIM inside the geometry input.
 *
 * - <b> capsDiscipline</b> This attribute is a requirement if doing aeroelastic analysis within Astros. capsDiscipline allows
 * the AIM to determine which bodies are meant for structural analysis and which are used for aerodynamics. Options
 * are:  Structure and Aerodynamic (case insensitive).
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
 * - <b> capsConnect</b> This is a name assigned to any geometric body where the user wishes to create
 * "fictitious" connections such as springs, dampers, and/or rigid body connections to. The user must manually specify
 * the connection between two <c>capsConnect</c> entities using the "Connect" tuple (see \ref aimInputsAstros).
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsConnect $springStart</c>.
 *
 * - <b> capsConnectLink</b> Similar to <c>capsConnect</c>, this is a name assigned to any geometric body where
 * the user wishes to create "fictitious" connections to. A connection is automatically made if a <c>capsConnectLink</c>
 * matches a <c>capsConnect</c> group. Again, further specifics of the connection are input using the "Connect"
 * tuple (see \ref aimInputsAstros). Recall that a string in ESP starts with a $.
 * For example, attribute <c>capsConnectLink $springEnd</c>.
 *
 * - <b> capsBound </b> This is used to mark surfaces on the structural grid in which data transfer with an external
 * solver will take place. See \ref dataTransferAstros for additional details.
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
#include "vlmSpanSpace.h" // Auto spanwise distribution for vortex lattice
#include "feaUtils.h"     // FEA utilities
#include "nastranUtils.h" // Nastran utilities
#include "astrosUtils.h"  // Astros utilities

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#define SLASH     '\\'
#else
#define SLASH     '/'
#endif

#define MXCHAR  255

//#define DEBUG

enum aimInputs
{
  Proj_Name = 1,               /* index is 1-based */
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
  File_Format,
  Mesh_File_Format,
  Design_Variable,
  Design_Variable_Relation,
  Design_Constraint,
  ObjectiveMinMax,
  ObjectiveResponseType,
  VLM_Surface,
  Support,
  Connect,
  Parameter,
  Mesh,
  NUMINPUT = Mesh              /* Total number of inputs */
};

enum aimOutputs
{
  EigenValue = 1,              /* index is 1-based */
  EigenRadian,
  EigenFrequency,
  EigenGeneralMass,
  EigenGeneralStiffness,
  Tmax,
  T1max,
  T2max,
  T3max,
  NUMOUTPUT = T3max            /* Total number of outputs */
};


typedef struct {

    // Project name
    const char *projectName; // Project name

    feaProblemStruct feaProblem;

    feaUnitsStruct units; // units system

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

} aimStorage;


static int initiate_aimStorage(aimStorage *astrosInstance)
{

    int status;

    // Set initial values for astrosInstance
    astrosInstance->projectName = NULL;

    /*
    // Check to make sure data transfer is ok
    astrosInstance->dataTransferCheck = (int) true;
     */

    status = initiate_feaUnitsStruct(&astrosInstance->units);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&astrosInstance->attrMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to constraint index map
    status = initiate_mapAttrToIndexStruct(&astrosInstance->constraintMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to load index map
    status = initiate_mapAttrToIndexStruct(&astrosInstance->loadMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for transfer to index map
    status = initiate_mapAttrToIndexStruct(&astrosInstance->transferMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for connect to index map
    status = initiate_mapAttrToIndexStruct(&astrosInstance->connectMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for response to index map
    status = initiate_mapAttrToIndexStruct(&astrosInstance->responseMap);
    if (status != CAPS_SUCCESS) return status;

    status = initiate_feaProblemStruct(&astrosInstance->feaProblem);
    if (status != CAPS_SUCCESS) return status;

    // Mesh holders
    astrosInstance->numMesh = 0;
    astrosInstance->feaMesh = NULL;

    return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *astrosInstance)
{

    int status;
    int i;

    status = destroy_feaUnitsStruct(&astrosInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&astrosInstance->attrMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to constraint index map
    status = destroy_mapAttrToIndexStruct(&astrosInstance->constraintMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to load index map
    status = destroy_mapAttrToIndexStruct(&astrosInstance->loadMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Transfer to index map
    status = destroy_mapAttrToIndexStruct(&astrosInstance->transferMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Connect to index map
    status = destroy_mapAttrToIndexStruct(&astrosInstance->connectMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Response to index map
    status = destroy_mapAttrToIndexStruct(&astrosInstance->responseMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (astrosInstance->feaMesh != NULL) {

        for (i = 0; i < astrosInstance->numMesh; i++) {
            status = destroy_meshStruct(&astrosInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        EG_free(astrosInstance->feaMesh);
    }

    astrosInstance->feaMesh = NULL;
    astrosInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&astrosInstance->feaProblem);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    // NULL projetName
    astrosInstance->projectName = NULL;

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(aimStorage *aimInfo, aimStorage *astrosInstance)
{
  // Function return flag
  int status;
  int i, remesh = (int)true;

  // Meshing related variables
  double tessParam[3] = {0.025, 0.001, 15};
  int edgePointMin = 2;
  int edgePointMax = 50;
  int quadMesh = (int) false;

  // analysis input values
  capsValue *TessParams = NULL;
  capsValue *EdgePoint_Min = NULL;
  capsValue *EdgePoint_Max = NULL;
  capsValue *QuadMesh = NULL;

  for (i = 0; i < astrosInstance->numMesh; i++) {
      remesh = remesh && (astrosInstance->feaMesh[i].egadsTess->oclass == EMPTY);
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
      tessParam[0] = TessParams->vals.reals[0]; // Gets multiplied by bounding box size
      tessParam[1] = TessParams->vals.reals[1]; // Gets multiplied by bounding box size
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

  status = fea_createMesh(aimInfo,
                          tessParam,
                          edgePointMin,
                          edgePointMax,
                          quadMesh,
                          &astrosInstance->attrMap,
                          &astrosInstance->constraintMap,
                          &astrosInstance->loadMap,
                          &astrosInstance->transferMap,
                          &astrosInstance->connectMap,
                          &astrosInstance->responseMap,
                          &astrosInstance->numMesh,
                          &astrosInstance->feaMesh,
                          &astrosInstance->feaProblem );
  if (status != CAPS_SUCCESS) return status;

  return CAPS_SUCCESS;
}


static int _combineVLM(char *type, int numfeaAero, feaAeroStruct feaAero[],
                       int combineID, feaAeroStruct *combine)
{

    int status;

    int i, j, k;
    int sectionIndex;
    int found = (int) false, skip;

    vlmSectionStruct *tempSection;

    for (i=0; i < numfeaAero; i++) {
        if (strcasecmp(feaAero[i].vlmSurface.surfaceType, type) == 0) break;
    }

    if (i >= numfeaAero) {
        printf("SurfaceType, %s, not found!\n", type);
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    combine->name = EG_strdup(feaAero[i].name);
    if (combine->name == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }
    combine->vlmSurface.surfaceType = EG_strdup(feaAero[i].vlmSurface.surfaceType);
    if (combine->vlmSurface.surfaceType == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    combine->surfaceID = combineID;

    // ADD something for coordinate systems

    // Populate vmlSurface structure
    sectionIndex = 0;
    for (i = 0; i < numfeaAero; i++) {

        if (strcasecmp(feaAero[i].vlmSurface.surfaceType, type) != 0) continue;

        if (found == (int) false) {
            combine->vlmSurface.Cspace = feaAero[i].vlmSurface.Cspace;
            combine->vlmSurface.Sspace = feaAero[i].vlmSurface.Sspace;
            combine->vlmSurface.Nchord = 0;
            combine->vlmSurface.NspanTotal = 0;
            found = (int) true;
        }

        if (combine->vlmSurface.Nchord < feaAero[i].vlmSurface.Nchord) {
            combine->vlmSurface.Nchord = feaAero[i].vlmSurface.Nchord;
        }

        combine->vlmSurface.NspanTotal += feaAero[i].vlmSurface.NspanTotal;

        // Get grids
        combine->numGridID += feaAero[i].numGridID;
        combine->gridIDSet = (int *) EG_reall(combine->gridIDSet,
                                                combine->numGridID*sizeof(int));
        if (combine->gridIDSet == NULL) {
            status = EGADS_MALLOC;
            goto cleanup;
        }

        memcpy(&combine->gridIDSet[combine->numGridID -feaAero[i].numGridID],
               feaAero[i].gridIDSet,
               feaAero[i].numGridID*sizeof(int));

        // Copy section information
        for (j = 0; j < feaAero[i].vlmSurface.numSection; j++) {

            skip = (int) false;
            for (k = 0; k < combine->vlmSurface.numSection; k++) {

                // Check geometry
                status = EG_isEquivalent(combine->vlmSurface.vlmSection[k].ebody,
                                         feaAero[i].vlmSurface.vlmSection[j].ebody);
                if (status == EGADS_SUCCESS) {
                    skip = (int) true;
                    break;
                }

                // Check geometry
                status = EG_isSame(combine->vlmSurface.vlmSection[k].ebody,
                                   feaAero[i].vlmSurface.vlmSection[j].ebody);
                if (status == EGADS_SUCCESS) {
                    skip = (int) true;
                    break;
                }
            }

            if (skip == (int) true) continue;

            combine->vlmSurface.numSection += 1;

            tempSection = (vlmSectionStruct *) EG_reall(combine->vlmSurface.vlmSection,
                                                        combine->vlmSurface.numSection*sizeof(vlmSectionStruct));

            if (tempSection == NULL) {
                combine->vlmSurface.numSection -= 1;
                status = EGADS_MALLOC;
                goto cleanup;
            }

            combine->vlmSurface.vlmSection = tempSection;

            status = initiate_vlmSectionStruct(&combine->vlmSurface.vlmSection[sectionIndex]);
            if (status != CAPS_SUCCESS) {
                combine->vlmSurface.numSection -= 1;
                goto cleanup;
            }

            // Copy the section data - This also copies the control data for the section
            status = copy_vlmSectionStruct( &feaAero[i].vlmSurface.vlmSection[j],
                                            &combine->vlmSurface.vlmSection[sectionIndex]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // Reset the sectionIndex that is keeping track of the section order.
            combine->vlmSurface.vlmSection[sectionIndex].sectionIndex = sectionIndex;

            sectionIndex += 1;
        }
    }

    // Order cross sections for the surface - just in case
    status = vlm_orderSections(combine->vlmSurface.numSection, combine->vlmSurface.vlmSection);
    if (status != CAPS_SUCCESS) goto cleanup;

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS)
      printf("\tPremature exit in _combineVLM, status = %d\n", status);
    return status;
 }


static int createVLMMesh(void *instStore, void *aimInfo, capsValue *aimInputs)
{

    int projectionMethod = (int) true;

    int status, status2; // Function return status

    int i, j, k, surfaceIndex = 0, sectionIndex, transferIndex; // Indexing

    char *analysisType = NULL;
    aimStorage *astrosInstance;

    // Bodies
    const char *intents;
    int   numBody; // Number of Bodies
    ego  *bodies;

    // Aeroelastic information
    int numVLMSurface = 0;
    vlmSurfaceStruct *vlmSurface = NULL;
    int numSpanWise = 0;

    int wingCheck = (int) false, finCheck = (int) false, canardCheck = (int) false;
    int feaAeroTempCombineCount = 0, type=-1;
    feaAeroStruct *feaAeroTemp = NULL, *feaAeroTempCombine = NULL;

    // Vector variables
    double A[3], B[3], C[3], D[3], P[3], p[3], N[3], n[3], d_proj[3];

    double *a, *b, *c, *d;
    double apbArea, apcArea, cpdArea, bpdArea, Area;

    feaMeshDataStruct *feaData;

    astrosInstance = (aimStorage *) instStore;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    if (status != CAPS_SUCCESS) return status;

#ifdef DEBUG
    printf(" astrosAIM/createVLMMesh nbody = %d!\n", numBody);
#endif

    if ((numBody <= 0) || (bodies == NULL)) {
#ifdef DEBUG
        printf(" astrosAIM/createVLMMesh No Bodies!\n");
#endif
        return CAPS_SOURCEERR;
    }

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    // Get aerodynamic reference quantities
    status = fea_retrieveAeroRef(numBody, bodies, &astrosInstance->feaProblem.feaAeroRef);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Cleanup Aero storage first
    if (astrosInstance->feaProblem.feaAero != NULL) {

        for (i = 0; i < astrosInstance->feaProblem.numAero; i++) {
            status = destroy_feaAeroStruct(&astrosInstance->feaProblem.feaAero[i]);
            if (status != CAPS_SUCCESS) goto cleanup;
        }

        EG_free(astrosInstance->feaProblem.feaAero);
    }

    astrosInstance->feaProblem.numAero = 0;

    // Get AVL surface information
    if (aimInputs[VLM_Surface-1].nullVal != IsNull) {

        status = get_vlmSurface(aimInputs[VLM_Surface-1].length,
                                aimInputs[VLM_Surface-1].vals.tuple,
                                &astrosInstance->attrMap,
                                0.0, // default Cspace
                                &numVLMSurface,
                                &vlmSurface);
        if (status != CAPS_SUCCESS) goto cleanup;

    } else {
        printf("An analysis type of Aeroelastic set but no VLM_Surface tuple specified\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    printf("\nGetting FEA vortex lattice mesh\n");

    status = vlm_getSections(aimInfo,
                             numBody,
                             bodies,
                             "Aerodynamic",
                             astrosInstance->attrMap,
                             vlmPLANEYZ,
                             numVLMSurface,
                             &vlmSurface);
    if (status != CAPS_SUCCESS) goto cleanup;
    if (vlmSurface == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    for (i = 0; i < numVLMSurface; i++) {

        // Compute auto spacing
        if (vlmSurface[i].NspanTotal > 0)
            numSpanWise = vlmSurface[i].NspanTotal;
        else if (vlmSurface[i].NspanSection > 0)
            numSpanWise = (vlmSurface[i].numSection-1)*vlmSurface[i].NspanSection;
        else {
            AIM_ERROR(aimInfo  , "Only one of numSpanTotal and numSpanPerSection can be non-zero!");
            AIM_ADDLINE(aimInfo, "    numSpanTotal      = %d", vlmSurface[i].NspanTotal);
            AIM_ADDLINE(aimInfo, "    numSpanPerSection = %d", vlmSurface[i].NspanSection);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        status = vlm_equalSpaceSpanPanels(aimInfo, numSpanWise,
                                          vlmSurface[i].numSection,
                                          vlmSurface[i].vlmSection);
        AIM_STATUS(aimInfo, status);
    }

    // Split the surfaces that have more than 2 sections into a new surface
    for (i = 0; i < numVLMSurface; i++) {

        if (vlmSurface->numSection < 2) {
            printf("Error: Surface '%s' has less than two-sections!\n", vlmSurface[i].name);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(&astrosInstance->transferMap,
                                         vlmSurface[i].name,
                                         &transferIndex);
        if (status == CAPS_NOTFOUND) {
            printf("\tA corresponding capsBound name not found for \"%s\". Surface will be ignored!\n",
                   vlmSurface[i].name);
            continue;
        } else if (status != CAPS_SUCCESS) goto cleanup;

        for (j = 0; j < vlmSurface[i].numSection-1; j++) {

            // Increment the number of Aero surfaces
            astrosInstance->feaProblem.numAero += 1;

            surfaceIndex = astrosInstance->feaProblem.numAero - 1;

            // Allocate
            feaAeroTemp = (feaAeroStruct *) EG_reall(astrosInstance->feaProblem.feaAero,
                      astrosInstance->feaProblem.numAero*sizeof(feaAeroStruct));


            if (feaAeroTemp == NULL) {
                astrosInstance->feaProblem.numAero -= 1;

                status = EGADS_MALLOC;
                goto cleanup;
            }

            astrosInstance->feaProblem.feaAero = feaAeroTemp;

            // Initiate feaAeroStruct
            status = initiate_feaAeroStruct(&astrosInstance->feaProblem.feaAero[surfaceIndex]);
            if (status != CAPS_SUCCESS) goto cleanup;

            // Get surface Name - copy from original surface
            astrosInstance->feaProblem.feaAero[surfaceIndex].name = EG_strdup(vlmSurface[i].name);
            if (astrosInstance->feaProblem.feaAero[surfaceIndex].name == NULL) {
              status = EGADS_MALLOC;
              goto cleanup;
            }

            // Get surface ID - Multiple by 1000 !!
            astrosInstance->feaProblem.feaAero[surfaceIndex].surfaceID =
                                        1000*astrosInstance->feaProblem.numAero;

            // ADD something for coordinate systems

            // Sections aren't necessarily stored in order coming out of vlm_getSections, however sectionIndex is!
            sectionIndex = vlmSurface[i].vlmSection[j].sectionIndex;

            // Populate vmlSurface structure
            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.Cspace = vlmSurface[i].Cspace;
            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.Sspace = vlmSurface[i].Sspace;

            // use the section span count for the sub-surface
            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.NspanTotal = vlmSurface[i].vlmSection[sectionIndex].Nspan;
            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.Nchord     = vlmSurface[i].Nchord;

            // Copy surface type
            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.surfaceType = EG_strdup(vlmSurface[i].surfaceType);

            // Copy section information
            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.numSection = 2;

            astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection = (vlmSectionStruct *) EG_alloc(2*sizeof(vlmSectionStruct));
            if (astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            for (k = 0; k < astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.numSection; k++) {

                // Add k to section indexing variable j to get j and j+1 during iterations

                // Sections aren't necessarily stored in order coming out of vlm_getSections, however sectionIndex is!
                sectionIndex = vlmSurface[i].vlmSection[j+k].sectionIndex;

                status = initiate_vlmSectionStruct(&astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k]);
                if (status != CAPS_SUCCESS) goto cleanup;

                // Copy the section data - This also copies the control data for the section
                status = copy_vlmSectionStruct( &vlmSurface[i].vlmSection[sectionIndex],
                                                &astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k]);
                if (status != CAPS_SUCCESS) goto cleanup;

                // Reset the sectionIndex that is keeping track of the section order.
                astrosInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k].sectionIndex = k;
            }
        }
    }

    // Determine which grid points are to be used for each spline
    for (i = 0; i < astrosInstance->feaProblem.numAero; i++) {

        // Debug
        //printf("\tDetermining grid points\n");

        // Get the transfer index for this surface - it has already been checked to make sure the name is in the
        // transfer index map
        AIM_NOTNULL(astrosInstance->feaProblem.feaAero, aimInfo, status);
        status = get_mapAttrToIndexIndex(&astrosInstance->transferMap,
                                         astrosInstance->feaProblem.feaAero[i].name,
                                         &transferIndex);
        if (status != CAPS_SUCCESS) goto cleanup;

        if (projectionMethod == (int) false) { // Look for attributes

            for (j = 0; j < astrosInstance->feaProblem.feaMesh.numNode; j++) {

                if (astrosInstance->feaProblem.feaMesh.node[j].analysisType == MeshStructure) {
                    feaData = (feaMeshDataStruct *) astrosInstance->feaProblem.feaMesh.node[j].analysisData;
                } else {
                    continue;
                }

                if (feaData->transferIndex != transferIndex) continue;
                if (feaData->transferIndex == CAPSMAGIC) continue;


                astrosInstance->feaProblem.feaAero[i].numGridID += 1;
                k = astrosInstance->feaProblem.feaAero[i].numGridID;

                astrosInstance->feaProblem.feaAero[i].gridIDSet = (int *)
                       EG_reall(astrosInstance->feaProblem.feaAero[i].gridIDSet,
                                k*sizeof(int));

                if (astrosInstance->feaProblem.feaAero[i].gridIDSet == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                astrosInstance->feaProblem.feaAero[i].gridIDSet[k-1] =
                    astrosInstance->feaProblem.feaMesh.node[j].nodeID;
            }

        } else { // Projection method


            /*
             *   n = A X B Create a normal vector/ plane between A and B
             *
             *   d_proj = C - (C * n)*n/ ||n||^2 , projection of point d on plane created by AxB
             *
             *   p = D - (D * n)*n/ ||n||^2 , projection of point p on plane created by AxB
             *
             *                         (section 2)
             *                     LE(c)---------------->TE(d)
             *   Grid Point       -^                   ^ -|
             *           |^      -            -         - |
             *           | -     A      -   C          - d_proj
             *           |  D   -    -                 -
             *           |   - - -     (section 1     -
             *           p    LE(a)----------B------->TE(b)
             */

            // Vector between section 2 and 1
            a = astrosInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[0].xyzLE;
            b = astrosInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[0].xyzTE;
            c = astrosInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[1].xyzLE;
            d = astrosInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[1].xyzTE;

            // Debug
            //printf("a = %f %f %f\n", a[0], a[1], a[2]);
            //printf("b = %f %f %f\n", b[0], b[1], b[2]);
            //printf("c = %f %f %f\n", c[0], c[1], c[2]);
            //printf("d = %f %f %f\n", d[0], d[1], d[2]);

            // Vector between LE of section 1 and LE of section 2
            A[0] = c[0] - a[0];
            A[1] = c[1] - a[1];
            A[2] = c[2] - a[2];

            // Vector between LE and TE of section 1
            B[0] = b[0] - a[0];
            B[1] = b[1] - a[1];
            B[2] = b[2] - a[2];

            // Vector between LE of section 1 and TE of section 2
            C[0] = d[0] - a[0];
            C[1] = d[1] - a[1];
            C[2] = d[2] - a[2];

            // Normal vector between A and B
            cross_DoubleVal(A, B, N);

            // Normalize normal vector
            n[0] = N[0]/sqrt(dot_DoubleVal(N, N));
            n[1] = N[1]/sqrt(dot_DoubleVal(N, N));
            n[2] = N[2]/sqrt(dot_DoubleVal(N, N));

            //printf("n = %f %f %f\n", n[0], n[1], n[2]);

            // Projection of vector C on plane created by AxB
            d_proj[0] = C[0] - dot_DoubleVal(C, n)*n[0] + a[0];
            d_proj[1] = C[1] - dot_DoubleVal(C, n)*n[1] + a[1];
            d_proj[2] = C[2] - dot_DoubleVal(C, n)*n[2] + a[2];

            //printf("d_proj = %f %f %f\n", d_proj[0], d_proj[1], d_proj[2]);

            // Vector between LE of section 1 and TE of section 2 where the TE as been projected on A x B plane
            C[0] = d_proj[0] - a[0];
            C[1] = d_proj[1] - a[1];
            C[2] = d_proj[2] - a[2];

            // Area of the rectangle (first triangle)
            cross_DoubleVal(A, C, N);
            Area = 0.5*sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);

            // Area of the rectangle (second triangle)
            cross_DoubleVal(C, B, N);
            Area = 0.5*sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2])  + Area;

            // Debug
            //printf("\tArea = %f\n",Area);

            for (j = 0; j < astrosInstance->feaProblem.feaMesh.numNode; j++) {

                if (astrosInstance->feaProblem.feaMesh.node[j].analysisType == MeshStructure) {
                    feaData = (feaMeshDataStruct *) astrosInstance->feaProblem.feaMesh.node[j].analysisData;
                } else {
                    continue;
                }

                if (feaData->transferIndex != transferIndex) continue;
                if (feaData->transferIndex == CAPSMAGIC) continue;

                D[0] = astrosInstance->feaProblem.feaMesh.node[j].xyz[0] - a[0];

                D[1] = astrosInstance->feaProblem.feaMesh.node[j].xyz[1] - a[1];

                D[2] = astrosInstance->feaProblem.feaMesh.node[j].xyz[2] - a[2];

                // Projection of vector D on plane created by AxB
                p[0] = D[0] - dot_DoubleVal(D, n)*n[0] + a[0];
                p[1] = D[1] - dot_DoubleVal(D, n)*n[1] + a[1];
                p[2] = D[2] - dot_DoubleVal(D, n)*n[2] + a[2];

                // First triangle
                A[0] = a[0] - p[0];
                A[1] = a[1] - p[1];
                A[2] = a[2] - p[2];

                B[0] = b[0] - p[0];
                B[1] = b[1] - p[1];
                B[2] = b[2] - p[2];
                cross_DoubleVal(A, B, P);
                apbArea = 0.5*sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2]);

                // Second triangle
                A[0] = a[0] - p[0];
                A[1] = a[1] - p[1];
                A[2] = a[2] - p[2];

                B[0] = c[0] - p[0];
                B[1] = c[1] - p[1];
                B[2] = c[2] - p[2];
                cross_DoubleVal(A, B, P);
                apcArea = 0.5*sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2]);

                // Third triangle
                A[0] = c[0] - p[0];
                A[1] = c[1] - p[1];
                A[2] = c[2] - p[2];

                B[0] = d_proj[0] - p[0];
                B[1] = d_proj[1] - p[1];
                B[2] = d_proj[2] - p[2];
                cross_DoubleVal(A, B, P);
                cpdArea = 0.5*sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2]);

                // Fourth triangle
                A[0] = b[0] - p[0];
                A[1] = b[1] - p[1];
                A[2] = b[2] - p[2];

                B[0] = d_proj[0] - p[0];
                B[1] = d_proj[1] - p[1];
                B[2] = d_proj[2] - p[2];
                cross_DoubleVal(A, B, P);
                bpdArea = 0.5*sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2]);

                // Debug
                //printf("Area of Triangle  = %f\n", (apbArea + apcArea + cpdArea + bpdArea));
                //if (fabs(apbArea + apcArea + cpdArea + bpdArea - Area) > 1E-5) {
                //  printf("Point index %d\n", j);
                //  printf("\tAreas don't match - %f vs %f!!\n", Area, (apbArea + apcArea + cpdArea + bpdArea));
                //  printf("\tPoint - %f %f %f\n", D[0] + a[0], D[1] + a[1], D[2] + a[2]);
                //  printf("\tPoint projection - %f %f %f\n", p[0], p[1], p[2]);
                //}

                if (fabs(apbArea + apcArea + cpdArea + bpdArea - Area) > 1E-5)  continue;

                astrosInstance->feaProblem.feaAero[i].numGridID += 1;

                if (astrosInstance->feaProblem.feaAero[i].numGridID == 1) {
                    astrosInstance->feaProblem.feaAero[i].gridIDSet = (int *)
                  EG_alloc(astrosInstance->feaProblem.feaAero[i].numGridID*sizeof(int));
                } else {
                    astrosInstance->feaProblem.feaAero[i].gridIDSet = (int *)
                  EG_reall(astrosInstance->feaProblem.feaAero[i].gridIDSet,
                           astrosInstance->feaProblem.feaAero[i].numGridID*sizeof(int));
                }

                if (astrosInstance->feaProblem.feaAero[i].gridIDSet == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                astrosInstance->feaProblem.feaAero[i].gridIDSet[
                    astrosInstance->feaProblem.feaAero[i].numGridID-1] =
              astrosInstance->feaProblem.feaMesh.node[j].nodeID;
            }
        }

        if (astrosInstance->feaProblem.feaAero[i].numGridID > 0) {
            printf("\tSurface %d: Number of points found for aero-spline = %d\n",
                   i+1, astrosInstance->feaProblem.feaAero[i].numGridID );
        }
        else {
            printf("\tError: No points found for aero-spline for surface %d\n", i+1);
            status = CAPS_NOTFOUND;
            goto cleanup;
        }
    }

    // Need to combine all aero surfaces into one for static, opt and trim analysis
    if (strcasecmp(analysisType, "Aeroelastic") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrimOpt") == 0 ) {

        printf("\t(Re-)Combining all aerodynamic surfaces into a 'Wing', 'Canard', and/or  'Fin' single surfaces !\n");

        if (astrosInstance->feaProblem.feaAero == NULL) {
            status = CAPS_NULLVALUE;
            goto cleanup;
        }

        for (i = 0; i < astrosInstance->feaProblem.numAero; i++) {
            if (astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType == NULL){
                printf("DEVELOPER ERROR: no surfaceType set (surfcae index %d)!\n", i);
                status = CAPS_BADVALUE;
                goto cleanup;
            }
//          printf("Surface Type = |%s|\n",
//                 astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType);

            if        (wingCheck == (int) false &&
                       strcasecmp(astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType,
                                  "Wing") == 0) {
                type = 0;
                wingCheck = (int) true;
            } else if (canardCheck == (int) false &&
                       strcasecmp(astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType,
                                  "Canard") == 0) {
                type = 1;
                canardCheck = (int) true;
            } else if (finCheck == (int) false &&
                       strcasecmp(astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType,
                                  "Fin") == 0) {
                type = 2;
                finCheck = (int) true;
            } else {
//                printf("Skipping, surface %d, %d %d %d,  %d %d %d\n", i, wingCheck, canardCheck, finCheck,
//                        strcasecmp(astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType, "Wing"),
//                        strcasecmp(astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType, "Canard"),
//                        strcasecmp(astrosInstance->feaProblem.feaAero[i].vlmSurface.surfaceType, "Fin"));
                continue;
            }

            feaAeroTemp = (feaAeroStruct *) EG_reall(feaAeroTempCombine,
                                                     (feaAeroTempCombineCount+1)*
                                                     sizeof(feaAeroStruct));
            if (feaAeroTemp == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            feaAeroTempCombine = feaAeroTemp;

            status = initiate_feaAeroStruct(&feaAeroTempCombine[feaAeroTempCombineCount]);
            if (status != CAPS_SUCCESS) goto cleanup;

//            printf("TYPE %d\n", type);

            switch (type) {
                case 0:
                    status = _combineVLM("Wing",
                                         astrosInstance->feaProblem.numAero,
                                         astrosInstance->feaProblem.feaAero,
                                         1000*(feaAeroTempCombineCount+1),
                                         &feaAeroTempCombine[feaAeroTempCombineCount]);
                    if (status != CAPS_SUCCESS) goto cleanup;
                    break;
                case 1:
                    status = _combineVLM("Canard",
                                         astrosInstance->feaProblem.numAero,
                                         astrosInstance->feaProblem.feaAero,
                                         1000*(feaAeroTempCombineCount+1),
                                         &feaAeroTempCombine[feaAeroTempCombineCount]);
                    if (status != CAPS_SUCCESS) goto cleanup;
                    break;

                case 2:
                    status = _combineVLM("Fin",
                                         astrosInstance->feaProblem.numAero,
                                         astrosInstance->feaProblem.feaAero,
                                         1000*(feaAeroTempCombineCount+1),
                                         &feaAeroTempCombine[feaAeroTempCombineCount]);
                    if (status != CAPS_SUCCESS) goto cleanup;
                    break;
            }

            feaAeroTempCombineCount += 1;
        }

        // Free old feaProblem Aero
        if (astrosInstance->feaProblem.feaAero != NULL) {

            for (i = 0; i < astrosInstance->feaProblem.numAero; i++) {
                status = destroy_feaAeroStruct(&astrosInstance->feaProblem.feaAero[i]);
                if (status != CAPS_SUCCESS)
                    printf("Status %d during destroy_feaAeroStruct\n", status);
            }
        }

        if (astrosInstance->feaProblem.feaAero != NULL)
            EG_free(astrosInstance->feaProblem.feaAero);
        astrosInstance->feaProblem.feaAero = NULL;
        astrosInstance->feaProblem.numAero = 0;

        // Point to new data
        astrosInstance->feaProblem.numAero = feaAeroTempCombineCount;
        astrosInstance->feaProblem.feaAero = feaAeroTempCombine;
    }

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS)
        printf("\tPremature exit in createVLMMesh, status = %d\n", status);

    if (status != CAPS_SUCCESS && feaAeroTempCombine != NULL) {

        for (i = 0; i < feaAeroTempCombineCount+1; i++) {
            (void) destroy_feaAeroStruct(&feaAeroTempCombine[i]);
        }
/*@-kepttrans@*/
        EG_free(feaAeroTempCombine);
/*@+kepttrans@*/
    }

    if (vlmSurface != NULL) {

        for (i = 0; i < numVLMSurface; i++) {
            status2 = destroy_vlmSurfaceStruct(&vlmSurface[i]);
            if (status2 != CAPS_SUCCESS)
                printf("\tPremature exit in destroy_vlmSurfaceStruct, status = %d\n",
                       status2);
        }
    }

    if (vlmSurface != NULL) EG_free(vlmSurface);
    numVLMSurface = 0;

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

    aimStorage *astrosInstance=NULL;

#ifdef DEBUG
    printf("astrosAIM/aimInitialize   instance = %d!\n", inst);
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

    // Allocate astrosInstance
    AIM_ALLOC(astrosInstance, 1, aimStorage, aimInfo, status);
    *instStore = astrosInstance;

    // Initialize instance storage
    (void) initiate_aimStorage(astrosInstance);

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

    /*! \page aimInputsAstros AIM Inputs
     * The following list outlines the Astros inputs along with their default value available
     * through the AIM interface. Unless noted these values will be not be linked to
     * any parent AIMs with variables of the same name.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" astrosAIM/aimInputs  index = %d!\n", index);
#endif

    *ainame = NULL;

    // Astros Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("astros_CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsAstros
         * - <B> Proj_Name = "astros_CAPS"</B> <br>
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

        /*! \page aimInputsAstros
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
         * zero ignores this phase.
         */

    } else if (index == Edge_Point_Min) {
        *ainame               = EG_strdup("Edge_Point_Min");
        defval->type          = Integer;
        defval->vals.integer  = 2;
        defval->lfixed        = Fixed;
        defval->nrow          = 1;
        defval->ncol          = 1;
        defval->nullVal       = NotNull;

        /*! \page aimInputsAstros
         * - <B> aimInputsAstros = 2</B> <br>
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

        /*! \page aimInputsAstros
         * - <B> Edge_Point_Max = 50</B> <br>
         * Maximum number of points on an edge including end points to use when creating a surface mesh (min 2).
         */


    } else if (index == Quad_Mesh) {
        *ainame               = EG_strdup("Quad_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = (int) false;

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
         * - <B> Analysis = NULL</B> <br>
         * Analysis tuple used to input analysis/case information for the model, see \ref feaAnalysis for additional details.
         */
    } else if (index == Analysis_Type) {
        *ainame              = EG_strdup("Analysis_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Modal");
        defval->lfixed       = Change;

        /*! \page aimInputsAstros
         * - <B> Analysis_Type = "Modal"</B> <br>
         * Type of analysis to generate files for, options include "Modal", "Static", "AeroelasticTrim", "AeroelasticTrimOpt", "AeroelasticFlutter", and "Optimization".
         * Note: "Aeroelastic" and "StaticOpt" are still supported and refer to "AeroelasticTrim" and "Optimization".
         */
    } else if (index == File_Format) {
        *ainame              = EG_strdup("File_Format");
        defval->type         = String;
        defval->vals.string  = EG_strdup("Small"); // Small, Large, Free
        defval->lfixed       = Change;

        /*! \page aimInputsAstros
         * - <B> File_Format = "Small"</B> <br>
         * Formatting type for the bulk file. Options: "Small", "Large", "Free".
         */

    } else if (index == Mesh_File_Format) {
        *ainame              = EG_strdup("Mesh_File_Format");
        defval->type         = String;
        defval->vals.string  = EG_strdup("Free"); // Small, Large, Free
        defval->lfixed       = Change;

        /*! \page aimInputsAstros
         * - <B> Mesh_File_Format = "Free"</B> <br>
         * Formatting type for the mesh file. Options: "Small", "Large", "Free".
         */

    } else if (index == Design_Variable) {
        *ainame              = EG_strdup("Design_Variable");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAstros
         * - <B> Design_Variable = NULL</B> <br>
         * The design variable tuple used to input design variable information for the model optimization, see \ref feaDesignVariable for additional details.
         */

    } else if (index == Design_Variable_Relation) {
        *ainame              = EG_strdup("Design_Variable_Relation");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAstros
         * - <B> Design_Variable_Relation = NULL</B> <br>
         * The design variable relation tuple is used to input design variable relation information for the model optimization, see \ref feaDesignVariableRelation for additional details.
         */

    }  else if (index == Design_Constraint) {
        *ainame              = EG_strdup("Design_Constraint");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAstros
         * - <B> Design_Constraint = NULL</B> <br>
         * The design constraint tuple used to input design constraint information for the model optimization, see \ref feaDesignConstraint for additional details.
         */

    } else if (index == ObjectiveMinMax) {
        *ainame              = EG_strdup("ObjectiveMinMax");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Max"); // Max, Min
        defval->lfixed       = Change;

        /*! \page aimInputsAstros
         * - <B> ObjectiveMinMax = "Max"</B> <br>
         * Maximize or minimize the design objective during an optimization. Option: "Max" or "Min".
         */

    } else if (index == ObjectiveResponseType) {
        *ainame              = EG_strdup("ObjectiveResponseType");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Weight"); // Weight
        defval->lfixed       = Change;

        /*! \page aimInputsAstros
         * - <B> ObjectiveResponseType = "Weight"</B> <br>
         * Object response type (see Astros manual).
         */
    } else if (index == VLM_Surface) {
        *ainame              = EG_strdup("VLM_Surface");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAstros
         * - <B>VLM_Surface = NULL </B> <br>
         * Vortex lattice method tuple input. See \ref vlmSurface for additional details.
         */
    } else if (index == Support) {
        *ainame              = EG_strdup("Support");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
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

        /*! \page aimInputsAstros
         * - <B> Parameter = NULL</B> <br>
         * Parameter tuple used to define user entries. This can be used to input things to ASTROS such as CONVERT or MFORM etc.
         * The input is in Tuple form ("DATACARD", "DATAVALUE"). All inputs are strings.  Example: ("CONVERT","MASS,  0.00254").
         * Note:  Inputs assume a "," delimited entry.  Notice the "," after MASS in the Example.
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

        /*! \page aimInputsAstros
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
    int found;

    int status; // Status return

    //int found; // Boolean operator

    int *tempIntegerArray = NULL; // Temporary array to store a list of integers
    char *noQuoteString = NULL;

    // Analysis information
    char *analysisType = NULL;
    int optFlag; // 0 - ANALYSIS, 1 - OPTIMIZATION Set based on analysisType char input

    // Optimization Information
    //char *objectiveResp = NULL;
    const char *geomInName;

    // File IO
    char *filename = NULL; // Output file name
    FILE *fp = NULL; // Output file pointer
    int addComma = (int) false; // Add comma between inputs

    feaLoadStruct *feaLoad;
    int numThermalLoad=0, numGravityLoad=0;
    int nGeomIn;
    capsValue *geomInVal;
    aimStorage *astrosInstance;

#ifdef DEBUG
    // Bodies
    const char *intents;
    int   numBody; // Number of Bodies
    ego  *bodies;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    if (status != CAPS_SUCCESS) return status;

    printf(" astrosAIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif
    astrosInstance = (aimStorage *) instStore;
    if (aimInputs == NULL) return CAPS_NULLVALUE;

    // Get project name
    astrosInstance->projectName = aimInputs[Proj_Name-1].vals.string;

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    // Get FEA mesh if we don't already have one
    if (aim_newGeometry(aimInfo) == CAPS_SUCCESS) {

        status = checkAndCreateMesh(aimInfo, astrosInstance);
        if (status != CAPS_SUCCESS) goto cleanup;

        // Get Aeroelastic mesh
        if( strcasecmp(analysisType, "Aeroelastic") == 0 ||
            strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
            strcasecmp(analysisType, "AeroelasticTrimOpt") == 0 ||
            strcasecmp(analysisType, "AeroelasticFlutter") == 0) {

            status = createVLMMesh(instStore, aimInfo, aimInputs);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    // Note: Setting order is important here.
    // 1. Materials should be set before properties.
    // 2. Coordinate system should be set before mesh and loads
    // 3. Mesh should be set before loads, constraints, and supports
    // 4. Constraints and loads should be set before analysis
    // 5. Optimization should be set after properties, but before analysis

    // Set material properties
    if (aimInputs[Material-1].nullVal == NotNull) {
        status = fea_getMaterial(aimInfo,
                                 aimInputs[Material-1].length,
                                 aimInputs[Material-1].vals.tuple,
                                 &astrosInstance->units,
                                 &astrosInstance->feaProblem.numMaterial,
                                 &astrosInstance->feaProblem.feaMaterial);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Material tuple is NULL - No materials set\n");

    // Set property properties
    if (aimInputs[Property-1].nullVal == NotNull) {
        status = fea_getProperty(aimInfo,
                                 aimInputs[Property-1].length,
                                 aimInputs[Property-1].vals.tuple,
                                 &astrosInstance->attrMap,
                                 &astrosInstance->units,
                                 &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;

        // Assign element "subtypes" based on properties set
        status = fea_assignElementSubType(astrosInstance->feaProblem.numProperty,
                                          astrosInstance->feaProblem.feaProperty,
                                          &astrosInstance->feaProblem.feaMesh);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Property tuple is NULL - No properties set\n");

    // Set constraint properties
    if (aimInputs[Constraint-1].nullVal == NotNull) {
        status = fea_getConstraint(aimInputs[Constraint-1].length,
                                   aimInputs[Constraint-1].vals.tuple,
                                   &astrosInstance->constraintMap,
                                   &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Constraint tuple is NULL - No constraints applied\n");

    // Set support properties
    if (aimInputs[Support-1].nullVal == NotNull) {
        status = fea_getSupport(aimInputs[Support-1].length,
                                aimInputs[Support-1].vals.tuple,
                                &astrosInstance->constraintMap,
                                &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Support tuple is NULL - No supports applied\n");

    // Set connection properties
    if (aimInputs[Connect-1].nullVal == NotNull) {
        status = fea_getConnection(aimInputs[Connect-1].length,
                                   aimInputs[Connect-1].vals.tuple,
                                   &astrosInstance->connectMap,
                                   &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;


        // Unify all connectionID's for RBE2 cards sake to be used for MPC in case control
        for (i = 0; i < astrosInstance->feaProblem.numConnect; i++) {

            if (astrosInstance->feaProblem.feaConnect[i].connectionType == RigidBody
                || astrosInstance->feaProblem.feaConnect[i].connectionType == RigidBodyInterpolate) {

                astrosInstance->feaProblem.feaConnect[i].connectionID = astrosInstance->feaProblem.numConnect+1;
            }
        }

    } else printf("Connect tuple is NULL - Using defaults\n");

    // Set load properties
    if (aimInputs[Load-1].nullVal == NotNull) {
        status = fea_getLoad(aimInputs[Load-1].length,
                             aimInputs[Load-1].vals.tuple,
                             &astrosInstance->loadMap,
                             &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;

        // Loop through loads to see if any of them are supposed to be from an external source
        for (i = 0; i < astrosInstance->feaProblem.numLoad; i++) {

            if (astrosInstance->feaProblem.feaLoad[i].loadType == PressureExternal) {

                // Transfer external pressures from the AIM discrObj
                status = fea_transferExternalPressure(aimInfo,
                                                      &astrosInstance->feaProblem.feaMesh,
                                                      &astrosInstance->feaProblem.feaLoad[i]);
                if (status != CAPS_SUCCESS) return status;

            } // End PressureExternal if
        } // End load for loop
    } else printf("Load tuple is NULL - No loads applied\n");

    // Set design variables
    if (aimInputs[Design_Variable-1].nullVal == NotNull) {
        status = fea_getDesignVariable(aimInfo,
                                       (int)true,
                                       aimInputs[Design_Variable-1].length,
                                       aimInputs[Design_Variable-1].vals.tuple,
                                       aimInputs[Design_Variable_Relation-1].length,
                                       aimInputs[Design_Variable_Relation-1].vals.tuple,
                                       &astrosInstance->attrMap,
                                       &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Design_Variable tuple is NULL - No design variables applied\n");

    // Set design constraints
    if (aimInputs[Design_Constraint-1].nullVal == NotNull) {
        status = fea_getDesignConstraint(aimInputs[Design_Constraint-1].length,
                                         aimInputs[Design_Constraint-1].vals.tuple,
                                         &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;
    } else printf("Design_Constraint tuple is NULL - No design constraints applied\n");

    // Set analysis settings
    if (aimInputs[Analysix-1].nullVal == NotNull) {
        status = fea_getAnalysis(aimInputs[Analysix-1].length,
                                 aimInputs[Analysix-1].vals.tuple,
                                 &astrosInstance->feaProblem);
        if (status != CAPS_SUCCESS) return status;
    } else {
        printf("Analysis tuple is NULL\n"); // Its ok to not have an analysis tuple we will just create one

        status = fea_createDefaultAnalysis(&astrosInstance->feaProblem, analysisType);
        if (status != CAPS_SUCCESS) return status;
    }


    // Set file format type
    if        (strcasecmp(aimInputs[File_Format-1].vals.string, "Small") == 0) {
        astrosInstance->feaProblem.feaFileFormat.fileType = SmallField;
    } else if (strcasecmp(aimInputs[File_Format-1].vals.string, "Large") == 0) {
        astrosInstance->feaProblem.feaFileFormat.fileType = LargeField;
    } else if (strcasecmp(aimInputs[File_Format-1].vals.string, "Free") == 0)  {
        astrosInstance->feaProblem.feaFileFormat.fileType = FreeField;
    } else {
        printf("Unrecognized \"File_Format\", valid choices are [Small, Large, or Free]. Reverting to default\n");
    }

    // Set grid file format type
    if        (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Small") == 0) {
        astrosInstance->feaProblem.feaFileFormat.gridFileType = SmallField;
    } else if (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Large") == 0) {
        astrosInstance->feaProblem.feaFileFormat.gridFileType = LargeField;
    } else if (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Free") == 0)  {
        astrosInstance->feaProblem.feaFileFormat.gridFileType = FreeField;
    } else {
        printf("Unrecognized \"Mesh_File_Format\", valid choices are [Small, Large, or Free]. Reverting to default\n");
    }

    // Write Astros Mesh
    filename = EG_alloc(MXCHAR +1);
    if (filename == NULL) return EGADS_MALLOC;
    strcpy(filename, astrosInstance->projectName);

    status = mesh_writeAstros(aimInfo,
                              filename,
                              1,
                              &astrosInstance->feaProblem.feaMesh,
                              astrosInstance->feaProblem.feaFileFormat.gridFileType,
                              astrosInstance->feaProblem.numDesignVariable,
                              astrosInstance->feaProblem.feaDesignVariable,
                              1.0);
    if (status != CAPS_SUCCESS) {
        EG_free(filename);
        return status;
    }

    // Write Astros subElement types not supported by astros_writeMesh
    strcat(filename, ".bdf");
    fp = aim_fopen(aimInfo, filename, "a");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        EG_free(filename);
        return CAPS_IOERR;
    }
    EG_free(filename);

    printf("Writing subElement types (if any) - appending mesh file\n");
    status = astros_writeSubElementCard(fp,
                                         &astrosInstance->feaProblem.feaMesh,
                                         astrosInstance->feaProblem.numProperty,
                                         astrosInstance->feaProblem.feaProperty,
                                         &astrosInstance->feaProblem.feaFileFormat);
    if (status != CAPS_SUCCESS) goto cleanup;

    // Connections
    for (i = 0; i < astrosInstance->feaProblem.numConnect; i++) {

        if (i == 0) {
            printf("Writing connection cards - appending mesh file\n");
        }

        status = astros_writeConnectionCard(fp,
                                            &astrosInstance->feaProblem.feaConnect[i],
                                            &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }
    if (fp != NULL) fclose(fp);
    fp = NULL;


    // Write astros input file
    filename = EG_alloc(MXCHAR +1);
    if (filename == NULL) return EGADS_MALLOC;
    strcpy(filename, astrosInstance->projectName);
    strcat(filename, ".dat");


    printf("\nWriting Astros instruction file....\n");
    fp = aim_fopen(aimInfo, filename, "w");
    EG_free(filename);
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        return CAPS_IOERR;
    }

    // define file format delimiter type
    /*
    if (astrosInstance->feaProblem.feaFileFormat.fileType == FreeField) {
        delimiter = ",";
    } else {
        delimiter = " ";
    }
    */

    //////////////// Executive control ////////////////
    fprintf(fp, "ASSIGN DATABASE CAPS PASS NEW\n");

    //////////////// Case control ////////////////
    fprintf(fp, "SOLUTION\n");
    fprintf(fp, "TITLE = %s\n", astrosInstance->projectName);
    // Analysis type
    if     (strcasecmp(analysisType, "Modal")              == 0 ||
            strcasecmp(analysisType, "Static")             == 0 ||
            strcasecmp(analysisType, "AeroelasticTrim")    == 0 ||
            strcasecmp(analysisType, "AeroelasticFlutter") == 0 ||
            strcasecmp(analysisType, "Aeroelastic")        == 0) {

        fprintf(fp, "ANALYZE\n");
        optFlag = 0;
    }
    else if(strcasecmp(analysisType, "StaticOpt")          == 0 ||
            strcasecmp(analysisType, "Optimization")       == 0 ||
            strcasecmp(analysisType, "AeroelasticTrimOpt") == 0) {

        fprintf(fp, "OPTIMIZE STRATEGY=((FSD,10), (MP,20)), MAXITER=30, NRFAC=1.5,\n");
        fprintf(fp, "EPS= 1.00, MOVLIM=1.5\n");
        fprintf(fp, "PRINT DCON=ALL, GDES=ALL, GPWG=ALL\n");
        optFlag = 1;
    }

    else {
        AIM_ERROR(aimInfo, "Unrecognized \"Analysis_Type\", %s\n", analysisType);
        status = CAPS_BADVALUE;
        goto cleanup;
    }


    // Set up the case information
    if (astrosInstance->feaProblem.numAnalysis == 0) {
        printf("Error: No analyses in the feaProblem! (this shouldn't be possible)\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }


    // Write sub-case information if multiple analysis tuples were provide - will always have at least 1
    for (i = 0; i < astrosInstance->feaProblem.numAnalysis; i++) {

        // Write boundary constraints/supports/etc.
        fprintf(fp, " BOUNDARY");

        addComma = (int) false;
        // Write support for sub-case
        if (astrosInstance->feaProblem.feaAnalysis[i].numSupport != 0) {

            if (astrosInstance->feaProblem.feaAnalysis[i].numSupport > 1) {
                printf("\tWARNING: More than 1 support is not supported at this time for a given case!\n");

            } else {
                fprintf(fp, " SUPPORT = %d ", astrosInstance->feaProblem.feaAnalysis[i].supportSetID[0]);
                addComma = (int) true;
            }
        }

        // Write constraint for sub-case - see warning statement below for behavior
        if (astrosInstance->feaProblem.numConstraint != 0) {
            if (addComma == (int) true) fprintf(fp,",");

            fprintf(fp, " SPC = %d ", astrosInstance->feaProblem.numConstraint+i+1); //TODO - change to i+1 to just i
            addComma = (int) true;
        }

        // Issue some warnings regarding constraints if necessary
        if (astrosInstance->feaProblem.feaAnalysis[i].numConstraint == 0 && astrosInstance->feaProblem.numConstraint !=0) {

            printf("\tWarning: No constraints specified for case %s, assuming "
                    "all constraints are applied!!!!\n", astrosInstance->feaProblem.feaAnalysis[i].name);

        } else if (astrosInstance->feaProblem.numConstraint == 0) {

            printf("\tWarning: No constraints specified for case %s!!!!\n", astrosInstance->feaProblem.feaAnalysis[i].name);
        }

        // Write MPC for sub-case - currently only supported when we have RBE2 elements - see above for unification - TODO - investigate
        for (j = 0; j < astrosInstance->feaProblem.numConnect; j++) {

            if (astrosInstance->feaProblem.feaConnect[j].connectionType == RigidBody
                || astrosInstance->feaProblem.feaConnect[i].connectionType == RigidBodyInterpolate) {

                if (addComma == (int) true) fprintf(fp,",");

                fprintf(fp, " MPC = %d ", astrosInstance->feaProblem.feaConnect[j].connectionID);
                addComma = (int) true;
                break;
            }
        }

        if (astrosInstance->feaProblem.feaAnalysis[i].analysisType == Modal ||
            astrosInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter ) {

            if (addComma == (int) true) fprintf(fp,",");
            fprintf(fp," METHOD = %d ", astrosInstance->feaProblem.feaAnalysis[i].analysisID);
            //addComma = (int) true;
        }

        fprintf(fp, "\n"); // End boundary line

        fprintf(fp, "    LABEL = %s\n", astrosInstance->feaProblem.feaAnalysis[i].name);

        // Write discipline
        if (astrosInstance->feaProblem.feaAnalysis[i].analysisType == Static) { // Static

            fprintf(fp, "    STATICS ");

            // Issue some warnings regarding loads if necessary
            if (astrosInstance->feaProblem.feaAnalysis[i].numLoad == 0 && astrosInstance->feaProblem.numLoad !=0) {
                printf("\tWarning: No loads specified for static case %s, assuming "
                        "all loads are applied!!!!\n", astrosInstance->feaProblem.feaAnalysis[i].name);

            } else if (astrosInstance->feaProblem.numLoad == 0) {
                printf("\tWarning: No loads specified for static case %s!!!!\n", astrosInstance->feaProblem.feaAnalysis[i].name);
            }

            addComma = (int) false;
            found = (int) false;
            numThermalLoad = 0;
            numGravityLoad = 0;

            if (astrosInstance->feaProblem.numLoad !=0) {
                fprintf(fp, "(");

                for (k = 0; k < astrosInstance->feaProblem.numLoad; k++) {

                    feaLoad = &astrosInstance->feaProblem.feaLoad[k];

                    if (astrosInstance->feaProblem.feaAnalysis[i].numLoad != 0) { // if loads specified in analysis

                        for (j = 0; j < astrosInstance->feaProblem.feaAnalysis[i].numLoad; j++) { // See if the load is in the loadSet

                            if (feaLoad->loadID == astrosInstance->feaProblem.feaAnalysis[i].loadSetID[j] ) break;
                        }

                        if (j >= astrosInstance->feaProblem.feaAnalysis[i].numLoad) continue; // If it isn't in the loadSet move on
                    } else {
                        //pass
                    }

                    if (feaLoad->loadType == Thermal && numThermalLoad == 0) {

                        if (addComma == (int) true) fprintf(fp, ",");
                        fprintf(fp, " THERMAL = %d", feaLoad->loadID);
                        addComma = (int) true;

                        numThermalLoad += 1;
                        if (numThermalLoad > 1) {
                            printf("More than 1 Thermal load found - astrosAIM does NOT currently doesn't support multiple thermal loads in a given case!\n");
                        }

                        continue;
                    }

                    if (feaLoad->loadType == Gravity && numGravityLoad == 0) {

                        if (addComma == (int) true) fprintf(fp, ",");
                        fprintf(fp, " GRAVITY = %d", feaLoad->loadID);
                        addComma = (int) true;

                        numGravityLoad += 1;
                        if (numGravityLoad > 1) {
                            printf("More than 1 Gravity load found - astrosAIM does NOT currently doesn't support multiple gravity loads in a given case!\n");
                        }

                        continue;
                    }

                    found = (int) true;
                }

                if (found == (int) true) {
                    if (addComma == (int) true) fprintf(fp, ",");
                    fprintf(fp, " MECH = %d", astrosInstance->feaProblem.numLoad+i+1);
                }

                fprintf(fp,")");
            }

            if (optFlag == 0) {
                fprintf(fp, "\n");
                fprintf(fp, "    PRINT DISP=ALL, STRESS=ALL\n");
            } else {
                fprintf(fp, ", CONST( STRESS = %d)\n", astrosInstance->feaProblem.numDesignConstraint+i+1);
                fprintf(fp, "    PRINT DISP(ITER=LAST)=ALL, STRESS(ITER=LAST)=ALL\n");
            }
        }

        if (astrosInstance->feaProblem.feaAnalysis[i].analysisType == Modal) {// Modal
            fprintf(fp, "    MODES\n");
            fprintf(fp, "    PRINT (MODES=ALL) DISP=ALL, ROOT=ALL\n");
        }
        if (astrosInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim) {// Trim
            fprintf(fp, "    SAERO SYMMETRIC (TRIM=%d)", astrosInstance->feaProblem.feaAnalysis[i].analysisID);

            if (optFlag == 0) {
                fprintf(fp, "\n");
                fprintf(fp, "    PRINT DISP=ALL, GPWG=ALL, TRIM, TPRE=ALL, STRESS=ALL\n");
            } else {
                fprintf(fp, ", CONST(STRESS = %d)\n", astrosInstance->feaProblem.numDesignConstraint+i+1);
                fprintf(fp, "    PRINT (ITER=LAST) DISP=ALL, GPWG=ALL, TRIM, TPRE=ALL, STRESS=ALL\n");
            }

        }

        if (astrosInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {// Flutter
            fprintf(fp, "    MODES\n");
            fprintf(fp, "    FLUTTER (FLCOND = %d)\n", astrosInstance->feaProblem.feaAnalysis[i].analysisID);
            fprintf(fp, "    PRINT (MODES=ALL) DISP=ALL, ROOT=ALL\n");
        }
    }

//    } else { // If no sub-cases
//
//        if (strcasecmp(analysisType, "Modal") == 0) {
//            printf("Warning: No eigenvalue analysis information specified in \"Analysis\" tuple, through "
//                    "AIM input \"Analysis_Type\" is set to \"Modal\"!!!!\n");
//            return CAPS_NOTFOUND;
//        }
//
//        fprintf(fp, "  BOUNDARY");
//
//        // Write support information
//        addComma = (int) false;
//        if (astrosInstance->feaProblem.numSupport != 0) {
//            if (astrosInstance->feaProblem.numSupport > 1) {
//                printf("\tWARNING: More than 1 support is not supported at this time for a given case!\n");
//            } else {
//                fprintf(fp, " SUPPORT = %d ", astrosInstance->feaProblem.numSupport+1);
//                addComma = (int) true;
//            }
//        }
//
//        // Write constraint information
//        if (astrosInstance->feaProblem.numConstraint != 0) {
//            if (addComma == (int) true) fprintf(fp, ",");
//            fprintf(fp, " SPC = %d ", astrosInstance->feaProblem.numConstraint+1);
//        } else {
//            printf("\tWarning: No constraints specified for job!!!!\n");
//        }
//
//        // Write MPC for sub-case - currently only supported when we have RBE2 elements - see above for unification
//        for (j = 0; j < astrosInstance->feaProblem.numConnect; j++) {
//
//            if (astrosInstance->feaProblem.feaConnect[j].connectionType == RigidBody) {
//
//                if (addComma == (int) true) fprintf(fp,",");
//
//                fprintf(fp, " MPC = %d ", astrosInstance->feaProblem.feaConnect[j].connectionID);
//                addComma = (int) true;
//                break;
//            }
//        }
//
//        fprintf(fp, "\n");
//
//        // Write discipline
//        if (strcasecmp(analysisType, "Static") == 0) { // Static loads
//
//            fprintf(fp, "    STATICS");
//
//            // Write loads for sub-case
//            if (astrosInstance->feaProblem.numLoad != 0) {
//                fprintf(fp, "(");
//
//                addComma = (int) false;
//
//                k = 0;
//                for (j = 0; j < astrosInstance->feaProblem.numLoad; j++) {
//                    if (astrosInstance->feaProblem.feaLoad[j].loadType == Thermal) {
//
//                        if (addComma == (int) true) fprintf(fp, ",");
//                        fprintf(fp, " THERMAL = %d", astrosInstance->feaProblem.feaLoad[j].loadID);
//                        addComma = (int) true;
//
//                        numThermalLoad += 1;
//                        if (numThermalLoad > 1) {
//                            printf("More than 1 Thermal load found - astrosAIM does NOT currently doesn't support multiple thermal loads!\n");
//                        }
//
//                        continue;
//                    }
//
//                    if (astrosInstance->feaProblem.feaLoad[j].loadType == Gravity) {
//                        if (addComma == (int) true) fprintf(fp, ",");
//                        fprintf(fp, " GRAVITY = %d", astrosInstance->feaProblem.feaLoad[j].loadID);
//                        addComma = (int) true;
//                        continue;
//                    }
//
//                    numGravityLoad += 1;
//                    if (numGravityLoad > 1) {
//                        printf("More than 1 Gravity load found - astrosAIM does NOT currently doesn't support multiple gravity loads!\n");
//                    }
//
//                    k += 1;
//                }
//
//                if (k != 0) {
//                    if (addComma == (int) true) fprintf(fp, ",");
//                    fprintf(fp, " MECH = %d", astrosInstance->feaProblem.numLoad+1);
//                }
//
//                fprintf(fp,")");
//            } else {
//                printf("\tWarning: No loads specified for job!!!!\n");
//            }
//
//            fprintf(fp, "\n");
//
//            fprintf(fp, "    PRINT DISP=ALL, STRESS=ALL\n");
//        }
//    }


    fprintf(fp, "END\n$\n"); // End Case control

    //////////////// Bulk data ////////////////
    fprintf(fp, "BEGIN BULK(SORT)\n");
    fprintf(fp, "$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|\n");

    //PRINT Parameter ENTRIES IN BULK DATA

    if (aimInputs[Parameter-1].nullVal == NotNull) {
        for (i = 0; i < aimInputs[Parameter-1].length; i++) {
            noQuoteString = string_removeQuotation(aimInputs[Parameter-1].vals.tuple[i].value);
            AIM_NOTNULL(noQuoteString, aimInfo, status);
            fprintf(fp, "%s, %s\n", aimInputs[Parameter-1].vals.tuple[i].name, noQuoteString);
            EG_free(noQuoteString);
        }
    }

    // Turn off auto SPC
    //fprintf(fp, "%-8s %7s %7s\n", "PARAM", "AUTOSPC", "N");

    // Optimization Objective Response Response, SOL 200 only
    if (strcasecmp(analysisType, "StaticOpt") == 0 || strcasecmp(analysisType, "Optimization") == 0) {
        /*
        objectiveResp = aimInputs[ObjectiveResponseType-1].vals.string;
        if     (strcasecmp(objectiveResp, "Weight") == 0) objectiveResp = "WEIGHT";
        else {
            AIM_ERROR(aimInfo, "\tUnrecognized \"ObjectiveResponseType\", %s\n", objectiveResp);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        fprintf(fp,"%-8s", "DRESP1");

        tempString = convert_integerToString(1, 7, 1);
        fprintf(fp, "%s%s", delimiter, tempString);
        EG_free(tempString);

        fprintf(fp, "%s%7s", delimiter, objectiveResp);
        fprintf(fp, "%s%7s", delimiter, objectiveResp);


        fprintf(fp, "\n");
         */
    }

    // Write AERO Card
    if (strcasecmp(analysisType, "AeroelasticFlutter") == 0 ) {

        printf("\tWriting aero card\n");
        status = astros_writeAEROCard(fp,
                                      &astrosInstance->feaProblem.feaAeroRef,
                                      &astrosInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Write AESTAT and AESURF cards
    if (strcasecmp(analysisType, "Aeroelastic") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrimOpt") == 0) {

        printf("\tWriting aeros card\n");
        status = astros_writeAEROSCard(fp,
                                       &astrosInstance->feaProblem.feaAeroRef,
                                       &astrosInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);

        // No AESTAT Cards in ASTROS


        // fprintf(fp,"\n");
    }

    // Analysis Cards - Eigenvalue and design objective included, as well as combined load, constraint, and design constraints
    for (i = 0; i < astrosInstance->feaProblem.numAnalysis; i++) {

        if (i == 0) printf("\tWriting analysis cards\n");

        status = astros_writeAnalysisCard(fp,
                                          &astrosInstance->feaProblem.feaAnalysis[i],
                                          &astrosInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);

        if (astrosInstance->feaProblem.feaAnalysis[i].numLoad != 0) {

            // Create a temporary list of load IDs
            tempIntegerArray = (int *) EG_alloc(astrosInstance->feaProblem.feaAnalysis[i].numLoad*sizeof(int));
            if (tempIntegerArray == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            k = 0;
            for (j = 0; j <  astrosInstance->feaProblem.feaAnalysis[i].numLoad; j++) {
                for (l = 0; l < astrosInstance->feaProblem.numLoad; l++) {
                    if (astrosInstance->feaProblem.feaAnalysis[i].loadSetID[j] == astrosInstance->feaProblem.feaLoad[l].loadID) break;
                }

                if (l >= astrosInstance->feaProblem.numLoad) continue;
                if (astrosInstance->feaProblem.feaLoad[l].loadType == Gravity) continue;
                if (astrosInstance->feaProblem.feaLoad[l].loadType == Thermal) continue;
                tempIntegerArray[k] = astrosInstance->feaProblem.feaLoad[l].loadID;
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
                                              astrosInstance->feaProblem.numLoad+i+1,
                                              k,
                                              tempIntegerArray,
                                              astrosInstance->feaProblem.feaLoad,
                                              &astrosInstance->feaProblem.feaFileFormat);

            if (status != CAPS_SUCCESS) goto cleanup;

            // Free temporary load ID list
            EG_free(tempIntegerArray);
            tempIntegerArray = NULL;

        } else { // If no loads for an individual analysis are specified assume that all loads should be applied

            if (astrosInstance->feaProblem.numLoad != 0) {

                // Create a temporary list of load IDs
                tempIntegerArray = (int *) EG_alloc(astrosInstance->feaProblem.numLoad*sizeof(int));
                if (tempIntegerArray == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                k = 0;
                for (j = 0; j < astrosInstance->feaProblem.numLoad; j++) {
                    if (astrosInstance->feaProblem.feaLoad[j].loadType == Gravity) continue;
                    if (astrosInstance->feaProblem.feaLoad[j].loadType == Thermal) continue;
                    tempIntegerArray[k] = astrosInstance->feaProblem.feaLoad[j].loadID;
                    k += 1;
                }

                tempIntegerArray = (int *) EG_reall(tempIntegerArray, k*sizeof(int));
                if (tempIntegerArray == NULL)  {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }

                //TOOO: eliminate load add card?
                // Write combined load card
                printf("\tWriting load ADD cards\n");
                status = nastran_writeLoadADDCard(fp,
                                                  astrosInstance->feaProblem.numLoad+i+1,
                                                  k,
                                                  tempIntegerArray,
                                                  astrosInstance->feaProblem.feaLoad,
                                                  &astrosInstance->feaProblem.feaFileFormat);

                if (status != CAPS_SUCCESS) goto cleanup;

                // Free temporary load ID list
                EG_free(tempIntegerArray);
                tempIntegerArray = NULL;
            }

        }

        if (astrosInstance->feaProblem.feaAnalysis[i].numConstraint != 0) {

            // Write combined constraint card
            printf("\tWriting constraint cards--each subcase individually\n");
            fprintf(fp,"$\n$ Constraint(s)\n");

            for (j = 0; j < astrosInstance->feaProblem.feaAnalysis[i].numConstraint; j++) {
                k = astrosInstance->feaProblem.feaAnalysis[i].constraintSetID[j] - 1;

                // one spc set per subcase, each different
                status = astros_writeConstraintCard(fp,
                                                    astrosInstance->feaProblem.numConstraint+i+1,
                                                    &astrosInstance->feaProblem.feaConstraint[k],
                                                    &astrosInstance->feaProblem.feaFileFormat);
                if (status != CAPS_SUCCESS) goto cleanup;
            }

//                printf("\tWriting constraint ADD cards\n");
//                status = nastran_writeConstraintADDCard(fp,
//                                                        astrosInstance->feaProblem.numConstraint+i+1,
//                                                        astrosInstance->feaProblem.feaAnalysis[i].numConstraint,
//                                                        astrosInstance->feaProblem.feaAnalysis[i].constraintSetID,
//                                                        &astrosInstance->feaProblem.feaFileFormat);
//                AIM_STATUS(aimInfo, status);

        } else { // If no constraints for an individual analysis are specified assume that all constraints should be applied

            if (astrosInstance->feaProblem.numConstraint != 0) {

                // Write combined constraint card
                printf("\tWriting constraint cards--all constraints for each subcase\n");
                fprintf(fp,"$\n$ Constraint(s)\n");

                for (j = 0; j < astrosInstance->feaProblem.numConstraint; j++) {

                    // one spc set per subcase, each the same
                    status = astros_writeConstraintCard(fp,
                                                        astrosInstance->feaProblem.numConstraint+i+1,
                                                        &astrosInstance->feaProblem.feaConstraint[j],
                                                        &astrosInstance->feaProblem.feaFileFormat);
                    if (status != CAPS_SUCCESS) goto cleanup;
                }
            }
        }

        if (astrosInstance->feaProblem.feaAnalysis[i].numDesignConstraint != 0) {

            printf("\tWriting design constraint cards--no subcases\n");
            fprintf(fp,"$\n$ Design constraint(s)\n");
            for( j = 0; j < astrosInstance->feaProblem.feaAnalysis[i].numDesignConstraint; j++) {
                k = astrosInstance->feaProblem.feaAnalysis[i].designConstraintSetID[j] - 1;

                // one design constraint set per subcase analysis, each may be different
                status = astros_writeDesignConstraintCard(fp,
                                                          astrosInstance->feaProblem.numDesignConstraint+i+1,
                                                          &astrosInstance->feaProblem.feaDesignConstraint[k],
                                                          astrosInstance->feaProblem.numMaterial,
                                                          astrosInstance->feaProblem.feaMaterial,
                                                          astrosInstance->feaProblem.numProperty,
                                                          astrosInstance->feaProblem.feaProperty,
                                                          &astrosInstance->feaProblem.feaFileFormat);
                if (status != CAPS_SUCCESS) goto cleanup;
            }

        } else { // If no design constraints for an individual analysis are specified assume that all design constraints should be applied

            if (astrosInstance->feaProblem.numDesignConstraint != 0) {

                printf("\tWriting design constraint cards\n");
                fprintf(fp,"$\n$ Design constraint(s)\n");
                for( j = 0; j < astrosInstance->feaProblem.numDesignConstraint; j++) {

                    // one design constraint set per subcase analysis, all the same
                    status = astros_writeDesignConstraintCard(fp,
                                                              astrosInstance->feaProblem.numDesignConstraint+i+1,
                                                              &astrosInstance->feaProblem.feaDesignConstraint[j],
                                                              astrosInstance->feaProblem.numMaterial,
                                                              astrosInstance->feaProblem.feaMaterial,
                                                              astrosInstance->feaProblem.numProperty,
                                                              astrosInstance->feaProblem.feaProperty,
                                                              &astrosInstance->feaProblem.feaFileFormat);
                    if (status != CAPS_SUCCESS) goto cleanup;
                }
            }

        }
    }

    // Loads
    for (i = 0; i < astrosInstance->feaProblem.numLoad; i++) {

        if (i == 0) {
            printf("\tWriting load cards\n");
            fprintf(fp,"$\n$ Load(s)\n");
        }

        status = astros_writeLoadCard(fp,
                                      &astrosInstance->feaProblem.feaMesh,
                                      &astrosInstance->feaProblem.feaLoad[i],
                                      &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Constraints
    // Move to subcase level because ASTROS does not support SPCADD card--DB 8 Mar 18
    /* for (i = 0; i < astrosInstance->feaProblem.numConstraint; i++) {

        if (i == 0) {
            printf("\tWriting constraint cards\n");
            fprintf(fp,"$\n$ Constraint(s)\n");
        }

        status = nastran_writeConstraintCard(fp,
                                             &astrosInstance->feaProblem.feaConstraint[i],
                                             &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    } */

    // Supports
    for (i = 0; i < astrosInstance->feaProblem.numSupport; i++) {

        if (i == 0) {
            printf("\tWriting support cards\n");
            fprintf(fp,"$\n$ Support(s)\n");
        }

        status = astros_writeSupportCard(fp,
                                         &astrosInstance->feaProblem.feaSupport[i],
                                         &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }


    // Materials
    for (i = 0; i < astrosInstance->feaProblem.numMaterial; i++) {

        if (i == 0) {
            printf("\tWriting material cards\n");
            fprintf(fp,"$\n$ Material(s)\n");
        }

        status = nastran_writeMaterialCard(fp,
                                           &astrosInstance->feaProblem.feaMaterial[i],
                                           &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Properties
    for (i = 0; i < astrosInstance->feaProblem.numProperty; i++) {

        if (i == 0) {
            printf("\tWriting property cards\n");
            fprintf(fp,"$\n$ Property(ies)\n");
        }

        status = astros_writePropertyCard(fp,
                                          &astrosInstance->feaProblem.feaProperty[i],
                                          &astrosInstance->feaProblem.feaFileFormat,
                                          astrosInstance->feaProblem.numDesignVariable,
                                          astrosInstance->feaProblem.feaDesignVariable);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Coordinate systems
    for (i = 0; i < astrosInstance->feaProblem.numCoordSystem; i++) {

        if (i == 0) {
            printf("\tWriting coordinate system cards\n");
            fprintf(fp,"$\n$ Coordinate system(s)\n");
        }

        status = nastran_writeCoordinateSystemCard(fp, &astrosInstance->feaProblem.feaCoordSystem[i], &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Optimization - design variables
    for (i = 0; i < astrosInstance->feaProblem.numDesignVariable; i++) {

        if (i == 0) {
            printf("\tWriting design variables and analysis - design variable relation cards\n");
            fprintf(fp,"$\n$ Design variable(s)\n");
        }

        status = astros_writeDesignVariableCard(fp,
                                                &astrosInstance->feaProblem.feaDesignVariable[i],
                                                astrosInstance->feaProblem.numProperty,
                                                astrosInstance->feaProblem.feaProperty,
                                                &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
    }

    // Optimization - design variables - geometry
    nGeomIn = aim_getIndex(aimInfo, NULL, GEOMETRYIN);
    if (nGeomIn > 0) {
      status = aim_getValue(aimInfo, 1, GEOMETRYIN, &geomInVal);
      if (status != CAPS_SUCCESS) {
        printf("Error: Cannot get Geometry In Value Structures\n");
        nGeomIn = 0;
      }
    }
    for (i = 0; i < astrosInstance->feaProblem.numDesignVariable; i++) {

        // Geometric parameterization - only if needed
        for (j = 0; j < nGeomIn; j++) {

            status = aim_getName(aimInfo, j+1, GEOMETRYIN, &geomInName);
            if (status != CAPS_SUCCESS) goto cleanup;

            if (strcmp(astrosInstance->feaProblem.feaDesignVariable[i].name, geomInName) == 0) break;
        }

        // If name isn't found in Geometry inputs skip write geometric design variables
        if (j >= nGeomIn) continue;

        if(aim_getGeomInType(aimInfo, j+1) != 0) {
            printf("Error: Geometric sensitivity not available for CFGPMTR = %s\n", geomInName);
            status = CAPS_NOSENSITVTY;
            goto cleanup;
        }

        printf(">>> Writing geometry parametrization\n");
        status = astros_writeGeomParametrization(fp,
                                                 aimInfo,
                                                 astrosInstance->feaProblem.numDesignVariable,
                                                 astrosInstance->feaProblem.feaDesignVariable,
                                                 nGeomIn, geomInVal,
                                                 &astrosInstance->feaProblem.feaMesh,
                                                 &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
        printf(">>> Done writing geometry parametrization\n");

        break; // Only need to call astros_writeGeomParametrization once!
    }

    // Optimization - design constraints
    // Move to subcase level because ASTROS does not support DCONADD card--DB 7 Mar 18
    /* for( i = 0; i < astrosInstance->feaProblem.numDesignConstraint; i++) {

          if (i == 0) {
              printf("\tWriting design constraints and responses cards\n");
              fprintf(fp,"$\n$ Design constraint(s)\n");
          }

          status = astros_writeDesignConstraintCard(fp,
                                                    &astrosInstance->feaProblem.feaDesignConstraint[i],
                                                    astrosInstance->feaProblem.numMaterial,
                                                    astrosInstance->feaProblem.feaMaterial,
                                                    astrosInstance->feaProblem.numProperty,
                                                    astrosInstance->feaProblem.feaProperty,
                                                    &astrosInstance->feaProblem.feaFileFormat);
        if (status != CAPS_SUCCESS) goto cleanup;
      } */

    // Aeroelastic
    if (strcasecmp(analysisType, "Aeroelastic") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrimOpt") == 0 ) {

        printf("\tWriting aeroelastic cards\n");
        for (i = 0; i < astrosInstance->feaProblem.numAero; i++){

            status = astros_writeCAeroCard(fp,
                                           &astrosInstance->feaProblem.feaAero[i],
                                           &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = astros_checkAirfoil(aimInfo,
                                         &astrosInstance->feaProblem.feaAero[i]);
            if (status == CAPS_SOURCEERR) {

                j = (int) false;
                printf("\tBody topology used in aerodynamic surface %d, isn't suitable for airfoil shape, switching to panel", i+1);

            } else if (status != CAPS_SUCCESS) {

                goto cleanup;

            } else {

                j = (int) true;
            }

            status = astros_writeAirfoilCard(fp,
                                             j, // useAirfoilShape
                                             &astrosInstance->feaProblem.feaAero[i],
                                             &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = astros_writeAeroData(aimInfo,
                                          fp,
                                          j, // useAirfoilShape
                                          &astrosInstance->feaProblem.feaAero[i],
                                          &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = astros_writeAeroSplineCard(fp,
                                                &astrosInstance->feaProblem.feaAero[i],
                                                &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = nastran_writeSet1Card(fp,
                                           &astrosInstance->feaProblem.feaAero[i],
                                           &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    // Aeroelastic
    if (strcasecmp(analysisType, "AeroelasticFlutter") == 0) {

        printf("\tWriting unsteady aeroelastic cards\n");
        for (i = 0; i < astrosInstance->feaProblem.numAero; i++){

            status = nastran_writeCAeroCard(fp,
                                            &astrosInstance->feaProblem.feaAero[i],
                                            &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = astros_writeAeroSplineCard(fp,
                                                &astrosInstance->feaProblem.feaAero[i],
                                                &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;

            status = nastran_writeSet1Card(fp,
                                           &astrosInstance->feaProblem.feaAero[i],
                                           &astrosInstance->feaProblem.feaFileFormat);
            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    // Include mesh file
    fprintf(fp,"$\nINCLUDE %s.bdf\n$\n", astrosInstance->projectName);

    // End bulk data
    fprintf(fp,"ENDDATA\n");

    fclose(fp);
    fp = NULL;
/*
////////////////////////////////////////
    printf("\n\n\nTESTING OUT READER\n\n");

    // FO6 data variables
    int numGridPoint = 0;
    int numEigenVector = 0;
    double **dataMatrix = NULL;
    filename = (char *) EG_alloc((strlen(astrosInstance->projectName) +
                                  strlen(".out") + 2)*sizeof(char));

    sprintf(filename, "%s%s", astrosInstance->projectName, ".out");

    // Open file
    fp = aim_fopen(aimInfo, filename, "r");
    if (filename != NULL) EG_free(filename);
    filename = NULL;

    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);

        return CAPS_IOERR;
    }

    status = astros_readOUTDisplacement(fp,
                                         -9,
                                         &numGridPoint,
                                         &dataMatrix);

    for ( i = 0 ; i < numGridPoint ; i++) printf("dataMatrix  = %f %f %f %f %f %f %f %f\n",
                                                 dataMatrix[i][0],
                                                 dataMatrix[i][1],
                                                 dataMatrix[i][2],
                                                 dataMatrix[i][3],
                                                 dataMatrix[i][4],
                                                 dataMatrix[i][5],
                                                 dataMatrix[i][6],
                                                 dataMatrix[i][7]);
    fclose(fp);
///////////////////////////////////
*/

    status = CAPS_SUCCESS;
    goto cleanup;

    cleanup:
        if (status != CAPS_SUCCESS) printf("\tPremature exit in astrosAIM preAnalysis, status = %d\n", status);

        if (fp != NULL) fclose(fp);

        EG_free(tempIntegerArray);
        tempIntegerArray = NULL;

        return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  /*! \page aimExecuteAstros AIM Execution
   *
   * If auto execution is enabled when creating an Astros AIM,
   * the AIM will execute Astros just-in-time with the command line:
   *
   * \code{.sh}
   * $ASTROS_ROOT/astros < $Proj_Name.dat > $Proj_Name.out
   * \endcode
   *
   * where preAnalysis generated the file Proj_Name + ".dat" which contains the input information.
   * The environemtn variable ASTROS_ROOT is assumed to point to the location where the
   * "astros.exe" executable and run files "ASTRO.D01" and "ASTRO.IDX" are located.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an Astros AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or Astros can be executed by calling preAnalysis, system call, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * astros.preAnalysis()
   *
   * print ("\n\nRunning......")
   * astros.system(ASTROS_ROOT + os.sep + "astros.exe < " + astros.input.Proj_Name + ".dat > " + astros.input.Proj_Name + ".out"); # Run via system call
   *
   * print ("\n\postAnalysis......")
   * astros.postAnalysis()
   * \endcode
   */

  int status = CAPS_SUCCESS;
  char command[PATH_MAX], *env;
  aimStorage *astrosInstance;
  *state = 0;

  astrosInstance = (aimStorage *) instStore;
  if (astrosInstance == NULL) return CAPS_NULLVALUE;

  env = getenv("ASTROS_ROOT");
  if (env == NULL) {
    AIM_ERROR(aimInfo, "ASTROS_ROOT environment variable is not set!");
    return CAPS_EXECERR;
  }

  snprintf(command, PATH_MAX, "%s%cASTRO.D01", env, SLASH);
  status = aim_cpFile(aimInfo, command, "");
  if (status != CAPS_SUCCESS) return status;

  snprintf(command, PATH_MAX, "%s%cASTRO.IDX", env, SLASH);
  status = aim_cpFile(aimInfo, command, "");
  if (status != CAPS_SUCCESS) return status;

  snprintf(command, PATH_MAX, "%s%castros.exe < %s.dat > %s.out",
           env, SLASH, astrosInstance->projectName, astrosInstance->projectName);

  return aim_system(aimInfo, NULL, command);
}


// ********************** AIM Function Break *****************************
// Check that astros ran without errors
int
aimPostAnalysis(void *instStore, /*@unused@*/ void *aimInfo,
                /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS;

  char *filename = NULL; // File to open
  char extOUT[] = ".out";
  FILE *fp = NULL; // File pointer
  aimStorage *astrosInstance;

  size_t linecap = 0;
  char *line = NULL; // Temporary line holder
  int withErrors = (int) false;
  int terminated = (int) false;

#ifdef DEBUG
  printf(" astrosAIM/aimPostAnalysis!\n");
#endif
  astrosInstance = (aimStorage *) instStore;

  filename = (char *) EG_alloc((strlen(astrosInstance->projectName) +
                                strlen(extOUT) +1)*sizeof(char));
  if (filename == NULL) return EGADS_MALLOC;

  sprintf(filename, "%s%s", astrosInstance->projectName, extOUT);

  fp = aim_fopen(aimInfo, filename, "r");

  EG_free(filename); // Free filename allocation

  if (fp == NULL) {
      AIM_ERROR(aimInfo, " astrosAIM/aimPostAnalysis Cannot open Output file!");

      return CAPS_IOERR;
  }

  // Scan the file for the string
/*@-nullpass@*/
  while( !feof(fp) ) {

      // Get line from file
      status = getline(&line, &linecap, fp);
      if (status < 0) break;

      if (terminated == (int) false) terminated = (int) (strstr(line, "A S T R O S  T E R M I N A T E D") != NULL);
      if (withErrors == (int) false) withErrors = (int) (strstr(line, "W I T H  E R R O R S") != NULL);
  }
/*@+nullpass@*/
  fclose(fp);
  EG_free(line);
  status = CAPS_SUCCESS;

  if (terminated == (int) false) {
    AIM_ERROR(aimInfo, "Astros did not run to termination!");
    status = CAPS_EXECERR;
  }

  if (withErrors == (int) true) {
    AIM_ERROR(aimInfo, "");
    AIM_ADDLINE(aimInfo, "****************************************");
    AIM_ADDLINE(aimInfo, "***                                  ***");
    AIM_ADDLINE(aimInfo, "*** A S T R O S  T E R M I N A T E D ***");
    AIM_ADDLINE(aimInfo, "***      W I T H  E R R O R S        ***");
    AIM_ADDLINE(aimInfo, "***                                  ***");
    AIM_ADDLINE(aimInfo, "****************************************");
    status = CAPS_EXECERR;
  }

  return status;
}


// Set Astros output variables
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsAstros AIM Outputs
     * The following list outlines the Astros outputs available through the AIM interface.
     */

    #ifdef DEBUG
        printf(" astrosAIM/aimOutputs instance = %d  index = %d!\n", iIndex, index);
    #endif

    /*! \page aimOutputsAstros AIM Outputs
     * - <B>EigenValue</B> = List of Eigen-Values (\f$ \lambda\f$) after a modal solve.
     * - <B>EigenRadian</B> = List of Eigen-Values in terms of radians (\f$ \omega = \sqrt{\lambda}\f$ ) after a modal solve.
     * - <B>EigenFrequency</B> = List of Eigen-Values in terms of frequencies (\f$ f = \frac{\omega}{2\pi}\f$) after a modal solve.
     * - <B>EigenGeneralMass</B> = List of generalized masses for the Eigen-Values.
     * - <B>EigenGeneralStiffness</B> = List of generalized stiffness for the Eigen-Values.
     * .
     */

    //printf("***** index=%d *****", index);

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
        *aoname = EG_strdup("Tmax");

    } else if (index == 7) {
        *aoname = EG_strdup("T1max");

    } else if (index == 8) {
        *aoname = EG_strdup("T2max");

    } else if (index == 9) {
        *aoname = EG_strdup("T3max");
    }

    //printf(" %s\n", *aoname);

    if (index <= 5) {
        form->type          = Double;
        form->units         = NULL;
        form->lfixed        = Change;
        form->sfixed        = Change;
        form->vals.reals = NULL;
        form->vals.real = 0;
    } else {
        form->type = Double;
        form->dim  = Vector;
        form->nrow  = 1;
        form->ncol  = 1;
        form->units  = NULL;
        form->vals.reals = NULL;
        form->vals.real = 0;
    }

    return CAPS_SUCCESS;
}


// Calculate Astros output
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int status = CAPS_SUCCESS; // Function return status
    aimStorage *astrosInstance;

    int i; //Indexing

    int numEigenVector;
    double **dataMatrix = NULL;

    char *filename = NULL; // File to open
    char extOUT[] = ".out";
    FILE *fp = NULL; // File pointer

    astrosInstance = (aimStorage *) instStore;

    filename = (char *) EG_alloc((strlen(astrosInstance->projectName) + strlen(extOUT) +1)*sizeof(char));
    if (filename == NULL) return EGADS_MALLOC;

    sprintf(filename, "%s%s", astrosInstance->projectName, extOUT);

    fp = aim_fopen(aimInfo, filename, "r");

    EG_free(filename); // Free filename allocation

    if (fp == NULL) {
#ifdef DEBUG
        printf(" astrosAIM/aimCalcOutput Cannot open Output file!\n");
#endif
        return CAPS_IOERR;
    }

    if (index <= 5) {

        status = astros_readOUTEigenValue(fp, &numEigenVector, &dataMatrix);
        if ((status == CAPS_SUCCESS) && (dataMatrix != NULL)) {

            val->nrow = numEigenVector;
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

        if (dataMatrix != NULL) {
            for (i = 0; i < numEigenVector; i++) {
                if (dataMatrix[i] != NULL) EG_free(dataMatrix[i]);
            }
            EG_free(dataMatrix);
        }

    } else if (index <= 9) {
        double T1max=0, T2max=0, T3max=0, Tmax=0, TT;
        int ipnt, numGridPoint;

        status = astros_readOUTDisplacement(fp, -1, &numGridPoint, &dataMatrix);
        if ((status == CAPS_SUCCESS) && (dataMatrix != NULL)) {
            val->dim  = Scalar;
            val->nrow = 1;
            val->ncol = 1;
            val->length = val->nrow * val->ncol;

            for (ipnt = 0; ipnt < numGridPoint; ipnt++) {
                TT = sqrt(pow(dataMatrix[ipnt][2], 2)
                        + pow(dataMatrix[ipnt][3], 2)
                        + pow(dataMatrix[ipnt][4], 2));

                if (fabs(dataMatrix[ipnt][2]) > T1max) T1max = fabs(dataMatrix[ipnt][2]);
                if (fabs(dataMatrix[ipnt][3]) > T2max) T2max = fabs(dataMatrix[ipnt][3]);
                if (fabs(dataMatrix[ipnt][4]) > T3max) T3max = fabs(dataMatrix[ipnt][4]);
                if (TT                        > Tmax ) Tmax  = TT;
            }

            if        (index == 6) {
                val->vals.real = Tmax;
            } else if (index == 7) {
                val->vals.real = T1max;
            } else if (index == 8) {
                val->vals.real = T2max;
            } else {
                val->vals.real = T3max;
            }
        }

        if (dataMatrix != NULL) {
            for (i = 0; i < numGridPoint; i++) {
                if (dataMatrix[i] != NULL) EG_free(dataMatrix[i]);
            }
            EG_free(dataMatrix);
        }
    }

    if (fp != NULL) fclose(fp);

    return status;
}


void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *astrosInstance;

#ifdef DEBUG
    printf(" astrosAIM/Cleanup  numInstance = %d!\n", numInstance);
#endif

    astrosInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(astrosInstance);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during clean up\n", status);

    EG_free(astrosInstance);
}


int aimDiscr(char *tname, capsDiscr *discr) {

    int status; // Function return status

    int numBody, i;
    aimStorage *astrosInstance;

    // EGADS objects
    ego *bodies = NULL, *tess = NULL;

    const char   *intents;

#ifdef DEBUG
    printf(" astrosAIM/aimDiscr: tname = %s, instance = %d!\n", tname);
#endif

    if (tname == NULL) return CAPS_NOTFOUND;

    astrosInstance = (aimStorage *) discr->instStore;

    /*if (astrosInstance->dataTransferCheck == (int) false) {
        printf("The volume is not suitable for data transfer - possibly the volume mesher "
                "added unaccounted for points\n");
        return CAPS_BADVALUE;
    }*/

    // Currently this ONLY works if the capsTranfer lives on single body!
    status = aim_getBodies(discr->aInfo, &intents, &numBody, &bodies);
    if (status != CAPS_SUCCESS) {
        printf(" astrosAIM/aimDiscr: aim_getBodies = %d!\n", status);
        return status;
    }
    if (bodies == NULL) {
        printf(" astrosAIM/aimDiscr: Null Bodies!\n");
        return CAPS_NULLOBJ;
    }

    // Check and generate/retrieve the mesh
    status = checkAndCreateMesh(discr->aInfo, astrosInstance);
    if (status != CAPS_SUCCESS) goto cleanup;

    AIM_ALLOC(tess, astrosInstance->numMesh, ego, discr->aInfo, status);
    for (i = 0; i < astrosInstance->numMesh; i++) {
      tess[i] = astrosInstance->feaMesh[i].egadsTess;
    }

    status = mesh_fillDiscr(tname, &astrosInstance->attrMap, astrosInstance->numMesh, tess, discr);
    if (status != CAPS_SUCCESS) goto cleanup;

#ifdef DEBUG
    printf(" astrosAIM/aimDiscr: Instance = %d, Finished!!\n", iIndex);
#endif

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit in astrosAIM aimDiscr, status = %d\n", status);

    AIM_FREE(tess);
    return status;
}


int aimTransfer(capsDiscr *discr, const char *dataName, int numPoint,
                int dataRank, double *dataVal, /*@unused@*/ char **units)
{

    /*! \page dataTransferAstros Astros Data Transfer
     *
     * The Astros AIM has the ability to transfer displacements and eigenvectors from the AIM and pressure
     * distributions to the AIM using the conservative and interpolative data transfer schemes in CAPS.
     *
     * \section dataFromAstros Data transfer from Astros (FieldOut)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements from the *.out file
     * </ul>
     *
     * <ul>
     *  <li> <B>"EigenVector_#"</B> </li> <br>
     *   Retrieves modal eigen-vectors from the *.out file, where "#" should be replaced by the
     *   corresponding mode number for the eigen-vector (e.g. EigenVector_3 would correspond to the third mode,
     *   while EigenVector_6 would be the sixth mode).
     * </ul>
     *
     * \section dataToAstros Data transfer to Astros (FieldIn)
     * <ul>
     *  <li> <B>"Pressure"</B> </li> <br>
     *  Writes appropriate load cards using the provided pressure distribution.
     * </ul>
     *
     */

    int status; // Function return status
    int i, j, dataPoint, bIndex; // Indexing

    char *extOUT = ".out";

    // FO6 data variables
    int numGridPoint = 0;
    int numEigenVector = 0;

    double **dataMatrix = NULL;
    aimStorage *astrosInstance;

    // Specific EigenVector to use
    int eigenVectorIndex = 0;

    // Variables used in global node mapping
    //int *storage;
    int globalNodeID;

    // Filename stuff
    char *filename = NULL;
    FILE *fp; // File pointer

#ifdef DEBUG
    printf(" astrosAIM/aimTransfer name = %s   npts = %d/%d!\n",
           dataName, numPoint, dataRank);
#endif
    astrosInstance = (aimStorage *) discr->instStore;

    //Get the appropriate parts of the tessellation to data
    //storage = (int *) discr->ptrm;
    //capsGroupList = &storage[1]; // List of boundary ID (attrMap) in the transfer

    if (strcasecmp(dataName, "Displacement") != 0 &&
        strncmp(dataName, "EigenVector", 11) != 0) {

        printf("Unrecognized data transfer variable - %s\n", dataName);
        return CAPS_NOTFOUND;
    }

    filename = (char *) EG_alloc((strlen(astrosInstance->projectName) +
                                  strlen(extOUT) + 1)*sizeof(char));
    if (filename == NULL) return EGADS_MALLOC;

    sprintf(filename, "%s%s", astrosInstance->projectName, extOUT);

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

            status = astros_readOUTDisplacement(fp,
                                                 -1,
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

            printf("Invalid rank for dataName \"%s\" - excepted a rank of 3!!!\n", dataName);
            status = CAPS_BADRANK;

        } else {

            status = astros_readOUTEigenVector(fp,
                                               &numEigenVector,
                                               &numGridPoint,
                                               &dataMatrix);
        }

        fclose(fp);

    } else {
        status = CAPS_NOTFOUND;
    }
    AIM_STATUS(discr->aInfo, status);
    if (dataMatrix == NULL) return CAPS_NULLVALUE;

    // Check EigenVector range
    if (strncmp(dataName, "EigenVector", 11) == 0)  {
        if (eigenVectorIndex > numEigenVector) {
            AIM_ERROR(discr->aInfo, "Only %d EigenVectors found but index %d requested!",
                   numEigenVector, eigenVectorIndex);
            status = CAPS_RANGEERR;
            goto cleanup;
        }

        if (eigenVectorIndex < 1) {
            AIM_ERROR(discr->aInfo, "For EigenVector_# notation, # must be >= 1, currently # = %d",
                   eigenVectorIndex);
            status = CAPS_RANGEERR;
            goto cleanup;
        }
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
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit in astrosAIM aimTransfer, status = %d\n", status);

    // Free data matrix
    if (dataMatrix != NULL) {
        j = 0;
        if      (strcasecmp(dataName, "Displacement") == 0) j = numGridPoint;
        else if (strncmp(dataName, "EigenVector", 11) == 0) j = numEigenVector;

        for (i = 0; i < j; i++) {
            AIM_FREE(dataMatrix[i]);
        }
        AIM_FREE(dataMatrix);
    }

    return status;
}


void aimFreeDiscrPtr(void *ptr)
{
    // Extra information to store into the discr void pointer - just a int array
    EG_free(ptr);
}


int aimLocateElement(capsDiscr *discr, double *params, double *param,
                     int *bIndex, int *eIndex, double *bary)
{
#ifdef DEBUG
    printf(" astrosAIM/aimLocateElement!\n");
#endif

    return aim_locateElement(discr, params, param, bIndex,eIndex, bary);
}


int aimInterpolation(capsDiscr *discr, const char *name,
                     int bIndex, int eIndex, double *bary,
                     int rank, double *data,
                     double *result)
{
#ifdef DEBUG
    printf(" astrosAIM/aimInterpolation  %s!\n", name);
#endif

    return aim_interpolation(discr, name, bIndex, eIndex, bary, rank, data,
                             result);

}


int aimInterpolateBar(capsDiscr *discr,  const char *name,
                      int bIndex, int eIndex,  double *bary,
                      int rank,  double *r_bar,
                      double *d_bar)
{
#ifdef DEBUG
    printf(" astrosAIM/aimInterpolateBar  %s!\n",  name);
#endif

    return aim_interpolateBar(discr, name, bIndex, eIndex, bary, rank, r_bar,
                              d_bar);
}


int aimIntegration(capsDiscr *discr,  const char *name,
                   int bIndex, int eIndex,  int rank,
                   double *data, double *result)
{
#ifdef DEBUG
    printf(" astrosAIM/aimIntegration  %s!\n", name);
#endif
    return aim_integration(discr, name, bIndex, eIndex, rank, data, result);
}


int aimIntegrateBar( capsDiscr *discr,  const char *name,
                     int bIndex, int eIndex,  int rank,
                     double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" astrosAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank, r_bar, d_bar);
}
