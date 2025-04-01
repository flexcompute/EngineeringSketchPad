/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Nastran AIM
 *
 *     Written by Dr. Ryan Durscher and Dr. Ed Alyanak AFRL/RQVC
 *
 *     This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

/*! \mainpage Introduction
 *
 * \section overviewNastran Nastran AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact (primarily
 * through input files) with the finite element structural solver Nastran \cite Nastran.
 *
 * Details on the use of units are outlined in \ref aimUnitsNastran.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsNastran and
 * \ref aimOutputsNastran and \ref attributeNastran, respectively.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferNastran
 *
 * \section nastranExamples Examples
 * Example problems using the Nastran AIM may be found at \ref examplesNastran .
 *  - \ref nastranSingleLoadEx
 *  - \ref nastranMultipleLoadEx
 *  - \ref nastranModalEx
 *  - \ref nastranOptEx
 *  - \ref nastranCompositeEx
 *  - \ref nastranCompOptimizationEx
 *
 * \section clearanceNastran Clearance Statement
 *  This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */

/*! \page attributeNastran AIM attributes
 * The following list of attributes are required for the Nastran AIM inside the geometry input.
 *
 * - <b> capsDiscipline</b> This attribute is a requirement if doing aeroelastic analysis within Nastran. capsDiscipline allows
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
 * - <b> capsIgnore</b> It is possible that there is a geometric body (or entity) that you do not want the Nastran AIM to pay attention to when creating
 * a finite element model. The capsIgnore attribute allows a body (or entity) to be in the geometry and ignored by the AIM.  For example,
 * because of limitations in OpenCASCADE a situation where two edges are overlapping may occur; capsIgnore allows the user to only
 *  pay attention to one of the overlapping edges.
 *
 * - <b> capsConnect</b> This is a name assigned to any geometric body where the user wishes to create
 * "fictitious" connections such as springs, dampers, and/or rigid body connections to. The user must manually specify
 * the connection between two <c>capsConnect</c> entities using the "Connect" tuple (see \ref aimInputsNastran).
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsConnect $springStart</c>.
 *
 * - <b> capsConnectLink</b> Similar to <c>capsConnect</c>, this is a name assigned to any geometric body where
 * the user wishes to create "fictitious" connections to. A connection is automatically made if a <c>capsConnectLink</c>
 * matches a <c>capsConnect</c> group. Again further specifics of the connection are input using the "Connect"
 * tuple (see \ref aimInputsNastran). Recall that a string in ESP starts with a $.
 * For example, attribute <c>capsConnectLink $springEnd</c>.
 *
 * - <b> capsResponse</b> This is a name assigned to any geometric body that will be used to define design sensitivity
 * responses for optimization. Specific information for the responses are input using the "Design_Response" tuple (see
 * \ref aimInputsNastran). Recall that a string in ESP starts with a $. For examples,
 * attribute <c>capsResponse $displacementNode</c>.
 *
 * - <b> capsBound </b> This is used to mark surfaces on the structural grid in which data transfer with an external
 * solver will take place. See \ref dataTransferNastran for additional details.
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
#else
#include <unistd.h>
#endif

#define MXCHAR     255

//#define DEBUG

enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Tess_Params,
  Edge_Point_Min,
  Edge_Point_Max,
  Quad_Mesh,
  Nastran,
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
  Design_Equation,
  Design_Table,
  Design_Response,
  Design_Equation_Response,
  Design_Opt_Param,
  Objective_Min_Max,
  Objective_Response_Type,
  Mass_Increment,
  VLM_Surface,
  VLM_Control,
  VLM_Camber_Twist,
  Visualize_Flutter,
  Support,
  Connect,
  Parameter,
  Model_Parameter,
  Aero_Reference,
  Mesh_Morph,
  Mesh,
  NUMINPUT = Mesh              /* Total number of inputs */
};

 enum aimOutputs
 {
   outEigenValue = 1,          /* index is 1-based */
   outEigenRadian,
   outEigenFrequency,
   outEigenGeneralMass,
   outEigenGeneralStiffness,
   outObjective,
   outObjectiveHistory,
   outMass,
   outCG,
   outIxx,
   outIyy,
   outIzz,
   outIxy,
   outIxz,
   outIyz,
   outI_Vector,
   outI_Lower,
   outI_Upper,
   outI_Tensor,
   outMassProp,
   outMassPropLink,
   outF06FileLink,
   NUMOUTPUT = outF06FileLink  /* Total number of outputs */
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

    // Attribute to reference map
    mapAttrToIndexStruct referenceMap;

    // Mesh holders
    int numMesh;
    meshStruct *feaMesh;

    // mass properties
    capsValue mass;
    capsValue Ixx, Iyy, Izz;
    capsValue Ixy, Ixz, Iyz;
    capsValue CG;

    feaMassPropStruct feaMassProp;
    feaSolFileStruct feaF06File;

} aimStorage;



static int initiate_aimStorage(void *aimInfo, aimStorage *nastranInstance)
{

    int i, status;

    // Set initial values for nastranInstance
    nastranInstance->projectName = NULL;

#define NUMVAR 8
    capsValue *values[NUMVAR];

/*
    // Check to make sure data transfer is ok
    nastranInstance->dataTransferCheck = (int) true;
*/

    status = initiate_feaUnitsStruct(&nastranInstance->units);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->attrMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to constraint index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->constraintMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to load index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->loadMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for transfer to index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->transferMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for connect to index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->connectMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for response to index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->responseMap);
    if (status != CAPS_SUCCESS) return status;

    // Container for reference to index map
    status = initiate_mapAttrToIndexStruct(&nastranInstance->referenceMap);
    if (status != CAPS_SUCCESS) return status;

    status = initiate_feaProblemStruct(&nastranInstance->feaProblem);
    if (status != CAPS_SUCCESS) return status;

    // Mesh holders
    nastranInstance->numMesh = 0;
    nastranInstance->feaMesh = NULL;

    // Mass properties
    values[0] = &nastranInstance->mass;
    values[1] = &nastranInstance->Ixx;
    values[2] = &nastranInstance->Iyy;
    values[3] = &nastranInstance->Izz;
    values[4] = &nastranInstance->Ixy;
    values[5] = &nastranInstance->Ixz;
    values[6] = &nastranInstance->Iyz;
    values[7] = &nastranInstance->CG;
#if NUMVAR != 8
#error "Developer error! Update values array!"
#endif
    for (i = 0; i < NUMVAR; i++) {
      aim_initValue(values[i]);
      values[i]->type = Double;
    }
#undef NUMVAR

    nastranInstance->CG.length = 3;
    AIM_ALLOC(nastranInstance->CG.vals.reals, 3, double, aimInfo, status);

    initiate_feaMassPropStruct(&nastranInstance->feaMassProp);

    initiate_feaSolFileStruct(&nastranInstance->feaF06File);

cleanup:
    return status;
}


static int destroy_aimStorage(aimStorage *nastranInstance)
{
    int status;
    int i;

    status = destroy_feaUnitsStruct(&nastranInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->attrMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to constraint index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->constraintMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Attribute to load index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->loadMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Transfer to index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->transferMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Connect to index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->connectMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Response to index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->responseMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Reference to index map
    status = destroy_mapAttrToIndexStruct(&nastranInstance->referenceMap);
    if (status != CAPS_SUCCESS) printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (nastranInstance->feaMesh != NULL) {

        for (i = 0; i < nastranInstance->numMesh; i++) {
            status = destroy_meshStruct(&nastranInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        AIM_FREE(nastranInstance->feaMesh);
    }

    nastranInstance->feaMesh = NULL;
    nastranInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&nastranInstance->feaProblem);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    // NULL projetName
    nastranInstance->projectName = NULL;

    // mass properties
    aim_freeValue(&nastranInstance->mass);
    aim_freeValue(&nastranInstance->Ixx);
    aim_freeValue(&nastranInstance->Iyy);
    aim_freeValue(&nastranInstance->Izz);
    aim_freeValue(&nastranInstance->Ixy);
    aim_freeValue(&nastranInstance->Ixz);
    aim_freeValue(&nastranInstance->Iyz);
    aim_freeValue(&nastranInstance->CG);

    destroy_feaMassPropStruct(&nastranInstance->feaMassProp);

    destroy_feaSolFileStruct(&nastranInstance->feaF06File);

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(void *aimInfo, aimStorage *nastranInstance)
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
  capsValue *TessParams;
  capsValue *EdgePoint_Min;
  capsValue *EdgePoint_Max;
  capsValue *QuadMesh;

  for (i = 0; i < nastranInstance->numMesh; i++) {
      remesh = remesh && (nastranInstance->feaMesh[i].egadsTess->oclass == EMPTY);
  }
  if (remesh == (int) false) return CAPS_SUCCESS;

  // retrieve or create the mesh from fea_createMesh
  status = aim_getValue(aimInfo, Tess_Params,    ANALYSISIN, &TessParams);
  AIM_STATUS(aimInfo, status);

  status = aim_getValue(aimInfo, Edge_Point_Min, ANALYSISIN, &EdgePoint_Min);
  AIM_STATUS(aimInfo, status);

  status = aim_getValue(aimInfo, Edge_Point_Max, ANALYSISIN, &EdgePoint_Max);
  AIM_STATUS(aimInfo, status);

  status = aim_getValue(aimInfo, Quad_Mesh,      ANALYSISIN, &QuadMesh);
  AIM_STATUS(aimInfo, status);

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

  if (QuadMesh != NULL) quadMesh = 2*QuadMesh->vals.integer;

  status = fea_createMesh(aimInfo,
                          tessParam,
                          edgePointMin,
                          edgePointMax,
                          quadMesh,
                          &nastranInstance->attrMap,
                          &nastranInstance->constraintMap,
                          &nastranInstance->loadMap,
                          &nastranInstance->transferMap,
                          &nastranInstance->connectMap,
                          &nastranInstance->responseMap,
                          &nastranInstance->referenceMap,
                          &nastranInstance->numMesh,
                          &nastranInstance->feaMesh,
                          &nastranInstance->feaProblem );
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


static int createVLMMesh(void *aimInfo, aimStorage *nastranInstance,
                         capsValue *aimInputs)
{

    int projectionMethod = (int) true;

    int status, status2; // Function return status

    int i, j, k, surfaceIndex = 0, sectionIndex, transferIndex; // Indexing

    // Bodies
    const char *intents;
    int   numBody; // Number of Bodies
    ego  *bodies;

    // Aeroelastic information
    int numVLMSurface = 0;
    vlmSurfaceStruct *vlmSurface = NULL;
    int numSpanWise;

    int numVLMControl = 0;
    vlmControlStruct *vlmControl = NULL;

    // Vector variables
    double A[3], B[3], C[3], D[3], P[3], p[3], N[3], n[3], d_proj[3];

    double *a, *b, *c, *d;
    double apbArea, apcArea, cpdArea, bpdArea, Area;

    feaMeshDataStruct *feaData;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    if (status != CAPS_SUCCESS) return status;

#ifdef DEBUG
    printf(" nastranAIM/createVLMMesh  nbody = %d!\n", numBody);
#endif

    if ((numBody <= 0) || (bodies == NULL)) {
#ifdef DEBUG
        printf(" nastranAIM/createVLMMesh No Bodies!\n");
#endif

        return CAPS_SOURCEERR;
    }

    // Get aerodynamic reference quantities
    status = fea_retrieveAeroRef(numBody, bodies,
                                 &nastranInstance->feaProblem.feaAeroRef);
    AIM_STATUS(aimInfo, status);

    // Cleanup Aero storage first
    if (nastranInstance->feaProblem.feaAero != NULL) {

        for (i = 0; i < nastranInstance->feaProblem.numAero; i++) {
            status = destroy_feaAeroStruct(&nastranInstance->feaProblem.feaAero[i]);
            AIM_STATUS(aimInfo, status);
        }

        AIM_FREE(nastranInstance->feaProblem.feaAero);
    }

    nastranInstance->feaProblem.numAero = 0;

    // Get VLM surface information
    if (aimInputs[VLM_Surface-1].nullVal != IsNull) {

        status = get_vlmSurface(aimInfo,
                                aimInputs[VLM_Surface-1].length,
                                aimInputs[VLM_Surface-1].vals.tuple,
                                &nastranInstance->attrMap,
                                0.0, // default Cspace
                                &numVLMSurface,
                                &vlmSurface);
        AIM_STATUS(aimInfo, status);

    } else {
        AIM_ERROR(aimInfo, "An analysis type of Aeroelastic set but no VLM_Surface tuple specified\n");
        status = CAPS_NOTFOUND;
        goto cleanup;
    }

    // Get VLM control surface information
    if (aimInputs[VLM_Control-1].nullVal == NotNull) {

        status = get_vlmControl(aimInfo,
                                aimInputs[VLM_Control-1].length,
                                aimInputs[VLM_Control-1].vals.tuple,
                                nastranInstance->units.length == NULL ? NULL : "degree",
                                &numVLMControl,
                                &vlmControl);

        AIM_STATUS(aimInfo, status);
    }

    printf("\nGetting FEA vortex lattice mesh\n");

    status = vlm_getSections(aimInfo,
                             numBody,
                             bodies,
                             "Aerodynamic",
                             nastranInstance->attrMap,
                             vlmGENERIC,
                             numVLMSurface,
                             &vlmSurface);
    AIM_STATUS(aimInfo, status);
    if (vlmSurface == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    for (i = 0; i < numVLMSurface; i++) {

        // Compute equal spacing
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
            AIM_ERROR(aimInfo, "Surface '%s' has less than two-sections!", vlmSurface[i].name);
            status = CAPS_BADVALUE;
            goto cleanup;
        }

        status = get_mapAttrToIndexIndex(&nastranInstance->transferMap,
                                         vlmSurface[i].name,
                                         &transferIndex);
        if (status == CAPS_NOTFOUND) {
            printf("\tA corresponding capsBound name not found for \"%s\". Surface will be ignored!\n",
                   vlmSurface[i].name);
            continue;
        } else AIM_STATUS(aimInfo, status);

        for (j = 0; j < vlmSurface[i].numSection-1; j++) {

            // Increment the number of Aero surfaces
            nastranInstance->feaProblem.numAero += 1;

            surfaceIndex = nastranInstance->feaProblem.numAero - 1;

            // Allocate
            AIM_REALL(nastranInstance->feaProblem.feaAero, nastranInstance->feaProblem.numAero, feaAeroStruct, aimInfo, status);

            // Initiate feaAeroStruct
            status = initiate_feaAeroStruct(&nastranInstance->feaProblem.feaAero[surfaceIndex]);
            AIM_STATUS(aimInfo, status);

            // Get surface Name - copy from original surface
            AIM_STRDUP(nastranInstance->feaProblem.feaAero[surfaceIndex].name, vlmSurface[i].name, aimInfo, status);

            // Get surface ID - Multiple by 10000 !!
            nastranInstance->feaProblem.feaAero[surfaceIndex].surfaceID =
                10000*nastranInstance->feaProblem.numAero;

            // ADD something for coordinate systems

            // Sections aren't necessarily stored in order coming out of vlm_getSections, however sectionIndex is!
            sectionIndex = vlmSurface[i].vlmSection[j].sectionIndex;

            // Populate vmlSurface structure
            nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.Cspace =
                     vlmSurface[i].Cspace;
            nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.Sspace =
                     vlmSurface[i].Sspace;

            // use the section span count for the sub-surface
            nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.NspanTotal =
                     vlmSurface[i].vlmSection[sectionIndex].Nspan;
            nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.Nchord     =
                     vlmSurface[i].Nchord;

            // Copy section information
            nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.numSection = 2;

            AIM_ALLOC(nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection, 2, vlmSectionStruct, aimInfo, status);

            for (k = 0; k < nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.numSection; k++) {

                // Add k to section indexing variable j to get j and j+1 during iterations

                // Sections aren't necessarily stored in order coming out of vlm_getSections, however sectionIndex is!
                sectionIndex = vlmSurface[i].vlmSection[j+k].sectionIndex;

                status = initiate_vlmSectionStruct(&nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k]);
                AIM_STATUS(aimInfo, status);

                // Copy the section data - This also copies the control data for the section
                status = copy_vlmSectionStruct(&vlmSurface[i].vlmSection[sectionIndex],
                                              &nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k]);
                AIM_STATUS(aimInfo, status);

                // Reset the sectionIndex that is keeping track of the section order.
                nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface.vlmSection[k].sectionIndex = k;
            }


            if (numVLMControl > 0) {
                AIM_NOTNULL(vlmControl, aimInfo, status);
                // transfer control surface data to sections
                status = get_ControlSurface(bodies,
                                            numVLMControl,
                                            vlmControl,
                                            &nastranInstance->feaProblem.feaAero[surfaceIndex].vlmSurface);
                AIM_STATUS(aimInfo, status);
            }
        }
    }

    if (nastranInstance->feaProblem.feaAero == NULL) {
        status = CAPS_NULLVALUE;
        goto cleanup;
    }

    // Determine which grid points are to be used for each spline
    for (i = 0; i < nastranInstance->feaProblem.numAero; i++) {

        // Debug
        //printf("\tDetermining grid points\n");

        // Get the transfer index for this surface - it has already been checked to make sure the name is in the
        // transfer index map
        status = get_mapAttrToIndexIndex(&nastranInstance->transferMap,
                                         nastranInstance->feaProblem.feaAero[i].name,
                                         &transferIndex);
        if (status != CAPS_SUCCESS ) goto cleanup;

        if (projectionMethod == (int) false) { // Look for attributes

            for (j = 0; j < nastranInstance->feaProblem.feaMesh.numNode; j++) {

                if (nastranInstance->feaProblem.feaMesh.node[j].analysisType ==
                    MeshStructure) {
                    feaData = (feaMeshDataStruct *)
                       nastranInstance->feaProblem.feaMesh.node[j].analysisData;
                } else {
                    continue;
                }

                if (feaData->transferIndex != transferIndex) continue;
                if (feaData->transferIndex == CAPSMAGIC) continue;


                nastranInstance->feaProblem.feaAero[i].numGridID += 1;
                k = nastranInstance->feaProblem.feaAero[i].numGridID;

                AIM_REALL(nastranInstance->feaProblem.feaAero[i].gridIDSet, k, int, aimInfo, status);

                nastranInstance->feaProblem.feaAero[i].gridIDSet[k-1] =
                    nastranInstance->feaProblem.feaMesh.node[j].nodeID;
            }

        } else { // Projection method


            /*
             *   n = A X B Create a normal vector/ plane between A and B
             *
             *   d_proj = C - (C * n)*n/ ||n||^2 , projection of point d on plane created by AxB
             *
             *   p = D - (D * n)*n/ ||n||^2 , projection of point p on plane created by AxB
             *
             *                            (section 2)
             *                     LE(c)---------------->TE(d)
             *   Grid Point       -^                   ^ -|
             *           |^      -            -         - |
             *           | -     A      -   C          - d_proj
             *           |  D   -    -                 -
             *           |   - - -     (section 1)    -
             *           p    LE(a)----------B------->TE(b)
             */

            // Vector between section 2 and 1
            a = nastranInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[0].xyzLE;
            b = nastranInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[0].xyzTE;
            c = nastranInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[1].xyzLE;
            d = nastranInstance->feaProblem.feaAero[i].vlmSurface.vlmSection[1].xyzTE;

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

            for (j = 0; j < nastranInstance->feaProblem.feaMesh.numNode; j++) {

                if (nastranInstance->feaProblem.feaMesh.node[j].analysisType ==
                    MeshStructure) {
                    feaData = (feaMeshDataStruct *)
                       nastranInstance->feaProblem.feaMesh.node[j].analysisData;
                } else {
                    continue;
                }

                if (feaData->transferIndex != transferIndex) continue;
                if (feaData->transferIndex == CAPSMAGIC) continue;

                D[0] = nastranInstance->feaProblem.feaMesh.node[j].xyz[0] - a[0];

                D[1] = nastranInstance->feaProblem.feaMesh.node[j].xyz[1] - a[1];

                D[2] = nastranInstance->feaProblem.feaMesh.node[j].xyz[2] - a[2];

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
                //	printf("Point index %d\n", j);
                //	printf("\tAreas don't match - %f vs %f!!\n", Area, (apbArea + apcArea + cpdArea + bpdArea));
                //	printf("\tPoint - %f %f %f\n", D[0] + a[0], D[1] + a[1], D[2] + a[2]);
                //	printf("\tPoint projection - %f %f %f\n", p[0], p[1], p[2]);
                //}

                if (fabs(apbArea + apcArea + cpdArea + bpdArea - Area) > 1E-5) 	continue;

                nastranInstance->feaProblem.feaAero[i].numGridID += 1;

                AIM_REALL(nastranInstance->feaProblem.feaAero[i].gridIDSet,
                          nastranInstance->feaProblem.feaAero[i].numGridID, int, aimInfo, status);

                nastranInstance->feaProblem.feaAero[i].gridIDSet[
                    nastranInstance->feaProblem.feaAero[i].numGridID-1] =
                    nastranInstance->feaProblem.feaMesh.node[j].nodeID;
            }
        }

        printf("\tSurface %d: Number of points found for aero-spline = %d\n",
               i+1, nastranInstance->feaProblem.feaAero[i].numGridID );
    }

    status = CAPS_SUCCESS;

cleanup:

    if (status != CAPS_SUCCESS)
        printf("\tcreateVLMMesh status = %d\n", status);

    if (vlmSurface != NULL) {

        for (i = 0; i < numVLMSurface; i++) {
            status2 = destroy_vlmSurfaceStruct(&vlmSurface[i]);
            if (status2 != CAPS_SUCCESS)
                printf("\tdestroy_vlmSurfaceStruct status = %d\n", status2);
        }
    }

    AIM_FREE(vlmSurface);
    numVLMSurface = 0;

    // Destroy Control
    if (vlmControl != NULL) {
        for (i = 0; i < numVLMControl; i++) {
            (void) destroy_vlmControlStruct(&vlmControl[i]);
        }
        AIM_FREE(vlmControl);
    }
    numVLMControl = 0;

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
    const char *keyWord;
    char *keyValue = NULL;
    double real = 1.0;
    feaUnitsStruct *units=NULL;

    aimStorage *nastranInstance=NULL;

#ifdef DEBUG
    printf("nastranAIM/aimInitialize   instance = %d!\n", inst);
#endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUTPUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 5;

    /* specify the name of each field variable */
    AIM_ALLOC(strs, *nFields, char *, aimInfo, status);

    strs[0]  = EG_strdup("Displacement");
    strs[1]  = EG_strdup("EigenVector");
    strs[2]  = EG_strdup("EigenVector_#");
    strs[3]  = EG_strdup("Pressure");
    strs[4]  = EG_strdup("Temperature");
    for (i = 0; i < *nFields; i++)
      if (strs[i] == NULL) { status = EGADS_MALLOC; goto cleanup; }
    *fnames  = strs;

    /* specify the dimension of each field variable */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);
    ints[0]  = 3;
    ints[1]  = 3;
    ints[2]  = 3;
    ints[3]  = 1;
    ints[4]  = 1;
    *franks  = ints;
    ints = NULL;

    /* specify if a field is an input field or output field */
    AIM_ALLOC(ints, *nFields, int, aimInfo, status);

    ints[0]  = FieldOut;
    ints[1]  = FieldOut;
    ints[2]  = FieldOut;
    ints[3]  = FieldIn;
    ints[4]  = FieldIn;
    *fInOut  = ints;
    ints = NULL;

    // Allocate nastranInstance
    AIM_ALLOC(nastranInstance, 1, aimStorage, aimInfo, status);
    *instStore = nastranInstance;

    // Initialize instance storage
    status = initiate_aimStorage(aimInfo, nastranInstance);
    AIM_STATUS(aimInfo, status);


    /*! \page aimUnitsNastran AIM Units
     *  A unit system may be optionally specified during AIM instance initiation. If
     *  a unit system is provided, all AIM  input values which have associated units must be specified as well.
     *  If no unit system is used, AIM inputs, which otherwise would require units, will be assumed
     *  unit consistent. A unit system may be specified via a JSON string dictionary for example (using pyCAPS):
     *  \code{.py}
     *  nastran = capsProblem.analysis.create(aim="nastranAIM", unitSys = {"mass": "kg", "length": "inch", "time": "hour", "temperature" : "Kelvin"})
     *  \endcode
     */
    if (unitSys != NULL) {
      units = &nastranInstance->units;

      // Do we have a json string?
      if (strncmp( unitSys, "{", 1) != 0) {
        AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
        return CAPS_BADVALUE;
      }

      /*! \page aimUnitsNastran
       *  \section jsonStringNastran JSON String Dictionary
       *  The key arguments of the dictionary are described in the following:
       *
       *  <ul>
       *  <li> <B>mass = "None"</B> </li> <br>
       *  Mass units - e.g. "kilogram", "k", "slug", ...
       *  </ul>
       */
      keyWord = "mass";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->mass = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->mass, &real, "kg", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->mass, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsNastran
       *  <ul>
       *  <li> <B>length = "None"</B> </li> <br>
       *  Length units - e.g. "meter", "m", "inch", "in", "mile", ...
       *  </ul>
       */
      keyWord = "length";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->length = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->length, &real, "m", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->length, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsNastran
       *  <ul>
       *  <li> <B>time = "None"</B> </li> <br>
       *  Time units - e.g. "s", "second", "hour", "day", "year", ...
       *  </ul>
       */
      keyWord = "time";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->time = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->time, &real, "s", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->time, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      /*! \page aimUnitsNastran
       *  <ul>
       *  <li> <B>temperature = "None"</B> </li> <br>
       *  Temperature units - e.g. "C", "F", "K", "R", ...
       *  </ul>
       */
      keyWord = "temperature";
      status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
      if (status == CAPS_SUCCESS) {
        units->temperature = string_removeQuotation(keyValue);
        AIM_FREE(keyValue);
        real = 1;
        status = aim_convert(aimInfo, 1, units->temperature, &real, "K", &real);
        AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->temperature, keyWord);
      } else {
        AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      // create derived units from the base units
      status = fea_feaDerivedUnits(aimInfo, units);
      AIM_STATUS(aimInfo, status);
    }

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

    /*! \page aimInputsNastran AIM Inputs
     * The following list outlines the Nastran inputs along with their default value available
     * through the AIM interface. Unless noted these values will be not be linked to
     * any parent AIMs with variables of the same name.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" nastranAIM/aimInputs index = %d!\n", index);
#endif

    *ainame = NULL;

    // Nastran Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("nastran_CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsNastran
         * - <B> Proj_Name = "nastran_CAPS"</B> <br>
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
        AIM_ALLOC(defval->vals.reals, defval->nrow, double, aimInfo, status);
        defval->vals.reals[0] = 0.025;
        defval->vals.reals[1] = 0.001;
        defval->vals.reals[2] = 15.00;

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
         * - <B> Edge_Point_Max = 50</B> <br>
         * Maximum number of points on an edge including end points to use when creating a surface mesh (min 2).
         */

    } else if (index == Quad_Mesh) {
        *ainame               = EG_strdup("Quad_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = (int) false;

        /*! \page aimInputsNastran
         * - <B> Quad_Mesh = False</B> <br>
         * Create a quadratic mesh on four edge faces when creating the boundary element model.
         */

    } else if (index == Nastran) {
        *ainame               = EG_strdup("Nastran");
        defval->type          = String;
        defval->dim           = Vector;
        defval->lfixed        = Change;
        defval->nullVal       = IsNull;
        defval->vals.string   = NULL;

        /*! \page aimInputsNastran
         * - <B> Nastran = NULL</B> <br>
         * List of strings to specify NASTRAN statements <br>
         *   The NASTRAN statement is used to specify values for certain Executive System operational parameters. <br>
         *   These parameters are also called system cells. The NASTRAN statement is used for exceptional circumstances <br>
         *   and is therefore not needed in most runs. The NASTRAN statement may also be specified in the <br>
         *   runtime configuration (RC) files at the system, user, and job level as described in the .
         */

    } else if (index == Property) {
        *ainame              = EG_strdup("Property");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
         * - <B> Analysis = NULL</B> <br>
         * Analysis tuple used to input analysis/case information for the model, see \ref feaAnalysis for additional details.
         */
    } else if (index == Analysis_Type) {
        *ainame              = EG_strdup("Analysis_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Modal");
        defval->lfixed       = Change;

        /*! \page aimInputsNastran
         * - <B> Analysis_Type = "Modal"</B> <br>
         * Type of analysis to generate files for, options include "Modal", "Static", "AeroelasticTrim", "AeroelasticFlutter", and "Optimization".
         * Note: "Aeroelastic" and "StaticOpt" are still supported and refer to "AeroelasticTrim" and "Optimization".
         */
    } else if (index == File_Format) {
        *ainame              = EG_strdup("File_Format");
        defval->type         = String;
        defval->vals.string  = EG_strdup("Small"); // Small, Large, Free
        defval->lfixed       = Change;

        /*! \page aimInputsNastran
         * - <B> File_Format = "Small"</B> <br>
         * Formatting type for the bulk file. Options: "Small", "Large", "Free".
         */

    } else if (index == Mesh_File_Format) {
        *ainame              = EG_strdup("Mesh_File_Format");
        defval->type         = String;
        defval->vals.string  = EG_strdup("Free"); // Small, Large, Free
        defval->lfixed       = Change;

        /*! \page aimInputsNastran
         * - <B> Mesh_File_Format = "Small"</B> <br>
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
         * - <B> Design_Opt_Param = NULL</B> <br>
         * The design optimization parameter tuple used to input parameters used in design optimization, see \ref feaDesignOptParam for additional details.
         */

    } else if (index == Objective_Min_Max) {
        *ainame              = EG_strdup("Objective_Min_Max");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Max"); // Max, Min
        defval->lfixed       = Change;

        /*! \page aimInputsNastran
         * - <B> Objective_Min_Max = "Max"</B> <br>
         * Maximize or minimize the design objective during an optimization. Option: "Max" or "Min".
         */

    } else if (index == Objective_Response_Type) {
        *ainame              = EG_strdup("Objective_Response_Type");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Weight"); // Weight
        defval->lfixed       = Change;

        /*! \page aimInputsNastran
         * - <B> Objective_Response_Type = "Weight"</B> <br>
         * Object response type (see Nastran manual).
         */

    } else if (index == Mass_Increment) {
        *ainame              = EG_strdup("Mass_Increment");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsNastran
         * - <B> Mass_Increment = NULL </B> <br>
         * FEA Mass Increments (see Nastran manual).
         */

    } else if (index == VLM_Surface) {
        *ainame              = EG_strdup("VLM_Surface");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsNastran
         * - <B>VLM_Surface = NULL </B> <br>
         * Vortex lattice method tuple input,  see \ref vlmSurface for additional details.
         */

    }  else if (index == VLM_Control) {
        *ainame              = EG_strdup("VLM_Control");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsNastran
         * - <B>VLM_Control = NULL </B> <br>
         * Vortex lattice method control surface tuple input, see \ref vlmControl for additional details.
         */
    } else if (index == Support) {
        *ainame              = EG_strdup("Support");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
         * - <B> Parameter = NULL</B> <br>
         * Parameter tuple used to define PARAM entries. Note, entries are output exactly as inputed, that is, if the PARAM entry
         * requires an integer entry the user must input an integer!
         */
    } else if (index == Model_Parameter) {
        *ainame              = EG_strdup("Model_Parameter");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsNastran
         * - <B> Model_Parameter = NULL</B> <br>
         * Model_Parameter tuple used to define MDLPRM entries. Note, entries are output exactly as inputed, that is, if the MDLPRM entry
         * requires an integer entry the user must input an integer!
         */

    } else if (index == Aero_Reference) {
        *ainame              = EG_strdup("Aero_Reference");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;
        defval->dim          = Vector;

        /*! \page aimInputsNastran
         * - <B> Aero_Reference = NULL</B> <br>
         * A JSON dictionary used to define aerodynamic reference parameters. see \ref feaAeroReference for additional details
         */

    } else if (index == VLM_Camber_Twist) {
        *ainame              = EG_strdup("VLM_Camber_Twist");
        defval->type         = Boolean;
        defval->nullVal      = NotNull;
        defval->vals.integer = (int) true;

        /*! \page aimInputsNastran
         * - <B> VLM_Camber_Twist = True</B> <br>
         * Apply camber & twist to VLM sections. Option: True or False.
         */

    } else if (index == Visualize_Flutter) {
        *ainame              = EG_strdup("Visualize_Flutter");
        defval->type         = Boolean;
        defval->nullVal      = NotNull;
        defval->vals.integer = (int) false;

        /*! \page aimInputsNastran
        * - <B> Visualize_Flutter = False</B> <br>
        * Determines if flutter cards are written for visualization. Option: True or False.
        */

    } else if (index == Mesh_Morph) {
        *ainame              = EG_strdup("Mesh_Morph");
        defval->type         = Boolean;
        defval->lfixed       = Fixed;
        defval->vals.integer = (int) false;
        defval->dim          = Scalar;
        defval->nullVal      = NotNull;

        /*! \page aimInputsNastran
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

        /*! \page aimInputsNastran
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

    int found = (int)false;
    int i;

    const char *analysisType = NULL;
    const char *discipline = NULL;

    const char *intents;
    int   numBody; // Number of Bodies
    ego  *bodies;

    //vlmControlStruct *dummyControlSurface=NULL;

    aimStorage *nastranInstance;

    nastranInstance = (aimStorage *) instStore;
    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Get project name
    nastranInstance->projectName = aimInputs[Proj_Name-1].vals.string;

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    // Get FEA mesh if we don't already have one
    status = checkAndCreateMesh(aimInfo, nastranInstance);
    AIM_STATUS(aimInfo, status);

    // Get Aeroelastic mesh
    if ((strcasecmp(analysisType, "Aeroelastic") == 0) ||
        (strcasecmp(analysisType, "AeroelasticTrim") == 0) ||
        (strcasecmp(analysisType, "AeroelasticFlutter") == 0) ||
        (strcasecmp(analysisType, "Gust") == 0) ||
        (strcasecmp(analysisType, "Optimization") == 0)) {

        found = (int) true;
        if (strcasecmp(analysisType, "Optimization") == 0) { // Is this aeroelastic optimization?
            found = (int) false;

            status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
            AIM_STATUS(aimInfo, status);

            for (i = 0; i < numBody; i++) {
                status = retrieve_CAPSDisciplineAttr(bodies[i], &discipline);
                if ((status == CAPS_SUCCESS) && (discipline != NULL)) {
                    if (strcasecmp(discipline, "Aerodynamic") == 0) {
                        found = (int) true; // We at least have aerodynamic bodies
                        break;
                    }
                }
            }
        }

        if (found == (int) true) {
            status = createVLMMesh(aimInfo, nastranInstance, aimInputs);
            AIM_STATUS(aimInfo, status);
        }
    }


    // Set aero reference parameters
    if (aimInputs[Aero_Reference-1].nullVal == NotNull) {
        status = fea_getAeroReference(aimInfo,
                                      aimInputs[Aero_Reference-1].length,
                                      aimInputs[Aero_Reference-1].vals.tuple,
                                      &nastranInstance->referenceMap,
                                      &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Aero_Reference value is NULL - No aero reference parameters set\n");

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
                                 &nastranInstance->units,
                                 &nastranInstance->feaProblem.numMaterial,
                                 &nastranInstance->feaProblem.feaMaterial);
        AIM_STATUS(aimInfo, status);
    } else printf("Material tuple is NULL - No materials set\n");

    // Set property properties
    if (aimInputs[Property-1].nullVal == NotNull) {
        status = fea_getProperty(aimInfo,
                                 aimInputs[Property-1].length,
                                 aimInputs[Property-1].vals.tuple,
                                 &nastranInstance->attrMap,
                                 &nastranInstance->units,
                                 &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);


        // Assign element "subtypes" based on properties set
        status = fea_assignElementSubType(nastranInstance->feaProblem.numProperty,
                                          nastranInstance->feaProblem.feaProperty,
                                          &nastranInstance->feaProblem.feaMesh);
        AIM_STATUS(aimInfo, status);

    } else printf("Property tuple is NULL - No properties set\n");

    // Set constraint properties
    if (aimInputs[Constraint-1].nullVal == NotNull) {
        status = fea_getConstraint(aimInfo,
                                   aimInputs[Constraint-1].length,
                                   aimInputs[Constraint-1].vals.tuple,
                                   &nastranInstance->constraintMap,
                                   &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Constraint tuple is NULL - No constraints applied\n");

    // Set support properties
    if (aimInputs[Support-1].nullVal == NotNull) {
        status = fea_getSupport(aimInputs[Support-1].length,
                                aimInputs[Support-1].vals.tuple,
                                &nastranInstance->constraintMap,
                                &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Support tuple is NULL - No supports applied\n");

    // Set connection properties
    if (aimInputs[Connect-1].nullVal == NotNull) {
        status = fea_getConnection(aimInfo,
                                   aimInputs[Connect-1].length,
                                   aimInputs[Connect-1].vals.tuple,
                                   &nastranInstance->connectMap,
                                   &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Connect tuple is NULL - Using defaults\n");

    // Set load properties
    if (aimInputs[Load-1].nullVal == NotNull) {
        status = fea_getLoad(aimInfo,
                             aimInputs[Load-1].length,
                             aimInputs[Load-1].vals.tuple,
                             &nastranInstance->loadMap,
                             &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Load tuple is NULL - No loads applied\n");

    // Set design variables
    if (aimInputs[Design_Variable-1].nullVal == NotNull) {
        status = fea_getDesignVariable(aimInfo,
                                       (int)true,
                                       aimInputs[Design_Variable-1].length,
                                       aimInputs[Design_Variable-1].vals.tuple,
                                       aimInputs[Design_Variable_Relation-1].length,
                                       aimInputs[Design_Variable_Relation-1].vals.tuple,
                                       &nastranInstance->attrMap,
                                       &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Variable tuple is NULL - No design variables applied\n");

    // Set design constraints
    if (aimInputs[Design_Constraint-1].nullVal == NotNull) {
        status = fea_getDesignConstraint(aimInfo,
                                         aimInputs[Design_Constraint-1].length,
                                         aimInputs[Design_Constraint-1].vals.tuple,
                                         &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Constraint tuple is NULL - No design constraints applied\n");

    // Set design equations
    if (aimInputs[Design_Equation-1].nullVal == NotNull) {
        status = fea_getDesignEquation(aimInputs[Design_Equation-1].length,
                                       aimInputs[Design_Equation-1].vals.tuple,
                                       &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Equation tuple is NULL - No design equations applied\n");

    // Set design table constants
    if (aimInputs[Design_Table-1].nullVal == NotNull) {
        status = fea_getDesignTable(aimInputs[Design_Table-1].length,
                                    aimInputs[Design_Table-1].vals.tuple,
                                    &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Table tuple is NULL - No design table constants applied\n");

    // Set design optimization parameters
    if (aimInputs[Design_Opt_Param-1].nullVal == NotNull) {
        status = fea_getDesignOptParam(aimInputs[Design_Opt_Param-1].length,
                                       aimInputs[Design_Opt_Param-1].vals.tuple,
                                       &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Opt_Param tuple is NULL - No design optimization parameters applied\n");

    // Set design constraints
    if (aimInputs[Mass_Increment-1].nullVal == NotNull) {
        status = fea_getMassIncrement(aimInputs[Mass_Increment-1].length,
                                      aimInputs[Mass_Increment-1].vals.tuple,
                                      &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Mass_Increment tuple is NULL - No mass increments applied\n");

    // Set design responses
    if (aimInputs[Design_Response-1].nullVal == NotNull) {
        status = fea_getDesignResponse(aimInfo,
                                       aimInputs[Design_Response-1].length,
                                       aimInputs[Design_Response-1].vals.tuple,
                                       &nastranInstance->responseMap,
                                       &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Response tuple is NULL - No design responses applied\n");

    // Set design equation responses
    if (aimInputs[Design_Equation_Response-1].nullVal == NotNull) {
        status = fea_getDesignEquationResponse(aimInputs[Design_Equation_Response-1].length,
                                               aimInputs[Design_Equation_Response-1].vals.tuple,
                                               &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status);
    } else printf("Design_Equation_Response tuple is NULL - No design equation responses applied\n");

    // Set analysis settings
    if (aimInputs[Analysix-1].nullVal == NotNull) {
        status = fea_getAnalysis(aimInfo,
                                 aimInputs[Analysix-1].length,
                                 aimInputs[Analysix-1].vals.tuple,
                                 &nastranInstance->feaProblem);
        AIM_STATUS(aimInfo, status); // It ok to not have an analysis tuple
    } else {
        printf("Analysis tuple is NULL\n");
        status = fea_createDefaultAnalysis(aimInfo, &nastranInstance->feaProblem, analysisType);
        AIM_STATUS(aimInfo, status);
    }


    // Set file format type
    if        (strcasecmp(aimInputs[File_Format-1].vals.string, "Small") == 0) {
      nastranInstance->feaProblem.feaFileFormat.fileType = SmallField;
    } else if (strcasecmp(aimInputs[File_Format-1].vals.string, "Large") == 0) {
      nastranInstance->feaProblem.feaFileFormat.fileType = LargeField;
    } else if (strcasecmp(aimInputs[File_Format-1].vals.string, "Free") == 0)  {
      nastranInstance->feaProblem.feaFileFormat.fileType = FreeField;
    } else {
      printf("Unrecognized \"File_Format\", valid choices are [Small, Large, or Free]. Reverting to default\n");
    }

    // Set grid file format type
    if        (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Small") == 0) {
        nastranInstance->feaProblem.feaFileFormat.gridFileType = SmallField;
    } else if (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Large") == 0) {
        nastranInstance->feaProblem.feaFileFormat.gridFileType = LargeField;
    } else if (strcasecmp(aimInputs[Mesh_File_Format-1].vals.string, "Free") == 0)  {
        nastranInstance->feaProblem.feaFileFormat.gridFileType = FreeField;
    } else {
        printf("Unrecognized \"Mesh_File_Format\", valid choices are [Small, Large, or Free]. Reverting to default\n");
    }

    status = CAPS_SUCCESS;
cleanup:
    return status;
}

int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{

    int i, j, k, l; // Indexing

    int status; // Status return

    int found; // Boolean operator

    int *tempIntegerArray = NULL; // Temporary array to store a list of integers

    // Analysis information
    int haveSubAeroelasticTrim = (int) false;
    int haveSubAeroelasticFlutter = (int) false;

    // Mesh length scaling
    const char *lengthUnits=NULL;
    double scaleFactor=1;

    // Optimization Information
    char *objectiveMinMax = NULL, *objectiveResp = NULL;

    // Aeroelastic Information
    int numAEStatSurf = 0;
    //char **aeStatSurf = NULL;

    // File format information
    char *tempString = NULL, *delimiter = NULL;
    char fieldString[16];
    char *NASTRAN=NULL;

    // File IO
    char *filename = NULL; // Output file name
    FILE *fp = NULL; // Output file pointer

    feaLoadStruct *feaLoad = NULL;
    int numThermalLoad=0;

    int numSetID;
    int tempID, *setID = NULL;

    const char *analysisType = NULL;

    vlmControlStruct *dummyControlSurface=NULL;


    const aimStorage *nastranInstance;

    nastranInstance = (const aimStorage *) instStore;
    AIM_NOTNULL(aimInputs, aimInfo, status);

    if ( aimInputs[Mesh_Morph-1].vals.integer == (int) true &&
         aimInputs[Mesh-1].nullVal == NotNull) { // If we are mighty morphing
        // store the current mesh for future iterations
        status = aim_storeMeshRef(aimInfo, (aimMeshRef *) aimInputs[Mesh-1].vals.AIMptr, NULL);
        AIM_STATUS(aimInfo, status);
    }

    // Analysis type
    analysisType = aimInputs[Analysis_Type-1].vals.string;

    // Write Nastran Mesh
    AIM_ALLOC(filename, (MXCHAR +1), char, aimInfo, status);
    strcpy(filename, nastranInstance->projectName);

    if (nastranInstance->feaProblem.numLoad > 0) {
        AIM_ALLOC(feaLoad, nastranInstance->feaProblem.numLoad, feaLoadStruct, aimInfo, status);
        for (i = 0; i < nastranInstance->feaProblem.numLoad; i++) initiate_feaLoadStruct(&feaLoad[i]);
        for (i = 0; i < nastranInstance->feaProblem.numLoad; i++) {
            status = copy_feaLoadStruct(aimInfo, &nastranInstance->feaProblem.feaLoad[i], &feaLoad[i]);
            AIM_STATUS(aimInfo, status);

            if (feaLoad[i].loadType == PressureExternal) {

                // Transfer external pressures from the AIM discrObj
                status = fea_transferExternalPressure(aimInfo,
                                                      &nastranInstance->feaProblem.feaMesh,
                                                      &nastranInstance->units,
                                                      &feaLoad[i]);
                AIM_STATUS(aimInfo, status);
            } else if (feaLoad[i].loadType == ThermalExternal) {

                // Transfer external temperature from the AIM discrObj
                status = fea_transferExternalTemperature(aimInfo,
                                                         &feaLoad[i]);
                AIM_STATUS(aimInfo, status);
            }
        }
    }

    if (nastranInstance->units.length != NULL) {
      status = aim_capsLength(aimInfo, &lengthUnits);
      AIM_NOTFOUND(aimInfo, status);
      if (status == CAPS_NOTFOUND) {
        AIM_ERROR(aimInfo, "capsLength attribute must be specified!");
        goto cleanup;
      }
      AIM_NOTNULL(lengthUnits, aimInfo, status);

      status = aim_convert(aimInfo, 1, lengthUnits, &scaleFactor, nastranInstance->units.length, &scaleFactor);
      AIM_STATUS(aimInfo, status);
    }

    status = mesh_writeNASTRAN(aimInfo,
                               filename,
                               1,
                               &nastranInstance->feaProblem.feaMesh,
                               nastranInstance->feaProblem.numProperty,
                               nastranInstance->feaProblem.feaProperty,
                               nastranInstance->feaProblem.feaFileFormat.gridFileType,
                               scaleFactor);
    AIM_STATUS(aimInfo, status);

    // Write Nastran subElement types not supported by mesh_writeNASTRAN
    strcat(filename, ".bdf");
    fp = aim_fopen(aimInfo, filename, "a");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file: %s", filename);
        status = CAPS_IOERR;
        goto cleanup;
    }
    AIM_FREE(filename);

    printf("Writing subElement types (if any) - appending mesh file\n");
    status = nastran_writeSubElementCard(aimInfo, fp,
                                         &nastranInstance->feaProblem.feaMesh,
                                         nastranInstance->feaProblem.numProperty,
                                         nastranInstance->feaProblem.feaProperty,
                                         &nastranInstance->feaProblem.feaFileFormat);
    AIM_STATUS(aimInfo, status);

    // Connections
    for (i = 0; i < nastranInstance->feaProblem.numConnect; i++) {

        if (i == 0) {
            printf("Writing connection cards - appending mesh file\n");
        }

        status = nastran_writeConnectionCard(fp,
                                             &nastranInstance->feaProblem.feaConnect[i],
                                             &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Multiple Mass Conditions Set
    for (i = 0; i < nastranInstance->feaProblem.numAnalysis; i++) {
        if (i == 0) printf("\tWriting mass increment set\n");

        if (nastranInstance->feaProblem.feaAnalysis[i].numMassIncrement > 0) {

          status = nastran_writeMassIncrementSet(fp,
                                                 i+1,
                                                 nastranInstance->feaProblem.feaAnalysis[i].numMassIncrement,
                                                 nastranInstance->feaProblem.feaAnalysis[i].massIncrementSetID,
                                                 nastranInstance->feaProblem.feaMassIncrement,
                                                 &nastranInstance->feaProblem.feaFileFormat);

            if (status != CAPS_SUCCESS) goto cleanup;
        }
    }

    if (fp != NULL) fclose(fp);
    fp = NULL;

    // Write nastran input file
    AIM_ALLOC(filename, (MXCHAR +1), char, aimInfo, status);
    strcpy(filename, nastranInstance->projectName);
    strcat(filename, ".dat");

    printf("\nWriting Nastran instruction file....\n");
    fp = aim_fopen(aimInfo, filename, "w");
    if (fp == NULL) {
        AIM_ERROR(aimInfo, "Unable to open file: %s\n", filename);
        AIM_FREE(filename);
        status = CAPS_IOERR;
        goto cleanup;
    }
    AIM_FREE(filename);

    // define file format delimiter type
    if (nastranInstance->feaProblem.feaFileFormat.fileType == FreeField) {
        delimiter = ",";
    } else {
        delimiter = " ";
    }

    //////////////// Executive control ////////////////
    if (aimInputs[Nastran-1].nullVal == NotNull) {
      NASTRAN = aimInputs[Nastran-1].vals.string;
      for (i = 0; i < aimInputs[Nastran-1].length; i++) {
        fprintf(fp, "NASTRAN %s\n", NASTRAN);
        NASTRAN += strlen(NASTRAN)+1;
      }
    }

    fprintf(fp, "ID CAPS generated Problem FOR Nastran\n");

    // Analysis type
    if     (strcasecmp(analysisType, "Modal")               == 0) fprintf(fp, "SOL 103\n");
    else if(strcasecmp(analysisType, "Static")              == 0) fprintf(fp, "SOL 101\n");
    else if(strcasecmp(analysisType, "Craig-Bampton")       == 0) fprintf(fp, "SOL 31\n");
    else if(strcasecmp(analysisType, "StaticOpt")           == 0) fprintf(fp, "SOL 200\n");
    else if(strcasecmp(analysisType, "Optimization")        == 0) fprintf(fp, "SOL 200\n");
    else if(strcasecmp(analysisType, "Aeroelastic")         == 0) fprintf(fp, "SOL 144\n");
    else if(strcasecmp(analysisType, "AeroelasticTrim")     == 0) fprintf(fp, "SOL 144\n");
    else if(strcasecmp(analysisType, "AeroelasticFlutter")  == 0) fprintf(fp, "SOL 145\n");
    else if(strcasecmp(analysisType, "Gust")                == 0) fprintf(fp, "SOL 146\n");
    else {
        AIM_ERROR(aimInfo, "Unrecognized \"Analysis_Type\", %s", analysisType);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    fprintf(fp, "CEND\n\n");

    if (nastranInstance->feaProblem.feaMesh.numNode> 10000) fprintf(fp, "LINE=%d\n", nastranInstance->feaProblem.feaMesh.numNode*10);
    else                                                           fprintf(fp, "LINE=10000\n");


    // Set up the case information
     if (nastranInstance->feaProblem.numAnalysis == 0) {
         printf("Error: No analyses in the feaProblem! (this shouldn't be possible)\n");
         status = CAPS_BADVALUE;
         goto cleanup;
     }

    //////////////// Case control ////////////////

    // Set flutter visualization
    nastranInstance->feaProblem.feaAnalysis->visualFlutter = aimInputs[Visualize_Flutter -1].vals.integer;

    // Write output request information
    if (nastranInstance->feaProblem.feaAnalysis->visualFlutter == (int) true) {
      printf("\tFLFACT VELOCITIES WILL BE NEGATIVE, FLUTTER VISUALIZATION IS TRUE\n");

      fprintf(fp, "DISP (PRINT,PUNCH,PHASE) = ALL\n"); // Output all displacements with phase
    }
    else{
      fprintf(fp, "DISP (PRINT,PUNCH) = ALL\n"); // Output all displacements
    }

    fprintf(fp, "STRE (PRINT,PUNCH) = ALL\n"); // Output all stress

    fprintf(fp, "STRA (PRINT,PUNCH) = ALL\n"); // Output all strain

    // Design objective information, SOL 200 only
    if ((strcasecmp(analysisType, "StaticOpt") == 0) ||
        (strcasecmp(analysisType, "Optimization") == 0)) {

        objectiveMinMax = aimInputs[Objective_Min_Max-1].vals.string;
        if     (strcasecmp(objectiveMinMax, "Min") == 0) fprintf(fp, "DESOBJ(MIN) = 1\n");
        else if(strcasecmp(objectiveMinMax, "Max") == 0) fprintf(fp, "DESOBJ(MAX) = 1\n");
        else {
            printf("Unrecognized \"Objective_Min_Max\", %s, defaulting to \"Min\"\n", objectiveMinMax);
            //objectiveMinMax = "Min";
            fprintf(fp, "DESOBJ(MIN) = 1\n");
        }

    }

    // Write sub-case information if multiple analysis tuples were provide - will always have at least 1
    for (i = 0; i < nastranInstance->feaProblem.numAnalysis; i++) {
        fprintf(fp, "SUBCASE %d\n", i+1);
        fprintf(fp, "    LABEL = %s\n", nastranInstance->feaProblem.feaAnalysis[i].name);

        if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == Static) {
            fprintf(fp,"    ANALYSIS = STATICS\n");
        } else if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == Modal) {
            fprintf(fp,"    ANALYSIS = MODES\n");
        } else if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim) {
            fprintf(fp,"    ANALYSIS = SAERO\n");
            haveSubAeroelasticTrim = (int) true;
        } else if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {
            fprintf(fp,"    ANALYSIS = FLUTTER\n");
            haveSubAeroelasticFlutter = (int) true;
        } else if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == Optimization) {
            printf("\t *** WARNING :: INPUT TO ANALYSIS CASE INPUT analysisType should NOT be Optimization or StaticOpt - Defaulting to Static\n");
            fprintf(fp,"    ANALYSIS = STATICS\n");
            }

        // Write support for sub-case
        if (nastranInstance->feaProblem.feaAnalysis[i].numSupport != 0) {
            if (nastranInstance->feaProblem.feaAnalysis[i].numSupport > 1) {
                printf("\tWARNING: More than 1 support is not supported at this time for sub-cases!\n");

            } else {
                fprintf(fp, "    SUPORT1 = %d\n", nastranInstance->feaProblem.feaAnalysis[i].supportSetID[0]);
            }
        }

        // Write constraint for sub-case
        if (nastranInstance->feaProblem.numConstraint != 0) {
            fprintf(fp, "    SPC = %d\n", nastranInstance->feaProblem.numConstraint+i+1); //TODO - change to i+1 to just i
        }

        // Issue some warnings regarding constraints if necessary
        if (nastranInstance->feaProblem.feaAnalysis[i].numConstraint == 0 && nastranInstance->feaProblem.numConstraint !=0) {

            printf("\tWarning: No constraints specified for case %s, assuming "
                    "all constraints are applied!!!!\n", nastranInstance->feaProblem.feaAnalysis[i].name);

        } else if (nastranInstance->feaProblem.numConstraint == 0) {

            printf("\tWarning: No constraints specified for case %s!!!!\n", nastranInstance->feaProblem.feaAnalysis[i].name);
        }

  //        // Write MPC for sub-case - currently only supported when we have RBE2 elements - see above for unification - TODO - investigate
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

        if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == Modal ||
            nastranInstance->feaProblem.feaAnalysis[i].analysisType == Gust) {
            fprintf(fp,"    METHOD = %d\n", nastranInstance->feaProblem.feaAnalysis[i].analysisID);
        }

        if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {
            fprintf(fp,"    METHOD = %d\n", nastranInstance->feaProblem.feaAnalysis[i].analysisID);
            fprintf(fp,"    FMETHOD = %d\n", 100+nastranInstance->feaProblem.feaAnalysis[i].analysisID);
        }

        if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim) {
            fprintf(fp,"    TRIM = %d\n", nastranInstance->feaProblem.feaAnalysis[i].analysisID);
        }

        if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticTrim ||
            nastranInstance->feaProblem.feaAnalysis[i].analysisType == AeroelasticFlutter) {

            if (nastranInstance->feaProblem.feaAnalysis[i].aeroSymmetryXY != NULL) {

                if(strcmp("SYM",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"    AESYMXY = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYM",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"    AESYMXY = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYM",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"    AESYMXY = %s\n","ASYMMETRIC");
                } else if(strcmp("SYMMETRIC",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"    AESYMXY = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYMMETRIC",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"    AESYMXY = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYMMETRIC",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY) == 0) {
                    fprintf(fp,"    AESYMXY = %s\n","ASYMMETRIC");
                } else {
                    printf("\t*** Warning *** aeroSymmetryXY Input %s to nastranAIM not understood!\n",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXY );
                }
            }

            if (nastranInstance->feaProblem.feaAnalysis[i].aeroSymmetryXZ != NULL) {

                if(strcmp("SYM",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"    AESYMXZ = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYM",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"    AESYMXZ = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYM",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"    AESYMXZ = %s\n","ASYMMETRIC");
                } else if(strcmp("SYMMETRIC",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"    AESYMXZ = %s\n","SYMMETRIC");
                } else if(strcmp("ANTISYMMETRIC",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"    AESYMXZ = %s\n","ANTISYMMETRIC");
                } else if(strcmp("ASYMMETRIC",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ) == 0) {
                    fprintf(fp,"    AESYMXZ = %s\n","ASYMMETRIC");
                } else {
                    printf("\t*** Warning *** aeroSymmetryXZ Input %s to nastranAIM not understood!\n",nastranInstance->feaProblem.feaAnalysis->aeroSymmetryXZ );
                }

            }
        }

        // Issue some warnings regarding loads if necessary
        if (nastranInstance->feaProblem.feaAnalysis[i].numLoad == 0 && nastranInstance->feaProblem.numLoad !=0) {
            printf("\tWarning: No loads specified for case %s, assuming "
                    "all loads are applied!!!!\n", nastranInstance->feaProblem.feaAnalysis[i].name);
        } else if (nastranInstance->feaProblem.numLoad == 0) {
            printf("\tWarning: No loads specified for case %s!!!!\n", nastranInstance->feaProblem.feaAnalysis[i].name);
        }

        // Write loads for sub-case
        if (feaLoad != NULL) {

            found = (int) false;

            for (k = 0; k < nastranInstance->feaProblem.numLoad; k++) {

                if (nastranInstance->feaProblem.feaAnalysis[i].numLoad != 0) { // if loads specified in analysis

                    for (j = 0; j < nastranInstance->feaProblem.feaAnalysis[i].numLoad; j++) { // See if the load is in the loadSet

                        if (feaLoad[k].loadID == nastranInstance->feaProblem.feaAnalysis[i].loadSetID[j]) break;
                    }

                    if (j >= nastranInstance->feaProblem.feaAnalysis[i].numLoad) continue; // If it isn't in the loadSet move on
                } else {
                    //pass
                }

                if (feaLoad[k].loadType == Thermal && numThermalLoad == 0) {

                    fprintf(fp, "    Temperature = %d\n", feaLoad[k].loadID);
                    numThermalLoad += 1;
                    if (numThermalLoad > 1) {
                        printf("More than 1 Thermal load found - nastranAIM does NOT currently doesn't support multiple thermal loads in a given case!\n");
                    }

                    continue;
                }

                found = (int) true;
            }

            if (found == (int) true) {
                fprintf(fp, "    LOAD = %d\n", nastranInstance->feaProblem.numLoad+i+1);
            }
        }

        //if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == Optimization) {
            // Write objective function
            //fprintf(fp, "    DESOBJ(%s) = %d\n", nastranInstance->feaProblem.feaAnalysis[i].objectiveMinMax,
            //                                   nastranInstance->feaProblem.feaAnalysis[i].analysisID);
            // Write optimization constraints
        if (nastranInstance->feaProblem.feaAnalysis[i].numDesignConstraint != 0) {
            fprintf(fp, "    DESSUB = %d\n", nastranInstance->feaProblem.numDesignConstraint+i+1);
        }
        //}

        // Write response spanning for sub-case
        if (nastranInstance->feaProblem.feaAnalysis[i].numDesignResponse != 0) {

            numSetID = nastranInstance->feaProblem.feaAnalysis[i].numDesignResponse;
            // setID = nastranInstance->feaProblem.feaAnalysis[i].designResponseSetID;
            AIM_ALLOC(setID, numSetID, int, aimInfo, status);

            for (j = 0; j < numSetID; j++) {
                tempID = nastranInstance->feaProblem.feaAnalysis[i].designResponseSetID[j];
                setID[j] = tempID + 100000;
            }

            tempID = i+1;
            status = nastran_writeSetCard(fp, tempID, numSetID, setID);

            AIM_FREE(setID);

            AIM_STATUS(aimInfo, status);
            fprintf(fp, "    DRSPAN = %d\n", tempID);
        }

        // Add MASSET information to sub-case case control
        if (nastranInstance->feaProblem.feaAnalysis[i].numMassIncrement != 0) {
            fprintf(fp, "\tMASSSET = %d\n", 100+i+1);
        }
    }

//
//    // Check thermal load - currently only a single thermal load is supported - also it only works for the global level - no subcase
//    found = (int) false;
//    for (i = 0; i < nastranInstance->feaProblem.numLoad; i++) {
//
//        if (feaLoad[i].loadType != Thermal) continue;
//
//        if (found == (int) true) {
//            printf("More than 1 Thermal load found - nastranAIM does NOT currently doesn't support multiple thermal loads!\n");
//        }
//
//        found = (int) true;
//
//        fprintf(fp, "TEMPERATURE = %d\n", feaLoad[i].loadID);
//    }
//
//    // Design objective information, SOL 200 only
//    if ((strcasecmp(analysisType, "StaticOpt") == 0) || (strcasecmp(analysisType, "Optimization") == 0)) {
//
//        objectiveMinMax = aimInputs[ObjectiveMinMax-1].vals.string;
//        if     (strcasecmp(objectiveMinMax, "Min") == 0) fprintf(fp, "DESOBJ(MIN) = 1\n");
//        else if(strcasecmp(objectiveMinMax, "Max") == 0) fprintf(fp, "DESOBJ(MAX) = 1\n");
//        else {
//            printf("Unrecognized \"ObjectiveMinMax\", %s, defaulting to \"Min\"\n", objectiveMinMax);
//            objectiveMinMax = "Min";
//            fprintf(fp, "DESOBJ(MIN) = 1\n");
//        }
//
//    }
//
//    // Modal analysis - only
//    // If modal - we are only going to use the first analysis structure we come across that has its type as modal
//    if (strcasecmp(analysisType, "Modal") == 0) {
//
//        // Look through analysis structures for a modal one
//        found = (int) false;
//        for (i = 0; i < nastranInstance->feaProblem.numAnalysis; i++) {
//            if (nastranInstance->feaProblem.feaAnalysis[i].analysisType == Modal) {
//                found = (int) true;
//                break;
//            }
//        }
//
//        // Write out analysis ID if a modal analysis structure was found
//        if (found == (int) true)  {
//            fprintf(fp, "METHOD = %d\n", nastranInstance->feaProblem.feaAnalysis[i].analysisID);
//        } else {
//            printf("Warning: No eigenvalue analysis information specified in \"Analysis\" tuple, through "
//            		"AIM input \"Analysis_Type\" is set to \"Modal\"!!!!\n");
//            status = CAPS_NOTFOUND;
//            goto cleanup;
//        }
//
//        // Write support for sub-case
//        if (nastranInstance->feaProblem.feaAnalysis[i].numSupport != 0) {
//            if (nastranInstance->feaProblem.feaAnalysis[i].numSupport > 1) {
//                printf("\tWARNING: More than 1 support is not supported at this time for sub-cases!\n");
//
//            } else {
//                fprintf(fp, "SUPORT1 = %d\n", nastranInstance->feaProblem.feaAnalysis[i].supportSetID[0]);
//            }
//        }
//
//        // Write constraint for sub-case
//        if (nastranInstance->feaProblem.numConstraint != 0) {
//            fprintf(fp, "SPC = %d\n", nastranInstance->feaProblem.numConstraint+i+1);
//        }
//
//        // Issue some warnings regarding constraints if necessary
//        if (nastranInstance->feaProblem.feaAnalysis[i].numConstraint == 0 && nastranInstance->feaProblem.numConstraint !=0) {
//            printf("\tWarning: No constraints specified for modal case %s, assuming "
//                    "all constraints are applied!!!!\n", nastranInstance->feaProblem.feaAnalysis[i].name);
//        } else if (nastranInstance->feaProblem.numConstraint == 0) {
//            printf("\tWarning: No constraints specified for modal case %s!!!!\n", nastranInstance->feaProblem.feaAnalysis[i].name);
//        }
//
//    }

//        } else { // If no sub-cases
//
//            if (nastranInstance->feaProblem.numSupport != 0) {
//                if (nastranInstance->feaProblem.numSupport > 1) {
//                    printf("\tWARNING: More than 1 support is not supported at this time for a given case!\n");
//                } else {
//                    fprintf(fp, "SUPORT1 = %d ", nastranInstance->feaProblem.numSupport+1);
//                }
//            }
//
//            // Write constraint information
//            if (nastranInstance->feaProblem.numConstraint != 0) {
//                fprintf(fp, "SPC = %d\n", nastranInstance->feaProblem.numConstraint+1);
//            } else {
//                printf("\tWarning: No constraints specified for job!!!!\n");
//            }
//
//            // Write load card
//            if (nastranInstance->feaProblem.numLoad > 0) {
//                fprintf(fp, "LOAD = %d\n", nastranInstance->feaProblem.numLoad+1);
//            } else {
//                printf("\tWarning: No loads specified for job!!!!\n");
//            }
//
//            // What about an objective function if no analysis tuple? Do we need to add a new capsValue?
//
//            // Write design constraints
//            if (nastranInstance->feaProblem.numDesignConstraint != 0) {
//                fprintf(fp, "    DESSUB = %d\n", nastranInstance->feaProblem.numDesignConstraint+1);
//            }
//        }
//    }


    //////////////// Bulk data ////////////////
    fprintf(fp, "\nBEGIN BULK\n");
    fprintf(fp, "$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---10--|\n");
    //PRINT PARAM ENTRIES IN BULK DATA

    if (aimInputs[Parameter-1].nullVal == NotNull) {
        for (i = 0; i < aimInputs[Parameter-1].length; i++) {

            tempString = string_removeQuotation(aimInputs[Parameter-1].vals.tuple[i].value);
            AIM_NOTNULL(tempString, aimInfo, status);

            fprintf(fp, "PARAM, %s, %s\n", aimInputs[Parameter-1].vals.tuple[i].name,
                                           tempString);
            AIM_FREE(tempString);
        }
    }

    fprintf(fp, "PARAM, %s\n", "POST, -1\n"); // Output OP2 file
    fprintf(fp, "PARAM, %s\n", "GRDPNT, 0\n"); // Write out mass properties to F06 file

    // Add MDLPRM for HDF5 Output
    if (aimInputs[Model_Parameter-1].nullVal == NotNull) {
        for (i = 0; i < aimInputs[Model_Parameter-1].length; i++) {

            tempString = string_removeQuotation(aimInputs[Model_Parameter-1].vals.tuple[i].value);
            AIM_NOTNULL(tempString, aimInfo, status);

            fprintf(fp, "MDLPRM, %s, %s\n", aimInputs[Model_Parameter-1].vals.tuple[i].name, tempString);
            AIM_FREE(tempString);
        }
    }

    // Mass increment
    for (i = 0; i < nastranInstance->feaProblem.numAnalysis; i++) {
        if (i == 0) printf("\tWriting mass increment cards\n");

        if (nastranInstance->feaProblem.feaAnalysis[i].numMassIncrement > 0) {

            status = nastran_writeMasssetCard(fp,
                                              i+1,
                                              &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);
        }
    }

    // Turn off auto SPC
    //fprintf(fp, "%-8s %7s %7s\n", "PARAM", "AUTOSPC", "N");

    // Optimization Objective Response Response, SOL 200 only
    if (strcasecmp(analysisType, "StaticOpt") == 0 || strcasecmp(analysisType, "Optimization") == 0) {

        objectiveResp = aimInputs[Objective_Response_Type-1].vals.string;
        if     (strcasecmp(objectiveResp, "Weight") == 0) objectiveResp = "WEIGHT";
        else {
            printf("\tUnrecognized \"Objective_Response_Type\", %s, defaulting to \"Weight\"\n", objectiveResp);
            objectiveResp = "WEIGHT";
        }

        fprintf(fp,"%-8s", "DRESP1");

        status = convert_integerToString(1, 7, 1, fieldString);
        AIM_STATUS(aimInfo, status);
        fprintf(fp, "%s%s", delimiter, fieldString);

        fprintf(fp, "%s%7s", delimiter, objectiveResp);
        fprintf(fp, "%s%7s", delimiter, objectiveResp);


        fprintf(fp, "\n");
    }

    // Write AEROS, AESTAT and AESURF cards
    if (strcasecmp(analysisType, "AeroelasticFlutter") == 0 ||
        strcasecmp(analysisType, "Gust") == 0 ||
        haveSubAeroelasticFlutter == (int) true) {

        printf("\tWriting aero card\n");
        status = nastran_writeAEROCard(fp,
                                       &nastranInstance->feaProblem.feaAeroRef,
                                       &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Write AEROS, AESTAT and AESURF cards
    if (strcasecmp(analysisType, "Aeroelastic") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
        haveSubAeroelasticTrim == (int) true) {

        printf("\tWriting aeros card\n");
        status = nastran_writeAEROSCard(fp,
                                        &nastranInstance->feaProblem.feaAeroRef,
                                        &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);

        numAEStatSurf = 0;
        for (i = 0; i < nastranInstance->feaProblem.numAnalysis; i++) {

            if (nastranInstance->feaProblem.feaAnalysis[i].analysisType != AeroelasticTrim) continue;

            if (i == 0) printf("\tWriting aestat cards\n");

            // Loop over rigid variables
            for (j = 0; j < nastranInstance->feaProblem.feaAnalysis[i].numRigidVariable; j++) {

                found = (int) false;

                // Loop over previous rigid variables
                for (k = 0; k < i; k++) {
                    for (l = 0; l < nastranInstance->feaProblem.feaAnalysis[k].numRigidVariable; l++) {

                        // If current rigid variable was previous written - mark as found
                        if (strcmp(nastranInstance->feaProblem.feaAnalysis[i].rigidVariable[j],
                                   nastranInstance->feaProblem.feaAnalysis[k].rigidVariable[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // Make sure constraint isn't already in rigid variables too!
                for (k = 0; k < i; k++) {
                    for (l = 0; l < nastranInstance->feaProblem.feaAnalysis[k].numRigidVariable; l++) {

                        // If current rigid constraint was previous written - mark as found
                        if (strcmp(nastranInstance->feaProblem.feaAnalysis[i].rigidConstraint[j],
                                   nastranInstance->feaProblem.feaAnalysis[k].rigidVariable[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // Check if is control surface (will be defined by AESURF)
                status = fea_findControlSurfaceByName(&nastranInstance->feaProblem,
                                                      nastranInstance->feaProblem.feaAnalysis[i].rigidVariable[j],
                                                      &dummyControlSurface);
                if (status == CAPS_SUCCESS) {
                    continue;
                }

                // If not write out an aestat card
                numAEStatSurf += 1;

                fprintf(fp,"%-8s", "AESTAT");

                status = convert_integerToString(numAEStatSurf, 7, 1, fieldString);
                AIM_STATUS(aimInfo, status);
                fprintf(fp, "%s%s", delimiter, fieldString);

                fprintf(fp, "%s%7s\n", delimiter, nastranInstance->feaProblem.feaAnalysis[i].rigidVariable[j]);
            }

            // Loop over rigid Constraints
            for (j = 0; j < nastranInstance->feaProblem.feaAnalysis[i].numRigidConstraint; j++) {

                found = (int) false;

                // Loop over previous rigid constraints
                for (k = 0; k < i; k++) {
                    for (l = 0; l < nastranInstance->feaProblem.feaAnalysis[k].numRigidConstraint; l++) {

                        // If current rigid constraint was previous written - mark as found
                        if (strcmp(nastranInstance->feaProblem.feaAnalysis[i].rigidConstraint[j],
                                   nastranInstance->feaProblem.feaAnalysis[k].rigidConstraint[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // Make sure constraint isn't already in rigid variables too!
                for (k = 0; k < i; k++) {
                    for (l = 0; l < nastranInstance->feaProblem.feaAnalysis[k].numRigidVariable; l++) {

                        // If current rigid constraint was previous written - mark as found
                        if (strcmp(nastranInstance->feaProblem.feaAnalysis[i].rigidConstraint[j],
                                   nastranInstance->feaProblem.feaAnalysis[k].rigidVariable[l]) == 0) {
                            found = (int) true;
                            break;
                        }
                    }
                }

                // If already found continue
                if (found == (int) true) continue;

                // Check if is control surface (will be defined by AESURF)
                status = fea_findControlSurfaceByName(&nastranInstance->feaProblem,
                                                      nastranInstance->feaProblem.feaAnalysis[i].rigidConstraint[j],
                                                      &dummyControlSurface);
                if (status == CAPS_SUCCESS) {
                    continue;
                }

                // If not write out an aestat card
                numAEStatSurf += 1;

                fprintf(fp,"%-8s", "AESTAT");

                status = convert_integerToString(numAEStatSurf, 7, 1, fieldString);
                AIM_STATUS(aimInfo, status);
                fprintf(fp, "%s%s", delimiter, fieldString);

                fprintf(fp, "%s%7s\n", delimiter, nastranInstance->feaProblem.feaAnalysis[i].rigidConstraint[j]);
            }
        }

        fprintf(fp,"\n");
    }

    // Analysis Cards - Eigenvalue and design objective included, as well as combined load, constraint, and design constraints
    for (i = 0; i < nastranInstance->feaProblem.numAnalysis; i++) {

        if (i == 0) printf("\tWriting analysis cards, %d\n", i);

        // TODO: FIX THIS
        nastranInstance->feaProblem.feaAnalysis[i].visualFlutter = aimInputs[Visualize_Flutter -1].vals.integer;

        status = nastran_writeAnalysisCard(fp,
                                           &nastranInstance->feaProblem.feaAnalysis[i],
                                           &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);

        if (nastranInstance->feaProblem.feaAnalysis[i].numLoad != 0) {
            AIM_NOTNULL(feaLoad, aimInfo, status);

            // Create a temporary list of load IDs
            AIM_ALLOC(tempIntegerArray, nastranInstance->feaProblem.feaAnalysis[i].numLoad, int, aimInfo, status);

            k = 0;
            for (j = 0; j <  nastranInstance->feaProblem.feaAnalysis[i].numLoad; j++) {
                for (l = 0; l < nastranInstance->feaProblem.numLoad; l++) {
                    if (nastranInstance->feaProblem.feaAnalysis[i].loadSetID[j] == feaLoad[l].loadID) break;
                }
                if (feaLoad[l].loadType == Thermal) continue;
                tempIntegerArray[j] = feaLoad[l].loadID;
                k += 1;
            }
            if (k == 0) { AIM_FREE(tempIntegerArray); continue; }

            AIM_REALL(tempIntegerArray, k, int, aimInfo, status);

            // Write combined load card
            printf("\tWriting load ADD cards\n");
            status = nastran_writeLoadADDCard(fp,
                                               nastranInstance->feaProblem.numLoad+i+1,
                                               k,
                                               tempIntegerArray,
                                               feaLoad,
                                               &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

            // Free temporary load ID list
            AIM_FREE(tempIntegerArray);

        } else { // If no loads for an individual analysis are specified assume that all loads should be applied

            if (feaLoad != NULL) {

                // Create a temporary list of load IDs
                AIM_ALLOC(tempIntegerArray, nastranInstance->feaProblem.numLoad, int, aimInfo, status);

                k = 0;
                for (j = 0; j < nastranInstance->feaProblem.numLoad; j++) {
                    if (feaLoad[j].loadType == Gravity) continue;
                    tempIntegerArray[k] = feaLoad[j].loadID;
                    k += 1;
                }

                AIM_REALL(tempIntegerArray, k, int, aimInfo, status);

                //TOOO: eliminate load add card?
                // Write combined load card
                printf("\tWriting load ADD cards\n");
                status = nastran_writeLoadADDCard(fp,
                                                  nastranInstance->feaProblem.numLoad+i+1,
                                                  k,
                                                  tempIntegerArray,
                                                  feaLoad,
                                                  &nastranInstance->feaProblem.feaFileFormat);
                AIM_STATUS(aimInfo, status);

                // Free temporary load ID list
                AIM_FREE(tempIntegerArray);
            }

        }

        if (nastranInstance->feaProblem.feaAnalysis[i].numConstraint != 0) {
            // Write combined constraint card
            printf("\tWriting constraint ADD cards\n");
            status = nastran_writeConstraintADDCard(fp,
                                                    nastranInstance->feaProblem.numConstraint+i+1,
                                                    nastranInstance->feaProblem.feaAnalysis[i].numConstraint,
                                                    nastranInstance->feaProblem.feaAnalysis[i].constraintSetID,
                                                    &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

        } else { // If no constraints for an individual analysis are specified assume that all constraints should be applied

            if (nastranInstance->feaProblem.numConstraint != 0) {

                printf("\tWriting combined constraint cards\n");

                // Create a temporary list of constraint IDs
                AIM_ALLOC(tempIntegerArray, nastranInstance->feaProblem.numConstraint, int, aimInfo, status);

                for (j = 0; j < nastranInstance->feaProblem.numConstraint; j++) {
                    tempIntegerArray[j] = nastranInstance->feaProblem.feaConstraint[j].constraintID;
                }

                // Write combined constraint card
                status = nastran_writeConstraintADDCard(fp,
                                                        nastranInstance->feaProblem.numConstraint+i+1,
                                                        nastranInstance->feaProblem.numConstraint,
                                                        tempIntegerArray,
                                                        &nastranInstance->feaProblem.feaFileFormat);
                AIM_STATUS(aimInfo, status);

                // Free temporary constraint ID list
                AIM_FREE(tempIntegerArray);
            }
        }

        if (nastranInstance->feaProblem.feaAnalysis[i].numDesignConstraint != 0) {

            // Write combined design constraint card
            printf("\tWriting design constraint ADD cards\n");
            status = nastran_writeDesignConstraintADDCard(fp,
                                                          nastranInstance->feaProblem.numDesignConstraint+i+1,
                                                          nastranInstance->feaProblem.feaAnalysis[i].numDesignConstraint,
                                                          nastranInstance->feaProblem.feaAnalysis[i].designConstraintSetID,
                                                          &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

        } else { // If no design constraints for an individual analysis are specified assume that all design constraints should be applied

            if (nastranInstance->feaProblem.numDesignConstraint != 0) {

                // Create a temporary list of design constraint IDs
                AIM_ALLOC(tempIntegerArray, nastranInstance->feaProblem.numDesignConstraint, int, aimInfo, status);

                for (j = 0; j < nastranInstance->feaProblem.numDesignConstraint; j++) {
                    tempIntegerArray[j] = nastranInstance->feaProblem.feaDesignConstraint[j].designConstraintID;
                }

                // Write combined design constraint card
                printf("\tWriting design constraint ADD cards\n");
                status = nastran_writeDesignConstraintADDCard(fp,
                                                              nastranInstance->feaProblem.numDesignConstraint+i+1,
                                                              nastranInstance->feaProblem.numDesignConstraint,
                                                              tempIntegerArray,
                                                              &nastranInstance->feaProblem.feaFileFormat);
                // Free temporary design constraint ID list
                AIM_FREE(tempIntegerArray);
                AIM_STATUS(aimInfo, status);
            }

        }
    }


    // Loads
    for (i = 0; i < nastranInstance->feaProblem.numLoad; i++) {
        AIM_NOTNULL(feaLoad, aimInfo, status);

        if (i == 0) printf("\tWriting load cards\n");

        status = nastran_writeLoadCard(fp,
                                       &feaLoad[i],
                                       &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Constraints
    for (i = 0; i < nastranInstance->feaProblem.numConstraint; i++) {

        if (i == 0) printf("\tWriting constraint cards\n");

        status = nastran_writeConstraintCard(fp,
                                             &nastranInstance->feaProblem.feaConstraint[i],
                                             &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Supports
    for (i = 0; i < nastranInstance->feaProblem.numSupport; i++) {

        if (i == 0) printf("\tWriting support cards\n");
        j = (int) true;
        status = nastran_writeSupportCard(fp,
                                          &nastranInstance->feaProblem.feaSupport[i],
                                          &nastranInstance->feaProblem.feaFileFormat,
                                          &j);
        AIM_STATUS(aimInfo, status);
    }


    // Materials
    for (i = 0; i < nastranInstance->feaProblem.numMaterial; i++) {

        if (i == 0) printf("\tWriting material cards\n");

        status = nastran_writeMaterialCard(fp,
                                           &nastranInstance->feaProblem.feaMaterial[i],
                                           &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Properties
    for (i = 0; i < nastranInstance->feaProblem.numProperty; i++) {

        if (i == 0) printf("\tWriting property cards\n");

        status = nastran_writePropertyCard(fp,
                                           &nastranInstance->feaProblem.feaProperty[i],
                                           &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Coordinate systems
    // printf("DEBUG: Number of coord systems: %d\n", nastranInstance->feaProblem.numCoordSystem);
    for (i = 0; i < nastranInstance->feaProblem.numCoordSystem; i++) {

        if (i == 0) printf("\tWriting coordinate system cards\n");

        status = nastran_writeCoordinateSystemCard(fp,
                                                   &nastranInstance->feaProblem.feaCoordSystem[i],
                                                   &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design variables
    for( i = 0; i < nastranInstance->feaProblem.numDesignVariable; i++) {

        if (i == 0) printf("\tWriting design variable cards\n");

        status = nastran_writeDesignVariableCard(aimInfo, fp,
                                                 &nastranInstance->feaProblem.feaDesignVariable[i],
                                                 &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design variable relations
    for( i = 0; i < nastranInstance->feaProblem.numDesignVariableRelation; i++) {

        if (i == 0) printf("\tWriting design variable relation cards\n");

        status = nastran_writeDesignVariableRelationCard(aimInfo,
                                                         fp,
                                                         &nastranInstance->feaProblem.feaDesignVariableRelation[i],
                                                         &nastranInstance->feaProblem,
                                                         &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design constraints
    for( i = 0; i < nastranInstance->feaProblem.numDesignConstraint; i++) {

        if (i == 0) printf("\tWriting design constraints and responses cards\n");

        status = nastran_writeDesignConstraintCard(fp,
                                                   &nastranInstance->feaProblem.feaDesignConstraint[i],
                                                   &nastranInstance->feaProblem,
                                                   &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design equations
    for( i = 0; i < nastranInstance->feaProblem.numEquation; i++) {

        if (i == 0) printf("\tWriting design equation cards\n");

        status = nastran_writeDesignEquationCard(fp,
                                                 &nastranInstance->feaProblem.feaEquation[i],
                                                 &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design table constants
    if (nastranInstance->feaProblem.feaDesignTable.numConstant > 0)
        printf("\tWriting design table card\n");
    status = nastran_writeDesignTableCard(fp,
                                          &nastranInstance->feaProblem.feaDesignTable,
                                          &nastranInstance->feaProblem.feaFileFormat);
    AIM_STATUS(aimInfo, status);

    // Optimization - design optimization parameters
    if (nastranInstance->feaProblem.feaDesignOptParam.numParam > 0)
        printf("\tWriting design optimization parameters card\n");
    status = nastran_writeDesignOptParamCard(fp,
                                             &nastranInstance->feaProblem.feaDesignOptParam,
                                             &nastranInstance->feaProblem.feaFileFormat);
    AIM_STATUS(aimInfo, status);

    // Optimization - design responses
    for( i = 0; i < nastranInstance->feaProblem.numDesignResponse; i++) {

        if (i == 0) printf("\tWriting design response cards\n");

        status = nastran_writeDesignResponseCard(fp,
                                                 &nastranInstance->feaProblem.feaDesignResponse[i],
                                                 &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Optimization - design equation responses
    for( i = 0; i < nastranInstance->feaProblem.numEquationResponse; i++) {

        if (i == 0) printf("\tWriting design equation response cards\n");

        status = nastran_writeDesignEquationResponseCard(fp,
                                                 &nastranInstance->feaProblem.feaEquationResponse[i],
                                                 &nastranInstance->feaProblem,
                                                 &nastranInstance->feaProblem.feaFileFormat);
        AIM_STATUS(aimInfo, status);
    }

    // Aeroelastic
    if (strcasecmp(analysisType, "Aeroelastic") == 0 ||
        strcasecmp(analysisType, "AeroelasticTrim") == 0 ||
        strcasecmp(analysisType, "AeroelasticFlutter") == 0 ||
        strcasecmp(analysisType, "Gust") == 0 ||
        haveSubAeroelasticTrim == (int) true ||
        haveSubAeroelasticFlutter == (int) true ) {

        printf("\tWriting aeroelastic cards\n");
        for (i = 0; i < nastranInstance->feaProblem.numAero; i++){
            status = nastran_writeCAeroCard(fp,
                                            &nastranInstance->feaProblem.feaAero[i],
                                            &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

            status = nastran_writeAeroSplineCard(fp,
                                                 &nastranInstance->feaProblem.feaAero[i],
                                                 &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);

            status = nastran_writeSet1Card(fp,
                                           &nastranInstance->feaProblem.feaAero[i],
                                           &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);
        }

        if (aimInputs[VLM_Camber_Twist-1].vals.integer == (int) true ) {


            status = nastran_writeAeroCamberTwist(aimInfo,
                                                  fp,
                                                  nastranInstance->feaProblem.numAero,
                                                  nastranInstance->feaProblem.feaAero,
                                                  &nastranInstance->feaProblem.feaFileFormat);
            AIM_STATUS(aimInfo, status);
        }
    }

    printf("\tDone writing aeroelastic cards\n");

    // Include mesh file
    fprintf(fp,"\nINCLUDE \'%s.bdf\'\n\n", nastranInstance->projectName);

    // End bulk data
    fprintf(fp,"ENDDATA\n");

    fclose(fp);
    fp = NULL;

    printf("\tDone writing files\n");

/*
////////////////////////////////////////
    printf("\n\n\nTESTING F06 READER\n\n");

    // FO6 data variables
    int numGridPoint = 0;
    int numEigenVector = 0;
    double **dataMatrix = NULL;
    filename = (char *) EG_alloc((strlen(nastranInstance->projectName) +
                                      strlen(".f06") + 2)*sizeof(char));

    snprintf(filename,strlen(nastranInstance->projectName) +
                                      strlen(".f06") + 2,"%s%s", nastranInstance->projectName, ".f06");

    // Open file
    fp = aim_fopen(aimInfo, filename, "r");
    AIM_FREE(filename);

    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);

        return CAPS_IOERR;
    }

    status = nastran_readF06Displacement(fp,
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

cleanup:
    if (status != CAPS_SUCCESS)
        printf("\tPremature exit in nastranAIM preAnalysis, status = %d\n",
               status);

    if (feaLoad != NULL) {
        for (i = 0; i < nastranInstance->feaProblem.numLoad; i++) {
            destroy_feaLoadStruct(&feaLoad[i]);
        }
        AIM_FREE(feaLoad);
    }

    AIM_FREE(filename);
    AIM_FREE(tempIntegerArray);

    if (fp != NULL) fclose(fp);

    return status;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS; // Function return status

  int i; //Indexing
  aimStorage *nastranInstance;

  char filename[PATH_MAX]; // File to open
  char extF06[] = ".f06";
  FILE *fp = NULL; // File pointer

  double mass = 0.0;
  double CG[3];
  double inertia[6];
  double inertiaAboutCG[3];
  double transformMatrix[9];

  nastranInstance = (aimStorage *) instStore;

  // initialize matrices to zero
  for (i = 0; i < 3; i++) {
    CG[i] = 0.0;
    inertiaAboutCG[i] = 0.0;
  }

  for (i = 0; i < 6; i++)
    inertia[i] = 0.0;

  for (i = 0; i < 9; i++)
    transformMatrix[i] = 0.0;

  snprintf(filename, PATH_MAX, "%s%s", nastranInstance->projectName, extF06);

  fp = aim_fopen(aimInfo, filename, "r");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Unable to open file: %s", filename);
    return CAPS_IOERR;
  }

  status = nastran_readF06GridPointWeightGeneratorOutput(
      fp, &mass, CG, inertia, inertiaAboutCG, transformMatrix);
  if (status == CAPS_NOTFOUND) {
    AIM_ERROR(aimInfo, "Unable to find Grid Point Weight Generator Output "
                       "section in f06.");
    goto cleanup;
  }
  else if (status != CAPS_SUCCESS) {
    AIM_ERROR(aimInfo, "Failed to extract mass properties from F06 file");
    goto cleanup;
  }

  nastranInstance->mass.vals.real   = mass;

  nastranInstance->CG.vals.reals[0] = CG[0];
  nastranInstance->CG.vals.reals[1] = CG[1];
  nastranInstance->CG.vals.reals[1] = CG[1];

  nastranInstance->Ixx.vals.real = inertia[I11];
  nastranInstance->Iyy.vals.real = inertia[I22];
  nastranInstance->Izz.vals.real = inertia[I33];
  nastranInstance->Ixy.vals.real = inertia[I21];
  nastranInstance->Ixz.vals.real = inertia[I31];
  nastranInstance->Iyz.vals.real = inertia[I32];

cleanup:
  if (fp != NULL) fclose(fp);

  return status;
}


// Set Nastran output variables
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               int index, char **aoname, capsValue *form)
{
    int status = CAPS_SUCCESS;

    /*! \page aimOutputsNastran AIM Outputs
     * The following list outlines the Nastran outputs available through the AIM interface.
     */

    #ifdef DEBUG
        printf(" nastranAIM/aimOutputs index = %d!\n", index);
    #endif

    /*! \page aimOutputsNastran AIM Outputs
     * - <B>EigenValue</B> = List of Eigen-Values (\f$ \lambda\f$) after a modal solve.
     * - <B>EigenRadian</B> = List of Eigen-Values in terms of radians (\f$ \omega = \sqrt{\lambda}\f$ ) after a modal solve.
     * - <B>EigenFrequency</B> = List of Eigen-Values in terms of frequencies (\f$ f = \frac{\omega}{2\pi}\f$) after a modal solve.
     * - <B>EigenGeneralMass</B> = List of generalized masses for the Eigen-Values.
     * - <B>EigenGeneralStiffness</B> = List of generalized stiffness for the Eigen-Values.
     * - <B>Objective</B> = Final objective value for a design optimization case.
     * - <B>ObjectiveHistory</B> = List of objective value for the history of a design optimization case.
     */

    if (index == outEigenValue) {
        *aoname = EG_strdup("EigenValue");

    } else if (index == outEigenRadian) {
        *aoname = EG_strdup("EigenRadian");

    } else if (index == outEigenFrequency) {
        *aoname = EG_strdup("EigenFrequency");

    } else if (index == outEigenGeneralMass) {
        *aoname = EG_strdup("EigenGeneralMass");

    } else if (index == outEigenGeneralStiffness) {
        *aoname = EG_strdup("EigenGeneralStiffness");

    } else if (index == outObjective) {
        *aoname = EG_strdup("Objective");

    } else if (index == outObjectiveHistory) {
        *aoname = EG_strdup("ObjectiveHistory");


    /*! \page aimOutputsNastran AIM Outputs
     * - <B>Mass</B> = Total mass of the model.
     * - <B>CG</B> = Center of gravity of the model.
     * - <B>Ixx</B> = Moment of inertia
     * - <B>Iyy</B> = Moment of inertia
     * - <B>Izz</B> = Moment of inertia
     * - <B>Ixy</B> = Moment of inertia
     * - <B>Izy</B> = Moment of inertia
     * - <B>Iyz</B> = Moment of inertia
     * - <B>I_Vector</B> = Moment of inertia vector
     *  \f[
     *    \vec{I} = \begin{bmatrix} I_{xx} & I_{yy} & I_{zz} & I_{xy} & I_{xz} & I_{yz} \end{bmatrix}
     *  \f]
     * - <B>I_Lower</B> = Moment of inertia lower triangular tensor
     * \f[
     *  \vec{I}_{lower} = \begin{bmatrix} I_{xx} & -I_{xy} & I_{yy} & -I_{xz} & -I_{yz} & I_{zz} \end{bmatrix},
     * \f]
     * - <B>I_Upper</B> = Moment of inertia upper triangular tensor
     * \f[
     *  \vec{I}_{upper} = \begin{bmatrix} I_{xx} & -I_{xy} & -I_{xz} & I_{yy} & -I_{yz} & I_{zz} \end{bmatrix},
     * \f]
     * - <B>I_Tensor</B> = Moment of inertia tensor
     * \f[
     * \bar{\bar{I}} =
     *  \begin{bmatrix}
     *    I_{xx} & -I_{xy} & -I_{xz} \\
     *   -I_{xy} &  I_{yy} & -I_{yz} \\
     *   -I_{xz} & -I_{yz} &  I_{zz}
     *  \end{bmatrix}
     * \f]
     * - <B>MassProp</B> = JSON String Mass Properties
     * - <B>MassPropLink</B> = Mass Properties
     * Mass properties that can be linked to analysis input MassPropLink
     */

    } else if (index == outMass) {
        *aoname = EG_strdup("Mass");
        form->dim  = Scalar;

    } else if (index == outCG) {
        *aoname = EG_strdup("CG");
        form->nrow = 3;
        form->dim  = Vector;

    } else if (index == outIxx) {
        *aoname = EG_strdup("Ixx");
        form->dim  = Scalar;

    } else if (index == outIyy) {
        *aoname = EG_strdup("Iyy");
        form->dim  = Scalar;

    } else if (index == outIzz) {
        *aoname = EG_strdup("Izz");
        form->dim  = Scalar;

    } else if (index == outIxy) {
        *aoname = EG_strdup("Ixy");
        form->dim  = Scalar;

    } else if (index == outIxz) {
        *aoname = EG_strdup("Ixz");
        form->dim  = Scalar;

    } else if (index == outIyz) {
        *aoname = EG_strdup("Iyz");
        form->dim  = Scalar;

    } else if (index == outI_Vector) {
        *aoname = EG_strdup("I_Vector");
        form->nrow = 6;
        form->dim  = Vector;

    } else if (index == outI_Lower) {
        *aoname = EG_strdup("I_Lower");
        form->nrow = 6;
        form->dim  = Vector;

    } else if (index == outI_Upper) {
        *aoname = EG_strdup("I_Upper");
        form->nrow = 6;
        form->dim  = Vector;

    } else if (index == outI_Tensor) {
        *aoname = EG_strdup("I_Tensor");
        form->nrow = 9;
        form->dim  = Array2D;

    } else if (index == outMassProp) {
        *aoname = EG_strdup("MassProp");
        form->type = String;
        form->nullVal = IsNull;
        goto cleanup;

    } else if (index == outMassPropLink) {
        *aoname = EG_strdup("MassPropLink");
        form->type = Pointer;
        form->nullVal = IsNull;
        AIM_STRDUP(form->units, "feaMassPropStruct", aimInfo, status);
        goto cleanup;

    } else if (index == outF06FileLink) {
        *aoname = EG_strdup("F06FileLink");
        form->type = Pointer;
        form->nullVal = IsNull;
        AIM_STRDUP(form->units, "feaSolFileStruct", aimInfo, status);
        goto cleanup;
    }

    if (*aoname == NULL) return EGADS_MALLOC;

    form->type       = Double;
    form->lfixed     = Fixed;
    form->sfixed     = Fixed;
    form->vals.reals = NULL;
    form->vals.real  = 0;

    if (form->nrow > 1) {
      AIM_ALLOC(form->vals.reals, form->nrow, double, aimInfo, status);
      memset(form->vals.reals, 0, form->nrow*sizeof(double));
    }

//    if (units->momentOfInertia != NULL &&
//        index >= outI_Vector && index <= outI_Tensor) {
//      AIM_STRDUP(form->units, units->momentOfInertia, aimInfo, status);
//    }

    if (index == outI_Tensor) {
        form->nrow     = 3;
        form->ncol     = 3;
    }

cleanup:
    return status;
}


// Calculate Nastran output
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int status = CAPS_SUCCESS; // Function return status

    int i; //Indexing

    double mass = 0, CGx = 0, CGy = 0, CGz = 0;
    double Ixx = 0;
    double Ixy = 0;
    double Izz = 0;
    double Ixz = 0;
    double Iyy = 0;
    double Iyz = 0;
    double *I;
    char massProp[512];

    int numData = 0;
    double *dataVector = NULL;
    double **dataMatrix = NULL;
    aimStorage *nastranInstance;

    char filename[PATH_MAX]; // File to open
    char aimFile[PATH_MAX];
    char extF06[] = ".f06", extOP2[] = ".op2";
    FILE *fp = NULL; // File pointer
    feaUnitsStruct *units=NULL;

    nastranInstance = (aimStorage *) instStore;
    units = &nastranInstance->units;

    if (index >= outEigenValue &&
        index <= outEigenGeneralStiffness) {

        snprintf(filename, PATH_MAX, "%s%s", nastranInstance->projectName, extF06);

        fp = aim_fopen(aimInfo, filename, "r");
        if (fp == NULL) {
#ifdef DEBUG
            printf(" nastranAIM/aimCalcOutput Cannot open Output file!\n");
#endif
            AIM_ERROR(aimInfo, "Cannot open %s", filename);
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

                AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
                if (val->vals.reals != NULL) {

                    for (i = 0; i < val->length; i++) {

                        val->vals.reals[i] = dataMatrix[i][index-1];
                    }

                } else status = EGADS_MALLOC;
            }
        }
    } else if (index == outObjective || index == outObjectiveHistory) {
        snprintf(filename, PATH_MAX, "%s%s", nastranInstance->projectName, extOP2);

        status = nastran_readOP2Objective(filename, &numData, &dataVector);
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

                AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
                if (val->vals.reals != NULL) {

                    for (i = 0; i < val->length; i++) {

                        val->vals.reals[i] = dataVector[i];
                    }

                } else status = EGADS_MALLOC;
            }
        }

    } else if (index == outMass) {
      status = aim_copyValue(&nastranInstance->mass, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outCG) {
      status = aim_copyValue(&nastranInstance->CG, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outIxx) {
      status = aim_copyValue(&nastranInstance->Ixx, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outIyy) {
      status = aim_copyValue(&nastranInstance->Iyy, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outIzz) {
      status = aim_copyValue(&nastranInstance->Izz, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outIxy) {
      status = aim_copyValue(&nastranInstance->Ixy, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outIxz) {
      status = aim_copyValue(&nastranInstance->Ixz, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outIyz) {
      status = aim_copyValue(&nastranInstance->Iyz, val);
      AIM_STATUS(aimInfo, status);

    } else if (index == outI_Vector) {
      val->type = Double;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

      val->vals.reals[0] = nastranInstance->Ixx.vals.real;
      val->vals.reals[1] = nastranInstance->Iyy.vals.real;
      val->vals.reals[2] = nastranInstance->Izz.vals.real;
      val->vals.reals[3] = nastranInstance->Ixy.vals.real;
      val->vals.reals[4] = nastranInstance->Ixz.vals.real;
      val->vals.reals[5] = nastranInstance->Iyz.vals.real;

    } else if (index == outI_Lower) {
      val->type = Double;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

      val->vals.reals[0] =  nastranInstance->Ixx.vals.real;
      val->vals.reals[1] = -nastranInstance->Ixy.vals.real;
      val->vals.reals[2] =  nastranInstance->Iyy.vals.real;
      val->vals.reals[3] = -nastranInstance->Ixz.vals.real;
      val->vals.reals[4] = -nastranInstance->Iyz.vals.real;
      val->vals.reals[5] =  nastranInstance->Izz.vals.real;

    } else if (index == outI_Upper) {
      val->type = Double;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

      val->vals.reals[0] =  nastranInstance->Ixx.vals.real;
      val->vals.reals[1] = -nastranInstance->Ixy.vals.real;
      val->vals.reals[2] = -nastranInstance->Ixz.vals.real;
      val->vals.reals[3] =  nastranInstance->Iyy.vals.real;
      val->vals.reals[4] = -nastranInstance->Iyz.vals.real;
      val->vals.reals[5] =  nastranInstance->Izz.vals.real;

    } else if (index == outI_Tensor) {
      val->type = Double;
      AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

      Ixx = nastranInstance->Ixx.vals.real;
      Iyy = nastranInstance->Iyy.vals.real;
      Iyz = nastranInstance->Iyz.vals.real;
      Ixy = nastranInstance->Ixy.vals.real;
      Ixz = nastranInstance->Ixz.vals.real;
      Izz = nastranInstance->Izz.vals.real;
      I = val->vals.reals;

      // populate the tensor
      I[0] =  Ixx; I[1] = -Ixy; I[2] = -Ixz;
      I[3] = -Ixy; I[4] =  Iyy; I[5] = -Iyz;
      I[6] = -Ixz; I[7] = -Iyz; I[8] =  Izz;

    } else if (index == outMassProp) {

      mass = nastranInstance->mass.vals.real;

      CGx = nastranInstance->CG.vals.reals[0];
      CGy = nastranInstance->CG.vals.reals[1];
      CGz = nastranInstance->CG.vals.reals[2];

      Ixx = nastranInstance->Ixx.vals.real;
      Iyy = nastranInstance->Iyy.vals.real;
      Iyz = nastranInstance->Iyz.vals.real;
      Ixy = nastranInstance->Ixy.vals.real;
      Ixz = nastranInstance->Ixz.vals.real;
      Izz = nastranInstance->Izz.vals.real;

      if (units->mass != NULL) {

        snprintf(massProp, 512, "{\"mass\":[%20.14le, %s], \"CG\":[[%20.14le,%20.14le,%20.14le], %s], \"massInertia\":[[%20.14le, %20.14le, %20.14le, %20.14le, %20.14le, %20.14le], %s]}",
            mass, units->mass, CGx, CGy, CGz, units->length, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, units->momentOfInertia);

      } else {

        snprintf(massProp, 512, "{\"mass\":%20.14le, \"CG\":[%20.14le,%20.14le,%20.14le], \"massInertia\":[%20.14le, %20.14le, %20.14le, %20.14le, %20.14le, %20.14le]}",
            mass, CGx, CGy, CGz, Ixx, Iyy, Izz, Ixy, Ixz, Iyz);

      }

      AIM_STRDUP(val->vals.string, massProp, aimInfo, status);
    } else if (index == outMassPropLink) {

      nastranInstance->feaMassProp.mass = nastranInstance->mass.vals.real;
      for (i = 0; i < 3; i++)
        nastranInstance->feaMassProp.CG[i] = nastranInstance->CG.vals.reals[i];

      nastranInstance->feaMassProp.massInertia[I11] = nastranInstance->Ixx.vals.real;
      nastranInstance->feaMassProp.massInertia[I22] = nastranInstance->Iyy.vals.real;
      nastranInstance->feaMassProp.massInertia[I32] = nastranInstance->Iyz.vals.real;
      nastranInstance->feaMassProp.massInertia[I21] = nastranInstance->Ixy.vals.real;
      nastranInstance->feaMassProp.massInertia[I31] = nastranInstance->Ixz.vals.real;
      nastranInstance->feaMassProp.massInertia[I33] = nastranInstance->Izz.vals.real;

      nastranInstance->feaMassProp.massUnit = units->mass;
      nastranInstance->feaMassProp.lengthUnit = units->length;
      nastranInstance->feaMassProp.momentOfInertiaUnit = units->momentOfInertia;

      /*@-immediatetrans@*/
      val->vals.AIMptr = (void*)&nastranInstance->feaMassProp;
      /*@+immediatetrans@*/

    } else if (index == outF06FileLink) {

      // Get the local filename in the analysis directory
      snprintf(filename, PATH_MAX, "%s%s", nastranInstance->projectName, extF06);

      // Get the absolute path to the file
      status = aim_file(aimInfo, filename, aimFile);
      AIM_STATUS(aimInfo, status);

      AIM_FREE(nastranInstance->feaF06File.filename);

      nastranInstance->feaF06File.fileForm = MSC_NASTRAN;
      AIM_STRDUP(nastranInstance->feaF06File.filename, aimFile, aimInfo, status);

      /*@-immediatetrans@*/
      val->vals.AIMptr = (void*)&nastranInstance->feaF06File;
      /*@+immediatetrans@*/
    }

cleanup:
    // Restore the path we came in with
    if (fp != NULL) fclose(fp);

    if (dataMatrix != NULL) {
        for (i = 0; i < numData; i++) {
            if (dataMatrix[i] != NULL) AIM_FREE(dataMatrix[i]);
        }
        AIM_FREE(dataMatrix);
    }

    if (dataVector != NULL) AIM_FREE(dataVector);

    return status;
}


void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *nastranInstance;

#ifdef DEBUG
    printf(" nastranAIM/aimCleanup!\n");
#endif
    nastranInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(nastranInstance);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during clean up of instance\n", status);

    AIM_FREE(nastranInstance);
}


int aimDiscr(char *tname, capsDiscr *discr)
{
    int        status = CAPS_SUCCESS; // Function return status
    int        i;
    ego        *tess = NULL;
    aimStorage *nastranInstance;

#ifdef DEBUG
    printf(" nastranAIM/aimDiscr: tname = %s!\n", tname);
#endif
    if (tname == NULL) return CAPS_NOTFOUND;

    nastranInstance = (aimStorage *) discr->instStore;

    // Check and generate/retrieve the mesh
    status = checkAndCreateMesh(discr->aInfo, nastranInstance);
    AIM_STATUS(discr->aInfo, status);

    AIM_ALLOC(tess, nastranInstance->numMesh, ego, discr->aInfo, status);
    for (i = 0; i < nastranInstance->numMesh; i++) {
      tess[i] = nastranInstance->feaMesh[i].egadsTess;
    }

    status = mesh_fillDiscr(tname, &nastranInstance->attrMap, nastranInstance->numMesh, tess, discr);
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

    /*! \page dataTransferNastran AIM Data Transfer
     *
     * The Nastran AIM has the ability to transfer displacements and eigenvectors from the AIM and pressure
     * distributions to the AIM using the conservative and interpolative data transfer schemes in CAPS.
     *
     * \section dataFromNastran Data transfer from Nastran (FieldOut)
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
     * \section dataToNastran Data transfer to Nastran (FieldIn)
     * <ul>
     *  <li> <B>"Pressure"</B> </li> <br>
     *  Writes appropriate load cards using the provided pressure distribution.
     * </ul>
     *
     */

    int status; // Function return status
    int i, j, dataPoint, bIndex; // Indexing
    aimStorage *nastranInstance;

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
    size_t stringLength;
    char *filename = NULL;
    FILE *fp=NULL; // File pointer

#ifdef DEBUG
    printf(" nastranAIM/aimTransfer name = %s  npts = %d/%d!\n",
           dataName, numPoint, dataRank);
#endif
    nastranInstance = (aimStorage *) discr->instStore;

    //Get the appropriate parts of the tessellation to data
    //storage = (int *) discr->ptrm;
    //capsGroupList = &storage[0]; // List of boundary ID (attrMap) in the transfer

    if (strcasecmp(dataName, "Displacement") != 0 &&
        strncmp(dataName, "EigenVector", 11) != 0) {

        printf("Unrecognized data transfer variable - %s\n", dataName);
        return CAPS_NOTFOUND;
    }

    stringLength = strlen(nastranInstance->projectName) + strlen(extF06) + 1;
    AIM_ALLOC(filename, stringLength, char, discr->aInfo, status);

    snprintf(filename,stringLength,"%s%s", nastranInstance->projectName, extF06);

    // Open file
    fp = aim_fopen(discr->aInfo, filename, "r");
    if (fp == NULL) {
        printf("Unable to open file: %s\n", filename);
        AIM_FREE(filename);
        return CAPS_IOERR;
    }

    AIM_FREE(filename);

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
            AIM_FREE(dataMatrix[i]);
        }
        AIM_FREE(dataMatrix);
    }

    return status;
}


void aimFreeDiscrPtr(void *ptr)
{
    // Extra information to store into the discr void pointer - just a int array
    AIM_FREE(ptr);
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
