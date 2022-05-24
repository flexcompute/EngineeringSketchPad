/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AFLR3 AIM
 *
 *      Written by Dr. Ryan Durscher AFRL/RQVC
 *
 *      This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 *
 *
 */

/*!\mainpage Introduction
 *
 * \section overviewAFLR3 AFLR3 AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact with the
 * unstructured, volumetric grid generator AFLR3 \cite Marcum1995 \cite Marcum1998.
 *
 * The AFLR3 AIM provides the CAPS users with the ability to generate "unstructured tetrahedral element grids" using an
 * "Advancing-Front/Local-Reconnection (AFLR) procedure." Additionally, an "Advancing-Normal Boundary-Layer (ANBL) procedure" may be
 * used "to generate a tetrahedral/pentahedral/hexahedral BL grid adjacent to" specified surfaces.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsAFLR3 and \ref aimOutputsAFLR3, respectively.
 * The complete AFLR documentation is available at the <a href="https://www.simcenter.msstate.edu/software/documentation/system/index.html">SimCenter</a>.
 *
 * Example volumes meshes:
 *  \image html multiBodyAFRL3.png "AFLR3 meshing example - Multiple Airfoils with Boundary Layer" width=500px
 *  \image latex multiBodyAFRL3.png "AFLR3 meshing example - Multiple Airfoils with Boundary Layer" width=5cm
 *
 * \section clearanceAFLR3 Clearance Statement
 *  This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */

#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"
#include "aimMesh.h"

#include "meshUtils.h"       // Collection of helper functions for meshing
#include "miscUtils.h"
#include "deprecateUtils.h"

#include "aflr3_Interface.h" // Bring in AFLR3 'interface' functions

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

#define MXCHAR 255

//#define DEBUG


enum aimOutputs
{
  NumberOfElement = 1,         /* index is 1-based */
  NumberOfNode,
  Volume_Mesh,
  NUMOUT = Volume_Mesh         /* Total number of outputs */
};

#define NODATATRANSFER "noDataTransfer.%d"

////// Global variables //////
typedef struct {

    int numSurfaceMesh;
    meshStruct *surfaceMesh;

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    mapAttrToIndexStruct meshMap;

    // Container for mesh input
    meshInputStruct meshInput;

    // Mesh references for link
    int numMeshRef;
    aimMeshRef *meshRef;

} aimStorage;


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


static int get_MeshIndex(mapAttrToIndexStruct *meshMap, ego ebody,
                         meshElementStruct *element, int *meshIndex)
{

    int status;

    int numFace;

    const char *meshName;
    ego *efaces=NULL;

    *meshIndex = -1;

    if (element->elementType != Triangle &&
        element->elementType != Quadrilateral) {
        printf("Error: Unsupported element type\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    status = EG_getBodyTopos(ebody, NULL, FACE, &numFace, &efaces);
    if (status < EGADS_SUCCESS) goto cleanup;
    if (efaces == NULL) {
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    // get the face index
    if (element->topoIndex < 1 || element->topoIndex > numFace) {
        printf("Error: Element '%d': Invalid face topological index: %d, [1-%d]\n",
               element->elementID, element->topoIndex, numFace);
        goto cleanup;
    }

    status = retrieve_CAPSMeshAttr(efaces[element->topoIndex-1], &meshName);
    if (status == EGADS_SUCCESS) {

        status = get_mapAttrToIndexIndex(meshMap, meshName, meshIndex);

        if (status != CAPS_SUCCESS)
          printf("\tNo capsMesh \"%s\" not found in attribute map\n", meshName);
    } else {
        // OK not to find capsMesh attribute
    }

    status = CAPS_SUCCESS;
    goto cleanup;

cleanup:
    if (status != CAPS_SUCCESS)
      printf("\tPremature exit in get_MeshIndex, status = %d\n", status);
    EG_free(efaces);
    return status;
}


static int destroy_aimStorage(aimStorage *aflr3Instance)
{

    int i; // Indexing

    int status; // Function return status

#if _MSC_VER >= 1900
    _iob[0] = *stdin;
    _iob[1] = *stdout;
    _iob[2] = *stderr;
#endif

    // Destroy meshInput
    status = destroy_meshInputStruct(&aflr3Instance->meshInput);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, aflr3AIM meshInput cleanup!!!\n", status);

    // Surface meshes are always referenced
    aflr3Instance->numSurfaceMesh = 0;
    aflr3Instance->surfaceMesh = NULL;

    // Destroy attribute to index map
    status = destroy_mapAttrToIndexStruct(&aflr3Instance->groupMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, aflr3AIM destroy_mapAttrToIndexStruct cleanup!!!\n",
             status);

    status = destroy_mapAttrToIndexStruct(&aflr3Instance->meshMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, aflr3AIM destroy_mapAttrToIndexStruct cleanup!!!\n",
             status);

    // Free the meshRef
    for (i = 0; i < aflr3Instance->numMeshRef; i++)
      aim_freeMeshRef(&aflr3Instance->meshRef[i]);
    AIM_FREE(aflr3Instance->meshRef);

    return status;
}


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS; // Function status return
    aimStorage *aflr3Instance=NULL;

#ifdef WIN32
    #if  _MSC_VER >= 1900
        _iob[0] = *stdin;
        _iob[1] = *stdout;
        _iob[2] = *stderr;
    #endif
#endif

#ifdef DEBUG
    printf("\n aflr3AIM/aimInitialize   instance = %d!\n", inst);
#endif

    /* specify the number of analysis input and out "parameters" */
    *nIn     = NUMINPUT;
    *nOut    = NUMOUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames  = NULL;
    *franks   = NULL;
    *fInOut   = NULL;

    // Allocate aflr3Instance
    AIM_ALLOC(aflr3Instance, 1, aimStorage, aimInfo, status);
    *instStore = aflr3Instance;

    // Set initial values for aflr3Instance

    // Container for surface meshes
    aflr3Instance->numSurfaceMesh = 0;
    aflr3Instance->surfaceMesh = NULL;

    // Mesh reference passed to solver
    aflr3Instance->numMeshRef = 0;
    aflr3Instance->meshRef = NULL;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&aflr3Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&aflr3Instance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Container for mesh input
    status = initiate_meshInputStruct(&aflr3Instance->meshInput);
    AIM_STATUS(aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);
    return status;
}


int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
    /*! \page aimInputsAFLR3 AIM Inputs
     * The following list outlines the AFLR3 meshing options along with their default value available
     * through the AIM interface.
     */

    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" aflr3AIM/aimInputs index = %d!\n", index);
#endif

    // Meshing Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name"); // If NULL a volume grid won't be written by the AIM
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->vals.string  = NULL;
        //defval->vals.string  = EG_strdup("CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsAFLR3
         * - <B> Proj_Name = NULL</B> <br>
         * This corresponds to the output name of the mesh. If left NULL, the mesh is not written to a file.
         */

    } else if (index == Mesh_Quiet_Flag) {
        *ainame               = EG_strdup("Mesh_Quiet_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsAFLR3
         * - <B> Mesh_Quiet_Flag = False</B> <br>
         * Suppression of mesh generator (not including errors)
         */

    } else if (index == Mesh_Format) {
        *ainame               = EG_strdup("Mesh_Format");
        defval->type          = String;
        defval->vals.string   = EG_strdup("AFLR3"); // AFLR3, SU2, VTK

        /*! \page aimInputsAFLR3
         * - <B> Mesh_Format = "AFLR3"</B> <br>
         * Mesh output format. Available format names  include: "AFLR3", "SU2", "Nastran", "Tecplot", and "VTK".
         */

    } else if (index == Mesh_ASCII_Flag) {
        *ainame               = EG_strdup("Mesh_ASCII_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = true;

        /*! \page aimInputsAFLR3
         * - <B> Mesh_ASCII_Flag = True</B> <br>
         * Output mesh in ASCII format, otherwise write a binary file, if applicable.
         */

    } else if (index == Mesh_Gen_Input_String) {
        *ainame               = EG_strdup("Mesh_Gen_Input_String");
        defval->type          = String;
        defval->nullVal       = IsNull;
        defval->vals.string   = NULL;

        /*! \page aimInputsAFLR3
         * - <B> Mesh_Gen_Input_String = NULL</B> <br>
         * Meshing program command line string (as if called in bash mode). Use this to specify more
         * complicated options/use features of the mesher not currently exposed through other AIM input
         * variables. Note that this is the exact string that will be provided to the volume mesher; no
         * modifications will be made. If left NULL an input string will be created based on default values
         * of the relevant AIM input variables.
         */

    } else if (index == Multiple_Mesh) {
        *ainame               = EG_strdup("Multiple_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = (int) false;

        /*! \page aimInputsAFLR3
         * - <B> Multiple_Mesh = False</B> <br>
         * If set to True a volume will be generated for each body. When set to False (default value) only a single volume
         * mesh will be created.
         */

    } else if (index == Mesh_Sizing) {
        *ainame              = EG_strdup("Mesh_Sizing");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAFLR3
         * - <B>Mesh_Sizing = NULL </B> <br>
         * See \ref meshSizingProp for additional details.
         */
    } else if (index == BL_Initial_Spacing) {
        *ainame               = EG_strdup("BL_Initial_Spacing");
        defval->type          = Double;
        defval->nullVal       = NotNull;
        defval->vals.real     = 0.0;

        /*! \page aimInputsAFLR3
         * - <B> BL_Initial_Spacing = 0.0</B> <br>
         * Initial mesh spacing when growing a boundary layer that is applied to all bodies. <br><br>
         * Note: Both "BL_Initial_Spacing" and "BL_Thickness" must be non-zero for values to be applied. If "Multiple_Mesh" is
         * False (default value) these values will not be applied to the largest body (if more than 1 body exist in the AIM),
         * as that body is assumed to be a bounding box (e.g. a farfield boundary in a CFD simulation). Boundary spacing
         * and thickness specified through the use of the "Mesh_Sizing" input (see \ref meshSizingProp for additional details)
         * will take precedence over the values specified for "BL_Initial_Spacing" and "BL_Thickness".
         */

    } else if (index == BL_Thickness) {
        *ainame               = EG_strdup("BL_Thickness");
        defval->type          = Double;
        defval->nullVal       = NotNull;
        defval->vals.real     = 0.0;

        /*! \page aimInputsAFLR3
         * - <B> BL_Thickness = 0.0</B> <br>
         * Total boundary layer thickness that is applied to all bodies. <br>
         * This is a lower bound on the desired thickness. The height can be limited with "nbl".<br><br>
         * Note: see "BL_Initial_Spacing" and "BL_Max_Layers" for additional details
         */

    } else if (index == BL_Max_Layers) {
        *ainame               = EG_strdup("BL_Max_Layers");
        defval->type          = Integer;
        defval->nullVal       = IsNull;
        defval->vals.integer  = 10000;

        /*! \page aimInputsAFLR3
         * - <B> BL_Max_Layers = 10000</B> <br>
         * Maximum BL grid layers to generate.
         */

    } else if (index == BL_Max_Layer_Diff) {
        *ainame               = EG_strdup("BL_Max_Layer_Diff");
        defval->type          = Integer;
        defval->nullVal       = IsNull;
        defval->vals.integer  = 0;

        /*! \page aimInputsAFLR3
         * - <B> BL_Max_Layer_Diff = 0</B> <br>
         * Maximum difference in BL levels.<br>
         * If BL_Max_Layer_Diff > 0 then the maximum difference between the number of BL levels for<br>
         * the BL nodes on a given BL boundary surface face is limited to BL_Max_Layer_Diff. <br>
         * Any active BL node that would allow the number of levels to be greater is <br>
         * terminated. <br>
         * If BL_Max_Layer_Diff = 0 then the difference in BL levels is ignored.
         */

    } else if (index == Surface_Mesh) {
        *ainame             = AIM_NAME(Surface_Mesh);
        defval->type        = Pointer;
        defval->dim         = Vector;
        defval->lfixed      = Change;
        defval->sfixed      = Change;
        defval->vals.AIMptr = NULL;
        defval->nullVal     = IsNull;
        AIM_STRDUP(defval->units, "meshStruct", aimInfo, status);

        /*! \page aimInputsAFLR3
         * - <B>Surface_Mesh = NULL</B> <br>
         * A Surface_Mesh link.
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
    int i, j, propIndex, ibody, pointIndex[4] = {-1,-1,-1,-1}, meshIndex; // Indexing

    int status; // Return status

    const char *intents;

    // Incoming bodies
    ego *bodies = NULL, body;
    int numBody = 0;

    int tempInt1, tempInt2, numElementCheck;

    // Container for volume mesh
    int numVolumeMesh=0;
    meshStruct *volumeMesh=NULL;

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;

    // Boundary Layer meshing related variables
    double globalBLThickness = 0.0, globalBLSpacing = 0.0, capsMeshLength = 0;
    int createBL = (int) false; // Indicate if boundary layer meshes should be generated
    int *blFlag = NULL ;        // Boundary layer flag [numElement] - Total elements for all bodies
    double *blThickness = NULL; // Boundary layer thickness [numNode] - Total nodes for all bodies
    double *blSpacing = NULL;   // Boundary layer spacing [numNode] - Total nodes for all bodies

    int elementOffSet = 0, nodeOffSet = 0;
    aimStorage *aflr3Instance;

    // Bounding box variables
    int boundingBoxIndex = CAPSMAGIC;
    double  box[6], boxMax[6] = {0,0,0,0,0,0};//, refLen = -1.0;

    // Output grid
    char *filename = NULL;
    char bodyNumberFile[42];
    char aimFile[PATH_MAX];
    FILE *fp;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" aflr3AIM/aimPreAnalysis numBody = %d!\n", numBody);
#endif

    if (numBody <= 0 || bodies == NULL) {
#ifdef DEBUG
        printf(" aflr3AIM/aimPreAnalysis No Bodies!\n");
#endif
        return CAPS_SOURCEERR;
    }

    aflr3Instance = (aimStorage *) instStore;
    if (aimInputs == NULL) return CAPS_NULLVALUE;

    // remove previous meshes
    for (ibody = 0; ibody < aflr3Instance->numMeshRef; ibody++) {
      status = aim_deleteMeshes(aimInfo, &aflr3Instance->meshRef[ibody]);
      AIM_STATUS(aimInfo, status);
    }

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(aflr3Instance);
    AIM_STATUS(aimInfo, status, "aflr3AIM aimStorage cleanup!!!");

    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            3, // Only search down to the node level of the EGADS bodyIndex
                                            &aflr3Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    // Get capsMesh name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSMeshAttrToIndexMap(numBody,
                                           bodies,
                                           3, // Only search down to the node level of the EGADS bodyIndex
                                           &aflr3Instance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Setup meshing input structure

    // Get Tessellation parameters -Tess_Params
    aflr3Instance->meshInput.paramTess[0] = 0;
    aflr3Instance->meshInput.paramTess[1] = 0;
    aflr3Instance->meshInput.paramTess[2] = 0;

    aflr3Instance->meshInput.quiet = aimInputs[Mesh_Quiet_Flag-1].vals.integer;
    aflr3Instance->meshInput.outputASCIIFlag = aimInputs[Mesh_ASCII_Flag-1].vals.integer;

    // Mesh Format
    AIM_STRDUP(aflr3Instance->meshInput.outputFormat, aimInputs[Mesh_Format-1].vals.string, aimInfo, status);

    // Project Name
    if (aimInputs[Proj_Name-1].nullVal != IsNull)
        AIM_STRDUP(aflr3Instance->meshInput.outputFileName, aimInputs[Proj_Name-1].vals.string, aimInfo, status);

    // Set aflr3 specific mesh aimInputs
    if (aimInputs[Mesh_Gen_Input_String-1].nullVal != IsNull)
       AIM_STRDUP(aflr3Instance->meshInput.aflr3Input.meshInputString, aimInputs[Mesh_Gen_Input_String-1].vals.string, aimInfo, status);

    status = populate_bndCondStruct_from_mapAttrToIndexStruct(&aflr3Instance->groupMap,
                                                              &aflr3Instance->meshInput.bndConds);
    AIM_STATUS(aimInfo, status);

    // Get mesh sizing parameters
    if (aimInputs[Mesh_Sizing-1].nullVal != IsNull) {

        status = deprecate_SizingAttr(aimInfo,
                                      aimInputs[Mesh_Sizing-1].length,
                                      aimInputs[Mesh_Sizing-1].vals.tuple,
                                      &aflr3Instance->meshMap,
                                      &aflr3Instance->groupMap);
        AIM_STATUS(aimInfo, status);

        status = mesh_getSizingProp(aimInfo,
                                    aimInputs[Mesh_Sizing-1].length,
                                    aimInputs[Mesh_Sizing-1].vals.tuple,
                                    &aflr3Instance->meshMap,
                                    &numMeshProp,
                                    &meshProp);
        AIM_STATUS(aimInfo, status);
    }

    // Get surface mesh
    if (aimInputs[Surface_Mesh-1].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Surface_Mesh, "'Surface_Mesh' input must be linked to an output 'Surface_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get mesh
    aflr3Instance->numSurfaceMesh = aimInputs[Surface_Mesh-1].length;
    aflr3Instance->surfaceMesh    = (meshStruct *)aimInputs[Surface_Mesh-1].vals.AIMptr;

    if (aflr3Instance->numSurfaceMesh != numBody) {
        AIM_ANALYSISIN_ERROR(aimInfo, Surface_Mesh, "Number of linked surface meshes (%d) does not match the number of bodies (%d)!",
                             aflr3Instance->numSurfaceMesh, numBody);
        return CAPS_SOURCEERR;
    }

    // If we are only going to have one volume determine the bounding box body so that boundary layer
    // parameters wont be applied to it
    if (aimInputs[Multiple_Mesh-1].vals.integer == (int) false) {

        // Determine which body is the bounding body based on size
        if (numBody > 1) {
            for (ibody = 0; ibody < numBody; ibody++) {

                // Get bounding box for the body
                status = EG_getBoundingBox(bodies[ibody], box);
                if (status != EGADS_SUCCESS) {
                    printf(" EG_getBoundingBox = %d\n\n", status);
                    return status;
                }

                // Just copy the box coordinates on the first go around
                if (ibody == 0) {

                    memcpy(boxMax, box, sizeof(box));

                    // Set body as the bounding box (ie. farfield)
                    boundingBoxIndex = ibody;

                    // Else compare with the "max" box size
                } else if (boxMax[0] >= box[0] &&
                           boxMax[1] >= box[1] &&
                           boxMax[2] >= box[2] &&
                           boxMax[3] <= box[3] &&
                           boxMax[4] <= box[4] &&
                           boxMax[5] <= box[5]) {

                    // If bigger copy coordinates
                    memcpy(boxMax, box, sizeof(box));

                    // Set body as the bounding box (ie. farfield)
                    boundingBoxIndex = ibody;
                }
            }
        }
    }

    // Get global
    globalBLThickness = aimInputs[BL_Thickness-1].vals.real;
    globalBLSpacing   = aimInputs[BL_Initial_Spacing-1].vals.real;

    // check that both
    //   globalBLThickness == 0 && globalBLSpacing == 0
    // or
    //   globalBLThickness != 0 && globalBLSpacing != 0
    if (!((globalBLThickness == 0.0 && globalBLSpacing == 0.0) ||
          (globalBLThickness != 0.0 && globalBLSpacing != 0.0))) {
        AIM_ERROR(aimInfo, "Both BL_Thickness = %le and BL_Initial_Spacing = %le",
                 globalBLThickness, globalBLSpacing);
        AIM_ADDLINE(aimInfo, "must be zero or non-zero.");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // check if boundary layer meshing is requested
    if (globalBLThickness != 0.0 &&
        globalBLSpacing   != 0.0) {

        // Turn on BL flags if globals are set
        createBL = (int) true;
    }

    // Loop through surface meshes to set boundary layer parameters
    for (ibody = 0; ibody < aflr3Instance->numSurfaceMesh && createBL == (int) false; ibody++) {

        // Loop through meshing properties and see if boundaryLayerThickness and boundaryLayerSpacing have been specified
        for (propIndex = 0; propIndex < numMeshProp && createBL == (int) false; propIndex++) {

            // If no boundary layer specified in meshProp continue
            if (meshProp != NULL)
                if (meshProp[propIndex].boundaryLayerThickness == 0 ||
                    meshProp[propIndex].boundaryLayerSpacing   == 0) continue;

//            // Loop through Elements and see if marker match
//            for (i = 0; i < aflr3Instance->surfaceMesh[bodyIndex].numElement && createBL == (int) false; i++) {
//
//                //If they don't match continue
//                if (aflr3Instance->surfaceMesh[bodyIndex].element[i].markerID != meshProp[propIndex].attrIndex) continue;

            // Set face "bl" flag
            createBL = (int) true;
//            }
        }
    }

    // Get the capsMeshLenght if boundary layer meshing has been requested
    if (createBL == (int) true) {
        status = check_CAPSMeshLength(numBody, bodies, &capsMeshLength);
        if (capsMeshLength <= 0 || status != CAPS_SUCCESS) {
            if (status != CAPS_SUCCESS) {
              AIM_ERROR(aimInfo, "capsMeshLength is not set on any body.\n");
            } else {
              AIM_ERROR(aimInfo, "capsMeshLength: %f\n", capsMeshLength);
            }
            AIM_ADDLINE(aimInfo,
                   "\nThe capsMeshLength attribute must present on at least\n"
                   "one body for boundary layer generation.\n"
                   "\n"
                   "capsMeshLength should be a a positive value representative\n"
                   "of a characteristic length of the geometry,\n"
                   "e.g. the MAC of a wing or diameter of a fuselage.\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    // Loop through surface meshes to set boundary layer parameters
    for (ibody = 0; ibody < aflr3Instance->numSurfaceMesh; ibody++) {
        // Allocate blThickness and blSpacing
        AIM_REALL(blThickness, aflr3Instance->surfaceMesh[ibody].numNode+nodeOffSet, double, aimInfo, status);
        AIM_REALL(blSpacing  , aflr3Instance->surfaceMesh[ibody].numNode+nodeOffSet, double, aimInfo, status);

        // Set default to globals
        for (i = 0; i < aflr3Instance->surfaceMesh[ibody].numNode; i++) {
            blThickness[i + nodeOffSet] = globalBLThickness*capsMeshLength;
            blSpacing  [i + nodeOffSet] = globalBLSpacing*capsMeshLength;
        }

        // Allocate blFlag
        AIM_REALL(blFlag, aflr3Instance->surfaceMesh[ibody].numElement+elementOffSet, int, aimInfo, status);

        // Set BL flag
        for (i = 0; i < aflr3Instance->surfaceMesh[ibody].numElement; i++) {

            // Set default value to non-viscous
            blFlag[i + elementOffSet] = 0;

            if (globalBLThickness != 0.0 &&
                globalBLSpacing   != 0.0) {

                if (ibody == boundingBoxIndex) {
                    if (i == 0) {
                        printf("Global boundary layer parameters will not be specified on surface %d as it is the bounding box!\n", ibody+1);
                    }
                    continue;
                }

                // Turn on BL flags if globals are set
                blFlag[i + elementOffSet] = -1;
            }
        }

        // Get the body for the surface mesh from the tessellation
        status = EG_statusTessBody(aflr3Instance->surfaceMesh[ibody].egadsTess,
                                   &body, &tempInt1, &tempInt2);
        AIM_STATUS(aimInfo, status);

        // Loop through meshing properties and see if boundaryLayerThickness and boundaryLayerSpacing have been specified
        for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

            // If no boundary layer specified in meshProp continue
            if (meshProp != NULL)
                if (meshProp[propIndex].boundaryLayerThickness == 0 ||
                    meshProp[propIndex].boundaryLayerSpacing   == 0) continue;

            // Loop through Elements and see if marker match
            for (i = 0; i < aflr3Instance->surfaceMesh[ibody].numElement; i++) {

                status = get_MeshIndex(&aflr3Instance->meshMap,
                                       body,
                                       &aflr3Instance->surfaceMesh[ibody].element[i],
                                       &meshIndex);
                AIM_STATUS(aimInfo, status);

                //If they don't match continue
                if (meshProp != NULL)
                    if (meshIndex!= meshProp[propIndex].attrIndex) continue;

                // Set face "bl" flag
                blFlag[i + elementOffSet] = -1;

                if (aflr3Instance->surfaceMesh[ibody].element[i].elementType != Triangle &&
                    aflr3Instance->surfaceMesh[ibody].element[i].elementType != Quadrilateral) continue;

                // Get face indexing for Triangles - 1 bias
                if (aflr3Instance->surfaceMesh[ibody].element[i].elementType == Triangle) {
                    pointIndex[0] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[0] -1;
                    pointIndex[1] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[1] -1;
                    pointIndex[2] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[2] -1;
                    pointIndex[3] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[2] -1; //Repeat last point for simplicity

                }

                if (aflr3Instance->surfaceMesh[ibody].element[i].elementType == Quadrilateral) {
                    pointIndex[0] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[0] -1;
                    pointIndex[1] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[1] -1;
                    pointIndex[2] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[2] -1;
                    pointIndex[3] = aflr3Instance->surfaceMesh[ibody].element[i].connectivity[3] -1;
                }

                // Set boundary layer spacing
                if (meshProp == NULL) continue;
                blSpacing[pointIndex[0] + nodeOffSet] = meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
                blSpacing[pointIndex[1] + nodeOffSet] = meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
                blSpacing[pointIndex[2] + nodeOffSet] = meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
                blSpacing[pointIndex[3] + nodeOffSet] = meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;

                // Set boundary layer thickness
                blThickness[pointIndex[0] + nodeOffSet] = meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
                blThickness[pointIndex[1] + nodeOffSet] = meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
                blThickness[pointIndex[2] + nodeOffSet] = meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
                blThickness[pointIndex[3] + nodeOffSet] = meshProp[propIndex].boundaryLayerThickness*capsMeshLength;
            }
        }

        elementOffSet += aflr3Instance->surfaceMesh[ibody].numElement;
        nodeOffSet    += aflr3Instance->surfaceMesh[ibody].numNode;
    }

    // Create/setup volume meshes
    if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {

        AIM_ALLOC(volumeMesh, numBody, meshStruct, aimInfo, status);
        numVolumeMesh = numBody;

        for (ibody = 0; ibody < numVolumeMesh; ibody++) {
            // Set reference mesh - One surface per body
            AIM_ALLOC(volumeMesh[ibody].referenceMesh, 1, meshStruct, aimInfo, status);
            volumeMesh[ibody].numReferenceMesh = 1;

            volumeMesh[ibody].referenceMesh[0] = aflr3Instance->surfaceMesh[ibody];
        }

    } else {

        AIM_ALLOC(volumeMesh, 1, meshStruct, aimInfo, status);
        numVolumeMesh = 1;

        status = initiate_meshStruct(&volumeMesh[0]);
        AIM_STATUS(aimInfo, status);

        // Combine mesh - temporary store the combined mesh in the volume mesh
        status = mesh_combineMeshStruct(aflr3Instance->numSurfaceMesh,
                                        aflr3Instance->surfaceMesh,
                                        &volumeMesh[0]);
        AIM_STATUS(aimInfo, status);

        // Set reference meshes - All surfaces
        volumeMesh[0].numReferenceMesh = aflr3Instance->numSurfaceMesh;
        AIM_ALLOC(volumeMesh[0].referenceMesh, volumeMesh[0].numReferenceMesh, meshStruct, aimInfo, status);

        for (ibody = 0; ibody < aflr3Instance->numSurfaceMesh; ibody++) {
            volumeMesh[0].referenceMesh[ibody] = aflr3Instance->surfaceMesh[ibody];
        }

        // Report surface mesh
        printf("Number of surface nodes - %d\n",
               volumeMesh[0].numNode);
        printf("Number of surface elements - %d\n",
               volumeMesh[0].numElement);

    }

    // Run AFLR3
    elementOffSet = 0;
    nodeOffSet    = 0;
    for (i = 0; i < numVolumeMesh; i++) {

        snprintf(bodyNumberFile, 42, "aflr3_%d", i);
        status = aim_file(aimInfo, bodyNumberFile, aimFile);
        AIM_STATUS(aimInfo, status);

        if (numVolumeMesh > 1) {
            AIM_NOTNULL(blFlag, aimInfo, status);
            AIM_NOTNULL(blSpacing, aimInfo, status);
            AIM_NOTNULL(blThickness, aimInfo, status);
            printf("Getting volume mesh for body %d (of %d)\n",
                   i+1, numVolumeMesh);
/*@-nullpass@*/
            status = aflr3_Volume_Mesh(aimInfo, aimInputs,
                                       aflr3Instance->meshInput,
                                       aimFile,
                                       createBL,
                                       &blFlag[elementOffSet],
                                       &blSpacing[nodeOffSet],
                                       &blThickness[nodeOffSet],
                                       capsMeshLength,
                                       numMeshProp,
                                       meshProp,
                                       &volumeMesh[i].referenceMesh[0],
                                       &volumeMesh[i]);
/*@+nullpass@*/
            elementOffSet += volumeMesh[i].referenceMesh[0].numElement;
            nodeOffSet    += volumeMesh[i].referenceMesh[0].numNode;
        } else {
            printf("Getting volume mesh\n");
/*@-nullpass@*/
            status = aflr3_Volume_Mesh(aimInfo, aimInputs,
                                       aflr3Instance->meshInput,
                                       aimFile,
                                       createBL,
                                       blFlag,
                                       blSpacing,
                                       blThickness,
                                       capsMeshLength,
                                       numMeshProp,
                                       meshProp,
                                       &volumeMesh[i],
                                       &volumeMesh[i]);
/*@+nullpass@*/
        }

        AIM_STATUS(aimInfo, status, "Problem during volume meshing of bodyIndex %d\n", i+1);

        if (numVolumeMesh > 1) {
            printf("Volume mesh for body %d (of %d):\n",
                   i+1, numVolumeMesh);
        } else {
            printf("Volume mesh:\n");
        }
        printf("\tNumber of nodes = %d\n",
               volumeMesh[i].numNode);
        printf("\tNumber of elements = %d\n",
               volumeMesh[i].numElement);
        printf("\tNumber of triangles = %d\n",
               volumeMesh[i].meshQuickRef.numTriangle);
        printf("\tNumber of quadrilatarals = %d\n",
               volumeMesh[i].meshQuickRef.numQuadrilateral);
        printf("\tNumber of tetrahedrals = %d\n",
               volumeMesh[i].meshQuickRef.numTetrahedral);
        printf("\tNumber of pyramids= %d\n",
               volumeMesh[i].meshQuickRef.numPyramid);
        printf("\tNumber of prisms= %d\n",
               volumeMesh[i].meshQuickRef.numPrism);
        printf("\tNumber of hexahedrals= %d\n",
               volumeMesh[i].meshQuickRef.numHexahedral);
    }

    for (i = 0; i < numVolumeMesh; i++) {

        // Check to make sure the volume mesher didn't add any unaccounted for points/faces
        numElementCheck = 0;
        for (j = 0; j < volumeMesh[i].numReferenceMesh; j++) {
          numElementCheck += volumeMesh[i].referenceMesh[j].numElement;
        }

        if (volumeMesh[i].meshQuickRef.useStartIndex == (int) false &&
            volumeMesh[i].meshQuickRef.useListIndex  == (int) false) {

          status = mesh_retrieveNumMeshElements(volumeMesh[i].numElement,
                                                volumeMesh[i].element,
                                                Triangle,
                                                &volumeMesh[i].meshQuickRef.numTriangle);
          AIM_STATUS(aimInfo, status);

          status = mesh_retrieveNumMeshElements(volumeMesh[i].numElement,
                                                volumeMesh[i].element,
                                                Quadrilateral,
                                                &volumeMesh[i].meshQuickRef.numQuadrilateral);
          AIM_STATUS(aimInfo, status);

        }

        snprintf(bodyNumberFile, 42, NODATATRANSFER, i);
        status = aim_rmFile(aimInfo, bodyNumberFile);
        AIM_STATUS(aimInfo, status);

        if (numElementCheck != volumeMesh[i].meshQuickRef.numTriangle +
                               volumeMesh[i].meshQuickRef.numQuadrilateral) {

            fp = aim_fopen(aimInfo, bodyNumberFile,"w");
            if (fp == NULL) {
                AIM_ERROR(aimInfo, "Failed to open '%s'", bodyNumberFile);
                status = CAPS_IOERR;
                goto cleanup;
            }
            fprintf(fp, "shucks...");
            fclose(fp);  fp = NULL;

            printf("Volume mesher added surface elements - data transfer will NOT be possible.\n");
        }
    }


    // If filename is not NULL write the mesh
    if (aflr3Instance->meshInput.outputFileName != NULL) {

        for (ibody = 0; ibody < numVolumeMesh; ibody++) {

            if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {
                sprintf(bodyNumberFile, "%d", ibody);
                AIM_ALLOC(filename, strlen(aflr3Instance->meshInput.outputFileName) +
                                    2+strlen("_Vol")+strlen(bodyNumberFile), char, aimInfo, status);
            } else {
                AIM_ALLOC(filename, strlen(aflr3Instance->meshInput.outputFileName) +
                                           2, char, aimInfo, status);
            }

            strcpy(filename, aflr3Instance->meshInput.outputFileName);

            if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {
                strcat(filename,"_Vol");
                strcat(filename,bodyNumberFile);
            }

            if (strcasecmp(aflr3Instance->meshInput.outputFormat, "AFLR3") == 0) {

                status = mesh_writeAFLR3(aimInfo,
                                         filename,
                                         aflr3Instance->meshInput.outputASCIIFlag,
                                         &volumeMesh[ibody],
                                         1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(aflr3Instance->meshInput.outputFormat, "VTK") == 0) {

                status = mesh_writeVTK(aimInfo,
                                       filename,
                                       aflr3Instance->meshInput.outputASCIIFlag,
                                       &volumeMesh[ibody],
                                       1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(aflr3Instance->meshInput.outputFormat, "SU2") == 0) {

                status = mesh_writeSU2(aimInfo,
                                       filename,
                                       aflr3Instance->meshInput.outputASCIIFlag,
                                       &volumeMesh[ibody],
                                       aflr3Instance->meshInput.bndConds.numBND,
                                       aflr3Instance->meshInput.bndConds.bndID,
                                       1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(aflr3Instance->meshInput.outputFormat, "Tecplot") == 0) {

                status = mesh_writeTecplot(aimInfo,
                                           filename,
                                           aflr3Instance->meshInput.outputASCIIFlag,
                                           &volumeMesh[ibody],
                                           1.0);
                AIM_STATUS(aimInfo, status);


            } else if (strcasecmp(aflr3Instance->meshInput.outputFormat, "Nastran") == 0) {

                status = mesh_writeNASTRAN(aimInfo,
                                           filename,
                                           aflr3Instance->meshInput.outputASCIIFlag,
                                           &volumeMesh[ibody],
                                           LargeField,
                                           1.0);
                AIM_STATUS(aimInfo, status);

            } else {
                printf("Unrecognized mesh format, \"%s\", the volume mesh will not be written out\n",
                       aflr3Instance->meshInput.outputFormat);
            }

            if (filename != NULL) EG_free(filename);
            filename = NULL;
        }
    } else {
        printf("No project name (\"Proj_Name\") provided - A volume mesh will not be written out\n");
    }

    status = CAPS_SUCCESS;

cleanup:

    // Destroy volumeMesh allocated arrays
    if (volumeMesh != NULL) {
        for (i = 0; i < numVolumeMesh; i++) {
            (void) destroy_meshStruct(&volumeMesh[i]);
        }
    }
    AIM_FREE(volumeMesh);

    // Clean up meshProps
    if (meshProp != NULL) {
        for (i = 0; i < numMeshProp; i++) {
            (void) destroy_meshSizingStruct(&meshProp[i]);
        }
        AIM_FREE(meshProp);
    }

    // Free filename
    AIM_FREE(filename);

    // Free boundary layer parameters
    AIM_FREE(blFlag);
    AIM_FREE(blThickness);
    AIM_FREE(blSpacing);  // Boundary layer spacing [numNode] - Total nodes for all bodies

    return status;
}


/* the execution code from above should be moved here */
int aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(void *instStore, void *aimInfo,
                    /*@unused@*/ int restart, capsValue *aimInputs)
{
    int status = CAPS_SUCCESS;

    // Incoming bodies
    const char *intents;
    ego *bodies = NULL;
    int numBody = 0;

    // Container for volume mesh
    int numSurfaceMesh=0;
    meshStruct *surfaceMesh=NULL;

    int i, ibody, nodeOffSet;
    int noDataTransfer = (int)false;

    char bodyNumberFile[42];
    char aimFile[PATH_MAX];

    aimStorage *aflr3Instance;

    aflr3Instance = (aimStorage *) instStore;

    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    // Get mesh
    numSurfaceMesh = aimInputs[Surface_Mesh-1].length;
    surfaceMesh    = (meshStruct *)aimInputs[Surface_Mesh-1].vals.AIMptr;
    AIM_NOTNULL(surfaceMesh, aimInfo, status);

    // Create/setup volume meshes
    if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {

        AIM_ALLOC(aflr3Instance->meshRef, numBody, aimMeshRef, aimInfo, status);
        aflr3Instance->numMeshRef = numBody;

        for (ibody = 0; ibody < numBody; ibody++) {
          status = aim_initMeshRef(aflr3Instance->meshRef);
          AIM_STATUS(aimInfo, status);
        }

        for (ibody = 0; ibody < numBody; ibody++) {
            snprintf(bodyNumberFile, 42, "aflr3_%d", ibody);
            status = aim_file(aimInfo, bodyNumberFile, aimFile);
            AIM_STATUS(aimInfo, status);
            AIM_STRDUP(aflr3Instance->meshRef[ibody].fileName, aimFile, aimInfo, status);

            aflr3Instance->meshRef[ibody].maps[0].tess = surfaceMesh[ibody].egadsTess;

            snprintf(bodyNumberFile, 42, NODATATRANSFER, ibody);
            if (aim_isFile(aimInfo, bodyNumberFile) == CAPS_SUCCESS) continue;

            aflr3Instance->meshRef[0].nmap = 1;
            AIM_ALLOC(aflr3Instance->meshRef[ibody].maps[0].map, surfaceMesh[ibody].numNode, int, aimInfo, status);
            for (i = 0; i < surfaceMesh[ibody].numNode; i++)
              aflr3Instance->meshRef[ibody].maps[0].map[i] = i+1;
        }

    } else {

        AIM_ALLOC(aflr3Instance->meshRef, 1, aimMeshRef, aimInfo, status);
        aflr3Instance->numMeshRef = 1;
        ibody = 0;

        status = aim_initMeshRef(aflr3Instance->meshRef);
        AIM_STATUS(aimInfo, status);

        // set the filename without extensions where the grid is written for solvers
        snprintf(bodyNumberFile, 42, "aflr3_%d", ibody);
        status = aim_file(aimInfo, bodyNumberFile, aimFile);
        AIM_STATUS(aimInfo, status);
        AIM_STRDUP(aflr3Instance->meshRef[0].fileName, aimFile, aimInfo, status);

        snprintf(bodyNumberFile, 42, NODATATRANSFER, ibody);
        if (aim_isFile(aimInfo, bodyNumberFile) == CAPS_SUCCESS) {
            noDataTransfer = (int)true;
        }

        AIM_ALLOC(aflr3Instance->meshRef[0].maps, numSurfaceMesh, aimMeshTessMap, aimInfo, status);
        aflr3Instance->meshRef[0].nmap = numSurfaceMesh;

        nodeOffSet = 0;
        for (ibody = 0; ibody < numSurfaceMesh; ibody++) {
            aflr3Instance->meshRef[0].maps[ibody].tess = surfaceMesh[ibody].egadsTess;
            aflr3Instance->meshRef[0].maps[ibody].map = NULL;

            if (noDataTransfer == (int)true) continue;

            AIM_ALLOC(aflr3Instance->meshRef[0].maps[ibody].map, surfaceMesh[ibody].numNode, int, aimInfo, status);
            for (i = 0; i < surfaceMesh[ibody].numNode; i++)
              aflr3Instance->meshRef[0].maps[ibody].map[i] = nodeOffSet + i+1;

            nodeOffSet += surfaceMesh[ibody].numNode;
        }
    }

cleanup:
    return status;
}


int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               /*@unused@*/ int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsAFLR3 AIM Outputs
     * The following list outlines the AFLR3 AIM outputs available through the AIM interface.
     */

    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" aflr3AIM/aimOutputs index = %d!\n", index);
#endif

    if (index == NumberOfElement) {
        *aoname = EG_strdup("NumberOfElement");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsAFLR3
         * - <B> NumberOfElement </B> <br>
         * Number of elements in the surface mesh
         */

    } else if (index == NumberOfNode) {
        *aoname = EG_strdup("NumberOfNode");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsAFLR3
         * - <B> NumberOfNode </B> <br>
         * Number of vertices in the surface mesh
         */

    } else if (index == Volume_Mesh) {

        *aoname           = AIM_NAME(Volume_Mesh);
        form->type        = PointerMesh;
        form->dim         = Vector;
        form->lfixed      = Change;
        form->sfixed      = Fixed;
        form->vals.AIMptr = NULL;
        form->nullVal     = IsNull;

        /*! \page aimOutputsAFLR3
         * - <B> Volume_Mesh </B> <br>
         * The volume mesh for a link
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimStruc, status, "Unknown output index %d!", index);
    }

    AIM_NOTNULL(*aoname, aimStruc, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*aoname);
    return status;
}


int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ int index, capsValue *val)
{
    int        i, status = CAPS_SUCCESS;
    int        numElement, numNodes;
    int        nVertex, nTri, nQuad, nTet, nPyramid, nPrism, nHex;
    aimStorage *aflr3Instance;
    aimMesh    mesh;

#ifdef DEBUG
    printf(" aflr3AIM/aimCalcOutput index = %d!\n", index);
#endif
    aflr3Instance = (aimStorage *) instStore;

    if (NumberOfElement == index) {

        // Count the total number of surface elements
        numElement = 0;
        for (i = 0; i < aflr3Instance->numMeshRef; i++) {
            status = aim_readBinaryUgridHeader(aimInfo, &aflr3Instance->meshRef[i],
                                               &nVertex, &nTri, &nQuad,
                                               &nTet, &nPyramid, &nPrism, &nHex);
            AIM_STATUS(aimInfo, status);

            numElement += nTri + nQuad + nTet + nPyramid + nPrism + nHex;
        }

        val->vals.integer = numElement;

    } else if (NumberOfNode == index) {

        // Count the total number of surface vertices
        numNodes = 0;
        for (i = 0; i < aflr3Instance->numMeshRef; i++) {
             status = aim_readBinaryUgridHeader(aimInfo, &aflr3Instance->meshRef[i],
                                                &nVertex, &nTri, &nQuad,
                                                &nTet, &nPyramid, &nPrism, &nHex);
             AIM_STATUS(aimInfo, status);

             numNodes += nVertex;
        }

        val->vals.integer = numNodes;

    } else if (Volume_Mesh == index) {

        for (i = 0; i < aflr3Instance->numMeshRef; i++) {
            status = aim_queryMeshes( aimInfo, Volume_Mesh, &aflr3Instance->meshRef[i] );
            if (status > 0) {
/*@-immediatetrans@*/
                mesh.meshData = NULL;
                mesh.meshRef = &aflr3Instance->meshRef[i];
/*@+immediatetrans@*/

                status = aim_readBinaryUgrid(aimInfo, &mesh);
                AIM_STATUS(aimInfo, status);

                status = aim_writeMeshes(aimInfo, Volume_Mesh, &mesh);
                AIM_STATUS(aimInfo, status);

                status = aim_freeMeshData(mesh.meshData);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(mesh.meshData);
            }
            else
              AIM_STATUS(aimInfo, status);
        }

/*@-immediatetrans@*/
        // Return the volume mesh references
        val->nrow        = aflr3Instance->numMeshRef;
        val->vals.AIMptr = aflr3Instance->meshRef;
/*@+immediatetrans@*/

    } else {

        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

    }

cleanup:

    return status;
}


void aimCleanup(void *instStore)
{
    int        status;
    aimStorage *aflr3Instance;

#ifdef DEBUG
    printf(" aflr3AIM/aimCleanup!\n");
#endif
    aflr3Instance = (aimStorage *) instStore;

    // Clean up aflr3Instance data
    status = destroy_aimStorage(aflr3Instance);
    if (status != CAPS_SUCCESS)
      printf(" Status = %d, aflr3AIM aimStorage cleanup!!!\n", status);

    EG_free(aflr3Instance);
}

