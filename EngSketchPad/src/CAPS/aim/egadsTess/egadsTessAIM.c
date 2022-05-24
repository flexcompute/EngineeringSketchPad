/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             EGADS Tessellation AIM
 *
 *      This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */

/*!\mainpage Introduction
 *
 * \section overviewEgadsTess EGADS Tessellation AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact internal meshing
 * capability for the EGADS library.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsEgadsTess and \ref aimOutputsEgadsTess, respectively.
 *
 * \section clearanceEgadsTess Clearance Statement
 *  This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 *
 */



#include <string.h>
#include <math.h>
#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"       // Collection of helper functions for meshing
#include "miscUtils.h"
#include "deprecateUtils.h"

#ifdef WIN32
#define strcasecmp  stricmp
#define strncasecmp _strnicmp
#endif


enum aimInputs
{
  Proj_Name = 1,               /* index is 1-based */
  Mesh_Quiet_Flag,
  Mesh_Length_Factor,
  Tess_Params,
  Mesh_Format,
  Mesh_ASCII_Flag,
  Edge_Point_Min,
  Edge_Point_Max,
  Mesh_Sizing,
  Mesh_Elements,
  Multiple_Mesh,
  TFI_Templates,
  NUMINPUT = TFI_Templates     /* Total number of inputs */
};

enum aimOutputs
{
  Done = 1,                    /* index is 1-based */
  NumberOfElement,
  NumberOfNode,
  Surface_Mesh,
  NUMOUT = Surface_Mesh        /* Total number of outputs */
};


#define MXCHAR  255
#define EGADSTESSFILE "egadsTess_%d.eto"

//#define DEBUG

typedef struct {

    // Container for surface mesh
    int numSurface;
    meshStruct *surfaceMesh;

    // Container for mesh input
    meshInputStruct meshInput;

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    mapAttrToIndexStruct meshMap;

} aimStorage;


static int destroy_aimStorage(aimStorage *egadsInstance)
{

    int i; // Indexing

    int status; // Function return status

    // Destroy meshInput
    status = destroy_meshInputStruct(&egadsInstance->meshInput);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, egadsTessAIM meshInput cleanup!!!\n", status);

    // Destroy surface mesh allocated arrays
    for (i = 0; i < egadsInstance->numSurface; i++) {
        status = destroy_meshStruct(&egadsInstance->surfaceMesh[i]);
        if (status != CAPS_SUCCESS)
          printf("Status = %d, egadsTessAIM surfaceMesh cleanup!!!\n", status);
    }
    egadsInstance->numSurface = 0;
    AIM_FREE(egadsInstance->surfaceMesh);

    // Destroy attribute to index map
    status = destroy_mapAttrToIndexStruct(&egadsInstance->groupMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, egadsTessAIM attributeMap cleanup!!!\n", status);

    status = destroy_mapAttrToIndexStruct(&egadsInstance->meshMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, egadsTessAIM attributeMap cleanup!!!\n", status);
    return CAPS_SUCCESS;
}


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status; // Function return status

    aimStorage *egadsInstance=NULL;

#ifdef DEBUG
    printf("\n egadsTessAIM/aimInitialize   inst = %d!\n", inst);
#endif

    // Specify the number of analysis input and out "parameters"
    *nIn     = NUMINPUT;
    *nOut    = NUMOUT;
    if (inst == -1) return CAPS_SUCCESS;

    /* specify the field variables this analysis can generate and consume */
    *nFields = 0;
    *fnames  = NULL;
    *franks  = NULL;
    *fInOut  = NULL;

    // Allocate egadsInstance
    AIM_ALLOC(egadsInstance, 1, aimStorage, aimInfo, status);
    *instStore = egadsInstance;

    // Set initial values for egadsInstance //

    // Container for surface meshes
    egadsInstance->surfaceMesh = NULL;
    egadsInstance->numSurface = 0;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&egadsInstance->groupMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&egadsInstance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Container for mesh input
    status = initiate_meshInputStruct(&egadsInstance->meshInput);
    AIM_STATUS(aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);

    return status;
}


int aimInputs(/*@unused@*/ void *aimStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{

    /*! \page aimInputsEgadsTess AIM Inputs
     * The following list outlines the EGADS Tessellation meshing options along with their default value available
     * through the AIM interface.
     */

#ifdef DEBUG
    printf(" egadsTessAIM/aimInputs index = %d!\n", index);
#endif

    // Meshing Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name"); // If NULL a volume grid won't be written by the AIM
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->vals.string  = NULL;
        //defval->vals.string  = EG_strdup("CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsEgadsTess
         * - <B> Proj_Name = NULL</B> <br>
         * This corresponds to the output name of the mesh. If left NULL, the mesh is not written to a file.
         */

    } else if (index == Mesh_Quiet_Flag) {
        *ainame               = AIM_NAME(Mesh_Quiet_Flag);
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsAFLR4
         * - <B> Mesh_Quiet_Flag = False</B> <br>
         * Complete suppression of mesh generator (not including errors)
         */

    } else if (index == Mesh_Length_Factor) {
        *ainame              = EG_strdup("Mesh_Length_Factor");
        defval->type         = Double;
        defval->dim          = Scalar;
        defval->vals.real    = 1;
        defval->nullVal      = NotNull;

        /*! \page aimInputsEgadsTess
         * - <B> Mesh_Length_Factor = 1</B> <br>
         * Scaling factor to compute a meshing Reference_Length via:<br>
         * <br>
         * Reference_Length = capsMeshLength*Mesh_Length_Factor<br>
         * <br>
         * Reference_Length scales Tess_Params[0] and Tess_Params[1] in both aimInputs and Mesh_Sizing
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
            defval->vals.reals[0] = 0.10;
            defval->vals.reals[1] = 0.01;
            defval->vals.reals[2] = 15.0;
        } else return EGADS_MALLOC;

        /*! \page aimInputsEgadsTess
         * - <B> Tess_Params = [0.1, 0.01, 15.0]</B> <br>
         * Body tessellation parameters. Tess_Params[0] and Tess_Params[1] get scaled by Reference_Length if
         * it is set, otherwise by the bounding box of the largest body. (From the EGADS manual)
         * A set of 3 parameters that drive the EDGE discretization
         * and the FACE triangulation. The first is the maximum length of an EDGE segment or triangle side
         * (in physical space). A zero is flag that allows for any length. The second is a curvature-based
         * value that looks locally at the deviation between the centroid of the discrete object and the
         * underlying geometry. Any deviation larger than the input value will cause the tessellation to
         * be enhanced in those regions. The third is the maximum interior dihedral angle (in degrees)
         * between triangle facets (or Edge segment tangents for a WIREBODY tessellation), note that a
         * zero ignores this phase
         */
    } else if (index == Mesh_Format) {
        *ainame               = EG_strdup("Mesh_Format");
        defval->type          = String;
        defval->vals.string   = EG_strdup("AFLR3"); // TECPLOT, VTK, AFLR3, STL, AF, FAST, NASTRAN
        defval->lfixed        = Change;

        /*! \page aimInputsEgadsTess
         * - <B> Mesh_Format = "AFLR3"</B> <br>
         * Mesh output format. Available format names include: "AFLR3", "VTK", "TECPLOT", "STL" (quadrilaterals will be
         * split into triangles), "Airfoil", "FAST", "Nastran".
         *
         *   - "Airfoil" corresponds to the following file format in which the nodal coordinates of the body's edges
         *  are written. Bodies should be face bodies, planar, and have no holes. A *.af suffix is used for the file:
         * <br>"Character Name"
         * <br>x[0] y[0] x y coordinates
         * <br>x[1] y[1]
         * <br>...  ...
         */
    } else if (index == Mesh_ASCII_Flag) {
        *ainame               = EG_strdup("Mesh_ASCII_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = true;

        /*! \page aimInputsEgadsTess
         * - <B> Mesh_ASCII_Flag = True</B> <br>
         * Output mesh in ASCII format, otherwise write a binary file if applicable.
         */

    } else if (index == Edge_Point_Min) {
        *ainame               = EG_strdup("Edge_Point_Min");
        defval->type          = Integer;
        defval->vals.integer  = 0;
        defval->lfixed        = Fixed;
        defval->nrow          = 1;
        defval->ncol          = 1;
        defval->nullVal       = IsNull;

        /*! \page aimInputsEgadsTess
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

        /*! \page aimInputsEgadsTess
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

        /*! \page aimInputsEgadsTess
         * - <B>Mesh_Sizing = NULL </B> <br>
         * See \ref meshSizingProp for additional details.
         */

    } else if (index == Mesh_Elements) {
        *ainame              = EG_strdup("Mesh_Elements");
        defval->type         = String;
        defval->nullVal      = NotNull;
        defval->vals.string  = EG_strdup("Tri");
        defval->lfixed       = Change;

        /*! \page aimInputsEgadsTess
         * - <B>Mesh_Elements = "Tri" </B> <br>
         * Element topology in the resulting mesh:
         *  - "Tri"   - All triangle elements
         *  - "Quad"  - All quadrilateral elements
         *  - "Mixed" - Quad elements for four-sided faces with TFI, triangle elements otherwise (deprecated)
         */

    } else if (index == Multiple_Mesh) {
        *ainame               = EG_strdup("Multiple_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = (int) true;

        /*! \page aimInputsEgadsTess
         * - <B> Multiple_Mesh = True</B> <br>
         * If set to True (default) a surface mesh will be generated and output (given "Proj_Name" is set) for each body.
         * When set to False only a single surface
         * mesh will be created. Note, this only affects the mesh when writing to a
         * file.
         */

    } else if (index == TFI_Templates) {
        *ainame               = EG_strdup("TFI_Templates");
        defval->type          = Boolean;
        defval->vals.integer  = (int) true;

        /*! \page aimInputsEgadsTess
         * - <B> TFI_Templates = True</B> <br>
         * Use Transfinite Interpolation and Templates to generate
         * structured triangulations on FACEs with 3 or 4 "sides" with similar opposing vertex counts.
         */
    }

    return CAPS_SUCCESS;
}


int aimPreAnalysis(void *aimStore, void *aimInfo, capsValue *aimInputs)
{
    int status; // Status return

    int i, bodyIndex; // Indexing

    aimStorage *egadsInstance;

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects
    ego etess = NULL;

    // Global settings
    int minEdgePoint = -1, maxEdgePoint = -1, quadMesh = 0;
    double refLen = 0, meshLenFac = 0, capsMeshLength = 0;

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;
    const char *MeshElements = NULL;

    // File output
    char bodyNumber[42];
    char aimFile[PATH_MAX];

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" egadsTessAIM/aimPreAnalysis numBody = %d!\n", numBody);
#endif

    if ((numBody <= 0) || (bodies == NULL)) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }

    if (aimInputs == NULL) {
#ifdef DEBUG
        printf(" egadsTessAIM/aimPreAnalysis No aimInputs!\n");
#endif
        return CAPS_SOURCEERR;
    }
    egadsInstance = (aimStorage *) aimStore;

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(egadsInstance);
    if (status != CAPS_SUCCESS) {
        AIM_ERROR(aimInfo, "Status = %d, egadsTessAIM aimStorage cleanup!!!", status);
        return status;
    }

    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            3,
                                            &egadsInstance->groupMap);
    AIM_STATUS(aimInfo, status);

    status = create_CAPSMeshAttrToIndexMap(numBody,
                                           bodies,
                                           3,
                                           &egadsInstance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Setup meshing input structure

    // Reference length for meshing
    meshLenFac = aimInputs[Mesh_Length_Factor-1].vals.real;

    status = check_CAPSMeshLength(numBody, bodies, &capsMeshLength);

    // TODO: Should capsMeshLength be optional?
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
    }
    */

    if (meshLenFac <= 0) {
        AIM_ERROR(aimInfo, "Mesh_Length_Factor is: %f\n", meshLenFac);
        AIM_ADDLINE(aimInfo, "Mesh_Length_Factor must be a positive number.");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    refLen = meshLenFac*capsMeshLength;

    // Get Tessellation parameters -Tess_Params
    egadsInstance->meshInput.paramTess[0] =
                 aimInputs[Tess_Params-1].vals.reals[0]; // Gets multiplied by bounding box size
    egadsInstance->meshInput.paramTess[1] =
                 aimInputs[Tess_Params-1].vals.reals[1]; // Gets multiplied by bounding box size
    egadsInstance->meshInput.paramTess[2] =
                 aimInputs[Tess_Params-1].vals.reals[2];

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
      AIM_ERROR(aimInfo, "Edge_Point_Max must be greater or equal Edge_Point_Min\n");
      AIM_ERROR(aimInfo, "Edge_Point_Max = %d, Edge_Point_Min = %d\n", maxEdgePoint, minEdgePoint);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // Get mesh sizing parameters
    if (aimInputs[Mesh_Sizing-1].nullVal != IsNull) {

        status = deprecate_SizingAttr(aimInfo,
                                      aimInputs[Mesh_Sizing-1].length,
                                      aimInputs[Mesh_Sizing-1].vals.tuple,
                                      &egadsInstance->meshMap,
                                      &egadsInstance->groupMap);
        AIM_STATUS(aimInfo, status);

        status = mesh_getSizingProp(aimInfo,
                                    aimInputs[Mesh_Sizing-1].length,
                                    aimInputs[Mesh_Sizing-1].vals.tuple,
                                    &egadsInstance->meshMap,
                                    &numMeshProp,
                                    &meshProp);
        AIM_STATUS(aimInfo, status);
    }

    // Get mesh element types
    MeshElements = aimInputs[Mesh_Elements-1].vals.string;

         if ( strncasecmp(MeshElements,"Tri",3)   == 0 ) { quadMesh = 0; }
    else if ( strncasecmp(MeshElements,"Quad",4)  == 0 ) { quadMesh = 1; }
    else if ( strncasecmp(MeshElements,"Mixed",3) == 0 ) { quadMesh = 2; }
    else {
        AIM_ERROR(  aimInfo, "Unknown Mesh_Elements = \"%s\"\n", MeshElements);
        AIM_ADDLINE(aimInfo, "       Should be one of \"Tri\", \"Quad\", or \"Mixed\"\n");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (aimInputs[TFI_Templates-1].vals.integer == (int) false) {
        for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++){
            // Disable TFI and templates
            status = EG_attributeAdd(bodies[bodyIndex], ".qParams", ATTRSTRING,
                                     0, NULL, NULL, "off");
            if (status != EGADS_SUCCESS) goto cleanup;
        }
    }

    // Modify the EGADS body tessellation based on given inputs
/*@-nullpass@*/
    status =  mesh_modifyBodyTess(numMeshProp,
                                  meshProp,
                                  minEdgePoint,
                                  maxEdgePoint,
                                  quadMesh,
                                  &refLen,
                                  egadsInstance->meshInput.paramTess,
                                  egadsInstance->meshMap,
                                  numBody,
                                  bodies);
/*@+nullpass@*/
    AIM_STATUS(aimInfo, status);

    // Run egadsTess for each body
    for (bodyIndex = 0 ; bodyIndex < numBody; bodyIndex++) {

        if (aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false)
          printf("Getting surface mesh for body %d (of %d)\n", bodyIndex+1, numBody);

        status = mesh_surfaceMeshEGADSBody(aimInfo,
                                           bodies[bodyIndex],
                                           refLen,
                                           egadsInstance->meshInput.paramTess,
                                           quadMesh,
                                           &etess);
        AIM_STATUS(aimInfo, status, "Problem during surface meshing of body %d", bodyIndex+1);
        AIM_NOTNULL(etess, aimInfo, status);

        // set the file name to write the egads file
        snprintf(bodyNumber, 42, EGADSTESSFILE, bodyIndex);
        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        remove(aimFile);
        status = EG_saveTess(etess, aimFile);
        AIM_STATUS(aimInfo, status);

        EG_deleteObject(etess);
    }

    status = CAPS_SUCCESS;

cleanup:

    // Clean up meshProps
    if (meshProp != NULL) {
        for (i = 0; i < numMeshProp; i++) {
            (void) destroy_meshSizingStruct(&meshProp[i]);
        }
        AIM_FREE(meshProp);
    }
    return status;
}


/* the execution code from above should be moved here */
int aimExecute(/*@unused@*/ void *aimStore, /*@unused@*/ void *aimStruc, int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *aimStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart,    /*@unused@*/ capsValue *aimInputs)
{

    int status = CAPS_SUCCESS;
    int bodyIndex;

    int numNodeTotal=0, numElemTotal=0;

    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    const char *intents;
    char bodyNumber[42];
    char aimFile[PATH_MAX];
    char *filename = NULL;

    aimStorage *egadsInstance;
    egadsInstance = (aimStorage *) aimStore;

    // Combined mesh
    meshStruct combineMesh;

    status = initiate_meshStruct(&combineMesh);
    AIM_STATUS(aimInfo, status);

    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodies, aimInfo, status);

    // Allocate surfaceMesh from number of bodies
    AIM_ALLOC(egadsInstance->surfaceMesh, numBody, meshStruct, aimInfo, status);
    egadsInstance->numSurface = numBody;

    // Initiate surface meshes
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++){
        status = initiate_meshStruct(&egadsInstance->surfaceMesh[bodyIndex]);
        AIM_STATUS(aimInfo, status);
    }

    if (egadsInstance->groupMap.mapName == NULL) {
        // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
        status = create_CAPSGroupAttrToIndexMap(numBody,
                                                bodies,
                                                3,
                                                &egadsInstance->groupMap);
        AIM_STATUS(aimInfo, status);
    }

    // Run egadsTess for each body
    for (bodyIndex = 0 ; bodyIndex < numBody; bodyIndex++) {

        status = copy_mapAttrToIndexStruct( &egadsInstance->groupMap,
                                            &egadsInstance->surfaceMesh[bodyIndex].groupMap );
        AIM_STATUS(aimInfo, status);

        // set the file name to read the egads file
        snprintf(bodyNumber, 42, EGADSTESSFILE, bodyIndex);
        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        status = EG_loadTess(bodies[bodyIndex], aimFile, &egadsInstance->surfaceMesh[bodyIndex].egadsTess);
        AIM_STATUS(aimInfo, status);

        status = mesh_surfaceMeshEGADSTess(aimInfo, &egadsInstance->surfaceMesh[bodyIndex]);
        AIM_STATUS(aimInfo, status);

        status = aim_newTess(aimInfo, egadsInstance->surfaceMesh[bodyIndex].egadsTess);
        AIM_STATUS(aimInfo, status);

        if (restart == 0 &&
            aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false) {
            printf("Body %d (of %d)\n", bodyIndex+1, numBody);

            printf("Number of nodes    = %d\n", egadsInstance->surfaceMesh[bodyIndex].numNode);
            printf("Number of elements = %d\n", egadsInstance->surfaceMesh[bodyIndex].numElement);

            if (egadsInstance->surfaceMesh[bodyIndex].meshQuickRef.useStartIndex == (int) true ||
                egadsInstance->surfaceMesh[bodyIndex].meshQuickRef.useListIndex  == (int) true) {

                printf("Number of node elements          = %d\n",
                       egadsInstance->surfaceMesh[bodyIndex].meshQuickRef.numNode);
                printf("Number of line elements          = %d\n",
                       egadsInstance->surfaceMesh[bodyIndex].meshQuickRef.numLine);
                printf("Number of triangle elements      = %d\n",
                       egadsInstance->surfaceMesh[bodyIndex].meshQuickRef.numTriangle);
                printf("Number of quadrilateral elements = %d\n",
                       egadsInstance->surfaceMesh[bodyIndex].meshQuickRef.numQuadrilateral);
            }

            numNodeTotal += egadsInstance->surfaceMesh[bodyIndex].numNode;
            numElemTotal += egadsInstance->surfaceMesh[bodyIndex].numElement;
        }
    }

    if (restart == 0 &&
        aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false) {
        printf("----------------------------\n");
        printf("Total number of nodes    = %d\n", numNodeTotal);
        printf("Total number of elements = %d\n", numElemTotal);
    }

    if (restart == 0) {
        // Project Name
        if (aimInputs[Proj_Name-1].nullVal != IsNull) {
            AIM_FREE(egadsInstance->meshInput.outputFileName);
            AIM_STRDUP(egadsInstance->meshInput.outputFileName, aimInputs[Proj_Name-1].vals.string, aimInfo, status);
        }

        // Mesh Format
        AIM_FREE(egadsInstance->meshInput.outputFormat);
        AIM_STRDUP(egadsInstance->meshInput.outputFormat, aimInputs[Mesh_Format-1].vals.string, aimInfo, status);

        // ASCII flag
        egadsInstance->meshInput.outputASCIIFlag = aimInputs[Mesh_ASCII_Flag-1].vals.integer;

        if (egadsInstance->meshInput.outputFileName != NULL) {

            // We need to combine the mesh
            if (aimInputs[Multiple_Mesh-1].vals.integer == (int) false) {

                status = mesh_combineMeshStruct(egadsInstance->numSurface,
                                                egadsInstance->surfaceMesh,
                                                &combineMesh);

                AIM_STATUS(aimInfo, status);

                filename = (char *) EG_alloc((strlen(egadsInstance->meshInput.outputFileName) +
                                              2)*sizeof(char));

                if (filename == NULL) goto cleanup;

                strcpy(filename, egadsInstance->meshInput.outputFileName);

                if (strcasecmp(egadsInstance->meshInput.outputFormat, "AFLR3") == 0) {

                    status = mesh_writeAFLR3(aimInfo, filename,
                                             egadsInstance->meshInput.outputASCIIFlag,
                                             &combineMesh,
                                             1.0);

                } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "VTK") == 0) {

                    status = mesh_writeVTK(aimInfo, filename,
                                           egadsInstance->meshInput.outputASCIIFlag,
                                           &combineMesh,
                                           1.0);

                } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "Tecplot") == 0) {

                    status = mesh_writeTecplot(aimInfo, filename,
                                               egadsInstance->meshInput.outputASCIIFlag,
                                               &combineMesh,
                                               1.0);

                } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "STL") == 0) {

                    status = mesh_writeSTL(aimInfo, filename,
                                           egadsInstance->meshInput.outputASCIIFlag,
                                           &combineMesh,
                                           1.0);

                } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "Airfoil") == 0) {

                    status = mesh_writeAirfoil(aimInfo, filename,
                                               egadsInstance->meshInput.outputASCIIFlag,
                                               &combineMesh,
                                               1.0);

                } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "FAST") == 0) {

                    status = mesh_writeFAST(aimInfo, filename,
                                            egadsInstance->meshInput.outputASCIIFlag,
                                            &combineMesh,
                                            1.0);

                } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "Nastran") == 0) {

                    status = mesh_writeNASTRAN(aimInfo, filename,
                                               egadsInstance->meshInput.outputASCIIFlag,
                                               &combineMesh,
                                               FreeField,
                                               1.0);
                } else {
                    printf("Unrecognized mesh format, \"%s\", the mesh will not be written out\n",
                           egadsInstance->meshInput.outputFormat);
                }

                AIM_FREE(filename);
                AIM_STATUS(aimInfo, status);

            } else {

                for (bodyIndex = 0; bodyIndex < egadsInstance->numSurface; bodyIndex++) {

                    if (egadsInstance->numSurface > 1) {
                        snprintf(bodyNumber, 42, "%d", bodyIndex);
                        filename = (char *) EG_alloc((strlen(egadsInstance->meshInput.outputFileName)  +
                                                      2 + strlen("_Surf_") + strlen(bodyNumber))*sizeof(char));
                    } else {
                        filename = (char *) EG_alloc((strlen(egadsInstance->meshInput.outputFileName) +
                                                      2)*sizeof(char));

                    }

                    if (filename == NULL) goto cleanup;

                    strcpy(filename, egadsInstance->meshInput.outputFileName);

                    if (egadsInstance->numSurface > 1) {
                        strcat(filename,"_Surf_");
                        strcat(filename, bodyNumber);
                    }

                    if (strcasecmp(egadsInstance->meshInput.outputFormat, "AFLR3") == 0) {

                        status = mesh_writeAFLR3(aimInfo, filename,
                                                 egadsInstance->meshInput.outputASCIIFlag,
                                                 &egadsInstance->surfaceMesh[bodyIndex],
                                                 1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "VTK") == 0) {

                        status = mesh_writeVTK(aimInfo, filename,
                                                egadsInstance->meshInput.outputASCIIFlag,
                                                &egadsInstance->surfaceMesh[bodyIndex],
                                                1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "Tecplot") == 0) {

                        status = mesh_writeTecplot(aimInfo, filename,
                                                   egadsInstance->meshInput.outputASCIIFlag,
                                                   &egadsInstance->surfaceMesh[bodyIndex],
                                                   1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "STL") == 0) {

                        status = mesh_writeSTL(aimInfo, filename,
                                               egadsInstance->meshInput.outputASCIIFlag,
                                               &egadsInstance->surfaceMesh[bodyIndex],
                                               1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "Airfoil") == 0) {

                        status = mesh_writeAirfoil(aimInfo, filename,
                                                   egadsInstance->meshInput.outputASCIIFlag,
                                                   &egadsInstance->surfaceMesh[bodyIndex],
                                                   1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "FAST") == 0) {

                        status = mesh_writeFAST(aimInfo, filename,
                                                egadsInstance->meshInput.outputASCIIFlag,
                                                &egadsInstance->surfaceMesh[bodyIndex],
                                                1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "Nastran") == 0) {

                        status = mesh_writeNASTRAN(aimInfo, filename,
                                                   egadsInstance->meshInput.outputASCIIFlag,
                                                   &egadsInstance->surfaceMesh[bodyIndex],
                                                   FreeField,
                                                   1.0);

                    } else if (strcasecmp(egadsInstance->meshInput.outputFormat, "ETO") == 0) {

                        filename = (char *) EG_reall(filename,(strlen(filename) + 5) *sizeof(char));
                        if (filename == NULL) {
                            status = EGADS_MALLOC;
                            goto cleanup;
                        }
                        strcat(filename,".eto");

                        status = EG_saveTess(egadsInstance->surfaceMesh[bodyIndex].egadsTess, filename);

                    } else {
                        printf("Unrecognized mesh format, \"%s\", the mesh will not be written out\n",
                               egadsInstance->meshInput.outputFormat);
                    }

                    AIM_FREE(filename);
                    AIM_STATUS(aimInfo, status);
                }
            }
        }
    }

    status = CAPS_SUCCESS;

cleanup:

    (void) destroy_meshStruct(&combineMesh);
    AIM_FREE(filename);

    return status;
}


int aimOutputs(/*@unused@*/ void *aimStore, /*@unused@*/ void *aimStruc,
                /*@unused@*/ int index, char **aoname, capsValue *form)
{
    /*! \page aimOutputsEgadsTess AIM Outputs
     * The following list outlines the EGADS Tessellation AIM outputs available through the AIM interface.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" egadsTessAIM/aimOutputs index = %d!\n", index);
#endif

    if (index == Done) {
        *aoname = EG_strdup("Done");
        form->type = Boolean;
        form->vals.integer = (int) false;

        /*! \page aimOutputsEgadsTess
         * - <B> Done </B> <br>
         * True if a surface mesh was created on all surfaces, False if not.
         */

    } else if (index == NumberOfElement) {
        *aoname = EG_strdup("NumberOfElement");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsEgadsTess
         * - <B> NumberOfElement </B> <br>
         * Number of elements in the surface mesh
         */

    } else if (index == NumberOfNode) {
        *aoname = EG_strdup("NumberOfNode");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsEgadsTess
         * - <B> NumberOfNode </B> <br>
         * Number of vertices in the surface mesh
         */

    } else if (index == Surface_Mesh) {
        *aoname           = AIM_NAME(Surface_Mesh);
        form->type        = Pointer;
        form->dim         = Vector;
        form->lfixed      = Change;
        form->sfixed      = Change;
        form->vals.AIMptr = NULL;
        form->nullVal     = IsNull;
        AIM_STRDUP(form->units, "meshStruct", aimStruc, status);

        /*! \page aimOutputsEgadsTess
         * - <B> Surface_Mesh </B> <br>
         * The surface mesh for a link.
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


// See if a surface mesh was generated for each body
int aimCalcOutput(void *aimStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int status = CAPS_SUCCESS;
    int surf; // Indexing
    int numElement, numNodes, nElem;
    aimStorage *egadsInstance;

#ifdef DEBUG
    printf(" egadsTessAIM/aimCalcOutput index = %d!\n", index);
#endif
    egadsInstance = (aimStorage *) aimStore;

    if (Done == index) {

        val->vals.integer = (int) false;

        // Check to see if surface meshes was generated
        for (surf = 0; surf < egadsInstance->numSurface; surf++ ) {

            if (egadsInstance->surfaceMesh[surf].numElement != 0) {

                val->vals.integer = (int) true;

            } else {

                val->vals.integer = (int) false;
                printf("No surface Tris and/or Quads were generated for surface - %d\n", surf);
                return CAPS_SUCCESS;
            }
        }

    } else if (NumberOfElement == index) {

        // Count the total number of surface elements
        numElement = 0;
        for (surf = 0; surf < egadsInstance->numSurface; surf++ ) {

            status = mesh_retrieveNumMeshElements(egadsInstance->surfaceMesh[surf].numElement,
                                                  egadsInstance->surfaceMesh[surf].element,
                                                  Triangle,
                                                  &nElem);
            AIM_STATUS(aimInfo, status);
            numElement += nElem;

            status = mesh_retrieveNumMeshElements(egadsInstance->surfaceMesh[surf].numElement,
                                                  egadsInstance->surfaceMesh[surf].element,
                                                  Quadrilateral,
                                                  &nElem);
            AIM_STATUS(aimInfo, status);
            numElement += nElem;
        }

        val->vals.integer = numElement;

    } else if (NumberOfNode == index) {

        // Count the total number of surface vertices
        numNodes = 0;
        for (surf = 0; surf < egadsInstance->numSurface; surf++ ) {
            numNodes += egadsInstance->surfaceMesh[surf].numNode;
        }

        val->vals.integer = numNodes;

    } else if (Surface_Mesh == index) {

         // Return the surface meshes
         val->nrow        = egadsInstance->numSurface;
         val->vals.AIMptr = egadsInstance->surfaceMesh;

     } else {

         status = CAPS_BADINDEX;
         AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

     }

 cleanup:

     return status;
}


void aimCleanup(void *aimStore)
{
    int        status;
    aimStorage *egadsInstance;

#ifdef DEBUG
    printf(" egadsTessAIM/aimClenup!\n");
#endif
    egadsInstance = (aimStorage *) aimStore;
    if (egadsInstance == NULL) return;

    status = destroy_aimStorage(egadsInstance);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, egadsTessAIM aimStorage cleanup!!!\n", status);

    EG_free(egadsInstance);
}
