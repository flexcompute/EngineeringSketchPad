/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             EGADS Tessellation AIM
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
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

    // quad meshing flag
    int quadMesh;

    // reference length
    double refLen;

    // Container for surface mesh
     int numMeshRef;
     aimMeshRef *meshRef;

    // Container for mesh input
    meshInputStruct meshInput;

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    mapAttrToIndexStruct meshMap;

    int numNodeTotal;
    int numElemTotal;

} aimStorage;


static int destroy_aimStorage(aimStorage *egadsInstance)
{
    int status; // Function return status
    int i;

    egadsInstance->quadMesh = 0;

    // Destroy meshInput
    status = destroy_meshInputStruct(&egadsInstance->meshInput);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, egadsTessAIM meshInput cleanup!!!\n", status);

    // Destroy surface mesh allocated arrays
    for (i = 0; i < egadsInstance->numMeshRef; i++) {
      status = aim_freeMeshRef(&egadsInstance->meshRef[i]);
      if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr4AIM surfaceMesh cleanup!!!\n", status);
    }
    AIM_FREE(egadsInstance->meshRef);
    egadsInstance->numMeshRef = 0;

    egadsInstance->numNodeTotal = 0;
    egadsInstance->numElemTotal = 0;

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

    egadsInstance->quadMesh = 0;
    egadsInstance->refLen = 0.0;

    // Container for surface meshes
    egadsInstance->numMeshRef = 0;
    egadsInstance->meshRef = NULL;

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


// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *aimStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{

    /*! \page aimInputsEgadsTess AIM Inputs
     * The following list outlines the EGADS Tessellation meshing options along with their default value available
     * through the AIM interface.
     */

    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" egadsTessAIM/aimInputs index = %d!\n", index);
#endif

    // Meshing Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotAllowed;
        defval->vals.string  = EG_strdup("egadsTess_CAPS");

        /*! \page aimInputsEgadsTess
         * - <B> Proj_Name = "egadsTess_CAPS"</B> <br>
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

        /*! \page aimInputsEgadsTess
         * \include{doc} Mesh_Format.dox
         */

    } else if (index == Mesh_Quiet_Flag) {
        *ainame               = AIM_NAME(Mesh_Quiet_Flag);
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsEgadsTess
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
         *  - "Mixed" - Quad elements for four-sided faces with TFI, triangle elements otherwise
         */

    } else if (index == Multiple_Mesh) {
        *ainame               = EG_strdup("Multiple_Mesh");
        defval->type          = String;
        AIM_STRDUP(defval->vals.string, "SingleDomain", aimInfo, status);

        /*! \page aimInputsEgadsTess
         * - <B> Multiple_Mesh = "SingleDomain"</B> <br>
         * If "SingleDomain": Generate a single surface mesh file is assuming multiple
         * bodies define a single computational domain (i.e. CFD)<br>
         * <br>
         * If "MultiFile": Generate a surface mesh file for each body.<br>
         * <br>
         * If "MultiDomain": Generate a single mesh file containing multiple surface meshes for each body.<br>
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

    } else {
        AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
        status = CAPS_BADINDEX;
        goto cleanup;
    }

    AIM_NOTNULL(*ainame, aimInfo, status);
    status = CAPS_SUCCESS;

cleanup:
    return status;
}


// ********************** AIM Function Break *****************************
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
    int status; // Function return status
    int i, j;

    int MultiMesh = -1;
    char bodyNumberFile[128];
    char aimFile[PATH_MAX];

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    // Global settings
    int minEdgePoint = -1, maxEdgePoint = -1;
    double meshLenFac = 0, capsMeshLength = 0;
    int bodyIndex;

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;
    const char *MeshElements = NULL;

    aimStorage *egadsInstance;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if ((numBody <= 0) || (bodies == NULL)) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }
    AIM_NOTNULL(aimInputs, aimInfo, status);

    egadsInstance = (aimStorage *) instStore;

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(egadsInstance);
    AIM_STATUS(aimInfo, status);

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


    // Get Tessellation parameters -Tess_Params
    egadsInstance->meshInput.paramTess[0] = aimInputs[Tess_Params-1].vals.reals[0]; // Gets multiplied by bounding box size
    egadsInstance->meshInput.paramTess[1] = aimInputs[Tess_Params-1].vals.reals[1]; // Gets multiplied by bounding box size
    egadsInstance->meshInput.paramTess[2] = aimInputs[Tess_Params-1].vals.reals[2];

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

         if ( strncasecmp(MeshElements,"Tri",3)   == 0 ) { egadsInstance->quadMesh = 0; }
    else if ( strncasecmp(MeshElements,"Quad",4)  == 0 ) { egadsInstance->quadMesh = 1; }
    else if ( strncasecmp(MeshElements,"Mixed",3) == 0 ) { egadsInstance->quadMesh = 2; }
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
            AIM_STATUS(aimInfo, status);
        }
    }

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

    egadsInstance->refLen = meshLenFac*capsMeshLength;

    // Modify the EGADS body tessellation based on given inputs
/*@-nullpass@*/
    status =  mesh_modifyBodyTess(aimInfo,
                                  numMeshProp,
                                  meshProp,
                                  minEdgePoint,
                                  maxEdgePoint,
                                  egadsInstance->quadMesh,
                                  &egadsInstance->refLen,
                                  egadsInstance->meshInput.paramTess,
                                  egadsInstance->meshMap,
                                  numBody,
                                  bodies);
/*@+nullpass@*/
    AIM_STATUS(aimInfo, status);


    if (strcasecmp(aimInputs[Multiple_Mesh-1].vals.string, "SingleDomain") == 0) {
      MultiMesh = 0;
    } else if (strcasecmp(aimInputs[Multiple_Mesh-1].vals.string, "MultiFile") == 0) {
      MultiMesh = 1;
    } else if (strcasecmp(aimInputs[Multiple_Mesh-1].vals.string, "MultiDomain") == 0) {
      MultiMesh = 2;
    } else {
      AIM_ERROR(aimInfo, "Multiple_Mesh = '%s' must be 'SingleDomain', 'MultiFile', or 'MultiDomain'", aimInputs[Multiple_Mesh-1].vals.string);
      status = CAPS_BADVALUE;
      goto cleanup;
    }


    if (MultiMesh == 0 || MultiMesh == 2) {

      AIM_ALLOC(egadsInstance->meshRef, 1, aimMeshRef, aimInfo, status);
      egadsInstance->numMeshRef = 1;

      status = aim_initMeshRef(egadsInstance->meshRef, aimSurfaceMesh);
      AIM_STATUS(aimInfo, status);

      AIM_ALLOC(egadsInstance->meshRef[0].maps, numBody, aimMeshTessMap, aimInfo, status);
      egadsInstance->meshRef[0].nmap = numBody;

       for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {
         egadsInstance->meshRef[0].maps[bodyIndex].tess = NULL;
         egadsInstance->meshRef[0].maps[bodyIndex].map = NULL;
       }

       // set the filename without extensions where the grid is written for solvers
       status = aim_file(aimInfo, aimInputs[Proj_Name-1].vals.string, aimFile);
       AIM_STATUS(aimInfo, status);
       AIM_STRDUP(egadsInstance->meshRef[0].fileName, aimFile, aimInfo, status);

    } else  if (MultiMesh == 1) {

      AIM_ALLOC(egadsInstance->meshRef, numBody, aimMeshRef, aimInfo, status);
      egadsInstance->numMeshRef = numBody;

      for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {
        status = aim_initMeshRef(&egadsInstance->meshRef[bodyIndex], aimSurfaceMesh);
        AIM_STATUS(aimInfo, status);

        AIM_ALLOC(egadsInstance->meshRef[bodyIndex].maps, 1, aimMeshTessMap, aimInfo, status);
        egadsInstance->meshRef[bodyIndex].nmap = 1;

        egadsInstance->meshRef[bodyIndex].maps[0].tess = NULL;
        egadsInstance->meshRef[bodyIndex].maps[0].map = NULL;

        // set the filename without extensions where the grid is written for solvers
        snprintf(bodyNumberFile, 128, "%s_%d", aimInputs[Proj_Name-1].vals.string, bodyIndex);
        status = aim_file(aimInfo, bodyNumberFile, aimFile);
        AIM_STATUS(aimInfo, status);
        AIM_STRDUP(egadsInstance->meshRef[bodyIndex].fileName, aimFile, aimInfo, status);
      }
    }

    for (i = 0; i < egadsInstance->numMeshRef; i++) {

      AIM_ALLOC(egadsInstance->meshRef[i].bnds, egadsInstance->groupMap.numAttribute, aimMeshBnd, aimInfo, status);
      egadsInstance->meshRef[i].nbnd = egadsInstance->groupMap.numAttribute;
      for (j = 0; j < egadsInstance->meshRef[i].nbnd; j++) {
        status = aim_initMeshBnd(egadsInstance->meshRef[i].bnds + j);
        AIM_STATUS(aimInfo, status);
      }

      for (j = 0; j < egadsInstance->meshRef[i].nbnd; j++) {
        AIM_STRDUP(egadsInstance->meshRef[i].bnds[j].groupName, egadsInstance->groupMap.attributeName[j], aimInfo, status);
        egadsInstance->meshRef[i].bnds[j].ID = egadsInstance->groupMap.attributeIndex[j];
      }
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


// ********************** AIM Function Break *****************************
int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{
    int status; // Status return

    int i, bodyIndex; // Indexing

    const aimStorage *egadsInstance;

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects
    ego etess = NULL;

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
    AIM_NOTNULL(aimInputs, aimInfo, status);

    egadsInstance = (const aimStorage *) instStore;

    // remove previous meshes
    for (i = 0; i < egadsInstance->numMeshRef; i++) {
      status = aim_deleteMeshes(aimInfo, &egadsInstance->meshRef[i]);
      AIM_STATUS(aimInfo, status);
    }

    // Run egadsTess for each body
    for (bodyIndex = 0 ; bodyIndex < numBody; bodyIndex++) {

        if (aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false)
          printf("Getting surface mesh for body %d (of %d)\n", bodyIndex+1, numBody);

        status = mesh_surfaceMeshEGADSBody(aimInfo,
                                           bodies[bodyIndex],
                                           egadsInstance->refLen,
                                           egadsInstance->meshInput.paramTess,
                                           egadsInstance->quadMesh,
                                           &etess);
        AIM_STATUS(aimInfo, status, "Problem during surface meshing of body %d", bodyIndex+1);
        AIM_NOTNULL(etess, aimInfo, status);

        // set the file name to write the egads file
        snprintf(bodyNumber, 42, EGADSTESSFILE, bodyIndex);
        status = aim_rmFile(aimInfo, bodyNumber);
        AIM_STATUS(aimInfo, status);

        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        status = EG_saveTess(etess, aimFile);
        AIM_STATUS(aimInfo, status);

        EG_deleteObject(etess);
    }

#ifdef DEBUG
    {
      ego *bodyCopy=NULL;
      AIM_ALLOC(bodyCopy, numBody, ego, aimInfo, status);
      for (bodyIndex = 0 ; bodyIndex < numBody; bodyIndex++) {
        EG_copyObject(bodies[bodyIndex], NULL, &bodyCopy[bodyIndex]);
      }
      ego context, model;
      EG_getContext(bodyCopy[0], &context);
      EG_makeTopology(context, NULL, MODEL, 0, NULL, numBody, bodyCopy, NULL, &model);
      remove("egadsTess_debug.egads");
      EG_saveModel(model, "egadsTess_debug.egads");
      EG_deleteObject(model);
      AIM_FREE(bodyCopy);
    }
#endif

    status = CAPS_SUCCESS;

cleanup:

    return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ const void *aimStore, /*@unused@*/ void *aimStruc, int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis(/*@unused@*/ void *aimStore, /*@unused@*/ void *aimInfo,
                    /*@unused@*/ int restart,    /*@unused@*/ capsValue *aimInputs)
{

    int status = CAPS_SUCCESS;
    int bodyIndex;
    int MultiMesh = -1;

    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    int state, nglobal, i, iglobal=1;
    int iface, numFace, numTri, numQuad;
    ego body, tess;

    // EGADS function return variables
    int           plen = 0, tlen = 0;
    const int    *tris = NULL, *tric = NULL, *ptype = NULL, *pindex = NULL;
    const double *points = NULL, *uv = NULL;

    const int *tessFaceQuadMap = NULL;
    int alen, atype;
    const double *reals = NULL;
    const char *string = NULL;

    const char *intents;
    char bodyNumber[42];
    char aimFile[PATH_MAX];
    aimMesh    mesh;

    aimStorage *egadsInstance;
    egadsInstance = (aimStorage *) aimStore;

    AIM_NOTNULL(aimInputs, aimInfo, status);

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodies, aimInfo, status);

    if (strcasecmp(aimInputs[Multiple_Mesh-1].vals.string, "SingleDomain") == 0) {
      MultiMesh = 0;
    } else if (strcasecmp(aimInputs[Multiple_Mesh-1].vals.string, "MultiFile") == 0) {
      MultiMesh = 1;
    } else if (strcasecmp(aimInputs[Multiple_Mesh-1].vals.string, "MultiDomain") == 0) {
      MultiMesh = 2;
    } else {
      AIM_ERROR(aimInfo, "Developer error! Unknown Multiple_Mesh %s", aimInputs[Multiple_Mesh-1].vals.string);
      status = CAPS_BADVALUE;
      goto cleanup;
    }


    // Read mesh for each body
    for (bodyIndex = 0 ; bodyIndex < numBody; bodyIndex++) {

      // set the file name to read the egads file
      snprintf(bodyNumber, 42, EGADSTESSFILE, bodyIndex);
      status = aim_file(aimInfo, bodyNumber, aimFile);
      AIM_STATUS(aimInfo, status);

      status = EG_loadTess(bodies[bodyIndex], aimFile, &tess);
      AIM_STATUS(aimInfo, status);

      status = aim_newTess(aimInfo, tess);
      AIM_STATUS(aimInfo, status);

      status = EG_statusTessBody(tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      if (MultiMesh == 0 || MultiMesh == 2) {

        egadsInstance->meshRef[0].maps[bodyIndex].tess = tess;

        AIM_ALLOC(egadsInstance->meshRef[0].maps[bodyIndex].map, nglobal, int, aimInfo, status);
        for (i = 0; i < nglobal; i++) egadsInstance->meshRef[0].maps[bodyIndex].map[i] = iglobal++;

      } else if (MultiMesh == 1) {

        egadsInstance->meshRef[bodyIndex].maps[0].tess = tess;

        AIM_ALLOC(egadsInstance->meshRef[bodyIndex].maps[0].map, nglobal, int, aimInfo, status);
        for (i = 0; i < nglobal; i++) egadsInstance->meshRef[bodyIndex].maps[0].map[i] = iglobal++;

      }

      // Get faces, edges, and nodes so we can check for attributes on them
      status = EG_getBodyTopos(bodies[bodyIndex], NULL, FACE, &numFace, NULL);
      AIM_STATUS(aimInfo, status);

      numQuad = numTri = 0;
      for (iface = 1; iface <= numFace; iface++) {
          status = EG_getTessFace(tess, iface, &plen, &points, &uv, &ptype, &pindex,
                                  &tlen, &tris, &tric);
          AIM_STATUS(aimInfo, status);

          numTri += tlen;
      }

      // check if the tessellation has a mixture of quad and tri
      status = EG_attributeRet(tess, ".mixed",
                               &atype, &alen, &tessFaceQuadMap, &reals, &string);
      AIM_NOTFOUND(aimInfo, status);

      // Do we have split quads?
      if (tessFaceQuadMap != NULL) {
          for (iface = 0; iface < numFace; iface++)
            numQuad += tessFaceQuadMap[iface];

          // subtract off the split quads from the total tri count
          numTri -= 2*numQuad;
      }

      if (restart == 0 &&
          aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false) {

        printf("Body %d (of %d)\n", bodyIndex+1, numBody);
        printf("Number of nodes    = %d\n", nglobal);
        printf("Number of elements = %d\n", numTri+numQuad);
        printf("Number of triangle elements      = %d\n", numTri);
        printf("Number of quadrilateral elements = %d\n", numQuad);
      }

      egadsInstance->numNodeTotal += nglobal;
      egadsInstance->numElemTotal += numTri+numQuad;
    }

    if (restart == 0 &&
        aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false) {
        printf("----------------------------\n");
        printf("Total number of nodes    = %d\n", egadsInstance->numNodeTotal);
        printf("Total number of elements = %d\n", egadsInstance->numElemTotal);
    }

    for (i = 0; i < egadsInstance->numMeshRef; i++) {
      status = aim_queryMeshes( aimInfo, Mesh_Format, ANALYSISIN, &egadsInstance->meshRef[i] );
      if (status > 0) {
  /*@-immediatetrans@*/
        mesh.meshData = NULL;
        mesh.meshRef = &egadsInstance->meshRef[i];
  /*@+immediatetrans@*/

        status = mesh_surfaceMeshData(aimInfo, &egadsInstance->groupMap, &mesh);
        AIM_STATUS(aimInfo, status);

        status = aim_writeMeshes(aimInfo, Mesh_Format, ANALYSISIN, &mesh);
        AIM_STATUS(aimInfo, status);

        status = aim_freeMeshData(mesh.meshData);
        AIM_STATUS(aimInfo, status);
        AIM_FREE(mesh.meshData);
      }
      else
        AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;

cleanup:
    return status;
}


// ********************** AIM Function Break *****************************
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
        form->type        = PointerMesh;
        form->dim         = Vector;
        form->lfixed      = Change;
        form->sfixed      = Fixed;
        form->vals.AIMptr = NULL;
        form->nullVal     = IsNull;

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


// ********************** AIM Function Break *****************************
int aimCalcOutput(void *aimStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int status = CAPS_SUCCESS;
    int i;
    aimStorage *egadsInstance;
    aimMesh    mesh;

#ifdef DEBUG
    printf(" egadsTessAIM/aimCalcOutput index = %d!\n", index);
#endif
    egadsInstance = (aimStorage *) aimStore;

    if (Done == index) {

      if (egadsInstance->numNodeTotal > 0 && egadsInstance->numElemTotal > 0)
        val->vals.integer = (int) true;
      else
        val->vals.integer = (int) false;

    } else if (NumberOfElement == index) {

      val->vals.integer = egadsInstance->numElemTotal;

    } else if (NumberOfNode == index) {

      val->vals.integer = egadsInstance->numNodeTotal;

    } else if (Surface_Mesh == index) {

      for (i = 0; i < egadsInstance->numMeshRef; i++) {
        status = aim_queryMeshes( aimInfo, Surface_Mesh, ANALYSISOUT, &egadsInstance->meshRef[i] );
        if (status > 0) {
/*@-immediatetrans@*/
          mesh.meshData = NULL;
          mesh.meshRef = &egadsInstance->meshRef[i];
/*@+immediatetrans@*/

          status = mesh_surfaceMeshData(aimInfo, &egadsInstance->groupMap, &mesh);
          AIM_STATUS(aimInfo, status);

          status = aim_writeMeshes(aimInfo, Surface_Mesh, ANALYSISOUT, &mesh);
          AIM_STATUS(aimInfo, status);

          status = aim_freeMeshData(mesh.meshData);
          AIM_STATUS(aimInfo, status);
          AIM_FREE(mesh.meshData);
        }
        else
          AIM_STATUS(aimInfo, status);
      }

      // Return the surface meshes
      /*@-immediatetrans@*/
      val->nrow        = egadsInstance->numMeshRef;
      val->vals.AIMptr = egadsInstance->meshRef;
      /*@+immediatetrans@*/

    } else {

      status = CAPS_BADINDEX;
      AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

    }

 cleanup:

     return status;
}


// ********************** AIM Function Break *****************************
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
