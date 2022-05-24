/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             TetGen AIM
 *
 *      Written by Dr. Ryan Durscher AFRL/RQVC
 *
 *     This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */


/*! \mainpage Introduction
 *
 * \section overviewTetGen TetGen AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact with
 * the open-source volume mesh generator, TetGen \cite Hang2015. TetGen is capable of generating exact constrained
 * Delaunay tetrahedralizations, boundary conforming Delaunay meshes, and Voronoi partitions.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsTetGen and \ref aimOutputsTetGen, respectively.
 *
 * Current issues include:
 *  - The holes or seed points provided to TetGen by creating an taking the cnetroid of a tetrahedron from an 'empty' mesh.
 *    This is guaranteed to work with solid bodies, but sheet bodies with multiple segregated regions where some
 *    regions are holes require manual seed points to indicate the hole.
 *  - (<b>Important</b>) If Tetgen is allowed to added Steiner points (see "Preserve_Surf_Mesh" in \ref aimInputsTetGen)
 *  discrete data transfer will <b>NOT</b> be possible.
 *
 * \section tetgenInterface TetGen Interface
 * In order to use TetGen, CAPS will automatically build the TetGen source code supplied in "library" mode. The directory
 * in which the source code exists is set in the ESP configuration script. The C++ API is interfaced within the AIM
 * through an interface function that takes the body tessellation and transfers the data to a "tetgenio" object
 * in PLCs format (Piecewise Linear Complexes). After volume meshing is complete the mesh can be output in various
 *  mesh formats (see \ref aimInputsTetGen for additional details).
 *
 *  \section clearanceTetgen Clearance Statement
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

#include "tetgen_Interface.hpp" // Bring in TetGen 'interface' function

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

//#define DEBUG


enum aimInputs
{
  Proj_Name = 1,                 /* index is 1-based */
  Preserve_Surf_Mesh,
  Mesh_Verbose_Flag,
  Mesh_Quiet_Flag,
  Quality_Rad_Edge,
  Quality_Angle,
  Mesh_Format,
  Mesh_ASCII_Flag,
  Mesh_Gen_Input_String,
  Ignore_Surface_Mesh_Extraction,
  Mesh_Tolerance,
  Multiple_Mesh,
  Regions,
  Holes,
  Surface_Mesh,
  NUMINPUT = Surface_Mesh        /* Total number of inputs */
};

enum aimOutputs
{
  NumberOfElement = 1,         /* index is 1-based */
  NumberOfNode,
  Volume_Mesh,
  NUMOUT = Volume_Mesh         /* Total number of outputs */
};

#define NODATATRANSFER "noDataTransfer.%d"


typedef struct {

    // Container for mesh input
    meshInputStruct meshInput;

    // Attribute to index map
    mapAttrToIndexStruct attrMap;

    // Mesh references for link
    int numMeshRef;
    aimMeshRef *meshRef;

} aimStorage;



static int destroy_aimStorage(aimStorage *tetgenInstance)
{

    int i; // Indexing

    int status; // Function return status

    // Destroy meshInput
    status = destroy_meshInputStruct(&tetgenInstance->meshInput);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, tetgenAIM  meshInput cleanup!!!\n", status);

    // Destroy attribute to index map
    status = destroy_mapAttrToIndexStruct(&tetgenInstance->attrMap);
    if (status != CAPS_SUCCESS)
      printf("Status = %d, tetgenAIM  attrMap cleanup!!!\n", status);

    // Free the meshRef
    for (i = 0; i < tetgenInstance->numMeshRef; i++)
      aim_freeMeshRef(&tetgenInstance->meshRef[i]);
    AIM_FREE(tetgenInstance->meshRef);

    return CAPS_SUCCESS;
}



/* ********************** Exposed AIM Functions ***************************** */

extern "C" int
aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
              /*@unused@*/ void **instStore, /*@unused@*/ int *major,
              /*@unused@*/ int *minor, int *nIn, int *nOut,
              int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int status = CAPS_SUCCESS; // Function status return

    aimStorage *tetgenInstance = NULL;

#ifdef DEBUG
    printf("\n tetgenAIM/aimInitialize   instance = %d!\n", inst);
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

    // Allocate tetgenInstance
    AIM_ALLOC(tetgenInstance, 1, aimStorage, aimInfo, status);
    *instStore = tetgenInstance;

    // Set initial values for tetgenInstance //

    // Mesh reference passed to solver
    tetgenInstance->numMeshRef = 0;
    tetgenInstance->meshRef = NULL;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&tetgenInstance->attrMap);
    AIM_STATUS(aimInfo, status);

    // Container for mesh input
    status = initiate_meshInputStruct(&tetgenInstance->meshInput);
    AIM_STATUS(aimInfo, status);


cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);
    return CAPS_SUCCESS;
}


extern "C" int
aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, int index,
          char **ainame, capsValue *defval)
{
    /*! \page aimInputsTetGen AIM Inputs
     * The following list outlines the TetGen meshing options along with their default value available
     * through the AIM interface.
     */
    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" tetgenAIM/aimInputs  index = %d!\n", index);
#endif

    // Meshing Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name"); // If NULL a volume grid won't be written by the AIM
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->vals.string  = NULL;
        //defval->vals.string  = EG_strdup("CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsTetGen
         * - <B> Proj_Name = NULL</B> <br>
         * This corresponds to the output name of the mesh. If left NULL, the mesh is not written to a file.
         */

    } else if (index == Preserve_Surf_Mesh) {
        *ainame               = EG_strdup("Preserve_Surf_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = true;

        /*! \page aimInputsTetGen
         * - <B> Preserve_Surf_Mesh = True</B> <br>
         * Tells TetGen to preserve the surface mesh provided (i.e. do not add Steiner points on the surface). Discrete data transfer
         * will <b>NOT</b> be possible if Steiner points are added.
         */
    } else if (index == Mesh_Verbose_Flag) {
        *ainame               = EG_strdup("Mesh_Verbose_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsTetGen
         * - <B> Mesh_Verbose_Flag = False</B> <br>
         * Verbose output from TetGen.
         */
    } else if (index == Mesh_Quiet_Flag) {
        *ainame               = EG_strdup("Mesh_Quiet_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsTetGen
         * - <B> Mesh_Quiet_Flag = False</B> <br>
         * Complete suppression of all TetGen output messages (not including errors).
         */
    } else if (index == Quality_Rad_Edge) {
        *ainame               = EG_strdup("Quality_Rad_Edge"); // TetGen specific parameters
        defval->type          = Double;
        defval->vals.real     = 1.5;

        /*! \page aimInputsTetGen
         * - <B> Quality_Rad_Edge = 1.5</B> <br>
         * TetGen maximum radius-edge ratio.
         */
    } else if (index == Quality_Angle) {
        *ainame               = EG_strdup("Quality_Angle"); // TetGen specific parameters
        defval->type          = Double;
        defval->vals.real     = 0.0;

        /*! \page aimInputsTetGen
         * - <B> Quality_Angle = 0.0</B> <br>
         * TetGen minimum dihedral angle (in degrees).
         */
    } else if (index == Mesh_Format) {
        *ainame               = EG_strdup("Mesh_Format");
        defval->type          = String;
        defval->vals.string   = EG_strdup("AFLR3"); // TECPLOT, VTK, AFLR3, SU2, NASTRAN
        defval->lfixed        = Change;

        /*! \page aimInputsTetGen
         * - <B> Mesh_Format = "AFLR3"</B> <br>
         * Mesh output format. Available format names  include: "AFLR3", "TECPLOT", "SU2", "VTK", and "NASTRAN".
         */
    } else if (index == Mesh_ASCII_Flag) {
        *ainame               = EG_strdup("Mesh_ASCII_Flag");
        defval->type          = Boolean;
        defval->vals.integer  = true;

        /*! \page aimInputsTetGen
         * - <B> Mesh_ASCII_Flag = True</B> <br>
         * Output mesh in ASCII format, otherwise write a binary file if applicable.
         */
    } else if (index == Mesh_Gen_Input_String) {
        *ainame               = EG_strdup("Mesh_Gen_Input_String");
        defval->type          = String;
        defval->nullVal       = IsNull;
        defval->vals.string   = NULL;
        defval->lfixed       = Change;

        /*! \page aimInputsTetGen
         * - <B> Mesh_Gen_Input_String = NULL</B> <br>
         * Meshing program command line string (as if called in bash mode). Use this to specify more
         * complicated options/use features of the mesher not currently exposed through other AIM input
         * variables. Note that this is the exact string that will be provided to the volume mesher; no
         * modifications will be made. If left NULL an input string will be created based on default values
         * of the relevant AIM input variables. See \ref tetgenCommandLineInput for options to include
         * in the input string.
         */
    } else if (index == Ignore_Surface_Mesh_Extraction) { // Deprecated
        *ainame               = EG_strdup("Ignore_Surface_Mesh_Extraction");
        defval->type          = Boolean;
        defval->vals.integer  = true;

//        /*! \page aimInputsTetGen
//         * - <B> Ignore_Surface_Mesh_Extraction = True</B> <br>
//         * If TetGen doesn't preserve the surface mesh provided (i.e. Steiner points are added) a simple search algorithm may be
//         * used to reconstruct a separate (not dependent on the volume mesh node numbering) representation of the surface mesh. In
//         * general, this has little use and can add a significant computational penalty. The default value of "True" is recommended.
//         */
    } else if (index == Mesh_Tolerance) {
        *ainame               = EG_strdup("Mesh_Tolerance"); // TetGen specific parameters
        defval->type          = Double;
        defval->vals.real     = 1E-16;

        /*! \page aimInputsTetGen
         * - <B> Mesh_Tolerance = 1E-16</B> <br>
         * Sets the tolerance for coplanar test in TetGen.
         */
    } else if (index == Multiple_Mesh) {
        *ainame               = EG_strdup("Multiple_Mesh");
        defval->type          = Boolean;
        defval->vals.integer  = (int) false;

        /*! \page aimInputsTetGen
         * - <B> Multiple_Mesh = False</B> <br>
         * If set to True a volume will be generated for each body. When set to False (default value) only a single volume
         * mesh will be created.
         */

    } else if (index == Regions) {
        *ainame              = EG_strdup("Regions");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsTetGen
         * - <B>Regions = NULL</B><br>
         * If this input is set, the volume mesh will be divided into regions,
         * each bounded by surface mesh and identified by an interior `seed`
         * point.  If a seed appears in the `Regions` input, then its region
         * will be meshed and the markers for all cells in the region will be
         * set to the region's `id`.  The markers for cells that fall outside
         * of any user-defined region or hole will be numbered automatically.
         * The input is a vector of tuples.  The tuple keys are ignored and the
         * tuple values are dictionaries; each requires an integer `id` entry
         * and a 3-vector `seed` point.  For example, from within a pyCAPS
         * script,
         *
         *  \code{.py}
         *  tetgen.input.Regions = {
         *    'A': { 'id': 10, 'seed': [0, 0,  1] },
         *    'B': { 'id': 20, 'seed': [0, 0, -1] }
         *  }
         *  \endcode
         *
         * Automatic hole detection will be disabled if one or both of the
         * `Regions` and `Holes` inputs is not NULL.
         */
    } else if (index == Holes) {
      *ainame              = EG_strdup("Holes");
      defval->type         = Tuple;
      defval->nullVal      = IsNull;
      defval->dim          = Vector;
      defval->lfixed       = Change;
      defval->vals.tuple   = NULL;

      /*! \page aimInputsTetGen
       * - <B>Holes = NULL</B><br>
       * If this input is set, the volume mesh will be divided into regions,
       * each bounded by surface mesh and identified by an interior `seed`
       * point.  If a seed appears in the `Holes` input, then its region will
       * not be meshed.  The input is a vector of tuples.  The tuple keys are
       * ignored and the tuple values are dictionaries; each requires a
       * 3-vector `seed` point.  For example, from within a pyCAPS script,
       *
       *  \code{.py}
       *  tetgen.input.Holes = {
       *    'A': { 'seed': [ 1, 0, 0] },
       *    'B': { 'seed': [-1, 0, 0] }
       *  }
       *  \endcode
       *
       * Automatic hole detection will be disabled if one or both of the
       * `Regions` and `Holes` inputs is not NULL.
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

        /*! \page aimInputsTetGen
         * - <B>Surface_Mesh = NULL</B> <br>
         * A Surface_Mesh link.
         */
    }

    AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
    return status;
}


extern "C" int
aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{
    int i, j, elem, ibody; // Indexing

    int status; // Return status

    aimStorage *tetgenInstance;

    // Incoming bodies
    const char *intents;
    ego *bodies = NULL;
    int numBody = 0;

    // Container for surface mesh
    int numSurfaceMesh = 0;
    meshStruct *surfaceMesh = NULL;

    // Container for volume mesh
    int numVolumeMesh = 0;
    meshStruct *volumeMesh = NULL;

    // Meshing related variables
    meshElementStruct *element;

    double  box[6], boxMax[6] = {0,0,0,0,0,0};
    int     bodyBoundingBox = 0;

    int numElementCheck;

    // File ouput
    char *filename=NULL;
    char bodyNumberFile[42];
    char aimFile[PATH_MAX];
    FILE *fp;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" tetgenAIM/aimPreAnalysis  numBody = %d!\n", numBody);
#endif

    if (numBody <= 0 || bodies == NULL) {
#ifdef DEBUG
        printf(" tetgenAIM/aimPreAnalysis No Bodies!\n");
#endif
        return CAPS_SOURCEERR;
    }
    tetgenInstance = (aimStorage *) instStore;

    // remove previous meshes
    for (ibody = 0; ibody < tetgenInstance->numMeshRef; ibody++) {
      status = aim_deleteMeshes(aimInfo, &tetgenInstance->meshRef[ibody]);
      AIM_STATUS(aimInfo, status);
    }

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(tetgenInstance);
    AIM_STATUS(aimInfo, status, "tetgenAIM aimStorage cleanup!!!");

    // Get capsGroup name and index mapping
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            2, // Only search down to the edge level of the EGADS body
                                            &tetgenInstance->attrMap);
    AIM_STATUS(aimInfo, status);

    // Get surface mesh
    if (aimInputs[Surface_Mesh-1].nullVal == IsNull) {
        AIM_ANALYSISIN_ERROR(aimInfo, Surface_Mesh, "'Surface_Mesh' input must be linked to an output 'Surface_Mesh'");
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    // Get mesh
    numSurfaceMesh = aimInputs[Surface_Mesh-1].length;
    surfaceMesh    = (meshStruct *)aimInputs[Surface_Mesh-1].vals.AIMptr;

    if (numSurfaceMesh != numBody) {
        AIM_ANALYSISIN_ERROR(aimInfo, Surface_Mesh, "Number of linked surface meshes (%d) does not match the number of bodies (%d)\n",
                             numSurfaceMesh, numBody);
        return CAPS_SOURCEERR;
    }

    // Create/setup volume meshes
    if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {

        AIM_ALLOC(volumeMesh, numBody, meshStruct, aimInfo, status);
        numVolumeMesh = numBody;

        for (ibody = 0; ibody < numVolumeMesh; ibody++) {
            status = initiate_meshStruct(&volumeMesh[ibody]);
            AIM_STATUS(aimInfo, status);

            // Set reference mesh - One surface per body
            volumeMesh[ibody].numReferenceMesh = 1;
            AIM_ALLOC(volumeMesh[ibody].referenceMesh, volumeMesh[ibody].numReferenceMesh, meshStruct, aimInfo, status);

            volumeMesh[ibody].referenceMesh[0] = surfaceMesh[ibody];
            printf("Tetgen MultiMesh TopoIndex = %d\n",
                   volumeMesh[0].referenceMesh[0].element[0].topoIndex);
        }

    } else {

        // Determine which body is the bounding body based on size
        bodyBoundingBox = 0;
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
                    bodyBoundingBox = ibody;

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
                    bodyBoundingBox = ibody;
                }
            }
        }

        // Swap the surface element orientation so normals point out of the computational domain
        for (ibody = 0; ibody < numBody; ibody++) {
            if (ibody != bodyBoundingBox) {
                // Swap two indices to reverse the normal vector of all elements on internal bodies
                // so they point out of the domain
                for (elem = 0; elem < surfaceMesh[ibody].numElement; elem++) {

                    element = surfaceMesh[ibody].element + elem;

                    // This should be valid for both Triangles and Quadrilaterals
                    i = element->connectivity[2];
                    element->connectivity[2] = element->connectivity[0];
                    element->connectivity[0] = i;
                }
            }
        }

        numVolumeMesh = 1;
        AIM_ALLOC(volumeMesh, numVolumeMesh, meshStruct, aimInfo, status);

        status = initiate_meshStruct(&volumeMesh[0]);
        AIM_STATUS(aimInfo, status);

        // Combine mesh - temporary store the combined mesh in the volume mesh
        status = mesh_combineMeshStruct(numSurfaceMesh,
                                        surfaceMesh,
                                        &volumeMesh[0]);
        AIM_STATUS(aimInfo, status);

        // Set reference meshes - All surfaces
        volumeMesh[0].numReferenceMesh = numSurfaceMesh;
        AIM_ALLOC(volumeMesh[0].referenceMesh, volumeMesh[0].numReferenceMesh, meshStruct, aimInfo, status);

        for (ibody = 0; ibody < numSurfaceMesh; ibody++) {
            volumeMesh[0].referenceMesh[ibody] = surfaceMesh[ibody];
        }

         // Report surface mesh
        printf("Number of surface nodes - %d\n",
               volumeMesh[0].numNode);
        printf("Number of surface elements - %d\n",
               volumeMesh[0].numElement);
    }

    // Setup meshing input structure - mesher specific input gets set below before entering interface
    tetgenInstance->meshInput.preserveSurfMesh = aimInputs[Preserve_Surf_Mesh-1].vals.integer;
    tetgenInstance->meshInput.quiet            = aimInputs[Mesh_Quiet_Flag-1].vals.integer;
    tetgenInstance->meshInput.outputASCIIFlag  = aimInputs[Mesh_ASCII_Flag-1].vals.integer;

    AIM_STRDUP(tetgenInstance->meshInput.outputFormat, aimInputs[Mesh_Format-1].vals.string, aimInfo, status);

    if (aimInputs[Proj_Name-1].nullVal != IsNull) {
        AIM_STRDUP(tetgenInstance->meshInput.outputFileName, aimInputs[Proj_Name-1].vals.string, aimInfo, status);
    }

    // Set tetgen specific mesh inputs
    tetgenInstance->meshInput.tetgenInput.meshQuality_rad_edge = aimInputs[Quality_Rad_Edge-1].vals.real;
    tetgenInstance->meshInput.tetgenInput.meshQuality_angle    = aimInputs[Quality_Angle-1].vals.real;
    tetgenInstance->meshInput.tetgenInput.verbose              = aimInputs[Mesh_Verbose_Flag-1].vals.integer;
    tetgenInstance->meshInput.tetgenInput.ignoreSurfaceExtract = aimInputs[Ignore_Surface_Mesh_Extraction-1].vals.integer;
    tetgenInstance->meshInput.tetgenInput.meshTolerance        = aimInputs[Mesh_Tolerance-1].vals.real;

    if (aimInputs[Regions-1].nullVal != IsNull) {
      populate_regions(&tetgenInstance->meshInput.tetgenInput.regions,
                       aimInputs[Regions-1].length,
                       aimInputs[Regions-1].vals.tuple);
    }
    if (aimInputs[Holes-1].nullVal != IsNull) {
      populate_holes(&tetgenInstance->meshInput.tetgenInput.holes,
                     aimInputs[Holes-1].length,
                     aimInputs[Holes-1].vals.tuple);
    }

    if (aimInputs[Mesh_Gen_Input_String-1].nullVal != IsNull) {

        AIM_STRDUP(tetgenInstance->meshInput.tetgenInput.meshInputString,
                   aimInputs[Mesh_Gen_Input_String-1].vals.string, aimInfo, status);
    }

    status = populate_bndCondStruct_from_mapAttrToIndexStruct(&tetgenInstance->attrMap,
                                                              &tetgenInstance->meshInput.bndConds);
    AIM_STATUS(aimInfo, status);

    // Call tetgen volume mesh interface
    for (ibody = 0; ibody < numVolumeMesh; ibody++) {

        snprintf(bodyNumberFile, 42, "tetgen_%d", ibody);
        status = aim_file(aimInfo, bodyNumberFile, aimFile);
        AIM_STATUS(aimInfo, status);

        // Call tetgen volume mesh interface for each body
        if (numVolumeMesh > 1) {
            printf("Getting volume mesh for body %d (of %d)\n", ibody+1, numBody);

            status = tetgen_VolumeMesh(aimInfo,
                                       tetgenInstance->meshInput,
                                       aimFile,
                                      &volumeMesh[ibody].referenceMesh[0],
                                      &volumeMesh[ibody]);
        } else {
            printf("Getting volume mesh\n");

            status = tetgen_VolumeMesh(aimInfo,
                                       tetgenInstance->meshInput,
                                       aimFile,
                                      &volumeMesh[ibody],
                                      &volumeMesh[ibody]);
        }

        if (status != CAPS_SUCCESS) {
            if (numVolumeMesh > 1) {
              AIM_ERROR(aimInfo, "TetGen volume mesh failed on body - %d!!!!", ibody+1);
            } else {
              AIM_ERROR(aimInfo, "TetGen volume mesh failed!!!!");
            }
            goto cleanup;
        }
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
            fclose(fp); fp = NULL;

            printf("Volume mesher added surface elements - data transfer will NOT be possible.\n");
        }
    }

    // If filename is not NULL write the mesh
    if (tetgenInstance->meshInput.outputFileName != NULL) {

        for (ibody = 0; ibody < numVolumeMesh; ibody++) {

            if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {
                sprintf(bodyNumberFile, "%d", ibody);
                AIM_ALLOC(filename, strlen(tetgenInstance->meshInput.outputFileName) + 2 +
                                    strlen("_Vol")+strlen(bodyNumberFile), char, aimInfo, status);
            } else {
                AIM_ALLOC(filename, strlen(tetgenInstance->meshInput.outputFileName)+2, char, aimInfo, status);
            }

            strcpy(filename, tetgenInstance->meshInput.outputFileName);

            if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {
                strcat(filename,"_Vol");
                strcat(filename,bodyNumberFile);
            }

            if (strcasecmp(tetgenInstance->meshInput.outputFormat, "AFLR3") == 0) {

                status = mesh_writeAFLR3(aimInfo,
                                         filename,
                                         tetgenInstance->meshInput.outputASCIIFlag,
                                         &volumeMesh[ibody],
                                         1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(tetgenInstance->meshInput.outputFormat, "VTK") == 0) {

                status = mesh_writeVTK(aimInfo,
                                       filename,
                                       tetgenInstance->meshInput.outputASCIIFlag,
                                       &volumeMesh[ibody],
                                       1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(tetgenInstance->meshInput.outputFormat, "SU2") == 0) {

                status = mesh_writeSU2(aimInfo,
                                       filename,
                                       tetgenInstance->meshInput.outputASCIIFlag,
                                       &volumeMesh[ibody],
                                       tetgenInstance->meshInput.bndConds.numBND,
                                       tetgenInstance->meshInput.bndConds.bndID,
                                       1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(tetgenInstance->meshInput.outputFormat, "Tecplot") == 0) {

                status = mesh_writeTecplot(aimInfo,
                                           filename,
                                           tetgenInstance->meshInput.outputASCIIFlag,
                                           &volumeMesh[ibody],
                                           1.0);
                AIM_STATUS(aimInfo, status);

            } else if (strcasecmp(tetgenInstance->meshInput.outputFormat, "Nastran") == 0) {

                status = mesh_writeNASTRAN(aimInfo,
                                           filename,
                                           tetgenInstance->meshInput.outputASCIIFlag,
                                           &volumeMesh[ibody],
                                           LargeField,
                                           1.0);
                AIM_STATUS(aimInfo, status);

            } else {
                AIM_ERROR(aimInfo, "Unrecognized mesh format, \"%s\"",
                          tetgenInstance->meshInput.outputFormat);
                status = CAPS_BADVALUE;
                goto cleanup;
            }
            AIM_FREE(filename);
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
        AIM_FREE(volumeMesh);
    }
    AIM_FREE(filename);

    return status;
}


/* the execution code from above should be moved here */
extern "C" int
aimExecute(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc, int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


/* no longer optional and needed for restart */
extern "C" int
aimPostAnalysis(void *instStore, void *aimInfo,
                int restart, capsValue *aimInputs)
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

    aimStorage *tetgenInstance;

    tetgenInstance = (aimStorage *) instStore;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    // Get mesh
    numSurfaceMesh = aimInputs[Surface_Mesh-1].length;
    surfaceMesh    = (meshStruct *)aimInputs[Surface_Mesh-1].vals.AIMptr;
    AIM_NOTNULL(surfaceMesh, aimInfo, status);

    // Create/setup volume meshes
    if (aimInputs[Multiple_Mesh-1].vals.integer == (int) true) {

        AIM_ALLOC(tetgenInstance->meshRef, numBody, aimMeshRef, aimInfo, status);
        tetgenInstance->numMeshRef = numBody;

        for (ibody = 0; ibody < numBody; ibody++) {
          status = aim_initMeshRef(tetgenInstance->meshRef);
          AIM_STATUS(aimInfo, status);
        }

        for (ibody = 0; ibody < numBody; ibody++) {
            snprintf(bodyNumberFile, 42, "tetgen_%d", ibody);
            status = aim_file(aimInfo, bodyNumberFile, aimFile);
            AIM_STATUS(aimInfo, status);
            AIM_STRDUP(tetgenInstance->meshRef[ibody].fileName, aimFile, aimInfo, status);

            tetgenInstance->meshRef[ibody].maps[0].tess = surfaceMesh[ibody].egadsTess;

            snprintf(bodyNumberFile, 42, NODATATRANSFER, ibody);
            if (aim_isFile(aimInfo, bodyNumberFile) == CAPS_SUCCESS) continue;

            tetgenInstance->meshRef[0].nmap = 1;
            AIM_ALLOC(tetgenInstance->meshRef[ibody].maps[0].map, surfaceMesh[ibody].numNode, int, aimInfo, status);
            for (i = 0; i < surfaceMesh[ibody].numNode; i++)
              tetgenInstance->meshRef[ibody].maps[0].map[i] = i+1;
        }

    } else {

        AIM_ALLOC(tetgenInstance->meshRef, 1, aimMeshRef, aimInfo, status);
        tetgenInstance->numMeshRef = 1;
        ibody = 0;

        status = aim_initMeshRef(tetgenInstance->meshRef);
        AIM_STATUS(aimInfo, status);

        // set the filename without extensions where the grid is written for solvers
        snprintf(bodyNumberFile, 42, "tetgen_%d", ibody);
        status = aim_file(aimInfo, bodyNumberFile, aimFile);
        AIM_STATUS(aimInfo, status);
        AIM_STRDUP(tetgenInstance->meshRef[0].fileName, aimFile, aimInfo, status);

        snprintf(bodyNumberFile, 42, NODATATRANSFER, ibody);
        if (aim_isFile(aimInfo, bodyNumberFile) == CAPS_SUCCESS) {
            noDataTransfer = (int)true;
        }

        AIM_ALLOC(tetgenInstance->meshRef[0].maps, numSurfaceMesh, aimMeshTessMap, aimInfo, status);
        tetgenInstance->meshRef[0].nmap = numSurfaceMesh;

        nodeOffSet = 0;
        for (ibody = 0; ibody < numSurfaceMesh; ibody++) {
            tetgenInstance->meshRef[0].maps[ibody].tess = surfaceMesh[ibody].egadsTess;
            tetgenInstance->meshRef[0].maps[ibody].map = NULL;

            if (noDataTransfer == (int)true) continue;

            AIM_ALLOC(tetgenInstance->meshRef[0].maps[ibody].map, surfaceMesh[ibody].numNode, int, aimInfo, status);
            for (i = 0; i < surfaceMesh[ibody].numNode; i++)
              tetgenInstance->meshRef[0].maps[ibody].map[i] = nodeOffSet + i+1;

            nodeOffSet += surfaceMesh[ibody].numNode;
        }
    }

cleanup:
    return status;
}


extern "C" int
aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,  int index,
           char **aoname, capsValue *form)
{
    /*! \page aimOutputsTetGen AIM Outputs
     * The following list outlines the TetGen AIM outputs available through the AIM interface.
     */

    int status = CAPS_SUCCESS;

#ifdef DEBUG
    printf(" tetgenAIM/aimOutputs  index = %d!\n", index);
#endif
    if (index == NumberOfElement) {
        *aoname = EG_strdup("NumberOfElement");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsTetGen
         * - <B> NumberOfElement </B> <br>
         * Number of elements in the surface mesh
         */

    } else if (index == NumberOfNode) {
        *aoname = EG_strdup("NumberOfNode");
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsTetGen
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

        /*! \page aimOutputsTetGen
         * - <B> Volume_Mesh </B> <br>
         * The volume mesh for a link
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);
    }

    AIM_NOTNULL(*aoname, aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*aoname);
    return status;
}


extern "C" int
aimCalcOutput(void *instStore, void *aimInfo, int index, capsValue *val)
{
    int        i, status = CAPS_SUCCESS;
    int        numElement, numNodes;
    int        nVertex, nTri, nQuad, nTet, nPyramid, nPrism, nHex;
    aimStorage *tetgenInstance;
    aimMesh    mesh;

#ifdef DEBUG
    printf(" tetgenAIM/aimCalcOutput  index = %d!\n", index);
#endif
    tetgenInstance = (aimStorage *) instStore;

    if (NumberOfElement == index) {

        // Count the total number of surface elements
        numElement = 0;
        for (i = 0; i < tetgenInstance->numMeshRef; i++) {
            status = aim_readBinaryUgridHeader(aimInfo, &tetgenInstance->meshRef[i],
                                               &nVertex, &nTri, &nQuad,
                                               &nTet, &nPyramid, &nPrism, &nHex);
            AIM_STATUS(aimInfo, status);

            numElement += nTri + nQuad + nTet + nPyramid + nPrism + nHex;
        }

        val->vals.integer = numElement;

    } else if (NumberOfNode == index) {

        // Count the total number of surface vertices
        numNodes = 0;
        for (i = 0; i < tetgenInstance->numMeshRef; i++) {
             status = aim_readBinaryUgridHeader(aimInfo, &tetgenInstance->meshRef[i],
                                                &nVertex, &nTri, &nQuad,
                                                &nTet, &nPyramid, &nPrism, &nHex);
             AIM_STATUS(aimInfo, status);

             numNodes += nVertex;
        }

        val->vals.integer = numNodes;

    } else if (Volume_Mesh == index) {

        for (i = 0; i < tetgenInstance->numMeshRef; i++) {
            status = aim_queryMeshes( aimInfo, Volume_Mesh, &tetgenInstance->meshRef[i] );
            if (status > 0) {
                mesh.meshData = NULL;
                mesh.meshRef = &tetgenInstance->meshRef[i];

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

        // Return the volume mesh references
        val->nrow        = tetgenInstance->numMeshRef;
        val->vals.AIMptr = tetgenInstance->meshRef;

    } else {

        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);

    }

cleanup:

    return status;
}

extern "C" void
aimCleanup(void *instStore)
{
    int status; // Function return status
    aimStorage *tetgenInstance;

#ifdef DEBUG
    printf(" tetgenAIM/aimCleanup!\n");
#endif
    tetgenInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(tetgenInstance);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, tetgenAIM  aimStorage cleanup!!!\n", status);

    EG_free(tetgenInstance);
}
