/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AFLR4 AIM
 *
 *      Modified from code written by Dr. Ryan Durscher AFRL/RQVC
 *
 *      This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 *
 */

/*!\mainpage Introduction
 *
 * \section overviewAFLR4 AFLR4 AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to interact with the
 * unstructured, surface grid generator AFLR4 \cite Marcum1995 \cite Marcum1998.
 *
 * The AFLR4 AIM provides the CAPS users with the ability to generate "unstructured, 3D surface grids" using an
 * "Advancing-Front/Local-Reconnection (AFLR) procedure." Both triangular and quadrilateral elements are supported.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsAFLR4 and
 * \ref aimOutputsAFLR4 and \ref attributeAFLR4, respectively.
 * The complete AFLR documentation is available at the <a href="https://www.simcenter.msstate.edu/software/documentation/system/index.html">SimCenter</a>.
 *
 * Example surface meshes:
 *  \image html wing2DAFRL4.png "AFLR4 meshing example - 2D Airfoil" width=500px
 *  \image latex wing2DAFRL4.png "AFLR4 meshing example - 2D Airfoil" width=5cm
 *
 *  \image html sphereAFRL4.png "AFLR4 meshing example - Sphere" width=500px
 *  \image latex sphereAFRL4.png "AFLR4 meshing example - Sphere" width=5cm
 *
 *  \image html wingAFRL4.png "AFLR4 meshing example - Wing" width=500px
 *  \image latex wingAFRL4.png "AFLR4 meshing example - Wing" width=5cm
 *
 *  \section clearanceAFLR4 Clearance Statement
 *  This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.
 */


/*! \page attributeAFLR4 AIM Attributes
 * The following list of attributes are available to guide the mesh generation with the AFLR4 AIM.
 *
 * - <b> capsMeshLength</b> [Required BODY attribute] This numeric <c>BODY</c> attribute
 * sets the AFLR4 <c>ref_len</c> input. <c>capsMeshLength</c> should be a a positive value representative
 * of a characteristic length of the geometry, e.g. the MAC of a wing or diameter of a fuselage.
 * The AIM input <c>Mesh_Length_Factor</c> may be used to apply a global scaling to increase or decrease the mesh resolution.<br>
 *
 * From the AFLR4 documentation:<br><br>
 * ref_len:<br>
 * Reference length for components/bodies in grid units. Reference length should<br>
 * be set to a physically relevant characteristic length for the configuration<br>
 * such as wing chord length or pipe diameter. If ref_len = 0 then it will be set<br>
 * to the bounding box for the largest component/body of interest.<br>
 * The parameters ref_len, max_scale, min_scale and abs_min_scale are all used to<br>
 * set spacing values on all component/body surfaces (those that are not on the<br>
 * farfield or symmetry plane-if any).<br>
 * <br>
 * max_spacing = max_scale * ref_len<br>
 * min_spacing = min_scale * ref_len<br>
 * abs_min_spacing = abs_min_scale * ref_len<br>
 *
 * - <b> AFLR_GBC</b> [Optional FACE attribute: Default STD_UG3_GBC] This string <c>FACE</c> attribute informs AFLR4 what BC
 * treatment should be employed for each geometric <c>FACE</c>. The BC
 * defaults to the string STD_UG3_GBC if none is specified.
 *
 * Predefined AFLR Grid BC string values are:
 *
 * |AFLR_GBC String        | Description                                                         |
 * |:----------------------|:--------------------------------------------------------------------|
 * |FARFIELD_UG3_GBC       | farfield surface<br>same as a standard surface except w/AFLR4       |
 * |STD_UG3_GBC            | standard surface                                                    |
 * |-STD_UG3_GBC           | standard BL generating surface                                      |
 * |BL_INT_UG3_GBC         | symmetry or standard surface<br>that intersects BL region           |
 * |TRANSP_SRC_UG3_GBC     | embedded/transparent surface<br>converted to source nodes by AFLR   |
 * |TRANSP_BL_INT_UG3_GBC  | embedded/transparent surface<br>that intersects BL region           |
 * |TRANSP_UG3_GBC         | embedded/transparent surface                                        |
 * |-TRANSP_UG3_GBC        | embedded/transparent BL generating surface                          |
 * |TRANSP_INTRNL_UG3_GBC  | embedded/transparent surface<br>converted to internal faces by AFLR |
 * |FIXED_BL_INT_UG3_GBC   | fixed surface with BL region<br>that intersects BL region           |
 *
 * Within AFLR4 the grid BC determines how automatic spacing is applied. Their are
 * four basic Grid BC types that are each treated differently.<br>
 * <br>
 * 1. Faces/surfaces that are part of the farfield should be given a
 * FARFIELD_UG3_GBC Grid BC. Farfield faces/surfaces are given a uniform spacing
 * independent of other faces/surfaces with different Grid BCs.<br>
 * <br>
 * 2. Faces/surfaces that represent standard solid surfaces should be given either
 * a STD_UG3_GBC or -STD_UG3_GBC (BL generating) Grid BC. Standard surfaces are
 * given a curvature dependent spacing that may be modified by proximity checking.<br>
 * <br>
 * 3. Faces/surfaces that intersect a BL region should be given either a
 * BL_INT_UG3_GBC (standard boundary surface) or TRANSP_BL_INT_UG3_GBC (embedded/
 * transparent surface with volume mesh on both sides) Grid BC. A common example
 * for the BL_INT_UG3_GBC Grid BC is a symmetry plane. Faces/surfaces set as BL
 * intersecting surfaces are excluded from auto spacing calculations within AFLR4
 * and use edge spacing derived from their neighbors.<br>
 * <br>
 * 4. Surfaces set as transparent surfaces will have a volume mesh on both sides.
 * They can have free edges and can have non-manifold connections to standard
 * solid surfaces and/or BL intersecting surfaces. Vertices in the final surface
 * mesh are not duplicated at non-manifold connections. Transparent surfaces use
 * curvature driven surface spacing as used on standard solid surfaces. However,
 * at non-manifold connections with standard solid surfaces they inherit the
 * surface spacing set on the solid surface they are attached to. They are also
 * excluded from proximity checking. Typical examples of transparent surfaces
 * include wake sheets or multi-material interface surfaces.<br>
 * <br>
 * - <b> AFLR4_Cmp_ID</b> [Optional FACE attribute]<br>
 * EGADS attribute AFLR4_Cmp_ID represents the component identifier for a given
 * face/surface. Component IDs are used for proximity checking. Proximity is only
 * checked between different components. A component is one or more CAD surfaces
 * that represent a component of the full configuration that should be treated
 * individually. For example, a wing-body-strut-nacelle configuration could be
 * considered as four components with wing surfaces set to component 1, body
 * surfaces set to component 2, nacelle surfaces set to 3, and store surfaces set
 * to 4. If each component is a topologically closed surface/body then there is
 * no need to set components. If component IDS are not specified then component
 * identifiers are set for each body defined in the EGADS model or topologically
 * closed surfaces/bodies of the overall configuration. Proximity checking is
 * disabled if there is only one component/body defined. Note that proximity
 * checking only applies to standard surfaces. Component identifiers are set by
 * one of three methods, chosen in the following order.<br>
 * <br>
 * 1. If defined by EGADS attribute AFLR4_Cmp_ID then attribute sets component
 *    identifier.<br>
 * 2. Else, if multiple bodies are defined in the EGADS model then bodies index is
 *    used to set component identifier.<br>
 * 3. Else, component identifiers are set an index based on topologically closed
 *    surfaces/bodies of the overall configuration.<br>
 * <br>
 * - <b> AFLR4_Isolated_Edge_Refinement_Flag </b> [Optional FACE attribute: Integer Range 0 to 2]<br>
 * Isolated edge refinement flag.
 * If Flag = 0 then do not refine isolated interior edges.
 * If Flag = 1 then refine isolated interior edges if the surface has local
 * curvature (as defined using cier).
 * If Flag = 2 then refine all isolated interior edges.
 * An isolated interior edges is connected only to boundary nodes. Isolated edges
 * are refined by placing a new node in the middle of of the edge.
 * Note that if not set then the isolated edge refinement flag is set to the global
 * value AFLR4_mier.
 * <br>
 * - <b> AFLR4_Edge_Refinement_Weight</b> [Optional FACE attribute: Default 0.0, Range 0 to 1]<br>
 * EGADS attribute AFLR4_Edge_Refinement_Weight represents the edge mesh spacing
 * scale factor weight for a given face/surface. Edge mesh spacing can be scaled
 * on a given face/surface based on the discontinuity level between adjacent
 * faces/surfaces on both sides of the edge. The edge mesh spacing scale factor
 * weight set with AFLR4_Edge_Refinement_Weight is used as an interpolation
 * weight between the unmodified spacing and the modified spacing. A value of one
 * applies the maximum modification and a value of zero applies no change in edge
 * spacing. Note that no modification is done to edges that belong to farfield or
 * BL intersecting face/surface.
 * <br>
 * - <b> AFLR4_Scale_Factor</b> [Optional FACE or EDGE attribute: Default 1.0]<br>
 * EGADS attribute AFLR4_Scale_Factor represents the AFLR4 mesh spacing
 * scale factor for a given face/edge. Curvature dependent spacing can be
 * scaled on the face/edge by the value of the scale factor set with
 * AFLR4_Scale_Factor.<br>
 * <br>
 * - <b> AFLR4_quad_local</b> [Optional FACE attribute: Default 0.0]<br>
 * Local quad-face combination flag.<br>
 * If AFLR4_quad_local is set to 0, then do not combine tria-faces to form
 * quad-faces.<br>
 * If AFLR4_quad_local is set to 1, then combine tria-faces to form quad-faces.<br>
 * This option also locally selects advancing-point point placement rather than
 * the default advancing-front point placement.
 * Note that if not set then the local quad-face combination flag is set to the
 * global quad-face combination flag AFLR4_quad.
 */

#ifdef WIN32
#define strcasecmp  stricmp
#define strncasecmp _strnicmp
typedef int         pid_t;
#endif

#include <string.h>
#include <math.h>

#include "aimUtil.h"

#include "meshUtils.h"       // Collection of helper functions for meshing
#include "miscUtils.h"
#include "deprecateUtils.h"

#include <ug/UG_LIB.h>
#include <aflr4/AFLR4_LIB.h> // Bring in AFLR4 API library

#include "aflr4_Interface.h" // Bring in AFLR4 'interface' functions

#define AFLR4FILE "aflr4_%d"

//#define DEBUG

enum aimOutputs
{
  Done = 1,                    /* index is 1-based */
  NumberOfElement,
  NumberOfNode,
  Surface_Mesh,
  NUMOUT = Surface_Mesh        /* Total number of outputs */
};


typedef struct {

    // Container for surface mesh
    int numMeshRef;
    aimMeshRef *meshRef;

    // Container for mesh input
    meshInputStruct meshInput;

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    mapAttrToIndexStruct meshMap;

    int numElemTotal;
    int numNodeTotal;

} aimStorage;


static int destroy_aimStorage(aimStorage *aflr4Instance, int inUpdate)
{
    int i, status; // Function return status

    // Destroy meshInput
    status = destroy_meshInputStruct(&aflr4Instance->meshInput);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr4AIM meshInput cleanup!!!\n", status);

    // Destroy surface mesh allocated arrays
    for (i = 0; i < aflr4Instance->numMeshRef; i++) {
      status = aim_freeMeshRef(&aflr4Instance->meshRef[i]);
      if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr4AIM surfaceMesh cleanup!!!\n", status);
    }
    AIM_FREE(aflr4Instance->meshRef);
    aflr4Instance->numMeshRef = 0;

    aflr4Instance->numElemTotal = 0;
    aflr4Instance->numNodeTotal = 0;

    if (inUpdate == (int)true) return status;

    // Destroy attribute to index map
    status = destroy_mapAttrToIndexStruct(&aflr4Instance->groupMap);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr4AIM attributeMap cleanup!!!\n", status);

    status = destroy_mapAttrToIndexStruct(&aflr4Instance->meshMap);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr4AIM attributeMap cleanup!!!\n", status);

    return CAPS_SUCCESS;
}


#ifdef WIN32
#if _MSC_VER >= 1900
    // FILE _iob[] = {*stdin, *stdout, *stderr};
    FILE _iob[3];

    extern FILE * __cdecl __iob_func(void)
    {
        return _iob;
    }
#endif
#endif


static int setAFLR4Attr(void *aimInfo,
                        ego body,
                        mapAttrToIndexStruct *meshMap,
                        int numMeshProp,
                        meshSizingStruct *meshProp)
{

    int status; // Function return status

    int iface, iedge; // Indexing

    int numFace, numEdge;
    ego *faces = NULL, *edges = NULL;

    const char *meshName = NULL;
    int attrIndex = 0;
    int propIndex = 0;
    const char* bcType = NULL;

    status = EG_getBodyTopos(body, NULL, FACE, &numFace, &faces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(faces, aimInfo, status);

    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &edges);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(edges, aimInfo, status);

    // Loop through the faces and set AFLR attributes
    for (iface = 0; iface < numFace; iface++) {

        status = retrieve_CAPSMeshAttr(faces[iface], &meshName);
        if ((status == EGADS_SUCCESS) && (meshName != NULL)) {

            status = get_mapAttrToIndexIndex(meshMap, meshName, &attrIndex);
            AIM_STATUS(aimInfo, status);

            for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

                // Check if the mesh property applies to the capsMesh of this face
                if (meshProp[propIndex].attrIndex != attrIndex) continue;

                // If bcType specified in meshProp
                if (meshProp[propIndex].bcType != NULL) {

                    bcType = meshProp[propIndex].bcType;
                    if      (strncasecmp(meshProp[propIndex].bcType, "Farfield"  ,  8) == 0)
                        bcType = "FARFIELD_UG3_GBC";
                    else if (strncasecmp(meshProp[propIndex].bcType, "Freestream", 10) == 0)
                        bcType = "FARFIELD_UG3_GBC";
                    else if (strncasecmp(meshProp[propIndex].bcType, "Viscous"   ,  7) == 0)
                        bcType = "-STD_UG3_GBC";
                    else if (strncasecmp(meshProp[propIndex].bcType, "Inviscid"  ,  8) == 0)
                        bcType = "STD_UG3_GBC";
                    else if (strncasecmp(meshProp[propIndex].bcType, "Symmetry"  ,  8) == 0 ||
                             strncasecmp(meshProp[propIndex].bcType, "BoundaryLayerIntersect",22) == 0)
                        bcType = "BL_INT_UG3_GBC";

                    // add the BC attribute
                    status = EG_attributeAdd(faces[iface], "AFLR_GBC", ATTRSTRING, 0, NULL, NULL, bcType);
                    AIM_STATUS(aimInfo, status);
                }

                // If scaleFactor specified in meshProp
                if (meshProp[propIndex].scaleFactor > 0) {

                    // add the scale factor attribute
                    status = EG_attributeAdd(faces[iface], "AFLR4_Scale_Factor", ATTRREAL, 1, NULL, &meshProp[propIndex].scaleFactor, NULL);
                    AIM_STATUS(aimInfo, status);
                }

                // If edgeWeight specified in meshProp
                if (meshProp[propIndex].edgeWeight >= 0) {

                    // add the edge scale factor weight attribute
                    status = EG_attributeAdd(faces[iface], "AFLR4_Edge_Refinement_Weight", ATTRREAL, 1, NULL, &meshProp[propIndex].edgeWeight, NULL);
                    AIM_STATUS(aimInfo, status);
                }
            }
        }
    } // Face loop


    // Loop through the edges and set AFLR attributes
    for (iedge = 0; iedge < numEdge; iedge++) {

        status = retrieve_CAPSMeshAttr(edges[iedge], &meshName);
        if ((status == EGADS_SUCCESS) && (meshName != NULL)) {

            status = get_mapAttrToIndexIndex(meshMap, meshName, &attrIndex);
            AIM_STATUS(aimInfo, status);

            for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

                // Check if the mesh property applies to the capsMesh of this face
                if (meshProp[propIndex].attrIndex != attrIndex) continue;

                // If scaleFactor specified in meshProp
                if (meshProp[propIndex].scaleFactor > 0) {

                    // add the scale factor attribute
                    status = EG_attributeAdd(edges[iedge], "AFLR4_Scale_Factor", ATTRREAL, 1, NULL, &meshProp[propIndex].scaleFactor, NULL);
                    AIM_STATUS(aimInfo, status);
                }
            }
        }
    } // Edge loop

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(faces);
    AIM_FREE(edges);

    return status;
}



/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
    int        status; // Function return
    aimStorage *aflr4Instance=NULL;

#ifdef WIN32
#if _MSC_VER >= 1900
    _iob[0] = *stdin;
    _iob[1] = *stdout;
    _iob[2] = *stderr;
#endif
#endif

#ifdef DEBUG
    printf("\n aflr4AIM/aimInitialize   instance = %d!\n", inst);
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
    AIM_ALLOC(aflr4Instance, 1, aimStorage, aimInfo, status);
    *instStore = aflr4Instance;

    // Set initial values for aflrInstance

    // Container for surface meshes
    aflr4Instance->numMeshRef = 0;
    aflr4Instance->meshRef = NULL;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&aflr4Instance->meshMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&aflr4Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    // Container for mesh input
    status = initiate_meshInputStruct(&aflr4Instance->meshInput);
    AIM_STATUS(aimInfo, status);

    aflr4Instance->numElemTotal = 0;
    aflr4Instance->numNodeTotal = 0;

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);
    return status;
}


// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
    /*! \page aimInputsAFLR4 AIM Inputs
     * The following list outlines the AFLR4 meshing options along with their default value available
     * through the AIM interface.
     *
     * Please consult the <a href="https://www.simcenter.msstate.edu/software/documentation/aflr4/index.html">AFLR4 documentation</a> for default values not present here.
     */

    int status; // error code
    UG_Param_Struct *AFLR4_Param_Struct_Ptr = NULL; // AFRL4 input structure used to get default values

    int min_ncell, mer_all;
    CHAR_UG_MAX no_prox;

    double ff_cdfr, abs_min_scale, BL_thickness, Re_l, curv_factor,
           max_scale, min_scale, erw_all; //, ref_len;


    status = ug_malloc_param(&AFLR4_Param_Struct_Ptr);
    if ((status != CAPS_SUCCESS) || (AFLR4_Param_Struct_Ptr == NULL)) {
        printf("ug_malloc_param status = %d\n", status);
        goto cleanup;
    }

    status = ug_initialize_param(4, AFLR4_Param_Struct_Ptr);
    if (status != CAPS_SUCCESS) {
        printf("ug_initialize_param status = %d\n", status);
        goto cleanup;
    }

    status = aflr4_initialize_param(AFLR4_Param_Struct_Ptr);
    if (status != CAPS_SUCCESS) {
        printf("aflr4_initialize_param status = %d\n", status);
        goto cleanup;
    }

#ifdef DEBUG
    printf(" aflr4AIM/aimInputs index = %d!\n", index);
#endif

    // Meshing Inputs
    if (index == Proj_Name) {
        *ainame              = EG_strdup("Proj_Name");
        defval->type         = String;
        defval->nullVal      = NotAllowed;
        AIM_STRDUP(defval->vals.string, "aflr4_CAPS", aimInfo, status);

        /*! \page aimInputsAFLR4
         * - <B> Proj_Name = "aflr4_CAPS"</B> <br>
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

        /*! \page aimInputsAFLR4
         * \include{doc} Mesh_Format.dox
         */

    } else if (index == Mesh_Quiet_Flag) {
        *ainame               = AIM_NAME(Mesh_Quiet_Flag);
        defval->type          = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsAFLR4
         * - <B> Mesh_Quiet_Flag = False</B> <br>
         * Complete suppression of mesh generator (not including errors)
         */

    } else if (index == Mesh_Gen_Input_String) {
        *ainame               = AIM_NAME(Mesh_Gen_Input_String);
        defval->type          = String;
        defval->nullVal       = IsNull;
        defval->vals.string   = NULL;

        /*! \page aimInputsAFLR4
         * - <B> Mesh_Gen_Input_String = NULL</B> <br>
         * Meshing program command line string (as if called in bash mode). Use this to specify more<br>
         * complicated options/use features of the mesher not currently exposed through other AIM input<br>
         * variables. Note that this is the exact string that will be provided to the volume mesher; no<br>
         * modifications will be made. If left NULL an input string will be created based on default values<br>
         * of the relevant AIM input variables.
         */

    } else if (index == Ff_cdfr) {

        status = ug_get_double_param((char *) "ff_cdfr", &ff_cdfr,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'ff_cdfr'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame              = EG_strdup("ff_cdfr");
        defval->type         = Double;
        defval->dim          = Scalar;
        defval->vals.real    = ff_cdfr;

        /*! \page aimInputsAFLR4
         * - <B>ff_cdfr</B> <br>
         * Farfield growth rate for field point spacing.<br>
         * The farfield spacing is set to a uniform value dependent upon the maximum size <br>
         * of the domain, maximum size of inner bodies, max and min body spacing, and <br>
         * farfield growth rate. ; <br>
         * ff_spacing = (ff_cdfr-1)*L+(min_spacing+max_spacing)/2 ; <br>
         * where L is the approximate distance between inner bodies and farfield. <br>
         */

    } else if (index == Min_ncell) {

        status = ug_get_int_param((char *) "min_ncell", &min_ncell,
                                  AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'min_ncell'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame              = EG_strdup("min_ncell");
        defval->type         = Integer;
        defval->dim          = Scalar;
        defval->vals.integer  = min_ncell;

        /*! \page aimInputsAFLR4
         * - <B>min_ncell</B> <br>
         * Minimum number of cells between two components/bodies.<br>
         * Proximity of components/bodies<br>
         * to each other is estimated and surface<br>
         * spacing is locally reduced if needed. Local<br>
         * surface spacing is selectively reduced when<br>
         * components/bodies are close and their<br>
         * existing local surface spacing would generate<br>
         * less than the minimum number of cells<br>
         * specified by min_ncell. Proximity checking is<br>
         * automatically disabled if min_ncell=1 or if<br>
         * there is only one component/body defined.
         */

    } else if (index == Mer_all) {

        status = ug_get_int_param((char *) "mer_all", &mer_all,
                                  AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'mer_all'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame              = EG_strdup("mer_all");
        defval->type         = Integer;
        defval->dim          = Scalar;
        defval->vals.integer = mer_all;

        /*! \page aimInputsAFLR4
         * - <B>mer_all</B> <br>
         * Global edge mesh spacing scale factor flag.<br>
         * Edge mesh spacing can be scaled on all surfaces based on discontinuity level<br>
         * between adjacent surfaces on both sides of the edge. For each surface the level<br>
         * of discontinuity (as defined by angerw1 and angerw2) determines the edge<br>
         * spacing scale factor for potentially reducing the edge spacing. See erw_ids and<br>
         * erw_list. This option is equivalent to setting erw_ids equal to the list of all<br>
         * surface IDS and the edge mesh spacing scale factor weight in erw_list equal to<br>
         * one. Note that no modification is done to edges that belong to surfaces with a<br>
         * grid BC of farfield (ff_ids) or BL intersecting (int_ids).
         */

    } else if (index == No_prox) {

        status = ug_get_char_param((char *) "-no_prox", no_prox,
                                   AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'no_prox'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame              = EG_strdup("no_prox");
        defval->type         = Boolean;
        defval->dim          = Scalar;
        defval->vals.integer = False;

        /*! \page aimInputsAFLR4
         * - <B>no_prox</B> <br>
         * Disable proximity check flag.<br>
         * If no_prox=False then proximity of components/bodies to each other is estimated <br>
         * and surface spacing is locally reduced if needed. <br>
         * If no_prox=True or if there is only one component/body defined then proximity <br>
         * checking is disabled. <br>
         *
         */

    } else if (index == Abs_min_scale) {

        status = ug_get_double_param((char *) "abs_min_scale", &abs_min_scale,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'abs_min_scale'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("abs_min_scale");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = abs_min_scale;

        /*! \page aimInputsAFLR4
         * - <B>abs_min_scale</B> <br>
         * Relative scale of absolute minimum spacing to<br>
         * reference length. The relative scale of<br>
         * absolute minimum spacing to reference length<br>
         * (ref_len) controls the absolute minimum<br>
         * spacing that can be set on any component/body<br>
         * surface by proximity checking (see min_ncell).<br>
         * Note that the value of abs_min_scale is limited to be less<br>
         * than or equal to min_scale.
         */

    } else if (index == BL_Thickness) {

        status = ug_get_double_param((char *) "BL_thickness", &BL_thickness,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'BL_thickness'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("BL_thickness");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = BL_thickness;

        /*! \page aimInputsAFLR4
         * - <B>BL_thickness</B> <br>
         * Boundary layer thickness for proximity checking. <br>
         * Proximity of components/bodies to each other is estimated and surface spacing <br>
         * is locally reduced if needed. Note that if the Reynolds Number, Re_l, is set <br>
         * then the BL_thickness value is set to an estimate for turbulent flow. If the <br>
         * set or calculated value of BL_thickness>0 then the boundary layer thickness is <br>
         * included in the calculation for the required surface spacing during proximity <br>
         * checking.
         */

    } else if (index == RE_l) {

        status = ug_get_double_param((char *) "Re_l", &Re_l,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'Re_l'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("Re_l");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = Re_l;

        /*! \page aimInputsAFLR4
         * - <B>Re_l</B> <br>
         * Reynolds Number for estimating BL thickness.<br>
         * The Reynolds Number based on reference length, Re_l, (if set) along with <br>
         * reference length, ref_len, are used to estimate the BL thickness, BL_thickness, <br>
         * for turbulent flow. If Re_l>0 then this estimated value is used to set <br>
         * BL_thickness.
         */

    } else if (index == Curv_factor) {

        status = ug_get_double_param((char *) "curv_factor", &curv_factor,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'curv_factor'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("curv_factor");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = curv_factor;

        /*! \page aimInputsAFLR4
         * - <B>curv_factor</B> <br>
         * Curvature factor<br>
         * For surface curvature the spacing is derived from the curvature factor divided<br>
         * by the curvature.<br>
         * Curvature = 1 / Curvature_Radius<br>
         * Spacing = curv_factor / Curvature<br>
         * The resulting spacing between is limited by the minimum and maximum spacing set<br>
         * by min_scale and max_scale. Note that if curv_factor=0 then surface curvature<br>
         * adjustment is not used.
         */

    } else if (index == Erw_all) {

        status = ug_get_double_param((char *) "erw_all", &erw_all,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'erw_all'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("erw_all");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = erw_all;

        /*! \page aimInputsAFLR4
         * - <B>erw_all</B> <br>
         * Global edge mesh spacing refinement weight.<br>
         * Edge mesh spacing can be scaled on all surfaces (if mer_all=1) based on<br>
         * discontinuity level between adjacent surfaces on both sides of the edge.<br>
         * For each surface the level of discontinuity (as defined by angerw1 and angerw2)<br>
         * determines the edge spacing scale factor for potentially reducing the edge<br>
         * spacing. The edge mesh spacing scale factor weight is then used as an<br>
         * interpolation weight between the unmodified spacing and the modified spacing.<br>
         * A value of one applies the maximum modification and a value of zero applies no<br>
         * change in edge spacing. If the global edge mesh spacing scale factor flag,<br>
         * mer_all, is set to 1 then that is equivalent to setting AFLR_Edge_Scale_Factor_Weight<br>
         * on all FACEs to the value erw_all. Note that no modification is<br>
         * done to edges that belong to surfaces with a grid BC of farfield (FARFIELD_UG3_GBC) or BL<br>
         * intersecting. Also, note that the global weight, erw_all, is not<br>
         * applicable if mer_all=0.
         */

    } else if (index == Max_scale) {

        status = ug_get_double_param((char *) "max_scale", &max_scale,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'max_scale'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("max_scale");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = max_scale;

        /*! \page aimInputsAFLR4
         * - <B>max_scale</B> <br>
         * Relative scale of maximum spacing to<br>
         * reference length. The relative scale of<br>
         * maximum spacing to reference length (ref_len)<br>
         * controls the maximum spacing that can be set<br>
         * on any component/body surface.
         */

    } else if (index == Min_scale) {

        status = ug_get_double_param((char *) "min_scale", &min_scale,
                                     AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'min_scale'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }

        *ainame           = EG_strdup("min_scale");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = min_scale;

        /*! \page aimInputsAFLR4
         * - <B>min_scale</B> <br>
         * Relative scale of minimum spacing to<br>
         * reference length. The relative scale of<br>
         * minimum spacing to reference length (ref_len)<br>
         * controls the minimum spacing that can be set<br>
         * on any component/body surface.
         */

    } else if (index == Mesh_Length_Factor) {

      /* There is no reasonable default for the ref_len parameter,
       * the user must always set it via capsMeshLength and Mesh_Length_Factor
       *
        status = ug_get_double_param ((char*)"ref_len", &ref_len, AFLR4_Param_Struct_Ptr);
        if (status == 1) status = CAPS_SUCCESS;
        if (status != CAPS_SUCCESS) {
            printf("Failed to retrieve default value for 'ref_len'\n");
            status = CAPS_NOTFOUND;
            goto cleanup;
        }
     */

        *ainame           = EG_strdup("Mesh_Length_Factor");
        defval->type      = Double;
        defval->dim       = Scalar;
        defval->vals.real = 1;
        defval->nullVal   = NotNull;

        /*! \page aimInputsAFLR4
         * - <B>Mesh_Length_Factor = 1</B> <br>
         * Global scaling factor to increase or decrease mesh resolution.<br>
         *
         * ref_len:<br>
         * Reference length for components/bodies in grid units. Reference length should<br>
         * be set to a physically relevant characteristic length for the configuration<br>
         * such as wing chord length or pipe diameter. If ref_len = 0 then it will be set<br>
         * to the bounding box for the largest component/body of interest.<br>
         * The parameters ref_len, max_scale, min_scale and abs_min_scale are all used to<br>
         * set spacing values on all component/body surfaces (those that are not on the<br>
         * farfield or symmetry plane-if any).<br>
         * <br>
         * max_spacing = max_scale * ref_len<br>
         * min_spacing = min_scale * ref_len<br>
         * abs_min_spacing = abs_min_scale * ref_len<br>
         */

    } else if (index == Mesh_Sizing) {
        *ainame              = EG_strdup("Mesh_Sizing");
        defval->type         = Tuple;
        defval->nullVal      = IsNull;
        //defval->units        = NULL;
        defval->dim          = Vector;
        defval->lfixed       = Change;
        defval->vals.tuple   = NULL;

        /*! \page aimInputsAFLR4
         * - <B> Mesh_Sizing = NULL </B> <br>
         * See \ref meshSizingProp for additional details.
         */

    } else if (index == Multiple_Mesh) {
        *ainame               = EG_strdup("Multiple_Mesh");
        defval->type          = String;
        AIM_STRDUP(defval->vals.string, "SingleDomain", aimInfo, status);

        /*! \page aimInputsAFLR4
         * - <B> Multiple_Mesh = "SingleDomain"</B> <br>
         * If "SingleDomain": Generate a single surface mesh file is assuming multiple
         * bodies define a single computational domain (i.e. CFD)<br>
         * <br>
         * If "MultiFile": Generate a surface mesh file for each body.<br>
         * <br>
         * If "MultiDomain": Generate a single mesh file containing multiple surface meshes for each body.<br>
         */

    } else if (index == EGADS_Quad) {
        *ainame              = EG_strdup("EGADS_Quad");
        defval->type         = Boolean;
        defval->vals.integer = false;

        /*! \page aimInputsAFLR4
         * - <B> EGADS_Quad = False </B> <br>
         * If true, apply EGADS quadding to the AFLR4 triangulation.
         */

    } else if (index == AFLR4_Quad) {
        *ainame              = EG_strdup("AFLR4_Quad");
        defval->type         = Boolean;
        defval->vals.integer = false;

        /*! \page aimInputsAFLR4
         * - <B> AFLR4_Quad = False </B> <br>
         * If true, apply -quad flag for AFLR4 quadding.
         */

    } else if (index == Skin) {
        *ainame              = EG_strdup("skin");
        defval->type         = Boolean;
        defval->vals.integer = false;

        /*! \page aimInputsAFLR4
         * - <B> skin = False </B> <br>
         * If true, apply -skin flag to automatically set the grid BCs for structural cases.
         */

    } else {
        AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
        status = CAPS_BADINDEX;
        goto cleanup;
    }

    AIM_NOTNULL(*ainame, aimInfo, status);
    status = CAPS_SUCCESS;

cleanup:
    ug_free_param(AFLR4_Param_Struct_Ptr);
    return status;
}


// ********************** AIM Function Break *****************************
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
    int status; // Function return status

    int i, j, bodyIndex;

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;

    int MultiMesh = -1;
    char bodyNumberFile[128];
    char aimFile[PATH_MAX];

    aimStorage *aflr4Instance;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" aflr4AIM/aimPreAnalysis numBody = %d!\n", numBody);
#endif

    if (numBody <= 0 || bodies == NULL) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }
    AIM_NOTNULL(aimInputs, aimInfo, status);

    aflr4Instance = (aimStorage *) instStore;

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(aflr4Instance, (int)true);
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

    if (aflr4Instance->groupMap.numAttribute == 0 ||
        aim_newGeometry(aimInfo) == CAPS_SUCCESS ) {
        // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
        status = create_CAPSGroupAttrToIndexMap(numBody,
                                                bodies,
                                                -1, // Faces only
                                                &aflr4Instance->groupMap);
        AIM_STATUS(aimInfo, status);
    }

    if (aflr4Instance->meshMap.numAttribute == 0 ||
        aim_newGeometry(aimInfo) == CAPS_SUCCESS ) {
        status = create_CAPSMeshAttrToIndexMap(numBody,
                                               bodies,
                                               2, // Body, Face, and Edge
                                               &aflr4Instance->meshMap);
        AIM_STATUS(aimInfo, status);
    }

    // Setup meshing input structure

    // Get Tessellation parameters -Tess_Params
    aflr4Instance->meshInput.paramTess[0] = 0;
    aflr4Instance->meshInput.paramTess[1] = 0;
    aflr4Instance->meshInput.paramTess[2] = 0;

    aflr4Instance->meshInput.quiet = aimInputs[Mesh_Quiet_Flag-1].vals.integer;

    // Set aflr4 specific mesh inputs
    if (aimInputs[Mesh_Gen_Input_String-1].nullVal != IsNull) {
        AIM_STRDUP(aflr4Instance->meshInput.aflr4Input.meshInputString, aimInputs[Mesh_Gen_Input_String-1].vals.string, aimInfo, status);
    }

    // Get mesh sizing parameters
    if (aimInputs[Mesh_Sizing-1].nullVal != IsNull) {

        status = deprecate_SizingAttr(aimInfo,
                                      aimInputs[Mesh_Sizing-1].length,
                                      aimInputs[Mesh_Sizing-1].vals.tuple,
                                      &aflr4Instance->meshMap,
                                      &aflr4Instance->groupMap);
        AIM_STATUS(aimInfo, status);

        status = mesh_getSizingProp(aimInfo,
                                    aimInputs[Mesh_Sizing-1].length,
                                    aimInputs[Mesh_Sizing-1].vals.tuple,
                                    &aflr4Instance->meshMap,
                                    &numMeshProp,
                                    &meshProp);
        AIM_STATUS(aimInfo, status);

        // apply the sizing attributes
        for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {
/*@-nullpass@*/
            status = setAFLR4Attr(aimInfo,
                                  bodies[bodyIndex],
                                  &aflr4Instance->meshMap,
                                  numMeshProp,
                                  meshProp);
/*@+nullpass@*/
            AIM_STATUS(aimInfo, status);
        }
    }


    if (MultiMesh == 0 || MultiMesh == 2) {

      AIM_ALLOC(aflr4Instance->meshRef, 1, aimMeshRef, aimInfo, status);
      aflr4Instance->numMeshRef = 1;

      status = aim_initMeshRef(aflr4Instance->meshRef, aimSurfaceMesh);
      AIM_STATUS(aimInfo, status);

      AIM_ALLOC(aflr4Instance->meshRef[0].maps, numBody, aimMeshTessMap, aimInfo, status);
      aflr4Instance->meshRef[0].nmap = numBody;

      for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {
        aflr4Instance->meshRef[0].maps[bodyIndex].tess = NULL;
        aflr4Instance->meshRef[0].maps[bodyIndex].map = NULL;
      }

      // set the filename without extensions where the grid is written for solvers
      bodyIndex = 0;
      if (aimInputs[Proj_Name-1].nullVal != IsNull)
        snprintf(bodyNumberFile, 128, "%s_%d", aimInputs[Proj_Name-1].vals.string, bodyIndex);
      else
        snprintf(bodyNumberFile, 128, AFLR4FILE, bodyIndex);
      status = aim_file(aimInfo, bodyNumberFile, aimFile);
      AIM_STATUS(aimInfo, status);
      AIM_STRDUP(aflr4Instance->meshRef[0].fileName, aimFile, aimInfo, status);

    } else  if (MultiMesh == 1) {

      AIM_ALLOC(aflr4Instance->meshRef, numBody, aimMeshRef, aimInfo, status);
      aflr4Instance->numMeshRef = numBody;

      for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {
        status = aim_initMeshRef(&aflr4Instance->meshRef[bodyIndex], aimSurfaceMesh);
        AIM_STATUS(aimInfo, status);

        AIM_ALLOC(aflr4Instance->meshRef[bodyIndex].maps, 1, aimMeshTessMap, aimInfo, status);
        aflr4Instance->meshRef[bodyIndex].nmap = 1;

        aflr4Instance->meshRef[bodyIndex].maps[0].tess = NULL;
        aflr4Instance->meshRef[bodyIndex].maps[0].map = NULL;

        // set the filename without extensions where the grid is written for solvers
        if (aimInputs[Proj_Name-1].nullVal != IsNull)
          snprintf(bodyNumberFile, 128, "%s_%d", aimInputs[Proj_Name-1].vals.string, bodyIndex);
        else
          snprintf(bodyNumberFile, 128, AFLR4FILE, bodyIndex);
        status = aim_file(aimInfo, bodyNumberFile, aimFile);
        AIM_STATUS(aimInfo, status);
        AIM_STRDUP(aflr4Instance->meshRef[bodyIndex].fileName, aimFile, aimInfo, status);
      }
    }

    for (i = 0; i < aflr4Instance->numMeshRef; i++) {

      AIM_ALLOC(aflr4Instance->meshRef[i].bnds, aflr4Instance->groupMap.numAttribute, aimMeshBnd, aimInfo, status);
      aflr4Instance->meshRef[i].nbnd = aflr4Instance->groupMap.numAttribute;
      for (j = 0; j < aflr4Instance->meshRef[i].nbnd; j++) {
        status = aim_initMeshBnd(aflr4Instance->meshRef[i].bnds + j);
        AIM_STATUS(aimInfo, status);
      }

      for (j = 0; j < aflr4Instance->meshRef[i].nbnd; j++) {
        AIM_STRDUP(aflr4Instance->meshRef[i].bnds[j].groupName, aflr4Instance->groupMap.attributeName[j], aimInfo, status);
        aflr4Instance->meshRef[i].bnds[j].ID = aflr4Instance->groupMap.attributeIndex[j];
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

    int i; // Indexing

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    // File output
    char *filename = NULL;

    const aimStorage *aflr4Instance;

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

#ifdef DEBUG
    printf(" aflr4AIM/aimPreAnalysis numBody = %d!\n", numBody);
#endif

    if (numBody <= 0 || bodies == NULL) {
        AIM_ERROR(aimInfo, "No Bodies!");
        return CAPS_SOURCEERR;
    }
    AIM_NOTNULL(aimInputs, aimInfo, status);

    aflr4Instance = (const aimStorage *) instStore;

    // remove previous meshes
    for (i = 0; i < aflr4Instance->numMeshRef; i++) {
      status = aim_deleteMeshes(aimInfo, &aflr4Instance->meshRef[i]);
      AIM_STATUS(aimInfo, status);
    }

    status = aflr4_Surface_Mesh(aimInfo,
                                aflr4Instance->meshInput.quiet,
                                numBody, bodies,
                                aimInputs,
                                aflr4Instance->meshInput);
    AIM_STATUS(aimInfo, status, "Problem during AFLR4 surface meshing");

    status = CAPS_SUCCESS;

cleanup:

    AIM_FREE(filename);
    return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ const void *instStore, /*@unused@*/ void *aimStruc,
               int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis( void *aimStore, void *aimInfo,
                     int restart,    capsValue *aimInputs)
{

    int status = CAPS_SUCCESS;
    int bodyIndex;
    int MultiMesh = 0;

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

    aimStorage *aflr4Instance;
    aflr4Instance = (aimStorage *) aimStore;

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
      snprintf(bodyNumber, 42, AFLR4TESSFILE, bodyIndex);
      status = aim_file(aimInfo, bodyNumber, aimFile);
      AIM_STATUS(aimInfo, status);

      status = EG_loadTess(bodies[bodyIndex], aimFile, &tess);
      AIM_STATUS(aimInfo, status);

      status = aim_newTess(aimInfo, tess);
      AIM_STATUS(aimInfo, status);

      status = EG_statusTessBody(tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      if (MultiMesh == 0 || MultiMesh == 2) {

        aflr4Instance->meshRef[0].maps[bodyIndex].tess = tess;

        AIM_ALLOC(aflr4Instance->meshRef[0].maps[bodyIndex].map, nglobal, int, aimInfo, status);
        for (i = 0; i < nglobal; i++) aflr4Instance->meshRef[0].maps[bodyIndex].map[i] = iglobal++;

      } else if (MultiMesh == 1) {

        aflr4Instance->meshRef[bodyIndex].maps[0].tess = tess;

        AIM_ALLOC(aflr4Instance->meshRef[bodyIndex].maps[0].map, nglobal, int, aimInfo, status);
        for (i = 0; i < nglobal; i++) aflr4Instance->meshRef[bodyIndex].maps[0].map[i] = iglobal++;

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

      aflr4Instance->numNodeTotal += nglobal;
      aflr4Instance->numElemTotal += numTri+numQuad;
    }

    if (restart == 0 &&
        aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false) {
        printf("----------------------------\n");
        printf("Total number of nodes    = %d\n", aflr4Instance->numNodeTotal);
        printf("Total number of elements = %d\n", aflr4Instance->numElemTotal);
    }

    /* Explicitly write out any requested meshes */
    for (i = 0; i < aflr4Instance->numMeshRef; i++) {
        status = aim_queryMeshes( aimInfo, Mesh_Format, ANALYSISIN, &aflr4Instance->meshRef[i] );
        if (status > 0) {
/*@-immediatetrans@*/
            mesh.meshData = NULL;
            mesh.meshRef = &aflr4Instance->meshRef[i];
/*@+immediatetrans@*/

            status = mesh_surfaceMeshData(aimInfo, &aflr4Instance->groupMap, &mesh);
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
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               /*@unused@*/ int index, char **aoname, capsValue *form)
{
    int status = CAPS_SUCCESS;

    /*! \page aimOutputsAFLR4 AIM Outputs
     * The following list outlines the AFLR4 AIM outputs available through the AIM interface.
     */

#ifdef DEBUG
    printf(" aflr4AIM/aimOutputs index = %d!\n", index);
#endif

    if (index == Node) {
        *aoname = AIM_NAME(Done);
        form->type = Boolean;
        form->vals.integer = (int) false;

        /*! \page aimOutputsAFLR4
         * - <B> Done </B> <br>
         * True if a surface mesh was created on all surfaces, False if not.
         */

    } else if (index == NumberOfElement) {
        *aoname = AIM_NAME(NumberOfElement);
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsAFLR4
         * - <B> NumberOfElement </B> <br>
         * Number of elements in the surface mesh
         */

    } else if (index == NumberOfNode) {
        *aoname = AIM_NAME(NumberOfNode);
        form->type = Integer;
        form->vals.integer = 0;

        /*! \page aimOutputsAFLR4
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

        /*! \page aimOutputsAFLR4
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
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int i, status = CAPS_SUCCESS;
    aimStorage *aflr4Instance;
    aimMesh    mesh;

#ifdef DEBUG
    printf(" aflr4AIM/aimCalcOutput  index = %d!\n", index);
#endif
    aflr4Instance = (aimStorage *) instStore;

    if (Done == index) {

      if (aflr4Instance->numNodeTotal > 0 && aflr4Instance->numElemTotal > 0)
        val->vals.integer = (int) true;
      else
        val->vals.integer = (int) false;

    } else if (NumberOfElement == index) {

      val->vals.integer = aflr4Instance->numElemTotal;

    } else if (NumberOfNode == index) {

      val->vals.integer = aflr4Instance->numNodeTotal;

    } else if (Surface_Mesh == index) {

      for (i = 0; i < aflr4Instance->numMeshRef; i++) {
        status = aim_queryMeshes( aimInfo, Surface_Mesh, ANALYSISOUT, &aflr4Instance->meshRef[i] );
        if (status > 0) {
/*@-immediatetrans@*/
          mesh.meshData = NULL;
          mesh.meshRef = &aflr4Instance->meshRef[i];
/*@+immediatetrans@*/

          status = mesh_surfaceMeshData(aimInfo, &aflr4Instance->groupMap, &mesh);
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
      val->nrow        = aflr4Instance->numMeshRef;
      val->vals.AIMptr = aflr4Instance->meshRef;

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
    aimStorage *aflr4Instance;

#ifdef DEBUG
    printf(" aflr4AIM/aimCleanup!\n");
#endif
    aflr4Instance = (aimStorage *) instStore;

    status = destroy_aimStorage(aflr4Instance, (int)false);
    if (status != CAPS_SUCCESS)
          printf(" Status = %d, aflr4AIM aimStorage cleanup!!!\n", status);

    EG_free(aflr4Instance);
}
