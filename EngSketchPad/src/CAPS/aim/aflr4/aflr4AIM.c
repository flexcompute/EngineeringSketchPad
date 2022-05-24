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
 * "Advancing-Front/Local-Reconnection (AFLR) procedure." Only triangular elements may be generated, with planned future support of quadrilateral elements.
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
 * - <b> AFLR4_Scale_Factor</b> [Optional FACE attribute: Default 1.0]<br>
 * EGADS attribute AFLR4_Scale_Factor represents the AFLR4 surface mesh spacing
 * scale factor for a given face/surface. Curvature dependent spacing can be
 * scaled on the face/surface by the value of the scale factor set with
 * AFLR4_Scale_Factor.<br>
 * <br>
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
    int numSurface;
    meshStruct *surfaceMesh;

    // Container for mesh input
    meshInputStruct meshInput;

    // Attribute to index map
    mapAttrToIndexStruct groupMap;

    mapAttrToIndexStruct meshMap;

} aimStorage;


static int destroy_aimStorage(aimStorage *aflr4Instance)
{

    int i; // Indexing

    int status; // Function return status

    // Destroy meshInput
    status = destroy_meshInputStruct(&aflr4Instance->meshInput);
    if (status != CAPS_SUCCESS)
        printf("Status = %d, aflr4AIM meshInput cleanup!!!\n", status);

    // Destroy surface mesh allocated arrays
    for (i = 0; i < aflr4Instance->numSurface; i++) {

        status = destroy_meshStruct(&aflr4Instance->surfaceMesh[i]);
        if (status != CAPS_SUCCESS)
            printf("Status = %d, aflr4AIM surfaceMesh cleanup!!!\n", status);

    }
    aflr4Instance->numSurface = 0;

    if (aflr4Instance->surfaceMesh != NULL) EG_free(aflr4Instance->surfaceMesh);
    aflr4Instance->surfaceMesh = NULL;

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

    status = EG_getBodyTopos(body, NULL, FACE, &numEdge, &edges);
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
                    else if (strncasecmp(meshProp[propIndex].bcType, "Symmetry"  ,  8) == 0)
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
    aflr4Instance->surfaceMesh = NULL;
    aflr4Instance->numSurface = 0;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&aflr4Instance->meshMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_mapAttrToIndexStruct(&aflr4Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    // Container for mesh input
    status = initiate_meshInputStruct(&aflr4Instance->meshInput);
    AIM_STATUS(aimInfo, status);

cleanup:
    if (status != CAPS_SUCCESS) AIM_FREE(*instStore);
    return status;
}


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
        *ainame              = AIM_NAME(Proj_Name); // If NULL a volume grid won't be written by the AIM
        defval->type         = String;
        defval->nullVal      = IsNull;
        defval->vals.string  = NULL;
        //defval->vals.string  = EG_strdup("CAPS");
        defval->lfixed       = Change;

        /*! \page aimInputsAFLR4
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

    } else if (index == Mesh_Format) {
        *ainame               = AIM_NAME(Mesh_Format);
        defval->type          = String;
        defval->vals.string   = EG_strdup("AFLR3"); // TECPLOT, VTK, AFLR3, STL, FAST
        defval->lfixed        = Change;

        /*! \page aimInputsAFLR4
         * - <B> Mesh_Format = "AFLR3"</B> <br>
         * Mesh output format. Available format names include: "AFLR3", "VTK", "TECPLOT", "STL" (quadrilaterals will be
         * split into triangles), "FAST", "ETO".
         */

    } else if (index == Mesh_ASCII_Flag) {
        *ainame               = AIM_NAME(Mesh_ASCII_Flag);
        defval->type          = Boolean;
        defval->vals.integer  = true;

        /*! \page aimInputsAFLR4
         * - <B> Mesh_ASCII_Flag = True</B> <br>
         * Output mesh in ASCII format, otherwise write a binary file if applicable.
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
       * the user must always set it via capsMeshLength andd Mesh_Length_Factor
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
         * Scaling factor to compute AFLR4 'ref_len' parameter via:<br>
         * <br>
         * ref_len = capsMeshLength * Mesh_Length_Factor<br>
         * <br>
         * where capsMeshLength is a numeric attribute that<br>
         * must be on at least one body and consistent if on multiple bodies.<br>
         * <br>
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
        defval->type          = Boolean;
        defval->vals.integer  = (int) true;

        /*! \page aimInputsAFLR4
         * - <B> Multiple_Mesh = True</B> <br>
         * If set to True (default) a surface mesh will be generated and output (given "Proj_Name" is set) for each body.
         * When set to False only a single surface
         * mesh will be created. Note, this only affects the mesh when writing to a
         * file.
         */

    } else if (index == EGADS_Quad) {
        *ainame              = EG_strdup("EGADS_Quad");
        defval->type         = Boolean;
        defval->vals.integer  = false;

        /*! \page aimInputsAFLR4
         * - <B> EGADS_Quad = False </B> <br>
         * Apply EGADS quadding to the AFLR4 triangulation.
         */

    } else {
        status = CAPS_BADINDEX;
        AIM_STATUS(aimInfo, status, "Unknown input index %d!", index);
    }

    status = CAPS_SUCCESS;

cleanup:
    ug_free_param(AFLR4_Param_Struct_Ptr);

    if (status != CAPS_SUCCESS) {
        printf("An error occurred creating aimInputs\n");
    }
    return status;
}


int aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{
    int status; // Status return

    int i, bodyIndex; // Indexing
    aimStorage *aflr4Instance;

    // Body parameters
    const char *intents;
    int numBody = 0; // Number of bodies
    ego *bodies = NULL; // EGADS body objects

    // File output
    char *filename = NULL;
    char bodyNumber[11];

    // Mesh attribute parameters
    int numMeshProp = 0;
    meshSizingStruct *meshProp = NULL;

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
    if (aimInputs == NULL) return CAPS_NULLVALUE;

    aflr4Instance = (aimStorage *) instStore;

    // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
    status = destroy_aimStorage(aflr4Instance);
    AIM_STATUS(aimInfo, status);

    // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(numBody,
                                            bodies,
                                            3, // Node level
                                            &aflr4Instance->groupMap);
    AIM_STATUS(aimInfo, status);

    status = create_CAPSMeshAttrToIndexMap(numBody,
                                           bodies,
                                           3, // Node level
                                           &aflr4Instance->meshMap);
    AIM_STATUS(aimInfo, status);

    // Setup meshing input structure

    // Get Tessellation parameters -Tess_Params
    aflr4Instance->meshInput.paramTess[0] = 0;
    aflr4Instance->meshInput.paramTess[1] = 0;
    aflr4Instance->meshInput.paramTess[2] = 0;

    aflr4Instance->meshInput.quiet            = aimInputs[Mesh_Quiet_Flag-1].vals.integer;
    aflr4Instance->meshInput.outputASCIIFlag  = aimInputs[Mesh_ASCII_Flag-1].vals.integer;

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


    status = aflr4_Surface_Mesh(aimInfo,
                                aflr4Instance->meshInput.quiet,
                                numBody, bodies,
                                aimInputs,
                                aflr4Instance->meshInput);
    AIM_STATUS(aimInfo, status, "Problem during AFLR4 surface meshing");

    if (aflr4Instance->meshInput.outputFileName != NULL) {

        for (bodyIndex = 0; bodyIndex < aflr4Instance->numSurface; bodyIndex++) {

            if (aflr4Instance->numSurface > 1) {
/*@-bufferoverflowhigh@*/
                sprintf(bodyNumber, "%d", bodyIndex);
/*@+bufferoverflowhigh@*/
                filename = (char *) EG_alloc((strlen(aflr4Instance->meshInput.outputFileName) +
                                              strlen("_Surf_") + 2 +
                                              strlen(bodyNumber))*sizeof(char));
            } else {
                filename = (char *) EG_alloc((strlen(aflr4Instance->meshInput.outputFileName) +
                                              2)*sizeof(char));
            }

            if (filename == NULL) {
                status = EGADS_MALLOC;
                goto cleanup;
            }

            strcpy(filename, aflr4Instance->meshInput.outputFileName);

            if (aflr4Instance->numSurface > 1) {
                strcat(filename,"_Surf_");
                strcat(filename, bodyNumber);
            }

            if (strcasecmp(aflr4Instance->meshInput.outputFormat, "AFLR3") == 0) {

                status = mesh_writeAFLR3(aimInfo,
                                         filename,
                                         aflr4Instance->meshInput.outputASCIIFlag,
                                         &aflr4Instance->surfaceMesh[bodyIndex],
                                         1.0);

            } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "VTK") == 0) {

                status = mesh_writeVTK(aimInfo,
                                       filename,
                                       aflr4Instance->meshInput.outputASCIIFlag,
                                       &aflr4Instance->surfaceMesh[bodyIndex],
                                       1.0);

            } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Tecplot") == 0) {

                status = mesh_writeTecplot(aimInfo,
                                           filename,
                                           aflr4Instance->meshInput.outputASCIIFlag,
                                           &aflr4Instance->surfaceMesh[bodyIndex],
                                           1.0);

            } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "STL") == 0) {

                status = mesh_writeSTL(aimInfo,
                                       filename,
                                       aflr4Instance->meshInput.outputASCIIFlag,
                                       &aflr4Instance->surfaceMesh[bodyIndex],
                                       1.0);

            } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "FAST") == 0) {

                status = mesh_writeFAST(aimInfo,
                                        filename,
                                        aflr4Instance->meshInput.outputASCIIFlag,
                                        &aflr4Instance->surfaceMesh[bodyIndex],
                                        1.0);

            } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "ETO") == 0) {

                filename = (char *) EG_reall(filename,(strlen(filename) + 5) *sizeof(char));
                if (filename == NULL) {
                    status = EGADS_MALLOC;
                    goto cleanup;
                }
                strcat(filename,".eto");

                status = EG_saveTess(aflr4Instance->surfaceMesh[bodyIndex].egadsTess, filename);

            } else {
                printf("Unrecognized mesh format, \"%s\", the volume mesh will not be written out\n", aflr4Instance->meshInput.outputFormat);
            }

            if (filename != NULL) EG_free(filename);
            filename = NULL;

            AIM_STATUS(aimInfo, status);
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

    if (status != CAPS_SUCCESS)
      printf("Error: aflr4AIM status %d\n", status);

    if (filename != NULL) EG_free(filename);
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
int aimPostAnalysis( void *aimStore, void *aimInfo,
                     int restart,    capsValue *aimInputs)
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

    aimStorage *aflr4Instance;
    aflr4Instance = (aimStorage *) aimStore;

    // Combined mesh
    meshStruct combineMesh;

    AIM_NOTNULL(aimInputs, aimInfo, status);

    status = initiate_meshStruct(&combineMesh);
    AIM_STATUS(aimInfo, status);

    // Get AIM bodies
    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodies, aimInfo, status);

    // Allocate surfaceMesh from number of bodies
    AIM_ALLOC(aflr4Instance->surfaceMesh, numBody, meshStruct, aimInfo, status);
    aflr4Instance->numSurface = numBody;

    // Initiate surface meshes
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++){
        status = initiate_meshStruct(&aflr4Instance->surfaceMesh[bodyIndex]);
        AIM_STATUS(aimInfo, status);
    }

    if (aflr4Instance->groupMap.mapName == NULL) {
        // Get capsGroup name and index mapping to make sure all faces have a capsGroup value
        status = create_CAPSGroupAttrToIndexMap(numBody,
                                                bodies,
                                                3,
                                                &aflr4Instance->groupMap);
        AIM_STATUS(aimInfo, status);
    }

    // Run egadsTess for each body
    for (bodyIndex = 0 ; bodyIndex < numBody; bodyIndex++) {

        status = copy_mapAttrToIndexStruct( &aflr4Instance->groupMap,
                                            &aflr4Instance->surfaceMesh[bodyIndex].groupMap );
        AIM_STATUS(aimInfo, status);

        // set the file name to read the egads file
        snprintf(bodyNumber, 42, AFLR4TESSFILE, bodyIndex);
        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        status = EG_loadTess(bodies[bodyIndex], aimFile, &aflr4Instance->surfaceMesh[bodyIndex].egadsTess);
        AIM_STATUS(aimInfo, status);

        status = mesh_surfaceMeshEGADSTess(aimInfo, &aflr4Instance->surfaceMesh[bodyIndex]);
        AIM_STATUS(aimInfo, status);

        status = aim_newTess(aimInfo, aflr4Instance->surfaceMesh[bodyIndex].egadsTess);
        AIM_STATUS(aimInfo, status);

        if (restart == 0 &&
            aimInputs[Mesh_Quiet_Flag-1].vals.integer == (int)false) {
            printf("Body %d (of %d)\n", bodyIndex+1, numBody);

            printf("Number of nodes    = %d\n", aflr4Instance->surfaceMesh[bodyIndex].numNode);
            printf("Number of elements = %d\n", aflr4Instance->surfaceMesh[bodyIndex].numElement);

            if (aflr4Instance->surfaceMesh[bodyIndex].meshQuickRef.useStartIndex == (int) true ||
                aflr4Instance->surfaceMesh[bodyIndex].meshQuickRef.useListIndex  == (int) true) {

                printf("Number of node elements          = %d\n",
                       aflr4Instance->surfaceMesh[bodyIndex].meshQuickRef.numNode);
                printf("Number of line elements          = %d\n",
                       aflr4Instance->surfaceMesh[bodyIndex].meshQuickRef.numLine);
                printf("Number of triangle elements      = %d\n",
                       aflr4Instance->surfaceMesh[bodyIndex].meshQuickRef.numTriangle);
                printf("Number of quadrilateral elements = %d\n",
                       aflr4Instance->surfaceMesh[bodyIndex].meshQuickRef.numQuadrilateral);
            }

            numNodeTotal += aflr4Instance->surfaceMesh[bodyIndex].numNode;
            numElemTotal += aflr4Instance->surfaceMesh[bodyIndex].numElement;
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
            AIM_FREE(aflr4Instance->meshInput.outputFileName);
            AIM_STRDUP(aflr4Instance->meshInput.outputFileName, aimInputs[Proj_Name-1].vals.string, aimInfo, status);
        }

        // Mesh Format
        AIM_FREE(aflr4Instance->meshInput.outputFormat);
        AIM_STRDUP(aflr4Instance->meshInput.outputFormat, aimInputs[Mesh_Format-1].vals.string, aimInfo, status);

        // ASCII flag
        aflr4Instance->meshInput.outputASCIIFlag = aimInputs[Mesh_ASCII_Flag-1].vals.integer;

        if (aflr4Instance->meshInput.outputFileName != NULL) {

            // We need to combine the mesh
            if (aimInputs[Multiple_Mesh-1].vals.integer == (int) false) {

                status = mesh_combineMeshStruct(aflr4Instance->numSurface,
                                                aflr4Instance->surfaceMesh,
                                                &combineMesh);

                AIM_STATUS(aimInfo, status);

                filename = (char *) EG_alloc((strlen(aflr4Instance->meshInput.outputFileName) +
                                              2)*sizeof(char));

                if (filename == NULL) goto cleanup;

                strcpy(filename, aflr4Instance->meshInput.outputFileName);

                if (strcasecmp(aflr4Instance->meshInput.outputFormat, "AFLR3") == 0) {

                    status = mesh_writeAFLR3(aimInfo, filename,
                                             aflr4Instance->meshInput.outputASCIIFlag,
                                             &combineMesh,
                                             1.0);

                } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "VTK") == 0) {

                    status = mesh_writeVTK(aimInfo, filename,
                                           aflr4Instance->meshInput.outputASCIIFlag,
                                           &combineMesh,
                                           1.0);

                } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Tecplot") == 0) {

                    status = mesh_writeTecplot(aimInfo, filename,
                                               aflr4Instance->meshInput.outputASCIIFlag,
                                               &combineMesh,
                                               1.0);

                } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "STL") == 0) {

                    status = mesh_writeSTL(aimInfo, filename,
                                           aflr4Instance->meshInput.outputASCIIFlag,
                                           &combineMesh,
                                           1.0);

                } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Airfoil") == 0) {

                    status = mesh_writeAirfoil(aimInfo, filename,
                                               aflr4Instance->meshInput.outputASCIIFlag,
                                               &combineMesh,
                                               1.0);

                } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "FAST") == 0) {

                    status = mesh_writeFAST(aimInfo, filename,
                                            aflr4Instance->meshInput.outputASCIIFlag,
                                            &combineMesh,
                                            1.0);

                } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Nastran") == 0) {

                    status = mesh_writeNASTRAN(aimInfo, filename,
                                               aflr4Instance->meshInput.outputASCIIFlag,
                                               &combineMesh,
                                               FreeField,
                                               1.0);
                } else {
                    printf("Unrecognized mesh format, \"%s\", the mesh will not be written out\n",
                           aflr4Instance->meshInput.outputFormat);
                }

                AIM_FREE(filename);
                AIM_STATUS(aimInfo, status);

            } else {

                for (bodyIndex = 0; bodyIndex < aflr4Instance->numSurface; bodyIndex++) {

                    if (aflr4Instance->numSurface > 1) {
                        snprintf(bodyNumber, 42, "%d", bodyIndex);
                        filename = (char *) EG_alloc((strlen(aflr4Instance->meshInput.outputFileName)  +
                                                      2 + strlen("_Surf_") + strlen(bodyNumber))*sizeof(char));
                    } else {
                        filename = (char *) EG_alloc((strlen(aflr4Instance->meshInput.outputFileName) +
                                                      2)*sizeof(char));

                    }

                    if (filename == NULL) goto cleanup;

                    strcpy(filename, aflr4Instance->meshInput.outputFileName);

                    if (aflr4Instance->numSurface > 1) {
                        strcat(filename,"_Surf_");
                        strcat(filename, bodyNumber);
                    }

                    if (strcasecmp(aflr4Instance->meshInput.outputFormat, "AFLR3") == 0) {

                        status = mesh_writeAFLR3(aimInfo, filename,
                                                 aflr4Instance->meshInput.outputASCIIFlag,
                                                 &aflr4Instance->surfaceMesh[bodyIndex],
                                                 1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "VTK") == 0) {

                        status = mesh_writeVTK(aimInfo, filename,
                                                aflr4Instance->meshInput.outputASCIIFlag,
                                                &aflr4Instance->surfaceMesh[bodyIndex],
                                                1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Tecplot") == 0) {

                        status = mesh_writeTecplot(aimInfo, filename,
                                                   aflr4Instance->meshInput.outputASCIIFlag,
                                                   &aflr4Instance->surfaceMesh[bodyIndex],
                                                   1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "STL") == 0) {

                        status = mesh_writeSTL(aimInfo, filename,
                                               aflr4Instance->meshInput.outputASCIIFlag,
                                               &aflr4Instance->surfaceMesh[bodyIndex],
                                               1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Airfoil") == 0) {

                        status = mesh_writeAirfoil(aimInfo, filename,
                                                   aflr4Instance->meshInput.outputASCIIFlag,
                                                   &aflr4Instance->surfaceMesh[bodyIndex],
                                                   1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "FAST") == 0) {

                        status = mesh_writeFAST(aimInfo, filename,
                                                aflr4Instance->meshInput.outputASCIIFlag,
                                                &aflr4Instance->surfaceMesh[bodyIndex],
                                                1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "Nastran") == 0) {

                        status = mesh_writeNASTRAN(aimInfo, filename,
                                                   aflr4Instance->meshInput.outputASCIIFlag,
                                                   &aflr4Instance->surfaceMesh[bodyIndex],
                                                   FreeField,
                                                   1.0);

                    } else if (strcasecmp(aflr4Instance->meshInput.outputFormat, "ETO") == 0) {

                        filename = (char *) EG_reall(filename,(strlen(filename) + 5) *sizeof(char));
                        if (filename == NULL) {
                            status = EGADS_MALLOC;
                            goto cleanup;
                        }
                        strcat(filename,".eto");

                        status = EG_saveTess(aflr4Instance->surfaceMesh[bodyIndex].egadsTess, filename);

                    } else {
                        printf("Unrecognized mesh format, \"%s\", the mesh will not be written out\n",
                               aflr4Instance->meshInput.outputFormat);
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
        form->type        = Pointer;
        form->dim         = Vector;
        form->lfixed      = Change;
        form->sfixed      = Change;
        form->vals.AIMptr = NULL;
        form->nullVal     = IsNull;
        AIM_STRDUP(form->units, "meshStruct", aimStruc, status);

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


// See if a surface mesh was generated for each body
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
    int status = CAPS_SUCCESS;
    int surf; // Indexing
    int numElement, numNodes, nElem;
    aimStorage *aflr4Instance;

#ifdef DEBUG
    printf(" aflr4AIM/aimCalcOutput  index = %d!\n", index);
#endif
    aflr4Instance = (aimStorage *) instStore;

    if (Done == index) {

        val->vals.integer = (int) false;

        // Check to see if surface meshes was generated
        for (surf = 0; surf < aflr4Instance->numSurface; surf++ ) {

            if (aflr4Instance->surfaceMesh[surf].numElement != 0) {

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
        for (surf = 0; surf < aflr4Instance->numSurface; surf++ ) {

            status = mesh_retrieveNumMeshElements(aflr4Instance->surfaceMesh[surf].numElement,
                                                  aflr4Instance->surfaceMesh[surf].element,
                                                  Triangle,
                                                  &nElem);
            AIM_STATUS(aimInfo, status);
            numElement += nElem;

            status = mesh_retrieveNumMeshElements(aflr4Instance->surfaceMesh[surf].numElement,
                                                  aflr4Instance->surfaceMesh[surf].element,
                                                  Quadrilateral,
                                                  &nElem);
            AIM_STATUS(aimInfo, status);
            numElement += nElem;
        }

        val->vals.integer = numElement;

    } else if (NumberOfNode == index) {

        // Count the total number of surface vertices
        numNodes = 0;
        for (surf = 0; surf < aflr4Instance->numSurface; surf++ ) {
            numNodes += aflr4Instance->surfaceMesh[surf].numNode;
        }

        val->vals.integer = numNodes;

    } else if (Surface_Mesh == index) {

        // Return the surface meshes
        val->nrow        = aflr4Instance->numSurface;
        val->vals.AIMptr = aflr4Instance->surfaceMesh;

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

    status = destroy_aimStorage(aflr4Instance);
    if (status != CAPS_SUCCESS)
          printf(" Status = %d, aflr4AIM aimStorage cleanup!!!\n", status);

    EG_free(aflr4Instance);
}
