// AFLR4 interface functions - Modified from functions provided with
// AFLR4_LIB source written by David L. Marcum

//  Modified by Dr. Ryan Durscher AFRL/RQVC
//  Modified by Dr. Marshall Galbraith MIT

#ifdef WIN32
#define strcasecmp stricmp
#define strtok_r   strtok_s
typedef int         pid_t;
#endif

#include <aflr4/AFLR4_LIB.h> // Bring in AFLR4 API library
#include <egads_aflr4/EGADS_AFLR4_LIB_INC.h>

#include "meshUtils.h" // Collection of helper functions for meshing
#include "miscUtils.h"
#include "aimUtil.h"

#include "aflr4_Interface.h" // Bring in AFLR4 'interface' functions

#ifndef S_SPLINT_S
#define AFLR_STATUS(aimInfo, statys, ...) \
if (status != 0) { status = CAPS_EXECERR; AIM_STATUS(aimInfo, status, ##__VA_ARGS__); }
#else
extern void AFLR_STATUS(void *aimInfo, int status, ...);
#endif


int aflr4_Surface_Mesh(void *aimInfo,
                       int quiet,
                       int numBody, ego *bodies,
                       capsValue *aimInputs,
                       meshInputStruct meshInput)
{
    int status; // Function return status

    int numFace = 0; // Number of faces on a body

    int faceIndex = 0, bodyIndex = 0;

    const char *pstring = NULL;
    const int *pints = NULL;
    const double *preals = NULL;
    int atype, n;

    double capsMeshLength = 0, meshLenFac = 0;

    ego *faces = NULL, *copy_bodies=NULL, context, tess, model=NULL;

    int min_ncell, mer_all, no_prox, quad, skin;
    int EGADSQuad = (int) false;

    double ff_cdfr, abs_min_scale, BL_thickness, Re_l, curv_factor,
           max_scale, min_scale, ref_len, sf_global, erw_all;

    const char *aflr4_debug = "aflr4_debug.egads";
    char aimFile[PATH_MAX];
    char bodyNumber[42];
    const char *attrtype = "";

     // Commandline inputs
    INT_ mmsg = 0;
    int  prog_argc   = 1;    // Number of arguments
    char **prog_argv = NULL; // String arrays
    char *meshInputString = NULL;

    UG_Param_Struct *AFLR4_Param_Struct_Ptr = NULL;

    // return values from egads_aflr4_tess
    ego *tessBodies = NULL;
    void *ext_cad_data = NULL;
    egads_struct *ptr = NULL;


    // Do some additional sanity checks on the attributes
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {

        status = EG_getBodyTopos (bodies[bodyIndex], NULL, FACE, &numFace,
                                  &faces);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(faces, aimInfo, status);

        for (faceIndex = 0; faceIndex < numFace ; faceIndex++) {

            // AFLR Grid BC values should be set for each CAD surface
            // check to make sure the values are not bad

            // check to see if AFLR_GBC is already set
            status = EG_attributeRet(faces[faceIndex], "AFLR_GBC", &atype, &n,
                                     &pints, &preals, &pstring);

            if (status == EGADS_SUCCESS) {
                if (atype == ATTRSTRING) {
/*@-nullpass@*/
                    if ( !( (strcasecmp (pstring, "STD_UG3_GBC") == 0)            ||
                            (strcasecmp (pstring, "-STD_UG3_GBC") == 0)           ||
                            (strcasecmp (pstring, "FARFIELD_UG3_GBC") == 0)       ||
                            (strcasecmp (pstring, "BL_INT_UG3_GBC") == 0)         ||
                            (strcasecmp (pstring, "TRANSP_SRC_UG3_GBC") == 0)     ||
                            (strcasecmp (pstring, "TRANSP_BL_INT_UG3_GBC") == 0)  ||
                            (strcasecmp (pstring, "TRANSP_UG3_GBC") == 0)         ||
                            (strcasecmp (pstring, "-TRANSP_UG3_GBC") == 0)        ||
                            (strcasecmp (pstring, "TRANSP_INTRNL_UG3_GBC") == 0)  ||
                            (strcasecmp (pstring, "-TRANSP_INTRNL_UG3_GBC") == 0) ) ) {
                        AIM_ERROR(  aimInfo, "Invalid AFLR_GBC on face %d of body %d: \"%s\"", faceIndex+1, bodyIndex+1, pstring);
                        AIM_ADDLINE(aimInfo, "Valid string values are:\n");
                        AIM_ADDLINE(aimInfo, "  FARFIELD_UG3_GBC       : farfield surface");
                        AIM_ADDLINE(aimInfo, "  STD_UG3_GBC            : standard surface");
                        AIM_ADDLINE(aimInfo, "  -STD_UG3_GBC           : standard surface");
                        AIM_ADDLINE(aimInfo, "                           BL generating surface");
                        AIM_ADDLINE(aimInfo, "  BL_INT_UG3_GBC         : symmetry or standard surface that intersects BL");
                        AIM_ADDLINE(aimInfo, "  TRANSP_SRC_UG3_GBC     : embedded/transparent surface");
                        AIM_ADDLINE(aimInfo, "                           converted to source nodes");
                        AIM_ADDLINE(aimInfo, "  TRANSP_BL_INT_UG3_GBC  : embedded/transparent surface that intersects BL");
                        AIM_ADDLINE(aimInfo, "  TRANSP_UG3_GBC         : embedded/transparent surface");
                        AIM_ADDLINE(aimInfo, "  -TRANSP_UG3_GBC        : embedded/transparent surface");
                        AIM_ADDLINE(aimInfo, "                           BL generating surface");
                        AIM_ADDLINE(aimInfo, "  TRANSP_INTRNL_UG3_GBC  : embedded/transparent surface");
                        AIM_ADDLINE(aimInfo, "                           converted to an internal surface");
                        AIM_ADDLINE(aimInfo, "                           coordinates are retained but connectivity is not");
                        AIM_ADDLINE(aimInfo, "  -TRANSP_INTRNL_UG3_GBC : embedded/transparent surface");
                        AIM_ADDLINE(aimInfo, "                           converted to an internal surface");
                        AIM_ADDLINE(aimInfo, "                           coordinates are retained but connectivity is not");
                        AIM_ADDLINE(aimInfo, "                           BL generating surface");
/*@+nullpass@*/
                        status = CAPS_BADVALUE;
                        goto cleanup;
                    }
                } else {
                    AIM_ERROR(  aimInfo, "AFLR_GBC on face %d of body %d has %d entries ",
                           faceIndex+1, bodyIndex+1, n);
                    if (atype == ATTRREAL)     { AIM_ADDLINE(aimInfo, "of reals");    }
                    else if (atype == ATTRINT) { AIM_ADDLINE(aimInfo, "of integers"); }
                    AIM_ADDLINE(aimInfo, "Should only contain a string!");
                    status = CAPS_BADVALUE;
                    goto cleanup;
                }
            }

            // check for deprecation
            status = EG_attributeRet(faces[faceIndex], "AFLR_Cmp_ID", &atype, &n,
                                     &pints, &preals, &pstring);
            if (status == EGADS_SUCCESS) {
              AIM_ERROR(aimInfo, "AFLR_Cmp_ID on face %d of body %d is deprecated",
                     faceIndex+1, bodyIndex+1);
              AIM_ADDLINE(aimInfo, "   use AFLR4_Cmp_ID instead!");
              status = CAPS_BADVALUE;
              goto cleanup;
            }

            // check to see if AFLR4_Cmp_ID is set correctly
            status = EG_attributeRet(faces[faceIndex], "AFLR4_Cmp_ID", &atype, &n,
                                     &pints, &preals, &pstring);

            if (status == EGADS_SUCCESS && (!(atype == ATTRREAL || atype == ATTRINT) || n != 1)) {
                //make sure it is only a single real
              if (atype == ATTRREAL)        { attrtype = "of reals";    }
              else if (atype == ATTRINT)    { attrtype = "of integers"; }
              else if (atype == ATTRSTRING) { attrtype = "of a string"; }
                AIM_ERROR(aimInfo, "AFLR4_Cmp_ID on face %d of body %d has %d entries %s",
                       faceIndex+1, bodyIndex+1, n, attrtype);
                AIM_ADDLINE(aimInfo, "Should only contain a single integer or real!");
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            // check to see if AFLR4_Isolated_Edge_Refinement_Flag is set correctly
            status = EG_attributeRet(faces[faceIndex], "AFLR4_Isolated_Edge_Refinement_Flag",
                                     &atype, &n, &pints, &preals, &pstring);

            if (status == EGADS_SUCCESS && (!(atype == ATTRREAL || atype == ATTRINT) || n != 1)) {
                //make sure it is only a single real
                if (atype == ATTRREAL)        { attrtype = "of reals";    }
                else if (atype == ATTRINT)    { attrtype = "of integers"; }
                else if (atype == ATTRSTRING) { attrtype = "of a string"; }
                AIM_ERROR(aimInfo, "AFLR4_Isolated_Edge_Refinement_Flag on face %d of body %d has %d entries %s",
                       faceIndex+1, bodyIndex+1, n, attrtype);
                AIM_ADDLINE(aimInfo, "Should only contain a single integer or real!");
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            // check for deprecation
            status = EG_attributeRet(faces[faceIndex], "AFLR_Scale_Factor", &atype,
                                     &n, &pints, &preals, &pstring);
            if (status == EGADS_SUCCESS) {
              AIM_ERROR(aimInfo, "AFLR_Scale_Factor on face %d of body %d is deprecated",
                     faceIndex+1, bodyIndex+1);
              AIM_ADDLINE(aimInfo, "   use AFLR4_Scale_Factor instead!");
              status = CAPS_BADVALUE;
              goto cleanup;
            }

            // check to see if AFLR4_Scale_Factor is already set
            status = EG_attributeRet(faces[faceIndex], "AFLR4_Scale_Factor", &atype,
                                     &n, &pints, &preals, &pstring);

            if (status == EGADS_SUCCESS && !(atype == ATTRREAL && n == 1)) {
                //make sure it is only a single real
                if (atype == ATTRREAL)        { attrtype = "of reals";    }
                else if (atype == ATTRINT)    { attrtype = "of integers"; }
                else if (atype == ATTRSTRING) { attrtype = "of a string"; }
                AIM_ERROR(aimInfo, "AFLR4_Scale_Factor on face %d of body %d has %d entries %s",
                       faceIndex+1, bodyIndex+1, n, attrtype);
                AIM_ADDLINE(aimInfo, "Should only contain a single real!");
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            // edge mesh spacing scale factor weight for each CAD definition
            // no need to set if scale factor weight is 0.0 (default)

            // check for deprecation
            status = EG_attributeRet(faces[faceIndex], "AFLR_Edge_Scale_Factor_Weight",
                                     &atype, &n, &pints, &preals, &pstring);
            if (status == EGADS_SUCCESS) {
              AIM_ERROR(aimInfo, "AFLR_Edge_Scale_Factor_Weight on face %d of body %d is deprecated",
                     faceIndex+1, bodyIndex+1);
              AIM_ADDLINE(aimInfo, "   use AFLR4_Edge_Refinement_Weight instead!");
              status = CAPS_BADVALUE;
              goto cleanup;
            }

            status = EG_attributeRet(faces[faceIndex], "AFLR4_Edge_Scale_Factor_Weight",
                                     &atype, &n, &pints, &preals, &pstring);
            if (status == EGADS_SUCCESS) {
              AIM_ERROR(aimInfo, "AFLR4_Edge_Scale_Factor_Weight on face %d of body %d is deprecated",
                     faceIndex+1, bodyIndex+1);
              AIM_ADDLINE(aimInfo, "   use AFLR4_Edge_Refinement_Weight instead!");
              status = CAPS_BADVALUE;
              goto cleanup;
            }

            // check to see if AFLR4_Edge_Refinement_Weight is already set
            status = EG_attributeRet(faces[faceIndex], "AFLR4_Edge_Refinement_Weight",
                                     &atype, &n, &pints, &preals, &pstring);

            if (status == EGADS_SUCCESS && !(atype == ATTRREAL && n == 1) ) {
                //make sure it is only a single real
                if (atype == ATTRREAL)        { attrtype = "of reals";    }
                else if (atype == ATTRINT)    { attrtype = "of integers"; }
                else if (atype == ATTRSTRING) { attrtype = "of a string"; }
                AIM_ERROR(aimInfo, "AFLR4_Edge_Refinement_Weight on face %d of body %d has %d entries %s",
                       faceIndex+1, bodyIndex+1, n, attrtype);
                AIM_ADDLINE(aimInfo, "Should only contain a single real!");
                status = CAPS_BADVALUE;
                goto cleanup;
            }
        }

        EG_free(faces); faces = NULL;
    }

    // get AFLR parameters from user inputs

    ff_cdfr       = aimInputs[Ff_cdfr-1].vals.real;
    min_ncell     = aimInputs[Min_ncell-1].vals.integer;
    mer_all       = aimInputs[Mer_all-1].vals.integer;
    no_prox       = aimInputs[No_prox-1].vals.integer;
    quad          = aimInputs[AFLR4_Quad-1].vals.integer;
    skin          = aimInputs[Skin-1].vals.integer;

    BL_thickness  = aimInputs[BL_Thickness-1].vals.real;
    Re_l          = aimInputs[RE_l-1].vals.real;
    curv_factor   = aimInputs[Curv_factor-1].vals.real;
    abs_min_scale = aimInputs[Abs_min_scale-1].vals.real;
    max_scale     = aimInputs[Max_scale-1].vals.real;
    min_scale     = aimInputs[Min_scale-1].vals.real;
    erw_all       = aimInputs[Erw_all-1].vals.real;
    meshLenFac    = aimInputs[Mesh_Length_Factor-1].vals.real;

    EGADSQuad     = aimInputs[EGADS_Quad-1].vals.integer;

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

    if (meshLenFac <= 0) {
      AIM_ERROR(aimInfo, "Mesh_Length_Factor is: %f\n", meshLenFac);
      AIM_ADDLINE(aimInfo, "Mesh_Length_Factor must be a positive number.\n");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // compute the reference length used by AFLR4
    ref_len = capsMeshLength;
    sf_global = meshLenFac;

    /*
    printf("ff_cdfr               = %f\n", ff_cdfr       );
    printf("min_ncell             = %d\n", min_ncell     );
    printf("mer_all               = %d\n", mer_all       );
    printf("no_prox               = %d\n", no_prox       );

    printf("BL_thickness          = %f\n", BL_thickness  );
    printf("Re_l                  = %f\n", Re_l          );
    printf("curv_factor           = %f\n", curv_factor   );
    printf("abs_min_scale         = %f\n", abs_min_scale );
    printf("max_scale             = %f\n", max_scale     );
    printf("min_scale             = %f\n", min_scale     );
    printf("ref_len               = %f\n", ref_len       );
    printf("erw_all               = %f\n", erw_all       );
    */

    // Allocate argument vector.

    status = ug_add_new_arg (&prog_argv, "allocate_and_initialize_argv");
    AFLR_STATUS(aimInfo, status);

    // set AFLR4 input parameters

    status = ug_add_flag_arg (  "min_ncell", &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_int_arg  (   min_ncell , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg (  "mer_all"  , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_int_arg  (   mer_all   , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    if (no_prox == True) {
      status = ug_add_flag_arg ("-no_prox"    , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    }
    if (quad == True) {
      status = ug_add_flag_arg ("-quad"       , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);

      // proximity currently does not work with -quad. Dave says it's a long story
      status = ug_add_flag_arg ("-no_prox"    , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    }
    if (skin == True) {
      status = ug_add_flag_arg ("-skin"       , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    }
    status = ug_add_flag_arg ( "ff_cdfr"      , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (ff_cdfr       , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "BL_thickness" , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (BL_thickness  , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "Re_l"         , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (Re_l          , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "curv_factor"  , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (curv_factor   , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "abs_min_scale", &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (abs_min_scale , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "max_scale"    , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (max_scale     , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "min_scale"    , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (min_scale     , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "ref_len"      , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (ref_len       , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "erw_all"      , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (erw_all       , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_flag_arg ( "sf_global"    , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
    status = ug_add_double_arg (sf_global     , &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);

    // Add meshInputString arguments (if any) to argument vector.
    // Note that since this comes after setting aimInputs the
    // meshInputString could override the value set by aimInputs.
    // Move this to before setting aimInputs arguments to reverse
    // that behavior.

    if (meshInput.aflr4Input.meshInputString != NULL) {
        meshInputString = EG_strdup(meshInput.aflr4Input.meshInputString);
        AIM_NOTNULL(meshInputString, aimInfo, status);

        status = ug_add_list_arg (meshInputString, &prog_argc, &prog_argv);
        AFLR_STATUS(aimInfo, status);
    }

    // Set AFLR4 case name.
    // Used for any requested output files.

    if (aimInputs[Proj_Name-1].nullVal != IsNull)
        ug_set_case_name (aimInputs[Proj_Name-1].vals.string);
    else
        ug_set_case_name ("AFLR4");

    // Set quiet message flag.

    if (quiet == (int)true ) {
        status = ug_add_flag_arg ("mmsg", &prog_argc, &prog_argv); AFLR_STATUS(aimInfo, status);
        status = ug_add_int_arg (0, &prog_argc, &prog_argv);       AFLR_STATUS(aimInfo, status);
    }

    // Register AFLR4-EGADS routines for CAD related setup & cleanup,
    // cad evaluation, cad bounds and generating boundary edge grids.

    // these calls are in aflr4_main_register - if that changes then these
    // need to change these
    aflr4_register_auto_cad_geom_setup (egads_auto_cad_geom_setup);
    aflr4_register_cad_geom_create_tess (egads_aflr4_create_tess);
    aflr4_register_cad_geom_data_cleanup (egads_cad_geom_data_cleanup);
    aflr4_register_cad_geom_file_read (egads_cad_geom_file_read);
    aflr4_register_cad_geom_file_write (egads_cad_geom_file_write);
    aflr4_register_cad_geom_reset_attr (egads_cad_geom_reset_attr);
    aflr4_register_cad_geom_setup (egads_cad_geom_setup);
    aflr4_register_cad_tess_to_dgeom (egads_aflr4_tess_to_dgeom);
    aflr4_register_get_mzone_bedge_data (egads_get_mzone_bedge_data);
    aflr4_register_merge_mzone_data (egads_merge_mzone_data);
    aflr4_register_set_ext_cad_data (egads_set_ext_cad_data);
    aflr4_register_set_mzone_data (egads_set_mzone_data);
    aflr4_register_set_mzone_edge_data (egads_set_mzone_edge_data);

    dgeom_register_cad_eval_curv_at_uv (egads_eval_curv_at_uv);
    dgeom_register_cad_eval_edge_arclen (egads_eval_edge_arclen);
    dgeom_register_cad_eval_uv_bounds (egads_eval_uv_bounds);
    dgeom_register_cad_eval_uv_at_t (egads_eval_uv_at_t);
    dgeom_register_cad_eval_uv_at_xyz (egads_eval_uv_at_xyz);
    dgeom_register_cad_eval_xyz_at_t (egads_eval_xyz_at_u);
    dgeom_register_cad_eval_xyz_at_uv (egads_eval_xyz_at_uv);
    dgeom_register_discrete_eval_xyz_at_t (surfgen_discrete_eval_xyz_at_t);

    // Register fork routines for parallel processing.
    // This is only required to use parallel processing in fork/shared memory
    // mode.

#ifndef WIN32
/*@-unrecog@*/
    ug_mp_register_fork_function (fork);
    ug_mp_register_mmap_function (mmap);
    ug_mp_register_pipe_function (pipe);
/*@+unrecog@*/
#endif

    // Malloc, initialize, and setup AFLR4 input parameter structure.
/*@-nullpass@*/
    status = aflr4_setup_param (mmsg, 0, prog_argc, prog_argv,
                                &AFLR4_Param_Struct_Ptr);
/*@+nullpass@*/
    AFLR_STATUS(aimInfo, status, "aflr4_setup_param failed!");

    // Allocate AFLR4-EGADS data structure, initialize, and link body data.

    AIM_ALLOC(copy_bodies, numBody, ego, aimInfo, status);
    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {
        status = EG_copyObject(bodies[bodyIndex], NULL, &copy_bodies[bodyIndex]);
        AIM_STATUS(aimInfo, status);
    }
    status = EG_getContext(bodies[0], &context);
    AIM_STATUS(aimInfo, status);
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, numBody,
                             copy_bodies, NULL, &model);
    AIM_STATUS(aimInfo, status);

    // Set CAD geometry data structure.
    // Note that memory management of the CAD geometry data structure is
    // controlled by DGEOM after this call.

    // required for setting input data
/*@-nullpass@*/
    (void) ug_set_int_param ("geom_type", 1, AFLR4_Param_Struct_Ptr);
/*@+nullpass@*/

    status = aflr4_set_ext_cad_data (&model);
    AFLR_STATUS(aimInfo, status);

    // Complete all tasks required for AFLR4 surface grid generation.

/*@-nullpass@*/
    status = aflr4_setup_and_grid_gen (1, AFLR4_Param_Struct_Ptr);
    /*@+nullpass@*/
    AFLR_STATUS(aimInfo, status);

//#define DUMP_TECPLOT_DEBUG_FILE
#ifdef DUMP_TECPLOT_DEBUG_FILE
    {
      int surf = 0;
      int numEdge = 0;
      int numSurface = 0;
      int numTriFace = 0;
      int numNodes = 0;
      int numQuadFace = 0;
      int glueId; // Id of glued composite surface

      // AFRL4 output arrays
      INT_1D *bcFlag = NULL;
      INT_1D *ieFlag = NULL;
      INT_1D *idFlag = NULL;
      INT_2D *edgeCon = NULL;
      INT_3D *triCon = NULL;
      INT_4D *quadCon = NULL;
      DOUBLE_2D *uv = NULL;
      DOUBLE_3D *xyz = NULL;

      // Get output id index (glue-only composite)
      dgeom_def_get_idef (0, &glueId);

      FILE *fp = aim_fopen(aimInfo, "aflr4_debug.tec", "w");
      fprintf(fp, "VARIABLES = X, Y, Z, u, v\n");

      numSurface = dgeom_get_ndef(); // Get number of surfaces meshed

      for (surf = 0; surf < numSurface ; surf++) {

        if (surf+1 == glueId) continue;

        status = aflr4_get_def (surf+1,
                                0,
                                &numEdge,
                                &numTriFace,
                                &numNodes,
                                &numQuadFace,
                                &bcFlag,
                                &ieFlag,
                                &idFlag,
                                &edgeCon,
                                &triCon,
                                &quadCon,
                                &uv,
                                &xyz);
        AFLR_STATUS(aimInfo, status);

        fprintf(fp, "ZONE T=\"def %d\" N=%d, E=%d, F=FEPOINT, ET=Quadrilateral, DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n",
                surf+1, numNodes, numTriFace+numQuadFace);
        for (int i = 0; i < numNodes; i++)
          fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e\n",
                  xyz[i+1][0], xyz[i+1][1], xyz[i+1][2], uv[i+1][0], uv[i+1][1]);
        for (int i = 0; i < numTriFace; i++)
          fprintf(fp, "%d %d %d %d\n", triCon[i+1][0], triCon[i+1][1], triCon[i+1][2], triCon[i+1][2]);
        for (int i = 0; i < numQuadFace; i++)
          fprintf(fp, "%d %d %d %d\n", quadCon[i+1][0], quadCon[i+1][1], quadCon[i+1][2], quadCon[i+1][3]);

        ug_free (bcFlag);  bcFlag = NULL;
        ug_free (ieFlag);  ieFlag = NULL;
        ug_free (idFlag);  idFlag = NULL;
        ug_free (edgeCon); edgeCon = NULL;
        ug_free (triCon);  triCon = NULL;
        ug_free (quadCon); quadCon = NULL;
        ug_free (uv);      uv = NULL;
        ug_free (xyz);     xyz = NULL;

      }
      fclose(fp); fp = NULL;

    }
#endif

    status = aflr4_cad_geom_create_tess(mmsg, 1, 0);
    if (status != EGADS_SUCCESS) goto cleanup;

    ext_cad_data = dgeom_get_ext_cad_data ();
    ptr = (egads_struct *) ext_cad_data;
    tessBodies = ptr->tess;

    if (tessBodies == NULL) {
        AIM_ERROR(aimInfo, "aflr4 did not produce EGADS tessellations");
        status = CAPS_NULLOBJ;
        goto cleanup;
    }
    if (ptr->ntess != numBody) {
        AIM_ERROR(aimInfo, "aflr4 body count missmatch!");
        status = CAPS_NULLOBJ;
        goto cleanup;
    }

    for (bodyIndex = 0; bodyIndex < numBody; bodyIndex++) {

        /* apply EGADS quading if requested */
        if (EGADSQuad == (int) true) {
            if (quiet == (int) false)
                printf("Creating EGADS quads for Body %d\n", bodyIndex+1);
            tess = tessBodies[bodyIndex];
            status = EG_quadTess(tess, &tessBodies[bodyIndex]);
            if (status < EGADS_SUCCESS) {
                printf(" EG_quadTess = %d  -- reverting...\n", status);
                tessBodies[bodyIndex] = tess;
            }
        }

        // set the file name to write the egads file
        snprintf(bodyNumber, 42, AFLR4TESSFILE, bodyIndex);
        status = aim_file(aimInfo, bodyNumber, aimFile);
        AIM_STATUS(aimInfo, status);

        remove(aimFile);
        status = EG_saveTess(tessBodies[bodyIndex], aimFile);
        AIM_STATUS(aimInfo, status);
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != 0 && dgeom_get_ext_cad_data () != NULL) {
        AIM_ERROR  (aimInfo, "AFLR4 mesh generation failed...");
        AIM_ADDLINE(aimInfo, "An EGADS file with all AFLR4 parameters");
        AIM_ADDLINE(aimInfo, "has been written to '%s'", aflr4_debug);

        remove(aflr4_debug);
    /*@-mustfreefresh@*/
        (void) EG_saveModel(egads_get_model (0), aflr4_debug);
    /*@+mustfreefresh@*/
    }

    if (ptr != NULL) {
        AIM_FREE (ptr->bodies);
        EG_deleteObject (ptr->model);
    }

    // Free all aflr4 data
    aflr4_free_all (0);

    // Free program arguements
/*@-nullpass*/
    ug_free_argv(prog_argv); prog_argv = NULL;
    ug_free_param (AFLR4_Param_Struct_Ptr);

    AIM_FREE(faces);

    AIM_FREE(meshInputString);

    AIM_FREE(copy_bodies);
    EG_deleteObject(model);
/*@+nullpass@*/

    return status;
}
