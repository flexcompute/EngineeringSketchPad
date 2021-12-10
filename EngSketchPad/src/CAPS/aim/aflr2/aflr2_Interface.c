// AFLR2 interface functions - Modified from functions provided with
//	AFLR2_LIB source (aflr2_cad_geom_setup.c) written by David L. Marcum

//  Modified by Dr. Ryan Durscher AFRL/RQVC

#include "aflr2c/AFLR2_LIB.h" // Bring in AFLR2 API library
#include "ug_gq/UG_GQ_LIB_INC.h"
#include "ug_io/UG_IO_LIB_INC.h"

#include "egads.h"     // Bring in EGADS

#include "aimUtil.h"
#include "aimMesh.h"

#include "meshUtils.h" // Collection of helper functions for meshing
#include "cfdUtils.h"  // Collection of helper functions for cfd analysis
#include "miscUtils.h"

#include "aflr2_Interface.h"

#ifdef WIN32
#define strtok_r   strtok_s
#endif

//TODO: xplt_Case_Name is needed for OSX at the moment to resolve a link error. Remove it as soon as you can...
char xplt_Case_Name[133];

//#define AFLR2_BOUNDARY_LAYER

#define CROSS(a,b,c)      a[0] = ((b)[1]*(c)[2]) - ((b)[2]*(c)[1]);\
                          a[1] = ((b)[2]*(c)[0]) - ((b)[0]*(c)[2]);\
                          a[2] = ((b)[0]*(c)[1]) - ((b)[1]*(c)[0])


static INT_
egads_eval_bedge(void *aimInfo, ego tess,
                 mapAttrToIndexStruct *groupMap, mapAttrToIndexStruct *meshMap,
                 INT_ *nbedge, INT_2D ** inibe, INT_1D **Bnd_Edge_ID_Flag,
                 INT_1D **Bnd_Edge_Mesh_ID_Flag, DOUBLE_2D ** x)
{

  /*
   * Generate a boundary edge grid for a given surface (idef) with a CAD geometry
   * definition. The boundary edge connectivity and coordinate arrays do
   * not need to be allocated prior to calling this routine. They should be set to
   * NULL if they are not allocated. In either case they will be reallocated for
   * the generated boundary edge grid size.
   *
   * Modifications are required for alternative CAD systems.
   *
   * Input Parameters
   *
   *  idef  : identification label for surface definition
   *  inibe : boundary edge node connectivity pointer
   *          must at least be set to NULL it will be reallocated
   *  u     : u,v-coordinates pointer
   *          must at least be set to NULL it will be reallocated
   *  x     : x,y,z-coordinates pointer
   *          must at least be set to NULL it will be reallocated
   *
   * Output Parameters
   *
   *  nbedge : number of boundary edges in boundary edge grid
   *  inibe  : boundary edge node connectivity for each boundary edge of
   *              boundary edge grid
   *  x      : x,y,z-coordinates of boundary edge grid
   *
   * Return Value
   *
   *  0     : no errors
   *  > 0   : there were errors
   *
   */

    INT_         ierr;
    int          cnt, i, iloop, iedge, last, mtype, oclass,
                 n, nloop, nedge, *sen, *senses, status;
    int          groupIndex, meshIndex, index;

    double       range[2], uvbox[4];

    const double *xyzs = NULL, *ts = NULL;
    double *uvs = NULL;

    ego          geom, *loops = NULL, *edges = NULL, *nodes = NULL;

    int tessStatus, numPoints, numEdge, numFace;
    ego body, *bodyFaces = NULL, *bodyEdges = NULL;
    const char *groupName = NULL, *meshName=NULL;

    // Get body from tessellation
    status = EG_statusTessBody(tess, &body, &tessStatus, &numPoints);
    AIM_STATUS(aimInfo, status);

    // Get edges
    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &bodyEdges);
    AIM_STATUS(aimInfo, status);

    // Get faces
    status = EG_getBodyTopos(body, NULL, FACE, &numFace, &bodyFaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodyFaces, aimInfo, status);

    status = EG_getTopology(bodyFaces[0], &geom, &oclass, &mtype, uvbox,
                            &nloop, &loops, &senses);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(loops, aimInfo, status);

    for (cnt = iloop = 0; iloop < nloop; iloop++) {
        status = EG_getTopology(loops[iloop], &geom, &oclass, &mtype,
                                NULL, &nedge, &edges, &senses);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(edges, aimInfo, status);

        for (iedge = 0; iedge < nedge; iedge++) {
            status = EG_getTopology(edges[iedge], &geom, &oclass,
                                    &mtype, range, &n, &nodes, &sen);
            AIM_STATUS(aimInfo, status);
            if (mtype == DEGENERATE) continue;

            index = EG_indexBodyTopo(body, edges[iedge]);
            if (index < EGADS_SUCCESS) {
                status = index;
                AIM_STATUS(aimInfo, status);
            }

            status = EG_getTessEdge(tess, index, &n, &xyzs, &ts);
            AIM_STATUS(aimInfo, status);

            
            cnt += n-1;
        }
    }

    ierr = 0;
    *nbedge = cnt;

    *inibe                 = (INT_2D *) ug_malloc (&ierr, ((*nbedge)+1) * sizeof (INT_2D));
    *Bnd_Edge_ID_Flag      = (INT_1D *) ug_malloc (&ierr, ((*nbedge)+1) * sizeof (INT_1D));
    *Bnd_Edge_Mesh_ID_Flag = (INT_1D *) ug_malloc (&ierr, ((*nbedge)+1) * sizeof (INT_1D));
    *x                  = (DOUBLE_2D *) ug_malloc (&ierr, ((*nbedge)+1) * sizeof (DOUBLE_2D));
    if (ierr) {
        ug_error_message("*** ERROR 104111 : unable to allocate required memory ***");
        return 104111;
    }

    for (cnt = iloop = 0; iloop < nloop; iloop++) {
        status = EG_getTopology(loops[iloop], &geom, &oclass, &mtype,
                                NULL, &nedge, &edges, &senses);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(edges, aimInfo, status);

        last = cnt;
        for (iedge = 0; iedge < nedge; iedge++) {
            status = EG_getTopology(edges[iedge], &geom, &oclass,
                                    &mtype, range, &n, &nodes, &sen);
            AIM_STATUS(aimInfo, status);
            if (mtype == DEGENERATE) continue;

            index = EG_indexBodyTopo(body, edges[iedge]);
            if (index < EGADS_SUCCESS) {
                AIM_STATUS(aimInfo, status);
                goto cleanup;
            }

            status = retrieve_CAPSGroupAttr(edges[iedge], &groupName);
            if (status == EGADS_SUCCESS) {
                AIM_NOTNULL(groupName, aimInfo, status);
                status = get_mapAttrToIndexIndex(groupMap, groupName, &groupIndex);
                if (status != CAPS_SUCCESS) {
                    AIM_ERROR(aimInfo, "No capsGroup \"%s\" not found in attribute map", groupName);
                    goto cleanup;
                }
            } else if (status == EGADS_NOTFOUND) {
                AIM_ERROR(aimInfo, "No capsGroup found on edge %d", index);
                print_AllAttr(aimInfo, edges[iedge] );
                goto cleanup;
            } else goto cleanup;

            status = retrieve_CAPSMeshAttr(edges[iedge], &meshName);
            if (status == EGADS_SUCCESS) {
                AIM_NOTNULL(meshName, aimInfo, status);
                status = get_mapAttrToIndexIndex(meshMap, meshName, &meshIndex);
                if (status != CAPS_SUCCESS) {
                    AIM_ERROR(aimInfo, "No capsMesh \"%s\" not found in attribute map", meshName);
                    goto cleanup;
                }
            } else {
                meshIndex = -1;
            }

            status = EG_getTessEdge(tess, index, &n, &xyzs, &ts);
            AIM_STATUS(aimInfo, status);
            AIM_NOTNULL(ts, aimInfo, status);

            AIM_REALL(uvs, 2*n, double, aimInfo, status);

            status = EG_getEdgeUVs(bodyFaces[0], edges[iedge], senses[iedge], n, ts, uvs);
            AIM_STATUS(aimInfo, status);

            for (i = 0; i < n-1; i++, cnt++) {

                (*Bnd_Edge_ID_Flag)[cnt+1] = groupIndex;
                (*Bnd_Edge_Mesh_ID_Flag)[cnt+1] = meshIndex;

                if (senses[iedge] == 1) {
                    (*x)[cnt+1][0] = uvs[2*i+0];
                    (*x)[cnt+1][1] = uvs[2*i+1];
                } else {
                    (*x)[cnt+1][0] = uvs[2*(n-i-1)+0];
                    (*x)[cnt+1][1] = uvs[2*(n-i-1)+1];
                }

                (*inibe)[cnt+1][0] = cnt+1;
                (*inibe)[cnt+1][1] = cnt+2;
            }
        }
        (*inibe)[cnt][1] = last+1;
    }

    status = EGADS_SUCCESS;

cleanup:

    AIM_FREE(bodyEdges);
    AIM_FREE(bodyFaces);
    AIM_FREE(uvs);

    return status;
}


static int
egads_xyz_bedge(void *aimInfo, ego tess, int ix, int iy,
                double *face_xyz, DOUBLE_2D* Coordiantes)
{
    int          cnt, i, iloop, iedge, mtype, oclass,
                 n, nloop, nedge, *sen, *senses, status;
    int          index;

    double       range[2], uvbox[4];

    const double *xyzs = NULL, *ts = NULL;

    ego          geom, *loops = NULL, *edges = NULL, *nodes = NULL;

    int tessStatus, numPoints, numEdge, numFace;
    ego body, *bodyFaces = NULL, *bodyEdges = NULL;

    // Get body from tessellation
    status = EG_statusTessBody(tess, &body, &tessStatus, &numPoints);
    AIM_STATUS(aimInfo, status);

    // Get edges
    status = EG_getBodyTopos(body, NULL, EDGE, &numEdge, &bodyEdges);
    AIM_STATUS(aimInfo, status);

    // Get faces
    status = EG_getBodyTopos(body, NULL, FACE, &numFace, &bodyFaces);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(bodyFaces, aimInfo, status);

    status = EG_getTopology(bodyFaces[0], &geom, &oclass, &mtype, uvbox,
                            &nloop, &loops, &senses);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(loops, aimInfo, status);

    for (cnt = iloop = 0; iloop < nloop; iloop++) {
        status = EG_getTopology(loops[iloop], &geom, &oclass, &mtype,
                                NULL, &nedge, &edges, &senses);
        AIM_STATUS(aimInfo, status);
        AIM_NOTNULL(edges, aimInfo, status);
        AIM_NOTNULL(senses, aimInfo, status);

        for (iedge = 0; iedge < nedge; iedge++) {
            status = EG_getTopology(edges[iedge], &geom, &oclass,
                                    &mtype, range, &n, &nodes, &sen);
            AIM_STATUS(aimInfo, status);
            if (mtype == DEGENERATE) continue;

            index = EG_indexBodyTopo(body, edges[iedge]);
            if (index < EGADS_SUCCESS) {
                AIM_STATUS(aimInfo, status);
                goto cleanup;
            }

            status = EG_getTessEdge(tess, index, &n, &xyzs, &ts);
            AIM_STATUS(aimInfo, status);
            AIM_NOTNULL(xyzs, aimInfo, status);

            for (i = 0; i < n-1; i++, cnt++) {

                if (senses[iedge] == 1) {
                  face_xyz[3*cnt+0] = xyzs[3*i+0];
                  face_xyz[3*cnt+1] = xyzs[3*i+1];
                  face_xyz[3*cnt+2] = xyzs[3*i+2];
                } else {
                  face_xyz[3*cnt+0] = xyzs[3*(n-i-1)+0];
                  face_xyz[3*cnt+1] = xyzs[3*(n-i-1)+1];
                  face_xyz[3*cnt+2] = xyzs[3*(n-i-1)+2];
                }

                Coordiantes[cnt+1][0] = face_xyz[3*cnt+ix];
                Coordiantes[cnt+1][1] = face_xyz[3*cnt+iy];
            }
        }
    }

    status = EGADS_SUCCESS;

cleanup:

    AIM_FREE(bodyEdges);
    AIM_FREE(bodyFaces);

    //printf("egads_eval_bedge- cleanup complete\n");
    return status;
}


int aflr2_Surface_Mesh(void *aimInfo,
                       int Message_Flag, ego bodyIn,
                       meshInputStruct *meshInput,
                       mapAttrToIndexStruct *groupMap,
                       mapAttrToIndexStruct *meshMap,
                       meshStruct *surfaceMesh,
                       aimMeshRef *meshRef)
{

    int status; // Function return status

    int i, j, elementIndex; // Indexing
    int faceAttr = 0;
    int numFace;

    double box[6], size, params[3]; // Bounding box variables

    ego tess = NULL, body, *bodyFaces = NULL;

#ifdef AFLR2_BOUNDARY_LAYER
    // Boundary Layer meshing related variables
    int propIndex, attrIndex = 0; // Indexing
    double capsMeshLength = 0;
    int createBL = (int) false; // Indicate if boundary layer meshes should be generated
#endif

    // AFRL2 arrays

    INT_1D *Bnd_Edge_Grid_BC_Flag = NULL;
    INT_1D *Bnd_Edge_ID_Flag = NULL;
    INT_1D *Bnd_Edge_Mesh_ID_Flag = NULL;
    INT_1D *Bnd_Edge_Error_Flag  =  NULL;
    INT_2D *Bnd_Edge_Connectivity = NULL;
    INT_3D *Tria_Connectivity = NULL;
    INT_4D *Quad_Connectivity = NULL;

    DOUBLE_2D *Coordinates = NULL;
    DOUBLE_1D *Initial_Normal_Spacing = NULL;
    DOUBLE_1D *BL_Thickness = NULL;

    INT_1D *BG_Bnd_Edge_Grid_BC_Flag = NULL;
    INT_1D *BG_Bnd_Edge_ID_Flag = NULL;
    INT_2D *BG_Bnd_Edge_Connectivity = NULL;
    INT_3D *BG_Tria_Neighbors = NULL;
    INT_3D *BG_Tria_Connectivity = NULL;

    DOUBLE_1D *BG_Spacing = NULL;
    DOUBLE_2D *BG_Coordinates = NULL;
    DOUBLE_3D *BG_Metric = NULL;

    DOUBLE_1D *Source_Spacing = NULL;
    DOUBLE_2D *Source_Coordinates = NULL;
    DOUBLE_3D *Source_Metric = NULL;

    INT_ Number_of_Bnd_Edges = 0;
    INT_ Number_of_Nodes = 0;
    INT_ Number_of_Quads = 0;
    INT_ Number_of_Trias = 0;

    INT_ Number_of_BG_Bnd_Edges = 0;
    INT_ Number_of_BG_Nodes = 0;
    INT_ Number_of_BG_Trias = 0;

    INT_ Number_of_Source_Nodes = 0;

    int ix = 0, iy = 1;

    int ntris, oclass, mtype, nloop, *ivec=NULL, *senses=NULL, *face_tris = NULL;
    ego geom, ref, *loops=NULL;
    double uvbox[4], result[18], *rvec=NULL, *face_xyz = NULL, *face_uv = NULL;
    char aimFile[PATH_MAX];

     // Commandline inputs
    int  prog_argc   = 1;    // Number of arguments
    char **prog_argv = NULL; // String arrays
    char *meshInputString = NULL;
    char *rest = NULL, *token = NULL;

    // Get bounding box for scaling
    status = EG_getBoundingBox(bodyIn, box);
    if (status != EGADS_SUCCESS) goto cleanup;

    size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
                (box[1]-box[4])*(box[1]-box[4]) +
                (box[2]-box[5])*(box[2]-box[5]));

    // Negating the first parameter triggers EGADS to only put vertexes on edges
    params[0] = -meshInput->paramTess[0]*size;
    params[1] =  meshInput->paramTess[1]*size;
    params[2] =  meshInput->paramTess[2];

    // Get Tessellation  - save to global variable
    status = EG_makeTessBody(bodyIn, params, &tess);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(tess, aimInfo, status);

    // Get bodyFaces
    status = EG_getBodyTopos(bodyIn, NULL, FACE, &numFace, &bodyFaces);
    AIM_STATUS(aimInfo, status);

    if (numFace != 1) {
      AIM_ERROR(aimInfo, "Body must have only one Face!!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }
    AIM_NOTNULL(bodyFaces, aimInfo, status);

    status = EG_getTopology(bodyFaces[0], &geom, &oclass, &mtype, uvbox,
                            &nloop, &loops, &senses);
    AIM_STATUS(aimInfo, status);

    if (geom->mtype != PLANE) {
      AIM_ERROR(aimInfo, "Body must be a PLANE surface!!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    status = EG_getGeometry(geom, &oclass, &mtype, &ref, &ivec, &rvec);
    AIM_STATUS(aimInfo, status);
    CROSS(result, rvec+3, rvec+6);

    if        (fabs(result[0])     < 1e-7 &&
               fabs(result[1])     < 1e-7 &&
               fabs(result[2] - 1) < 1e-7) { // z-constant plane
      ix = 0;
      iy = 1;
    } else if (fabs(result[0])     < 1e-7 &&
               fabs(result[1] - 1) < 1e-7 &&
               fabs(result[2])     < 1e-7) { // y-constant plane
      ix = 0;
      iy = 2;
    } else if (fabs(result[0] - 1) < 1e-7 &&
               fabs(result[1])     < 1e-7 &&
               fabs(result[2])     < 1e-7) { // x-constant plane
      ix = 2;
      iy = 1;
    } else {
      AIM_ERROR(aimInfo, "Body must be a PLANE surface aligned in a Cartesian plane!!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // extract the boundary tessellation from the egads body in uv-space
    status = egads_eval_bedge (aimInfo, tess, groupMap, meshMap,
                               &Number_of_Bnd_Edges, &Bnd_Edge_Connectivity,
                               &Bnd_Edge_ID_Flag, &Bnd_Edge_Mesh_ID_Flag, &Coordinates);
    AIM_STATUS(aimInfo, status);

    Number_of_Nodes = Number_of_Bnd_Edges;

//#define DUMP_TECPLOT_DEBUG_FILE
#ifdef DUMP_TECPLOT_DEBUG_FILE
    {
        FILE *fp = aim_fopen(aimInfo, "aflr2_debug.dat", "w");
        fprintf(fp, "VARIABLES = U, V, ID\n");

        fprintf(fp, "ZONE N=%d, E=%d, F=FEPOINT, ET=LINESEG\n", Number_of_Nodes, Number_of_Bnd_Edges);
        for (int i = 0; i < Number_of_Nodes; i++)
            fprintf(fp, "%22.15e %22.15e %d\n", Coordinates[i+1][0], Coordinates[i+1][1], Bnd_Edge_ID_Flag[i+1]);
        for (int i = 0; i < Number_of_Bnd_Edges; i++)
            fprintf(fp, "%d %d\n", Bnd_Edge_Connectivity[i+1][0], Bnd_Edge_Connectivity[i+1][1]);

        fclose(fp); fp = NULL;
    }
#endif


#ifdef AFLR2_BOUNDARY_LAYER
    // Loop through meshing properties and see if boundaryLayerThickness and boundaryLayerSpacing have been specified
    for (propIndex = 0; propIndex < numMeshProp && createBL == (int) false; propIndex++) {

      // If no boundary layer specified in meshProp continue
      if (meshProp[propIndex].boundaryLayerSpacing == 0) continue;

      // Loop through Elements and see if marker match
      for (i = 1; i <= Number_of_Bnd_Edges && createBL == (int) false; i++) {

        //If they don't match continue
        if (Bnd_Edge_ID_Flag[i] != meshProp[propIndex].attrIndex) continue;

        // Set face "bl" flag
        createBL = (int) true;
      }
    }

    // Get the capsMeshLenght if boundary layer meshing has been requested
    if (createBL == (int) true) {
        status = check_CAPSMeshLength(1, &bodyIn, &capsMeshLength);

        if (capsMeshLength <= 0 || status != CAPS_SUCCESS) {
            printf("**********************************************************\n");
            if (status != CAPS_SUCCESS)
              printf("capsMeshLength is not set on any body.\n");
            else
              printf("capsMeshLength: %f\n", capsMeshLength);
            printf("\n");
            printf("The capsMeshLength attribute must present on at least\n"
                   "one body for boundary layer generation.\n"
                   "\n"
                   "capsMeshLength should be a a positive value representative\n"
                   "of a characteristic length of the geometry,\n"
                   "e.g. the MAC of a wing or diameter of a fuselage.\n");
            printf("**********************************************************\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    // Loop through surface meshes to set boundary layer parameters

    // Allocate Initial_Normal_Spacing
    Initial_Normal_Spacing = (DOUBLE_1D *) ug_malloc (&status, (Number_of_Nodes+1) * sizeof (DOUBLE_1D));
    if (Initial_Normal_Spacing == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

    // Set default to none
    for (i = 1; i <= Number_of_Nodes; i++) {
        Initial_Normal_Spacing[i] = 0.0;
    }

    // Loop through meshing properties and see if boundaryLayerThickness and boundaryLayerSpacing have been specified
    for (propIndex = 0; propIndex < numMeshProp; propIndex++) {

        // If no boundary layer specified in meshProp continue
        if (meshProp[propIndex].boundaryLayerSpacing == 0) continue;

        // Loop through Elements and see if marker match
        for (i = 1; i <= Number_of_Bnd_Edges; i++) {

            //If they don't match continue
            if (Bnd_Edge_Mesh_ID_Flag[i] != meshProp[propIndex].attrIndex) continue;

            // Set boundary layer spacing
            Initial_Normal_Spacing[Bnd_Edge_Connectivity[i][0]] =
                        meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
            Initial_Normal_Spacing[Bnd_Edge_Connectivity[i][1]] =
                        meshProp[propIndex].boundaryLayerSpacing*capsMeshLength;
        }
    }
#endif //AFLR2_BOUNDARY_LAYER


    // Initialize and setup AFLR2 input parameter structure

    ug_set_prog_param_code (2);

    ug_set_prog_param_function1 (ug_initialize_aflr_param);
    ug_set_prog_param_function1 (ug_gq_initialize_param); // optional
    ug_set_prog_param_function2 (aflr2_initialize_param);
    ug_set_prog_param_function2 (ice2_initialize_param);

    status = ug_add_new_arg (&prog_argv, (char*)"allocate_and_initialize_argv");
    AIM_STATUS(aimInfo, status);

    // Parse input string
    if (meshInput->aflr4Input.meshInputString != NULL) {

        AIM_STRDUP(meshInputString, meshInput->aflr4Input.meshInputString, aimInfo, status);
        rest = meshInputString;
        while ((token = strtok_r(rest, " ", &rest))) {
            status = ug_add_flag_arg (token, &prog_argc, &prog_argv);
            if (status != CAPS_SUCCESS) {
                AIM_ERROR(aimInfo, "Failed to parse input string: %s", token);
                AIM_ADDLINE(aimInfo, "Complete input string: %s\n", meshInputString);
                goto cleanup;
            }
        }

        /* printf("Number of args = %d\n",prog_argc);
        for (i = 0; i <prog_argc ; i++) printf("Arg %d = %s\n", i, prog_argv[i]);
         */
    } else {
        // Set other command options
      //if (createBL == (int) true) {
      //    status = ug_add_flag_arg ((char*)"bl", &prog_argc, &prog_argv);
      //}
      if (status != 0) goto cleanup;
    }

    // set the last argument to 1 to print what arguments have been set
/*@-nullpass*/
    status = ug_check_prog_param (prog_argv, prog_argc, Message_Flag);
    AIM_STATUS(aimInfo, status); // Add error code

    status = aflr2_grid_generator (prog_argc, prog_argv,
                                   Message_Flag,
                                   &Number_of_Bnd_Edges, &Number_of_Nodes,
                                   &Number_of_Quads, &Number_of_Trias,
                                   &Number_of_BG_Bnd_Edges,
                                   &Number_of_BG_Nodes, &Number_of_BG_Trias,
                                   &Number_of_Source_Nodes,
                                   &Bnd_Edge_Connectivity,
                                   &Bnd_Edge_Error_Flag,
                                   &Bnd_Edge_Grid_BC_Flag,
                                   &Bnd_Edge_ID_Flag,
                                   &Quad_Connectivity, &Tria_Connectivity,
                                   &BG_Bnd_Edge_Connectivity,
                                   &BG_Bnd_Edge_Grid_BC_Flag,
                                   &BG_Bnd_Edge_ID_Flag,
                                   &BG_Tria_Neighbors,
                                   &BG_Tria_Connectivity,
                                   &Coordinates,
                                   Initial_Normal_Spacing,
                                   &BG_Coordinates,
                                   &BG_Spacing,
                                   &BG_Metric,
                                   &Source_Coordinates,
                                   &Source_Spacing,
                                   &Source_Metric);
    AIM_STATUS(aimInfo, status); // Add error code
/*@+nullpass*/
    AIM_NOTNULL(Coordinates          , aimInfo, status);
    AIM_NOTNULL(Tria_Connectivity    , aimInfo, status);
    AIM_NOTNULL(Quad_Connectivity    , aimInfo, status);
    AIM_NOTNULL(Bnd_Edge_ID_Flag     , aimInfo, status);
    AIM_NOTNULL(Bnd_Edge_Connectivity, aimInfo, status);

    ntris = Number_of_Trias + 2*Number_of_Quads;

    AIM_ALLOC(face_xyz , 3*Number_of_Nodes, double, aimInfo, status);
    AIM_ALLOC(face_uv  , 2*Number_of_Nodes, double, aimInfo, status);
    AIM_ALLOC(face_tris, 3*ntris          , int   , aimInfo, status);

    for (i = 0; i < Number_of_Nodes; i++) {
      status = EG_evaluate(bodyFaces[0], Coordinates[i+1], result);
      AIM_STATUS(aimInfo, status);

      face_uv[2*i+0] = Coordinates[i+1][0];
      face_uv[2*i+1] = Coordinates[i+1][1];

      for (j = 0; j < 3; j++)
        face_xyz[3*i+j] = result[j];

      // replace uv-coordinates with Cartesian coordinates
      Coordinates[i+1][0] = result[ix];
      Coordinates[i+1][1] = result[iy];
    }

    // get the xyz coordinates on the boundary loops
    status = egads_xyz_bedge(aimInfo, tess, ix, iy, face_xyz, Coordinates);
    AIM_STATUS(aimInfo, status);

    // Initiate triangle elements
    for (i = 0; i < Number_of_Trias; i++) {
      face_tris[3*i+0] = Tria_Connectivity[i+1][0];
      face_tris[3*i+1] = Tria_Connectivity[i+1][1];
      face_tris[3*i+2] = Tria_Connectivity[i+1][2];
    }

    // Initiate triangle elements
    for (i = 0; i < Number_of_Quads; i++) {
      face_tris[3*Number_of_Trias + 6*i+0] = Quad_Connectivity[i+1][0];
      face_tris[3*Number_of_Trias + 6*i+1] = Quad_Connectivity[i+1][1];
      face_tris[3*Number_of_Trias + 6*i+2] = Quad_Connectivity[i+1][2];

      face_tris[3*Number_of_Trias + 6*i+3] = Quad_Connectivity[i+1][0];
      face_tris[3*Number_of_Trias + 6*i+4] = Quad_Connectivity[i+1][2];
      face_tris[3*Number_of_Trias + 6*i+5] = Quad_Connectivity[i+1][3];
    }

    /* open the tessellation and set the face */
    status = EG_openTessBody( tess );
    AIM_STATUS(aimInfo, status);

    status = EG_setTessFace(tess, 1,
                            Number_of_Nodes,
                            face_xyz,
                            face_uv,
                            ntris,
                            face_tris);
    AIM_STATUS(aimInfo, status);

    if (Number_of_Quads > 0) {
      status = EG_attributeAdd(tess, ".mixed", ATTRINT, 1, &Number_of_Quads, NULL, NULL);
      AIM_STATUS(aimInfo, status);
    }

    status = EG_statusTessBody(tess, &body, &i, &j);
    AIM_STATUS(aimInfo, status, "Tessellation object was not built correctly!!!");

    // register the new tessellation with CAPS
    status = aim_newTess(aimInfo, tess);
    AIM_STATUS(aimInfo, status);

    snprintf(aimFile, PATH_MAX, "%s.lb8.ugrid", meshRef->fileName);

/*@-nullpass*/
    status = ug_io_write_2d_grid_file(aimFile,
                                      Message_Flag,
                                      Number_of_Bnd_Edges,
                                      Number_of_Nodes,
                                      Number_of_Quads,
                                      Number_of_Trias,
                                      Bnd_Edge_Connectivity,
                                      Bnd_Edge_Grid_BC_Flag,
                                      Bnd_Edge_ID_Flag,
                                      Quad_Connectivity,
                                      Tria_Connectivity,
                                      Coordinates,
                                      Initial_Normal_Spacing,
                                      BL_Thickness);
    AIM_STATUS(aimInfo, status);
/*@+nullpass*/

    surfaceMesh->meshType = Surface2DMesh;
    surfaceMesh->meshQuickRef.numTriangle      = Number_of_Trias;
    surfaceMesh->meshQuickRef.numQuadrilateral = Number_of_Quads;
    surfaceMesh->meshQuickRef.numLine          = Number_of_Bnd_Edges;

    surfaceMesh->meshQuickRef.startIndexTriangle      = 0;
    surfaceMesh->meshQuickRef.startIndexQuadrilateral = Number_of_Trias;
    surfaceMesh->meshQuickRef.startIndexLine          = Number_of_Trias + Number_of_Quads;

    surfaceMesh->meshQuickRef.useStartIndex = (int) true;

    // Nodes
    surfaceMesh->numNode = Number_of_Nodes;

    AIM_ALLOC(surfaceMesh->node, surfaceMesh->numNode, meshNodeStruct, aimInfo, status);

    // Initiate nodes
    for (i = 0; i < surfaceMesh->numNode; i++) {
        status = initiate_meshNodeStruct(&surfaceMesh->node[i], UnknownMeshAnalysis);
        AIM_STATUS(aimInfo, status);

        surfaceMesh->node[i].xyz[ix] = Coordinates[i+1][0];
        surfaceMesh->node[i].xyz[iy] = Coordinates[i+1][1];
    }

    // Allocate and initialize all elements
    surfaceMesh->numElement = Number_of_Trias + Number_of_Quads + Number_of_Bnd_Edges;
    AIM_ALLOC(surfaceMesh->element, surfaceMesh->numElement, meshElementStruct, aimInfo, status);

    // Initiate all elements
    for (i = 0; i < surfaceMesh->numElement; i++) {
        status = initiate_meshElementStruct(surfaceMesh->element + i, UnknownMeshAnalysis);
        AIM_STATUS(aimInfo, status);
    }

    // Initiate triangle elements
    for (i = 0; i < Number_of_Trias; i++) {
        elementIndex = i;
        surfaceMesh->element[elementIndex].elementType = Triangle;

        status = mesh_allocMeshElementConnectivity(&surfaceMesh->element[elementIndex]);
        AIM_STATUS(aimInfo, status);

        surfaceMesh->element[elementIndex].connectivity[0] = Tria_Connectivity[i+1][0];
        surfaceMesh->element[elementIndex].connectivity[1] = Tria_Connectivity[i+1][1];
        surfaceMesh->element[elementIndex].connectivity[2] = Tria_Connectivity[i+1][2];

        surfaceMesh->element[elementIndex].elementID = elementIndex+1;
        surfaceMesh->element[elementIndex].markerID = faceAttr;
    }

    // Initiate quad elements
    for (i = 0; i < Number_of_Quads; i++) {

        elementIndex = Number_of_Trias + i;

        surfaceMesh->element[elementIndex].elementType = Quadrilateral;

        status = mesh_allocMeshElementConnectivity(&surfaceMesh->element[elementIndex]);
        AIM_STATUS(aimInfo, status);

        surfaceMesh->element[elementIndex].connectivity[0] = Quad_Connectivity[i+1][0];
        surfaceMesh->element[elementIndex].connectivity[1] = Quad_Connectivity[i+1][1];
        surfaceMesh->element[elementIndex].connectivity[2] = Quad_Connectivity[i+1][2];
        surfaceMesh->element[elementIndex].connectivity[3] = Quad_Connectivity[i+1][3];

        surfaceMesh->element[elementIndex].elementID = elementIndex+1;
        surfaceMesh->element[elementIndex].markerID = faceAttr;

    }

    for (i = 0; i < Number_of_Bnd_Edges; i++) {

        elementIndex = Number_of_Trias + Number_of_Quads + i;

        surfaceMesh->element[elementIndex].elementType = Line;

        status = mesh_allocMeshElementConnectivity(&surfaceMesh->element[elementIndex]);
        AIM_STATUS(aimInfo, status);

        surfaceMesh->element[elementIndex].connectivity[0] = Bnd_Edge_Connectivity[i+1][0];
        surfaceMesh->element[elementIndex].connectivity[1] = Bnd_Edge_Connectivity[i+1][1];

        surfaceMesh->element[elementIndex].elementID = elementIndex+1;
        surfaceMesh->element[elementIndex].markerID = Bnd_Edge_ID_Flag[i+1];
    }

    status = CAPS_SUCCESS;

cleanup:
    if (status != CAPS_SUCCESS)
      (void) EG_deleteObject(tess);

    AIM_FREE(meshInputString);
    AIM_FREE(bodyFaces);
    AIM_FREE(ivec);
    AIM_FREE(rvec);
    AIM_FREE(face_uv);
    AIM_FREE(face_xyz);
    AIM_FREE(face_tris);

    //Free the arguments
    if (prog_argv != NULL)
      ug_free_argv (prog_argv);

    ug_free(Bnd_Edge_Grid_BC_Flag);
    ug_free(Bnd_Edge_ID_Flag);
    ug_free(Bnd_Edge_Mesh_ID_Flag);
    ug_free(Bnd_Edge_Error_Flag);
    ug_free(Bnd_Edge_Connectivity);
    ug_free(Tria_Connectivity);
    ug_free(Quad_Connectivity);

    ug_free(Coordinates);
    ug_free(Initial_Normal_Spacing);

    ug_free(BG_Bnd_Edge_Grid_BC_Flag);
    ug_free(BG_Bnd_Edge_ID_Flag);
    ug_free(BG_Bnd_Edge_Connectivity);
    ug_free(BG_Tria_Neighbors);
    ug_free(BG_Tria_Connectivity);

    ug_free(BG_Spacing);
    ug_free(BG_Coordinates);
    ug_free(BG_Metric);

    ug_free(Source_Spacing);
    ug_free(Source_Coordinates);
    ug_free(Source_Metric);

    return status;
}
