/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Tessellation Element Order Elevation AIM
 *
 *      Written by Marshall Galbraith, MIT.
 *        Based on that of Ryan Durscher and Dean E. Bryson, AFRL/RQVC
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

/*!\mainpage Introduction
 *
 * \section overviewCurveTess Curve Tessellation AIM Overview
 * This AIM provides the ability to elevate a triangular linear surfeace mesh to a high-order, "curved", surface meshes.
 * The algorithm only inserts high-order vertexes on element interior and edges and ensures they are on the geometry.
 * However, the original linear mesh vertexes are not modified. Hence, this algirthm is only suitable for isotropic meshes as the vertex insertion may produce negative Jacobian's for anisotpric elements.
 *
 * An outline of the AIM's inputs and outputs are provided in \ref aimInputsCurveTess and \ref aimOutputsCurveTess, respectively.
 *
 * \section clearanceCurveTess Clearance Statement
 *  This software has been cleared for public release.
 *
 */

#include <string.h>
#include <math.h>

#include "capsTypes.h"
#include "aimUtil.h"

#include "meshUtils.h"
#include "miscUtils.h"

#include "egads.h"

#ifdef WIN32
#define strcasecmp stricmp
#endif

//#define DEBUG

/* external functions not in headers */
extern int
EG_tessHOverts(const ego tess,
               int nstx, const double *st,
               int nItrix, const int *iTris,
               ego *nTess);

#define CURVETESSFILE "curveTess_%d_%d.eto"

enum aimInputs
{
    inProj_Name = 1,                /* index is 1-based */
    inElement_Order,
    inElement_Class,
    inMesh_Quiet_Flag,
    inMesh_Format,
    inMesh_Morph,
    inSurface_Mesh,
    NUMINPUT = inSurface_Mesh       /* Total number of inputs */
};

enum aimOutputs
{
  outNumberOfNode = 1,           /* index is 1-based */
  outNumberOfElement,
  outSurface_Mesh,
  NUMOUT = outSurface_Mesh       /* Total number of outputs */
};


typedef struct {

  // Mesh reference obtained from meshing AIM
  int numMeshRefIn;
  aimMeshRef *meshRefIn;

  // container for elevated surface mesh(es)
  int numMeshRefOut;
  aimMeshRef *meshRefOut;

  int numNodeTotal;
  int numElemTotal;

  // Attribute to index map
  mapAttrToIndexStruct groupMap;

} aimStorage;


static int initiate_aimStorage(aimStorage *curveTessInstance)
{
  // Mesh reference obtained from meshing AIM
  curveTessInstance->numMeshRefIn = 0;
  curveTessInstance->meshRefIn = NULL;

  // container for elevated surface mesh(es)
  curveTessInstance->numMeshRefOut = 0;
  curveTessInstance->meshRefOut = NULL;

  curveTessInstance->numNodeTotal = 0;
  curveTessInstance->numElemTotal = 0;

  initiate_mapAttrToIndexStruct(&curveTessInstance->groupMap);

  return CAPS_SUCCESS;

}


static int destroy_aimStorage(aimStorage *curveTessInstance)
{
  int i;

  // input meshes
  curveTessInstance->numMeshRefIn = 0;
  curveTessInstance->meshRefIn = NULL;

  // Destroy elevated mesh
  for (i = 0; i < curveTessInstance->numMeshRefOut; i++) {
    aim_freeMeshRef(&curveTessInstance->meshRefOut[i]);
  }
  AIM_FREE(curveTessInstance->meshRefOut);
  curveTessInstance->numMeshRefOut = 0;

  curveTessInstance->numNodeTotal = 0;
  curveTessInstance->numElemTotal = 0;

  destroy_mapAttrToIndexStruct(&curveTessInstance->groupMap);

  return CAPS_SUCCESS;
}


/* ********************** Exposed AIM Functions ***************************** */
int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, /*@unused@*/char ***fnames, /*@unused@*/int **franks, /*@unused@*/int **fInOut)
{
  int  status = CAPS_SUCCESS;

  aimStorage *curveTessInstance=NULL;

#ifdef DEBUG
  printf("curveTessAIM/aimInitialize   instance = %d!\n", inst);
#endif

  /* specify the number of analysis input and out "parameters" */
  *nIn     = NUMINPUT;
  *nOut    = NUMOUT;
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 0;

  // Allocate nastranInstance
  AIM_ALLOC(curveTessInstance, 1, aimStorage, aimInfo, status);
  *instStore = curveTessInstance;

  // Initialize instance storage
  (void) initiate_aimStorage(curveTessInstance);

cleanup:

  if (status != CAPS_SUCCESS) {
    AIM_FREE(*instStore);
  }

  return status;
}


// ********************** AIM Function Break *****************************
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{
  /*! \page aimInputsCurveTess AIM Inputs
   * The following list outlines the curveTess inputs along with their default value available
   * through the AIM interface.
   */

  int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" curveTessAIM/aimInputs index  = %d!\n", index);
#endif

  *ainame = NULL;


  if (index == inProj_Name) {
    *ainame              = EG_strdup("Proj_Name");
    defval->type         = String;
    defval->nullVal      = NotAllowed;
    AIM_STRDUP(defval->vals.string, "curveTess_CAPS", aimInfo, status);

    /*! \page aimInputsCurveTess
     * - <B> Proj_Name = "curveTess_CAPS"</B> <br>
     * Output name prefix for meshes to be written in formats specified by Mesh_Format.
     * These meshes are not linked to any analysis, but may be useful exploring meshing parameters.
     */

  } else if (index == inElement_Order) {
    *ainame              = EG_strdup("Element_Order");
    defval->type         = Integer;
    defval->vals.integer = 2;

    /*! \page aimInputsCurveTess
     * - <B> Element_Order = 2</B> <br>
     * Polynomial order for the elevated elements.<br>
     *  1 - liner<br>
     *  2 - quadratic<br>
     *  3 - qubic<br>
     *  etc.<br>
     */

  } else if (index == inElement_Class) {
      *ainame              = EG_strdup("Element_Class");
      defval->type         = Integer;
      defval->vals.integer = 1;

      /*! \page aimInputsCurveTess
       * - <B> Element_Class = 1</B> <br>
       * Class for the elevation for different Element_Order.<br>
       * (images from Exodus manual: https://sandialabs.github.io/seacas-docs/html/element_types.html)
       * <br>
       *  <table>
       *  <caption id="linear">Element_Order = 1 - liner</caption>
       *  <tr><th>Element_Class = 1  <th>Element_Class = 2
       *  <tr><td>
       *      \image html  tri3.png
       *      \image latex tri3.png  width=0.3\textwidth
       *      <td>
       *      \image html  tri4.png
       *      \image latex tri4.png  width=0.3\textwidth
       *  \if false
       *  <tr><td>
       *      \image html  quad4.png
       *      \image latex quad4.png width=0.3\textwidth
       *      <td>
       *      \image html  quad5.png
       *      \image latex quad5.png width=0.3\textwidth
       *  \endif
       *  </table>
       *<br>
       *  <table>
       *  <caption id="quadratic">Element_Order = 2 - quadratic</caption>
       *  <tr><th>Element_Class = 1  <th>Element_Class = 2
       *  <tr><td>
       *      \image html  tri6.png
       *      \image latex tri6.png  width=0.3\textwidth
       *      <td>
       *      \image html  tri7.png
       *      \image latex tri7.png  width=0.3\textwidth
       *  \if false
       *      <td>
       *      \image html  tri6.png
       *      \image latex tri6.png  width=0.3\textwidth
       *  <tr><td>
       *      \image html  quad9.png
       *      \image latex quad9.png width=0.3\textwidth
       *      <td>
       *      \image html  quad9.png
       *      \image latex quad9.png width=0.3\textwidth
       *      <td>
       *      \image html  quad8.png
       *      \image latex quad8.png width=0.3\textwidth
       *  \endif
       *  </table>
       *
       *  3 - qubic<br>
       *  etc.<br>
       */

  } else if (index == inMesh_Quiet_Flag) {
    *ainame               = EG_strdup("Mesh_Quiet_Flag");
    defval->type          = Boolean;
    defval->vals.integer  = false;

    /*! \page aimInputsCurveTess
     * - <B> inMesh_Quiet_Flag = False</B> <br>
     * Complete suppression of mesh generator (not including errors)
     */

  } else if (index == inMesh_Format) {
    *ainame               = EG_strdup("Mesh_Format");
    defval->type          = String;
    defval->vals.string   = NULL;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;
    defval->lfixed        = Change;

    /*! \page aimInputsCurveTess
     * \include{doc} Mesh_Format.dox
     */

  } else if (index == inMesh_Morph) {
    *ainame              = EG_strdup("Mesh_Morph");
    defval->type         = Boolean;
    defval->lfixed       = Fixed;
    defval->vals.integer = (int) false;
    defval->dim          = Scalar;
    defval->nullVal      = NotNull;

    /*! \page aimInputsCurveTess
     * - <B> Mesh_Morph = False</B> <br>
     * Project previous surface mesh onto new geometry.
     */

  } else if (index == inSurface_Mesh) {
    *ainame             = EG_strdup("Surface_Mesh");
    defval->type        = PointerMesh;
    defval->dim         = Vector;
    defval->lfixed      = Change;
    defval->sfixed      = Change;
    defval->vals.AIMptr = NULL;
    defval->nullVal     = IsNull;

    /*! \page aimInputsCurveTess
     * - <B>Mesh = NULL</B> <br>
     * A Surface_Mesh link.
     */
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

  cleanup:
  if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
  return status;
}


// ********************** AIM Function Break *****************************
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
               int index, char **aoname, capsValue *form)
{
  /*! \page aimOutputsCurveTess AIM Outputs
   * The following list outlines the curveTess outputs available through the AIM interface.
   */

  int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" curveTessAIM/aimOutputs index = %d!\n", index);
#endif


  if (index == outNumberOfNode) {
    *aoname = EG_strdup("NumberOfNode");
    form->type = Integer;
    form->vals.integer = 0;

    /*! \page aimOutputsCurveTess
     * - <B> NumberOfNode </B> <br>
     * Number of vertices in the surface mesh
     */

  } else if (index == outNumberOfElement) {
    *aoname = EG_strdup("NumberOfElement");
    form->type = Integer;
    form->vals.integer = 0;

    /*! \page aimOutputsCurveTess
     * - <B> NumberOfElement </B> <br>
     * Number of elements in the surface mesh
     */

  } else if (index == outSurface_Mesh) {

    *aoname           = AIM_NAME(Surface_Mesh);
    form->type        = PointerMesh;
    form->dim         = Vector;
    form->lfixed      = Change;
    form->sfixed      = Change;
    form->vals.AIMptr = NULL;
    form->nullVal     = IsNull;

    /*! \page aimOutputsCurveTess
     * - <B> Surface_Mesh </B> <br>
     * The elevated surface mesh for a link
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
int aimUpdateState(void *instStore, void *aimInfo,
                   capsValue *aimInputs)
{
  // Function return flag
  int status = CAPS_SUCCESS;

  int i, meshIndex;
  int imap;

  int atype, alen;

  // AIM input bodies
  const char *intents;
  int  numBody;
  ego *bodies = NULL;

  int state, nglobal;
  ego body;

  char aimFile[PATH_MAX];
  char filename[PATH_MAX];

  aimStorage *curveTessInstance = (aimStorage *)instStore;
  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(bodies, aimInfo, status);

#ifdef DEBUG
  printf(" curveTessAIM/aimUpdateState numBody = %d!\n", numBody);
#endif

  if (aimInputs[inSurface_Mesh-1].nullVal == IsNull &&
      aimInputs[inMesh_Morph-1].vals.integer == (int) false) {
    AIM_ANALYSISIN_ERROR(aimInfo, inSurface_Mesh, "'Surface_Mesh' input must be linked to an output 'Surface_Mesh'");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Cleanup previous aimStorage for the instance in case this is the second time through preAnalysis for the same instance
  status = destroy_aimStorage(curveTessInstance);
  AIM_STATUS(aimInfo, status);

  status = create_CAPSGroupAttrToIndexMap(numBody,
                                          bodies,
                                          3, // Body, Faces, and Edges
                                          &curveTessInstance->groupMap);
  AIM_STATUS(aimInfo, status);


  // get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  if (status != CAPS_SUCCESS) return status;

#ifdef DEBUG
  printf(" curveTessAIM/aimPreAnalysis numBody = %d!\n", numBody);
#endif

  if (numBody <= 0 || bodies == NULL) {
#ifdef DEBUG
    printf(" curveTessAIM/aimPreAnalysis No Bodies!\n");
#endif
    return CAPS_SOURCEERR;
  }

  // Get mesh
  if (aimInputs[inSurface_Mesh-1].nullVal == IsNull) {
    AIM_ANALYSISIN_ERROR(aimInfo, inSurface_Mesh, "'Surface_Mesh' input must be linked to an output 'Surface_Mesh'");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  curveTessInstance->numMeshRefIn = aimInputs[inSurface_Mesh-1].length;
  curveTessInstance->meshRefIn    = (aimMeshRef *)aimInputs[inSurface_Mesh-1].vals.AIMptr;

  for ( meshIndex = 0; meshIndex < curveTessInstance->numMeshRefIn; meshIndex++ ) {
    for ( imap = 0; imap < curveTessInstance->meshRefIn[meshIndex].nmap; imap++ ) {

      // check if the tessellation has a mixture of quad and tri
      status = EG_attributeRet(curveTessInstance->meshRefIn[meshIndex].maps[imap].tess, ".mixed",
                               &atype, &alen, NULL, NULL, NULL);
      if (status == EGADS_SUCCESS) {
        AIM_ERROR(aimInfo, "Only triangle elements are currently supported");
        status = CAPS_SOURCEERR;
        goto cleanup;
      }

      status = EG_statusTessBody(curveTessInstance->meshRefIn[meshIndex].maps[imap].tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      for (i = 0; i < numBody; i++) {
        if (bodies[i] == body) break;
      }

      if (i == numBody) {
        AIM_ERROR(aimInfo, "Failed to match all bodies between input Surface_Mesh and curveTess bodies!");
        status = CAPS_SOURCEERR;
        goto cleanup;
      }

    }
  }

  AIM_ALLOC(curveTessInstance->meshRefOut, curveTessInstance->numMeshRefIn, aimMeshRef, aimInfo, status);
  for ( meshIndex = 0; meshIndex < curveTessInstance->numMeshRefIn; meshIndex++ ) {
    aim_initMeshRef(&curveTessInstance->meshRefOut[meshIndex], aimSurfaceMesh);
  }
  curveTessInstance->numMeshRefOut = curveTessInstance->numMeshRefIn;


  for ( meshIndex = 0; meshIndex < curveTessInstance->numMeshRefIn; meshIndex++ ) {

    // copy ID/groupName information
    AIM_ALLOC(curveTessInstance->meshRefOut[meshIndex].bnds, curveTessInstance->meshRefIn[meshIndex].nbnd, aimMeshBnd, aimInfo, status);
    curveTessInstance->meshRefOut[meshIndex].nbnd = curveTessInstance->meshRefIn[meshIndex].nbnd;
    for (i = 0; i < curveTessInstance->meshRefOut[meshIndex].nbnd; i++) {
      curveTessInstance->meshRefOut[meshIndex].bnds[i].groupName = NULL;
      curveTessInstance->meshRefOut[meshIndex].bnds[i].ID = 0;
    }

    for (i = 0; i < curveTessInstance->meshRefOut[meshIndex].nbnd; i++) {
      AIM_STRDUP(curveTessInstance->meshRefOut[meshIndex].bnds[i].groupName, curveTessInstance->meshRefIn[meshIndex].bnds[i].groupName, aimInfo, status);
      curveTessInstance->meshRefOut[meshIndex].bnds[i].ID = curveTessInstance->meshRefIn[meshIndex].bnds[i].ID;
    }

    // allocate map information
    AIM_ALLOC(curveTessInstance->meshRefOut[meshIndex].maps, curveTessInstance->meshRefIn[meshIndex].nmap, aimMeshTessMap, aimInfo, status);
    for ( imap = 0; imap < curveTessInstance->meshRefIn[meshIndex].nmap; imap++ ) {
      curveTessInstance->meshRefOut[meshIndex].maps[imap].map = NULL;
      curveTessInstance->meshRefOut[meshIndex].maps[imap].tess = NULL;
    }
    curveTessInstance->meshRefOut[meshIndex].nmap = curveTessInstance->meshRefIn[meshIndex].nmap;

    // set the filename without extensions where the grid is written for solvers
    snprintf(filename, PATH_MAX, "%s_%d", aimInputs[inProj_Name-1].vals.string, meshIndex);
    status = aim_file(aimInfo, filename, aimFile);
    AIM_STATUS(aimInfo, status);

    AIM_STRDUP(curveTessInstance->meshRefOut[meshIndex].fileName, aimFile, aimInfo, status);
  }


  status = CAPS_SUCCESS;

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
int aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{

  // general
  int meshIndex, status;
  int imap;

  char curveTessFile[42];
  char aimFile[PATH_MAX];

  // mesh elevation
  int nElevatedTris, nElevatedVerts;
  int *elevatedTris=NULL;
  double *elevatedVerts=NULL;

  const aimStorage *curveTessInstance;

  // AIM input bodies
  const char *intents;
  int  numBody;
  ego *bodies = NULL;

  int bodyIndex;
  ego tess = NULL;

  curveTessInstance = (const aimStorage *) instStore;
  if (aimInputs == NULL) return CAPS_NULLVALUE;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

  // element order
  if( aimInputs[inElement_Order-1].vals.integer == 1 ) {

    // element class
    if( aimInputs[inElement_Class-1].vals.integer == 1 ) {

      /*   Tri3
       *                     1: (0.0, 0.0)
       *                     2: (1.0, 0.0)
       *                     3: (0.0, 1.0)
       *
       *                3
       *               / \
       *              /   \
       *             /     \
       *            /       \
       *           /         \                   --i-- : local vertex index (total element scope)
       *          /           \                    |\
       *         /             \
       *        /               \                 ---
       *       /                 \               \ j   : internal vertex index (internal tri scope)
       *      /                   \
       *     /                     \              (k)  : internal tri index (total element scope)
       *    1-----------------------2
       */

      nElevatedTris = 1;

      AIM_ALLOC(elevatedTris, 3*nElevatedTris, int, aimInfo, status);

      elevatedTris[0]  = 1; // Tri 1
      elevatedTris[1]  = 2;
      elevatedTris[2]  = 3;

      nElevatedVerts = 3;

      AIM_ALLOC(elevatedVerts, 2*nElevatedVerts, double, aimInfo, status);

      elevatedVerts[0]  = 0.0; // u1, v1
      elevatedVerts[1]  = 0.0;
      elevatedVerts[2]  = 1.0; // u2, v2
      elevatedVerts[3]  = 0.0;
      elevatedVerts[4]  = 0.0; // u3, v3
      elevatedVerts[5]  = 1.0;

    } else if( aimInputs[inElement_Class-1].vals.integer == 2 ) {

      /*   Tri4
       *                     1: (0.0, 0.0)
       *                     2: (1.0, 0.0)
       *                     3: (0.0, 1.0)
       *                     4: (1/3, 1/3)
       *
       *                3
       *               /|\
       *              /1|2\
       *             /  |  \
       *            /   |   \
       *           /    |    \                   --i-- : local vertex index (total element scope)
       *          /     |     \                    |\
       *         / (3)  4  (2) \
       *        /     3/3\3     \                 ---
       *       /    _/     \_    \               \ j   : internal vertex index (internal tri scope)
       *      / __/    (1)    \__ \
       *     /2/1               2\1\              (k)  : internal tri index (total element scope)
       *    1-----------------------2
       */

      nElevatedTris = 3;

      AIM_ALLOC(elevatedTris, 3*nElevatedTris, int, aimInfo, status);

      elevatedTris[0]  = 1; // Tri 1
      elevatedTris[1]  = 2;
      elevatedTris[2]  = 4;
      elevatedTris[3]  = 2; // Tri 2
      elevatedTris[4]  = 3;
      elevatedTris[5]  = 4;
      elevatedTris[6]  = 3; // Tri 3
      elevatedTris[7]  = 1;
      elevatedTris[8]  = 4;

      nElevatedVerts = 4;

      AIM_ALLOC(elevatedVerts, 2*nElevatedVerts, double, aimInfo, status);

      elevatedVerts[0]  = 0.0; // u1, v1
      elevatedVerts[1]  = 0.0;
      elevatedVerts[2]  = 1.0; // u2, v2
      elevatedVerts[3]  = 0.0;
      elevatedVerts[4]  = 0.0; // u3, v3
      elevatedVerts[5]  = 1.0;
      elevatedVerts[6]  = 1./3.; // u4, v4
      elevatedVerts[7]  = 1./3.;

    } else {
      AIM_ERROR(aimInfo, "Unsupported Element_Order = %d Element_Class = %d",
                aimInputs[inElement_Order-1].vals.integer,
                aimInputs[inElement_Class-1].vals.integer);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

  } else if( aimInputs[inElement_Order-1].vals.integer == 2 ) {

    // element class
    if( aimInputs[inElement_Class-1].vals.integer == 1 ) {

      /*  Tri6
       *
       *    3                1: (0.0, 0.0)
       *    |\               2: (1.0, 0.0)
       *    |3\              3: (0.0, 1.0)
       *    |  \             4: (0.5, 0.0)
       *    |(3)\            5: (0.5, 0.5)
       *    |    \           6: (0.0, 0.5)
       *    |1   2\
       *    6------5         --i-- : local vertex index (total element scope)
       *    |\3   2|\          |\
       *    |3\    |3\
       *    |  \(4)|  \       ---
       *    |(1)\  |(2)\     \ j   : internal vertex index (internal tri scope)
       *    |    \1|    \
       *    |1   2\|1   2\    (k)  : internal tri index (total element scope)
       *    1------4------2
       */

      nElevatedTris = 4;

      AIM_ALLOC(elevatedTris, 3*nElevatedTris, int, aimInfo, status);

      elevatedTris[0]  = 1; // Tri 1
      elevatedTris[1]  = 4;
      elevatedTris[2]  = 6;
      elevatedTris[3]  = 4; // Tri 2
      elevatedTris[4]  = 2;
      elevatedTris[5]  = 5;
      elevatedTris[6]  = 6; // Tri 3
      elevatedTris[7]  = 5;
      elevatedTris[8]  = 3;
      elevatedTris[9]  = 4; // Tri 4
      elevatedTris[10] = 5;
      elevatedTris[11] = 6;

      nElevatedVerts = 6;

      AIM_ALLOC(elevatedVerts, 2*nElevatedVerts, double, aimInfo, status);

      elevatedVerts[0]  = 0.0; // u1, v1
      elevatedVerts[1]  = 0.0;
      elevatedVerts[2]  = 1.0; // u2, v2
      elevatedVerts[3]  = 0.0;
      elevatedVerts[4]  = 0.0; // u3, v3
      elevatedVerts[5]  = 1.0;
      elevatedVerts[6]  = 0.5; // u4, v4
      elevatedVerts[7]  = 0.0;
      elevatedVerts[8]  = 0.5; // u5, v5
      elevatedVerts[9]  = 0.5;
      elevatedVerts[10] = 0.0; // u6, v6
      elevatedVerts[11] = 0.5;

    } else if( aimInputs[inElement_Class-1].vals.integer == 2 ) {

      /*  Tri7
       *    3                1: (0.0, 0.0)
       *    |\               2: (1.0, 0.0)
       *    |3\              3: (0.0, 1.0)
       *    |  \             4: (0.5, 0.0)
       *    |(3)\            5: (0.5, 0.5)
       *    |    \           6: (0.0, 0.5)
       *    |1   2\          7: (1/3, 1/3)
       *    6------5
       *    |\     |\        --i-- : local vertex index (total element scope)
       *    |3\  7 |3\         |\
       *    |  \   |  \
       *    |(1)\  |(2)\      ---
       *    |    \ |    \    \ j   : internal vertex index (internal tri scope)
       *    |1   2\|1   2\
       *    1------4------2   (k)  : internal tri index (total element scope)
       */

      nElevatedTris = 6;

      AIM_ALLOC(elevatedTris, 3*nElevatedTris, int, aimInfo, status);

      elevatedTris[0]  = 1; // Tri 1
      elevatedTris[1]  = 4;
      elevatedTris[2]  = 6;
      elevatedTris[3]  = 4; // Tri 2
      elevatedTris[4]  = 2;
      elevatedTris[5]  = 5;
      elevatedTris[6]  = 6; // Tri 3
      elevatedTris[7]  = 5;
      elevatedTris[8]  = 3;
      elevatedTris[9]  = 4; // Tri 4
      elevatedTris[10] = 5;
      elevatedTris[11] = 7;
      elevatedTris[12] = 5; // Tri 5
      elevatedTris[13] = 6;
      elevatedTris[14] = 7;
      elevatedTris[15] = 6; // Tri 6
      elevatedTris[16] = 4;
      elevatedTris[17] = 7;

      nElevatedVerts = 7;

      AIM_ALLOC(elevatedVerts, 2*nElevatedVerts, double, aimInfo, status);

      elevatedVerts[0]  = 0.0; // u1, v1
      elevatedVerts[1]  = 0.0;
      elevatedVerts[2]  = 1.0; // u2, v2
      elevatedVerts[3]  = 0.0;
      elevatedVerts[4]  = 0.0; // u3, v3
      elevatedVerts[5]  = 1.0;
      elevatedVerts[6]  = 0.5; // u4, v4
      elevatedVerts[7]  = 0.0;
      elevatedVerts[8]  = 0.5; // u5, v5
      elevatedVerts[9]  = 0.5;
      elevatedVerts[10] = 0.0; // u6, v6
      elevatedVerts[11] = 0.5;
      elevatedVerts[12] = 1./3.; // u7, v7
      elevatedVerts[13] = 1./3.;

    } else {
      AIM_ERROR(aimInfo, "Unsupported Element_Order = %d Element_Class = %d",
                aimInputs[inElement_Order-1].vals.integer,
                aimInputs[inElement_Class-1].vals.integer);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

  } else {
    AIM_ERROR(aimInfo, "Unsupported Element_Order = %d", aimInputs[inElement_Order-1].vals.integer );
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // set the filename without extensions where the grid is written for solvers
  status = aim_file(aimInfo, aimInputs[inProj_Name-1].vals.string, aimFile);
  AIM_STATUS(aimInfo, status);

  /*
   * finally, do the tessellation order elevation and put the result back in a mesh struct
   */

  bodyIndex = 0;
  for ( meshIndex = 0; meshIndex < curveTessInstance->numMeshRefIn; meshIndex++ ) {
    for ( imap = 0; imap < curveTessInstance->meshRefIn[meshIndex].nmap; imap++ ) {

      if (aimInputs[inMesh_Quiet_Flag-1].vals.integer == (int)false)
        printf("Curving surface mesh for body %d (of %d)\n", bodyIndex+1, numBody);

      // make the call to EG_tessHOverts
      status = EG_tessHOverts(curveTessInstance->meshRefIn[meshIndex].maps[imap].tess,
                              nElevatedVerts,
                              elevatedVerts,
                              nElevatedTris,
                              elevatedTris,
                              &tess);
      AIM_STATUS(aimInfo, status);
      AIM_NOTNULL(tess, aimInfo, status);

      // set the file name to write the egads file
      snprintf(curveTessFile, 42, CURVETESSFILE, meshIndex, imap);
      status = aim_rmFile(aimInfo, curveTessFile);
      AIM_STATUS(aimInfo, status);

      status = aim_file(aimInfo, curveTessFile, aimFile);
      AIM_STATUS(aimInfo, status);

      status = EG_saveTess(tess, aimFile);
      AIM_STATUS(aimInfo, status);

      EG_deleteObject(tess);
      bodyIndex++;
    }
  }

  status = CAPS_SUCCESS;

cleanup:

  AIM_FREE(elevatedTris);
  AIM_FREE(elevatedVerts);

  return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ const void *instStore, /*@unused@*/ void *aimInfo,
               int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int aimPostAnalysis(void *instStore, void *aimInfo,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *aimInputs)
{
  int i, meshIndex, status;
  int imap, bodyIndex;

  // AIM input bodies
  const char *intents;
  int  numBody;
  ego *bodies = NULL;

  char curveTessFile[42];
  char aimFile[PATH_MAX];
  aimMesh    mesh;

  int state, nglobal;
  ego body;

  int numFace, numQuad, numTri;
  int iface;
  int plen, tlen;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;

  int atype, nHOtris;

  aimStorage *curveTessInstance;

  curveTessInstance = (aimStorage *) instStore;
  if (aimInputs == NULL) return CAPS_NULLVALUE;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
  AIM_STATUS(aimInfo, status);

  bodyIndex = 0;
  for ( meshIndex = 0; meshIndex < curveTessInstance->numMeshRefIn; meshIndex++ ) {
    for ( imap = 0; imap < curveTessInstance->meshRefIn[meshIndex].nmap; imap++ ) {

      // set the file name to write the egads file
      snprintf(curveTessFile, 42, CURVETESSFILE, meshIndex, imap);

      status = aim_file(aimInfo, curveTessFile, aimFile);
      AIM_STATUS(aimInfo, status);

      status = EG_statusTessBody(curveTessInstance->meshRefIn[meshIndex].maps[imap].tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      status = EG_loadTess(body, aimFile, &curveTessInstance->meshRefOut[meshIndex].maps[imap].tess);
      AIM_STATUS(aimInfo, status);

      status = aim_newTess(aimInfo, curveTessInstance->meshRefOut[meshIndex].maps[imap].tess);
      AIM_STATUS(aimInfo, status);

      status = EG_attributeRet(curveTessInstance->meshRefOut[meshIndex].maps[imap].tess, ".HOtris",
                               &atype, &nHOtris, NULL, NULL, NULL);
      AIM_STATUS(aimInfo, status);

      // Get faces, edges, and nodes so we can check for attributes on them
      status = EG_getBodyTopos(body, NULL, FACE, &numFace, NULL);
      AIM_STATUS(aimInfo, status);

      numQuad = numTri = 0;
      for (iface = 1; iface <= numFace; iface++) {
        status = EG_getTessFace(curveTessInstance->meshRefOut[meshIndex].maps[imap].tess,
                                iface, &plen, &points, &uv, &ptype, &pindex,
                                &tlen, &tris, &tric);
        AIM_STATUS(aimInfo, status);

        numTri += tlen/nHOtris;
      }

      status = EG_statusTessBody(curveTessInstance->meshRefOut[meshIndex].maps[imap].tess, &body, &state, &nglobal);
      AIM_STATUS(aimInfo, status);

      if (restart == 0 &&
          aimInputs[inMesh_Quiet_Flag-1].vals.integer == (int)false) {

        printf("Body %d (of %d)\n", bodyIndex+1, numBody);
        printf("Number of nodes    = %d\n", nglobal);
        printf("Number of elements = %d\n", numTri+numQuad);
        //printf("Number of triangle elements      = %d\n", numTri);
        //printf("Number of quadrilateral elements = %d\n", numQuad);
      }

      curveTessInstance->numNodeTotal += nglobal;
      curveTessInstance->numElemTotal += numTri+numQuad;
      bodyIndex++;
    }
  }

  if (restart == 0 &&
      aimInputs[inMesh_Quiet_Flag-1].vals.integer == (int)false) {
    printf("----------------------------\n");
    printf("Total number of nodes    = %d\n", curveTessInstance->numNodeTotal);
    printf("Total number of elements = %d\n", curveTessInstance->numElemTotal);
  }


  for (i = 0; i < curveTessInstance->numMeshRefOut; i++) {
    status = aim_queryMeshes( aimInfo, inMesh_Format, ANALYSISIN, &curveTessInstance->meshRefOut[i] );
    if (status > 0) {
/*@-immediatetrans@*/
      mesh.meshData = NULL;
      mesh.meshRef = &curveTessInstance->meshRefOut[i];
/*@+immediatetrans@*/

      status = mesh_surfaceMeshData(aimInfo, &curveTessInstance->groupMap, &mesh);
      AIM_STATUS(aimInfo, status);

      status = aim_writeMeshes(aimInfo, inMesh_Format, ANALYSISIN, &mesh);
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
/* aimCalcOutput: Calculate/Retrieve Output Information */
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo, int index,
                  capsValue *val)
{
  int i;
  int status = CAPS_SUCCESS;

  aimMesh    mesh;
  aimStorage *curveTessInstance;

  curveTessInstance = (aimStorage *) instStore;


  if (outNumberOfNode == index) {

    val->vals.integer = curveTessInstance->numNodeTotal;

  } else if (outNumberOfElement == index) {

    val->vals.integer = curveTessInstance->numElemTotal;

  } else if (index == outSurface_Mesh) {

    for (i = 0; i < curveTessInstance->numMeshRefOut; i++) {
      status = aim_queryMeshes( aimInfo, outSurface_Mesh, ANALYSISOUT, &curveTessInstance->meshRefOut[i] );
      if (status > 0) {
        /*@-immediatetrans@*/
        mesh.meshData = NULL;
        mesh.meshRef = &curveTessInstance->meshRefOut[i];
        /*@+immediatetrans@*/

        status = mesh_surfaceMeshData(aimInfo, &curveTessInstance->groupMap, &mesh);
        AIM_STATUS(aimInfo, status);

        status = aim_writeMeshes(aimInfo, outSurface_Mesh, ANALYSISOUT, &mesh);
        AIM_STATUS(aimInfo, status);

        status = aim_freeMeshData(mesh.meshData);
        AIM_STATUS(aimInfo, status);
        AIM_FREE(mesh.meshData);
      }
      else
        AIM_STATUS(aimInfo, status);
    }

    // Return the surface meshes
    val->nrow        = curveTessInstance->numMeshRefOut;
    val->vals.AIMptr = curveTessInstance->meshRefOut;

  } else {
    status = CAPS_BADINDEX;
    AIM_STATUS(aimInfo, status, "Unknown output index %d!", index);
  }

  cleanup:

  return status;
}


// ********************** AIM Function Break *****************************
void aimCleanup(void *instStore)
{
    int status; // Returning status
    aimStorage *curveTessInstance;

#ifdef DEBUG
    printf(" curveTessAIM/aimCleanup!\n");
#endif
    curveTessInstance = (aimStorage *) instStore;

    status = destroy_aimStorage(curveTessInstance);
    if (status != CAPS_SUCCESS)
        printf("Error: Status %d during clean up of instance\n", status);

    EG_free(curveTessInstance);
}

