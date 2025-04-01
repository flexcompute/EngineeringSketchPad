/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Cart3D AIM
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

/* libCart3d defaults to WORD_BIT 32 unless told otherwise
 * this makes that assumption explicit and suppresses the warning
 */
#define WORD_BIT 32

#ifdef S_SPLINT_S
#define __INT64_H_
#  define INT64     unsigned long long int
#endif
#include <geomStructures.h>
#undef DOT

int writeSurfTrix(const p_tsTriangulation p_config, const int nComps,
                  const char *const p_fileName, const int options);

#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>          /* Needed in some systems for FLT_MAX definition */
#include <sys/stat.h> // chmod

#include "xddm.h"
#include "xddmUtils.h"

#include "writeTrix.h"
#include "aimUtil.h"
#include "meshUtils.h"
#include "miscUtils.h"
#include "cfdUtils.h"

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define CART3D_BBOX "cart3d_BBox"
#define CART3D_SENSOR "cart3d_Sensor"
#define CART3D_EQAREA "cart3d_EqArea"


//#define DEBUG

enum aimInputs
{
  Tess_Params = 1,                   /* index is 1-based */
  mesh2d,
  outer_box,
  nDiv,
  maxR,
  preSpec,
  symmX,
  symmY,
  symmZ,
  halfBody,
  Boundary_Condition,
  Mach,
  alpha,
  beta,
  Gamma,
  TargetCL,
  TargetCL_Tol,
  TargetCL_Start_Iter,
  TargetCL_Freq,
  TargetCL_NominalAlphaStep,
  Pressure_Scale_Factor,
  maxCycles,
  nMultiGridLevels,
  MultiGridCycleType,
  MultiGridPreSmoothing,
  MultiGridPostSmoothing,
  RK,
  CFL,
  Limiter,
  FluxFun,
  iForce,
  iHist,
  nOrders,
  nAdaptCycles,
  Adapt_Functional,
  Design_Variable,
  Design_Functional,
  Design_Sensitivity,
  Design_Adapt,
  Design_Run_Config,
  Design_Gradient_Memory_Budget,
  Xslices,
  Yslices,
  Zslices,
  y_is_spanwise,
  Model_X_axis,
  Model_Y_axis,
  Model_Z_axis,
  Restart,
  aerocsh,
  aerocsh_arg,
  NUMINPUT = aerocsh_arg             /* Total number of inputs */
};

#define NUMOUT    13                 /* number of outputs */


typedef struct {

  // Functional used for mesh adaptation
  cfdDesignFunctionalStruct *adaptFunctional;

  // Design information
  cfdDesignStruct design;

  // capsGroup to index mapping
  mapAttrToIndexStruct groupMap;

  // Boundary/surface properties
  cfdBoundaryConditionStruct bcProps;

  // Key indicating if aero.csh should be restarted
  const char *aero_start;

  int ntess;
  ego *tess;
} aimStorage;

/*!\mainpage Introduction
 *
 * \section overviewCART3D Cart3D AIM Overview
 *
 * This documentation contains four sections to document the use of the CART3D AIM.
 * \ref examplesCART3D contains example *.csm input files and pyCAPS scripts designed to make use of the CART3D AIM.  These example
 * scripts make extensive use of the \ref attributeCART3D and Cart3D \ref aimInputsCART3D and Cart3D \ref aimOutputsCART3D.
 *
 * The Cart3D AIM can automatically execute Cart3D, with details provided in \ref aimExecuteCART3D.
 *
 * Details of the AIM's automated data transfer capabilities are outlined in \ref dataTransferCART3D
 *
 * \section dependeciesCART3D Dependencies
 * ESP client of libxddm. For XDDM documentation, see $CART3D/doc/xddm/xddm.html. The library uses XML Path Language (XPath) to
 * navigate the elements of XDDM documents. For XPath tutorials, see the web,
 * e.g. www.developer.com/net/net/article.php/3383961/NET-and-XML-XPath-Queries.htm
 *
 * Dependency: libxml2, www.xmlsoft.org. This library is usually present on
 * most systems, check existence of 'xml2-config' script
 */

/*! \page examplesCART3D Cart3D Examples
 *
 * This example contains a set of *.csm and pyCAPS (*.py) inputs that uses the Cart3D AIM.  A user should have knowledge
 * on the generation of parametric geometry in Engineering Sketch Pad (ESP) before attempting to integrate with any AIM.
 * Specifically, this example makes use of Design Parameters, Set Parameters, User Defined Primitive (UDP) and attributes
 * in ESP.
 *
 */

/*! \page attributeCART3D Cart3D attributes
 * The following list of attributes drives the Cart3D geometric definition.
 *
 * \subsection bodyattr Body level attributes
 * - <b> capsAIM</b> This attribute is a CAPS requirement to indicate the analysis the geometry
 * representation supports.
 *
 * - <b> capsReferenceArea</b>  This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the Reference_Area entry in the Cart3D input.
 *
 * - <b> capsReferenceChord</b>  This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the Reference_Length entry in the Cart3D input.
 *
 * - <b> capsReferenceX</b>  This attribute may exist on any <em> Body</em>.  Its
 * value will be used in the Moment_Point entry in the Cart3D input.
 *
 * - <b> capsReferenceY</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used in the Moment_Point entry in the Cart3D input.
 *
 * - <b> capsReferenceZ</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 *  value will be used in the Moment_Point entry in the Cart3D input.
 *
 * - <b> cart3d_BBox </b> [Optional] This attribute may exist on any <em> Body</em>.
 *  If present, the bounding box of the body is used to define a BBox region in the preSpec.c3d.cntl file (the body it self is ignored).
 *  The BBox refinement level is given by the cart3d_BBox value. Note: This will replace the preSpec.c3d.cntl
 *  file generated by autoInputs.
 *
 * - <b> cart3d_Sensor </b> [Optional] This attribute may exist on any Line Wire or Node <em> Body</em>.
 *  Indicates the body is used to construct a Line or Point sensor. The string value of the attribute is the name of the sensor.<br>
 *  Note: The sensor name must be unique from all other Sensor names and capsGroup names.
 *  - <b> cart3d_EqArea </b> [Optional] Additional attribute to convert a line sensor to an equivalent area. The value is the radius.
 *
 * \subsection faceattr Face level attributes
 * - <b> capsGroup</b> The capsGroup attribute defines the geometric components. The loads for each component are available as Dynamic Outputs.
 */

static const char* const runDir = "design";


// Write Cart3D data transfer files
static int
cart3d_dataTransfer(void *aimInfo,
                    ego tess,
                    double **dxyz_out)
{

    /*! \page dataTransferCART3D Cart3D Data Transfer
     *
     * \section dataToCart3D Data transfer to Cart3D (FieldIn)
     *
     * <ul>
     *  <li> <B>"Displacement"</B> </li> <br>
     *   Retrieves nodal displacements (as from a structural solver)
     *   and updates Cart3D's surface mesh.
     * </ul>
     */

    int status; // Function return status
    int i, ibound, ibody, iglobal; // Indexing

    // Discrete data transfer variables
    capsDiscr *discr;
    char **boundNames = NULL;
    int numBoundName = 0;
    enum capsdMethod dataTransferMethod;
    int numDataTransferPoint;
    int dataTransferRank;
    double *dataTransferData;
    char *units;

    int state, nGlobal;
    ego body;

    // Data transfer Out variable
    double *dxyz = NULL;

    int foundDisplacement = (int) false;

    *dxyz_out = NULL;

    status = aim_getBounds(aimInfo, &numBoundName, &boundNames);
    if (status == CAPS_NOTFOUND) return CAPS_SUCCESS;
    AIM_STATUS(aimInfo, status);

    foundDisplacement = (int) false;
    for (ibound = 0; ibound < numBoundName; ibound++) {
      AIM_NOTNULL(boundNames, aimInfo, status);

      status = aim_getDiscr(aimInfo, boundNames[ibound], &discr);
      if (status != CAPS_SUCCESS) continue;

      status = aim_getDataSet(discr,
                              "Displacement",
                              &dataTransferMethod,
                              &numDataTransferPoint,
                              &dataTransferRank,
                              &dataTransferData,
                              &units);
      if (status != CAPS_SUCCESS) continue;

      foundDisplacement = (int) true;

      // Is the rank correct?
      if (dataTransferRank != 3) {
        AIM_ERROR(aimInfo, "Displacement transfer data found however rank is %d not 3!!!!", dataTransferRank);
        status = CAPS_BADRANK;
        goto cleanup;
      }

    } // Loop through bound names

    if (foundDisplacement != (int) true ) {
        status = CAPS_SUCCESS;
        goto cleanup;
    }

    status  = EG_statusTessBody(tess, &body, &state, &nGlobal);
    AIM_STATUS(aimInfo, status);

    AIM_ALLOC(dxyz, 3*nGlobal, double, aimInfo, status);
    memset(dxyz, 0, 3*nGlobal*sizeof(double));

    // now apply the displacements
    for (ibound = 0; ibound < numBoundName; ibound++) {
      AIM_NOTNULL(boundNames, aimInfo, status);

      status = aim_getDiscr(aimInfo, boundNames[ibound], &discr);
      if (status != CAPS_SUCCESS) continue;

      status = aim_getDataSet(discr,
                              "Displacement",
                              &dataTransferMethod,
                              &numDataTransferPoint,
                              &dataTransferRank,
                              &dataTransferData,
                              &units);
      if (status != CAPS_SUCCESS) continue;

      if (numDataTransferPoint != discr->nPoints &&
          numDataTransferPoint > 1) {
        AIM_ERROR(aimInfo, "Developer error!! %d != %d", numDataTransferPoint, discr->nPoints);
        status = CAPS_MISMATCH;
        goto cleanup;
      }

      // Find the disc tessellation in the original list of tessellations
      for (ibody = 0; ibody < discr->nBodys; ibody++) {
        if (discr->bodys[ibody].tess == tess) {
          break;
        }
      }
      if (ibody == discr->nBodys && numDataTransferPoint > 1) {
        AIM_ERROR(aimInfo, "Could not find matching tessellation!");
        status = CAPS_MISMATCH;
        goto cleanup;
      }

      for (i = 0; i < discr->nPoints; i++) {

        if (ibody+1 != discr->tessGlobal[2*i+0]) {
          AIM_ERROR(aimInfo, "Delveoper Error! Inconsistent body index %d != %d", ibody+1, discr->tessGlobal[2*i+0]);
          status = CAPS_SOURCEERR;
          goto cleanup;
        }
        iglobal = discr->tessGlobal[2*i+1];

        if (numDataTransferPoint == 1) {
          // A single point means this is an initialization phase
          dxyz[3*(iglobal-1)  ] = dataTransferData[0];
          dxyz[3*(iglobal-1)+1] = dataTransferData[1];
          dxyz[3*(iglobal-1)+2] = dataTransferData[2];
        } else {
          // Apply delta displacements
          dxyz[3*(iglobal-1)+0] = dataTransferData[3*i+0];
          dxyz[3*(iglobal-1)+1] = dataTransferData[3*i+1];
          dxyz[3*(iglobal-1)+2] = dataTransferData[3*i+2];
        }
      }
    } // Loop through bound names

    *dxyz_out = dxyz;
    status = CAPS_SUCCESS;

// Clean-up
cleanup:
/*@-dependenttrans@*/
    if (status != CAPS_SUCCESS) AIM_FREE(dxyz);
/*@+dependenttrans@*/
    AIM_FREE(boundNames);

    return status;
}


// Write out all variants of Components.i.tri and
// retrieve geometric sensitivities w.r.t. design variable
// to write runDir/builder.xml
static int
cart3d_Components(void *aimInfo,
                  int numDesignVariable,
       /*@null@*/ cfdDesignVariableStruct designVariable[],
                  const mapAttrToIndexStruct *groupMap,
                  int nBody,
                  ego *tess)
{
  int status; // Function return status

  int i, j, k, m, irow, icol; // Indexing
  int ibody;

  ego body;
  int state, nface, plen;
  const int *tris, *tric, *ptype, *pindex;
  const double *points, *uv;

  int iov, iot;

  int ndvar = 0, filesExist = (int)true;

  char aimFile[PATH_MAX];

  int nPoint;
  double *dxyz = NULL;

  int geomIndex;

  char tmp[128], name[128], geomDir[128];

  int nvert=0, opts = 0;
  int nnode, nedge, ntri;
  double *xyzs=NULL;
  int *triFaceConn = NULL;
  int *triFaceCompID = NULL;
  int *triFaceTopoID = NULL;
  int *bndEdgeConn = NULL;
  int *bndEdgeCompID = NULL;
  int *bndEdgeTopoID = NULL;
  int twoDMesh = (int)false;
  int numQuadFace;
  int *quadFaceConn = NULL;
  int *quadFaceCompID = NULL;
  int *quadFaceTopoID = NULL;

  int *nodeConn = NULL;
  int *nodeCompID = NULL;
  int *nodeTopoID = NULL;

  capsValue *geomInVal;

  p_tsXddm p_xddm = NULL;
  p_tsTriangulation p_surf=NULL;

  // Allocate xddm
  p_xddm = xddm_new(1);
  AIM_NOTNULL(p_xddm, aimInfo, status);

  // Set the parent "Model" name
  AIM_NOTNULL(p_xddm->p_e, aimInfo, status);
  AIM_STRDUP(p_xddm->p_e->p_nn, "Model", aimInfo, status);

  xddm_addAttribute("ID", "CAPSmodel", &p_xddm->p_e->nAttr, &p_xddm->p_e->p_attr);

  // Determine number of geometry input variables
  for (i = 0; i < numDesignVariable; i++) {

    AIM_NOTNULL(designVariable, aimInfo, status);
    if (designVariable[i].type != DesignVariableGeometry) continue;

    // Get the value index
    geomIndex = aim_getIndex(aimInfo, designVariable[i].name, GEOMETRYIN);

    if(aim_getGeomInType(aimInfo, geomIndex) != 0) {
      AIM_ERROR(aimInfo, "'%s' is not a DESPMTR - can't get sensitivity",
                designVariable[i].name);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    status = aim_getValue(aimInfo, geomIndex, GEOMETRYIN, &geomInVal);
    AIM_STATUS(aimInfo, status);

    // account for the length of the DESPMTR
    ndvar += geomInVal->length;
  }

  status = c3d_newTriangulation(&p_surf, 0, 1);
  if (status != 0 || p_surf == NULL) {
    AIM_ERROR(aimInfo, "c3d_newTriangulation failed\n");
    status = EGADS_MALLOC;
    goto cleanup;
  }


  // Allocate variables
  p_xddm->a_v = xddm_newVariable(ndvar);
  p_xddm->nv = ndvar;
  for (i = j = 0; i < numDesignVariable; i++) {

    AIM_NOTNULL(designVariable, aimInfo, status);
    if (designVariable[i].type != DesignVariableGeometry) continue;

    geomIndex = aim_getIndex(aimInfo, designVariable[i].name, GEOMETRYIN);

    status = aim_getValue(aimInfo, geomIndex, GEOMETRYIN, &geomInVal);
    AIM_STATUS(aimInfo, status);

    // Set the variable name
    for (irow = 0; irow < geomInVal->nrow; irow++) {
      for (icol = 0; icol < geomInVal->ncol; icol++) {

        snprintf(tmp, 128, "%s_%d_%d_", designVariable[i].name, irow+1, icol+1);
        AIM_STRDUP(p_xddm->a_v[j].p_e->p_id, tmp, aimInfo, status);
        AIM_STRDUP(p_xddm->a_v[j].p_e->p_nn, "Variable", aimInfo, status);

        if (geomInVal->length == 1)
          p_xddm->a_v[j].val = geomInVal->vals.real;
        else
          p_xddm->a_v[j].val = geomInVal->vals.reals[irow*geomInVal->ncol + icol];

        if (geomInVal->limits.dlims[0] != geomInVal->limits.dlims[1]) {
          p_xddm->a_v[j].minVal = geomInVal->limits.dlims[0];
          p_xddm->a_v[j].maxVal = geomInVal->limits.dlims[1];
        }
        j++;
      }
    }
  }


  for (i = 0; i < p_xddm->p_e->nAttr; i++) {
    if (strcmp(p_xddm->p_e->p_attr[i].p_name, "ID") == 0) {
      strcpy(p_surf->geomName, p_xddm->p_e->p_attr[i].p_value);
      break;
    }
  }

  /* linearization lives at the verts */
  status = c3d_allocVertData(&p_surf, 1);
  if (status != 0) {
    AIM_ERROR(aimInfo, "c3d_allocVertData failed");
    status = EGADS_MALLOC;
    goto cleanup;
  }
  p_surf->p_vertData[0].dim    = 3;
  p_surf->p_vertData[0].offset = i*3;
  p_surf->p_vertData[0].type   = VTK_Float64;
  p_surf->p_vertData[0].info   = TRIX_shapeLinearization;
  strcpy(p_surf->p_vertData[0].name, "None");

  status = c3d_allocTriData(&p_surf, 2);             /* 2 tags; body and capsGroup IDs */
  if (status != 0) {
    AIM_ERROR(aimInfo, "c3d_allocTriData failed\n");
    status = EGADS_MALLOC;
    goto cleanup;
  }

  strcpy(p_surf->p_triData[0].name, "IntersectComponents");
  p_surf->p_triData[0].dim    = 1;
  p_surf->p_triData[0].offset = 0;
  p_surf->p_triData[0].type   = VTK_Int16;
  p_surf->p_triData[0].info   = TRIX_componentTag;
  strcpy(p_surf->p_triData[1].name, "GMPtags");
  p_surf->p_triData[1].dim    = 1;
  p_surf->p_triData[1].offset = 1;
  p_surf->p_triData[1].type   = VTK_Int16;
  p_surf->p_triData[1].info   = TRIX_componentTag;

  //p_surf->infoCode += DP_VERTS_CODE;

  p_surf->nVerts = 0;
  p_surf->nTris  = 0;

  for (ibody = 0; ibody < nBody; ibody++) {
    status  = EG_statusTessBody(tess[ibody], &body, &state, &nvert);
    if (status != EGADS_SUCCESS) goto cleanup;
    p_surf->nVerts += nvert;

    status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) goto cleanup;

    for (i = 1; i <= nface; i++) {
      status = EG_getTessFace(tess[ibody], i, &plen, &points, &uv, &ptype, &pindex,
                              &ntri, &tris, &tric);
      if (status != EGADS_SUCCESS) goto cleanup;
      p_surf->nTris += ntri;
    }
  }

  status = c3d_allocTriangulation(&p_surf);
  if (status != 0) {
    AIM_ERROR(aimInfo, "c3d_allocTriangulation failed\n");
    status = EGADS_MALLOC;
    goto cleanup;
  }

  iov = 0;
  iot = 0;
  for (ibody = 0; ibody < nBody; ibody++) {

    // Get tessellation
    status = mesh_bodyTessellation(aimInfo, tess[ibody], groupMap,
                                   &nvert, &xyzs,
                                   &ntri, &triFaceConn, &triFaceCompID, &triFaceTopoID,
                                   &nedge, &bndEdgeConn, &bndEdgeCompID, &bndEdgeTopoID,
                                   &nnode, &nodeConn, &nodeCompID, &nodeTopoID,
                                   twoDMesh,
                                   NULL,
                                   &numQuadFace, &quadFaceConn, &quadFaceCompID, &quadFaceTopoID);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(xyzs, aimInfo, status);
    AIM_NOTNULL(triFaceConn, aimInfo, status);
    AIM_NOTNULL(triFaceCompID, aimInfo, status);


    //See if we have data transfer information
    status = cart3d_dataTransfer(aimInfo,
                                 tess[ibody],
                                 &dxyz);
    AIM_STATUS(aimInfo, status);


    for (i = 0; i < nvert; i++) {              /* save vert locations */
      p_surf->a_Verts[iov+i].x[0] = (float) xyzs[3*i+0];
      p_surf->a_Verts[iov+i].x[1] = (float) xyzs[3*i+1];
      p_surf->a_Verts[iov+i].x[2] = (float) xyzs[3*i+2];
      //surf->a_dpVerts[i].x[0] = xyzs[3*i  ];
      //surf->a_dpVerts[i].x[1] = xyzs[3*i+1];
      //surf->a_dpVerts[i].x[2] = xyzs[3*i+2];

      if (dxyz != NULL) {
        p_surf->a_Verts[iov+i].x[0] += dxyz[3*i+0];
        p_surf->a_Verts[iov+i].x[1] += dxyz[3*i+1];
        p_surf->a_Verts[iov+i].x[2] += dxyz[3*i+2];
      }

      p_surf->a_scalar0[3*iov+3*i+0] = 0;
      p_surf->a_scalar0[3*iov+3*i+1] = 0;
      p_surf->a_scalar0[3*iov+3*i+2] = 0;
    }

    for (i = 0; i < ntri; i++) {
      p_surf->a_Tris[iot+i].vtx[0] = iov + triFaceConn[3*i+0] - 1;
      p_surf->a_Tris[iot+i].vtx[1] = iov + triFaceConn[3*i+1] - 1;
      p_surf->a_Tris[iot+i].vtx[2] = iov + triFaceConn[3*i+2] - 1;

      /* -- component numbers -- */
      p_surf->a_scalar0_t[iot+i]               = ibody;
      p_surf->a_scalar0_t[iot+p_surf->nTris+i] = triFaceCompID[i];
    }

    iov += nvert;
    iot += ntri;

    AIM_FREE(xyzs);
    AIM_FREE(dxyz);
    AIM_FREE(triFaceConn);
    AIM_FREE(triFaceCompID);
    AIM_FREE(triFaceTopoID);
    AIM_FREE(bndEdgeConn);
    AIM_FREE(bndEdgeCompID);
    AIM_FREE(bndEdgeTopoID);
    AIM_FREE(quadFaceConn);
    AIM_FREE(quadFaceCompID);
    AIM_FREE(quadFaceTopoID);

    AIM_FREE(nodeConn);
    AIM_FREE(nodeCompID);
    AIM_FREE(nodeTopoID);
  }

  status = aim_file(aimInfo, "inputs/Components.i.tri", aimFile);
  AIM_STATUS(aimInfo, status);

  AIM_NOTNULL(p_surf, aimInfo, status);
  status = writeSurfTrix(p_surf, 1, aimFile, opts);
  if (status != 0) {
    AIM_ERROR(aimInfo, "io_writeSurfTrix failed: %d", status);
    status = EGADS_WRITERR;
    goto cleanup;
  }

  if (aim_isDir(aimInfo, runDir) == CAPS_SUCCESS) {
    // Write runDir/builder.xml
    snprintf(tmp, 128, "%s/builder.xml", runDir);
    status = aim_file(aimInfo, tmp, aimFile);
    AIM_STATUS(aimInfo, status);

    status = xddm_writeFile(aimFile, p_xddm, 0);
    AIM_STATUS(aimInfo, status);
  }

  // Only write out geometric sensitivities if necessary
  if (ndvar == 0) {
    status = CAPS_SUCCESS;
    goto cleanup;
  }

  snprintf(geomDir, 128, "%s/geometry_CAPSmodel", runDir);
  status = aim_mkDir(aimInfo, geomDir);
  AIM_STATUS(aimInfo, status);

  strcat(geomDir, "/geometry");
  status = aim_mkDir(aimInfo, geomDir);
  AIM_STATUS(aimInfo, status);

  /* check if all sensitivity files exist and geometry is old
   * (no need to re-create sensitivities) */
  filesExist = (int)true;
  for (i = m = 0; i < numDesignVariable && filesExist == (int)true; i++) {

    AIM_NOTNULL(designVariable, aimInfo, status);
    if (designVariable[i].type != DesignVariableGeometry) continue;

    // Get the value index
    geomIndex = aim_getIndex(aimInfo, designVariable[i].name, GEOMETRYIN);
    status = aim_getValue(aimInfo, geomIndex, GEOMETRYIN, &geomInVal);
    AIM_STATUS(aimInfo, status);

    for (irow = 0; irow < geomInVal->nrow; irow++) {
      for (icol = 0; icol < geomInVal->ncol; icol++) {

        snprintf(name, 128, "%s%s", "Model__CAPSmodel__Variable__", p_xddm->a_v[m].p_e->p_id);

        snprintf(tmp, 128, "%s/%s/Components.i.tri", geomDir, name);
        if (aim_isFile(aimInfo, tmp) != CAPS_SUCCESS)
          filesExist = (int)false;
        m++;
      }
    }
  }

  if (filesExist == (int)true &&
      aim_newGeometry(aimInfo) != CAPS_SUCCESS) {
    status = CAPS_SUCCESS;
    goto cleanup;
  }

  // Loop over the geometry in values
  for (i = m = 0; i < numDesignVariable; i++) {

    AIM_NOTNULL(designVariable, aimInfo, status);
    if (designVariable[i].type != DesignVariableGeometry) continue;

    geomIndex = aim_getIndex(aimInfo, designVariable[i].name, GEOMETRYIN);

    status = aim_getValue(aimInfo, geomIndex, GEOMETRYIN, &geomInVal);
    AIM_STATUS(aimInfo, status);

    for (irow = 0; irow < geomInVal->nrow; irow++) {
      for (icol = 0; icol < geomInVal->ncol; icol++) {

        iov = 0;
        for (ibody = 0; ibody < nBody; ibody++) {

          status = aim_tessSensitivity(aimInfo,
                                       designVariable[i].name,
                                       irow+1, icol+1, // row, col
                                       tess[ibody],
                                       &nPoint, &dxyz);
          AIM_STATUS(aimInfo, status, "Sensitivity for: %s", designVariable[ibody].name);
          AIM_NOTNULL(dxyz, aimInfo, status);

          // Set the variable name this this body
          snprintf(p_surf->p_vertData[0].name, STRING_LEN, "%s%s", "Model__CAPSmodel__Variable__", p_xddm->a_v[m].p_e->p_id);

          /* actual sensitivities */
          for (k = 0; k < nPoint; k++) {
            p_surf->a_scalar0[3*iov+3*k  ] = dxyz[3*k  ];
            p_surf->a_scalar0[3*iov+3*k+1] = dxyz[3*k+1];
            p_surf->a_scalar0[3*iov+3*k+2] = dxyz[3*k+2];
          }

          iov += nPoint;
          AIM_FREE(dxyz);
        }

        // Create the directory
        snprintf(tmp, 128, "%s/%s", geomDir, p_surf->p_vertData[0].name);
        status = aim_mkDir(aimInfo, tmp);
        AIM_STATUS(aimInfo, status);

        // and make the file name
        strcat(tmp, "/Components.i.tri");
        status = aim_file(aimInfo, tmp, aimFile);
        AIM_STATUS(aimInfo, status);

        // write the surface trix file
        status = writeSurfTrix(p_surf, 1, aimFile, opts);
        if (status != 0) {
          AIM_ERROR(aimInfo, "io_writeSurfTrix failed: %d", status);
          status = EGADS_WRITERR;
          goto cleanup;
        }
        m++;
      }
    }
  }

  status = CAPS_SUCCESS;

cleanup:
  if (p_surf != NULL) {
    c3d_freeTriangulation(p_surf, 0);
    free(p_surf); /* must use free */
  }

  xddm_free(p_xddm, opts);

  AIM_FREE(dxyz);

  AIM_FREE(xyzs);
  AIM_FREE(triFaceTopoID);
  AIM_FREE(bndEdgeConn);
  AIM_FREE(bndEdgeCompID);
  AIM_FREE(bndEdgeTopoID);
  AIM_FREE(quadFaceConn);
  AIM_FREE(quadFaceCompID);
  AIM_FREE(quadFaceTopoID);

  AIM_FREE(nodeConn);
  AIM_FREE(nodeCompID);
  AIM_FREE(nodeTopoID);

  return status;
}

// Write runDir/design.xml
static int
cart3d_AeroFun(void *aimInfo,
               const char *funName,
               cfdDesignFunctionalStruct *designFunctional,
               p_tsXddmAFun p_afun )
{
  int status; // Function return status

  int i, j; // Indexing

  int force=0, moment=0, frame, funType;

  const char *name;
  const char *names[] = {"C_A"  , "C_Y"  , "C_N",
                         "C_D"  , "C_S"  , "C_L",
                         "C_l"  , "C_m"  , "C_n",
                         "C_M_x", "C_M_y", "C_M_z",
                         "LoD"};

  AIM_STRDUP(p_afun->p_e->p_nn, "AeroFun", aimInfo, status);

  funType = 0;
  for (j = 0; j < designFunctional->numComponent; j++) {

    name = designFunctional->component[j].name;

    if (strcmp("C_A", name) == 0 ||
        strcmp("C_Y", name) == 0 ||
        strcmp("C_N", name) == 0 ||
        strcmp("C_D", name) == 0 ||
        strcmp("C_S", name) == 0 ||
        strcmp("C_L", name) == 0) {

      if (funType != 0 && funType != 1) {
        AIM_ERROR(aimInfo, "Functional components for '%s' must all be Forces!", funName);
        status = CAPS_BADVALUE;
        goto cleanup;
      }
      funType = 1;

      if (       strcmp("C_A", name) == 0 ||
                 strcmp("C_D", name) == 0) {
        force = 0;
      } else if (strcmp("C_Y", name) == 0 ||
                 strcmp("C_S", name) == 0) {
        force = 1;
      } else if (strcmp("C_N", name) == 0 ||
                 strcmp("C_L", name) == 0) {
        force = 2;
      }

      if (strcmp("C_D", name) == 0 ||
          strcmp("C_S", name) == 0 ||
          strcmp("C_L", name) == 0) {
        frame = 0; // Aerodynamic frame
      } else {
        frame = 1; // Aircraft (body)
      }

      status = xddm_addAeroFunForce(p_afun,
                                    name,
                                    force,
                                    frame,
                                    designFunctional->component[j].form,
                                    designFunctional->component[j].power,
                                    designFunctional->component[j].target,
                                    designFunctional->component[j].weight,
                                    0,
                                    designFunctional->component[j].boundaryName);
      AIM_STATUS(aimInfo, status);

    } else if (strcmp("C_l", name) == 0 ||
               strcmp("C_m", name) == 0 ||
               strcmp("C_n", name) == 0 ||
               strcmp("C_M_x", name) == 0 ||
               strcmp("C_M_y", name) == 0 ||
               strcmp("C_M_z", name) == 0) {

      if (funType != 0 && funType != 2) {
        AIM_ERROR(aimInfo, "Functional components for '%s' must all be Moments!", funName);
        status = CAPS_BADVALUE;
        goto cleanup;
      }
      funType = 2;

      if (       strcmp("C_l", name) == 0 ||
                 strcmp("C_M_x", name) == 0) {
        moment = 0;
      } else if (strcmp("C_m", name) == 0 ||
                 strcmp("C_M_y", name) == 0) {
        moment = 1;
      } else if (strcmp("C_n", name) == 0 ||
                 strcmp("C_M_z", name) == 0) {
        moment = 2;
      }

      if (strcmp("C_M_x", name) == 0 ||
          strcmp("C_M_y", name) == 0 ||
          strcmp("C_M_z", name) == 0) {
        frame = 0; // Aerodynamic frame
      } else {
        frame = 1; // Aircraft (body)
      }

      status = xddm_addAeroFunMoment_Point(p_afun,
                                           name,
                                           0, // index
                                           moment,
                                           frame,
                                           designFunctional->component[j].form,
                                           designFunctional->component[j].power,
                                           designFunctional->component[j].target,
                                           designFunctional->component[j].weight,
                                           0,
                                           designFunctional->component[j].boundaryName);
      AIM_STATUS(aimInfo, status);

    } else if (strcmp(name, "LoD") == 0) {

      if (funType != 0 && funType != 3) {
        AIM_ERROR(aimInfo, "Functional components for '%s' must all be LoD!", funName);
        status = CAPS_BADVALUE;
        goto cleanup;
      }
      funType = 3;

      frame = 0; // Aerodynamic frame
      status = xddm_addAeroFunLoD(p_afun,
                                  name,
                                  frame,
                                  designFunctional->component[j].form,
                                  1,
                                  designFunctional->component[j].power,
                                  designFunctional->component[j].bias,
                                  designFunctional->component[j].target,
                                  designFunctional->component[j].weight,
                                  0,
                                  designFunctional->component[j].boundaryName);
      AIM_STATUS(aimInfo, status);

    } else {
      AIM_ERROR(aimInfo, "Unknown function: '%s'", name);
      AIM_ADDLINE(aimInfo, "Case sensitive available functions:");
      for (i = 0; i < sizeof(names)/sizeof(char*); i++)
        AIM_ADDLINE(aimInfo, "%s", names[i]);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  }

  status = CAPS_SUCCESS;

cleanup:
  return status;
}

// Write runDir/design.xml
static int
cart3d_designxml(void *aimInfo,
                 const cfdDesignStruct *design,
                 const int sensitivity,
                 const char *adapt)
{
  int status; // Function return status

  int i, n; // Indexing

  int nAnalysisVar;

  p_tsXddm p_xddm = NULL;

  char aimFile[PATH_MAX], tmp[128];

  const char *name;
  unsigned options = 0;

  p_xddm = xddm_new(1);
  AIM_NOTNULL(p_xddm, aimInfo, status);

  /* Set the parent "Optimize" name */
  AIM_NOTNULL(p_xddm->p_e, aimInfo, status);
  AIM_STRDUP(p_xddm->p_e->p_nn, "Optimize", aimInfo, status);

  /* Use configuration to enable/disable sensitivities */
  if (sensitivity == (int)true) {
    options = XDDM_LIN;
  } else {
    options = XDDM_NOLIN;
  }

  /* Set the geometry intersect */
  p_xddm->a_geo = xddm_newIntersect(1);
  AIM_NOTNULL(p_xddm->a_geo, aimInfo, status);
  p_xddm->ng = 1;
  AIM_STRDUP(p_xddm->a_geo->p_e->p_nn, "Intersect", aimInfo, status);
  AIM_STRDUP(p_xddm->a_geo->p_e->p_id, "CAPSmodel", aimInfo, status);
  xddm_addAttribute("Parts", "builder.xml", &p_xddm->a_geo->p_e->nAttr, &p_xddm->a_geo->p_e->p_attr);

  /* Create a design point */
  // set the design point ID (required)

  p_xddm->a_kids = xddm_new(1);
  AIM_NOTNULL(p_xddm->a_kids, aimInfo, status);
  p_xddm->nk = 1;

  AIM_STRDUP(p_xddm->a_kids->p_e->p_id, "DP1", aimInfo, status);
  AIM_STRDUP(p_xddm->a_kids->p_e->p_nn, "DesignPoint", aimInfo, status);
  xddm_addAttribute("Geometry", "CAPSmodel", &p_xddm->a_kids->p_e->nAttr, &p_xddm->a_kids->p_e->p_attr);

  p_xddm->a_kids->a_ap = xddm_newAnalysis(design->numDesignFunctional);
  AIM_NOTNULL(p_xddm->a_kids->a_ap, aimInfo, status);
  p_xddm->a_kids->na = design->numDesignFunctional;

  for (i = 0; i < design->numDesignFunctional; i++) {
    name = design->designFunctional[i].name;
    AIM_STRDUP(p_xddm->a_kids->a_ap[i].p_e->p_id, name, aimInfo, status);
    AIM_STRDUP(p_xddm->a_kids->a_ap[i].p_e->p_nn, "Analysis", aimInfo, status);

    p_xddm->a_kids->a_ap[i].p_afun = xddm_newAeroFun(1);
    AIM_NOTNULL(p_xddm->a_kids->a_ap[i].p_afun, aimInfo, status);

    if ((adapt != NULL) && (strcasecmp(adapt, name) == 0)) {
      xddm_addAttribute("Options", "Adapt", &p_xddm->a_kids->a_ap[i].p_afun->p_e->nAttr, &p_xddm->a_kids->a_ap[i].p_afun->p_e->p_attr);
    }

    status = cart3d_AeroFun(aimInfo, name,
                            &design->designFunctional[i],
                            p_xddm->a_kids->a_ap[i].p_afun );
    AIM_STATUS(aimInfo, status);
  }

  // Add AnalysisIn design variables

  nAnalysisVar = 0;
  for (i = 0; i < design->numDesignVariable; i++) {
    if (design->designVariable[i].type != DesignVariableAnalysis) continue;
    nAnalysisVar += design->designVariable[i].var->length;
    if (design->designVariable[i].var->length != 1) {
      AIM_ERROR(aimInfo, "Developer error: Design variable '%s' length != 1", design->designVariable[i].name);
      status = CAPS_BADVALUE;
      goto cleanup;
    }
  }

  p_xddm->a_kids->a_v = xddm_newVariable(nAnalysisVar);
  AIM_NOTNULL(p_xddm->a_kids->a_v, aimInfo, status);
  p_xddm->a_kids->nv = nAnalysisVar;

  n = 0;
  for (i = 0; i < design->numDesignVariable; i++) {
    if (design->designVariable[i].type != DesignVariableAnalysis) continue;

    name = design->designVariable[i].name;
    AIM_STRDUP(p_xddm->a_kids->a_v[n].p_e->p_id, name, aimInfo, status);
    AIM_STRDUP(p_xddm->a_kids->a_v[n].p_e->p_nn, "Variable", aimInfo, status);

    if (design->designVariable[i].var->length != 1) {
      AIM_ERROR(aimInfo, "AnalysisIn design variables can only be lenght 1!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    p_xddm->a_kids->a_v[n].val = design->designVariable[i].var->vals.real;

    if (design->designVariable[i].lowerBound[0] != -FLT_MAX)
      p_xddm->a_kids->a_v[n].minVal = design->designVariable[i].lowerBound[0];

    if (design->designVariable[i].upperBound[0] !=  FLT_MAX)
      p_xddm->a_kids->a_v[n].maxVal = design->designVariable[i].upperBound[0];

    if (design->designVariable[i].typicalSize[0] != 0)
      p_xddm->a_kids->a_v[n].typicalSize = design->designVariable[i].typicalSize[0];

    n++;
  }

  // Write design.xml

  snprintf(tmp, 128, "%s/design.xml", runDir);
  status = aim_file(aimInfo, tmp, aimFile);
  AIM_STATUS(aimInfo, status);

  status = xddm_writeFile(aimFile, p_xddm, options);
  AIM_STATUS(aimInfo, status);

cleanup:

  xddm_free(p_xddm, options);

  return status;
}


// Write inputs/Config.xml
static int
cart3d_configxml(void *aimInfo,
                 const mapAttrToIndexStruct *groupMap)
{
  int status; // Function return status

  int i; // Indexing

  unsigned options = 0;
  p_tsXddm p_xddm = NULL;
  char aimFile[PATH_MAX];
  char tmp[128];
  const char *configFile = "inputs/Config.xml";

  if (aim_newGeometry(aimInfo) != CAPS_SUCCESS &&
      aim_isFile(aimInfo, configFile) == CAPS_SUCCESS) return CAPS_SUCCESS;

  p_xddm = xddm_new(1);
  AIM_NOTNULL(p_xddm->p_e, aimInfo, status);

  /* Set the parent "Configuration" name */
  AIM_STRDUP(p_xddm->p_e->p_nn, "Configuration", aimInfo, status);
  xddm_addAttribute("Name", "CAPSmodel", &p_xddm->p_e->nAttr, &p_xddm->p_e->p_attr);
  xddm_addAttribute("Source", "Components.i.tri", &p_xddm->p_e->nAttr, &p_xddm->p_e->p_attr);


  /* Create components (hijack Analysis structure for now...) */
  p_xddm->a_ap = xddm_newAnalysis(groupMap->numAttribute);
  AIM_NOTNULL(p_xddm->a_ap, aimInfo, status);
  p_xddm->na = groupMap->numAttribute;

  for (i = 0; i < groupMap->numAttribute; i++) {

    AIM_STRDUP(p_xddm->a_ap[i].p_e->p_nn, "Component", aimInfo, status);

    xddm_addAttribute("Name", groupMap->attributeName[i], &p_xddm->a_ap[i].p_e->nAttr, &p_xddm->a_ap[i].p_e->p_attr);
    xddm_addAttribute("Type", "tri"                     , &p_xddm->a_ap[i].p_e->nAttr, &p_xddm->a_ap[i].p_e->p_attr);

    p_xddm->a_ap[i].p_afun = xddm_newAeroFun(1);
    AIM_NOTNULL(p_xddm->a_ap[i].p_afun, aimInfo, status);
    AIM_STRDUP(p_xddm->a_ap[i].p_afun->p_e->p_nn, "Data", aimInfo, status);

    snprintf(tmp, 128, "Face Label=%d", groupMap->attributeIndex[i]);
    AIM_STRDUP(p_xddm->a_ap[i].p_afun->p_text, tmp, aimInfo, status);
  }

  // Write Config.xml

  status = aim_file(aimInfo, configFile, aimFile);
  AIM_STATUS(aimInfo, status);

  status = xddm_writeFile(aimFile, p_xddm, options);
  AIM_STATUS(aimInfo, status);

cleanup:

  xddm_free(p_xddm, options);

  return status;
}


// Write inputs/Functionals.xml
static int
cart3d_functionalsxml(void *aimInfo,
                      cfdDesignFunctionalStruct *adaptFunctional)
{
  int status; // Function return status

  p_tsXddm p_xddm = NULL;
  char aimFile[PATH_MAX];

  unsigned options = 0;

  p_xddm = xddm_new(1);

  /* Set the parent "Functionals" name */
  AIM_NOTNULL(p_xddm->p_e, aimInfo, status);
  AIM_STRDUP(p_xddm->p_e->p_nn, "Functionals", aimInfo, status);

  /* Create one AeroFun */
  p_xddm->a_afun = xddm_newAeroFun(1);
  AIM_NOTNULL(p_xddm->a_afun, aimInfo, status);
  p_xddm->naf = 1;

  /* Specify this functional is for adaptation */
  //AIM_STRDUP(p_xddm->a_afun->p_e, "Adapt", aimInfo, status);

  AIM_STRDUP(p_xddm->a_afun->p_e->p_id, adaptFunctional->name, aimInfo, status);
  xddm_addAttribute("Options", "Adapt", &p_xddm->a_afun->p_e->nAttr, &p_xddm->a_afun->p_e->p_attr);

  status = cart3d_AeroFun(aimInfo, adaptFunctional->name,
                          adaptFunctional, p_xddm->a_afun);
  AIM_STATUS(aimInfo, status);

  // Write Functionals.xml

  status = aim_file(aimInfo, "inputs/Functionals.xml", aimFile);
  AIM_STATUS(aimInfo, status);

  status = xddm_writeFile(aimFile, p_xddm, options);
  AIM_STATUS(aimInfo, status);

cleanup:

  xddm_free(p_xddm, options);

  return status;
}



// Write inputs/aero.csh
static int
cart3d_aerocsh(void *aimInfo,
               int limiter,
               int mesh2d,
               int nAdapt,
               int maxR,
               int preSpec,
               int y_is_spanwise,
               int aerocsh_len,
               const char *aerocsh_str)
{
  int status = CAPS_SUCCESS; // Function return status

  int i; // Indexing

  char aimFile[PATH_MAX];

  size_t linecap=0;
  char *line=NULL;

  const char *CART3D;
  const char *header = "# STOP: no user specified";
  FILE *fin=NULL, *fout=NULL;

  CART3D = getenv("CART3D");
  if (CART3D == NULL) {
    AIM_ERROR(aimInfo, "CART3D environment variable must point to the Cart3D installation!");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Input aero.csh
  snprintf(aimFile, PATH_MAX, "%s/bin/aero.csh", CART3D);
  fin = fopen(aimFile, "r");
  if (fin == NULL) {
    AIM_ERROR(aimInfo, "Failed to open: %s", aimFile);
    status = CAPS_IOERR;
    goto cleanup;
  }

  // Output aero.csh

  fout = aim_fopen(aimInfo, "inputs/aero.csh", "w");
  if (fout == NULL) {
    AIM_ERROR(aimInfo, "Failed to open: %s", aimFile);
    status = CAPS_IOERR;
    goto cleanup;
  }

  // Insert CAPS inputs after: '# STOP: no user specified'
  while (getline(&line, &linecap, fin) != -1) {
      AIM_NOTNULL(line, aimInfo, status);
      i = 0;
      while (line[i] == ' ' && line[i] != '\0') i++;
      if (strncasecmp(header, &line[i], strlen(header)) == 0) {

        fprintf(fout, "# CAPS aerocsh inputs\n");

        for (i = 0; i < aerocsh_len; i++) {
          fprintf(fout, "%s\n", aerocsh_str);
          aerocsh_str += strlen(aerocsh_str)+1;
        }

        fprintf(fout, "set limiter = %d\n", limiter);
        fprintf(fout, "set mesh2d = %d\n", mesh2d);
        fprintf(fout, "set n_adapt_cycles = %d\n", nAdapt);
        fprintf(fout, "set maxR = %d\n", maxR);
        fprintf(fout, "set use_preSpec = %d\n", preSpec);
        fprintf(fout, "set y_is_spanwise = %d\n", y_is_spanwise);
        fprintf(fout, "# ---------------------------------------------\n");

      }

      // write the line
      fprintf(fout, "%s", line);
  }

  fclose(fout); fout = NULL;

  // Make sure aero.csh has execute permission
  status = aim_file(aimInfo, "inputs/aero.csh", aimFile);
  AIM_STATUS(aimInfo, status);
  chmod(aimFile, S_IXUSR | S_IRUSR | S_IWUSR);

  status = CAPS_SUCCESS;

cleanup:
/*@-dependenttrans*/
  if (fin  != NULL) fclose(fin);
  if (fout != NULL) fclose(fout);
/*@+dependenttrans*/

  if (line != NULL) free(line);

  return status;
}


// Fill bcProps in a cfdBCStruct format with boundary condition information from incoming BC Tuple
int cart3d_getBoundaryCondition(void *aimInfo,
                                int numTuple,
                                capsTuple bcTuple[],
                                mapAttrToIndexStruct *attrMap,
                                cfdBoundaryConditionStruct *bcProps) {

    /*! \page cart3dBoundaryConditions CFD Boundary Conditions
     * Structure for the boundary condition tuple  = ("CAPS Group Name", "Value").
     * "CAPS Group Name" defines the capsGroup on which the boundary condition should be applied.
     *  The "Value" can either be a JSON String dictionary (see Section \ref jsonStringCFDBoundary) or a single string keyword string
     *  (see Section \ref keyStringCFDBoundary)
     */

    int status; //Function return

    int i, j, iprop; // Indexing

    int bcIndex;

    char *keyValue = NULL;
    char *keyWord = NULL;

    printf("Getting Cart3D boundary conditions\n");

    //printf("Number of tuple pairs = %d\n", numTuple);

    for (i = 0; i < numTuple; i++) {

        printf("\tBoundary condition name - %s\n", bcTuple[i].name);

        status = get_mapAttrToIndexIndex(attrMap, bcTuple[i].name, &bcIndex);
        if (status == CAPS_NOTFOUND) {
            AIM_ANALYSISIN_ERROR(aimInfo, Boundary_Condition, "BC name \"%s\" not found in capsGroup attributes", bcTuple[i].name);
            AIM_ADDLINE(aimInfo, "Available capsGroup attributes:");
            for (j = 0; j < attrMap->numAttribute; j++) {
                AIM_ADDLINE(aimInfo, "\t%s", attrMap->attributeName[j]);
            }
            return status;
        }

        // Find the boundary
        iprop = 0;
        while (bcProps->surfaceProp[iprop].bcID != bcIndex) {
          iprop++;
          if (iprop == bcProps->numSurfaceProp) {
            AIM_ERROR(aimInfo, "Developer error!");
            status = CAPS_NOTIMPLEMENT;
            goto cleanup;
          }
        }

        // bcIndex is 1 bias coming from attribute mapper
        //bcIndex = bcIndex -1;

        // Do we have a json string?
        if (strncmp(bcTuple[i].value, "{", 1) == 0) {

            //printf("\t\tJSON String - %s\n",bcTuple[i].value);

            /*! \page cart3dBoundaryConditions
             * \section jsonStringCFDBoundary JSON String Dictionary
             *
             * If "Value" is a JSON string dictionary
             * \if (FUN3D || SU2)
             *  (eg. "Value" = {"bcType": "Viscous", "wallTemperature": 1.1})
             * \endif
             *  the following keywords ( = default values) may be used:
             *
             * <ul>
             * <li> <B>bcType = ""</B> </li> <br>
             *      Boundary condition type. Options:
             *
             *  - Riemann             (Generic INLET/EXHAUST BC -- 2004)
             *  - InletPressRatio
             *  - InletVelocity
             *  - PowerStag
             *  - PowerMFR
             * </ul>
             */

            // Get property Type
            keyWord = "bcType";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {
                AIM_NOTNULL(keyValue, aimInfo, status);

                //{UnknownBoundary, Inviscid, Viscous, Farfield, Extrapolate, Freestream,
                // BackPressure, Symmetry, SubsonicInflow, SubsonicOutflow,
                // MassflowIn, MassflowOut, FixedInflow, FixedOutflow, MachOutflow}
                if      (strcasecmp(keyValue, "\"Riemann\"")         == 0) bcProps->surfaceProp[iprop].surfaceType = Farfield;
                else if (strcasecmp(keyValue, "\"InletPressRatio\"") == 0) bcProps->surfaceProp[iprop].surfaceType = BackPressure;
                else if (strcasecmp(keyValue, "\"InletVelocity\"")   == 0) bcProps->surfaceProp[iprop].surfaceType = SubsonicOutflow;
                else if (strcasecmp(keyValue, "\"PowerStag\"")       == 0) bcProps->surfaceProp[iprop].surfaceType = SubsonicInflow;
                else if (strcasecmp(keyValue, "\"PowerMFR\"")        == 0) bcProps->surfaceProp[iprop].surfaceType = MassflowIn;
                else {

                    AIM_ERROR(aimInfo, "Unrecognized \"%s\" specified (%s) for Boundary_Condition tuple %s, current options:"
                            "\" Riemann, InletPressRatio, InletVelocity, PowerStag, PowerMFR"
                            "\"",
                            keyWord, keyValue, bcTuple[i].name);
                    status = CAPS_NOTFOUND;
                    goto cleanup;
                }

            } else {
                AIM_ERROR(aimInfo, "\tNo \"%s\" specified for tuple %s in json string", keyWord,
                        bcTuple[i].name);
                status = CAPS_BADVALUE;
                goto cleanup;
            }

            AIM_FREE(keyValue);

            // Wall specific properties

            // Stagnation quantities

            /*! \page cart3dBoundaryConditions
             *
             * <ul>
             * <li>  <B>totalPressure = 0.0</B> </li> <br>
             *  Total pressure (Pt/Ptinf) for <B>PowerStag</B>
             * </ul>
             *
             */
            keyWord = "totalPressure";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].totalPressure);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cart3dBoundaryConditions
             *
             * <ul>
             * <li>  <B>totalTemperature = 0.0</B> </li> <br>
             *  Total temperature (Tt/Ttinf) for <B>PowerStag</B> and <B>PowerMFR</B>
             * </ul>
             */
            keyWord = "totalTemperature";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].totalTemperature);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            // Static quantities

            /*! \page cart3dBoundaryConditions
             *
             *  <ul>
             *  <li>  <B>staticDensity = 0.0</B> </li> <br>
             *  Static density for <B>Riemann</B>
             *  </ul>
             */
            keyWord = "staticDensity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].staticDensity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cart3dBoundaryConditions
             *
             *  <ul>
             *  <li> <B>staticPressure = 0.0</B> </li> <br>
             *  Static pressure for <B>Riemann</B>, or P/P_inf for <B>InletPressRatio</B>
             *  </ul>
             */
            keyWord = "staticPressure";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].staticPressure);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }


            // Velocity components
            /*! \page cart3dBoundaryConditions
             *
             * <ul>
             *  <li>  <B>uVelocity = 0.0</B> </li> <br>
             *  X-velocity component for <B>Riemann</B>
             * </ul>
             */
            keyWord = "uVelocity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].uVelocity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cart3dBoundaryConditions
             *
             *  <ul>
             *  <li>  <B>vVelocity = 0.0</B> </li> <br>
             *  Y-velocity component for <B>Riemann</B>
             *  </ul>
             *
             */
            keyWord = "vVelocity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].vVelocity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cart3dBoundaryConditions
             *
             *  <ul>
             *  <li>  <B>wVelocity = 0.0</B> </li> <br>
             *  Z-velocity component for <B>Riemann</B>
             * </ul>
             */
            keyWord = "wVelocity";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].wVelocity);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            /*! \page cart3dBoundaryConditions
             *
             *  <ul>
             *  <li> <B>machNumber = 0.0</B> </li> <br>
             *  Mach number (Vnorm/c_inf) for <B>InletVelocity</B>
             * </ul>
             */
            keyWord = "machNumber";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].machNumber);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

            // Massflow
            /*! \page cart3dBoundaryConditions
             *
             *  <ul>
             *  <li> <B>massflow = 0.0</B> </li> <br>
             *  Massflow for <B>PowerMFR</B>
             *  </ul>
             */
            keyWord = "massflow";
            status = search_jsonDictionary( bcTuple[i].value, keyWord, &keyValue);
            if (status == CAPS_SUCCESS) {

                status = string_toDouble(keyValue, &bcProps->surfaceProp[iprop].massflow);
                AIM_STATUS(aimInfo, status);
                AIM_FREE(keyValue);
            }

        } else {
            AIM_ERROR(aimInfo, "BC data for '%s' must be a json string!\n",
                      bcProps->surfaceProp[i].name);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    printf("\tDone getting Cart3D boundary conditions\n");

    status = CAPS_SUCCESS;
cleanup:
    AIM_FREE(keyValue);
    return status;
}


/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int    *ints = NULL, i, status = CAPS_SUCCESS;
  char   **strs = NULL;
  aimStorage *cartInstance = NULL;

#ifdef DEBUG
  printf("\n Cart3DAIM/aimInitialize   instance = %d!\n", inst);
#endif

  /* specify the number of analysis input and out "parameters" */
  *nIn     = NUMINPUT;
  *nOut    = NUMOUT;
  if (inst == -1) return CAPS_SUCCESS;

  /* get the storage */
  AIM_ALLOC(cartInstance, 1, aimStorage, aimInfo, status);
  *instStore = cartInstance;

  cartInstance->adaptFunctional = 0;

  initiate_cfdDesignStruct(&cartInstance->design);
  cartInstance->ntess = 0;
  cartInstance->tess  = NULL;

  initiate_mapAttrToIndexStruct(&cartInstance->groupMap);

  initiate_cfdBoundaryConditionStruct(&cartInstance->bcProps);

  cartInstance->aero_start = NULL;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 5;

  /* specify the name of each field variable */
  AIM_ALLOC(strs, *nFields, char *, aimInfo, status);
  strs[0]  = EG_strdup("Cp");
  strs[1]  = EG_strdup("Density");
  strs[2]  = EG_strdup("Velocity");
  strs[3]  = EG_strdup("Pressure");
  strs[4]  = EG_strdup("Displacement");
  for (i = 0; i < *nFields; i++)
    if (strs[i] == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
  *fnames  = strs;

  /* specify the dimension of each field variable */
  AIM_ALLOC(ints, *nFields, int, aimInfo, status);

  ints[0]  = 1;
  ints[1]  = 1;
  ints[2]  = 3;
  ints[3]  = 1;
  ints[4]  = 3;
  *franks   = ints;
  ints = NULL;

  /* specify if a field is an input field or output field */
  AIM_ALLOC(ints, *nFields, int, aimInfo, status);

  ints[0]  = FieldOut;
  ints[1]  = FieldOut;
  ints[2]  = FieldOut;
  ints[3]  = FieldOut;
  ints[4]  = FieldIn;
  *fInOut  = ints;
  ints = NULL;


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


/* ********************** Exposed AIM Functions ***************************** */
int
aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, int index,
          char **name, capsValue *defval)
{
  int status = CAPS_SUCCESS;

    /*! \page aimInputsCART3D AIM Inputs
     * The following list outlines the Cart3D inputs along with their default value available
     * through the AIM interface.
     */

#ifdef DEBUG
  printf(" Cart3DAIM/aimInputs   index = %d!\n", index);
#endif

  if (index == Tess_Params) {
    *name              = EG_strdup("Tess_Params");
    defval->type       = Double;
    defval->dim        = Vector;
    defval->nrow       = 3;
    defval->ncol       = 1;
    defval->units      = NULL;
    defval->vals.reals = (double *) EG_alloc(defval->nrow*sizeof(double));
    if (defval->vals.reals == NULL) {
      status = EGADS_MALLOC;
    } else {
      defval->vals.reals[0] = 0.025;
      defval->vals.reals[1] = 0.001;
      defval->vals.reals[2] = 15.00;
    }

      /*! \page aimInputsCART3D
     * - <B> Tess_Params = [double, double, double]</B>. <Default [0.025, 0.001, 15.00]> <br>
     * These parameters are used to create the surface mesh for Cart3D. Their order definition is as follows.
     *  1. Max Edge Length (0 is any length)
     *  2. Max Sag or distance from mesh segment and actual curved geometry
     *  3. Max angle in degrees between triangle facets
     */

  } else if (index == mesh2d) {
    *name                = EG_strdup("mesh2d");
    defval->type         = Boolean;
    defval->vals.integer = (int)false;
      /*! \page aimInputsCART3D
     * - <B> mesh2d = bool</B>. <Default false> <br>
     * Specify 2D analysis
     */
  } else if (index == outer_box) {
    *name             = EG_strdup("outer_box");
    defval->type      = Double;
    defval->vals.real = 30.;
      /*! \page aimInputsCART3D
     * - <B> outer_box = double</B>. <Default 30> <br>
     * Factor of outer boundary box based on geometry length scale defined by the  diagonal of the
     * 3D tightly fitting bounding box around body being modeled.
     */
  } else if (index == nDiv) {
    *name                = EG_strdup("nDiv");
    defval->type         = Integer;
    defval->vals.integer = 5;
      /*! \page aimInputsCART3D
     * - <B> nDiv = int</B>. <Default 5> <br>
     * nominal # of divisions in backgrnd mesh
     */

  } else if (index == maxR) {
    *name                = EG_strdup("maxR");
    defval->type         = Integer;
    defval->vals.integer = 7;
      /*! \page aimInputsCART3D
     * - <B> maxR = int</B>. <Default 7> <br>
     * Max Num of cell refinements to perform
     */
  } else if (index == preSpec) {
    *name                = EG_strdup("preSpec");
    defval->type         = Boolean;
    defval->vals.integer = (int)false;
      /*! \page aimInputsCART3D
     * - <B> preSpec = Bool</B>. <Default false> <br>
     * Use preSpec.c3d.cntl to guide mesh generation.
     * Enabled automatically if cart3d_preSpec bodies are provided.
     */
  } else if (index == symmX) {
    *name                = EG_strdup("symmX");
    defval->type         = Boolean;
    defval->vals.integer = (int)false;
      /*! \page aimInputsCART3D
     * - <B> symmX = bool</B>. <Default false> <br>
     * Symmetry on x-min boundary
     */
  } else if (index == symmY) {
    *name                = EG_strdup("symmY");
    defval->type         = Boolean;
    defval->vals.integer = (int)false;
      /*! \page aimInputsCART3D
     * - <B> symmY = int</B>. <Default false> <br>
     * Symmetry on y-min boundary
     */
  } else if (index == symmZ) {
    *name                = EG_strdup("symmZ");
    defval->type         = Boolean;
    defval->vals.integer = (int)false;
      /*! \page aimInputsCART3D
     * - <B> symmZ = int</B>. <Default false> <br>
     * Symmetry on z-min boundary
     */
  } else if (index == halfBody) {
    *name                = EG_strdup("halfBody");
    defval->type         = Boolean;
    defval->vals.integer = (int)false;
      /*! \page aimInputsCART3D
     * - <B> halfBody = bool</B>. <Default false> <br>
     * Input geometry is a half-body
     */
  } else if (index == Boundary_Condition) {
    *name              = EG_strdup("Boundary_Condition");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    //defval->units        = NULL;
    defval->dim          = Vector;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;

    /*! \page aimInputsCART3D
     * - <B>Boundary_Condition = NULL </B> <br>
     * See \ref cart3dBoundaryConditions for additional details.
     */

  } else if (index == Mach) {
    *name             = EG_strdup("Mach");
    defval->type      = Double;
    defval->vals.real = 0.0;
    defval->nullVal   = IsNull;
      /*! \page aimInputsCART3D
     * - <B> Mach = NULL</B>. <Default none> <br>
     */
  } else if (index == alpha) {
    *name             = EG_strdup("alpha");
    defval->type      = Double;
    defval->vals.real = 0.0;
      /*! \page aimInputsCART3D
     * - <B> alpha = double</B>. <Default 0.0> <br>
     * Angle of Attach in Degrees
     */
  } else if (index == beta) {
    *name             = EG_strdup("beta");
    defval->type      = Double;
    defval->vals.real = 0.0;
      /*! \page aimInputsCART3D
     * - <B> beta = double</B>. <Default 0.0> <br>
     * Side Slip Angle in Degrees
     */
  } else if (index == Gamma) {
    *name             = EG_strdup("gamma");
    defval->type      = Double;
    defval->vals.real = 1.4;
      /*! \page aimInputsCART3D
     * - <B> gamma = double</B>. <Default 1.4> <br>
     * Ratio of specific heats (default is air)
     */
  } else if (index == Pressure_Scale_Factor) {
    *name              = EG_strdup("Pressure_Scale_Factor");
    defval->type         = Double;
    defval->vals.real    = 1.0;

    /*! \page aimInputsCART3D
     * - <B>Pressure_Scale_Factor = 1.0</B> <br>
     * Value to scale Pressure data when transferring data. Data is scaled based on Pressure = Pressure_Scale_Factor*Pressure.
     */
  } else if (index == TargetCL) {
    *name             = EG_strdup("TargetCL");
    defval->type      = Double;
    defval->nullVal   = IsNull;
      /*! \page aimInputsCART3D
     * - <B> TargetCL = NULL</B>.<br>
     * Target lift coefficient. If set, Cart3D will adjust alpha to such that CL matches TargetCL.<br>
     * See TargetCL_ inputs for additional controls.
     */
  } else if (index == TargetCL_Tol) {
    *name             = EG_strdup("TargetCL_Tol");
    defval->type      = Double;
    defval->vals.real = 0.01;
      /*! \page aimInputsCART3D
     * - <B> TargetCL_Tol = double</B>. <Default 0.01> <br>
     * Target lift coefficient tolerance
     */
  } else if (index == TargetCL_Start_Iter) {
    *name                = EG_strdup("TargetCL_Start_Iter");
    defval->type         = Integer;
    defval->vals.integer = 60;
      /*! \page aimInputsCART3D
     * - <B> TargetCL_Start_Iter = int</B>. <Default 60> <br>
     * Iteration for first Target lift coefficient alpha adjustment
     */
  } else if (index == TargetCL_Freq) {
    *name                = EG_strdup("TargetCL_Freq");
    defval->type         = Integer;
    defval->vals.integer = 5;
      /*! \page aimInputsCART3D
     * - <B> TargetCL_Freq = int</B>. <Default 5> <br>
     * Iteration frequency for Target lift coefficient alpha adjustment
     */
  } else if (index == TargetCL_NominalAlphaStep) {
    *name             = EG_strdup("TargetCL_NominalAlphaStep");
    defval->type      = Double;
    defval->vals.real = 0.2;
      /*! \page aimInputsCART3D
     * - <B> TargetCL_NominalAlphaStep = double</B>. <Default 0.2> <br>
     * Initial alpha step and max step size in Degrees
     */
  } else if (index == maxCycles) {
    *name                = EG_strdup("maxCycles");
    defval->type         = Integer;
    defval->vals.integer = 1000;
      /*! \page aimInputsCART3D
     * - <B> maxCycles = int</B>. <Default 1000> <br>
     * Number of iterations
     */
  } else if (index == nMultiGridLevels) {
    *name                = EG_strdup("nMultiGridLevels");
    defval->type         = Integer;
    defval->vals.integer = 1;
      /*! \page aimInputsCART3D
     * - <B> nMultiGridLevels = int</B>. <Default 1> <br>
     * number of multigrid levels in the mesh (1 is a single mesh)
     */
  } else if (index == MultiGridCycleType) {
    *name                = EG_strdup("MultiGridCycleType");
    defval->type         = Integer;
    defval->vals.integer = 2;
    /*! \page aimInputsCART3D
     * - <B> MultiGridCycleType = int</B>. <Default 2> <br>
     * MultiGrid cycletype: 1 = "V-cycle", 2 = "W-cycle"
     * 'sawtooth' cycle is: V-cycle with MultiGridPreSmoothing = 1, MultiGridPostSmoothing = 0
     */
  }  else if (index == MultiGridPreSmoothing) {
    *name                = EG_strdup("MultiGridPreSmoothing");
    defval->type         = Integer;
    defval->vals.integer = 1;
    /*! \page aimInputsCART3D
     * - <B> MultiGridPreSmoothing = int</B>. <Default 1> <br>
     * number of pre-smoothing  passes in multigrid
     */
  } else if (index == MultiGridPostSmoothing) {
    *name                = EG_strdup("MultiGridPostSmoothing");
    defval->type         = Integer;
    defval->vals.integer = 1;
    /*! \page aimInputsCART3D
     * - <B> MultiGridPostSmoothing = int</B>. <Default 1> <br>
     * number of post-smoothing  passes in multigrid
     */
  } else if (index == RK) {
    *name             = EG_strdup("RK");
    defval->type      = Double;
    defval->dim       = Array2D;
    defval->nrow      = 5;
    defval->ncol      = 2;
    defval->sfixed    = Fixed;
    defval->lfixed    = Change;
    AIM_ALLOC(defval->vals.reals, defval->nrow*defval->ncol, double, aimInfo, status);
    defval->vals.reals[2*0+0] = 0.0695;
    defval->vals.reals[2*0+1] = 1;
    defval->vals.reals[2*1+0] = 0.1602;
    defval->vals.reals[2*1+1] = 0;
    defval->vals.reals[2*2+0] = 0.2898;
    defval->vals.reals[2*2+1] = 0;
    defval->vals.reals[2*3+0] = 0.5060;
    defval->vals.reals[2*3+1] = 0;
    defval->vals.reals[2*4+0] = 1.0;
    defval->vals.reals[2*4+1] = 0;
    /*! \page aimInputsCART3D
     * - <B> RK = double matrix, </B>. <Default [[0.0695, 1],[0.1602,0], [0.2898,0], [0.5060,0], [1.0,0]> <br>
     * Runge-Kutta Stage Coefficients. May be any number of rows, must be 2 columns with the 2nd column
     * indicating if gradients are evaluated (0 = no new evaluation at this stage, 1 = Yes, re-evaluate at this stage)
     */
  } else if (index == CFL) {
    *name             = EG_strdup("CFL");
    defval->type      = Double;
    defval->vals.real = 1.2;
    /*! \page aimInputsCART3D
     * - <B> CFL = double</B>. <Default 1.2> <br>
     * CFL number typically between 0.9 and 1.4
     */
  } else if (index == Limiter) {
    *name                = EG_strdup("Limiter");
    defval->type         = Integer;
    defval->vals.integer = 2;
    /*! \page aimInputsCART3D
     * - <B> Limiter = int</B>. <Default 2> <br>
     * organized in order of increasing dissipation. <br>
     * 0 = no Limiter, 1 = Barth-Jespersen, 2 = van Leer, 3 = sin limiter, 4 = van Albada, 5 = MinMod
     */
  } else if (index == FluxFun) {
    *name                = EG_strdup("FluxFun");
    defval->type         = Integer;
    defval->vals.integer = 0;
    /*! \page aimInputsCART3D
     * - <B> FluxFun = int</B>. <Default 0> <br>
     * 0 = van Leer, 1 = van Leer-Hanel, 2 = Colella 1998, 3 = HLLC (alpha test)
     */
  } else if (index == iForce) {
    *name                = EG_strdup("iForce");
    defval->type         = Integer;
    defval->vals.integer = 1;
    /*! \page aimInputsCART3D
     * - <B> iForce = int</B>. <Default 10> <br>
     * Report force & mom. information every iForce cycles
     */
  } else if (index == iHist) {
    *name                = EG_strdup("iHist");
    defval->type         = Integer;
    defval->vals.integer = 1;
    /*! \page aimInputsCART3D
     * - <B> iHist = int</B>. <Default 1> <br>
     * Update 'history.dat' every iHist cycles
     */
  } else if (index == nOrders) {
    *name                = EG_strdup("nOrders");
    defval->type         = Integer;
    defval->vals.integer = 8;
    /*! \page aimInputsCART3D
     * - <B> nOrders = int</B>. <Default 8> <br>
     * Num of orders of Magnitude reduction in residual
     */

  } else if (index == nAdaptCycles) {
    *name                = EG_strdup("nAdaptCycles");
    defval->type         = Integer;
    defval->nullVal      = NotAllowed;
    defval->vals.integer = 0;

    /*! \page aimInputsCART3D
     * - <B> nAdaptCycles = 0</B> <br>
     * Number of adaptation cycles.
     */

  } else if (index == Adapt_Functional) {
    *name                = EG_strdup("Adapt_Functional");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

   /*! \page aimInputsCART3D
    * - <B> Adapt_Functional = NULL</B> <br>
    * Single valued tuple that defines the functional used to drive mesh adaptation, see \ref cfdDesignFunctional for additional details on functionals.
    */

  } else if (index == Design_Variable) {
    *name                 = EG_strdup("Design_Variable");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsCART3D
     * - <B> Design_Variable = NULL</B> <br>
     * The design variable tuple is used to input design variable information for optimization, see \ref cfdDesignVariable for additional details.
     */
  } else if (index == Design_Functional) {
    *name                = EG_strdup("Design_Functional");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

   /*! \page aimInputsCART3D
    * - <B> Design_Functional = NULL</B> <br>
    * The design objective tuple is used to input objective information for optimization, see \ref cfdDesignFunctional for additional details.
    * The value of the design functionals become available as Dynamic Output Value Objects using the "name" of the functionals.
    */

  } else if (index == Design_Sensitivity) {
    *name                = EG_strdup("Design_Sensitivity");
    defval->type         = Boolean;
    defval->lfixed       = Fixed;
    defval->vals.integer = (int)false;
    defval->dim          = Scalar;

    /*! \page aimInputsCART3D
     * - <B> Design_Sensitivity = False</B> <br>
     * Create geometric sensitivities Cart3D input files needed to compute Design_Functional sensitivities w.r.t Design_Variable.
     */

  } else if (index == Design_Adapt) {
    *name                = EG_strdup("Design_Adapt");
    defval->type         = String;
    defval->nullVal      = IsNull;
    defval->lfixed       = Fixed;
    defval->vals.string  = NULL;
    defval->dim          = Scalar;

   /*! \page aimInputsCART3D
    * - <B> Design_Adapt = NULL</B> <br>
    * String name of a Design_Functional to be used for adjoint based mesh adaptation.
    */

  } else if (index == Design_Run_Config) {
    *name                = EG_strdup("Design_Run_Config");
    defval->type         = String;
    defval->nullVal      = NotAllowed;
    defval->dim          = Scalar;
    AIM_STRDUP(defval->vals.string, "production", aimInfo, status);

   /*! \page aimInputsCART3D
    * - <B> Design_Run_Config = "production"</B> <br>
    * run_config = debug || archive || standard || production <br>
    * debug:      no files are deleted; easily traceable but needs large disk-space<br>
    * archive:    compresses and tars critical files and deletes temp files;<br>
    *             becomes slower for cases with large number of design variables<br>
    *             but fully traceable<br>
    * standard:   similar to production, but keeps more files, in particular adjoint<br>
    *             solution(s) on the finest mesh<br>
    * production: keeps critical files; reasonable storage and max speed
    */

  } else if (index == Design_Gradient_Memory_Budget) {
    *name                = EG_strdup("Design_Gradient_Memory_Budget");
    defval->type         = Integer;
    defval->nullVal      = NotAllowed;
    defval->dim          = Scalar;
    defval->vals.integer = 32;

   /*! \page aimInputsCART3D
    * - <B> Design_Gradient_Memory_Budget = 32</B> <br>
    * This flag controls the parallel evaluation of components of the gradient. For
    * example, gradient of objective function J with respect to design variables x,
    * y, z is [ dJ/dx dJ/dy dJ/dz ]^T. The gradient components are independent and
    * hence, can be evaluated in parallel. Evaluating all the components
    * simultaneously is ideal, but you are limited by the number of design
    * variables and the size of your mesh (memory limit).  The framework can
    * compute efficient job partitioning automatically. To do this, you just need
    * to specify the memory limit for the run (in GB) via the flag
    */

  } else if (index == Xslices) {
    *name               = EG_strdup("Xslices");
    defval->type          = Double;
    defval->lfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;
    /*! \page aimInputsCART3D
     * - <B> Xslices = double </B> or <B> [double, ... , double] </B> <br>
     * X slice locations created in output.
     */
  } else if (index == Yslices) {
    *name               = EG_strdup("Yslices");
    defval->type          = Double;
    defval->lfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;

     /*! \page aimInputsCART3D
    * - <B> Yslices = double </B> or <B> [double, ... , double] </B> <br>
    * Y slice locations created in output.
      */
  } else if (index == Zslices) {
    *name               = EG_strdup("Zslices");
    defval->type          = Double;
    defval->lfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;

     /*! \page aimInputsCART3D
    * - <B> Zslices = double </B> or <B> [double, ... , double] </B> <br>
    * Z slice locations created in output.
    */

  } else if (index == y_is_spanwise) {
    *name = EG_strdup("y_is_spanwise");
    defval->type = Boolean;
    defval->vals.integer = (int)false;

    /*! \page aimInputsCART3D
     * - <B> y_is_spanwise = bool </B> <Default false> <br>
     * If false, then alpha is defined in the x-y plane, <br>
     * otherwise alpha is in the x-z plane
     */
  } else if (index == Model_X_axis) {
    *name = EG_strdup("Model_X_axis");
    defval->type = String;
    defval->vals.string = EG_strdup("-Xb");

    /*! \page aimInputsCART3D
     * - <B> Model_X_axis = "-Xb" </B> <br>
     * Model_X_axis defines x-axis orientation.
     */
  } else if (index == Model_Y_axis) {
    *name = EG_strdup("Model_Y_axis");
    defval->type = String;
    defval->vals.string = EG_strdup("Yb");

    /*! \page aimInputsCART3D
     * - <B> Model_Y_axis = "Yb" </B> <br>
     * Model_Y_axis defines y-axis orientation.
     */
  } else if (index == Model_Z_axis) {
    *name = EG_strdup("Model_Z_axis");
    defval->type = String;
    defval->vals.string = EG_strdup("-Zb");

    /*! \page aimInputsCART3D
     * - <B> Model_Z_axis = "-Zb" </B> <br>
     * Model_Z_axis defines z-axis orientation.
     */
  } else if (index == Restart) {
    *name = EG_strdup("Restart");
    defval->type = Boolean;
    defval->vals.integer = (int)false;

    /*! \page aimInputsCART3D
     * - <B> Restart = False </B> <br>
     * Use the "restart" option for aero.csh or not
     */
  } else if (index == aerocsh) {
    *name                = EG_strdup("aerocsh");
    defval->type         = String;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->nrow         = 0;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

   /*! \page aimInputsCART3D
    * - <B> aerocsh = NULL</B> <br>
    * List of strings that can be used to override defaults in the Cart3D aero.csh script.
    * Please refer to the Cart3D documentation for all available aero.csh inputs.
    */

  } else if (index == aerocsh_arg) {
    *name                = EG_strdup("aerocsh_arg");
    defval->type         = String;
    AIM_STRDUP(defval->vals.string, "", aimInfo, status);

   /*! \page aimInputsCART3D
    * - <B> aerocsh_arg = NULL</B> <br>
    * String of arguments used when executing aero.csh
    */

  } else {
    AIM_ERROR(aimInfo, "Uknown index %d", index);
    status = CAPS_BADINDEX;
    goto cleanup;
  }

  AIM_NOTNULL(*name, aimInfo, status);
cleanup:
  if (status != CAPS_SUCCESS) AIM_FREE(*name);
  return status;
}

// ********************** AIM Function Break *****************************
int
aimUpdateState(void *instStore, void *aimInfo,
               capsValue *aimInputs)
{
  int status; // Function return status

  int i, j;

  double params[3], box[6], size;

  // Body parameters
  const char *intents;
  int ntess = 0;
  int nBody = 0; // Number of bodies
  ego *bodies = NULL; // EGADS body objects

  int atype, alen;
  const int    *ints;
  const char   *string;
  const double *reals;
  double xyz[3];

  int oclass, mtype, nnode, *senses;
  double trange[4];
  int nedge;
  ego *eedges=NULL, *enodes=NULL, geom;

  aimStorage   *cartInstance;

  cartInstance = (aimStorage *) instStore;
  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
  AIM_STATUS(aimInfo, status);

  if (nBody < 1 || bodies == NULL) {
    AIM_ERROR(aimInfo, "Cart3D AIM requires at least 1 body");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if (aimInputs[RK-1].ncol != 2) {
    AIM_ERROR(aimInfo, "RK input must have 2 columns! Input has %d columns.", aimInputs[RK-1].ncol);
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, Boundary_Condition) == CAPS_SUCCESS) {

    // Get capsGroup name and index mapping to make sure all bodies have a capsGroup value
    status = create_CAPSGroupAttrToIndexMap(nBody,
                                            bodies,
                                            1, // Down to face
                                            &cartInstance->groupMap);
    AIM_STATUS(aimInfo, status);

    status = destroy_cfdBoundaryConditionStruct(&cartInstance->bcProps);
    AIM_STATUS(aimInfo, status);

    if (cartInstance->groupMap.numAttribute > 0) {
      AIM_ALLOC(cartInstance->bcProps.surfaceProp, cartInstance->groupMap.numAttribute, cfdSurfaceStruct, aimInfo, status);
      cartInstance->bcProps.numSurfaceProp = cartInstance->groupMap.numAttribute;
      for (i = 0; i < cartInstance->bcProps.numSurfaceProp; i++) {
        status = initiate_cfdSurfaceStruct(&cartInstance->bcProps.surfaceProp[i]);
        AIM_STATUS(aimInfo, status);
      }

      for (i = 0; i < cartInstance->bcProps.numSurfaceProp; i++) {
        cartInstance->bcProps.surfaceProp[i].bcID = cartInstance->groupMap.attributeIndex[i];
        AIM_STRDUP(cartInstance->bcProps.surfaceProp[i].name, cartInstance->groupMap.attributeName[i], aimInfo, status);
      }
    }

    status = cart3d_getBoundaryCondition(aimInfo,
                                         aimInputs[Boundary_Condition-1].length,
                                         aimInputs[Boundary_Condition-1].vals.tuple,
                                         &cartInstance->groupMap,
                                         &cartInstance->bcProps);
    AIM_STATUS(aimInfo, status);
  }

  // Get design variables
  if (aimInputs[Design_Variable-1].nullVal == NotNull &&
      (aim_newAnalysisIn(aimInfo, Design_Variable) == CAPS_SUCCESS ||
       cartInstance->design.numDesignVariable == 0) ) {

    if (aimInputs[Design_Functional-1].nullVal == IsNull) {
      AIM_ERROR(aimInfo, "\"Design_Variable\" has been set, but no values have been provided for \"Design_Functional\"!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    status = cfd_getDesignVariable(aimInfo,
                                   aimInputs[Design_Variable-1].length,
                                   aimInputs[Design_Variable-1].vals.tuple,
                                   &cartInstance->design.numDesignVariable,
                                   &cartInstance->design.designVariable);
    AIM_STATUS(aimInfo, status);
  }


  // Check for Design_Functional
  if (aimInputs[Design_Functional-1].nullVal == NotNull &&
      (cartInstance->design.numDesignFunctional == 0 ||
       aim_newAnalysisIn(aimInfo, Design_Functional ) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, Design_Variable   ) == CAPS_SUCCESS) ) {

      status = cfd_getDesignFunctional(aimInfo,
                                       aimInputs[Design_Functional-1].length,
                                       aimInputs[Design_Functional-1].vals.tuple,
                                       &cartInstance->bcProps,
                                       cartInstance->design.numDesignVariable,
                                       cartInstance->design.designVariable,
                                       &cartInstance->design.numDesignFunctional,
                                       &cartInstance->design.designFunctional);
      AIM_STATUS(aimInfo, status);
  }

  // Check if adaptation is requested
  if (aimInputs[Design_Functional-1].nullVal == IsNull &&
      aimInputs[Adapt_Functional-1].nullVal == NotNull &&
      aim_newAnalysisIn(aimInfo, Adapt_Functional) == CAPS_SUCCESS) {

    i = cartInstance->adaptFunctional == NULL ? 0 : 1;

    status = cfd_getDesignFunctional(aimInfo,
                                     aimInputs[Adapt_Functional-1].length,
                                     aimInputs[Adapt_Functional-1].vals.tuple,
                                     &cartInstance->bcProps,
                                     0,
                                     NULL,
                                     &i,
                                     &cartInstance->adaptFunctional);
    AIM_STATUS(aimInfo, status);
  }


  // Only do surface meshing if the geometry or Tess_Params have changed
  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, Tess_Params) == CAPS_SUCCESS) {

    // skip any preSpec and Sensor bodies
    for (i = 0; i < nBody; i++) {
      status = EG_attributeRet(bodies[i], CART3D_SENSOR, &atype, &alen,
                               &ints, &reals, &string);
      if (status == EGADS_SUCCESS) {
        if ( atype != ATTRSTRING ) {
            AIM_ERROR(aimInfo, "%s should be a string value!\n", CART3D_SENSOR);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
        if (aim_isNodeBody(bodies[i], xyz) == CAPS_SUCCESS) continue;

        status = EG_getBodyTopos(bodies[i], NULL, EDGE, &nedge, &eedges);
        AIM_STATUS(aimInfo, status);
        if (nedge != 1) {
          AIM_ERROR(aimInfo, "%s %s may only have 1 edge. Found %d!\n", CART3D_SENSOR, string, nedge);
          status = CAPS_BADVALUE;
          goto cleanup;
        }
        AIM_NOTNULL(eedges, aimInfo, status);

        status = EG_getTopology(eedges[0], &geom, &oclass, &mtype, trange, &nnode, &enodes, &senses);
        AIM_STATUS(aimInfo, status);
        if (geom == NULL || geom->mtype != LINE) {
          AIM_ERROR(aimInfo, "%s %s must be a LINE!\n", CART3D_SENSOR, string);
          status = CAPS_BADVALUE;
          goto cleanup;
        }

        status = EG_attributeRet(bodies[i], CART3D_EQAREA, &atype, &alen,
                                 &ints, &reals, &string);
        if (status == EGADS_SUCCESS) {
          if ( !(atype == ATTRREAL && alen == 1) ) {
            AIM_ERROR(aimInfo, "%s should be followed by a single real value!\n", CART3D_EQAREA);
            status = CAPS_BADVALUE;
            goto cleanup;
          }
        }
        AIM_FREE(eedges);
        continue;
      }

      status = EG_attributeRet(bodies[i], CART3D_BBOX, &atype, &alen,
                               &ints, &reals, &string);
      if (status == EGADS_SUCCESS) {
        if ( !((atype == ATTRINT || atype == ATTRREAL) && alen == 1) ) {
            AIM_ERROR(aimInfo, "%s should be followed by a single int value!\n", CART3D_BBOX);
            status = CAPS_BADVALUE;
            goto cleanup;
        }
        continue;
      }

      ntess++;
    }

    AIM_REALL(cartInstance->tess, ntess, ego, aimInfo, status);
    cartInstance->ntess = ntess;

    for (i = j = 0; i < nBody; i++) {
      status = EG_attributeRet(bodies[i], CART3D_SENSOR, &atype, &alen,
                               &ints, &reals, &string);
      if (status == EGADS_SUCCESS) continue;

      status = EG_attributeRet(bodies[i], CART3D_BBOX, &atype, &alen,
                               &ints, &reals, &string);
      if (status == EGADS_SUCCESS) continue;

      status = EG_getBoundingBox(bodies[i], box);
      AIM_STATUS(aimInfo, status);

      size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                  (box[4]-box[1])*(box[4]-box[1]) +
                  (box[5]-box[2])*(box[5]-box[2]));
      printf(" Body size = %lf\n", size);

      params[0] = aimInputs[Tess_Params-1].vals.reals[0]*size;
      params[1] = aimInputs[Tess_Params-1].vals.reals[1]*size;
      params[2] = aimInputs[Tess_Params-1].vals.reals[2];
      printf(" Tessellating body %d with  MaxEdge = %lf  Sag = %lf  Angle = %lf\n",
             i+1, params[0], params[1], params[2]);

      status = EG_makeTessBody(bodies[i], params, &cartInstance->tess[j]);
      AIM_STATUS(aimInfo, status);

      /* store away the tessellation */
      status = aim_newTess(aimInfo, cartInstance->tess[j]);
      AIM_STATUS(aimInfo, status);

      j++;
    }
  }

  status = CAPS_SUCCESS;
cleanup:
  AIM_FREE(eedges);

  return status;
}

// ********************** AIM Function Break *****************************
int
aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{
  int          i, j, varid, status, nBody;
  int          atype, alen, cleanDir, level;
  int          oclass, mtype, nnode, nedge;
  int          foundSref=(int)false, foundCref=(int)false, foundXref=(int)false;
  int          makePreSpec=(int)false, *senses;
  double       Sref, Cref, Xref, Yref, Zref, box[6], xyz[3], xyz0[9], xyz1[9], trange[4];
  char         line[128], aimFile[PATH_MAX], designxml[42];
  unsigned     options = 0;
  const int    *ints;
  const char   *string, *eqstring, *intents;
  const double *reals;
  const char   *inputsDir="inputs";
  char         args[128];
  p_tsXddm     p_xddm = NULL;
  ego          *bodies, *eedges=NULL, *enodes=NULL, geom;
  FILE         *fp=NULL;
  const aimStorage *cartInstance;

#ifdef DEBUG
  printf(" Cart3DAIM/aimPreAnalysis!\n");
#endif

  if (aimInputs == NULL) {
#ifdef DEBUG
    printf(" Cart3DAIM/aimPreAnalysis -- NULL inputs!\n");
#endif
    return CAPS_NULLOBJ;
  }

  if (aimInputs[Mach-1].nullVal == IsNull) {
    AIM_ERROR(aimInfo, "'Mach' is not set");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  cartInstance = (const aimStorage *) instStore;

  snprintf(designxml, 42, "%s/design.xml", runDir);

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
  AIM_STATUS(aimInfo, status);

  if (nBody < 1 || bodies == NULL) {
    AIM_ERROR(aimInfo, "Cart3D AIM requires at least 1 body");
    status = CAPS_BADVALUE;
    goto cleanup;
  }

  // Initialize reference values
  Sref = Cref = 1.0;
  Xref = Yref = Zref = 0.0;

  // Loop over bodies and look for reference quantity attributes
  for (i = 0; i < nBody; i++) {
    status = EG_attributeRet(bodies[i], "capsReferenceArea", &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
        if (atype == ATTRREAL && alen == 1) {
            Sref = (double) reals[0];
            foundSref = (int)true;
        } else {
            AIM_ERROR(aimInfo, "capsReferenceArea should be followed by a single real value!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceChord", &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
        if (atype == ATTRREAL && alen == 1) {
            Cref = (double) reals[0];
            foundCref = (int)true;
        } else {
            AIM_ERROR(aimInfo, "capsReferenceChord should be followed by a single real value!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceX", &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
        if (atype == ATTRREAL && alen == 1) {
            Xref = (double) reals[0];
            foundXref = (int)true;
        } else {
            AIM_ERROR(aimInfo, "capsReferenceX should be followed by a single real value!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceY", &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
        if (atype == ATTRREAL && alen == 1) {
            Yref = (double) reals[0];
        } else {
            AIM_ERROR(aimInfo, "capsReferenceY should be followed by a single real value!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    status = EG_attributeRet(bodies[i], "capsReferenceZ", &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_SUCCESS){
        if (atype == ATTRREAL && alen == 1) {
            Zref = (double) reals[0];
        } else {
            AIM_ERROR(aimInfo, "capsReferenceZ should be followed by a single real value!\n");
            status = CAPS_BADVALUE;
            goto cleanup;
        }
    }

    status = EG_attributeRet(bodies[i], CART3D_BBOX, &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_SUCCESS) {
      makePreSpec = (int)true;
      aimInputs[preSpec-1].vals.integer = (int)true;
    }
  }

  if (foundSref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceArea is not set on any body!");
      status = CAPS_BADVALUE;
      goto cleanup;
  }
  if (foundCref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceChord is not set on any body!");
      status = CAPS_BADVALUE;
      goto cleanup;
  }
  if (foundXref == (int)false) {
      AIM_ERROR(aimInfo, "capsReferenceX is not set on any body!");
      status = CAPS_BADVALUE;
      goto cleanup;
  }


  // Create the inputs directory
  status = aim_mkDir(aimInfo, inputsDir);
  AIM_STATUS(aimInfo, status);

  // remove any previous adapt directories
  status = aim_rmDir(aimInfo, "inputs/adapt??");
  AIM_STATUS(aimInfo, status);

  if ( aimInputs[Design_Functional-1].nullVal == NotNull ) {
    // Create the runDir directory
    status = aim_mkDir(aimInfo, runDir);
    AIM_STATUS(aimInfo, status);

    // remove the previous tbl file
    snprintf(line, 128, "%s/entire.DP1.tbl", runDir);
    status = aim_rmFile(aimInfo, line);
    AIM_STATUS(aimInfo, status);

    if (aim_newAnalysisIn(aimInfo, Design_Run_Config) == CAPS_SUCCESS ||
        aim_newAnalysisIn(aimInfo, Design_Gradient_Memory_Budget) == CAPS_SUCCESS) {

      snprintf(line, 128, "%s/OPTIONS", runDir);
      fp = aim_fopen(aimInfo, line, "w");
      if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open: %s", line);
        status = CAPS_IOERR;
        goto cleanup;
      }

      fprintf(fp, "run_config = %s\n", aimInputs[Design_Run_Config-1].vals.string);
      fprintf(fp, "gradient_memory_budget = %d\n", aimInputs[Design_Gradient_Memory_Budget-1].vals.integer);
      fclose(fp); fp = NULL;
    }
  }

  // Write out inputs/Config.xml
  status = cart3d_configxml(aimInfo, &cartInstance->groupMap);
  AIM_STATUS(aimInfo, status);

  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS) {
    // Remove old geometric sensitivity files
    snprintf(line, 128, "%s/geometry_*", runDir);
    status = aim_rmDir(aimInfo, line);
    AIM_STATUS(aimInfo, status);
  }

  // Copy over aero.csh with modified inputs
  if (aim_isFile(aimInfo, "inputs/aero.csh") == CAPS_NOTFOUND ||
      aim_newAnalysisIn(aimInfo, Limiter      ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, mesh2d       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, nAdaptCycles ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, maxR         ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, preSpec      ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, y_is_spanwise) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, aerocsh      ) == CAPS_SUCCESS) {
    status = cart3d_aerocsh(aimInfo,
                            aimInputs[Limiter-1].vals.integer,
                            aimInputs[mesh2d-1].vals.integer,
                            aimInputs[nAdaptCycles-1].vals.integer,
                            aimInputs[maxR-1].vals.integer,
                            aimInputs[preSpec-1].vals.integer,
                            aimInputs[y_is_spanwise-1].vals.integer,
                            aimInputs[aerocsh-1].length,
                            aimInputs[aerocsh-1].vals.string);
    AIM_STATUS(aimInfo, status);
  }

  /* check to see if previous M*A*B*_DP1 need to be removed
   * it should be kept only to restart a gradient calculation
   * */
  if (aimInputs[Design_Functional-1].nullVal == NotNull &&
      aim_isFile(aimInfo, designxml) == CAPS_SUCCESS) {

    /* check if only Design_Sensitivity is new input */
    cleanDir = aim_newGeometry(aimInfo) == CAPS_SUCCESS;
    for (i = 0; i < NUMINPUT && cleanDir == (int)false; i++) {
      if (i+1 == Design_Sensitivity) continue;
      if (aim_newAnalysisIn(aimInfo, i+1) == CAPS_SUCCESS)
        cleanDir = (int)true;
    }

    /* the desing.xml already has derivatives, then the directory needs to be removed */
    if (cleanDir == (int)false) {
      // Read design.xml
      status = aim_file(aimInfo, designxml, aimFile);
      AIM_STATUS(aimInfo, status);

      p_xddm = xddm_readFile(aimFile, "/Optimize", &options);
      if (p_xddm == NULL) {
        AIM_ERROR(aimInfo, "xddm_readFile failed to read: %s", aimFile);
        status = CAPS_IOERR;
        goto cleanup;
      }

      if (p_xddm->nk != 1) {
        AIM_ERROR(aimInfo, "'%s' should only have 1 design point, found %d", aimFile, p_xddm->nk);
        status = CAPS_IOERR;
        goto cleanup;
      }

      /* Loop over analysis params and look for derivatives */
      for (i = 0; i < p_xddm->a_kids->na; i++) {
        if (p_xddm->a_kids->a_ap[i].ndvs > 0) {
          cleanDir = (int)true;
          break;
        }
      }
    }

    if (cleanDir == (int)true) {
      // remove the previous design run directoy
      snprintf(line, 128, "%s/M*A*B*_DP1", runDir);
      //printf("****** Removing %s\n", line);
      status = aim_rmDir(aimInfo, line);
      AIM_STATUS(aimInfo, status);
    }
  }

  // Get design functionals
  // alpha, beta, Mach might be design variables that must be updated
  if ( aimInputs[Design_Functional-1].nullVal == NotNull &&
      ((int)true == (int)true || // right now design.xml must always be overwritten
       aim_newAnalysisIn(aimInfo, Design_Functional ) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, Design_Variable   ) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, Design_Sensitivity) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, Design_Adapt      ) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, alpha             ) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, beta              ) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, Mach              ) == CAPS_SUCCESS ||
       aim_isFile(aimInfo, "design/design.xml"      ) == CAPS_NOTFOUND)) {

      // write design.xml
      status = cart3d_designxml(aimInfo, &cartInstance->design,
                                aimInputs[Design_Sensitivity-1].vals.integer,
                                aimInputs[Design_Adapt-1].vals.string);
      AIM_STATUS(aimInfo, status);
  }

  // Check if adaptation is requested
  if (aimInputs[Design_Functional-1].nullVal == IsNull &&
      aimInputs[Adapt_Functional-1 ].nullVal == NotNull &&
      (aim_newAnalysisIn(aimInfo, Adapt_Functional ) == CAPS_SUCCESS ||
       aim_isFile(aimInfo, "inputs/Functionals.xml") == CAPS_NOTFOUND)) {

    AIM_NOTNULL(cartInstance->adaptFunctional, aimInfo, status);
    status = cart3d_functionalsxml(aimInfo, cartInstance->adaptFunctional);
    AIM_STATUS(aimInfo, status);
  }


  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, Tess_Params       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, Design_Variable   ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, Design_Sensitivity) == CAPS_SUCCESS ||
      aim_isFile(aimInfo, "design/builder.xml"     ) == CAPS_NOTFOUND) {

    // Only compute sensitivities if needed
    if (aimInputs[Design_Sensitivity-1].vals.integer == (int)true) {

      status = cart3d_Components(aimInfo,
                                 cartInstance->design.numDesignVariable,
                                 cartInstance->design.designVariable,
                                 &cartInstance->groupMap,
                                 cartInstance->ntess,
                                 cartInstance->tess);
      AIM_STATUS(aimInfo, status);

    } else {

      status = cart3d_Components(aimInfo,
                                 0,
                                 NULL,
                                 &cartInstance->groupMap,
                                 cartInstance->ntess,
                                 cartInstance->tess);
      AIM_STATUS(aimInfo, status);
    }
  }

  // Only do volume meshing if the geometry has been changed or a meshing input variable has been changed
  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, Tess_Params) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, mesh2d     ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, outer_box  ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, nDiv       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, maxR       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, symmX       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, symmY       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, symmZ       ) == CAPS_SUCCESS ||
      aim_newAnalysisIn(aimInfo, halfBody   ) == CAPS_SUCCESS ||
      aim_isFile(aimInfo, "inputs/input.c3d") == CAPS_NOTFOUND ||
      aim_isFile(aimInfo, "inputs/preSpec.c3d.cntl") == CAPS_NOTFOUND) {

    args[0] = '\0';
    if (aimInputs[mesh2d-1].vals.integer == (int)true) {
      strcat(args, " -mesh2d");
    }
    if (aimInputs[symmX-1].vals.integer == (int)true) {
      strcat(args, " -symmX");
    }
    if (aimInputs[symmY-1].vals.integer == (int)true) {
      strcat(args, " -symmY");
    }
    if (aimInputs[symmZ-1].vals.integer == (int)true) {
      strcat(args, " -symmZ");
    }
    if (aimInputs[halfBody-1].vals.integer == (int)true) {
      strcat(args, " -halfBody");
    }

    snprintf(line, 128, "autoInputs -r %lf -nDiv %d -maxR %d %s > autoInputs.out",
             aimInputs[outer_box-1].vals.real, aimInputs[nDiv-1].vals.integer,
             aimInputs[maxR-1].vals.integer, args);
    printf(" Executing: %s\n", line);
    status = aim_system(aimInfo, inputsDir, line);
    AIM_STATUS(aimInfo, status);
  }

  /* create input.cntl */
  fp = aim_fopen(aimInfo, "inputs/input.cntl", "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Cannot open: inputs/input.cntl!");
    status = CAPS_IOERR;
    goto cleanup;
  }

  fprintf(fp, "$__Case_Information:\n\n");
  fprintf(fp, "Mach     %lf\n", aimInputs[Mach-1].vals.real);
  fprintf(fp, "alpha    %lf\n", aimInputs[alpha-1].vals.real);
  fprintf(fp, "beta     %lf\n", aimInputs[beta-1].vals.real);
  fprintf(fp, "gamma    %lf\n", aimInputs[Gamma-1].vals.real);

  fprintf(fp, "\n$__File_Name_Information:\n\n");
  fprintf(fp, "MeshInfo           Mesh.c3d.Info\n");
  fprintf(fp, "MeshFile           Mesh.mg.c3d\n\n");

  fprintf(fp, "$__Solver_Control_Information:\n\n");
  for (i = 0; i < aimInputs[RK-1].nrow; i++)
    fprintf(fp, "RK        %lf     %d\n", aimInputs[RK-1].vals.reals[2*i+0], (int)aimInputs[RK-1].vals.reals[2*i+1]);
  fprintf(fp, "\n");
  fprintf(fp, "CFL           %lf\n", aimInputs[CFL-1].vals.real);
  fprintf(fp, "Limiter       %d\n",  aimInputs[Limiter-1].vals.integer);
  fprintf(fp, "FluxFun       %d\n",  aimInputs[FluxFun-1].vals.integer);
  fprintf(fp, "maxCycles     %d\n",  aimInputs[maxCycles-1].vals.integer);
  fprintf(fp, "Precon        0\n");
  fprintf(fp, "wallBCtype    0\n");
  fprintf(fp, "nMGlev        %d\n",   aimInputs[nMultiGridLevels-1].vals.integer);
  fprintf(fp, "MG_cycleType  %d\n",   aimInputs[MultiGridCycleType-1].vals.integer);
  fprintf(fp, "MG_nPre       %d\n",   aimInputs[MultiGridPreSmoothing-1].vals.integer);
  fprintf(fp, "MG_nPost      %d\n\n", aimInputs[MultiGridPostSmoothing-1].vals.integer);

  fprintf(fp, "$__Boundary_Conditions:\n\n");
/*                   # BC types: 0 = FAR FIELD
                                 1 = SYMMETRY
                                 2 = INFLOW  (specify all)
                                 3 = OUTFLOW (simple extrap)  */
                       /* (0/1/2) direction - Low BC - Hi BC */
  if (aimInputs[symmX-1].vals.integer == (int)true) {
    fprintf(fp, "Dir_Lo_Hi     0   1 0\n");
  } else {
    fprintf(fp, "Dir_Lo_Hi     0   0 0\n");
  }
  if (aimInputs[symmY-1].vals.integer == (int)true) {
    fprintf(fp, "Dir_Lo_Hi     1   1 0\n");
  } else {
    fprintf(fp, "Dir_Lo_Hi     1   0 0\n");
  }
  if (aimInputs[symmZ-1].vals.integer == (int)true) {
    fprintf(fp, "Dir_Lo_Hi     2   1 0\n\n");
  } else if (aimInputs[mesh2d-1].vals.integer == (int)true) {
    fprintf(fp, "Dir_Lo_Hi     2   1 1\n\n");
  } else {
    fprintf(fp, "Dir_Lo_Hi     2   0 0\n\n");
  }

  for (i = 0; i < cartInstance->bcProps.numSurfaceProp; i++) {
    if (cartInstance->bcProps.surfaceProp[i].surfaceType == Farfield) {

      // Generic INLET/EXHAUST BC -- 2004
      fprintf(fp, "SurfBC %d  %12.16e  %12.16e %12.16e %12.16e  %12.16e  # compID Rho Xvel Yvel Zvel Press\n",
              cartInstance->bcProps.surfaceProp[i].bcID,
              cartInstance->bcProps.surfaceProp[i].staticDensity,
              cartInstance->bcProps.surfaceProp[i].uVelocity,
              cartInstance->bcProps.surfaceProp[i].vVelocity,
              cartInstance->bcProps.surfaceProp[i].wVelocity,
              cartInstance->bcProps.surfaceProp[i].staticPressure);

    } else if (cartInstance->bcProps.surfaceProp[i].surfaceType == BackPressure) {

      //InletPressRatio -- 2018
      fprintf(fp, "InletPressRatioBC %d  %12.16e  #  compID  P/P_inf\n",
              cartInstance->bcProps.surfaceProp[i].bcID,
              cartInstance->bcProps.surfaceProp[i].staticPressure);

    } else if (cartInstance->bcProps.surfaceProp[i].surfaceType == SubsonicOutflow) {

      //InletVelocity
      fprintf(fp, "InletVelocityBC %d  %12.16e #  compID  Vnorm/c_inf\n",
              cartInstance->bcProps.surfaceProp[i].bcID,
              cartInstance->bcProps.surfaceProp[i].machNumber);


    } else if (cartInstance->bcProps.surfaceProp[i].surfaceType == SubsonicInflow) {

      //PowerStag
      fprintf(fp, "PowerStagBC %d  %12.16e  %12.16e  # compID Pt/Ptinf Tt/Ttinf\n",
              cartInstance->bcProps.surfaceProp[i].bcID,
              cartInstance->bcProps.surfaceProp[i].totalPressure,
              cartInstance->bcProps.surfaceProp[i].totalTemperature);

    } else if (cartInstance->bcProps.surfaceProp[i].surfaceType == MassflowIn) {

      //PowerMFR
      fprintf(fp, "PowerMFRBC %d  %12.16e  %12.16e  # compID  MFR Tt/Ttinf\n",
              cartInstance->bcProps.surfaceProp[i].bcID,
              cartInstance->bcProps.surfaceProp[i].massflow,
              cartInstance->bcProps.surfaceProp[i].totalTemperature);

    }
  }
  if (cartInstance->bcProps.numSurfaceProp > 0)
    fprintf(fp, "\n");

  fprintf(fp, "$__Convergence_History_reporting:\n\n");
  fprintf(fp, "iForce     %d\n", aimInputs[iForce-1].vals.integer);
  fprintf(fp, "iHist      %d\n", aimInputs[iHist-1].vals.integer);
  fprintf(fp, "nOrders    %d\n", aimInputs[nOrders-1].vals.integer);

  fprintf(fp, "\n$__Partition_Information:\n\n");
  fprintf(fp, "nPart      1\n");
  fprintf(fp, "type       1\n");

  fprintf(fp, "\n$__Post_Processing:\n\n");
  varid = Xslices-1;
  if (aimInputs[varid].nullVal != IsNull) {
    if (aimInputs[varid].length == 1) {
      fprintf(fp,"Xslices %12.16e\n",aimInputs[varid].vals.real);
    } else {
      fprintf(fp,"Xslices");
      for (j = 0; j < aimInputs[varid].length; j++) {
        fprintf(fp," %lf",aimInputs[varid].vals.reals[j]);
      }
      fprintf(fp,"\n");
    }
  }

  varid = Yslices-1;
  if (aimInputs[varid].nullVal != IsNull) {
    if (aimInputs[varid].length == 1) {
      fprintf(fp,"Yslices %12.16e\n",aimInputs[varid].vals.real);
    } else {
      fprintf(fp,"Yslices");
      for (j = 0; j < aimInputs[varid].length; j++) {
        fprintf(fp," %lf",aimInputs[varid].vals.reals[j]);
      }
      fprintf(fp,"\n");
    }
  }

  varid = Zslices-1;
  if (aimInputs[varid].nullVal != IsNull) {
    if (aimInputs[varid].length == 1) {
      fprintf(fp,"Zslices %12.16e\n",aimInputs[varid].vals.real);
    } else {
      fprintf(fp,"Zslices");
      for (j = 0; j < aimInputs[varid].length; j++) {
        fprintf(fp," %lf",aimInputs[varid].vals.reals[j]);
      }
      fprintf(fp,"\n");
    }
  }

  fprintf(fp, "#\n");
  fprintf(fp, "#lineSensor Name(%%s) Xorig(%%g) Yorig(%%g) Zorig(%%g)   Xdest(%%g) Ydest(%%g) Zdest(%%g)\n");
  fprintf(fp, "#lineSensor (Name) (float) (float) (float)   (float) (float) (float)\n");

  for (i = 0; i < nBody; i++) {
    status = EG_attributeRet(bodies[i], CART3D_SENSOR, &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_NOTFOUND) continue;
    status = EG_attributeRet(bodies[i], CART3D_EQAREA, &atype, &alen,
                             &ints, &reals, &eqstring);
    if (status == EGADS_SUCCESS) continue;
    if (aim_isNodeBody(bodies[i], xyz) == CAPS_SUCCESS) continue;
    status = EG_attributeRet(bodies[i], CART3D_SENSOR, &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_NOTFOUND) continue;

    status = EG_getBodyTopos(bodies[i], NULL, EDGE, &nedge, &eedges);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(eedges, aimInfo, status);

    status = EG_getTopology(eedges[0], &geom, &oclass, &mtype, trange, &nnode, &enodes, &senses);
    AIM_STATUS(aimInfo, status);

    status = EG_evaluate(eedges[0], &trange[0], xyz0);
    AIM_STATUS(aimInfo, status);

    status = EG_evaluate(eedges[0], &trange[1], xyz1);
    AIM_STATUS(aimInfo, status);

    AIM_FREE(eedges);

    fprintf(fp, "lineSensor  %s  %12.16e %12.16e  %12.16e   %12.16e %12.16e %12.16e\n",
            string, xyz0[0], xyz0[1], xyz0[2], xyz1[0], xyz1[1], xyz1[2]);
  }

  fprintf(fp, "#\n");
  fprintf(fp, "#                             ...Equivalent Area line sensor\n");
  fprintf(fp, "#eaSensor Name(%%s) Xorig Yorig Zorig        Xdest  Ydest  Zdest    Radius\n");
  fprintf(fp, "#eaSensor (Name) (float) (float) (float)   (float) (float) (float) (float)\n");

  for (i = 0; i < nBody; i++) {
    status = EG_attributeRet(bodies[i], CART3D_SENSOR, &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_NOTFOUND) continue;
    status = EG_attributeRet(bodies[i], CART3D_EQAREA, &atype, &alen,
                             &ints, &reals, &eqstring);
    if (status == EGADS_NOTFOUND) continue;
    if (aim_isNodeBody(bodies[i], xyz) == CAPS_SUCCESS) continue;

    status = EG_getBodyTopos(bodies[i], NULL, EDGE, &nedge, &eedges);
    AIM_STATUS(aimInfo, status);
    AIM_NOTNULL(eedges, aimInfo, status);

    status = EG_getTopology(eedges[0], &geom, &oclass, &mtype, trange, &nnode, &enodes, &senses);
    AIM_STATUS(aimInfo, status);

    status = EG_evaluate(eedges[0], &trange[0], xyz0);
    AIM_STATUS(aimInfo, status);

    status = EG_evaluate(eedges[0], &trange[1], xyz1);
    AIM_STATUS(aimInfo, status);

    AIM_FREE(eedges);

    fprintf(fp, "eaSensor  %s  %12.16e %12.16e  %12.16e   %12.16e %12.16e %12.16e    %12.16e\n",
            string, xyz0[0], xyz0[1], xyz0[2], xyz1[0], xyz1[1], xyz1[2], reals[0]);
  }

  fprintf(fp, "#\n");
  fprintf(fp, "#pointSensor Name(%%s) Xlocation(%%g) Ylocation(%%g) Zlocation(%%g)\n");

  for (i = 0; i < nBody; i++) {
    status = EG_attributeRet(bodies[i], CART3D_SENSOR, &atype, &alen,
                             &ints, &reals, &string);
    if (status == EGADS_NOTFOUND) continue;
    if (aim_isNodeBody(bodies[i], xyz) != CAPS_SUCCESS) continue;

    fprintf(fp, "pointSensor  %s  %12.16e %12.16e %12.16e\n",
            string, xyz[0], xyz[1], xyz[2]);
  }


  fprintf(fp, "\n$__Force_Moment_Processing:\n\n");
/* ... Axis definitions (with respect to body axis directions (Xb,Yb,Zb)
                         w/ usual stability and control orientation)  */
  fprintf(fp, "Model_X_axis  %s\n", aimInputs[Model_X_axis-1].vals.string);
  fprintf(fp, "Model_Y_axis  %s\n", aimInputs[Model_Y_axis-1].vals.string);
  fprintf(fp, "Model_Z_axis  %s\n", aimInputs[Model_Z_axis-1].vals.string);
  fprintf(fp, "Reference_Area   %12.16e all\n", Sref);
  fprintf(fp, "Reference_Length %12.16e all\n", Cref);
  fprintf(fp, "Force all\n\n");
  fprintf(fp, "Moment_Point %12.16e %12.16e %12.16e all\n",Xref ,Yref ,Zref );


  if (aimInputs[TargetCL-1].nullVal != IsNull) {
    //                          ...live steering this section is optional
    //                             if it exists it will get re-parsed every
    //                             iCLfreq iterations
    fprintf(fp,"\n");
    fprintf(fp, "$__Steering_Information:\n");
    //TargetCL  (CL target value) (CL tolerance)  (iCLStart iter) (iCLfrequency)
    //              (float)           (float)         (int)            (int)
    fprintf(fp, "TargetCL  %12.16e  %12.16e  %d %d\n", aimInputs[TargetCL-1].vals.real,
                                                       aimInputs[TargetCL_Tol-1].vals.real,
                                                       aimInputs[TargetCL_Start_Iter-1].vals.integer,
                                                       aimInputs[TargetCL_Freq-1].vals.integer);
    // NominalAlphaStep  %f  # (OPTIONAL) choose Initial alpha step, also max step size OPTIONAL, default <0.2>
    fprintf(fp, "NominalAlphaStep %12.16e\n", aimInputs[TargetCL_NominalAlphaStep-1].vals.real);
  }

  fclose(fp); fp = NULL;

  // overwrite preSpec.c3d.cntl generated by autoInputs
  if (makePreSpec == (int)true) {
    fp = aim_fopen(aimInfo, "inputs/preSpec.c3d.cntl", "w");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open: inputs/preSpec.c3d.cntl!");
      status = CAPS_IOERR;
      goto cleanup;
    }
/*
# FORMAT:
#    (min and maxes are defined in the same space as the mesh is generated)
#
# BBox: level   Xmin   Xmax      Ymin   Ymax      Zmin    Zmax
#      (int)  (float) (float)   (float) (float)  (float)  (float)
#
# XLev: NumXtraLevels Comp  Comp  Comp ...
# XLev:  (int)        (int) (int) (int) ...
*/

    fprintf(fp, "$__Prespecified_Adaptation_Regions:\n");
    for (i = 0; i < nBody; i++) {
      status = EG_attributeRet(bodies[i], CART3D_BBOX, &atype, &alen,
                               &ints, &reals, &string);
      if (status == EGADS_NOTFOUND) continue;

      level = atype == ATTRINT ? ints[0] : reals[0] + 0.00001;

      status = EG_getBoundingBox(bodies[i], box);
      AIM_STATUS(aimInfo, status);

      fprintf(fp, "BBox  %d  %12.16e %12.16e   %12.16e %12.16e   %12.16e %12.16e\n",
              level, box[0], box[3], box[1], box[4], box[2], box[5]);

    }
  }

  status = CAPS_SUCCESS;

cleanup:

  if (fp != NULL) fclose(fp);
  xddm_free(p_xddm, options);
  AIM_FREE(eedges);

  return status;
}


// ********************** AIM Function Break *****************************
int aimExecute(/*@unused@*/ const void *instStore, void *aimInfo,
               int *state)
{
  /*! \page aimExecuteCART3D AIM Execution
   *
   * If auto execution is enabled when creating an Cart3D AIM,
   * the AIM will execute Cart3D just-in-time with the command line:
   *
   * \code{.sh}
   * ./aero.csh > aero.out
   * \endcode
   *
   * in the "inputs" analysis directory or
   *
   * \code{.sh}
   * c3d_objGrad.csh
   * \endcode
   *
   * in the "design" analysis directory when computing functional sensitvities w.r.t. design variables.
   * The AIM preAnalysis generates the inputs/aero.csh script required to execute Cart3D by
   * copying it from $CART3D/bin/aero.csh and inserting various analysis input modifications.
   *
   * The analysis can be also be explicitly executed with caps_execute in the C-API
   * or via Analysis.runAnalysis in the pyCAPS API.
   *
   * Calling preAnalysis and postAnalysis is NOT allowed when auto execution is enabled.
   *
   * Auto execution can also be disabled when creating an Cart3D AIM object.
   * In this mode, caps_execute and Analysis.runAnalysis can be used to run the analysis,
   * or Cart3D can be executed by calling preAnalysis, system, and posAnalysis as demonstrated
   * below with a pyCAPS example:
   *
   * \code{.py}
   * print ("\n\preAnalysis......")
   * cart.preAnalysis()
   *
   * print ("\n\nRunning......")
   * cart.system("./aero.csh", "inputs"); # Run via system call in inputs analysis directory
   *
   * print ("\n\postAnalysis......")
   * cart.postAnalysis()
   * \endcode
   */
  int status = CAPS_SUCCESS;
  capsValue *design_functional = NULL;
  capsValue *aerocsh_args = NULL;
  const char* const objgcmd = "c3d_objGrad.csh";
  char aerocmd[256];

  *state = 0;

  aim_getValue(aimInfo, Design_Functional, ANALYSISIN, &design_functional);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(design_functional, aimInfo, status);

  aim_getValue(aimInfo, aerocsh_arg, ANALYSISIN, &aerocsh_args);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(aerocsh_args, aimInfo, status);

  if (design_functional->nullVal == NotNull) {
    printf(" Executing: %s\n", objgcmd);
    status = aim_system(aimInfo, runDir, objgcmd);
  } else {
    snprintf(aerocmd, 256, "./aero.csh %s > aero.out", aerocsh_args->vals.string);
    printf(" Executing: %s\n", aerocmd);
    status = aim_system(aimInfo, "inputs", aerocmd);
    if (status != CAPS_SUCCESS) {
      aim_system(aimInfo, "inputs", "cat aero.out");
    }
  }

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************

/* no longer optional and needed for restart */
int aimPostAnalysis(void *instStore, void *aimInfo,
       /*@unused@*/ int restart, capsValue *inputs)
{
  int status = CAPS_SUCCESS; // Function return status

  int i, j, idv, irow, icol; // Indexing
  int index, count, indexed;
  size_t len;

  unsigned options = 0;
  p_tsXddm p_xddm = NULL;
  char aimFile[PATH_MAX], tmp[128];

  const char *name;
  char *str;

  FILE *fp=NULL;
  size_t     linecap = 0;
  const char comp = ':';
  char       *valstr, *line = NULL;

  capsValue value, *inVal;
  aimStorage *cartInstance;

#define NUMCOEF 12
  const char *startLoads[NUMCOEF]= { "   Axial Force (C_A):",
                                     " Lateral Force (C_Y):",
                                     "  Normal Force (C_N):",
                                     "    Drag Force (C_D):",
                                     "    Side Force (C_S):",
                                     "    Lift Force (C_L):",
                                     "  Rolling Moment",
                                     " Pitching Moment",
                                     "   Yawing Moment",
                                     "   X Aero Moment",
                                     "   Y Aero Moment",
                                     "   Z Aero Moment" };

  const char *names[NUMCOEF] = {"C_A"  , "C_Y"  , "C_N",
                                "C_D"  , "C_S"  , "C_L",
                                "C_l"  , "C_m"  , "C_n",
                                "C_M_x", "C_M_y", "C_M_z"};


  cartInstance = (aimStorage*)instStore;

  AIM_NOTNULL(inputs, aimInfo, status);
  if (inputs[Design_Functional-1].nullVal == NotNull) {

    aim_initValue(&value);

    // Read design.xml
    snprintf(tmp, 128, "%s/design.xml", runDir);
    status = aim_file(aimInfo, tmp, aimFile);
    AIM_STATUS(aimInfo, status);

    p_xddm = xddm_readFile(aimFile, "/Optimize", &options);
    if (p_xddm == NULL) {
      AIM_ERROR(aimInfo, "xddm_readFile failed to read: %s", aimFile);
      status = CAPS_IOERR;
      goto cleanup;
    }

    if (p_xddm->nk != 1) {
      AIM_ERROR(aimInfo, "'%s' should only have 1 design point, found %d", aimFile, p_xddm->nk);
      status = CAPS_IOERR;
      goto cleanup;
    }

    /* Loop over analysis params */
    for (i = 0; i < p_xddm->a_kids[0].na; i++) {

      /* get the name and value */
      name = p_xddm->a_kids[0].a_ap[i].p_e->p_id;
      value.vals.real = p_xddm->a_kids[0].a_ap[i].val;

      /* look for derivatives */
      if (p_xddm->a_kids[0].a_ap[i].ndvs > 0) {
        value.type = DoubleDeriv;

        /* allocate derivatives */
        AIM_ALLOC(value.derivs, cartInstance->design.numDesignVariable, capsDeriv, aimInfo, status);
        for (idv = 0; idv < cartInstance->design.numDesignVariable; idv++) {
          value.derivs[idv].name  = NULL;
          value.derivs[idv].deriv = NULL;
          value.derivs[idv].len_wrt  = 0;
        }
        value.nderiv = cartInstance->design.numDesignVariable;

        /* set derivatives */
        for (j = 0; j < p_xddm->a_kids[0].a_ap[i].ndvs; j++) {
          str = strstr(p_xddm->a_kids[0].a_ap[i].pa_dvs[j], "__Variable__");
          AIM_NOTNULL(str, aimInfo, status);

          str += strlen("__Variable__");

          // find the variable and check if it has an index
          indexed = (int)false;
          for (idv = 0; idv < cartInstance->design.numDesignVariable; idv++) {
            len = strlen(cartInstance->design.designVariable[idv].name);
            if (strncmp(cartInstance->design.designVariable[idv].name, str, len) == 0) {
              if (strlen(str) > len)
                indexed = (int)true;
              break;
            }
          }

          AIM_FREE(value.derivs[idv].name);
          AIM_STRDUP(value.derivs[idv].name, str, aimInfo, status);

          irow = icol = 0;

          /* Extract indexing of the variable */
          if (indexed == (int)true) {
            str = value.derivs[idv].name + strlen(value.derivs[idv].name);
            count = 0;
            while(count < 3) {
              if (str[0] == '_') {
                str[0] = ' ';
                count++;
              }
              str--;
            }

            status = sscanf(value.derivs[idv].name, "%s %d %d", tmp, &irow, &icol);
            if (status != 3) {
              AIM_ERROR(aimInfo, "Failed to parse: %s", value.derivs[idv].name);
              status = CAPS_IOERR;
              goto cleanup;
            }
            irow--; // 0-based indexing
            icol--;
            value.derivs[idv].name[strlen(tmp)] = '\0';
          }

          // Get the value object so the dimension can be extracted
          if (cartInstance->design.designVariable[idv].type == DesignVariableGeometry) {
            index = aim_getIndex(aimInfo, value.derivs[idv].name, GEOMETRYIN);
            status = aim_getValue(aimInfo, index, GEOMETRYIN, &inVal);
            AIM_STATUS(aimInfo, status);
          } else {
            index = aim_getIndex(aimInfo, value.derivs[idv].name, ANALYSISIN);
            status = aim_getValue(aimInfo, index, ANALYSISIN, &inVal);
            AIM_STATUS(aimInfo, status);
          }

          if (value.derivs[idv].deriv == NULL) {
            AIM_ALLOC(value.derivs[idv].deriv, inVal->length, double, aimInfo, status);
            value.derivs[idv].len_wrt  = inVal->length;
          }

          value.derivs[idv].deriv[inVal->ncol*irow + icol] = p_xddm->a_kids[0].a_ap[i].a_lin[j];
        }

      } else {
        value.type = Double;
      }

      /* create the dynamic output */
      status = aim_makeDynamicOutput(aimInfo, name, &value);
      AIM_STATUS(aimInfo, status);
    }
  }


  if (inputs[Design_Functional-1].nullVal == IsNull) {

    /* open the Cart3D loads file */
    fp = aim_fopen(aimInfo, "inputs/BEST/FLOW/loadsCC.dat", "r");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open: inputs/BEST/FLOW/loadsCC.dat!");
      status = CAPS_IOERR;
      goto cleanup;
    }

    aim_initValue(&value);

    for (i = 0; i < cartInstance->groupMap.numAttribute; i++) {

      value.type = Tuple;
      value.dim  = Vector;
      AIM_ALLOC(value.vals.tuple, NUMCOEF, capsTuple, aimInfo, status);
      value.length = value.nrow = NUMCOEF;
      for (j = 0; j < value.nrow; j++) {
        value.vals.tuple[j].name = NULL;
        value.vals.tuple[j].value = NULL;
      }

      for (index = 0; index < NUMCOEF; index++) {
        /* scan the file for the string */
        valstr = NULL;
        while (getline(&line, &linecap, fp) >= 0) {
          if (line == NULL) continue;
          if (line[0] == '#') continue;
          if (strncmp(line, cartInstance->groupMap.attributeName[i], linecap) == 0) continue;

          valstr = strstr(line, startLoads[index]);
          if (valstr != NULL) {
            valstr = strchr(line,comp);

#ifdef DEBUG
            printf("valstr > |%s|\n",valstr);
#endif
            /* found it -- get the value */
            AIM_STRDUP(value.vals.tuple[index].name , names[index], aimInfo, status);
            AIM_STRDUP(value.vals.tuple[index].value, valstr+1    , aimInfo, status);
            break;
          }
        }

        if (valstr == NULL) {
          AIM_ERROR(aimInfo, "Cannot find '%s %s' in inputs/BEST/FLOW/loadsCC.dat file!",
                    cartInstance->groupMap.attributeName[i], startLoads[index]);
          status = CAPS_NOTFOUND;
          goto cleanup;
        }
      }

      //rewind the file
      rewind(fp);
      free(line); line = NULL;

      /* create the dynamic output */
      status = aim_makeDynamicOutput(aimInfo, cartInstance->groupMap.attributeName[i], &value);
      AIM_STATUS(aimInfo, status);
    }
  }


cleanup:

  xddm_free(p_xddm, options);

  if (fp != NULL) fclose(fp);
  if (line != NULL) free(line);

  return status;
}


// ********************** AIM Function Break *****************************
int
aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo, int index,
           char **aoname, capsValue *form)
{
  int status = CAPS_SUCCESS;
  const char *names[NUMOUT] = {"C_A"  , "C_Y"  , "C_N",
                               "C_D"  , "C_S"  , "C_L",
                               "C_l"  , "C_m"  , "C_n",
                               "C_M_x", "C_M_y", "C_M_z",
                               "alpha"};
#ifdef DEBUG
  printf(" Cart3DAIM/aimOutputs index = %d!\n", index);
#endif

  /*! \page aimOutputsCART3D AIM Outputs
  * Integrated force outputs on the entire body are available as outputs from the loadsCC.dat output file.
  * - <B> C_A</B>. entire Axial Force
  * - <B> C_Y</B>. entire Lateral Force
  * - <B> C_N</B>. entire Normal Force
  * - <B> C_D</B>. entire Drag Force
  * - <B> C_S</B>. entire Side Force
  * - <B> C_L</B>. entire Lift Force
  * - <B> C_l</B>. entire  Rolling Moment
  * - <B> C_m</B>. Pitching Moment
  * - <B> C_n</B>. Yawing Moment
  * - <B> C_M_x</B>. X Aero Moment
  * - <B> C_M_y</B>. Y Aero Moment
  * - <B> C_M_z</B>. Z Aero Moment
  * - <B> alpha</B>. Angle of attach (may change using TargetCL)
  */

  AIM_STRDUP(*aoname, names[index-1], aimInfo, status);
  form->type = Double;

cleanup:
  return status;
}


// ********************** AIM Function Break *****************************
int
aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, capsValue *val)
{
  int        status, i, n;
  size_t     linecap = 0, headercap = 0;
  const char comp = ':';
  char       *valstr, *restHeader, *restLine, *token, *line = NULL, *header=NULL, tmp[128];
  FILE       *fp=NULL;
  const char *startLoads[NUMOUT]= { "entire   Axial Force (C_A):",
                                    "entire Lateral Force (C_Y):",
                                    "entire  Normal Force (C_N):",
                                    "entire    Drag Force (C_D):",
                                    "entire    Side Force (C_S):",
                                    "entire    Lift Force (C_L):",
                                    "entire  Rolling Moment",
                                    "entire Pitching Moment",
                                    "entire   Yawing Moment",
                                    "entire   X Aero Moment",
                                    "entire   Y Aero Moment",
                                    "entire   Z Aero Moment",
                                    "alpha"};

  const char *startTbl[NUMOUT]= { "CAxial",
                                  "CY(Lateral)",
                                  "CNormal",
                                  "CDrag",
                                  "CSide",
                                  "CLift",
                                  "Cl(Roll)",
                                  "Cm(Pitch)",
                                  "Cn(Yaw)",
                                  "C_M_x",
                                  "C_M_y",
                                  "C_M_z",
                                  "Alpha"};

  capsValue *design_functional=NULL;

  status = aim_getValue(aimInfo, Design_Functional, ANALYSISIN, &design_functional);
  AIM_STATUS(aimInfo, status);
  AIM_NOTNULL(design_functional, aimInfo, status);

  if (design_functional->nullVal == NotNull) {

    // open the tbl file
    snprintf(tmp, 128, "%s/entire.DP1.tbl", runDir);
    fp = aim_fopen(aimInfo, tmp, "r");
    if (fp == NULL) {
      AIM_ERROR(aimInfo, "Cannot open: %s!", tmp);
      return CAPS_IOERR;
    }

    // The table file looks like:
    //
    //# Fri Oct 8 11:13:39 EDT 2021
    //# 1.DesignIter 2.Beta 3.Mach 4.Alpha 5.CAxial 6.CY(Lateral) 7.CNormal 8.CDrag 9.CSide 10.CLift 11.L/D 12.Cl(Roll) 13.Cm(Pitch) 14.Cn(Yaw) 15.nCells 16.maxRef
    //#
    //  0   0.00000000   0.50000000   1.00000000   0.02693920   0.00062873  -0.01373210   0.02669544   0.00062873  -0.01420016  -0.53193196   0.04208490  -0.01011480   0.04414300    10424     7

    if (getline(&line, &linecap, fp) < 0) {
      status = CAPS_IOERR;
      AIM_STATUS(aimInfo, status);
    }

    if (getline(&header, &headercap, fp) < 0) {
      status = CAPS_IOERR;
      AIM_STATUS(aimInfo, status);
    }

    for (i = 0; i < 2; i++)
      if (getline(&line, &linecap, fp) < 0) {
        status = CAPS_IOERR;
        AIM_STATUS(aimInfo, status);
      }

    restHeader = header+2; // skip '# '
    restLine = line;
    while ((token = strtok_r(restHeader, " ", &restHeader))) {
      sscanf(token, "%d.%s", &n, tmp);
      valstr = strtok_r(restLine, " ", &restLine);
      if (valstr == NULL) continue;
      if (strncasecmp(startTbl[index-1], tmp, strlen(startTbl[index-1])) == 0) {
        /* found it -- get the value */
        sscanf(valstr, "%lf", &val->vals.real);
        status = CAPS_SUCCESS;
        goto cleanup;
      }
    }
    AIM_ERROR(aimInfo, "Cannot find '%s' in entire.DP1.tbl file!", startTbl[index-1]);
    status = CAPS_NOTFOUND;

  } else {

    if (index == NUMOUT) { // alpha
      /* open the Cart3D input.cntl file */
      fp = aim_fopen(aimInfo, "inputs/BEST/FLOW/input.cntl", "r");
      if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open: inputs/BEST/FLOW/input.cntl!");
        return CAPS_IOERR;
      }
    } else {
      /* open the Cart3D loads file */
      fp = aim_fopen(aimInfo, "inputs/BEST/FLOW/loadsCC.dat", "r");
      if (fp == NULL) {
        AIM_ERROR(aimInfo, "Cannot open: inputs/BEST/FLOW/loadsCC.dat!");
        return CAPS_IOERR;
      }
    }

    /* scan the file for the string */
    valstr = NULL;
    while (getline(&line, &linecap, fp) >= 0) {
      if (line == NULL) continue;
      if (line[0] == '#') continue;
      valstr = strstr(line, startLoads[index-1]);
      if (valstr != NULL) {
        if (index == NUMOUT)
          valstr += strlen(startLoads[index-1]);
        else
          valstr = strchr(line,comp);


#ifdef DEBUG
        printf("valstr > |%s|\n",valstr);
#endif
        /* found it -- get the value */
        sscanf(valstr+1, "%lf", &val->vals.real);
        break;
      }
    }
    free(line); line = NULL;

    if (valstr == NULL) {
      if (index == NUMOUT) { // alpha
        AIM_ERROR(aimInfo, "Cannot find '%s' in inputs/BEST/FLOW/input.cntl file!", startLoads[index-1]);
      } else {
        AIM_ERROR(aimInfo, "Cannot find '%s' in inputs/BEST/FLOW/loadsCC.dat file!", startLoads[index-1]);
      }
      status = CAPS_NOTFOUND;
      goto cleanup;
    }
  }


cleanup:
  if (line != NULL) free(line);
  if (header != NULL) free(header);
  if (fp != NULL) fclose(fp);

  return status;
}


// ********************** AIM Function Break *****************************
void
aimCleanup(void *instStore)
{
#ifdef DEBUG
  printf(" Cart3DAIM/aimCleanup!\n");
#endif

  int i;
  aimStorage *cartInstance = NULL;

  cartInstance = (aimStorage*)instStore;

  destroy_cfdDesignFunctionalStruct(cartInstance->adaptFunctional);
  AIM_FREE(cartInstance->adaptFunctional);

  // Design information
  destroy_cfdDesignStruct(&cartInstance->design);

  destroy_mapAttrToIndexStruct(&cartInstance->groupMap);

  // BC properties
  destroy_cfdBoundaryConditionStruct(&cartInstance->bcProps);

  for (i = 0; i < cartInstance->ntess; i++)
    EG_deleteObject(cartInstance->tess[i]);
  AIM_FREE(cartInstance->tess);

  AIM_FREE(instStore);
}


// ********************** AIM Function Break *****************************
int
aimDiscr(char *tname, capsDiscr *discr)
{
  int           ibody, iface, i, vID=0, itri, ielem, status, found;
  int           nFace, atype, alen, tlen, nBodyDisc;
  int           ntris, nBody, state, nGlobal, global;
  int           *vid = NULL;
  const int     *ints, *ptype, *pindex, *tris, *nei;
  const double  *reals, *xyz, *uv;
  const char    *string, *intents;
  ego           body, *bodies, *faces=NULL, *tess = NULL;
  capsBodyDiscr *discBody;
  aimStorage    *cartInstance = NULL;

#ifdef DEBUG
  printf(" capsAIM/aimDiscr: tname = %s, instance = %d!\n",
         tname, aim_getInstance(discr->aInfo));
#endif

  /* create the discretization structure for this capsBound */

  status = aim_getBodies(discr->aInfo, &intents, &nBody, &bodies);
  AIM_STATUS(discr->aInfo, status);

  cartInstance = (aimStorage*)discr->instStore;
  tess = cartInstance->tess;

  /* find any faces with our boundary marker */
  for (nBodyDisc = ibody = 0; ibody < nBody; ibody++) {
    if (tess[ibody] == NULL) continue;

    status = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nFace, &faces);
    AIM_STATUS(discr->aInfo, status, "getBodyTopos (Face) for Body %d",
               ibody+1);
    AIM_NOTNULL(faces, discr->aInfo, status);

    found = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;
#ifdef DEBUG
      printf(" cart3dAIM/aimDiscr: Body %d/Face %d matches %s!\n",
             ibody+1, iface+1, tname);
#endif
      /* count the number of face with capsBound */
      found = 1;
    }
    if (found == 1) nBodyDisc++;
    AIM_FREE(faces);
  }
  if (nBodyDisc == 0) {
    printf(" cart3dAIM/aimDiscr: No Faces match %s!\n", tname);
    return CAPS_SUCCESS;
  }

  /* specify our single triangle element type */
  discr->nTypes = 1;
  AIM_ALLOC(discr->types, discr->nTypes, capsEleType, discr->aInfo, status);

  /* define triangle element type */
  status = aim_nodalTriangleType( &discr->types[0]);
  AIM_STATUS(discr->aInfo, status);

  /* allocate the body discretizations */
  discr->nBodys = nBodyDisc;
  AIM_ALLOC(discr->bodys, discr->nBodys, capsBodyDiscr, discr->aInfo, status);

  /* get the tessellation and
   make up a linear continuous triangle discretization */
  vID = nBodyDisc = 0;
  for (ibody = 0; ibody < nBody; ibody++) {
    if (tess[ibody] == NULL) continue;

    status = EG_statusTessBody(tess[ibody], &body, &state, &nGlobal);
    AIM_STATUS(discr->aInfo, status);

    ntris = 0;
    AIM_FREE(faces);
    status = EG_getBodyTopos(body, NULL, FACE, &nFace, &faces);
    AIM_STATUS(discr->aInfo, status, "getBodyTopos (Face) for Body %d", ibody+1);
    AIM_NOTNULL(faces, discr->aInfo, status);

    found = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;

      status = EG_getTessFace(tess[ibody], iface+1, &alen, &xyz, &uv,
                              &ptype, &pindex, &tlen, &tris, &nei);
      if (status != EGADS_SUCCESS) {
        printf(" cart3dAIM: EG_getTessFace %d = %d for Body %d!\n",
               iface+1, status, ibody+1);
        continue;
      }
      ntris += tlen;
      found = 1;
    }
    if (found == 0) continue;
    if (ntris == 0) {
      AIM_ERROR(discr->aInfo, "No faces with capsBound = %s", tname);
      status = CAPS_SOURCEERR;
      goto cleanup;
    }

    discBody = &discr->bodys[nBodyDisc];
    aim_initBodyDiscr(discBody);

    discBody->tess = tess[ibody];
    discBody->nElems = ntris;

    AIM_ALLOC(discBody->elems   ,   ntris, capsElement, discr->aInfo, status);
    AIM_ALLOC(discBody->gIndices, 6*ntris, int        , discr->aInfo, status);

    AIM_FREE(vid);
    AIM_ALLOC(vid, nGlobal, int, discr->aInfo, status);
    for (i = 0; i < nGlobal; i++) vid[i] = 0;

    ielem = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;

      status = EG_getTessFace(tess[ibody], iface+1, &alen, &xyz, &uv,
                            &ptype, &pindex, &tlen, &tris, &nei);
      AIM_STATUS(discr->aInfo, status);

      /* construct global vertex indices */
      for (i = 0; i < alen; i++) {
        status = EG_localToGlobal(tess[ibody], iface+1, i+1, &global);
        AIM_STATUS(discr->aInfo, status);
        if (vid[global-1] != 0) continue;
        vid[global-1] = vID+1;
        vID++;
      }

      /* fill elements */
      for (itri = 0; itri < tlen; itri++, ielem++) {
        discBody->elems[ielem].tIndex      = 1;
        discBody->elems[ielem].eIndex      = iface+1;
/*@-immediatetrans@*/
        discBody->elems[ielem].gIndices    = &discBody->gIndices[6*ielem];
/*@+immediatetrans@*/
        discBody->elems[ielem].dIndices    = NULL;
        discBody->elems[ielem].eTris.tq[0] = itri+1;

        status = EG_localToGlobal(tess[ibody], iface+1, tris[3*itri  ],
                                  &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[0] = vid[global-1];
        discBody->elems[ielem].gIndices[1] = tris[3*itri  ];

        status = EG_localToGlobal(tess[ibody], iface+1, tris[3*itri+1],
                                  &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[2] = vid[global-1];
        discBody->elems[ielem].gIndices[3] = tris[3*itri+1];

        status = EG_localToGlobal(tess[ibody], iface+1, tris[3*itri+2],
                                  &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[4] = vid[global-1];
        discBody->elems[ielem].gIndices[5] = tris[3*itri+2];
      }
    }
    nBodyDisc++;
  }

  /* set the total number of points */
  discr->nPoints = vID;

  status = CAPS_SUCCESS;

cleanup:

  /* free up our stuff */
  AIM_FREE(faces);
  AIM_FREE(vid);

  return status;
}


// ********************** AIM Function Break *****************************
int
aimLocateElement(capsDiscr *discr, double *params, double *param,
                 int *bIndex, int *eIndex, double *bary)
{
    return aim_locateElement(discr, params, param, bIndex, eIndex, bary);
}


// ********************** AIM Function Break *****************************
int
aimTransfer(capsDiscr *discr, const char *name, int npts, int rank, double *data,
            /*@unused@*/ char **units)
{
  /*! \page dataTransferCART3D AIM Data Transfer
   *
   * The Cart3D AIM has the ability to transfer surface data (e.g. pressure distributions) to and from the AIM
   * using the conservative and interpolative data transfer schemes in CAPS. Currently these transfers may only
   * take place on triangular meshes.
   *
   * \section dataFromCart3D Data transfer from Cart3D (FieldOut)
   *
   * <ul>
   * <li> <B>"Pressure" </B> </li> <br>
   *  Loads the pressure distribution from Components.i.trix file.
   *  This distribution may be scaled based on
   *  Pressure = Pressure_Scale_Factor*Pressure, where "Pressure_Scale_Factor"
   *  is an AIM input (\ref aimInputsCART3D)
   * </ul>
   *
   */
  int    i, j, global, status, bIndex, dim, nBody;
  double **rvec=NULL, scale = 1.0;
  capsValue *Pressure_Scale_Factor_Value=NULL;
  char aimFile[PATH_MAX];
  const char *intents;
  ego        *bodies;

#ifdef DEBUG
  printf(" Cart3DAIM/aimTransfer name = %s  npts = %d/%d!\n",
         name, npts, len_wrt);
#endif

  status = aim_getBodies(discr->aInfo, &intents, &nBody, &bodies);
  AIM_STATUS(discr->aInfo, status);

  status = aim_file(discr->aInfo, "inputs/BEST/FLOW/Components.i.trix", aimFile);
  AIM_STATUS(discr->aInfo, status);

  /* try and read the trix file */
  status = readTrix(aimFile, name, &dim, &rvec);
  if (status != 0) {
    return CAPS_IOERR;
  }
  if (dim != rank) {
    AIM_ERROR(discr->aInfo, "Rank/Dim missmatch!");
    return CAPS_IOERR;
  }
  AIM_NOTNULL(rvec, discr->aInfo, status);

  if (strcmp(name, "Pressure") == 0) {
    status = aim_getValue(discr->aInfo, Pressure_Scale_Factor, ANALYSISIN, &Pressure_Scale_Factor_Value);
    AIM_STATUS(discr->aInfo, status);
    AIM_NOTNULL(Pressure_Scale_Factor_Value, discr->aInfo, status);
    scale = Pressure_Scale_Factor_Value->vals.real;
  }

  /* move the appropriate parts of the tessellation to data */
  for (i = 0; i < npts; i++) {
    /* points might span multiple bodies */
    bIndex = discr->tessGlobal[2*i  ];
    global = discr->tessGlobal[2*i+1];
    for (j = 0; j < rank; j++) data[rank*i+j] = rvec[bIndex-1][rank*(global-1)+j] * scale;
  }

  status = CAPS_SUCCESS;
cleanup:
  if (rvec != NULL) {
    for (i = 0; i < nBody; ++i) {
      AIM_FREE(rvec[i]);
    }
  }
  AIM_FREE(rvec);
  return status;
}


// ********************** AIM Function Break *****************************
int
aimInterpolation(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                 int eIndex, double *bary, int rank, double *data,
                 double *result)
{
#ifdef DEBUG
    printf(" Cart3DAIM/aimInterpolation  %s!\n", name);
#endif

    return  aim_interpolation(discr, name, bIndex, eIndex,
                              bary, rank, data, result);
}


// ********************** AIM Function Break *****************************
int
aimInterpolateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                  int eIndex, double *bary, int rank, double *r_bar,
                  double *d_bar)
{
#ifdef DEBUG
    printf(" Cart3DAIM/aimInterpolateBar  %s!\n", name);
#endif

    return  aim_interpolateBar(discr, name, bIndex, eIndex,
                               bary, rank, r_bar, d_bar);
}


// ********************** AIM Function Break *****************************
int
aimIntegration(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
               int eIndex, int rank, double *data, double *result)
{
#ifdef DEBUG
    printf(" Cart3DAIM/aimIntegration  %s!\n", name);
#endif

    return aim_integration(discr, name, bIndex, eIndex, rank,
                           data, result);
}


// ********************** AIM Function Break *****************************
int
aimIntegrateBar(capsDiscr *discr, /*@unused@*/ const char *name, int bIndex,
                int eIndex, int rank, double *r_bar, double *d_bar)
{
#ifdef DEBUG
    printf(" Cart3DAIM/aimIntegrateBar  %s!\n", name);
#endif

    return aim_integrateBar(discr, name, bIndex, eIndex, rank,
                            r_bar, d_bar);
}
