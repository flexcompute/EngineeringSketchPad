/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Cart3D AIM
 *
 *      Copyright 2014-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <string.h>
#include <math.h>

#include "bodyTess.h"
#include "writeTrix.h"
#include "capsTypes.h"
#include "aimUtil.h"
#include "miscUtils.h"

#ifdef WIN32
#define snprintf   _snprintf
#define strcasecmp stricmp
#endif

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


//#define DEBUG

enum aimInputs
{
  Tess_Params = 1,                   /* index is 1-based */
  outer_box,
  nDiv,
  maxR,
  Mach,
  alpha,
  beta,
  Gamma,
  maxCycles,
  SharpFeatureDivisions,
  nMultiGridLevels,
  MultiGridCycleType,
  MultiGridPreSmoothing,
  MultiGridPostSmoothing,
  CFL,
  Limiter,
  FluxFun,
  iForce,
  iHist,
  nOrders,
  Xslices,
  Yslices,
  Zslices,
  Model_X_axis,
  Model_Y_axis,
  Model_Z_axis,
  NUMINPUT = Model_Z_axis            /* Total number of inputs */
};

#define NUMOUT    12                 /* number of outputs */


/* currently setup for a single body */
typedef struct {
  int  nface;
  int  nvert;
  int  ntris;

  // Meshing parameters
  double tessParam[3];
  double outerBox;
  int nDiv, maxR, sharpFeatureDivisions, nMultiGridLevels;

} c3dAIM;

/*!\mainpage Introduction
 *
 * \section overviewCART3D CART3D AIM Overview
 * CART3D
 *
 * \section assumptionsCART3D Assumptions
 *
 * This documentation contains four sections to document the use of the CART3D AIM.
 * \ref examplesCART3D contains example *.csm input files and pyCAPS scripts designed to make use of the CART3D AIM.  These example
 * scripts make extensive use of the \ref attributeCART3D and CART3D \ref aimInputsCART3D and CART3D \ref aimOutputsCART3D.
 *
 * \section dependeciesCART3D Dependencies
 * ESP client of libxddm. For XDDM documentation, see $CART3D/doc/xddm/xddm.html. The library uses XML Path Language (XPath) to
 * navigate the elements of XDDM documents. For XPath tutorials, see the web,
 * e.g. www.developer.com/net/net/article.php/3383961/NET-and-XML-XPath-Queries.htm
 *
 * Dependency: libxml2, www.xmlsoft.org. This library is usually present on
 * most systems, check existence of 'xml2-config' script
 */

/*! \page examplesCART3D CART3D Examples
 *
 * This example contains a set of *.csm and pyCAPS (*.py) inputs that uses the CART3D AIM.  A user should have knowledge
 * on the generation of parametric geometry in Engineering Sketch Pad (ESP) before attempting to integrate with any AIM.
 * Specifically, this example makes use of Design Parameters, Set Parameters, User Defined Primitive (UDP) and attributes
 * in ESP.
 *
 */

/*! \page attributeCART3D CART3D attributes
 * The following list of attributes drives the CART3D geometric definition.
 *
 * - <b> capsAIM</b> This attribute is a CAPS requirement to indicate the analysis the geometry
 * representation supports.
 *
 *  - <b> capsReferenceArea</b>  [Optional: Default 1.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the Reference_Area entry in the CART3D input.
 *
 *  - <b> capsReferenceChord</b>  [Optional: Default 1.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used as the Reference_Length entry in the CART3D input.
 *
 *  - <b> capsReferenceX</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used in the Moment_Point entry in the CART3D input.
 *
 *  - <b> capsReferenceY</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 * value will be used in the Moment_Point entry in the CART3D input.
 *
 *  - <b> capsReferenceZ</b>  [Optional: Default 0.0] This attribute may exist on any <em> Body</em>.  Its
 *  value will be used in the Moment_Point entry in the CART3D input.
 *
 */




/* ********************** Exposed AIM Functions ***************************** */

int aimInitialize(int inst, /*@unused@*/ const char *unitSys, void *aimInfo,
                  /*@unused@*/ void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int    *ints = NULL, i, status = CAPS_SUCCESS;
  char   **strs = NULL;
  c3dAIM *cartInstance = NULL;

#ifdef DEBUG
  printf("\n Cart3DAIM/aimInitialize   instance = %d!\n", inst);
#endif

  /* specify the number of analysis input and out "parameters" */
  *nIn     = NUMINPUT;
  *nOut    = NUMOUT;
  if (inst == -1) return CAPS_SUCCESS;

  /* get the storage */
  AIM_ALLOC(cartInstance, 1, c3dAIM, aimInfo, status);
  *instStore = cartInstance;

  cartInstance->nface = 0;
  cartInstance->nvert = 0;
  cartInstance->ntris = 0;

  /* Set meshing parameters  */
  cartInstance->tessParam[0] = 0.025;
  cartInstance->tessParam[1] = 0.001;
  cartInstance->tessParam[2] = 15.00;
  cartInstance->outerBox = 30;
  cartInstance->nDiv = 5;
  cartInstance->maxR = 11;
  cartInstance->sharpFeatureDivisions = 2;
  cartInstance->nMultiGridLevels = 1;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 4;

  /* specify the name of each field variable */
  AIM_ALLOC(strs, *nFields, char *, aimInfo, status);
  strs[0]  = EG_strdup("Cp");
  strs[1]  = EG_strdup("Density");
  strs[2]  = EG_strdup("Velocity");
  strs[3]  = EG_strdup("Pressure");
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
  *franks   = ints;
  ints = NULL;

  /* specify if a field is an input field or output field */
  AIM_ALLOC(ints, *nFields, int, aimInfo, status);

  ints[0]  = FieldOut;
  ints[1]  = FieldOut;
  ints[2]  = FieldOut;
  ints[3]  = FieldOut;
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
     * The following list outlines the CART3D inputs along with their default value available
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
     * These parameters are used to create the surface mesh for CART3D. Their order definition is as follows.
     *  1. Max Edge Length (0 is any length)
     *  2. Max Sag or distance from mesh segment and actual curved geometry
     *  3. Max angle in degrees between triangle facets
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
    defval->vals.integer = 11;
      /*! \page aimInputsCART3D
     * - <B> maxR = int</B>. <Default 11> <br>
     * Max Num of cell refinements to perform
     */
  } else if (index == Mach) {
    *name             = EG_strdup("Mach");
    defval->type      = Double;
    defval->vals.real = 0.76;
      /*! \page aimInputsCART3D
     * - <B> Mach = double</B>. <Default 0.76> <br>
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
  } else if (index == maxCycles) {
    *name                = EG_strdup("maxCycles");
    defval->type         = Integer;
    defval->vals.integer = 1000;
      /*! \page aimInputsCART3D
     * - <B> maxCycles = int</B>. <Default 1000> <br>
     * Number of iterations
     */
  } else if (index == SharpFeatureDivisions) {
    *name                = EG_strdup("SharpFeatureDivisions");
    defval->type         = Integer;
    defval->vals.integer = 2;
      /*! \page aimInputsCART3D
     * - <B> SharpFeatureDivisions = int</B>. <Default 2> <br>
     * nominal # of ADDITIONAL divisions in backgrnd mesh around sharp features
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
  defval->vals.integer = 10;
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

  } else if (index == Xslices) {
    *name               = EG_strdup("Xslices");
    defval->type          = Double;
    defval->lfixed        = Change;
    defval->sfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;
      /*! \page aimInputsCART3D
     * - <B> Xslices = double </B> <br> OR
     * - <B> Xslices = [double, ... , double] </B> <br>
     * X slice locations created in output.
     */
  } else if (index == Yslices) {
    *name               = EG_strdup("Yslices");
    defval->type          = Double;
    defval->lfixed        = Change;
    defval->sfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;

     /*! \page aimInputsCART3D
    * - <B> Yslices = double </B> <br> OR
    * - <B> Yslices = [double, ... , double] </B> <br>
    * Y slice locations created in output.
      */
  } else if (index == Zslices) {
    *name               = EG_strdup("Zslices");
    defval->type          = Double;
    defval->lfixed        = Change;
    defval->sfixed        = Change;
    defval->nullVal       = IsNull;
    defval->dim           = Vector;

     /*! \page aimInputsCART3D
    * - <B> Zslices = double </B> <br> OR
    * - <B> Zslices = [double, ... , double] </B> <br>
    * Z slice locations created in output.
    */
  }
  else if (index == Model_X_axis) {
    *name = EG_strdup("Model_X_axis");
    defval->type = String;
    defval->units = NULL;
    defval->vals.string = EG_strdup("-Xb");
    
    /*! \page aimInputsCART3D
     * - <B> Model_X_axis = string </B> <br>
     * Model_X_axis defines x-axis orientation.
     */
  }
  else if (index == Model_Y_axis) {
    *name = EG_strdup("Model_Y_axis");
    defval->type = String;
    defval->units = NULL;
    defval->vals.string = EG_strdup("Yb");
    
    /*! \page aimInputsCART3D
     * - <B> Model_Y_axis = string </B> <br>
     * Model_Y_axis defines y-axis orientation.
     */
  }
  else if (index == Model_Z_axis) {
    *name = EG_strdup("Model_Z_axis");
    defval->type = String;
    defval->units = NULL;
    defval->vals.string = EG_strdup("-Zb");
    
    /*! \page aimInputsCART3D
     * - <B> Model_Z_axis = string </B> <br>
     * Model_Z_axis defines z-axis orientation.
     */
  }
  return status;
}


/* ********************** Exposed AIM Functions ***************************** */
int
aimPreAnalysis(void *instStore, void *aimInfo, capsValue *inputs)
{
  int          i, j, varid, status, nBody, nface, nedge, nvert, ntri, *tris;
  int          atype, alen;
  double       params[3], box[6], size, *xyzs, Sref, Cref, Xref, Yref, Zref;
  char         line[128], aimFile[PATH_MAX];
  const int    *ints;
  const char   *string, *intents;
  const double *reals;
  c3dAIM       *cartInstance;
  ego          *bodies, tess;
  verTags      *vtags;
  FILE         *fp;

#ifdef DEBUG
  printf(" Cart3DAIM/aimPreAnalysis!\n");
#endif
  
  if (inputs == NULL) {
#ifdef DEBUG
    printf(" Cart3DAIM/aimPreAnalysis -- NULL inputs!\n");
#endif
    return CAPS_NULLOBJ;
  }
  
  cartInstance = (c3dAIM *) instStore;

  // Get AIM bodies
  status = aim_getBodies(aimInfo, &intents, &nBody, &bodies);
  if (status != CAPS_SUCCESS) {
    return status;
  }

  if (nBody != 1) {
    printf(" Cart3DAIM/aimPreAnalysis nBody = %d!\n", nBody);
    return CAPS_SOURCEERR;
  }

  // Initialize reference values
  Sref = Cref = 1.0;
  Xref = Yref = Zref = 0.0;

  // Loop over bodies and look for reference quantity attributes
  for (i = 0; i < nBody; i++) {
    status = EG_attributeRet(bodies[i], "capsReferenceArea", &atype, &alen,
                             &ints, &reals, &string);
    if ((status == EGADS_SUCCESS) && (atype == ATTRREAL)) Sref = (double) reals[0];
    
    status = EG_attributeRet(bodies[i], "capsReferenceChord", &atype, &alen,
                             &ints, &reals, &string);
    if ((status == EGADS_SUCCESS) && (atype == ATTRREAL)) Cref = (double) reals[0];
    
    status = EG_attributeRet(bodies[i], "capsReferenceX", &atype, &alen,
                             &ints, &reals, &string);
    if ((status == EGADS_SUCCESS) && (atype == ATTRREAL)) Xref = (double) reals[0];
    
    status = EG_attributeRet(bodies[i], "capsReferenceY", &atype, &alen,
                             &ints, &reals, &string);
    if ((status == EGADS_SUCCESS) && (atype == ATTRREAL)) Yref = (double) reals[0];
    
    status = EG_attributeRet(bodies[i], "capsReferenceZ", &atype, &alen,
                             &ints, &reals, &string);
    if ((status == EGADS_SUCCESS) && (atype == ATTRREAL)) Zref = (double) reals[0];
  }

  // Only do the meshing if the geometry has been changed or a meshing input variable has been changed
  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS ||
      cartInstance->tessParam[0]          != inputs[Tess_Params-1].vals.reals[0] ||
      cartInstance->tessParam[1]          != inputs[Tess_Params-1].vals.reals[1] ||
      cartInstance->tessParam[2]          != inputs[Tess_Params-1].vals.reals[2] ||
      cartInstance->outerBox              != inputs[outer_box-1].vals.real ||
      cartInstance->nDiv                  != inputs[nDiv-1].vals.integer ||
      cartInstance->maxR                  != inputs[maxR-1].vals.integer ||
      cartInstance->sharpFeatureDivisions != inputs[SharpFeatureDivisions-1].vals.integer ||
      cartInstance->nMultiGridLevels      != inputs[nMultiGridLevels-1].vals.integer) {

    // Set new meshing parameters
    cartInstance->tessParam[0] = inputs[Tess_Params-1].vals.reals[0];
    cartInstance->tessParam[1] = inputs[Tess_Params-1].vals.reals[1];
    cartInstance->tessParam[2] = inputs[Tess_Params-1].vals.reals[2];

    cartInstance->outerBox = inputs[outer_box-1].vals.real;
    cartInstance->nDiv     = inputs[nDiv-1].vals.integer;
    cartInstance->maxR     = inputs[maxR-1].vals.integer;
    cartInstance->sharpFeatureDivisions = inputs[SharpFeatureDivisions-1].vals.integer;
    cartInstance->nMultiGridLevels      = inputs[nMultiGridLevels-1].vals.integer;

    status = EG_getBoundingBox(bodies[0], box);
    if (status != EGADS_SUCCESS) return status;

    size = sqrt((box[3]-box[0])*(box[3]-box[0]) +
          (box[4]-box[1])*(box[4]-box[1]) +
          (box[5]-box[2])*(box[5]-box[2]));
    printf(" Body size = %lf\n", size);

    params[0] = inputs[0].vals.reals[0]*size;
    params[1] = inputs[0].vals.reals[1]*size;
    params[2] = inputs[0].vals.reals[2];
    printf(" Tessellating with  MaxEdge = %lf  Sag = %lf  Angle = %lf\n",
           params[0], params[1], params[2]);

    status = EG_makeTessBody(bodies[0], params, &tess);
    if (status != EGADS_SUCCESS) {
      return status;
    }

    /* cleanup old data and get new tessellation */
    status = bodyTess(tess, &nface, &nedge, &nvert, &xyzs, &vtags, &ntri, &tris);
    if (status != EGADS_SUCCESS) {
      EG_deleteObject(tess);
      return status;
    }
    cartInstance->nface = nface;
    cartInstance->nvert = nvert;
    cartInstance->ntris = ntri;

    status = aim_file(aimInfo, "Components.i.tri", aimFile);
    if (status != EGADS_SUCCESS) {
      EG_deleteObject(tess);
      return status;
    }

    /* write the tri file */
    status = writeTrix(aimFile, NULL, 1, nvert, xyzs, 0, NULL,
             ntri, tris);
    EG_free(tris);
    EG_free(vtags);
    EG_free(xyzs);
    if (status != 0) {
      printf(" writeTrix return = %d\n", status);
      EG_deleteObject(tess);
      return CAPS_IOERR;
    }

    /* store away the tessellation */
    status = aim_newTess(aimInfo, tess);
    if (status != 0) {
      printf(" aim_setTess return = %d\n", status);
      EG_deleteObject(tess);
      return status;
    }

    snprintf(line, 128, "autoInputs -r %lf -nDiv %d -maxR %d",
         inputs[1].vals.real, inputs[2].vals.integer,
         inputs[3].vals.integer);
    printf(" Executing: %s\n", line);
    status = aim_system(aimInfo, NULL, line);
    if (status != 0) {
      printf(" ERROR: autoInputs return = %d\n", status);
      return CAPS_EXECERR;
    }

    snprintf(line, 128, "cubes -reorder -sf %d",
        inputs[9].vals.integer);
    printf(" Executing: %s\n", line);
    status = aim_system(aimInfo, NULL, line);
    if (status != 0) {
      printf(" ERROR: cubes return = %d\n", status);
      return CAPS_EXECERR;
    }

    snprintf(line, 128, "mgPrep -n %d",
        inputs[10].vals.integer);
    printf(" Executing: %s\n", line);
    status = aim_system(aimInfo, NULL, line);
    if (status != 0) {
      printf(" ERROR: mgPrep return = %d\n", status);
      return CAPS_EXECERR;
    }
  }

  /* create and output input.cntl */
  fp = aim_fopen(aimInfo, "input.cntl", "w");
  if (fp == NULL) {
#ifdef DEBUG
    printf(" Cart3DAIM/aimCalcOutput Cannot open input.cntl!\n");
#endif
    return CAPS_DIRERR;
  }

  fprintf(fp, "$__Case_Information:\n\n");
  fprintf(fp, "Mach     %lf\n", inputs[4].vals.real);
  fprintf(fp, "alpha    %lf\n", inputs[5].vals.real);
  fprintf(fp, "beta     %lf\n", inputs[6].vals.real);
  fprintf(fp, "gamma    %lf\n", inputs[7].vals.real);

  fprintf(fp, "\n$__File_Name_Information:\n\n");
  fprintf(fp, "MeshInfo           Mesh.c3d.Info\n");
  fprintf(fp, "MeshFile           Mesh.mg.c3d\n\n");

  fprintf(fp, "$__Solver_Control_Information:\n\n");
  fprintf(fp, "RK        0.0695     1\n");
  fprintf(fp, "RK        0.1602     0\n");
  fprintf(fp, "RK        0.2898     0\n");
  fprintf(fp, "RK        0.5060     0\n");
  fprintf(fp, "RK        1.0        0\n\n");
  fprintf(fp, "CFL           %lf\n", inputs[14].vals.real);
  fprintf(fp, "Limiter       %d\n",  inputs[15].vals.integer);
  fprintf(fp, "FluxFun       %d\n",  inputs[16].vals.integer);
  fprintf(fp, "maxCycles     %d\n",  inputs[8].vals.integer);
  fprintf(fp, "Precon        0\n");
  fprintf(fp, "wallBCtype    0\n");
  fprintf(fp, "nMGlev        %d\n",   inputs[10].vals.integer);
  fprintf(fp, "MG_cycleType  %d\n",   inputs[11].vals.integer);
  fprintf(fp, "MG_nPre       %d\n",   inputs[12].vals.integer);
  fprintf(fp, "MG_nPost      %d\n\n", inputs[13].vals.integer);

  fprintf(fp, "$__Boundary_Conditions:\n\n");
/*                   # BC types: 0 = FAR FIELD
                                 1 = SYMMETRY
                                 2 = INFLOW  (specify all)
                                 3 = OUTFLOW (simple extrap)  */
  fprintf(fp, "Dir_Lo_Hi     0   0 0\n"); /* (0/1/2) direction - Low BC - Hi BC */
  fprintf(fp, "Dir_Lo_Hi     1   0 0\n");
  fprintf(fp, "Dir_Lo_Hi     2   0 0\n\n");

  fprintf(fp, "$__Convergence_History_reporting:\n\n");
  fprintf(fp, "iForce     %d\n", inputs[17].vals.integer);
  fprintf(fp, "iHist      %d\n", inputs[18].vals.integer);
  fprintf(fp, "nOrders    %d\n", inputs[19].vals.integer);
  fprintf(fp, "refArea    %lf\n", Sref);
  fprintf(fp, "refLength  %lf\n", Cref);

  fprintf(fp, "\n$__Partition_Information:\n\n");
  fprintf(fp, "nPart      1\n");
  fprintf(fp, "type       1\n");

  fprintf(fp, "\n$__Post_Processing:\n\n");
  varid = 20;
  if (inputs[varid].nullVal != IsNull) {
    if (inputs[varid].length == 1) {
      fprintf(fp,"Xslices %lf\n",inputs[varid].vals.real);
    } else {
      fprintf(fp,"Xslices");
      for (j = 0; j < inputs[varid].length; j++) {
        fprintf(fp," %lf",inputs[varid].vals.reals[j]);
      }
      fprintf(fp,"\n");
    }
  }

  varid = 21;
  if (inputs[varid].nullVal != IsNull) {
    if (inputs[varid].length == 1) {
      fprintf(fp,"Yslices %lf\n",inputs[varid].vals.real);
    } else {
      fprintf(fp,"Yslices");
      for (j = 0; j < inputs[varid].length; j++) {
        fprintf(fp," %lf",inputs[varid].vals.reals[j]);
      }
      fprintf(fp,"\n");
    }
  }

  varid = 22;
    if (inputs[varid].nullVal != IsNull) {
      if (inputs[varid].length == 1) {
        fprintf(fp,"Zslices %lf\n",inputs[varid].vals.real);
      } else {
        fprintf(fp,"Zslices");
        for (j = 0; j < inputs[varid].length; j++) {
          fprintf(fp," %lf",inputs[varid].vals.reals[j]);
        }
        fprintf(fp,"\n");
      }
    }

  fprintf(fp, "\n$__Force_Moment_Processing:\n\n");
/* ... Axis definitions (with respect to body axis directions (Xb,Yb,Zb)
                         w/ usual stability and control orientation)  */
  fprintf(fp, "Model_X_axis  %s\n", inputs[23].vals.string);
  fprintf(fp, "Model_Y_axis  %s\n", inputs[24].vals.string);
  fprintf(fp, "Model_Z_axis  %s\n", inputs[25].vals.string);
  fprintf(fp, "Reference_Area   %lf all\n", Sref);
  fprintf(fp, "Reference_Length %lf all\n", Cref);
  fprintf(fp, "Force entire\n\n");
  fprintf(fp, "Moment_Point %lf %lf %lf entire\n",Xref ,Yref ,Zref );

  fclose(fp);

  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************

/* no longer optional and needed for restart */
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int
aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc, int index,
           char **aoname, capsValue *form)
{
  char *names[12] = {"C_A", "C_Y", "C_N",   "C_D",   "C_S", "C_L", "C_l",
                     "C_m", "C_n", "C_M_x", "C_M_y", "C_M_z"};
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
  */

  *aoname    = EG_strdup(names[index-1]);
  form->type = Double;

  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
int
aimCalcOutput(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, capsValue *val)
{
  int        offset;
  size_t     linecap = 0;
  const char comp = ':';
  char       *valstr, *line = NULL;
  FILE       *fp;
  char       *start[12]= { "entire   Axial Force (C_A):",
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
                           "entire   Z Aero Moment" };

  if (index < 7) {
    offset = 27;
  } else {
    offset = 1;
  }

  /* open the Cart3D loads file */
  fp = aim_fopen(aimInfo, "loadsCC.dat", "r");
  if (fp == NULL) {
#ifdef DEBUG
    printf(" Cart3DAIM/aimCalcOutput Cannot open Output file!\n");
#endif
    return CAPS_DIRERR;
  }

  /* scan the file for the string */
  valstr = NULL;
  while (getline(&line, &linecap, fp) >= 0) {
    if (line == NULL) continue;
    valstr = strstr(line, start[index-1]);
    if (valstr != NULL) {
      if (index > 6) {
        valstr = strchr(line,comp);
      }
      
#ifdef DEBUG
      printf("valstr > |%s|\n",valstr);
#endif
      /* found it -- get the value */
      sscanf(&valstr[offset], "%lf", &val->vals.real);
      break;
    }
  }
  fclose(fp);
  if (line != NULL) EG_free(line);

  if (valstr == NULL) {
#ifdef DEBUG
    printf(" Cart3DAIM/aimCalcOutput Cannot find %s in Output file!\n",
           start[index-1]);
#endif
    return CAPS_NOTFOUND;
  }

  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
void
aimCleanup(/*@null@*/ void *instStore)
{
#ifdef DEBUG
  printf(" Cart3DAIM/aimCleanup!\n");
#endif

  EG_free(instStore);
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

#ifdef DEBUG
  printf(" capsAIM/aimDiscr: tname = %s, instance = %d!\n",
         tname, aim_getInstance(discr->aInfo));
#endif

  /* create the discretization structure for this capsBound */

  status = aim_getBodies(discr->aInfo, &intents, &nBody, &bodies);
  if (status != CAPS_SUCCESS) goto cleanup;

  /* bodies is 2*nBodies long where the last nBodies
     are the tess objects for the first nBodies */
  tess = bodies + nBody;

  /* find any faces with our boundary marker */
  for (nBodyDisc = ibody = 0; ibody < nBody; ibody++) {
    if (tess[ibody] == NULL) continue;
    status = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nFace, &faces);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (faces == NULL) {
      status = EGADS_TOPOERR;
      goto cleanup;
    }
    found = 0;
    for (iface = 0; iface < nFace; iface++) {
      status = EG_attributeRet(faces[iface], "capsBound", &atype, &alen,
                               &ints, &reals, &string);
      if (status != EGADS_SUCCESS)    continue;
      if (atype  != ATTRSTRING)       continue;
      if (strcmp(string, tname) != 0) continue;
#ifdef DEBUG
      printf(" skeletonAIM/aimDiscr: Body %d/Face %d matches %s!\n",
             ibody+1, iface+1, tname);
#endif
      /* count the number of face with capsBound */
      found = 1;
    }
    if (found == 1) nBodyDisc++;
    AIM_FREE(faces);
  }
  if (nBodyDisc == 0) {
    printf(" skeletonAIM/aimDiscr: No Faces match %s!\n", tname);
    return CAPS_SUCCESS;
  }

  /* specify our single triangle element type */
  discr->nTypes = 1;
  AIM_ALLOC(discr->types, discr->nTypes, capsEleType, discr->aInfo, status);

  discr->types[0].nref  = 3;
  discr->types[0].ndata = 0;         /* data at geom reference positions
                                        (i.e. vertex centered/iso-parametric) */
  discr->types[0].ntri  = 1;
  discr->types[0].nmat  = 0;         /* match points at geom ref positions */
  discr->types[0].tris  = NULL;
  discr->types[0].gst   = NULL;
  discr->types[0].dst   = NULL;
  discr->types[0].matst = NULL;

  /* specify the numbering for the points on the triangle */
  AIM_ALLOC(discr->types[0].tris, discr->types[0].nref, int, discr->aInfo, status);

  discr->types[0].tris[0] = 1;
  discr->types[0].tris[1] = 2;
  discr->types[0].tris[2] = 3;

  /* specify the reference coordinates for each point on the triangle */
  AIM_ALLOC(discr->types[0].gst, 2*discr->types[0].nref, double, discr->aInfo, status);

  discr->types[0].gst[0] = 0.0;   /* s = 0, t = 0 */
  discr->types[0].gst[1] = 0.0;
  discr->types[0].gst[2] = 1.0;   /* s = 1, t = 0 */
  discr->types[0].gst[3] = 0.0;
  discr->types[0].gst[4] = 0.0;   /* s = 0, t = 1 */
  discr->types[0].gst[5] = 1.0;

  /* allocate the body discretizations */
  AIM_ALLOC(discr->bodys, discr->nBodys, capsBodyDiscr, discr->aInfo, status);

  /* get the tessellation and
   make up a linear continuous triangle discretization */
  vID = nBodyDisc = 0;
  for (ibody = 0; ibody < nBody; ibody++) {
    if (tess[ibody] == NULL) continue;
    ntris = 0;
    AIM_FREE(faces);
    status = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nFace, &faces);
    if ((status != EGADS_SUCCESS) || (faces == NULL)) {
      printf(" skeletonAIM: getBodyTopos (Face) = %d for Body %d!\n",
             status, ibody+1);
      status = EGADS_TOPOERR;
      goto cleanup;
    }
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
        printf(" skeletonAIM: EG_getTessFace %d = %d for Body %d!\n",
               iface+1, status, ibody+1);
        continue;
      }
      ntris += tlen;
      found = 1;
    }
    if (found == 0) continue;
    if (ntris == 0) {
#ifdef DEBUG
      printf(" skeletonAIM/aimDiscr: ntris = %d!\n", ntris);
#endif
      status = CAPS_SOURCEERR;
      goto cleanup;
    }

    discBody = &discr->bodys[nBodyDisc];
    aim_initBodyDiscr(discBody);

    discBody->tess = tess[ibody-1];
    discBody->nElems = ntris;

    AIM_ALLOC(discBody->elems   ,   ntris, capsElement, discr->aInfo, status);
    AIM_ALLOC(discBody->gIndices, 6*ntris, int        , discr->aInfo, status);

    status = EG_statusTessBody(tess[ibody-1], &body, &state, &nGlobal);
    AIM_STATUS(discr->aInfo, status);

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

      status = EG_getTessFace(tess[ibody-1], iface, &alen, &xyz, &uv,
                            &ptype, &pindex, &tlen, &tris, &nei);
      AIM_STATUS(discr->aInfo, status);

      /* construct global vertex indices */
      for (i = 0; i < alen; i++) {
        status = EG_localToGlobal(tess[ibody-1], iface, i+1, &global);
        AIM_STATUS(discr->aInfo, status);
        if (vid[global-1] != 0) continue;
        vid[global-1] = vID+1;
        vID++;
      }

      /* fill elements */
      for (itri = 0; itri < tlen; itri++, ielem++) {
        discBody->elems[ielem].tIndex      = 1;
        discBody->elems[ielem].eIndex      = iface;
/*@-immediatetrans@*/
        discBody->elems[ielem].gIndices    = &discBody->gIndices[6*ielem];
/*@+immediatetrans@*/
        discBody->elems[ielem].dIndices    = NULL;
        discBody->elems[ielem].eTris.tq[0] = itri+1;

        status = EG_localToGlobal(tess[ibody-1], iface, tris[3*itri  ],
                                  &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[0] = vid[global-1];
        discBody->elems[ielem].gIndices[1] = tris[3*itri  ];

        status = EG_localToGlobal(tess[ibody-1], iface, tris[3*itri+1],
                                  &global);
        AIM_STATUS(discr->aInfo, status);
        discBody->elems[ielem].gIndices[2] = vid[global-1];
        discBody->elems[ielem].gIndices[3] = tris[3*itri+1];

        status = EG_localToGlobal(tess[ibody-1], iface, tris[3*itri+2],
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
  int    i, j, global, stat, bIndex;
  double *rvec;
  c3dAIM *cartInstance;

#ifdef DEBUG
  printf(" Cart3DAIM/aimTransfer name = %s  npts = %d/%d!\n",
         name, npts, rank);
#endif
  cartInstance = (c3dAIM *) discr->instStore;

  rvec = (double *) EG_alloc(rank*cartInstance->nvert*sizeof(double));
  if (rvec == NULL) {
#ifdef DEBUG
    printf(" Cart3DAIM/aimTransfer Cannot allocate %dx%d!\n",
           rank, cartInstance->nvert);
#endif
    return EGADS_MALLOC;
  }

  /* try and read the trix file */
  stat = readTrix("Components.i.trix", name, cartInstance->nvert, rank, rvec);
  if (stat != 0) {
    EG_free(rvec);
    return CAPS_IOERR;
  }

  /* move the appropriate parts of the tessellation to data */
  for (i = 0; i < npts; i++) {
    /* points might span multiple bodies */
    bIndex = discr->tessGlobal[2*i  ];
    global = discr->tessGlobal[2*i+1] +
             discr->bodys[bIndex-1].globalOffset;
    for (j = 0; j < rank; j++) data[rank*i+j] = rvec[rank*global+j];
  }

  EG_free(rvec);
  return CAPS_SUCCESS;
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
