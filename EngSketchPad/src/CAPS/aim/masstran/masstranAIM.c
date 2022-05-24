/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             masstran AIM
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


/*! \mainpage Introduction
 * \tableofcontents
 *
 * \section overviewMasstran Masstran AIM Overview
 * A module in the Computational Aircraft Prototype Syntheses (CAPS) has been developed to compute
 * mass properties using attributions for finite element structural solvers.
 *
 * An outline of the AIM's inputs, outputs and attributes are provided in \ref aimInputsMasstran and
 * \ref aimOutputsMasstran and \ref attributeMasstran, respectively.
 *
 * Details on the use of units are outlined in \ref aimUnitsMasstran.
 *
 * The mass properties are computed via the formulas:
 *
 *  \f{eqnarray*}{
 *  m &=& \sum_i m_i \\
 *  x_{cg} &=& \frac{1}{m}\sum_i m_i x_i \\
 *  y_{cg} &=& \frac{1}{m}\sum_i m_i y_i \\
 *  z_{cg} &=& \frac{1}{m}\sum_i m_i z_i \\
 * \big(I_{xx}\big)_{cg} &=& \sum_i m_i \big(y_i^2 + z_i^2\big) - m \big(y_{cg}^2 + z_{cg}^2\big) \\
 * \big(I_{yy}\big)_{cg} &=& \sum_i m_i \big(x_i^2 + z_i^2\big) - m \big(x_{cg}^2 + z_{cg}^2\big) \\
 * \big(I_{zz}\big)_{cg} &=& \sum_i m_i \big(x_i^2 + y_i^2\big) - m \big(x_{cg}^2 + y_{cg}^2\big) \\
 * \big(I_{xy}\big)_{cg} &=& \sum_i m_i \big(x_i     y_i  \big) - m \big(x_{cg}     y_{cg}  \big) \\
 * \big(I_{xz}\big)_{cg} &=& \sum_i m_i \big(x_i     z_i  \big) - m \big(x_{cg}     z_{cg}  \big) \\
 * \big(I_{yz}\big)_{cg} &=& \sum_i m_i \big(y_i     z_i  \big) - m \big(y_{cg}     z_{cg}  \big),
 * \f}
 *
 * where i represents an element index in the mesh, and the mass \f$m_i\f$ is computed from the density, thickness, and area of the element.
 *
 * The moment of inertias are accessible individually, in vector form as
 *
 * \f[
 *  \vec{I} = \begin{bmatrix} I_{xx} & I_{yy} & I_{zz} & I_{xy} & I_{xz} & I_{yz} \end{bmatrix},
 * \f]
 *
 * as lower/upper triangular form
 *
 * \f[
 *  \vec{I}_{lower} = \begin{bmatrix} I_{xx} & -I_{xy} & I_{yy} & -I_{xz} & -I_{yz} & I_{zz} \end{bmatrix},
 * \f]
 *
 * \f[
 *  \vec{I}_{upper} = \begin{bmatrix} I_{xx} & -I_{xy} & -I_{xz} & I_{yy} & -I_{yz} & I_{zz} \end{bmatrix},
 * \f]
 *
 * or in full tensor form as
 *
 * \f[
 * \bar{\bar{I}} =
 *  \begin{bmatrix}
 *    I_{xx} & -I_{xy} & -I_{xz} \\
 *   -I_{xy} &  I_{yy} & -I_{yz} \\
 *   -I_{xz} & -I_{yz} &  I_{zz}
 *  \end{bmatrix}.
 * \f]
 *
 *\section masstranExamples Examples
 * An example problem using the Masstran AIM may be found at \ref masstranExample.
 */


/*! \page attributeMasstran Masstran AIM attributes
 * The following list of attributes are required for the MYSTRAN AIM inside the geometry input.
 *
 * - <b> capsAIM</b> This attribute is a CAPS requirement to indicate the analysis the geometry
 * representation supports.
 *
 * - <b> capsGroup</b> This is a name assigned to any geometric body.  This body could be a solid, surface, face, wire, edge or node.
 * Recall that a string in ESP starts with a $.  For example, attribute <c>capsGroup $Wing</c>.
 *
 * - <b> capsIgnore</b> It is possible that there is a geometric body (or entity) that you do not want the Masstran AIM to pay attention to when creating
 * a finite element model. The capsIgnore attribute allows a body (or entity) to be in the geometry and ignored by the AIM.  For example,
 * because of limitations in OpenCASCADE a situation where two edges are overlapping may occur; capsIgnore allows the user to only
 *  pay attention to one of the overlapping edges.
 *
 */

#include <string.h>
#include <math.h>
#include "aimUtil.h"

#include "meshUtils.h"
#include "miscUtils.h"
#include "feaUtils.h"

//#define DEBUG

#define NINT(A)   (((A) < 0)   ? (int)(A-0.5) : (int)(A+0.5))
#define MIN(A,B)  (((A) < (B)) ? (A) : (B))
#define MAX(A,B)  (((A) < (B)) ? (B) : (A))

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


enum aimInputs
{
  Tess_Params = 1,               /* index is 1-based */
  Edge_Point_Min,
  Edge_Point_Max,
  Quad_Mesh,
  Property,
  Material,
  Surface_Mesh,
  NUMINPUT = Surface_Mesh        /* Total number of inputs */
};

#define NUMOUTPUT  15


typedef struct {

  double area;
  double mass;
  double Ixx, Iyy, Izz;
  double Ixy, Ixz, Iyz;
  double Cx, Cy, Cz;
  double CGx, CGy, CGz;

} massProperties;

/* AIM "local" per instance storage
   needed data should be added here & cleaned up in aimCleanup */
typedef struct {

  feaUnitsStruct units; // units system

  feaProblemStruct feaProblem;

  // Attribute to index map
  mapAttrToIndexStruct attrMap;

  // Mesh holders
  int numMesh;
  meshStruct *feaMesh;

  // mass properties
  massProperties massProp;

} aimStorage;



static int initiate_aimStorage(aimStorage *masstranInstance)
{

    int status;

    // Mesh holders
    masstranInstance->numMesh = 0;
    masstranInstance->feaMesh = NULL;

    status = initiate_feaUnitsStruct(&masstranInstance->units);
    if (status != CAPS_SUCCESS) return status;

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&masstranInstance->attrMap);
    if (status != CAPS_SUCCESS) return status;

    status = initiate_feaProblemStruct(&masstranInstance->feaProblem);
    if (status != CAPS_SUCCESS) return status;

    return CAPS_SUCCESS;
}


static int destroy_aimStorage(aimStorage *masstranInstance)
{

    int status;
    int i;
  
    status = destroy_feaUnitsStruct(&masstranInstance->units);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaUnitsStruct!\n", status);

    // Attribute to index map
    status = destroy_mapAttrToIndexStruct(&masstranInstance->attrMap);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_mapAttrToIndexStruct!\n", status);

    // Cleanup meshes
    if (masstranInstance->feaMesh != NULL) {

        for (i = 0; i < masstranInstance->numMesh; i++) {
            status = destroy_meshStruct(&masstranInstance->feaMesh[i]);
            if (status != CAPS_SUCCESS)
              printf("Error: Status %d during destroy_meshStruct!\n", status);
        }

        AIM_FREE(masstranInstance->feaMesh);
    }

    masstranInstance->feaMesh = NULL;
    masstranInstance->numMesh = 0;

    // Destroy FEA problem structure
    status = destroy_feaProblemStruct(&masstranInstance->feaProblem);
    if (status != CAPS_SUCCESS)
      printf("Error: Status %d during destroy_feaProblemStruct!\n", status);

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(void *aimInfo, aimStorage *masstranInstance)
{
  // Function return flag
  int status = CAPS_SUCCESS;
  int i, remesh = (int)true;

  // Meshing related variables
  double tessParam[3] = {0.025, 0.001, 15};
  int edgePointMin = 2;
  int edgePointMax = 50;
  int quadMesh = (int) false;

  // Dummy attribute to maps
  mapAttrToIndexStruct constraintMap;
  mapAttrToIndexStruct loadMap;
  mapAttrToIndexStruct transferMap;
  mapAttrToIndexStruct connectMap;

  // analysis input values
  capsValue *TessParams = NULL;
  capsValue *EdgePoint_Min = NULL;
  capsValue *EdgePoint_Max = NULL;
  capsValue *QuadMesh = NULL;

  for (i = 0; i < masstranInstance->numMesh; i++) {
      remesh = remesh && (masstranInstance->feaMesh[i].egadsTess->oclass == EMPTY);
  }
  if (remesh == (int) false) return CAPS_SUCCESS;

  // retrieve or create the mesh from fea_createMesh
  status = aim_getValue(aimInfo, Tess_Params, ANALYSISIN, &TessParams);
  AIM_STATUS(aimInfo, status);

  status = aim_getValue(aimInfo, Edge_Point_Min, ANALYSISIN, &EdgePoint_Min);
  AIM_STATUS(aimInfo, status);

  status = aim_getValue(aimInfo, Edge_Point_Max, ANALYSISIN, &EdgePoint_Max);
  AIM_STATUS(aimInfo, status);

  status = aim_getValue(aimInfo, Quad_Mesh, ANALYSISIN, &QuadMesh);
  AIM_STATUS(aimInfo, status);

  if (TessParams != NULL) {
      tessParam[0] = TessParams->vals.reals[0];
      tessParam[1] = TessParams->vals.reals[1];
      tessParam[2] = TessParams->vals.reals[2];
  }

  // Max and min number of points
  if (EdgePoint_Min != NULL && EdgePoint_Min->nullVal != IsNull) {
      edgePointMin = EdgePoint_Min->vals.integer;
      if (edgePointMin < 2) {
        AIM_ANALYSISIN_ERROR(aimInfo, Edge_Point_Min, "Edge_Point_Min = %d must be greater or equal to 2\n", edgePointMin);
        return CAPS_BADVALUE;
      }
  }

  if (EdgePoint_Max != NULL && EdgePoint_Max->nullVal != IsNull) {
      edgePointMax = EdgePoint_Max->vals.integer;
      if (edgePointMax < 2) {
        AIM_ANALYSISIN_ERROR(aimInfo, Edge_Point_Max, "Edge_Point_Max = %d must be greater or equal to 2\n", edgePointMax);
        return CAPS_BADVALUE;
      }
  }

  if (edgePointMin >= 2 && edgePointMax >= 2 && edgePointMin > edgePointMax) {
    AIM_ERROR  (aimInfo, "Edge_Point_Max must be greater or equal Edge_Point_Min");
    AIM_ADDLINE(aimInfo, "Edge_Point_Max = %d, Edge_Point_Min = %d\n",edgePointMax,edgePointMin);
    return CAPS_BADVALUE;
  }

  if (QuadMesh != NULL) quadMesh = QuadMesh->vals.integer;
  
  status = initiate_mapAttrToIndexStruct(&constraintMap);
  AIM_STATUS(aimInfo, status);

  status = initiate_mapAttrToIndexStruct(&loadMap);
  AIM_STATUS(aimInfo, status);

  status = initiate_mapAttrToIndexStruct(&transferMap);
  AIM_STATUS(aimInfo, status);

  status = initiate_mapAttrToIndexStruct(&connectMap);
  AIM_STATUS(aimInfo, status);

/*@-nullpass@*/
  status = fea_createMesh(aimInfo,
                          tessParam,
                          edgePointMin,
                          edgePointMax,
                          quadMesh,
                          &masstranInstance->attrMap,
                          &constraintMap,
                          &loadMap,
                          &transferMap,
                          &connectMap,
                          NULL,
                          &masstranInstance->numMesh,
                          &masstranInstance->feaMesh,
                          &masstranInstance->feaProblem);
/*@-nullpass@*/
  AIM_STATUS(aimInfo, status);

  status = destroy_mapAttrToIndexStruct(&constraintMap);
  AIM_STATUS(aimInfo, status);

  status = destroy_mapAttrToIndexStruct(&loadMap);
  AIM_STATUS(aimInfo, status);

  status = destroy_mapAttrToIndexStruct(&transferMap);
  AIM_STATUS(aimInfo, status);

  status = destroy_mapAttrToIndexStruct(&connectMap);
  AIM_STATUS(aimInfo, status);

cleanup:
  return status;
}


/****************** exposed AIM entry points -- Analysis **********************/

/* aimInitialize: Initialization Information for the AIM */
int aimInitialize(int inst, const char *unitSys, void *aimInfo,
                  void **instStore, /*@unused@*/ int *major,
                  /*@unused@*/ int *minor, int *nIn, int *nOut,
                  int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int status = CAPS_SUCCESS;
  aimStorage *masstranInstance;
  const char *keyWord;
  char *keyValue = NULL, *tmpUnits = NULL;
  double real = 1.0;
  feaUnitsStruct *units=NULL;

#ifdef DEBUG
  printf("\n masstranAIM/aimInitialize  Instance!\n", inst);
#endif

  /* specify the number of analysis inputs  defined in aimInputs
   *     and the number of analysis outputs defined in aimOutputs */
  *nIn    = NUMINPUT;
  *nOut   = NUMOUTPUT;

  /* return if "query" only */
  if (inst == -1) return CAPS_SUCCESS;

  /* specify the field variables this analysis can generate and consume */
  *nFields = 0;
  *fnames  = NULL;
  *franks  = NULL;
  *fInOut  = NULL;

  /* create our "local" storage for anything that needs to be persistent */
  masstranInstance = (aimStorage *) EG_alloc((inst+1)*sizeof(aimStorage));
  if (masstranInstance == NULL) return EGADS_MALLOC;

  initiate_aimStorage(masstranInstance);
  *instStore = masstranInstance;

  /*! \page aimUnitsMasstran AIM Units
   *  A unit system may be optionally specified during AIM instance initiation. If
   *  a unit system is provided, all AIM  input values which have associated units must be specified as well.
   *  If no unit system is used, AIM inputs, which otherwise would require units, will be assumed
   *  unit consistent. A unit system may be specified via a JSON string dictionary for example:
   *  unitSys = "{"mass": "kg", "length": "m"}"
   */
  if (unitSys != NULL) {
    units = &masstranInstance->units;

    // Do we have a json string?
    if (strncmp( unitSys, "{", 1) != 0) {
      AIM_ERROR(aimInfo, "unitSys ('%s') is expected to be a JSON string dictionary", unitSys);
      return CAPS_BADVALUE;
    }

    /*! \page aimUnitsMasstran
     *  \section jsonStringMasstran JSON String Dictionary
     *  The key arguments of the dictionary are described in the following:
     *
     *  <ul>
     *  <li> <B>mass = "None"</B> </li> <br>
     *  Mass units - e.g. "kilogram", "k", "slug", ...
     *  </ul>
     */
    keyWord = "mass";
    status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
    if (status == CAPS_SUCCESS) {
      units->mass = string_removeQuotation(keyValue);
      AIM_FREE(keyValue);
      real = 1;
      status = aim_convert(aimInfo, 1, units->mass, &real, "kg", &real);
      AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->mass, keyWord);
    } else {
      AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    /*! \page aimUnitsMasstran
     *  <ul>
     *  <li> <B>length = "None"</B> </li> <br>
     *  Length units - e.g. "meter", "m", "inch", "in", "mile", ...
     *  </ul>
     */
    keyWord = "length";
    status  = search_jsonDictionary(unitSys, keyWord, &keyValue);
    if (status == CAPS_SUCCESS) {
      units->length = string_removeQuotation(keyValue);
      AIM_FREE(keyValue);
      real = 1;
      status = aim_convert(aimInfo, 1, units->length, &real, "m", &real);
      AIM_STATUS(aimInfo, status, "unitSys ('%s'): %s is not a %s unit", unitSys, units->length, keyWord);
    } else {
      AIM_ERROR(aimInfo, "unitSys ('%s') does not contain '%s'", unitSys, keyWord);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // construct density volume unit
    status = aim_unitRaise(aimInfo, units->length, -3, &tmpUnits); // 1/length^3
    AIM_STATUS(aimInfo, status);
    status = aim_unitMultiply(aimInfo, units->mass, tmpUnits, &units->densityVol); // mass/length^3, e.g volume density
    AIM_STATUS(aimInfo, status);
    AIM_FREE(tmpUnits);

    // construct density area unit
    status = aim_unitRaise(aimInfo, units->length, -2, &tmpUnits); // 1/length^2
    AIM_STATUS(aimInfo, status);
    status = aim_unitMultiply(aimInfo, units->mass, tmpUnits, &units->densityArea); // mass/length^2, e.g area density
    AIM_STATUS(aimInfo, status);
    AIM_FREE(tmpUnits);

    status = aim_unitRaise(aimInfo, units->length, 2, &tmpUnits ); // length^2
    AIM_STATUS(aimInfo, status);
    status = aim_unitMultiply(aimInfo, units->mass, tmpUnits, &units->momentOfInertia ); // mass*length^2, e.g moment of inertia
    AIM_STATUS(aimInfo, status);
    AIM_FREE(tmpUnits);
  }

cleanup:
  AIM_FREE(keyValue);
  AIM_FREE(tmpUnits);

  return status;
}


/* aimInputs: Input Information for the AIM */
int aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
              int index, char **ainame, capsValue *defval)
{

  /*! \page aimInputsMasstran AIM Inputs
   * The following list outlines the Masstran inputs along with their default value available
   * through the AIM interface.
   */
  int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" nastranAIM/aimInputs  index = %d!\n", index);
#endif

  *ainame = NULL;

  // Masstran Inputs
  if (index == Tess_Params) {
    *ainame               = EG_strdup("Tess_Params");
    defval->type          = Double;
    defval->dim           = Vector;
    defval->nrow          = 3;
    defval->ncol          = 1;
    defval->units         = NULL;
    defval->lfixed        = Fixed;
    defval->vals.reals    = (double *) EG_alloc(defval->nrow*sizeof(double));
    if (defval->vals.reals != NULL) {
      defval->vals.reals[0] = 0.025;
      defval->vals.reals[1] = 0.001;
      defval->vals.reals[2] = 15.00;
    } else return EGADS_MALLOC;

    /*! \page aimInputsMasstran
     * - <B> Tess_Params = [0.025, 0.001, 15.0]</B> <br>
     * Body tessellation parameters used when creating a boundary element model.
     * Tess_Params[0] and Tess_Params[1] get scaled by the bounding
     * box of the body. (From the EGADS manual) A set of 3 parameters that drive the EDGE discretization
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
    defval->vals.integer  = 2;
    defval->lfixed        = Fixed;
    defval->nrow          = 1;
    defval->ncol          = 1;
    defval->nullVal       = NotNull;

    /*! \page aimInputsMasstran
     * - <B> Edge_Point_Min = 2</B> <br>
     * Minimum number of points on an edge including end points to use when creating a surface mesh (min 2).
     */

  } else if (index == Edge_Point_Max) {
    *ainame               = EG_strdup("Edge_Point_Max");
    defval->type          = Integer;
    defval->vals.integer  = 50;
    defval->lfixed        = Fixed;
    defval->nrow          = 1;
    defval->ncol          = 1;
    defval->nullVal       = NotNull;

    /*! \page aimInputsMasstran
     * - <B> Edge_Point_Max = 50</B> <br>
     * Maximum number of points on an edge including end points to use when creating a surface mesh (min 2).
     */

  } else if (index == Quad_Mesh) {
    *ainame               = EG_strdup("Quad_Mesh");
    defval->type          = Boolean;
    defval->vals.integer  = (int) false;

    /*! \page aimInputsMasstran
     * - <B> Quad_Mesh = False</B> <br>
     * Create a quadratic mesh on four edge faces when creating the boundary element model.
     */

  } else if (index == Property) {
    *ainame              = EG_strdup("Property");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsMasstran
     * - <B> Property = NULL</B> <br>
     * Property tuple used to input property information for the model, see \ref feaProperty for additional details.
     */
  } else if (index == Material) {
    *ainame              = EG_strdup("Material");
    defval->type         = Tuple;
    defval->nullVal      = IsNull;
    defval->lfixed       = Change;
    defval->vals.tuple   = NULL;
    defval->dim          = Vector;

    /*! \page aimInputsMasstran
     * - <B> Material = NULL</B> <br>
     * Material tuple used to input material information for the model, see \ref feaMaterial for additional details.
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

      /*! \page aimInputsMasstran
       * - <B>Surface_Mesh = NULL</B> <br>
       * A Surface_Mesh link.
       */
  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
  return status;
}

/* aimPreAnalysis: Parse Inputs, Generate Input File(s) & Possibly Tessellate */
int
aimPreAnalysis(void *instStore, void *aimInfo, capsValue *aimInputs)
{
  int i, j;   // Indexing
  int status; // Status return
  int found;

  double area   = 0;
  double cgxmom = 0;
  double cgymom = 0;
  double cgzmom = 0;
  double cxmom  = 0;
  double cymom  = 0;
  double czmom  = 0;
  double mass   = 0;
  double Ixx    = 0;
  double Ixy    = 0;
  double Izz    = 0;
  double Ixz    = 0;
  double Iyy    = 0;
  double Iyz    = 0;
  double Cx, Cy, Cz, CGx, CGy, CGz;
  double xelem[4], yelem[4], zelem[4], elemArea, thick, density, elemWeight;
  double xcent, ycent, zcent;
  double dx1[3], dx2[3], n[3];
  double Lscale = 1;

  int elem[4];
  int *n2a = NULL;

  int numBody = 0;
  ego *bodies=NULL;

  const char *bodyLunits=NULL;
  const char *intents=NULL;

  feaMeshDataStruct *feaData = NULL;
  int propertyID, materialID;
  //meshElementSubTypeEnum elementSubType;
  feaPropertyStruct *feaProperty;
  feaMaterialStruct *feaMaterial;
  feaUnitsStruct *units=NULL;

  meshStruct     *nasMesh;
  massProperties *massProp;
  aimStorage     *masstranInstance;

  masstranInstance = (aimStorage *) instStore;

  if (aimInputs == NULL) return CAPS_NULLVALUE;

  // Get FEA mesh if we don't already have one
  if (aim_newGeometry(aimInfo) == CAPS_SUCCESS) {
    status = checkAndCreateMesh(aimInfo, masstranInstance);
    if (status != CAPS_SUCCESS) goto cleanup;
  }

  if (masstranInstance->units.length != NULL) {
    units = &masstranInstance->units;

    status = aim_getBodies(aimInfo, &intents, &numBody, &bodies);
    AIM_STATUS(aimInfo, status);

    if (numBody == 0 || bodies == NULL) {
      AIM_ERROR(aimInfo, "No Bodies!");
      status = CAPS_SOURCEERR;
      goto cleanup;
    }

    // Get length units
    status = check_CAPSLength(numBody, bodies, &bodyLunits);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "capsLength is not set in *.csm file!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // conversion of the csm model units into units of Lunits
    Lscale = 1.0;
    status = aim_convert(aimInfo, 1, bodyLunits, &Lscale, units->length, &Lscale);
    AIM_STATUS(aimInfo, status);
  }

  // Note: Setting order is important here.
  // 1. Materials should be set before properties.
  // 2. Coordinate system should be set before mesh and loads
  // 3. Mesh should be set before loads, constraints, supports, and connections

  // Set material properties
  if (aimInputs[Material-1].nullVal == NotNull) {
    status = fea_getMaterial(aimInfo,
                             aimInputs[Material-1].length,
                             aimInputs[Material-1].vals.tuple,
                             &masstranInstance->units,
                             &masstranInstance->feaProblem.numMaterial,
                             &masstranInstance->feaProblem.feaMaterial);
    AIM_STATUS(aimInfo, status);
  } else printf("Material tuple is NULL - No materials set\n");

  // Set property properties
  if (aimInputs[Property-1].nullVal == NotNull) {
    status = fea_getProperty(aimInfo,
                             aimInputs[Property-1].length,
                             aimInputs[Property-1].vals.tuple,
                             &masstranInstance->attrMap,
                             &masstranInstance->units,
                             &masstranInstance->feaProblem);
    AIM_STATUS(aimInfo, status);


    // Assign element "subtypes" based on properties set
    status = fea_assignElementSubType(masstranInstance->feaProblem.numProperty,
                                      masstranInstance->feaProblem.feaProperty,
                                      &masstranInstance->feaProblem.feaMesh);
    AIM_STATUS(aimInfo, status);

  } else printf("Property tuple is NULL - No properties set\n");

  nasMesh = &masstranInstance->feaProblem.feaMesh;
  massProp = &masstranInstance->massProp;

  // maps from nodeID to mesh->node index
  mesh_nodeID2Array(nasMesh, &n2a);
  AIM_NOTNULL(n2a, aimInfo, status);

  for (i = 0; i < nasMesh->numElement; i++) {

    // Get the property and material
    if (nasMesh->element[i].analysisType == MeshStructure) {
        feaData = (feaMeshDataStruct *) nasMesh->element[i].analysisData;
        propertyID = feaData->propertyID;
        //coordID = feaData->coordID;
        //elementSubType = feaData->elementSubType;
    } else {
        feaData = NULL;
        propertyID = nasMesh->element[i].markerID;
        //coordID = 0;
        //elementSubType = UnknownMeshSubElement;
    }
    elemWeight = 0;

    found = (int) false;
    for (j = 0; j < masstranInstance->feaProblem.numProperty; j++) {
        if (propertyID == masstranInstance->feaProblem.feaProperty[j].propertyID) {
            found = (int) true;
            break;
        }
    }

    if (found == (int) false) {
        printf("No property information found for element %d!\n",
               nasMesh->element[i].elementID);
        continue;
    }

    feaProperty = masstranInstance->feaProblem.feaProperty + j;

    if (nasMesh->element[i].elementType == Node ) {

      if (feaData != NULL && feaData->elementSubType == ConcentratedMassElement) {

        elem[0] = n2a[nasMesh->element[i].connectivity[0]];

        xcent = nasMesh->node[elem[0]].xyz[0] * Lscale + feaProperty->massOffset[0];
        ycent = nasMesh->node[elem[0]].xyz[1] * Lscale + feaProperty->massOffset[1];
        zcent = nasMesh->node[elem[0]].xyz[2] * Lscale + feaProperty->massOffset[2];

        elemArea = 0;

        elemWeight = feaProperty->mass;

        // add the inertia at the point
        Ixx  += feaProperty->massInertia[I11];
        Iyy  += feaProperty->massInertia[I22];
        Izz  += feaProperty->massInertia[I33];
        Ixy  -= feaProperty->massInertia[I21];
        Ixz  -= feaProperty->massInertia[I31];
        Iyz  -= feaProperty->massInertia[I32];
      } else {
        continue; // nothing to do if the node is not a concentrated mass
      }

    } else if ( nasMesh->element[i].elementType == Triangle) {

      elem[0] = n2a[nasMesh->element[i].connectivity[0]];
      elem[1] = n2a[nasMesh->element[i].connectivity[1]];
      elem[2] = n2a[nasMesh->element[i].connectivity[2]];

      xelem[0] = nasMesh->node[elem[0]].xyz[0] * Lscale;
      yelem[0] = nasMesh->node[elem[0]].xyz[1] * Lscale;
      zelem[0] = nasMesh->node[elem[0]].xyz[2] * Lscale;
      xelem[1] = nasMesh->node[elem[1]].xyz[0] * Lscale;
      yelem[1] = nasMesh->node[elem[1]].xyz[1] * Lscale;
      zelem[1] = nasMesh->node[elem[1]].xyz[2] * Lscale;
      xelem[2] = nasMesh->node[elem[2]].xyz[0] * Lscale;
      yelem[2] = nasMesh->node[elem[2]].xyz[1] * Lscale;
      zelem[2] = nasMesh->node[elem[2]].xyz[2] * Lscale;

      xcent  = (xelem[0] + xelem[1] + xelem[2]) / 3;
      ycent  = (yelem[0] + yelem[1] + yelem[2]) / 3;
      zcent  = (zelem[0] + zelem[1] + zelem[2]) / 3;

      dx1[0] = xelem[1] - xelem[0];
      dx1[1] = yelem[1] - yelem[0];
      dx1[2] = zelem[1] - zelem[0];
      dx2[0] = xelem[2] - xelem[0];
      dx2[1] = yelem[2] - yelem[0];
      dx2[2] = zelem[2] - zelem[0];

      CROSS(n, dx1, dx2);

      elemArea  = sqrt(DOT(n, n))/2.0;      /* 1/2 for area */

    } else if ( nasMesh->element[i].elementType == Quadrilateral) {

      elem[0] = n2a[nasMesh->element[i].connectivity[0]];
      elem[1] = n2a[nasMesh->element[i].connectivity[1]];
      elem[2] = n2a[nasMesh->element[i].connectivity[2]];
      elem[3] = n2a[nasMesh->element[i].connectivity[3]];

      xelem[0] = nasMesh->node[elem[0]].xyz[0] * Lscale;
      yelem[0] = nasMesh->node[elem[0]].xyz[1] * Lscale;
      zelem[0] = nasMesh->node[elem[0]].xyz[2] * Lscale;
      xelem[1] = nasMesh->node[elem[1]].xyz[0] * Lscale;
      yelem[1] = nasMesh->node[elem[1]].xyz[1] * Lscale;
      zelem[1] = nasMesh->node[elem[1]].xyz[2] * Lscale;
      xelem[2] = nasMesh->node[elem[2]].xyz[0] * Lscale;
      yelem[2] = nasMesh->node[elem[2]].xyz[1] * Lscale;
      zelem[2] = nasMesh->node[elem[2]].xyz[2] * Lscale;
      xelem[3] = nasMesh->node[elem[3]].xyz[0] * Lscale;
      yelem[3] = nasMesh->node[elem[3]].xyz[1] * Lscale;
      zelem[3] = nasMesh->node[elem[3]].xyz[2] * Lscale;

      xcent  = (xelem[0] + xelem[1] + xelem[2] + xelem[3]) / 4;
      ycent  = (yelem[0] + yelem[1] + yelem[2] + yelem[3]) / 4;
      zcent  = (zelem[0] + zelem[1] + zelem[2] + zelem[3]) / 4;

      dx1[0] = xelem[2] - xelem[0];
      dx1[1] = yelem[2] - yelem[0];
      dx1[2] = zelem[2] - zelem[0];
      dx2[0] = xelem[3] - xelem[1];
      dx2[1] = yelem[3] - yelem[1];
      dx2[2] = zelem[3] - zelem[1];

      CROSS(n, dx1, dx2);

      elemArea  = sqrt(DOT(n, n))/2.0;      /* 1/2 for area */

    } else {
      AIM_ERROR(aimInfo, "Unknown element type %d", nasMesh->element[i].elementType);
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    if (feaData == NULL || feaData->elementSubType != ConcentratedMassElement) {

      materialID = feaProperty->materialID;

      for (j = 0; j < masstranInstance->feaProblem.numMaterial; j++) {
        if (materialID == masstranInstance->feaProblem.feaMaterial[j].materialID) {
          break;
        }
      }

      if (found == (int) false) {
        printf("No material information found for element %d!\n", nasMesh->element[i].elementID);
        continue;
      }

      feaMaterial = masstranInstance->feaProblem.feaMaterial + j;

      if (feaMaterial->density > 0 && feaProperty->massPerArea) {
        AIM_ERROR(aimInfo, "Cannot specify both Material 'density' and Property 'massPerArea");
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      density = 1;
      if (feaMaterial->density > 0)
        density = feaMaterial->density;

      thick = 1;
      if (feaProperty->membraneThickness > 0) {
        thick   = feaProperty->membraneThickness;
      }

      if (feaProperty->massPerArea > 0) {
        density = feaProperty->massPerArea;
        thick = 1;
      }

      elemWeight = elemArea * density * thick;
    }

    area   += elemArea;
    mass   += elemWeight;

    cxmom  += xcent * elemArea;
    cymom  += ycent * elemArea;
    czmom  += zcent * elemArea;

    cgxmom += xcent * elemWeight;
    cgymom += ycent * elemWeight;
    cgzmom += zcent * elemWeight;

    Ixx    += (ycent * ycent + zcent * zcent) * elemWeight;
    Iyy    += (xcent * xcent + zcent * zcent) * elemWeight;
    Izz    += (xcent * xcent + ycent * ycent) * elemWeight;
    Ixy    += (xcent         * ycent        ) * elemWeight;
    Ixz    += (xcent         * zcent        ) * elemWeight;
    Iyz    += (ycent         * zcent        ) * elemWeight;
  }

  /* compute statistics for whole Body */
  Cx  = cxmom / area;
  Cy  = cymom / area;
  Cz  = czmom / area;

  CGx = cgxmom / mass;
  CGy = cgymom / mass;
  CGz = cgzmom / mass;

  Ixx -= mass * (CGy * CGy + CGz * CGz);
  Iyy -= mass * (CGx * CGx + CGz * CGz);
  Izz -= mass * (CGx * CGx + CGy * CGy);
  Ixy -= mass *  CGx * CGy;
  Ixz -= mass *  CGx * CGz;
  Iyz -= mass *  CGy * CGz;

  // store the mass properties
  massProp->area = area;
  massProp->mass = mass;

  massProp->Cx = Cx;
  massProp->Cy = Cy;
  massProp->Cz = Cz;

  massProp->CGx = CGx;
  massProp->CGy = CGy;
  massProp->CGz = CGz;

  massProp->Ixx = Ixx;
  massProp->Iyy = Iyy;
  massProp->Izz = Izz;
  massProp->Ixy = Ixy;
  massProp->Ixz = Ixz;
  massProp->Iyz = Iyz;

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(n2a);

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
int aimPostAnalysis(/*@unused@*/ void *instStore, /*@unused@*/ void *aimStruc,
                    /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  return CAPS_SUCCESS;
}


/* aimOutputs: Output Information for the AIM */
int aimOutputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
               /*@unused@*/ int index, char **aoname, capsValue *form)
{
  int status = CAPS_SUCCESS;

#ifdef DEBUG
  printf(" masstranAIM/aimOutputs  index  = %d!\n", index);
#endif
  feaUnitsStruct *units=NULL;
  aimStorage     *masstranInstance;
  char *tmpUnits = NULL;

  masstranInstance = (aimStorage *) instStore;
  AIM_NOTNULL(masstranInstance, aimInfo, status);

  units = &masstranInstance->units;

  /*! \page aimOutputsMasstran AIM Outputs
   * The following list outlines the Masstran outputs available through the AIM interface.
   */
  /*! \page aimOutputsMasstran AIM Outputs
   * - <B>Area</B> = Total area of the mesh.
   * - <B>Mass</B> = Total mass of the model.
   * - <B>Centroid</B> = Centroid of the model.
   * - <B>CG</B> = Center of gravity of the model.
   * - <B>Ixx</B> = Moment of inertia
   * - <B>Iyy</B> = Moment of inertia
   * - <B>Izz</B> = Moment of inertia
   * - <B>Ixy</B> = Moment of inertia
   * - <B>Izy</B> = Moment of inertia
   * - <B>Iyz</B> = Moment of inertia
   * - <B>I_Vector</B> = Moment of inertia vector
   *  \f[
   *    \vec{I} = \begin{bmatrix} I_{xx} & I_{yy} & I_{zz} & I_{xy} & I_{xz} & I_{yz} \end{bmatrix}
   *  \f]
   * - <B>I_Lower</B> = Moment of inertia lower triangular tensor
   * \f[
   *  \vec{I}_{lower} = \begin{bmatrix} I_{xx} & -I_{xy} & I_{yy} & -I_{xz} & -I_{yz} & I_{zz} \end{bmatrix},
   * \f]
   * - <B>I_Upper</B> = Moment of inertia upper triangular tensor
   * \f[
   *  \vec{I}_{upper} = \begin{bmatrix} I_{xx} & -I_{xy} & -I_{xz} & I_{yy} & -I_{yz} & I_{zz} \end{bmatrix},
   * \f]
   * - <B>I_Tensor</B> = Moment of inertia tensor
   * \f[
   * \bar{\bar{I}} =
   *  \begin{bmatrix}
   *    I_{xx} & -I_{xy} & -I_{xz} \\
   *   -I_{xy} &  I_{yy} & -I_{yz} \\
   *   -I_{xz} & -I_{yz} &  I_{zz}
   *  \end{bmatrix}
   * \f]
   */

  if (index == 1) {
      *aoname = EG_strdup("Area");
      form->dim  = Scalar;

      if (units->length != NULL) {
        status = aim_unitRaise(aimInfo, units->length, 2, &form->units); // length^2
        AIM_STATUS(aimInfo, status);
      }

  } else if (index == 2) {
      *aoname = EG_strdup("Mass");
      form->dim  = Scalar;

      if (units->mass != NULL) {
        AIM_STRDUP(form->units, units->mass, aimInfo, status);
      }

  } else if (index == 3) {
      *aoname = EG_strdup("Centroid");
      form->nrow = 3;
      form->dim  = Vector;

      if (units->length != NULL) {
        AIM_STRDUP(form->units, units->length, aimInfo, status);
      }

  } else if (index == 4) {
      *aoname = EG_strdup("CG");
      form->nrow = 3;
      form->dim  = Vector;

      if (units->length != NULL) {
        AIM_STRDUP(form->units, units->length, aimInfo, status);
      }

  } else if (index == 5) {
      *aoname = EG_strdup("Ixx");
      form->dim  = Scalar;

  } else if (index == 6) {
      *aoname = EG_strdup("Iyy");
      form->dim  = Scalar;

  } else if (index == 7) {
      *aoname = EG_strdup("Izz");
      form->dim  = Scalar;

  } else if (index == 8) {
      *aoname = EG_strdup("Ixy");
      form->dim  = Scalar;

  } else if (index == 9) {
      *aoname = EG_strdup("Ixz");
      form->dim  = Scalar;

  } else if (index == 10) {
      *aoname = EG_strdup("Iyz");
      form->dim  = Scalar;

  } else if (index == 11) {
      *aoname = EG_strdup("I_Vector");
      form->nrow = 6;
      form->dim  = Vector;

  } else if (index == 12) {
      *aoname = EG_strdup("I_Lower");
      form->nrow = 6;
      form->dim  = Vector;

  } else if (index == 13) {
      *aoname = EG_strdup("I_Upper");
      form->nrow = 6;
      form->dim  = Vector;

  } else if (index == 14) {
      *aoname = EG_strdup("I_Tensor");
      form->nrow = 9;
      form->dim  = Array2D;

  } else if (index == 15) {
      *aoname = EG_strdup("MassProp");
      form->type = String;
      form->nullVal = IsNull;
      goto cleanup;

  }

#if NUMOUTPUT != 15
#error "NUMOUTPUT is inconsistent with the list of inputs"
#endif

  if (*aoname == NULL) return EGADS_MALLOC;

  form->type       = Double;
  form->lfixed     = Fixed;
  form->sfixed     = Fixed;
  form->vals.reals = NULL;
  form->vals.real  = 0;

  if (form->nrow > 1) {
    AIM_ALLOC(form->vals.reals, form->nrow, double, aimInfo, status);
    memset(form->vals.reals, 0, form->nrow*sizeof(double));
  }

  if (index >= 5 && units->momentOfInertia != NULL) {
    AIM_STRDUP(form->units, units->momentOfInertia, aimInfo, status);
  }

  if (index == 14) {
      form->nrow     = 3;
      form->ncol     = 3;
  }

cleanup:
  AIM_FREE(tmpUnits);

  return status;
}


/* aimCalcOutput: Calculate/Retrieve Output Information */
int aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo,
                  /*@unused@*/ int index, capsValue *val)
{
  int status = CAPS_SUCCESS;
  double mass = 0, CGx = 0, CGy = 0, CGz = 0;
  double Ixx = 0;
  double Ixy = 0;
  double Izz = 0;
  double Ixz = 0;
  double Iyy = 0;
  double Iyz = 0;
  double *I;
  char massProp[512];
  feaUnitsStruct *units=NULL;
  aimStorage *masstranInstance;

#ifdef DEBUG
  const char *name;

  status = aim_getName(aimInfo, index, ANALYSISOUT, &name);
  printf(" masstranAIM/aimCalcOutput instance = %d  index = %d %s %d!\n",
         aim_getInstance(aimInfo), index, name, status);
#endif
  
  masstranInstance = (aimStorage *) instStore;

  units = &masstranInstance->units;

  if (val->length > 1) {
    AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);
  }

  if (index == 1) {
    val->vals.real = masstranInstance->massProp.area;

  } else if (index == 2) {
    val->vals.real = masstranInstance->massProp.mass;

  } else if (index == 3) {
    val->vals.reals[0] = masstranInstance->massProp.Cx;
    val->vals.reals[1] = masstranInstance->massProp.Cy;
    val->vals.reals[2] = masstranInstance->massProp.Cz;

  } else if (index == 4) {
    val->vals.reals[0] = masstranInstance->massProp.CGx;
    val->vals.reals[1] = masstranInstance->massProp.CGy;
    val->vals.reals[2] = masstranInstance->massProp.CGz;

  } else if (index == 5) {
    val->vals.real = masstranInstance->massProp.Ixx;

  } else if (index == 6) {
    val->vals.real = masstranInstance->massProp.Iyy;

  } else if (index == 7) {
    val->vals.real = masstranInstance->massProp.Izz;

  } else if (index == 8) {
    val->vals.real = masstranInstance->massProp.Ixy;

  } else if (index == 9) {
    val->vals.real = masstranInstance->massProp.Ixz;

  } else if (index == 10) {
    val->vals.real = masstranInstance->massProp.Iyz;

  } else if (index == 11) {
    val->vals.reals[0] = masstranInstance->massProp.Ixx;
    val->vals.reals[1] = masstranInstance->massProp.Iyy;
    val->vals.reals[2] = masstranInstance->massProp.Izz;
    val->vals.reals[3] = masstranInstance->massProp.Ixy;
    val->vals.reals[4] = masstranInstance->massProp.Ixz;
    val->vals.reals[5] = masstranInstance->massProp.Iyz;

  } else if (index == 12) {
    val->vals.reals[0] =  masstranInstance->massProp.Ixx;
    val->vals.reals[1] = -masstranInstance->massProp.Ixy;
    val->vals.reals[2] =  masstranInstance->massProp.Iyy;
    val->vals.reals[3] = -masstranInstance->massProp.Ixz;
    val->vals.reals[4] = -masstranInstance->massProp.Iyz;
    val->vals.reals[5] =  masstranInstance->massProp.Izz;

  } else if (index == 13) {
    val->vals.reals[0] =  masstranInstance->massProp.Ixx;
    val->vals.reals[1] = -masstranInstance->massProp.Ixy;
    val->vals.reals[2] = -masstranInstance->massProp.Ixz;
    val->vals.reals[3] =  masstranInstance->massProp.Iyy;
    val->vals.reals[4] = -masstranInstance->massProp.Iyz;
    val->vals.reals[5] =  masstranInstance->massProp.Izz;

  } else if (index == 14) {

    Ixx = masstranInstance->massProp.Ixx;
    Iyy = masstranInstance->massProp.Iyy;
    Iyz = masstranInstance->massProp.Iyz;
    Ixy = masstranInstance->massProp.Ixy;
    Ixz = masstranInstance->massProp.Ixz;
    Izz = masstranInstance->massProp.Izz;
    I = val->vals.reals;

    // populate the tensor
    I[0] =  Ixx; I[1] = -Ixy; I[2] = -Ixz;
    I[3] = -Ixy; I[4] =  Iyy; I[5] = -Iyz;
    I[6] = -Ixz; I[7] = -Iyz; I[8] =  Izz;

  } else if (index == 15) {

    mass = masstranInstance->massProp.mass;

    CGx = masstranInstance->massProp.CGx;
    CGy = masstranInstance->massProp.CGy;
    CGz = masstranInstance->massProp.CGz;

    Ixx = masstranInstance->massProp.Ixx;
    Iyy = masstranInstance->massProp.Iyy;
    Iyz = masstranInstance->massProp.Iyz;
    Ixy = masstranInstance->massProp.Ixy;
    Ixz = masstranInstance->massProp.Ixz;
    Izz = masstranInstance->massProp.Izz;

    if (units->mass != NULL) {

      snprintf(massProp, 512, "{\"mass\":[%20.14le, %s], \"CG\":[[%20.14le,%20.14le,%20.14le], %s], \"massInertia\":[[%20.14le, %20.14le, %20.14le, %20.14le, %20.14le, %20.14le], %s]}",
          mass, units->mass, CGx, CGy, CGz, units->length, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, units->momentOfInertia);

    } else {

      snprintf(massProp, 512, "{\"mass\":%20.14le, \"CG\":[%20.14le,%20.14le,%20.14le], \"massInertia\":[%20.14le, %20.14le, %20.14le, %20.14le, %20.14le, %20.14le]}",
          mass, CGx, CGy, CGz, Ixx, Iyy, Izz, Ixy, Ixz, Iyz);

    }

    AIM_STRDUP(val->vals.string, massProp, aimInfo, status);
  }

cleanup:
  return status;
}


/* aimCleanup: Free up the AIMs storage */
void aimCleanup(void *instStore)
{
  aimStorage *masstranInstance;
#ifdef DEBUG
  printf(" masstranAIM/aimCleanup!\n");
#endif
  
  masstranInstance = (aimStorage *) instStore;
  destroy_aimStorage(masstranInstance);
  EG_free(masstranInstance);
}
