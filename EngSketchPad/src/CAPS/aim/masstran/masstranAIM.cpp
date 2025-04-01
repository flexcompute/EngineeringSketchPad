/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             masstran AIM
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
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

#include "Surreal/SurrealS.h"

#include <string.h>
#include <math.h>
#include "aimUtil.h"

#include "meshUtils.h"
#include "miscUtils.h"
#include "feaUtils.h"

#ifdef WIN32
#define strcasecmp stricmp
#define strncasecmp _strnicmp
#define strtok_r   strtok_s
#endif

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
  inTess_Params = 1,                   /* index is 1-based */
  inEdge_Point_Min,
  inEdge_Point_Max,
  inQuad_Mesh,
  inProperty,
  inMaterial,
  inMesh_Morph,
  inSurface_Mesh,
  inDesign_Variable,
  inDesign_Variable_Relation,
  NUMINPUT = inDesign_Variable_Relation /* Total number of inputs */
};

enum aimOutputs
{
  outArea = 1,               /* index is 1-based */
  outMass,
  outCentroid,
  outCG,
  outIxx,
  outIyy,
  outIzz,
  outIxy,
  outIxz,
  outIyz,
  outI_Vector,
  outI_Lower,
  outI_Upper,
  outI_Tensor,
  outMassProp,
  outMassPropLink,
  NUMOUTPUT = outMassPropLink /* Total number of inputs */
};


template<class T>
struct sensNode {
  T xyz[3];
};

template<class T>
struct sensPropertyStruct {
  propertyTypeEnum propertyType;

  int propertyID; // ID number of property

  int materialID; // ID number of material

  // Shell
  T membraneThickness;  // Membrane thickness
  T massPerArea;        // Mass per unit area or Non-structural mass per unit area

  // Concentrated Mass
  T mass; // Mass value
  T massOffset[3]; // Offset distance from the grid point to the center of gravity
  T massInertia[6]; // Mass moment of inertia measured at the mass center of gravity
};

template<class T>
struct sensMaterialStruct {
  materialTypeEnum materialType;

  int materialID; // ID number of material

  T density;         // Rho - material mass density
};

template<class T>
struct massProperties {
  T area;
  T mass;
  T Ixx, Iyy, Izz;
  T Ixy, Ixz, Iyz;
  T Cx, Cy, Cz;
  T CGx, CGy, CGz;
};

/* AIM "local" per instance storage
   needed data should be added here & cleaned up in aimCleanup */
typedef struct {

  feaUnitsStruct units; // units system
  double Lscale;        // length scale

  feaProblemStruct feaProblem;

  // Attribute to index map
  mapAttrToIndexStruct attrMap;

  // Mesh holders
  int numMesh;
  meshStruct *feaMesh;

  // mass properties
  capsValue area;
  capsValue mass;
  capsValue Ixx, Iyy, Izz;
  capsValue Ixy, Ixz, Iyz;
  capsValue C;
  capsValue CG;

#define NUMVAR 10
  capsValue *values[NUMVAR];

  feaMassPropStruct feaMassProp;

} aimStorage;

static int initiate_aimStorage(void *aimInfo, aimStorage *masstranInstance)
{
    int status = CAPS_SUCCESS;
    int i;

    // Mesh holders
    masstranInstance->numMesh = 0;
    masstranInstance->feaMesh = NULL;

    masstranInstance->Lscale = 1.0;
    status = initiate_feaUnitsStruct(&masstranInstance->units);
    AIM_STATUS(aimInfo, status);

    // Container for attribute to index map
    status = initiate_mapAttrToIndexStruct(&masstranInstance->attrMap);
    AIM_STATUS(aimInfo, status);

    status = initiate_feaProblemStruct(&masstranInstance->feaProblem);
    AIM_STATUS(aimInfo, status);

    masstranInstance->values[0] = &masstranInstance->area;
    masstranInstance->values[1] = &masstranInstance->mass;
    masstranInstance->values[2] = &masstranInstance->Ixx;
    masstranInstance->values[3] = &masstranInstance->Iyy;
    masstranInstance->values[4] = &masstranInstance->Izz;
    masstranInstance->values[5] = &masstranInstance->Ixy;
    masstranInstance->values[6] = &masstranInstance->Ixz;
    masstranInstance->values[7] = &masstranInstance->Iyz;
    masstranInstance->values[8] = &masstranInstance->C;
    masstranInstance->values[9] = &masstranInstance->CG;
#if NUMVAR != 10
#error "Developer error! Update values array!"
#endif
    for (i = 0; i < NUMVAR; i++) {
        aim_initValue(masstranInstance->values[i]);
        masstranInstance->values[i]->type = Double;
    }

    AIM_ALLOC(masstranInstance->C.vals.reals, 3, double, aimInfo, status);
    masstranInstance->C.dim = Vector;
    masstranInstance->C.nrow = 3;
    masstranInstance->C.length = 3;

    AIM_ALLOC(masstranInstance->CG.vals.reals, 3, double, aimInfo, status);
    masstranInstance->CG.dim = Vector;
    masstranInstance->CG.nrow = 3;
    masstranInstance->CG.length = 3;

    initiate_feaMassPropStruct(&masstranInstance->feaMassProp);

cleanup:
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

    aim_freeValue(&masstranInstance->area);
    aim_freeValue(&masstranInstance->mass);
    aim_freeValue(&masstranInstance->Ixx);
    aim_freeValue(&masstranInstance->Iyy);
    aim_freeValue(&masstranInstance->Izz);
    aim_freeValue(&masstranInstance->Ixy);
    aim_freeValue(&masstranInstance->Ixz);
    aim_freeValue(&masstranInstance->Iyz);
    aim_freeValue(&masstranInstance->C);
    aim_freeValue(&masstranInstance->CG);

    // units are borrowed
    masstranInstance->feaMassProp.massUnit = NULL;
    masstranInstance->feaMassProp.lengthUnit = NULL;
    masstranInstance->feaMassProp.momentOfInertiaUnit = NULL;
    destroy_feaMassPropStruct(&masstranInstance->feaMassProp);

    return CAPS_SUCCESS;
}


static int checkAndCreateMesh(void *aimInfo, aimStorage *masstranInstance)
{
  // Function return flag
  int status = CAPS_SUCCESS;
  int i, remesh = (int)false;

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
      remesh = remesh || (masstranInstance->feaMesh[i].egadsTess->oclass != TESSELLATION);
  }

  // retrieve or create the mesh from fea_createMesh
  status = aim_getValue(aimInfo, inTess_Params, ANALYSISIN, &TessParams);
  AIM_STATUS(aimInfo, status);

  remesh = remesh || (aim_newAnalysisIn(aimInfo, inTess_Params) == CAPS_SUCCESS);

  status = aim_getValue(aimInfo, inEdge_Point_Min, ANALYSISIN, &EdgePoint_Min);
  AIM_STATUS(aimInfo, status);

  remesh = remesh || (aim_newAnalysisIn(aimInfo, inEdge_Point_Min) == CAPS_SUCCESS);

  status = aim_getValue(aimInfo, inEdge_Point_Max, ANALYSISIN, &EdgePoint_Max);
  AIM_STATUS(aimInfo, status);

  remesh = remesh || (aim_newAnalysisIn(aimInfo, inEdge_Point_Max) == CAPS_SUCCESS);

  status = aim_getValue(aimInfo, inQuad_Mesh, ANALYSISIN, &QuadMesh);
  AIM_STATUS(aimInfo, status);

  remesh = remesh || (aim_newAnalysisIn(aimInfo, inQuad_Mesh) == CAPS_SUCCESS);

  remesh = remesh || (aim_newGeometry(aimInfo) == CAPS_SUCCESS);
  remesh = remesh || (masstranInstance->feaProblem.feaMesh.numNode == 0);

  if (remesh == (int) false) return CAPS_SUCCESS;

  if (TessParams != NULL) {
      tessParam[0] = TessParams->vals.reals[0];
      tessParam[1] = TessParams->vals.reals[1];
      tessParam[2] = TessParams->vals.reals[2];
  }

  // Max and min number of points
  if (EdgePoint_Min != NULL && EdgePoint_Min->nullVal != IsNull) {
      edgePointMin = EdgePoint_Min->vals.integer;
      if (edgePointMin < 2) {
        AIM_ANALYSISIN_ERROR(aimInfo, inEdge_Point_Min, "Edge_Point_Min = %d must be greater or equal to 2\n", edgePointMin);
        return CAPS_BADVALUE;
      }
  }

  if (EdgePoint_Max != NULL && EdgePoint_Max->nullVal != IsNull) {
      edgePointMax = EdgePoint_Max->vals.integer;
      if (edgePointMax < 2) {
        AIM_ANALYSISIN_ERROR(aimInfo, inEdge_Point_Max, "Edge_Point_Max = %d must be greater or equal to 2\n", edgePointMax);
        return CAPS_BADVALUE;
      }
  }

  if (edgePointMin >= 2 && edgePointMax >= 2 && edgePointMin > edgePointMax) {
    AIM_ERROR  (aimInfo, "Edge_Point_Max must be greater or equal Edge_Point_Min");
    AIM_ADDLINE(aimInfo, "Edge_Point_Max = %d, Edge_Point_Min = %d\n",edgePointMax,edgePointMin);
    return CAPS_BADVALUE;
  }

  if (QuadMesh != NULL) quadMesh = 2*QuadMesh->vals.integer;

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
                          NULL, NULL,
                          &masstranInstance->numMesh,
                          &masstranInstance->feaMesh,
                          &masstranInstance->feaProblem );
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
extern "C" int
aimInitialize(int inst, const char *unitSys, void *aimInfo,
              void **instStore, /*@unused@*/ int *major,
              /*@unused@*/ int *minor, int *nIn, int *nOut,
              int *nFields, char ***fnames, int **franks, int **fInOut)
{
  int status = CAPS_SUCCESS;
  aimStorage *masstranInstance=NULL;
  const char *keyWord;
  char *keyValue = NULL;
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
  AIM_ALLOC(masstranInstance, 1, aimStorage, aimInfo, status);

  initiate_aimStorage(aimInfo, masstranInstance);
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

    // Base units not used by this AIM
    AIM_STRDUP(units->time, "s", aimInfo, status);
    AIM_STRDUP(units->temperature, "Kelvin", aimInfo, status);

    status = fea_feaDerivedUnits(aimInfo, units);
    AIM_STATUS(aimInfo, status);

    // set units
    status = aim_unitRaise(aimInfo, units->length, 2, &masstranInstance->area.units); // length^2
    AIM_STATUS(aimInfo, status);
    AIM_STRDUP(masstranInstance->mass.units, units->mass           , aimInfo, status); // mass
    AIM_STRDUP(masstranInstance->C.units   , units->length         , aimInfo, status); // lenth
    AIM_STRDUP(masstranInstance->CG.units  , units->length         , aimInfo, status); // length
    AIM_STRDUP(masstranInstance->Ixx.units , units->momentOfInertia, aimInfo, status); // moment of inertia
    AIM_STRDUP(masstranInstance->Iyy.units , units->momentOfInertia, aimInfo, status); // moment of inertia
    AIM_STRDUP(masstranInstance->Izz.units , units->momentOfInertia, aimInfo, status); // moment of inertia
    AIM_STRDUP(masstranInstance->Ixy.units , units->momentOfInertia, aimInfo, status); // moment of inertia
    AIM_STRDUP(masstranInstance->Ixz.units , units->momentOfInertia, aimInfo, status); // moment of inertia
    AIM_STRDUP(masstranInstance->Iyz.units , units->momentOfInertia, aimInfo, status); // moment of inertia
  }

cleanup:
  AIM_FREE(keyValue);

  return status;
}


// ********************** AIM Function Break *****************************
/* aimInputs: Input Information for the AIM */
extern "C" int
aimInputs(/*@unused@*/ void *instStore, /*@unused@*/ void *aimInfo,
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
  if (index == inTess_Params) {
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

  } else if (index == inEdge_Point_Min) {
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

  } else if (index == inEdge_Point_Max) {
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

  } else if (index == inQuad_Mesh) {
    *ainame               = EG_strdup("Quad_Mesh");
    defval->type          = Boolean;
    defval->vals.integer  = (int) false;

    /*! \page aimInputsMasstran
     * - <B> Quad_Mesh = False</B> <br>
     * Create a quadratic mesh on four edge faces when creating the boundary element model.
     */

  } else if (index == inProperty) {
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
  } else if (index == inMaterial) {
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

  } else if (index == inMesh_Morph) {
      *ainame              = EG_strdup("Mesh_Morph");
      defval->type         = Boolean;
      defval->lfixed       = Fixed;
      defval->vals.integer = (int) false;
      defval->dim          = Scalar;
      defval->nullVal      = NotNull;

      /*! \page aimInputsMasstran
       * - <B> Mesh_Morph = False</B> <br>
       * Project previous surface mesh onto new geometry.
       */

  } else if (index == inSurface_Mesh) {
      *ainame             = EG_strdup("Surface_Mesh");
      defval->type        = PointerMesh;
      defval->dim         = Vector;
      defval->lfixed      = Change;
      defval->sfixed      = Fixed;
      defval->vals.AIMptr = NULL;
      defval->nullVal     = IsNull;

      /*! \page aimInputsMasstran
       * - <B>Surface_Mesh = NULL</B> <br>
       * A Surface_Mesh link.
       */

  } else if (index == inDesign_Variable) {
      *ainame              = EG_strdup("Design_Variable");
      defval->type         = Tuple;
      defval->nullVal      = IsNull;
      defval->lfixed       = Change;
      defval->vals.tuple   = NULL;
      defval->dim          = Vector;

      /*! \page aimInputsMasstran
       * - <B> Design_Variable = NULL</B> <br>
       * The design variable tuple is used to input design variable information for the model optimization, see \ref feaDesignVariable for additional details.
       */

  } else if (index == inDesign_Variable_Relation) {
      *ainame              = EG_strdup("Design_Variable_Relation");
      defval->type         = Tuple;
      defval->nullVal      = IsNull;
      //defval->units        = NULL;
      defval->lfixed       = Change;
      defval->vals.tuple   = NULL;
      defval->dim          = Vector;

      /*! \page aimInputsMasstran
       * - <B> Design_Variable_Relation = NULL</B> <br>
       * The design variable relation tuple is used to input design variable relation information for the model optimization, see \ref feaDesignVariableRelation for additional details.
       */

  }

  AIM_NOTNULL(*ainame, aimInfo, status);

cleanup:
  if (status != CAPS_SUCCESS) AIM_FREE(*ainame);
  return status;
}


// ********************** AIM Function Break *****************************
extern "C" int
aimUpdateState(void *instStore, void *aimInfo,
               capsValue *aimInputs)
{
  int status; // Function return status
  int i, j, idv, len, index, len_wrt;

  const feaUnitsStruct *units=NULL;

  int numBody = 0;
  ego *bodies=NULL;
  const char *intents=NULL;

  const char *bodyLunits=NULL, *name=NULL;
  capsValue *geomInVal, **values;

  aimStorage *masstranInstance;

  masstranInstance = (aimStorage *) instStore;
  AIM_NOTNULL(aimInputs, aimInfo, status);

  // Get FEA mesh if we don't already have one
  status = checkAndCreateMesh(aimInfo, masstranInstance);
  AIM_STATUS(aimInfo, status);

  // Note: Setting order is important here.
  // 1. Materials should be set before properties.
  // 2. Coordinate system should be set before mesh and loads
  // 3. Mesh should be set before loads, constraints, supports, and connections

  // Set material properties
  if (aimInputs[inMaterial-1].nullVal == NotNull) {
    if (aim_newAnalysisIn(aimInfo, inMaterial) == CAPS_SUCCESS ||
        masstranInstance->feaProblem.numMaterial == 0) {
      status = fea_getMaterial(aimInfo,
                               aimInputs[inMaterial-1].length,
                               aimInputs[inMaterial-1].vals.tuple,
                               &masstranInstance->units,
                               &masstranInstance->feaProblem.numMaterial,
                               &masstranInstance->feaProblem.feaMaterial);
      AIM_STATUS(aimInfo, status);
    }
  } else {
    printf("Material tuple is NULL - No materials set\n");

    if (masstranInstance->feaProblem.feaMaterial != NULL) {
        for (i = 0; i < masstranInstance->feaProblem.numMaterial; i++) {
            status = destroy_feaMaterialStruct(&masstranInstance->feaProblem.feaMaterial[i]);
            if (status != CAPS_SUCCESS) printf("Status %d during destroy_feaMaterialStruct\n", status);
        }
    }
    AIM_FREE(masstranInstance->feaProblem.feaMaterial);
    masstranInstance->feaProblem.numMaterial = 0;
  }

  // Set property properties
  if (aimInputs[inProperty-1].nullVal == NotNull) {

    if (aim_newAnalysisIn(aimInfo, inProperty) == CAPS_SUCCESS ||
        masstranInstance->feaProblem.numProperty == 0) {
      status = fea_getProperty(aimInfo,
                               aimInputs[inProperty-1].length,
                               aimInputs[inProperty-1].vals.tuple,
                               &masstranInstance->attrMap,
                               &masstranInstance->units,
                               &masstranInstance->feaProblem);
      AIM_STATUS(aimInfo, status);
    }

    // Assign element "subtypes" based on properties set
    status = fea_assignElementSubType(masstranInstance->feaProblem.numProperty,
                                      masstranInstance->feaProblem.feaProperty,
                                      &masstranInstance->feaProblem.feaMesh);
    AIM_STATUS(aimInfo, status);

  } else {
    printf("Property tuple is NULL - No properties set\n");
    for (i = 0; i < masstranInstance->feaProblem.numProperty; i++) {
      status = destroy_feaPropertyStruct(&masstranInstance->feaProblem.feaProperty[i]);
      AIM_STATUS(aimInfo, status);
    }
    AIM_FREE(masstranInstance->feaProblem.feaProperty);
    masstranInstance->feaProblem.numProperty = 0;
  }

  // Set design variables
  if (aimInputs[inDesign_Variable-1].nullVal == NotNull &&
      (aim_newAnalysisIn(aimInfo, inDesign_Variable) == CAPS_SUCCESS ||
       aim_newAnalysisIn(aimInfo, inDesign_Variable_Relation) == CAPS_SUCCESS ||
       masstranInstance->feaProblem.numDesignVariable == 0 ||
       masstranInstance->feaProblem.numDesignVariableRelation == 0)) {
    status = fea_getDesignVariable(aimInfo,
                                   (int)false,
                                   aimInputs[inDesign_Variable-1].length,
                                   aimInputs[inDesign_Variable-1].vals.tuple,
                                   aimInputs[inDesign_Variable_Relation-1].length,
                                   aimInputs[inDesign_Variable_Relation-1].vals.tuple,
                                   &masstranInstance->attrMap,
                                   &masstranInstance->feaProblem);
    AIM_STATUS(aimInfo, status);

    /* allocate derivatives */
    values = masstranInstance->values;
    for (i = 0; i < NUMVAR; i++) {
      values[i]->type = DoubleDeriv;

      if (values[i]->derivs != NULL) {
        for (idv = 0; idv < values[i]->nderiv; idv++) {
          AIM_FREE(values[i]->derivs[idv].name);
          AIM_FREE(values[i]->derivs[idv].deriv);
        }
        AIM_FREE(values[i]->derivs);
        values[i]->nderiv = 0;
      }

      AIM_ALLOC(values[i]->derivs, masstranInstance->feaProblem.numDesignVariable, capsDeriv, aimInfo, status);
      for (idv = 0; idv < masstranInstance->feaProblem.numDesignVariable; idv++) {
        values[i]->derivs[idv].name  = NULL;
        values[i]->derivs[idv].deriv = NULL;
        values[i]->derivs[idv].len_wrt = 0;
      }
      values[i]->nderiv = masstranInstance->feaProblem.numDesignVariable;
    }

    for (idv = 0; idv < masstranInstance->feaProblem.numDesignVariable; idv++) {

      name = masstranInstance->feaProblem.feaDesignVariable[idv].name;

      // Loop over the geometry in values and compute sensitivities for all bodies
      index = aim_getIndex(aimInfo, name, GEOMETRYIN);
      status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
      if (status == CAPS_BADINDEX) {
        len_wrt = 1;
      } else {
        len_wrt = geomInVal->length;
      }

      // setup derivative memory
      for (i = 0; i < NUMVAR; i++) {
        AIM_FREE(values[i]->derivs[idv].name);
        AIM_STRDUP(values[i]->derivs[idv].name, name, aimInfo, status);

        values[i]->derivs[idv].len_wrt = len_wrt;

        len = values[i]->length*values[i]->derivs[idv].len_wrt;
        AIM_ALLOC(values[i]->derivs[idv].deriv, len, double, aimInfo, status);
        for (j = 0; j < len; j++)
          values[i]->derivs[idv].deriv[j] = 0;
      }
    }

  } else if (aimInputs[inDesign_Variable-1].nullVal == IsNull){
    for (i = 0; i < masstranInstance->feaProblem.numDesignVariable; i++) {
      status = destroy_feaDesignVariableStruct(&masstranInstance->feaProblem.feaDesignVariable[i]);
      AIM_STATUS(aimInfo, status);
    }
    AIM_FREE(masstranInstance->feaProblem.feaDesignVariable);
    masstranInstance->feaProblem.numDesignVariable = 0;
    for (i = 0; i < masstranInstance->feaProblem.numDesignVariableRelation; i++) {
      status = destroy_feaDesignVariableRelationStruct(&masstranInstance->feaProblem.feaDesignVariableRelation[i]);
      AIM_STATUS(aimInfo, status);
    }
    AIM_FREE(masstranInstance->feaProblem.feaDesignVariableRelation);
    masstranInstance->feaProblem.numDesignVariableRelation = 0;

    /* remove derivatives */
    values = masstranInstance->values;
    for (i = 0; i < NUMVAR; i++) {
      values[i]->type = Double;

      for (idv = 0; idv < values[i]->nderiv; idv++) {
        AIM_FREE(values[i]->derivs[idv].name);
        AIM_FREE(values[i]->derivs[idv].deriv);
      }
      AIM_FREE(values[i]->derivs);
      values[i]->nderiv = 0;
    }
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
    status = aim_capsLength(aimInfo, &bodyLunits);
    if (status != CAPS_SUCCESS) {
      AIM_ERROR(aimInfo, "capsLength is not set in *.csm file!");
      status = CAPS_BADVALUE;
      goto cleanup;
    }

    // conversion of the csm model units into units of Lunits
    masstranInstance->Lscale = 1.0;
    status = aim_convert(aimInfo, 1, bodyLunits, &masstranInstance->Lscale, units->length, &masstranInstance->Lscale);
    AIM_STATUS(aimInfo, status);
  } else {
    status = destroy_feaUnitsStruct(&masstranInstance->units);
    AIM_STATUS(aimInfo, status);
  }


  status = CAPS_SUCCESS;
cleanup:

  return status;
}


// ********************** AIM Function Break *****************************
template<class T, class NodeType>
int
calcMassProps(void *aimInfo,
              const int numProperty,
              const sensPropertyStruct<T> *feaProperties,
              const int numMaterial,
              const sensMaterialStruct<T> *feaMaterials,
              const int *n2a,
              const meshStruct *nasMesh,
              const double Lscale,
              const NodeType *node,
              massProperties<T>& prop)
{
  int i, j;   // Indexing
  int status = CAPS_SUCCESS; // Status return

  T area;
  T mass;
  T Ixx, Ixy, Izz;
  T Ixz, Iyy, Iyz;
  T CGx, CGy, CGz;

  T cgxmom = 0;
  T cgymom = 0;
  T cgzmom = 0;
  T cxmom  = 0;
  T cymom  = 0;
  T czmom  = 0;
  T xelem[4], yelem[4], zelem[4], elemArea;
  T xcent, ycent, zcent;
  T dx1[3], dx2[3], n[3];
  T thick, density;
  T elemWeight;

  int elem[4];

  feaMeshDataStruct *feaData = NULL;
  int propertyID, materialID;
  const sensPropertyStruct<T> *feaProperty;
  const sensMaterialStruct<T> *feaMaterial;


  area = 0;
  mass = 0;
  Ixx = Ixy = Izz = 0;
  Ixz = Iyy = Iyz = 0;
  CGx = CGy = CGz = 0;

  for (i = 0; i < nasMesh->numElement; i++) {

    // Get the property and material
    if (nasMesh->element[i].analysisType == MeshStructure) {
        feaData = (feaMeshDataStruct *) nasMesh->element[i].analysisData;
        propertyID = feaData->propertyID;
    } else {
        feaData = NULL;
        propertyID = nasMesh->element[i].markerID;
    }
    elemWeight = 0;

    feaProperty = NULL;
    for (j = 0; j < numProperty; j++) {
        if (propertyID == feaProperties[j].propertyID) {
            feaProperty = &feaProperties[j];
            break;
        }
    }

    if (feaProperty == NULL) {
        AIM_ERROR(aimInfo, "No property ID %d information found for element %d!\n",
                  propertyID, nasMesh->element[i].elementID);
        status = CAPS_BADVALUE;
        goto cleanup;
    }

    if (nasMesh->element[i].elementType == Node ) {

      if (feaData != NULL && feaData->elementSubType == ConcentratedMassElement) {

        elem[0] = n2a[nasMesh->element[i].connectivity[0]];

        xcent = (node[elem[0]].xyz[0] + feaProperty->massOffset[0]) * Lscale;
        ycent = (node[elem[0]].xyz[1] + feaProperty->massOffset[1]) * Lscale;
        zcent = (node[elem[0]].xyz[2] + feaProperty->massOffset[2]) * Lscale;

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

      xelem[0] = node[elem[0]].xyz[0] * Lscale;
      yelem[0] = node[elem[0]].xyz[1] * Lscale;
      zelem[0] = node[elem[0]].xyz[2] * Lscale;
      xelem[1] = node[elem[1]].xyz[0] * Lscale;
      yelem[1] = node[elem[1]].xyz[1] * Lscale;
      zelem[1] = node[elem[1]].xyz[2] * Lscale;
      xelem[2] = node[elem[2]].xyz[0] * Lscale;
      yelem[2] = node[elem[2]].xyz[1] * Lscale;
      zelem[2] = node[elem[2]].xyz[2] * Lscale;

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

      xelem[0] = node[elem[0]].xyz[0] * Lscale;
      yelem[0] = node[elem[0]].xyz[1] * Lscale;
      zelem[0] = node[elem[0]].xyz[2] * Lscale;
      xelem[1] = node[elem[1]].xyz[0] * Lscale;
      yelem[1] = node[elem[1]].xyz[1] * Lscale;
      zelem[1] = node[elem[1]].xyz[2] * Lscale;
      xelem[2] = node[elem[2]].xyz[0] * Lscale;
      yelem[2] = node[elem[2]].xyz[1] * Lscale;
      zelem[2] = node[elem[2]].xyz[2] * Lscale;
      xelem[3] = node[elem[3]].xyz[0] * Lscale;
      yelem[3] = node[elem[3]].xyz[1] * Lscale;
      zelem[3] = node[elem[3]].xyz[2] * Lscale;

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

      feaMaterial = NULL;
      for (j = 0; j < numMaterial; j++) {
        if (materialID == feaMaterials[j].materialID) {
          feaMaterial = &feaMaterials[j];
          break;
        }
      }

      if (feaMaterial == NULL) {
        AIM_ERROR(aimInfo, "No material information found for element %d!\n", nasMesh->element[i].elementID);
        status = CAPS_BADVALUE;
        goto cleanup;
      }

      if (feaMaterial->density > 0 && feaProperty->massPerArea > 0) {
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
  prop.area = area;
  prop.mass = mass;

  prop.Cx  = cxmom / area;
  prop.Cy  = cymom / area;
  prop.Cz  = czmom / area;

  prop.CGx = CGx = cgxmom / mass;
  prop.CGy = CGy = cgymom / mass;
  prop.CGz = CGz = cgzmom / mass;

  prop.Ixx = Ixx - mass * (CGy * CGy + CGz * CGz);
  prop.Iyy = Iyy - mass * (CGx * CGx + CGz * CGz);
  prop.Izz = Izz - mass * (CGx * CGx + CGy * CGy);
  prop.Ixy = Ixy - mass *  CGx * CGy;
  prop.Ixz = Ixz - mass *  CGx * CGz;
  prop.Iyz = Iyz - mass *  CGy * CGz;

cleanup:
  return status;
}

// ********************** AIM Function Break *****************************
template<class T>
void setDesignVariableRelation(feaDesignVariableRelationStruct *designRelation, T *var)
{
  int idv;
  double val;

  *var = designRelation->constantRelationCoeff;
  for (idv = 0; idv < designRelation->numDesignVariable; idv++) {
    val = designRelation->designVariableSet[idv]->initialValue;
    *var += designRelation->linearRelationCoeff[idv]*val;
  }
}

// ********************** AIM Function Break *****************************
template<class T>
int getComponents(void *aimInfo, const aimStorage *masstranInstance,
                  sensPropertyStruct<T> **feaProperty_out,
                  sensMaterialStruct<T> **feaMaterial_out)
{
  int status = CAPS_SUCCESS;
  int i, j, idr, ip, im;

  sensPropertyStruct<T> *feaProperty = NULL;
  sensMaterialStruct<T> *feaMaterial = NULL;

  feaDesignVariableRelationStruct *designRelation = NULL;

  feaProperty = new sensPropertyStruct<T>[masstranInstance->feaProblem.numProperty];

  for (i = 0; i < masstranInstance->feaProblem.numProperty; i++) {
    feaProperty[i].propertyType      = masstranInstance->feaProblem.feaProperty[i].propertyType;
    feaProperty[i].propertyID        = masstranInstance->feaProblem.feaProperty[i].propertyID;
    feaProperty[i].materialID        = masstranInstance->feaProblem.feaProperty[i].materialID;
    feaProperty[i].membraneThickness = masstranInstance->feaProblem.feaProperty[i].membraneThickness;
    feaProperty[i].massPerArea       = masstranInstance->feaProblem.feaProperty[i].massPerArea;
    feaProperty[i].mass              = masstranInstance->feaProblem.feaProperty[i].mass;
    for (j = 0; j < 3; j++) {
      feaProperty[i].massOffset[j] = masstranInstance->feaProblem.feaProperty[i].massOffset[j];
    }
    for (j = 0; j < 6; j++) {
      feaProperty[i].massInertia[j] = masstranInstance->feaProblem.feaProperty[i].massInertia[j];
    }
  }

  feaMaterial = new sensMaterialStruct<T>[masstranInstance->feaProblem.numMaterial];

  for (i = 0; i < masstranInstance->feaProblem.numMaterial; i++) {
    feaMaterial[i].materialType = masstranInstance->feaProblem.feaMaterial[i].materialType;
    feaMaterial[i].materialID   = masstranInstance->feaProblem.feaMaterial[i].materialID;
    feaMaterial[i].density      = masstranInstance->feaProblem.feaMaterial[i].density;
  }

  /* set Design Variable Relation values */
  for (idr = 0; idr < masstranInstance->feaProblem.numDesignVariableRelation; idr++) {

    designRelation = &masstranInstance->feaProblem.feaDesignVariableRelation[idr];

    for (ip = 0; ip < designRelation->numPropertyID; ip++) {
      for (j = 0; j < masstranInstance->feaProblem.numProperty; j++) {
        if (designRelation->propertySetID[ip] == feaProperty[j].propertyID) {
          if (strncasecmp(designRelation->fieldName, "membraneThickness", 17) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].membraneThickness);
          } else if (strncasecmp(designRelation->fieldName, "massPerArea", 11) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massPerArea);
          } else if (strncasecmp(designRelation->fieldName, "mass", 4) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].mass);
          } else if (strncasecmp(designRelation->fieldName, "massOffset1", 11) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massOffset[0]);
          } else if (strncasecmp(designRelation->fieldName, "massOffset2", 11) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massOffset[1]);
          } else if (strncasecmp(designRelation->fieldName, "massOffset3", 11) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massOffset[2]);
          } else if (strncasecmp(designRelation->fieldName, "Ixx", 3) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massInertia[0]);
          } else if (strncasecmp(designRelation->fieldName, "Iyy", 3) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massInertia[1]);
          } else if (strncasecmp(designRelation->fieldName, "Izz", 3) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massInertia[2]);
          } else if (strncasecmp(designRelation->fieldName, "Ixy", 3) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massInertia[3]);
          } else if (strncasecmp(designRelation->fieldName, "Ixz", 3) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massInertia[4]);
          } else if (strncasecmp(designRelation->fieldName, "Iyz", 3) == 0) {
            setDesignVariableRelation(designRelation, &feaProperty[j].massInertia[5]);
          } else {
            AIM_ERROR(aimInfo, "Unknown Design_Variable_Relation fieldName = %s", designRelation->fieldName);
            status = CAPS_BADVALUE;
            goto cleanup;
          }
          break;
        }
      }
    }

    for (im = 0; im < designRelation->numMaterialID; im++) {
      for (j = 0; j < masstranInstance->feaProblem.numMaterial; j++) {
        if (designRelation->materialSetID[im] == feaMaterial[j].materialID) {
          if (strncasecmp(designRelation->fieldName, "density", 7) == 0) {
            setDesignVariableRelation(designRelation, &feaMaterial[j].density);
          } else {
            AIM_ERROR(aimInfo, "Unknown Design_Variable_Relation fieldName = %s", designRelation->fieldName);
            status = CAPS_BADVALUE;
            goto cleanup;
          }
          break;
        }
      }
    }
  }

  *feaProperty_out = feaProperty;
  *feaMaterial_out = feaMaterial;

cleanup:
  if (status != CAPS_SUCCESS) {
    delete [] feaProperty;
    delete [] feaMaterial;
    *feaProperty_out = NULL;
    *feaMaterial_out = NULL;
  }
  return status;
}

// ********************** AIM Function Break *****************************
extern "C" int
aimPreAnalysis(const void *instStore, void *aimInfo, capsValue *aimInputs)
{
  int status; // Status return

  int *n2a = NULL;

  FILE *fp=NULL;

  const meshStruct *nasMesh;
  const aimStorage *masstranInstance;

  massProperties<double> massProp;

  sensPropertyStruct<double> *feaProperty = NULL;
  sensMaterialStruct<double> *feaMaterial = NULL;

  masstranInstance = (const aimStorage *) instStore;
  AIM_NOTNULL(aimInputs, aimInfo, status);

  if ( aimInputs[inMesh_Morph-1].vals.integer == (int) true &&
       aimInputs[inSurface_Mesh-1].nullVal == NotNull) { // If we are mighty morphing
      // store the current mesh for future iterations
      status = aim_storeMeshRef(aimInfo, (aimMeshRef *) aimInputs[inSurface_Mesh-1].vals.AIMptr, NULL);
      AIM_STATUS(aimInfo, status);
  }

  status = getComponents(aimInfo, masstranInstance, &feaProperty, &feaMaterial);
  AIM_STATUS(aimInfo, status);

  nasMesh = &masstranInstance->feaProblem.feaMesh;

  // maps from nodeID to mesh->node index
  mesh_nodeID2Array(nasMesh, &n2a);
  AIM_NOTNULL(n2a, aimInfo, status);

  status = calcMassProps(aimInfo,
                         masstranInstance->feaProblem.numProperty,
                         feaProperty,
                         masstranInstance->feaProblem.numMaterial,
                         feaMaterial,
                         n2a,
                         nasMesh,
                         masstranInstance->Lscale,
                         nasMesh->node,
                         massProp);
  AIM_STATUS(aimInfo, status);

  // store the mass properties
  fp = aim_fopen(aimInfo, "masstran.out", "w");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Failed to open masstran.out");
    status = CAPS_IOERR;
    goto cleanup;
  }

  fprintf(fp, "%16.12e\n", massProp.area);
  fprintf(fp, "%16.12e\n", massProp.mass);

  fprintf(fp, "%16.12e\n", massProp.Cx);
  fprintf(fp, "%16.12e\n", massProp.Cy);
  fprintf(fp, "%16.12e\n", massProp.Cz);

  fprintf(fp, "%16.12e\n", massProp.CGx);
  fprintf(fp, "%16.12e\n", massProp.CGy);
  fprintf(fp, "%16.12e\n", massProp.CGz);

  fprintf(fp, "%16.12e\n", massProp.Ixx);
  fprintf(fp, "%16.12e\n", massProp.Iyy);
  fprintf(fp, "%16.12e\n", massProp.Izz);
  fprintf(fp, "%16.12e\n", massProp.Ixy);
  fprintf(fp, "%16.12e\n", massProp.Ixz);
  fprintf(fp, "%16.12e\n", massProp.Iyz);

  fclose(fp); fp = NULL;

  status = CAPS_SUCCESS;

cleanup:
  AIM_FREE(n2a);

  delete [] feaProperty;
  delete [] feaMaterial;

  if (fp != NULL) fclose(fp);

  return status;
}


// ********************** AIM Function Break *****************************
extern "C" int
aimExecute(/*@unused@*/ const void *instStore, /*@unused@*/ void *aimInfo,
           int *state)
{
  *state = 0;
  return CAPS_SUCCESS;
}


// ********************** AIM Function Break *****************************
extern "C" int
aimPostAnalysis( void *instStore, /*@unused@*/ void *aimInfo,
                 /*@unused@*/ int restart, /*@unused@*/ capsValue *inputs)
{
  int status = CAPS_SUCCESS;

  int i, j, d, idv, idr, irow, icol, ibody, im, ip; // Indexing
  int index;
  int *n2a = NULL;

  int numNode = 0;

  const char *name;
  double *dxyz = NULL;

  feaDesignVariableStruct *designVariable = NULL;
  feaDesignVariableRelationStruct *designRelation = NULL;

  capsValue *geomInVal;
  sensNode<SurrealS<1>> *sensNodes=NULL;

  FILE *fp = NULL;
  massProperties<SurrealS<1>> massPropS;

  sensPropertyStruct<SurrealS<1>> *feaProperty = NULL;
  sensMaterialStruct<SurrealS<1>> *feaMaterial = NULL;

  aimStorage *masstranInstance;

  masstranInstance = (aimStorage *) instStore;

  fp = aim_fopen(aimInfo, "masstran.out", "r");
  if (fp == NULL) {
    AIM_ERROR(aimInfo, "Failed to open masstran.out");
    return CAPS_IOERR;
  }

  // read the mass properties
  fscanf(fp, "%lf", &masstranInstance->area.vals.real);
  fscanf(fp, "%lf", &masstranInstance->mass.vals.real);

  fscanf(fp, "%lf", &masstranInstance->C.vals.reals[0]);
  fscanf(fp, "%lf", &masstranInstance->C.vals.reals[1]);
  fscanf(fp, "%lf", &masstranInstance->C.vals.reals[2]);

  fscanf(fp, "%lf", &masstranInstance->CG.vals.reals[0]);
  fscanf(fp, "%lf", &masstranInstance->CG.vals.reals[1]);
  fscanf(fp, "%lf", &masstranInstance->CG.vals.reals[2]);

  fscanf(fp, "%lf", &masstranInstance->Ixx.vals.real);
  fscanf(fp, "%lf", &masstranInstance->Iyy.vals.real);
  fscanf(fp, "%lf", &masstranInstance->Izz.vals.real);
  fscanf(fp, "%lf", &masstranInstance->Ixy.vals.real);
  fscanf(fp, "%lf", &masstranInstance->Ixz.vals.real);
  fscanf(fp, "%lf", &masstranInstance->Iyz.vals.real);

  fclose(fp); fp = NULL;

  if (inputs[inDesign_Variable-1].nullVal == NotNull) {

    // maps from nodeID to mesh->node index
    mesh_nodeID2Array(&masstranInstance->feaProblem.feaMesh, &n2a);
    AIM_NOTNULL(n2a, aimInfo, status);

    numNode = masstranInstance->feaProblem.feaMesh.numNode;

    /* copy over the mesh coordinates */
    sensNodes = new sensNode<SurrealS<1>>[numNode];

    numNode = 0; j = 0;
    for (ibody = 0; ibody < masstranInstance->numMesh; ibody++) {
      for (i = 0; i < masstranInstance->feaMesh[ibody].numNode; i++, j++) {
        sensNodes[j].xyz[0] = masstranInstance->feaMesh[ibody].node[i].xyz[0];
        sensNodes[j].xyz[1] = masstranInstance->feaMesh[ibody].node[i].xyz[1];
        sensNodes[j].xyz[2] = masstranInstance->feaMesh[ibody].node[i].xyz[2];
      }
      numNode += masstranInstance->feaMesh[ibody].numNode;
    }

    status = getComponents(aimInfo, masstranInstance, &feaProperty, &feaMaterial);
    AIM_STATUS(aimInfo, status);

    /* set derivatives for AnalysisIn variables */
    for (idv = 0; idv < masstranInstance->feaProblem.numDesignVariable; idv++) {

      designVariable = &masstranInstance->feaProblem.feaDesignVariable[idv];

      for (idr = 0; idr < designVariable->numRelation; idr++) {

        designRelation = designVariable->relationSet[idr];

        for (i = 0; i < designRelation->numDesignVariable; i++) {

          if (designVariable == designRelation->designVariableSet[i]) {

            for (ip = 0; ip < designRelation->numPropertyID; ip++) {
              for (j = 0; j < masstranInstance->feaProblem.numProperty; j++) {
                if (designRelation->propertySetID[ip] == feaProperty[j].propertyID) {
                  if (strncasecmp(designRelation->fieldName, "membraneThickness", 17) == 0) {
                    feaProperty[j].membraneThickness.deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "massPerArea", 11) == 0) {
                    feaProperty[j].massPerArea.deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "mass", 4) == 0) {
                    feaProperty[j].mass.deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "massOffset1", 11) == 0) {
                    feaProperty[j].massOffset[0].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "massOffset2", 11) == 0) {
                    feaProperty[j].massOffset[1].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "massOffset3", 11) == 0) {
                    feaProperty[j].massOffset[2].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "Ixx", 3) == 0) {
                    feaProperty[j].massInertia[0].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "Iyy", 3) == 0) {
                    feaProperty[j].massInertia[1].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "Izz", 3) == 0) {
                    feaProperty[j].massInertia[2].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "Ixy", 3) == 0) {
                    feaProperty[j].massInertia[3].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "Ixz", 3) == 0) {
                    feaProperty[j].massInertia[4].deriv() = designRelation->linearRelationCoeff[i];
                  } else if (strncasecmp(designRelation->fieldName, "Iyz", 3) == 0) {
                    feaProperty[j].massInertia[5].deriv() = designRelation->linearRelationCoeff[i];
                  } else {
                    AIM_ERROR(aimInfo, "Unknown Design_Variable_Relation fieldName = %s", designRelation->fieldName);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                  }

                  break;
                }
              }
            }

            for (im = 0; im < designRelation->numMaterialID; im++) {
              for (j = 0; j < masstranInstance->feaProblem.numMaterial; j++) {
                if (designRelation->materialSetID[im] == feaMaterial[j].materialID) {
                  if (strncasecmp(designRelation->fieldName, "density", 4) == 0) {
                    feaMaterial[j].density.deriv() = designRelation->linearRelationCoeff[i];
                  } else {
                    AIM_ERROR(aimInfo, "Unknown Design_Variable_Relation fieldName = %s", designRelation->fieldName);
                    status = CAPS_BADVALUE;
                    goto cleanup;
                  }
                  break;
                }
              }
            }
          }
        }
      }

      status = calcMassProps(aimInfo,
                             masstranInstance->feaProblem.numProperty,
                             feaProperty,
                             masstranInstance->feaProblem.numMaterial,
                             feaMaterial,
                             n2a,
                             &masstranInstance->feaProblem.feaMesh,
                             masstranInstance->Lscale,
                             sensNodes,
                             massPropS);
      AIM_STATUS(aimInfo, status);


      masstranInstance->area.derivs[idv].deriv[0] = massPropS.area.deriv();
      masstranInstance->mass.derivs[idv].deriv[0] = massPropS.mass.deriv();
      masstranInstance->Ixx .derivs[idv].deriv[0] = massPropS.Ixx .deriv();
      masstranInstance->Iyy .derivs[idv].deriv[0] = massPropS.Iyy .deriv();
      masstranInstance->Izz .derivs[idv].deriv[0] = massPropS.Izz .deriv();
      masstranInstance->Ixy .derivs[idv].deriv[0] = massPropS.Ixy .deriv();
      masstranInstance->Ixz .derivs[idv].deriv[0] = massPropS.Ixz .deriv();
      masstranInstance->Iyz .derivs[idv].deriv[0] = massPropS.Iyz .deriv();

      masstranInstance->C.derivs[idv].deriv[0] = massPropS.Cx.deriv();
      masstranInstance->C.derivs[idv].deriv[1] = massPropS.Cy.deriv();
      masstranInstance->C.derivs[idv].deriv[2] = massPropS.Cz.deriv();

      masstranInstance->CG.derivs[idv].deriv[0] = massPropS.CGx.deriv();
      masstranInstance->CG.derivs[idv].deriv[1] = massPropS.CGy.deriv();
      masstranInstance->CG.derivs[idv].deriv[2] = massPropS.CGz.deriv();

      for (j = 0; j < masstranInstance->feaProblem.numProperty; j++) {
        feaProperty[j].membraneThickness.deriv() = 0;
        feaProperty[j].massPerArea.deriv() = 0;
        feaProperty[j].mass.deriv() = 0;
        for (d = 0; d < 3; d++) {
          feaProperty[j].massOffset[d].deriv() = 0;
        }
        for (d = 0; d < 6; d++) {
          feaProperty[j].massInertia[d].deriv() = 0;
        }
      }
      for (j = 0; j < masstranInstance->feaProblem.numMaterial; j++) {
        feaMaterial[j].density.deriv() = 0;
      }
    }

    /* set derivatives for GeomegtryIn variables */
    for (idv = 0; idv < masstranInstance->feaProblem.numDesignVariable; idv++) {

      name = masstranInstance->feaProblem.feaDesignVariable[idv].name;

      // Loop over the geometry in values and compute sensitivities for all bodies
      index = aim_getIndex(aimInfo, name, GEOMETRYIN);
      status = aim_getValue(aimInfo, index, GEOMETRYIN, &geomInVal);
      if (status == CAPS_BADINDEX) continue;
      AIM_STATUS(aimInfo, status);

      if(aim_getGeomInType(aimInfo, index) != 0) {
          AIM_ERROR(aimInfo, "GeometryIn value %s is a configuration parameter and not a valid design parameter - can't get sensitivity\n",
                    name);
          status = CAPS_BADVALUE;
          goto cleanup;
      }

      for (irow = 0; irow < geomInVal->nrow; irow++) {
        for (icol = 0; icol < geomInVal->ncol; icol++) {

          j = 1; // 1-based
          for (ibody = 0; ibody < masstranInstance->numMesh; ibody++) {
            status = aim_tessSensitivity(aimInfo,
                                         name,
                                         irow+1, icol+1, // row, col
                                         masstranInstance->feaMesh[ibody].egadsTess,
                                         &i, &dxyz);
            AIM_STATUS(aimInfo, status, "Sensitivity for: %s\n", name);
            AIM_NOTNULL(dxyz, aimInfo, status);


            for (i = 0; i < masstranInstance->feaMesh[ibody].numNode; i++, j++) {
              if (n2a[j] < 0) continue;
              sensNodes[n2a[j]].xyz[0].deriv() = dxyz[3*i+0];  // dx/dGeomIn
              sensNodes[n2a[j]].xyz[1].deriv() = dxyz[3*i+1];  // dy/dGeomIn
              sensNodes[n2a[j]].xyz[2].deriv() = dxyz[3*i+2];  // dz/dGeomIn
            }
            AIM_FREE(dxyz);
          }

          status = calcMassProps(aimInfo,
                                 masstranInstance->feaProblem.numProperty,
                                 feaProperty,
                                 masstranInstance->feaProblem.numMaterial,
                                 feaMaterial,
                                 n2a,
                                 &masstranInstance->feaProblem.feaMesh,
                                 masstranInstance->Lscale,
                                 sensNodes,
                                 massPropS);
          AIM_STATUS(aimInfo, status);

          masstranInstance->area.derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.area.deriv();
          masstranInstance->mass.derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.mass.deriv();
          masstranInstance->Ixx .derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.Ixx .deriv();
          masstranInstance->Iyy .derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.Iyy .deriv();
          masstranInstance->Izz .derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.Izz .deriv();
          masstranInstance->Ixy .derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.Ixy .deriv();
          masstranInstance->Ixz .derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.Ixz .deriv();
          masstranInstance->Iyz .derivs[idv].deriv[geomInVal->ncol*irow + icol] = massPropS.Iyz .deriv();

          masstranInstance->C.derivs[idv].deriv[0*geomInVal->length + geomInVal->ncol*irow + icol] = massPropS.Cx.deriv();
          masstranInstance->C.derivs[idv].deriv[1*geomInVal->length + geomInVal->ncol*irow + icol] = massPropS.Cy.deriv();
          masstranInstance->C.derivs[idv].deriv[2*geomInVal->length + geomInVal->ncol*irow + icol] = massPropS.Cz.deriv();

          masstranInstance->CG.derivs[idv].deriv[0*geomInVal->length + geomInVal->ncol*irow + icol] = massPropS.CGx.deriv();
          masstranInstance->CG.derivs[idv].deriv[1*geomInVal->length + geomInVal->ncol*irow + icol] = massPropS.CGy.deriv();
          masstranInstance->CG.derivs[idv].deriv[2*geomInVal->length + geomInVal->ncol*irow + icol] = massPropS.CGz.deriv();
        }
      }
    }
  }

  status = CAPS_SUCCESS;
cleanup:

  if (fp != NULL) fclose(fp);

  AIM_FREE(n2a);
  AIM_FREE(dxyz);

  delete [] sensNodes;
  delete [] feaProperty;
  delete [] feaMaterial;

  return status;
}


// ********************** AIM Function Break *****************************
extern "C" int
aimOutputs( void *instStore, void *aimInfo,
            int index, char **aoname, capsValue *form)
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
   * - <B>MassProp</B> = JSON String Mass Properties
   * - <B>MassPropLink</B> = Mass Properties
   * Mass properties that can be linked to analysis input MassPropLink
   */

  if (index == outArea) {
      *aoname = EG_strdup("Area");
      form->dim  = Scalar;

  } else if (index == outMass) {
      *aoname = EG_strdup("Mass");
      form->dim  = Scalar;

  } else if (index == outCentroid) {
      *aoname = EG_strdup("Centroid");
      form->nrow = 3;
      form->dim  = Vector;

  } else if (index == outCG) {
      *aoname = EG_strdup("CG");
      form->nrow = 3;
      form->dim  = Vector;

  } else if (index == outIxx) {
      *aoname = EG_strdup("Ixx");
      form->dim  = Scalar;

  } else if (index == outIyy) {
      *aoname = EG_strdup("Iyy");
      form->dim  = Scalar;

  } else if (index == outIzz) {
      *aoname = EG_strdup("Izz");
      form->dim  = Scalar;

  } else if (index == outIxy) {
      *aoname = EG_strdup("Ixy");
      form->dim  = Scalar;

  } else if (index == outIxz) {
      *aoname = EG_strdup("Ixz");
      form->dim  = Scalar;

  } else if (index == outIyz) {
      *aoname = EG_strdup("Iyz");
      form->dim  = Scalar;

  } else if (index == outI_Vector) {
      *aoname = EG_strdup("I_Vector");
      form->nrow = 6;
      form->dim  = Vector;

  } else if (index == outI_Lower) {
      *aoname = EG_strdup("I_Lower");
      form->nrow = 6;
      form->dim  = Vector;

  } else if (index == outI_Upper) {
      *aoname = EG_strdup("I_Upper");
      form->nrow = 6;
      form->dim  = Vector;

  } else if (index == outI_Tensor) {
      *aoname = EG_strdup("I_Tensor");
      form->nrow = 9;
      form->dim  = Array2D;

  } else if (index == outMassProp) {
      *aoname = EG_strdup("MassProp");
      form->type = String;
      form->nullVal = IsNull;
      goto cleanup;

  } else if (index == outMassPropLink) {
      *aoname = EG_strdup("MassPropLink");
      form->type = Pointer;
      form->nullVal = IsNull;
      AIM_STRDUP(form->units, "feaMassPropStruct", aimInfo, status);
      goto cleanup;
  }

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

  if (units->momentOfInertia != NULL &&
      index >= outI_Vector && index <= outI_Tensor) {
    AIM_STRDUP(form->units, units->momentOfInertia, aimInfo, status);
  }

  if (index == outI_Tensor) {
      form->nrow     = 3;
      form->ncol     = 3;
  }

cleanup:
  AIM_FREE(tmpUnits);

  return status;
}


/* aimCalcOutput: Calculate/Retrieve Output Information */
extern "C" int
aimCalcOutput(void *instStore, /*@unused@*/ void *aimInfo,
              /*@unused@*/ int index, capsValue *val)
{
  int status = CAPS_SUCCESS;
  int idv, j, len_wrt;
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

  if (index < outI_Vector)
    aim_freeValue(val);

  if (index == outArea) {
    status = aim_copyValue(&masstranInstance->area, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outMass) {
    status = aim_copyValue(&masstranInstance->mass, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outCentroid) {
    status = aim_copyValue(&masstranInstance->C, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outCG) {
    status = aim_copyValue(&masstranInstance->CG, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outIxx) {
    status = aim_copyValue(&masstranInstance->Ixx, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outIyy) {
    status = aim_copyValue(&masstranInstance->Iyy, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outIzz) {
    status = aim_copyValue(&masstranInstance->Izz, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outIxy) {
    status = aim_copyValue(&masstranInstance->Ixy, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outIxz) {
    status = aim_copyValue(&masstranInstance->Ixz, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outIyz) {
    status = aim_copyValue(&masstranInstance->Iyz, val);
    AIM_STATUS(aimInfo, status);

  } else if (index == outI_Vector) {
    val->type = Double;
    AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

    val->vals.reals[0] = masstranInstance->Ixx.vals.real;
    val->vals.reals[1] = masstranInstance->Iyy.vals.real;
    val->vals.reals[2] = masstranInstance->Izz.vals.real;
    val->vals.reals[3] = masstranInstance->Ixy.vals.real;
    val->vals.reals[4] = masstranInstance->Ixz.vals.real;
    val->vals.reals[5] = masstranInstance->Iyz.vals.real;

    if (masstranInstance->Ixx.nderiv > 0) {
      val->type = DoubleDeriv;

      /* allocate derivatives */
      AIM_ALLOC(val->derivs, masstranInstance->Ixx.nderiv, capsDeriv, aimInfo, status);
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {
        val->derivs[idv].name  = NULL;
        val->derivs[idv].deriv = NULL;
        val->derivs[idv].len_wrt = 0;
      }
      val->nderiv = masstranInstance->Ixx.nderiv;

      /* set derivatives */
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {

        AIM_STRDUP(val->derivs[idv].name, masstranInstance->Ixx.derivs[idv].name, aimInfo, status);

        val->derivs[idv].len_wrt = masstranInstance->Ixx.derivs[idv].len_wrt;

        len_wrt = val->derivs[idv].len_wrt;
        AIM_ALLOC(val->derivs[idv].deriv, val->length*len_wrt, double, aimInfo, status);
        for (j = 0; j < len_wrt; j++) {
          val->derivs[idv].deriv[0*len_wrt + j] = masstranInstance->Ixx.derivs[idv].deriv[j];
          val->derivs[idv].deriv[1*len_wrt + j] = masstranInstance->Iyy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[2*len_wrt + j] = masstranInstance->Izz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[3*len_wrt + j] = masstranInstance->Ixy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[4*len_wrt + j] = masstranInstance->Ixz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[5*len_wrt + j] = masstranInstance->Iyz.derivs[idv].deriv[j];
        }
      }
    }

  } else if (index == outI_Lower) {
    val->type = Double;
    AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

    val->vals.reals[0] =  masstranInstance->Ixx.vals.real;
    val->vals.reals[1] = -masstranInstance->Ixy.vals.real;
    val->vals.reals[2] =  masstranInstance->Iyy.vals.real;
    val->vals.reals[3] = -masstranInstance->Ixz.vals.real;
    val->vals.reals[4] = -masstranInstance->Iyz.vals.real;
    val->vals.reals[5] =  masstranInstance->Izz.vals.real;

    if (masstranInstance->Ixx.nderiv > 0) {
      val->type = DoubleDeriv;

      /* allocate derivatives */
      AIM_ALLOC(val->derivs, masstranInstance->Ixx.nderiv, capsDeriv, aimInfo, status);
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {
        val->derivs[idv].name  = NULL;
        val->derivs[idv].deriv = NULL;
        val->derivs[idv].len_wrt = 0;
      }
      val->nderiv = masstranInstance->Ixx.nderiv;

      /* set derivatives */
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {

        AIM_STRDUP(val->derivs[idv].name, masstranInstance->Ixx.derivs[idv].name, aimInfo, status);

        val->derivs[idv].len_wrt = masstranInstance->Ixx.derivs[idv].len_wrt;

        len_wrt = val->derivs[idv].len_wrt;
        AIM_ALLOC(val->derivs[idv].deriv, val->length*len_wrt, double, aimInfo, status);
        for (j = 0; j < len_wrt; j++) {
          val->derivs[idv].deriv[0*len_wrt + j] =  masstranInstance->Ixx.derivs[idv].deriv[j];
          val->derivs[idv].deriv[1*len_wrt + j] = -masstranInstance->Ixy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[2*len_wrt + j] =  masstranInstance->Iyy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[3*len_wrt + j] = -masstranInstance->Ixz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[4*len_wrt + j] = -masstranInstance->Iyz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[5*len_wrt + j] =  masstranInstance->Izz.derivs[idv].deriv[j];
        }
      }
    }

  } else if (index == outI_Upper) {
    val->type = Double;
    AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

    val->vals.reals[0] =  masstranInstance->Ixx.vals.real;
    val->vals.reals[1] = -masstranInstance->Ixy.vals.real;
    val->vals.reals[2] = -masstranInstance->Ixz.vals.real;
    val->vals.reals[3] =  masstranInstance->Iyy.vals.real;
    val->vals.reals[4] = -masstranInstance->Iyz.vals.real;
    val->vals.reals[5] =  masstranInstance->Izz.vals.real;

    if (masstranInstance->Ixx.nderiv > 0) {
      val->type = DoubleDeriv;

      /* allocate derivatives */
      AIM_ALLOC(val->derivs, masstranInstance->Ixx.nderiv, capsDeriv, aimInfo, status);
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {
        val->derivs[idv].name  = NULL;
        val->derivs[idv].deriv = NULL;
        val->derivs[idv].len_wrt = 0;
      }
      val->nderiv = masstranInstance->Ixx.nderiv;

      /* set derivatives */
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {

        AIM_STRDUP(val->derivs[idv].name, masstranInstance->Ixx.derivs[idv].name, aimInfo, status);

        val->derivs[idv].len_wrt = masstranInstance->Ixx.derivs[idv].len_wrt;

        len_wrt = val->derivs[idv].len_wrt;
        AIM_ALLOC(val->derivs[idv].deriv, val->length*len_wrt, double, aimInfo, status);
        for (j = 0; j < len_wrt; j++) {
          val->derivs[idv].deriv[0*len_wrt + j] =  masstranInstance->Ixx.derivs[idv].deriv[j];
          val->derivs[idv].deriv[1*len_wrt + j] = -masstranInstance->Ixy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[2*len_wrt + j] = -masstranInstance->Ixz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[3*len_wrt + j] =  masstranInstance->Iyy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[4*len_wrt + j] = -masstranInstance->Iyz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[5*len_wrt + j] =  masstranInstance->Izz.derivs[idv].deriv[j];
        }
      }
    }

  } else if (index == outI_Tensor) {
    val->type = Double;
    AIM_ALLOC(val->vals.reals, val->length, double, aimInfo, status);

    Ixx = masstranInstance->Ixx.vals.real;
    Iyy = masstranInstance->Iyy.vals.real;
    Iyz = masstranInstance->Iyz.vals.real;
    Ixy = masstranInstance->Ixy.vals.real;
    Ixz = masstranInstance->Ixz.vals.real;
    Izz = masstranInstance->Izz.vals.real;
    I = val->vals.reals;

    // populate the tensor
    I[0] =  Ixx; I[1] = -Ixy; I[2] = -Ixz;
    I[3] = -Ixy; I[4] =  Iyy; I[5] = -Iyz;
    I[6] = -Ixz; I[7] = -Iyz; I[8] =  Izz;

    if (masstranInstance->Ixx.nderiv > 0) {
      val->type = DoubleDeriv;

      /* allocate derivatives */
      AIM_ALLOC(val->derivs, masstranInstance->Ixx.nderiv, capsDeriv, aimInfo, status);
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {
        val->derivs[idv].name  = NULL;
        val->derivs[idv].deriv = NULL;
        val->derivs[idv].len_wrt = 0;
      }
      val->nderiv = masstranInstance->Ixx.nderiv;

      /* set derivatives */
      for (idv = 0; idv < masstranInstance->Ixx.nderiv; idv++) {

        AIM_STRDUP(val->derivs[idv].name, masstranInstance->Ixx.derivs[idv].name, aimInfo, status);

        val->derivs[idv].len_wrt = masstranInstance->Ixx.derivs[idv].len_wrt;

        len_wrt = val->derivs[idv].len_wrt;
        AIM_ALLOC(val->derivs[idv].deriv, val->length*len_wrt, double, aimInfo, status);
        for (j = 0; j < len_wrt; j++) {
          val->derivs[idv].deriv[0*len_wrt + j] =  masstranInstance->Ixx.derivs[idv].deriv[j];
          val->derivs[idv].deriv[1*len_wrt + j] = -masstranInstance->Ixy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[2*len_wrt + j] = -masstranInstance->Ixz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[3*len_wrt + j] = -masstranInstance->Ixy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[4*len_wrt + j] =  masstranInstance->Iyy.derivs[idv].deriv[j];
          val->derivs[idv].deriv[5*len_wrt + j] = -masstranInstance->Iyz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[6*len_wrt + j] = -masstranInstance->Ixz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[7*len_wrt + j] = -masstranInstance->Iyz.derivs[idv].deriv[j];
          val->derivs[idv].deriv[8*len_wrt + j] =  masstranInstance->Izz.derivs[idv].deriv[j];
        }
      }
    }
  } else if (index == outMassProp) {

    mass = masstranInstance->mass.vals.real;

    CGx = masstranInstance->CG.vals.reals[0];
    CGy = masstranInstance->CG.vals.reals[1];
    CGz = masstranInstance->CG.vals.reals[2];

    Ixx = masstranInstance->Ixx.vals.real;
    Iyy = masstranInstance->Iyy.vals.real;
    Iyz = masstranInstance->Iyz.vals.real;
    Ixy = masstranInstance->Ixy.vals.real;
    Ixz = masstranInstance->Ixz.vals.real;
    Izz = masstranInstance->Izz.vals.real;

    if (units->mass != NULL) {

      snprintf(massProp, 512, "{\"mass\":[%20.14le, %s], \"CG\":[[%20.14le,%20.14le,%20.14le], %s], \"massInertia\":[[%20.14le, %20.14le, %20.14le, %20.14le, %20.14le, %20.14le], %s]}",
          mass, units->mass, CGx, CGy, CGz, units->length, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, units->momentOfInertia);

    } else {

      snprintf(massProp, 512, "{\"mass\":%20.14le, \"CG\":[%20.14le,%20.14le,%20.14le], \"massInertia\":[%20.14le, %20.14le, %20.14le, %20.14le, %20.14le, %20.14le]}",
          mass, CGx, CGy, CGz, Ixx, Iyy, Izz, Ixy, Ixz, Iyz);

    }

    AIM_STRDUP(val->vals.string, massProp, aimInfo, status);
  } else if (index == outMassPropLink) {

    masstranInstance->feaMassProp.mass = masstranInstance->mass.vals.real;
    for (j = 0; j < 3; j++)
      masstranInstance->feaMassProp.CG[j] = masstranInstance->CG.vals.reals[j];

    masstranInstance->feaMassProp.massInertia[I11] = masstranInstance->Ixx.vals.real;
    masstranInstance->feaMassProp.massInertia[I22] = masstranInstance->Iyy.vals.real;
    masstranInstance->feaMassProp.massInertia[I32] = masstranInstance->Iyz.vals.real;
    masstranInstance->feaMassProp.massInertia[I21] = masstranInstance->Ixy.vals.real;
    masstranInstance->feaMassProp.massInertia[I31] = masstranInstance->Ixz.vals.real;
    masstranInstance->feaMassProp.massInertia[I33] = masstranInstance->Izz.vals.real;

    masstranInstance->feaMassProp.massUnit = units->mass;
    masstranInstance->feaMassProp.lengthUnit = units->length;
    masstranInstance->feaMassProp.momentOfInertiaUnit = units->momentOfInertia;

    val->vals.AIMptr = (void*)&masstranInstance->feaMassProp;
  }

cleanup:
  return status;
}


/* aimCleanup: Free up the AIMs storage */
extern "C" void
aimCleanup(void *instStore)
{
  aimStorage *masstranInstance;
#ifdef DEBUG
  printf(" masstranAIM/aimCleanup!\n");
#endif

  masstranInstance = (aimStorage *) instStore;
  destroy_aimStorage(masstranInstance);
  EG_free(masstranInstance);
}
