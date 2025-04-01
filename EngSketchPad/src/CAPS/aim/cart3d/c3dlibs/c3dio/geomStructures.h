
/**
 * Support Libraries for Cart3D I/O Functions and Extensible Design
 * Description Markup
 * ================================================================
 *
 *
 * COPYRIGHT
 *
 * Copyright © 2022 United States Government as represented by the
 * Administrator of the National Aeronautics and Space Administration.
 * All Rights Reserved.
 *
 *
 * DISCLAIMERS
 *
 * No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
 * WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
 * INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT
 * SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
 * INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
 * FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM
 * TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
 * CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR
 * RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE
 * PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE
 * SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL
 * WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
 * PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
 *
 * Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
 * AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS,
 * DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY
 * DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE
 * OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD
 * HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
 * SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT
 * PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER
 * SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
 */

/*
 * $Id: geomStructures.h,v 1.14 2022/11/07 18:57:37 mnemec Exp $
 */

/* open source */

#ifndef __GEOMSTRUCTURES_H_
#define __GEOMSTRUCTURES_H_

#include "int64.h"
#include "basicTypes.h"
#include "stateVector.h"
#include "geomTypes.h"

/* Types to support extended triangulations: we re-use Visualization Toolkit's
 * VTK types, see http://www.vtk.org/VTK/img/file-formats.pdf
 *
 * 2014.05  added VTK_unset for initialization when no match found
 */
typedef enum{VTK_Int8,    VTK_UInt8,   VTK_Int16, VTK_UInt16,
  VTK_Int32,   VTK_UInt32,  VTK_Int64, VTK_UInt64,
  VTK_Float32, VTK_Float64, VTK_unset} VTKtype;

/* Types to identify the category of extended data-sets: flow solution, shape
 * linearization and component tags
 * 2014.05  these need to stay in sync with strings in c3d_vtk_trix.c
 */
typedef enum{TRIX_unset, TRIX_flowVariable, TRIX_shapeLinearization,
  TRIX_componentTag,  TRIX_other} TRIXtype;

/* Tri and vert extended meta data, e.g. sensitivities, solution, and
 * components
 */
typedef struct surfTriExtData tsSurfTriX;
typedef tsSurfTriX * p_tsSurfTriX;
struct surfTriExtData {
  char     name[STRING_LEN];
  int      dim;    /* vector dimension */
  int      offset; /* offset into a_scalar0 or a_scalar0_t in p_surf */
  VTKtype  type;   /* datatype for VTK output */
  TRIXtype info;   /* trix type categories to help select data-sets */
};

#define DP_VERTS_CODE 8 /* trigger allocation of double precision verts */

/* --- an Indexed Triangulation ----*/
typedef struct IndexedTriangulationStructure tsTriangulation;
typedef tsTriangulation * p_tsTriangulation;

struct IndexedTriangulationStructure {
  char         geomName[FILENAME_LEN];
  int          nVerts, nTris;          /* 080423: replaced nScalars         */
  int          nVertScalars;           /* ...size of a_scalar0              */
  int          nTriScalars;            /* ...size of a_scalar0_t            */
  unsigned int infoCode;               /* ...Codes for tsTriangulation      */
  /* defined for annotations: CP_CODE 1 (0001), MACH_CODE 2 (0010),
   * STATEVEC_CODE 4 geom(0100), DP_VERTS_CODE 8 (1000)
   */
  double       configBBox[DIM*2];      /* ...bounding box of the whole thing*/
  p_tsTri      a_Tris;                 /* ...Array of triangles             */
  p_tsVertex   a_Verts;                /* ...Array of Vertices              */
  p_tsDPVertex a_dpVerts;              /* ...Double precision verts         */
  double      *a_scalar0;              /* ...place to put data (node based) */
  double      *a_scalar1;              /* ...place to put data (node based) */
  double      *a_scalar2;              /* ...place to put data (node based) */
  p_tsState    a_U;                    /* ...place to put data (node based) */
  double      *a_scalar0_t;            /* ...place to put data (tri  based) */
  double      *a_scalar1_t;            /* ...place to put data (tri  based) */
  p_tsState    a_U_t;                  /* ...place to put data (tri  based) */
  int          nVertData, nTriData;    /* ...extended tri data size         */
  p_tsSurfTriX p_vertData;             /* ...meta info for a_scalar0        */
  p_tsSurfTriX p_triData;              /* ...meta info for a_scalar0_t      */
  double      *a_area;                 /* ...areas of the triangles         */
  bool         outwardNormals;         /* ...T if CCW norm points into flow */
};                                     /*    (will be false if meshInternal)*/

/* -- public prototypes -- */

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * -- c3d_newTriangulation(): Create or resize, one or more
   *                            triangulation structs.  When creating a
   *                            single tri struct, set startComp=0 and
   *                            nComps=1, when creating an array set
   *                            startComp=0 and nComps to number of
   *                            elements, and when resizing the array,
   *                            set an startComp to the last original
   *                            element and nComps to the new total
   *                            size.
   */
  int c3d_newTriangulation(p_tsTriangulation *pp_surf, const int startComp,
                           const int nComps);

  /**
   * -- c3d_allocTriangulation(): allocate or re-size space for verts,
   *                              tris and any data associated with the
   *                              triangulation.  The function allocs
   *                              space based on the internal nVerts,
   *                              nTris, nVertData and nTriData
   *                              values.
   */
  int c3d_allocTriangulation(p_tsTriangulation *pp_surf);

  /**
   * -- c3d_allocTriData(): allocate triData info array 
   */
  int c3d_allocTriData(p_tsTriangulation *pp_surf, const int nTriData);

  /**
   * -- c3d_allocVertData(): allocate vertData info array
   */
  int c3d_allocVertData(p_tsTriangulation *pp_surf, const int nVertData);

  int c3d_resizeScalars(p_tsTriangulation *pp_surf);

  /**
   * -- c3d_freeTriangulation(): destructor for c3d_allocTriangulation()
   */
  void c3d_freeTriangulation(p_tsTriangulation p_surf, const int verbose);

  int resizeTriangulation(p_tsTriangulation p_surf, const int *const p_nVerts,
                          const int *const p_nTris);
  p_tsTriangulation deepCopyTriangulation(const p_tsTriangulation p_surf);

  int c3d_readTrixHeader(const char *const p_name, int *const p_nVerts,
                         int *const p_nTris, int *const p_nVertScalars,
                         int *const p_nTriScalars, const bool verbose);
  int c3d_readTrixGeom(p_tsTriangulation p_surf);
  int c3d_writeTrix(FILE *p_strm, const char *const p_comment, const p_tsTriangulation p_surf,
                    const bool isTrixFile);

#ifdef __cplusplus
}
#endif
  
typedef struct NamedState tsNamState;  /* -- a Named State ------------------*/
typedef tsNamState *p_tsNamState;
struct NamedState{
  INT64    name;
  tsState  U;
  tsState  Ux;  /* extra state vector for vorticity, equilibrium sound speed */
  int      info;                       /* ...for keeping any node based info */
};

typedef struct QuadMeshStructure tsQuadMesh;
typedef tsQuadMesh *p_tsQuadMesh;
struct QuadMeshStructure{
  int          iDir;         /* ...dir normal to quad mesh (if there is one) */
  float        loc;          /* ...a reference location in iDir              */
  int          nQuads;       /* ...number of quads in the mesh               */
  int          niVerts;      /* ...number of iVerts in the mesh              */
  int         *a_foundHexes; /* ...pointer back into a tsHex array           */
  int         *a_foundState; /* ...pointer back into a state array           */
  p_iquad      a_quads;      /* ...connectvty of quad (indexed into a_iVerts */
  p_tsNamState a_iVerts;     /* ...array of int64 names and an assoc state   */
  int          nx;           /* ...number of states used in Ux array         */
};

typedef struct TrianglePolyArrayStructure tstPolys;
typedef tstPolys *p_tstPolys;

struct TrianglePolyArrayStructure{
  int      *p_tpEntry;            /* ...entry into the TrianglePolyList     */
  int      *p_tpIntTriList;       /* ...Intersected Triangle poly list      */
  p_dpoint3 p_tpCentroids;        /* ...Intersected Triangle poly centroids */
  double   *p_tpAreas;            /* ...Intersected Triangle poly areas     */
  int       nTotTriPolys;         /* ...how many of triPolys are there?     */
};

#endif /* __GEOMSTRUCTURES_H_ */
