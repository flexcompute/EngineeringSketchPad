#ifndef AIMMESH_H
#define AIMMESH_H
/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Mesh Function Prototypes
 *
 *      Copyright 2014-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


/* for AIM Mesh handling */

typedef struct {
  ego  tess;       /* the EGADS Tessellation Objects (contains Body) */
  int  *map;       /* the mapping between Tessellation vertices and
                      mesh vertices -- tess verts in length */
} aimMeshTessMap;

typedef struct {
  int             nmap;     /* number of EGADS Tessellation Objects */
  aimMeshTessMap *maps;     /* the EGADS Tess Object and map to mesh verts */
  char           *fileName; /* full path name (no extension) for grids */
} aimMeshRef;

typedef int    (*wrDLLFunc) (void);
typedef double aimMeshCoords[3];
typedef int    aimMeshIndices[2];

enum aimMeshElem {aimUnknownElem, aimLine, aimTri, aimQuad, aimTet, aimPyramid,
                  aimPrism, aimHex};

typedef struct {
  char             *groupName;  /* name of group or NULL */
  int              ID;          /* Group ID */
  enum aimMeshElem elementTopo; /* Element topology */
  int              order;       /* order of the element (1 - Linear) */
  int              nPoint;      /* number of points defining an element */
  int              nElems;      /* number of elements in the group */
  int              *elements;   /* Element-to-vertex connectivity
                                   nElem*nPoint in length */
} aimMeshElemGroup;

typedef struct {
  int              dim;         /* Physical dimension: 2D or 3D */
  int              nVertex;     /* total number of vertices in the mesh */
  aimMeshCoords    *verts;      /* the xyz coordinates of the vertices
                                   nVertex in length */
  int              nElemGroup;  /* number of element groups */
  aimMeshElemGroup *elemGroups; /* element groups -- nElemGroup in length */
  int              nTotalElems; /* total number of elements */
  aimMeshIndices   *elemMap;    /* group,elem map in original element ordering
                                   nTotalElems in length -- can be NULL */
} aimMeshData;

typedef struct {
  aimMeshData *meshData;
  aimMeshRef  *meshRef;
} aimMesh;


#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif


__ProtoExt__ int
  aim_deleteMeshes( void *aimInfo, aimMeshRef *meshRef );

__ProtoExt__ int
  aim_queryMeshes( void *aimInfo, int index, aimMeshRef *meshRef );

__ProtoExt__ int
  aim_writeMeshes( void *aimInfo, int index, aimMesh *mesh );

__ProtoExt__ int
  aim_initMeshRef( aimMeshRef *meshRef );

__ProtoExt__ int
  aim_freeMeshRef( /*@null@*/ aimMeshRef *meshRef );

__ProtoExt__ int
  aim_initMeshData( aimMeshData *meshData );

__ProtoExt__ int
  aim_freeMeshData( /*@null@*/ aimMeshData *meshData );

__ProtoExt__ int
  aim_addMeshElemGroup( void *aimStruc, /*@null@*/ const char *groupName,
                        int ID, enum aimMeshElem elementTopo, int order,
                        int nPoint, aimMeshData *meshData );

__ProtoExt__ int
  aim_addMeshElem( void *aimStruc, int nElems, aimMeshElemGroup *elemGroup );

__ProtoExt__ int
  aim_readBinaryUgridHeader(void *aimStruc, aimMeshRef *meshRef,
                          int *nVertex, int *nTri, int *nQuad,
                          int *nTet, int *nPyramid, int *nPrism, int *nHex);

__ProtoExt__ int
  aim_readBinaryUgrid( void *aimStruc, aimMesh *mesh );

/******************** meshWriter Dynamic Interface ***********************/

int         meshWrite( void *aimInfo, aimMesh *mesh );
const char *meshExtension( void );

/*************************************************************************/

#ifdef __cplusplus
}
#endif

#endif /* AIMMESH_H */
