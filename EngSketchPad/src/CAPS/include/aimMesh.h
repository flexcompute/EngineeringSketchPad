#ifndef AIMMESH_H
#define AIMMESH_H
/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             AIM Mesh Function Prototypes
 *
 *      Copyright 2014-2024, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <egads.h>
#include "capsTypes.h"

/* for AIM Mesh handling */

typedef struct {
  ego  tess;       /* the EGADS Tessellation Objects (contains Body) */
  int  *map;       /* the mapping between Tessellation vertices and
                      mesh vertices -- tess verts in length */
} aimMeshTessMap;

typedef struct {
  char             *groupName;  /* name of group or NULL */
  int              ID;          /* Group ID */
} aimMeshBnd;

enum aimMeshType {aimUnknownMeshType,
                  aimAreaMesh,
                  aimSurfaceMesh,
                  aimVolumeMesh};

typedef struct {
  enum aimMeshType type;     /* type of mesh referenced */
  int              nmap;     /* number of EGADS Tessellation Objects */
  aimMeshTessMap  *maps;     /* the EGADS Tess Object and map to mesh verts */
  int              nbnd;     /* number of boundary groups */
  aimMeshBnd      *bnds;     /* boundary group info */
  char            *fileName; /* full path name (no extension) for grids */
  int              _delTess; /* internal use only, whether tess/body ego are deleted */
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
  int              *elements;   /* Element-to-vertex connectivity (1-based)
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

/******************** meshWriter Helper Functions ************************/

__ProtoExt__ int
  aim_deleteMeshes( void *aimInfo, const aimMeshRef *meshRef );

__ProtoExt__ int
  aim_queryMeshes( void *aimInfo, int index, enum capssType subtype, aimMeshRef *meshRef );

__ProtoExt__ int
  aim_writeMeshes( void *aimInfo, int index, enum capssType subtype, aimMesh *mesh );

__ProtoExt__ int
  aim_initMeshBnd( aimMeshBnd *meshBnd );

__ProtoExt__ int
  aim_initMeshRef( aimMeshRef *meshRef, const enum aimMeshType type );

__ProtoExt__ int
  aim_freeMeshRef( /*@null@*/ aimMeshRef *meshRef );

__ProtoExt__ int
  aim_initMeshData( aimMeshData *meshData );

__ProtoExt__ int
  aim_freeMeshData( /*@null@*/ aimMeshData *meshData );

__ProtoExt__ int
  aim_elemTopoDim( enum aimMeshElem topo );

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

__ProtoExt__ int
  aim_storeMeshRef( void *aimStruc, const aimMeshRef *meshRef,
                    /*@null@*/ const char *meshextension );

__ProtoExt__ int
  aim_loadMeshRef( void *aimStruc, aimMeshRef *meshRef );

__ProtoExt__ int
  aim_morphMeshUpdate(void *aimInfo, aimMeshRef *meshRef, int numBody, ego *bodies);

/******************** meshWriter Dynamic Interface ***********************/

int         meshWrite( void *aimInfo, aimMesh *mesh );
const char *meshExtension( void );

/*************************************************************************/

#ifdef __cplusplus
}
#endif

#endif /* AIMMESH_H */
