/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             header for interference functions
 *
 *      Copyright 2020-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


#ifndef __CLOUD_H__
#define __CLOUD_H__

#include "egads.h"


/* data structures */
/* --------------- */
   
  /* the body-based cloud of vertices */

  typedef struct {
    ego      body;              /* solid body */
    ego      tess;              /* tessellation object */
    ego      xform;             /* current transform or NULL */
    double   tParams[3];        /* tessellation parameters */
    double   bbox[6];           /* bounding box */
    int      OML;               /* 1 - outer container */
    int      nFace;             /* number of Selected Faces -- 0 is all */
    int      *facep;            /* the body's Face pairs -- indices for
                                   selected Faces & Face IDs */
    int      *vBeg;             /* global Vertex start for each Face */
    int      nVert;             /* number of vertices */
    double   *xyzs;             /* the vertex coordinates (3*nVert in length) */
  } Cloud;

  /* the minimum distance for each vertex to target */

  typedef struct {
    Cloud  source;              /* source Cloud */
    Cloud  target;              /* target Cloud */
    int    *sVert;              /* source vertex index (nVert in target len) */
    double *sMin;               /* minimum distance (nVert target in len) */
    int    *tVert;              /* target vertex index (nVert source in len) */
    double *tMin;               /* minimum distance (nVert source in len) */
  } cloud2;

  typedef struct {
    Cloud  *source;             /* source Cloud */
    Cloud  *target;             /* target Cloud */
    int    *tVert;              /* target vertex index (nVert source in len) */
    double *min;                /* minimum distance (nVert source in len) */
    int    type;                /* the interference type */
    int    filled;              /* 1 - completed classify, 2 - minimize */
    ego    model;               /* the model returned from SBO */
    int    nBody;               /* number of Bodies in the Model */
    cloud2 *bodies;             /* the Cloud pair structure for each Body */
  } cloudPair;


/* API definition */
/* -------------- */

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
initializeCloud(ego    body,               /* (in)  body pointer */
                double tParams[],          /* (in)  tessellation parameters */
                int    OML,                /* (in)  outer container */
                Cloud  *cloud);            /* (out) pointer to Cloud struct */

__ProtoExt__ int
transformCloud(Cloud *cloud,               /* (in)  the Cloud to transform */
    /*@null@*/ ego   xform);               /* (in)  the transform (or NULL) */

__ProtoExt__ int
classifyClouds(Cloud     *source,          /* (in)  Source Cloud */
               Cloud     *target,          /* (in)  Target Cloud */
               cloudPair *pair);           /* (out) pointer to cloudPair */

__ProtoExt__ int
minimizeClouds(cloudPair *pair);

__ProtoExt__ int
endpointClouds(cloudPair *pair,            /* (in)  pointer to cloudPair */
               double    *dist,            /* (out) distance */
               double    xyzs[],           /* (out) coordinates on source */
               double    xyzt[]);          /* (out) coordinates on target */
  
__ProtoExt__ void
freeCloud(Cloud *cloud);

__ProtoExt__ void
freeCloudPair(cloudPair *pair);
 

#ifdef __cplusplus
}
#endif

#endif
