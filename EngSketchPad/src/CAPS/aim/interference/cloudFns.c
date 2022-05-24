/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             functions that find the interference between two bodies
 *
 *      Copyright 2020-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include "cloud.h"
#include "emp.h"


//#define DEBUG

#define MIN(A,B)         (((A) < (B)) ? (A) : (B))
#define MAX(A,B)         (((A) < (B)) ? (B) : (A))


  typedef struct {
    void      *mutex;                 /* the mutex or NULL for single thread */
    long      master;                 /* master thread ID */
    int       end;                    /* end of loop */
    int       index;                  /* current loop index */
    Cloud     *source;                /* the source cloud */
    Cloud     *target;                /* the target cloud */
    int       *tVert;                 /* the target vertex index */
    double    *min;                   /* minimum distance */
  } EMPcloud;



static void calcMinDist(void *struc)
{
  int      index, i, m;
  long     ID;
  double   dist, d, xyzs[3];
  Cloud    *source, *target;
  EMPcloud *tcloud;
  
  tcloud = (EMPcloud *) struc;
  source = tcloud->source;
  target = tcloud->target;
  
  /* get our identifier */
  ID = EMP_ThreadID();
  
  /* look for work */
  for (;;) {
    
    /* only one thread at a time here -- controlled by a mutex! */
    if (tcloud->mutex != NULL) EMP_LockSet(tcloud->mutex);
    index = tcloud->index;
    tcloud->index = index+1;
    if (tcloud->mutex != NULL) EMP_LockRelease(tcloud->mutex);
    if (index >= tcloud->end) break;

    xyzs[0] = source->xyzs[3*index  ];
    xyzs[1] = source->xyzs[3*index+1];
    xyzs[2] = source->xyzs[3*index+2];
    
    dist = 1.e200;
    for (m = i = 0; i < target->nVert; i++) {
      d = (xyzs[0]-target->xyzs[3*i  ])*(xyzs[0]-target->xyzs[3*i  ]) +
          (xyzs[1]-target->xyzs[3*i+1])*(xyzs[1]-target->xyzs[3*i+1]) +
          (xyzs[2]-target->xyzs[3*i+2])*(xyzs[2]-target->xyzs[3*i+2]);
      if (d >= dist) continue;
      dist = d;
      m    = i+1;
      if (dist == 0.0) break;
    }
    tcloud->tVert[index] = m;
    tcloud->min[index]   = sqrt(dist);
  }
  
  /* exhausted all work -- exit */
  if (ID != tcloud->master) EMP_ThreadExit();
}


static void minimizeCloud(Cloud *source, Cloud *target, int *tVert, double *min)
{
  int      i, np;
  long     start;
  void     **threads = NULL;
  EMPcloud tcloud;

  tcloud.mutex  = NULL;
  tcloud.master = EMP_ThreadID();
  tcloud.index  = 0;
  tcloud.end    = source->nVert;
  tcloud.source = source;
  tcloud.target = target;
  tcloud.tVert  = tVert;
  tcloud.min    = min;
  
  np = EMP_Init(&start);

  if (np > 1) {
    /* create the mutex to handle list synchronization */
    tcloud.mutex = EMP_LockCreate();
    if (tcloud.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(tcloud.mutex);
        np = 1;
      }
    }
  }
  
  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(calcMinDist, &tcloud);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  calcMinDist(&tcloud);
  
  /* wait for all others to return & cleanup */
  if (threads != NULL) {
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  }
  if (tcloud.mutex != NULL) EMP_LockDestroy(tcloud.mutex);
  if (threads != NULL) free(threads);

  (void) EMP_Done(&start);
}


static int annotate(Cloud *cloud, const char *name, ego *body)
{
  int i, nface, stat;
  ego cbody, *faces;

  stat = EG_copyObject(cloud->body, cloud->xform, &cbody);
  if (stat != EGADS_SUCCESS) return stat;
  
  stat = EG_getBodyTopos(cbody, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    EG_deleteObject(cbody);
    return stat;
  }
  
  for (i = 1; i <= nface; i++) {
    EG_attributeDel(faces[i-1], ".source");
    EG_attributeDel(faces[i-1], ".target");
    stat = EG_attributeAdd(faces[i-1], name, ATTRINT, 1, &i, NULL, NULL);
    if (stat != EGADS_SUCCESS) {
      printf(" annotate: %s %d = %d\n", name, i, stat);
      EG_free(faces);
      EG_deleteObject(cbody);
      return stat;
    }
  }
  EG_free(faces);
  
  *body = cbody;
  return EGADS_SUCCESS;
}


void
freeCloud(Cloud *cloud)
{
  if (cloud->facep != NULL) EG_free(cloud->facep);
  if (cloud->vBeg  != NULL) EG_free(cloud->vBeg);
  if (cloud->xyzs  != NULL) EG_free(cloud->xyzs);
  if (cloud->tess  != NULL) EG_deleteObject(cloud->tess);
  cloud->body = cloud->xform = cloud->tess   = NULL;
  cloud->OML  = cloud->nFace = cloud->nVert  = 0;
  cloud->vBeg = cloud->facep = NULL;
  cloud->xyzs = NULL;
}


static int
initCloud( ego    body,                    /* (in)  body pointer */
           double tParams[],               /* (in)  tessellation parameters */
           int    OML,                     /* (in)  outer container */
           int    nFace,                   /* (in)  number of Selected Faces
                                                    0 is all */
/*@null@*/ int    *faces,                  /* (in)  body Face pairs; indices
                                                    for selected Faces & IDs */
           Cloud  *cloud )                 /* (out) pointer to Cloud struct */
{
  int          i, j, len, n, stat, nfaces, nvert, pt, pi;
  ego          dum;
  const int    *tris, *tric, *ptype, *pindex;
  const double *xyzs, *uvs;

  cloud->body = cloud->xform = cloud->tess   = NULL;
  cloud->OML  = cloud->nFace = cloud->nVert  = 0;
  cloud->vBeg = cloud->facep = NULL;
  cloud->xyzs = NULL;
  
  stat = EG_getBodyTopos(body, NULL, FACE, &nfaces, NULL);
  if (stat != EGADS_SUCCESS) {
    printf(" initializeCloud: EG_getBodyTess = %d!\n", stat);
    return stat;
  }
  if (nFace != 0) {
    if (faces == NULL) {
      printf(" initializeCloud: nFace = %d with NULL faces!\n", nFace);
      return EGADS_NODATA;
    }
    for (i = 0; i < nFace; i++) {
      if ((faces[2*i] < 1) || (faces[2*i] > nfaces)) {
        printf(" initializeCloud: faces[%d] = %d [1-%d]!\n",
               2*i, faces[2*i], nfaces);
        return EGADS_INDEXERR;
      }
      for (j = i+1; j < nFace; j++)
        if (faces[2*i] == faces[2*j]) {
          printf(" initializeCloud: faces[%d] = faces[%d]!\n", 2*i, 2*j);
          return EGADS_INDEXERR;
        }
    }
    cloud->facep = (int *) EG_alloc(2*nFace*sizeof(int));
    if (cloud->facep == NULL) {
      printf(" initializeCloud: Malloc Error on %d Faces\n!", nFace);
      return EGADS_MALLOC;
    }
    for (i = 0; i < 2*nFace; i++) cloud->facep[i] = faces[i];
    cloud->nFace = nFace;
  }
  
  /* tessellate */
  stat = EG_makeTessBody(body, tParams, &cloud->tess);
  if ((stat != EGADS_SUCCESS) || (cloud->tess == NULL)) {
    printf(" initializeCloud: EG_makeTessBody = %d!\n", stat);
    freeCloud(cloud);
    return stat;
  }
  cloud->body       = body;
  cloud->tParams[0] = tParams[0];
  cloud->tParams[1] = tParams[1];
  cloud->tParams[2] = tParams[2];
  if ((nFace == 0) || (faces == NULL)) {
    stat = EG_statusTessBody(cloud->tess, &dum, &j, &nvert);
    if (stat != EGADS_SUCCESS) {
      printf(" initializeCloud: EG_statusTessBody = %d!\n", stat);
      freeCloud(cloud);
      return stat;
    }
  } else {
    for (nvert = i = 0; i < nFace; i++) {
      stat = EG_getTessFace(cloud->tess, faces[2*i], &len, &xyzs, &uvs, &ptype,
                            &pindex, &n, &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        printf(" initializeCloud: EG_getTessFace0 %d = %d!\n", i+1, stat);
        freeCloud(cloud);
        return stat;
      }
      nvert += len;
    }
  }
#ifdef DEBUG
  printf("    nvert = %d\n", nvert);
#endif
  if (nvert == 0) {
    printf(" initializeCloud: nVerts = 0!\n");
    freeCloud(cloud);
    return EGADS_NODATA;
  }

  cloud->xyzs = (double *) EG_alloc(3*nvert*sizeof(double));
  if (cloud->xyzs == NULL) {
    printf(" initializeCloud: Malloc Error on %d Coordinates\n!", nvert);
    freeCloud(cloud);
    return EGADS_MALLOC;
  }
  if ((nFace == 0) || (faces == NULL)) {
    for (i = 0; i < nvert; i++) {
      stat = EG_getGlobal(cloud->tess, i+1, &pt, &pi, &cloud->xyzs[3*i]);
      if (stat != EGADS_SUCCESS) {
        printf(" initializeCloud: EG_getGlobal %d = %d!\n", i+1, stat);
        freeCloud(cloud);
        return stat;
      }
    }
  } else {
    cloud->vBeg = (int *) EG_alloc(nFace*sizeof(int));
    if (cloud->vBeg == NULL) {
      printf(" initializeCloud: Malloc Error on %d Faces\n!", nFace);
      freeCloud(cloud);
      return EGADS_MALLOC;
    }
    for (nvert = i = 0; i < nFace; i++) {
      cloud->vBeg[i] = nvert;
      stat = EG_getTessFace(cloud->tess, faces[2*i], &len, &xyzs, &uvs, &ptype,
                            &pindex, &n, &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        printf(" initializeCloud: EG_getTessFace1 %d = %d!\n", i+1, stat);
        freeCloud(cloud);
        return stat;
      }
      for (j = 0; j < len; j++, nvert++) {
        cloud->xyzs[3*nvert  ] = xyzs[3*j  ];
        cloud->xyzs[3*nvert+1] = xyzs[3*j+1];
        cloud->xyzs[3*nvert+2] = xyzs[3*j+2];
      }
    }
  }
  
  cloud->bbox[0] = cloud->bbox[1] = cloud->bbox[2] =  1.e100;
  cloud->bbox[3] = cloud->bbox[4] = cloud->bbox[5] = -1.e100;
  for (i = 0; i < nvert; i++) {
    if (cloud->xyzs[3*i  ] < cloud->bbox[0]) cloud->bbox[0] = cloud->xyzs[3*i  ];
    if (cloud->xyzs[3*i  ] > cloud->bbox[3]) cloud->bbox[3] = cloud->xyzs[3*i  ];
    if (cloud->xyzs[3*i+1] < cloud->bbox[1]) cloud->bbox[1] = cloud->xyzs[3*i+1];
    if (cloud->xyzs[3*i+1] > cloud->bbox[4]) cloud->bbox[4] = cloud->xyzs[3*i+1];
    if (cloud->xyzs[3*i+2] < cloud->bbox[2]) cloud->bbox[2] = cloud->xyzs[3*i+2];
    if (cloud->xyzs[3*i+2] > cloud->bbox[5]) cloud->bbox[5] = cloud->xyzs[3*i+2];
  }
  cloud->nVert = nvert;
  
  if (OML == 1) cloud->OML = OML;
  
  return EGADS_SUCCESS;
}


int
initializeCloud(ego    body,               /* (in)  body pointer */
                double tParams[],          /* (in)  tessellation parameters */
                int    OML,                /* (in)  outer container */
                Cloud  *cloud)             /* (out) pointer to Cloud struct */
{
  return initCloud(body, tParams, OML, 0, NULL, cloud);
}


int
transformCloud(Cloud *cloud,               /* (in)  the Cloud to transform */
    /*@null@*/ ego   xform)                /* (in)  the transform (or NULL) */
{
  int    i, stat, pi, pt;
  double xyz[3], tform[12];
  
  if (cloud->nFace != 0) {
    printf(" transformCloud: Not a full Solid!\n");
    return EGADS_TESSTATE;
  }

  if (xform == NULL) {
    for (i = 1; i < 11; i++) tform[i] = 0.0;
    tform[0]  = tform[5] = tform[10]  = 1.0;
    tform[11] = 0.0;
  } else {
    stat = EG_getTransformation(xform, tform);
    if (stat != EGADS_SUCCESS) {
      printf(" transformCloud: EG_getTransform = %d\n", stat);
      return stat;
    }
  }

  cloud->bbox[0] = cloud->bbox[1] = cloud->bbox[2] =  1.e100;
  cloud->bbox[3] = cloud->bbox[4] = cloud->bbox[5] = -1.e100;
  for (i = 0; i < cloud->nVert; i++) {
    stat = EG_getGlobal(cloud->tess, i+1, &pt, &pi, xyz);
    if (stat != EGADS_SUCCESS) {
      printf(" transformCloud: EG_getGlobal %d = %d!\n", i+1, stat);
      freeCloud(cloud);
      return stat;
    }
    cloud->xyzs[3*i  ] = xyz[0]*tform[ 0] + xyz[1]*tform[ 1] +
                         xyz[2]*tform[ 2] +        tform[ 3];
    cloud->xyzs[3*i+1] = xyz[0]*tform[ 4] + xyz[1]*tform[ 5] +
                         xyz[2]*tform[ 6] +        tform[ 7];
    cloud->xyzs[3*i+2] = xyz[0]*tform[ 8] + xyz[1]*tform[ 9] +
                         xyz[2]*tform[10] +        tform[11];
    if (cloud->xyzs[3*i  ] < cloud->bbox[0]) cloud->bbox[0] = cloud->xyzs[3*i  ];
    if (cloud->xyzs[3*i  ] > cloud->bbox[3]) cloud->bbox[3] = cloud->xyzs[3*i  ];
    if (cloud->xyzs[3*i+1] < cloud->bbox[1]) cloud->bbox[1] = cloud->xyzs[3*i+1];
    if (cloud->xyzs[3*i+1] > cloud->bbox[4]) cloud->bbox[4] = cloud->xyzs[3*i+1];
    if (cloud->xyzs[3*i+2] < cloud->bbox[2]) cloud->bbox[2] = cloud->xyzs[3*i+2];
    if (cloud->xyzs[3*i+2] > cloud->bbox[5]) cloud->bbox[5] = cloud->xyzs[3*i+2];
  }
  cloud->xform = xform;
  
  return EGADS_SUCCESS;
}


void
freeCloudPair(cloudPair *pair)
{
  int i;
  
  for (i = 0; i < pair->nBody; i++) {
    if (pair->bodies == NULL) continue;
    if (pair->bodies[i].sVert != NULL) EG_free(pair->bodies[i].sVert);
    if (pair->bodies[i].sMin  != NULL) EG_free(pair->bodies[i].sMin);
    if (pair->bodies[i].tVert != NULL) EG_free(pair->bodies[i].tVert);
    if (pair->bodies[i].tMin  != NULL) EG_free(pair->bodies[i].tMin);
    if (pair->bodies[i].source.body != NULL) freeCloud(&pair->bodies[i].source);
    if (pair->bodies[i].target.body != NULL) freeCloud(&pair->bodies[i].target);
  }
  if (pair->tVert  != NULL) EG_free(pair->tVert);
  if (pair->min    != NULL) EG_free(pair->min);
  if (pair->bodies != NULL) EG_free(pair->bodies);
  if (pair->model  != NULL) EG_deleteObject(pair->model);
  pair->type   = 0;
  pair->nBody  = 0;
  pair->filled = 0;
  pair->source = NULL;
  pair->target = NULL;
  pair->tVert  = NULL;
  pair->min    = NULL;
  pair->model  = NULL;
  pair->bodies = NULL;
}


int
classifyClouds(Cloud     *source,          /* (in)  Source Cloud */
               Cloud     *target,          /* (in)  Target Cloud */
               cloudPair *pair)            /* (out) pointer to cloudPair */
{
  int          i, i0, i1, j, k, stat, nVert, outLevel;
  int          oclass, mtype, nbody, nfaces, atype, alen, *sf, *tf, *senses;
  double       tParams[3];
  ego          context, sbody, tbody, geom, *bodies, *faces;
  const char   *string;
  const int    *ints;
  const double *reals;
  
  pair->type   = 0;
  pair->nBody  = 0;
  pair->filled = 0;
  pair->source = source;
  pair->target = target;
  pair->tVert  = NULL;
  pair->min    = NULL;
  pair->model  = NULL;
  pair->bodies = NULL;
  
  if (target->OML == 1) {
    printf(" classifyClouds: target is outer!\n");
    freeCloudPair(pair);
    return EGADS_TOPOERR;
  }
  stat = EG_getContext(source->body, &context);
  if (stat != EGADS_SUCCESS) {
    printf(" classifyClouds: EG_getContext = %d!\n", stat);
    freeCloudPair(pair);
    return stat;
  }
  
  /* overlapping geometry? */
  outLevel = EG_setOutLevel(context, 0);
  stat = annotate(source, ".source", &sbody);
  if (stat != EGADS_SUCCESS) {
    printf(" classifyClouds: annotate source = %d!\n", stat);
    freeCloudPair(pair);
    return stat;
  }
  stat = annotate(target, ".target", &tbody);
  if (stat != EGADS_SUCCESS) {
    printf(" classifyClouds: annotate target = %d!\n", stat);
    freeCloudPair(pair);
    return stat;
  }
                        k = INTERSECTION;
  if (source->OML == 1) k = SUBTRACTION;
  stat = EG_generalBoolean(tbody, sbody, k, 0.0, &pair->model);
  EG_setOutLevel(context, outLevel);
  EG_deleteObject(sbody);
  EG_deleteObject(tbody);
  if ((stat != EGADS_SUCCESS) && (stat != EGADS_ATTRERR)) {
    printf(" classifyClouds: EG_generalBoolean = %d!\n", stat);
    freeCloudPair(pair);
    return stat;
  }
    
  /* set the interference type */
  if (stat == EGADS_SUCCESS) {
    pair->type = 2;
  } else {
    if ((source->bbox[0] >= target->bbox[0]) &&
        (source->bbox[3] <= target->bbox[3]) &&
        (source->bbox[1] >= target->bbox[1]) &&
        (source->bbox[4] <= target->bbox[4]) &&
        (source->bbox[2] >= target->bbox[2]) &&
        (source->bbox[5] <= target->bbox[5])) {
      pair->type = 3;
    } else if ((target->bbox[0] >= source->bbox[0]) &&
               (target->bbox[3] <= source->bbox[3]) &&
               (target->bbox[1] >= source->bbox[1]) &&
               (target->bbox[4] <= source->bbox[4]) &&
               (target->bbox[2] >= source->bbox[2]) &&
               (target->bbox[5] <= source->bbox[5])) {
      pair->type = 4;
    } else {
      pair->type = 1;
    }
  }
  
  if (pair->type != 2) {

    pair->tVert = (int *)    EG_alloc(source->nVert*sizeof(int));
    pair->min   = (double *) EG_alloc(source->nVert*sizeof(double));
    if ((pair->tVert == NULL) || (pair->min == NULL)) {
      printf(" classifyClouds: Error allocating %d tri storage!\n",
             source->nVert);
      freeCloudPair(pair);
      return EGADS_MALLOC;
    }
    for (i = 0; i < source->nVert; i++) {
      pair->tVert[i] = 0;
      pair->min[i]   = 1.e100;
    }
    
  } else {
  
/*@-nullpass@*/
    stat = EG_getTopology(pair->model, &geom, &oclass, &mtype, NULL, &nbody,
                          &bodies, &senses);
/*@+nullpass@*/
    if (stat != EGADS_SUCCESS) {
      printf(" classifyClouds: EG_getTopology = %d!\n", stat);
      freeCloudPair(pair);
      return stat;
    }
#ifdef DEBUG
    printf("    nBody = %d\n", nbody);
#endif
    pair->bodies = (cloud2 *) EG_alloc(nbody*sizeof(cloud2));
    if (pair->bodies == NULL) {
      printf(" classifyClouds: Malloc of %d Body storage!\n", nbody);
      freeCloudPair(pair);
      return stat;
    }
    for (i = 0; i < nbody; i++) {
      pair->bodies[i].source.body = NULL;
      pair->bodies[i].target.body = NULL;
      pair->bodies[i].sVert       = NULL;
      pair->bodies[i].sMin        = NULL;
      pair->bodies[i].tVert       = NULL;
      pair->bodies[i].tMin        = NULL;
    }
    pair->nBody = nbody;
    
    /* make the scaled paramneters smaller */
    tParams[0] = MIN(source->tParams[0], target->tParams[0])/2.0;
    tParams[1] = MIN(source->tParams[1], target->tParams[1])/2.0;
    tParams[2] = MIN(source->tParams[2], target->tParams[2]);
    
    for (i = 0; i < nbody; i++) {
      stat = EG_getBodyTopos(bodies[i], NULL, FACE, &nfaces, &faces);
      if (stat != EGADS_SUCCESS) {
        printf(" classifyClouds: EG_getBodyTopos %d = %d!\n", i+1, stat);
        freeCloudPair(pair);
        return stat;
      }
      for (i0 = i1 = j = 0; j < nfaces; j++) {
        stat = EG_attributeRet(faces[j], ".source", &atype, &alen,
                               &ints, &reals, &string);
        if (stat == EGADS_SUCCESS) {
#ifdef DEBUG
          printf("  Body %d Face %2d: source = %d\n", i+1, j+1, ints[0]);
#endif
          i0++;
        } else {
          stat = EG_attributeRet(faces[j], ".target", &atype, &alen,
                                 &ints, &reals, &string);
          if (stat == EGADS_SUCCESS) {
#ifdef DEBUG
            printf("  Body %d Face %2d: target = %d\n", i+1, j+1, ints[0]);
#endif
            i1++;
          } else {
            printf("  Body %d Face %2d: No marker!\n", i+1, j+1);
          }
        }
      }
      if ((i0 == 0) || (i1 == 0)) {
        printf(" classifyClouds: %d Source nFaces = %d  Target nFaces = %d!\n",
               i+1, i0, i1);
        freeCloudPair(pair);
        EG_free(faces);
        return EGADS_ATTRERR;
      }
      sf = (int *) EG_alloc(2*i0*sizeof(int));
      tf = (int *) EG_alloc(2*i1*sizeof(int));
      if ((sf == NULL) || (tf == NULL)) {
        printf(" classifyClouds: Malloc Face arrays = %d %d!\n", i0, i1);
        if (sf != NULL) EG_free(sf);
        if (tf != NULL) EG_free(tf);
        freeCloudPair(pair);
        EG_free(faces);
        return EGADS_MALLOC;
      }
      for (i0 = i1 = j = 0; j < nfaces; j++) {
        stat = EG_attributeRet(faces[j], ".source", &atype, &alen,
                               &ints, &reals, &string);
        if (stat == EGADS_SUCCESS) {
          sf[i0] = j+1;
          i0++;
          sf[i0] = ints[0];
          i0++;
        } else {
          stat = EG_attributeRet(faces[j], ".target", &atype, &alen,
                                 &ints, &reals, &string);
          if (stat == EGADS_SUCCESS) {
            tf[i1] = j+1;
            i1++;
            tf[i1] = ints[0];
            i1++;
          }
        }
      }
      EG_free(faces);

      stat = initCloud(bodies[i], tParams, 0, i0/2, sf, &pair->bodies[i].source);
      if (stat != EGADS_SUCCESS) {
        printf(" classifyClouds: initCloud %d source = %d!\n", i+1, stat);
        EG_free(sf);
        EG_free(tf);
        freeCloudPair(pair);
        EG_free(faces);
        return stat;
      }
      stat = initCloud(bodies[i], tParams, 0, i1/2, tf, &pair->bodies[i].target);
      if (stat != EGADS_SUCCESS) {
        printf(" classifyClouds: initCloud %d target = %d!\n", i+1, stat);
        EG_free(sf);
        EG_free(tf);
        freeCloudPair(pair);
        EG_free(faces);
        return stat;
      }
      EG_free(sf);
      EG_free(tf);
      
      nVert = pair->bodies[i].source.nVert;
      pair->bodies[i].tVert = (int *)    EG_alloc(nVert*sizeof(int));
      pair->bodies[i].tMin  = (double *) EG_alloc(nVert*sizeof(double));
      if ((pair->bodies[i].tVert == NULL) || (pair->bodies[i].tMin == NULL)) {
        printf(" classifyClouds: Error %d allocating %d source verts!\n",
               i+1, nVert);
        freeCloudPair(pair);
        return EGADS_MALLOC;
      }
      for (j = 0; j < nVert; j++) {
        pair->bodies[i].tVert[j] = 0;
        pair->bodies[i].tMin[j]  = 1.e100;
      }
      
      nVert = pair->bodies[i].target.nVert;
      pair->bodies[i].sVert = (int *)    EG_alloc(nVert*sizeof(int));
      pair->bodies[i].sMin  = (double *) EG_alloc(nVert*sizeof(double));
      if ((pair->bodies[i].sVert == NULL) || (pair->bodies[i].sMin == NULL)) {
        printf(" classifyClouds: Error %d allocating %d target verts!\n",
               i+1, nVert);
        freeCloudPair(pair);
        return EGADS_MALLOC;
      }
      for (j = 0; j < nVert; j++) {
        pair->bodies[i].sVert[j] = 0;
        pair->bodies[i].sMin[j]  = 1.e100;
      }
    }
  }
  
  pair->filled = 1;
         
  return EGADS_SUCCESS;
}


int
minimizeClouds(cloudPair *pair)
{
  int i;
 
  if (pair->filled == 0) {
    printf(" minimizeClouds: cloudPair not completely classified!\n");
    return EGADS_TOPOCNT;
  } else if (pair->filled == 2) {
    printf(" minimizeClouds: cloudPair already minimnized!\n");
    return EGADS_EXISTS;
  }

  if (pair->nBody == 0) {
    minimizeCloud(pair->source, pair->target, pair->tVert, pair->min);
  } else {
    for (i = 0; i < pair->nBody; i++) {
      minimizeCloud(&pair->bodies[i].source, &pair->bodies[i].target,
                     pair->bodies[i].tVert,   pair->bodies[i].tMin);
      minimizeCloud(&pair->bodies[i].target, &pair->bodies[i].source,
                     pair->bodies[i].sVert,   pair->bodies[i].sMin);
    }
  }
  
  /* mark as completed */
  pair->filled = 2;
  
  return EGADS_SUCCESS;
}


int
endpointClouds(cloudPair *pair,            /* (in)  pointer to cloudPair */
               double    *dist,            /* (out) distance */
               double    xyzs[],           /* (out) coordinates on source */
               double    xyzt[])           /* (out) coordinates on target */
{
  int    i, ibody, k, ib = 0, is = 0;
  double val = 0.0;
  Cloud  *source, *target;
  int    *Vert;
  
  if (pair->filled == 0) {
    printf(" endpointClouds: cloudPair not completely classified!\n");
    return EGADS_TOPOCNT;
  } else if (pair->filled == 1) {
    printf(" endpointClouds: cloudPair not minimnized!\n");
    return EGADS_EXISTS;
  }
  
  /* not intersecting */
  if (pair->nBody == 0) {

    source = pair->source;
    target = pair->target;
    Vert   = pair->tVert;

    if (source->OML != 1) {
      if (pair->type == 1) {
        val = pair->min[0];
        i   = 0;
        for (k = 1; k < source->nVert; k++)
          if (pair->min[k] < val) {
            val = pair->min[k];
            i   = k;
          }
      } else {
        printf(" endpointClouds Error: Type = %d!\n", pair->type);
        return EGADS_INDEXERR;
      }
    } else {
      if (pair->type == 4) {
        val = pair->min[0];
        i   = 0;
        for (k = 1; k < source->nVert; k++)
          if (pair->min[k] < val) {
            val = pair->min[k];
            i   = k;
          }
      } else {
        printf(" endpointClouds Error: OML Type = %d!\n", pair->type);
        return EGADS_INDEXERR;
      }
    }
    *dist   = val;
    
    /* source position */
    xyzs[0] = source->xyzs[3*i  ];
    xyzs[1] = source->xyzs[3*i+1];
    xyzs[2] = source->xyzs[3*i+2];
    
    /* target position */
    i       = Vert[i]-1;
    xyzt[0] = target->xyzs[3*i  ];
    xyzt[1] = target->xyzs[3*i+1];
    xyzt[2] = target->xyzs[3*i+2];
    
    return EGADS_SUCCESS;
  
  }
    
  /* intersecting */
  val = pair->bodies[0].tMin[0];
  for (ib = i = ibody = 0; ibody < pair->nBody; ibody++) {
    for (k = 0; k < pair->bodies[ibody].source.nVert; k++)
      if (pair->bodies[ibody].tMin[k] > val) {
        val = pair->bodies[ibody].tMin[k];
        ib  = ibody;
        i   = k;
        is  = 0;
      }
    for (k = 0; k < pair->bodies[ibody].target.nVert; k++)
      if (pair->bodies[ibody].sMin[k] > val) {
        val = pair->bodies[ibody].sMin[k];
        ib  = ibody;
        i   = k;
        is  = 1;
      }
  }
  *dist  = -val;
  source = &pair->bodies[ib].source;
  target = &pair->bodies[ib].target;
  
  if (is == 0) {
    Vert  =  pair->bodies[ib].tVert;
    
    /* source position */
    xyzs[0] = source->xyzs[3*i  ];
    xyzs[1] = source->xyzs[3*i+1];
    xyzs[2] = source->xyzs[3*i+2];
    
    /* target position */
    i       = Vert[i]-1;
    xyzt[0] = target->xyzs[3*i  ];
    xyzt[1] = target->xyzs[3*i+1];
    xyzt[2] = target->xyzs[3*i+2];
  } else {
    Vert  =  pair->bodies[ib].sVert;
    
    /* target position */
    xyzt[0] = target->xyzs[3*i  ];
    xyzt[1] = target->xyzs[3*i+1];
    xyzt[2] = target->xyzs[3*i+2];
    
    /* source position */
    i       = Vert[i]-1;
    xyzs[0] = source->xyzs[3*i  ];
    xyzs[1] = source->xyzs[3*i+1];
    xyzs[2] = source->xyzs[3*i+2];
  }
  
  return EGADS_SUCCESS;
}
