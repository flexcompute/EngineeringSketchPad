/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Report the Tolerances in an Model
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"


int main(int argc, char *argv[])
{
  int    i, j, n, stat, oclass, mtype, nbodies, *senses;
  double val, tol, size, bbox[6];
  ego    context, model, geom, ref, *objs, *bodies;
  
  if (argc != 2) {
    printf("\n Usage: tolerance modelFile\n\n");
    return 1;
  }
  
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" Error: EG_open = %d!\n", stat);
    return 1;
  }
  stat = EG_loadModel(context, 0, argv[1], &model);
  if (stat != EGADS_SUCCESS) {
    printf(" Error: EG_loadModel = %d!\n", stat);
    EG_close(context);
    return 1;
  }
  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" Error: EG_getTopology  = %d!\n", stat);
    EG_deleteObject(model);
    EG_close(context);
    return 1;
  }

  for (i = 0; i < nbodies; i++) {
    stat = EG_getBoundingBox(bodies[i], bbox);
    if (stat != EGADS_SUCCESS) continue;
                                  size = bbox[3] - bbox[0];
    if (size < bbox[4] - bbox[1]) size = bbox[4] - bbox[1];
    if (size < bbox[5] - bbox[2]) size = bbox[5] - bbox[2];
    printf("\n Body %d: Ref Size = %lf\n", i+1, size);
    
    stat = EG_getBodyTopos(bodies[i], NULL, NODE, &n, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" Error: EG_getBodyTopos NODE = %d!\n", stat);
      continue;
    }
    val = 0.0;
    for (j = 0; j < n; j++) {
      stat = EG_getTolerance(objs[j], &tol);
      if (stat != EGADS_SUCCESS) continue;
      if (tol > val) val = tol;
      stat = EG_getBody(objs[j], &ref);
      if (stat == EGADS_SUCCESS) {
        if (bodies[i] != ref) printf("  Wrong Body for NODE %d!\n", j+1);
      } else {
        printf(" Error: EG_getBody NODE = %d!\n", stat);
      }
    }
    EG_free(objs);
    printf("         Max Node tolerance = %lf  %le   %d\n",val, val/size, n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, EDGE, &n, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" Error: EG_getBodyTopos EDGE = %d!\n", stat);
      continue;
    }
    val = 0.0;
    for (j = 0; j < n; j++) {
      stat = EG_getTolerance(objs[j], &tol);
      if (stat != EGADS_SUCCESS) continue;
      if (tol > val) val = tol;
      stat = EG_getBody(objs[j], &ref);
      if (stat == EGADS_SUCCESS) {
        if (bodies[i] != ref) printf("  Wrong Body for EDGE %d!\n", j+1);
      } else {
        printf(" Error: EG_getBody EDGE = %d!\n", stat);
      }
    }
    EG_free(objs);
    printf("         Max Edge tolerance = %lf  %le   %d\n",val, val/size, n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, FACE, &n, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" Error: EG_getBodyTopos FACE = %d!\n", stat);
      continue;
    }
    val = 0.0;
    for (j = 0; j < n; j++) {
      stat = EG_getTolerance(objs[j], &tol);
      if (stat != EGADS_SUCCESS) continue;
      if (tol > val) val = tol;
      stat = EG_getBody(objs[j], &ref);
      if (stat == EGADS_SUCCESS) {
        if (bodies[i] != ref) printf("  Wrong Body for FACE %d!\n", j+1);
      } else {
        printf(" Error: EG_getBody FACE = %d!\n", stat);
      }
    }
    EG_free(objs);
    printf("         Max Face tolerance = %lf  %le   %d\n",val, val/size, n);
  }
  
  printf("\n");
  EG_close(context);
  return 0;
}
