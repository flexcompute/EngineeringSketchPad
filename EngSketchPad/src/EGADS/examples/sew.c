/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Sew 2 Bodies or sew up an IGES like input
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include "egads.h"



int main(int argc, char *argv[])
{
  int i, j, jj, k, n, mtype, oclass, nbody, nface1, nface2, *senses;
  ego context, model1, model2, newModel, geom, *bodies1, *bodies2;
  ego *faces, *faces1, *faces2;
  
  if ((argc != 2) && (argc != 3)) {
    printf("\n Usage: sew filename1 [filename2]\n\n");
    return 1;
  }

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  printf(" EG_loadModel 1    = %d\n", EG_loadModel(context, 0, argv[1], &model1));
  if (model1 == NULL) {
    printf(" EG_close          = %d\n", EG_close(context));
    return 1;
  }
  /* get all bodies */
  printf(" EG_getTopology 1  = %d\n", EG_getTopology(model1, &geom, &oclass,
                                                     &mtype, NULL, &nbody,
                                                     &bodies1, &senses));
  faces = bodies1;
  k     = nbody;
  
  if (argc != 2) {
    /* get all faces in body 1 */
    printf(" EG_getBodyTopos 1 = %d\n", EG_getBodyTopos(bodies1[0], NULL, FACE,
                                                        &nface1, &faces1));
  
    printf(" EG_loadModel 2    = %d\n", EG_loadModel(context, 0, argv[2], &model2));
    if (model2 == NULL) {
      EG_free(faces1);
      printf(" EG_deleteObject   = %d\n", EG_deleteObject(model1));
      printf(" EG_close          = %d\n", EG_close(context));
      return 1;
    }
    /* get all bodies */
    printf(" EG_getTopology 2  = %d\n", EG_getTopology(model2, &geom, &oclass,
                                                       &mtype, NULL, &nbody,
                                                       &bodies2, &senses));
    /* get all faces in body 2 */
    printf(" EG_getBodyTopos 2 = %d\n", EG_getBodyTopos(bodies2[0], NULL, FACE,
                                                        &nface2, &faces2));
    
    EG_setOutLevel(context, 2);
    /* find the match */
    printf(" EG_matchBodyFaces = %d\n", EG_matchBodyFaces(bodies1[0], bodies2[0],
                                                          0.0, &n, &senses));
    EG_setOutLevel(context, 1);
    if (n != 1) {
      printf("\n ** nMatch = %d **\n", n);
      EG_free(faces1);
      EG_free(faces2);
      printf(" EG_deleteObject   = %d\n", EG_deleteObject(model1));
      printf(" EG_deleteObject   = %d\n", EG_deleteObject(model2));
      printf(" EG_close          = %d\n", EG_close(context));
      return 1;
    }
    j  = senses[0];
    jj = senses[1];
    EG_free(senses);

    faces = (ego *) malloc((nface1+nface2)*sizeof(ego));
    if (faces == NULL) {
      printf(" ERROR: malloc on %d egos!\n", nface1+nface2-2);
      EG_free(faces1);
      EG_free(faces2);
      printf(" EG_deleteObject   = %d\n", EG_deleteObject(model1));
      printf(" EG_deleteObject   = %d\n", EG_deleteObject(model2));
      printf(" EG_close          = %d\n", EG_close(context));
      return 1;
    }
    for (k = i = 0; i < nface1; i++) {
      if (i == j-1) continue;
      faces[k] = faces1[i];
      k++;
    }
    for (i = 0; i < nface2; i++) {
      if (i == jj-1) continue;
      faces[k] = faces2[i];
      k++;
    }
  }
  
  /* merge each of them */
  printf(" EG_sewFaces       = %d\n", EG_sewFaces(k, faces, 0.0, 0, &newModel));
  if (argc != 2) {
    free(faces);
    EG_free(faces1);
    EG_free(faces2);
  }

  /* write it out */
  if (newModel != NULL) {
    printf(" EG_saveModel      = %d\n", EG_saveModel(newModel, "sew.egads"));
    printf("\n");
    printf(" EG_deleteObject   = %d\n", EG_deleteObject(newModel));
  }
  
  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model1));
  if (argc != 2) printf(" EG_deleteObject   = %d\n", EG_deleteObject(model2));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
