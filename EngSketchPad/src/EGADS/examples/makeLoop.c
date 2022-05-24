/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Test makeLoop from Edges
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"


int main(int argc, char *argv[])
{
  int i, j, n, stat, mtype, oclass, nbody, nedge, nloop, nface, *senses;
  ego context, model, geom, loop, *bodies, *edges, *faces, *loops;
  
  if (argc != 2) {
    printf("\n Usage: makeLoop modelFile\n\n");
    return 1;
  }

  /* initialize */
  printf(" EG_open          = %d\n", EG_open(&context));
  printf(" EG_loadModel     = %d\n", EG_loadModel(context, 0, argv[1], &model));
  /* get all bodies */
  printf(" EG_getTopology   = %d\n", EG_getTopology(model, &geom, &oclass, 
                                                    &mtype, NULL, &nbody,
                                                    &bodies, &senses));
  printf("\n");
  for (i = 0; i < nbody; i++) {
    stat = EG_getBodyTopos(bodies[i], NULL, FACE, &nface, &faces);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getBodyToposF = %d for Body %d\n", stat, i+1);
      continue;
    }
    for (j = 0; j < nface; j++) {
      stat = EG_getTopology(faces[j], &geom, &oclass, &mtype, NULL, &nloop,
                            &loops, &senses);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getTopology F = %d for Body %d Face %d\n", stat, i+1, j+1);
        continue;
      }
      stat = EG_getBodyTopos(bodies[i], faces[j], EDGE, &nedge, &edges);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getBodyToposE = %d for Body %d Face %d\n", stat, i+1, j+1);
        continue;
      }
      n = 0;
      do {
        stat = EG_makeLoop(nedge, edges, geom, 0.0, &loop);
        if (stat >= EGADS_SUCCESS) {
          if (loop->mtype == OPEN) printf("        Face %d: loop %d is Open\n",
                                          j+1, n+1);
          EG_deleteObject(loop);
          n++;
        }
      } while (stat > EGADS_SUCCESS);
      EG_free(edges);
      if (stat < EGADS_SUCCESS) {
        printf(" EG_makeLoop      = %d for Body %d Face %d (%d)\n",
               stat, i+1, j+1, nloop);
      } else {
        printf(" Body %d/Face %d: remade %d Loops from %d Edges\n",
               i+1, j+1, n, nedge);
      }
    }
    EG_free(faces);
  }

  printf(" \n");
  printf(" EG_close         = %d\n", EG_close(context));
  return 0;
}
