#include <stdio.h>
#include <stdlib.h>

#include "egads.h"


int main(int argc, char *argv[])
{
  int          i, j, k, stat, oclass, mtype, nbody, ibody, nface, *senses;
  int          ntris, npts, per;
  const int    *tris, *tric, *ptype, *pindex;
  double       box[6], params[3], lims[4], size, tol;
  const double *xyzs, *uvs;
  const char   *OCCrev;
  ego          context, model, geom, tess, *bodies, face, surf;
  ego          *fbodies, newmodel;
  
  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d %s\n\n", i, j, OCCrev);

  if (argc != 2) {
    printf("\n Usage: fitTri modelFile\n\n");
    return 1;
  }

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  printf(" EG_loadModel      = %d\n", EG_loadModel(context, 0, argv[1], &model));
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);  
  printf(" \n");
  
                            size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];
  tol = 1.e-6*size;

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n\n", nbody);

  params[0] =  0.015*size;
  params[1] =  0.001*size;
  params[2] = 12.0;
  
  for (ibody = 0; ibody < nbody; ibody++) {
    stat = EG_makeTessBody(bodies[ibody], params, &tess);
    if (stat != EGADS_SUCCESS) continue;
    printf(" EG_getBodyTopos   = %d\n", EG_getBodyTopos(bodies[ibody], NULL,
                                                        FACE, &nface, NULL));
    fbodies = (ego *) EG_alloc(nface*sizeof(ego));
    if (fbodies == NULL) {
      printf(" Fatal Error on fbodies!\n");
      EG_deleteObject(tess);
      EG_deleteObject(model);
      EG_close(context);
      return 1;
    }
    
    EG_setOutLevel(context, 2);		/* get some Debug output */
    k = 0;
    for (i = 1; i <= nface; i++) {
      fbodies[k] = NULL;
      stat = EG_getTessFace(tess, i, &npts, &xyzs, &uvs, &ptype, &pindex,
                            &ntris, &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      printf(" Face %d:\n", i);
      stat = EG_fitTriangles(context, npts, (double *) xyzs, ntris, tris, tric,
                             tol, &surf);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_fitTriangles   = %d\n", stat);
      } else {
        stat = EG_getRange(surf, lims, &per);
        if (stat != EGADS_SUCCESS) printf("  getRange = %d!\n", stat);
        stat = EG_makeFace(surf, SFORWARD, lims, &face);
        if (stat != EGADS_SUCCESS) printf("  makeFace = %d!\n", stat);
        stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &face,
                               NULL, &fbodies[k]);
        if (stat == EGADS_SUCCESS) {
          k++;
        } else {
          printf("  makeTopology = %d!\n", stat);
        }
        EG_deleteObject(face);
        EG_deleteObject(surf);
      }
    }
    stat = EG_makeTopology(context, NULL, MODEL, 0, NULL, k, fbodies,
                           NULL, &newmodel);
    if (stat != EGADS_SUCCESS) printf("  makeTopology on model = %d!\n", stat);
    stat = EG_saveModel(newmodel, "fitTri.egads");
    if (stat != EGADS_SUCCESS) printf("  saveModel = %d!\n", stat);
    EG_deleteObject(newmodel);
    EG_free(fbodies);
    EG_setOutLevel(context, 1);
    printf(" EG_deleteObject T = %d\n\n", EG_deleteObject(tess));
  }

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));

  return 0;
}
