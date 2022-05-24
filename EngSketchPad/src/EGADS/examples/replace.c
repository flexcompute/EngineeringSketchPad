/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Make a SolidBody and replace Faces
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

//#define INVALID

#include <math.h>
#include "egads.h"


int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  int       i, j, stat, nface, atype, len;
  double    data[6];
  ego       context, newModel, body, newBody, newFace, *faces, repl[8];
  const int *id;

  /* initialize */
  printf(" EG_open            = %d\n", EG_open(&context));
  
  /* make the box */
  data[0] = data[1] = data[2] = -1.0;
  data[3] = data[4] = data[5] =  2.0;
  stat = EG_makeSolidBody(context, BOX, data, &body);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody = %d\n", stat);
    return 1;
  }
  stat = EG_getBodyTopos(body, NULL, FACE, &nface, &faces);
  printf(" EG_getBodyTopos    = %d\n", stat);
  
  /* attribute the Faces */
  for (i = 0; i < nface; i++) {
    j = i+1;
    EG_attributeAdd(faces[i], "Face#", ATTRINT, 1, &j, NULL, NULL);
  }
  
  /* make a new Face */
  data[0] = 0.30;
  data[1] = 0.05;
  stat    = EG_makeFace(faces[4], 0, data, &newFace);
  printf(" EG_makeFace        = %d  %d\n", stat, faces[1]->mtype);
  j = -5;
  EG_attributeAdd(newFace, "Face#", ATTRINT, 1, &j, NULL, NULL);

  repl[0] = faces[3];
  repl[1] = NULL;
#ifdef INVALID
  repl[2] = faces[4];
  repl[3] = NULL;
  repl[4] = faces[5];
  repl[5] = NULL;
  repl[6] = faces[2];
  repl[7] = NULL;
  stat    = EG_replaceFaces(body, 4, repl, &newBody);
#else
  repl[2] = faces[4];
  repl[3] = newFace;
  stat    = EG_replaceFaces(body, 2, repl, &newBody);
#endif
  printf(" EG_replaceFaces    = %d\n", stat);
  EG_free(faces);
  EG_deleteObject(body);
  EG_deleteObject(newFace);
  
  /* look at the attributes */
  stat = EG_getBodyTopos(newBody, NULL, FACE, &nface, &faces);
  printf(" EG_getBodyTopos    = %d\n", stat);
  for (i = 0; i < nface; i++) {
    stat = EG_attributeRet(faces[i], "Face#", &atype, &len, &id, NULL, NULL);
    if (stat != EGADS_SUCCESS) continue;
    printf("  Face %d/%d:  Old ID = %d\n", i+1, nface, id[0]);
  }
  EG_free(faces);

  /* make a model and write it out */
  if (newBody != NULL) {
    printf(" EG_makeTopology   = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                        NULL, 1, &newBody, NULL,
                                                        &newModel));
    printf(" EG_saveModel      = %d\n", EG_saveModel(newModel, "replace.egads"));
    printf("\n");
    printf(" EG_deleteObject   = %d\n", EG_deleteObject(newModel));
  }

  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
