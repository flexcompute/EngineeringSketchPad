/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             make an offset FaceBody
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include "egads.h"


int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  int    i, stat;
  double info[14];
  ego    context, surface, face, body, newbody, model;

    
  /* initialize */
  stat = EG_open(&context);
  printf(" EG_open           = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;

#ifdef UNDOCUMENTED
  /* make a conical surface */
  for (i = 0; i < 12; i++) info[i] = 0.0;
  info[3]  = info[10] = 1.0;
  info[8]  = -1.0;
  info[12] = -0.463648;
  info[13] =  1.0;

  stat = EG_makeGeometry(context, SURFACE, CONICAL, NULL, NULL, info, &surface);
  printf(" EG_makeGeometry   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;

  stat = EG_getRange(surface, info, &i);
  printf(" EG_getRange       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  printf("                 U = %lf - %lf,  V = %le - %le,  per = %d\n",
         info[0], info[1], info[2], info[3], i);
  info[1] /= 2.0;
  info[2]  = 0.0;
  info[3]  = 2.0;
  
  /* make the facebody */
  stat = EG_makeFace(surface, SFORWARD, info, &face);
  printf(" EG_makeFace       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  
  stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &face, NULL,
                         &body);
  printf(" EG_makeTopology   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  EG_deleteObject(face);
  EG_deleteObject(surface);
  
  stat = EG_hollowBody(body, 0, NULL, 0.1, 0, &newbody, NULL);
  printf(" EG_hollowBody     = %d\n", stat);
  EG_deleteObject(body);
  
#else
  
  /* make a planar surface */
  for (i = 0; i < 9; i++) info[i] = 0.0;
  info[3] = info[7] = 1.0;

  stat = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, info, &surface);
  printf(" EG_makeGeometry   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  
  stat = EG_getRange(surface, info, &i);
  printf(" EG_getRange       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  info[0] = -1.0;
  info[1] =  1.0;
  info[2] = -1.0;
  info[3] =  1.0;
  
  /* make the face */
  stat = EG_makeFace(surface, SFORWARD, info, &face);
  printf(" EG_makeFace       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;

  /* hollow it out */
  info[0] = 0.2;
  info[1] = 0.1;
  stat = EG_makeFace(face, SFORWARD, info, &body);
  printf(" EG_makeFace       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  EG_deleteObject(face);
  
  /* make the FaceBody */
  stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &body, NULL,
                         &newbody);
  printf(" EG_makeTopology   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  EG_deleteObject(body);
  EG_deleteObject(surface);

#endif

  /* make a model and write it out */
  if (newbody != NULL) {
    printf(" EG_makeTopology   = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                        NULL, 1, &newbody, NULL,
                                                        &model));
    printf(" EG_saveModel      = %d\n", EG_saveModel(model, "offset.egads"));
    printf("\n");
    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  }
                                                
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
