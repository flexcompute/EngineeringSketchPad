/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             fit a cubic spline to a cone
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
  int    i, mDeg, stat, sizes[2], oclass, mtype, *splInfo;
  double info[14], *xyzs;
  ego    context, cone, tess, bspline, ref, face, body, model;

  if (argc != 2) {
    printf("\n usage: fit mDeg [0-8]!\n\n");
    return 1;
  }
  mDeg = -1;
  sscanf(argv[1], "%d", &mDeg);
  if ((mDeg < 0) || (mDeg > 8)) {
    printf("\n usage: fit mDeg [0-8]!\n\n");
    return 1;
  }
    
  /* initialize */
  stat = EG_open(&context);
  printf(" EG_open           = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  
  /* make a conical surface */
  for (i = 0; i < 12; i++) info[i] = 0.0;
  info[3]  = info[10] = 1.0;
  info[8]  = -1.0;
  info[12] = -0.463648;
  info[13] =  1.0;

  stat = EG_makeGeometry(context, SURFACE, CONICAL, NULL, NULL, info, &cone);
  printf(" EG_makeGeometry   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;

  stat = EG_getRange(cone, info, &i);
  printf(" EG_getRange       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;

  printf("                 U = %lf - %lf,  V = %le - %le,  per = %d\n",
         info[0], info[1], info[2], info[3], i);
  info[2]  = 0.0;
  info[3]  = 2.236068;
  sizes[0] = 32;
  sizes[1] = 16;
  stat     = EG_makeTessGeom(cone, info, sizes, &tess);
  printf(" EG_makeTessGeom   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  
  stat = EG_getTessGeom(tess, sizes, &xyzs);
  printf(" EG_getRange       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  
  /* make the apex degenerate at fp precision */
  for (i = sizes[0]*sizes[1]-32; i < sizes[0]*sizes[1]; i++) {
/*  printf("  %d:  %lf %lf %lf\n", i, xyzs[3*i  ], xyzs[3*i+1], xyzs[3*i+2]); */
    xyzs[3*i  ] = xyzs[3*i+2] = 0.0;
    xyzs[3*i+1] = 2.0;
  }
  
  /* fit the tessellated cone */
  stat = EG_approximate(context, mDeg, 1.e-7, sizes, xyzs, &bspline);
  printf(" EG_approximate    = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  EG_deleteObject(tess);
  EG_deleteObject(cone);
  
  stat = EG_getGeometry(bspline, &oclass, &mtype, &ref, &splInfo, &xyzs);
  printf(" EG_getGeometry    = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  stat = EG_getRange(bspline, info, &i);
  printf(" EG_getRange       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  printf("                 U = %lf - %lf,  V = %lf - %lf,  per = %d\n",
         info[0], info[1], info[2], info[3], i);
  printf("                     U Deg, nCp, nKnot = %d %d %d, V = %d %d %d\n",
         splInfo[1], splInfo[2], splInfo[3], splInfo[4], splInfo[5], splInfo[6]);
  
  /* make the facebody */
  stat = EG_makeFace(bspline, SFORWARD, info, &face);
  printf(" EG_makeFace       = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  
  stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &face, NULL,
                         &body);
  printf(" EG_makeTopology   = %d\n", stat);
  if (stat != EGADS_SUCCESS) return 1;
  EG_deleteObject(face);
  EG_deleteObject(bspline);

  /* make a model and write it out */
  if (body != NULL) {
    printf(" EG_makeTopology   = %d\n", EG_makeTopology(context, NULL, MODEL, 0,
                                                        NULL, 1, &body, NULL,
                                                        &model));
    printf(" EG_saveModel      = %d\n", EG_saveModel(model, "fit.egads"));
    printf("\n");
    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  }
                                                
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
