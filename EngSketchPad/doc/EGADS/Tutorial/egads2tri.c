/*
 *      EGADS: Engineering Geometry Aircraft Design System
 *
 *             Cart3D Tri File Export Example
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"
#include <math.h>
#include <string.h>

#ifdef WIN32
#define snprintf _snprintf
#endif

int main(int argc, char *argv[])
{
  int          i, j, k, status, oclass, mtype, nbody, nvert, ntriang, nface;
  int          plen, tlen, pty, pin, tin[3], *senses;
  const int    *ptype, *pindex, *tris, *tric;
  char         filename[20];
  const char   *OCCrev;
  float        arg;
  double       params[3], box[6], size, verts[3];
  const double *points, *uv;
  FILE         *fp;
  ego          context, model, geom, solid, *bodies, tess, *dum, *faces;

  if ((argc != 2) && (argc != 5)) {
    printf(" Usage: egads2tri Model [angle relSide relSag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  status = EG_open(&context);
  if (status != EGADS_SUCCESS) {
    printf(" EG_open = %d!\n\n", status);
    return 1;
  }
  status = EG_loadModel(context, 0, argv[1], &model);
  if (status != EGADS_SUCCESS) {
    printf(" EG_loadModel = %d\n\n", status);
    return 1;
  }
  status = EG_getBoundingBox(model, box);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getBoundingBox = %d\n\n", status);
    return 1;
  }
  size = sqrt((box[0]-box[3])*(box[0]-box[3]) + (box[1]-box[4])*(box[1]-box[4]) +
              (box[2]-box[5])*(box[2]-box[5]));

  /* get all bodies */
  status = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                          &bodies, &senses);
  if (status != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n\n", status);
    return 1;
  }

  params[0] =  0.025*size;
  params[1] =  0.001*size;
  params[2] = 15.0;
  if (argc == 5) {
    sscanf(argv[2], "%f", &arg);
    params[2] = arg;
    sscanf(argv[3], "%f", &arg);
    params[0] = arg;
    sscanf(argv[4], "%f", &arg);
    params[1] = arg;
    printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
           params[2], params[0], params[1]);
    params[0] *= size;
    params[1] *= size;
  }

  printf(" Number of Bodies = %d\n\n", nbody);

  /* write out each body as a different Cart3D ASCII tri file */

  for (i = 0; i < nbody; i++) {
    snprintf(filename, 20, "egads.%3.3d.a.tri", i+1);
    solid = bodies[i];
    
    mtype = 0;
    EG_getTopology(bodies[i], &geom, &oclass, &mtype, NULL, &j, &dum, &senses);
    if (mtype == SHEETBODY) {
      status = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, j, dum,
                               NULL, &solid);
      if (status == EGADS_SUCCESS) {
        printf(" SheetBody %d promoted to SolidBody\n", i);
        mtype = SOLIDBODY;
      } else {
        printf(" SheetBody %d cannot be promoted to SolidBody\n", i);
      }
    }
    if (mtype != SOLIDBODY) continue;   /* only Solid Bodies! */

    status = EG_makeTessBody(solid, params, &tess);
    if (status != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", i, status);
      if (solid != bodies[i]) EG_deleteObject(solid);
      continue;
    }
    status = EG_getBodyTopos(solid, NULL, FACE, &nface, &faces);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos %d = %d\n", i, status);
      if (solid != bodies[i]) EG_deleteObject(solid);
      EG_deleteObject(tess);
      continue;
    }
    EG_free(faces);

    /* get counts */
    status = EG_statusTessBody(tess, &geom, &j, &nvert);
    printf(" statusTessBody = %d %d  npts = %d\n", status, j, nvert);
    if (status != EGADS_SUCCESS) continue;
    ntriang = 0;
    for (j = 0; j < nface; j++) {
      status = EG_getTessFace(tess, j+1, &plen, &points, &uv, &ptype, &pindex,
                              &tlen, &tris, &tric);
      if (status != EGADS_SUCCESS) {
        printf(" Error: EG_getTessFace %d/%d = %d\n", j+1, nvert, status);
        continue;
      }
      ntriang += tlen;
    }

    /* write it out */

    fp = fopen(filename, "w");
    if (fp == NULL) {
      printf(" Can not Open file %s! NO FILE WRITTEN\n", filename);
      if (solid != bodies[i]) EG_deleteObject(solid);
      continue;
    }
    printf("\nWriting Cart3D component tri file %s\n", filename);
    /* header */
    fprintf(fp, "%d  %d\n", nvert, ntriang);
    /* ...vertList     */
    for (j = 0; j < nvert; j++) {
      status = EG_getGlobal(tess, j+1, &pty, &pin, verts);
      if (status != EGADS_SUCCESS)
        printf(" Error: EG_getGlobal %d/%d = %d\n", j+1, nvert, status);
      fprintf(fp, " %20.13le %20.13le %20.13le\n", verts[0],verts[1],verts[2]);
    }
    /* ...Connectivity */
    for (j = 0; j < nface; j++) {
      status = EG_getTessFace(tess, j+1, &plen, &points, &uv, &ptype, &pindex,
                              &tlen, &tris, &tric);
      if (status != EGADS_SUCCESS) continue;
      for (k = 0; k < tlen; k++) {
        status = EG_localToGlobal(tess, j+1, tris[3*k  ], &tin[0]);
        if (status != EGADS_SUCCESS)
          printf(" Error: EG_localToGlobal %d/%d = %d\n", j+1, tris[3*k  ],
                 status);
        status = EG_localToGlobal(tess, j+1, tris[3*k+1], &tin[1]);
        if (status != EGADS_SUCCESS)
          printf(" Error: EG_localToGlobal %d/%d = %d\n", j+1, tris[3*k+1],
                 status);
        status = EG_localToGlobal(tess, j+1, tris[3*k+2], &tin[2]);
        if (status != EGADS_SUCCESS)
          printf(" Error: EG_localToGlobal %d/%d = %d\n", j+1, tris[3*k+2],
                 status);
        fprintf(fp, "%6d %6d %6d\n",tin[0], tin[1], tin[2]);
      }
    }
    /* ...Component list*/
    for (j = 0; j < ntriang; j++) fprintf(fp, "%6d\n", 1);
    fclose(fp);
    
    if (solid != bodies[i]) EG_deleteObject(solid);
  }

  status = EG_deleteObject(tess);
  if (status != EGADS_SUCCESS) printf(" EG_deleteObject tess  = %d\n", status);
  status = EG_deleteObject(model);
  if (status != EGADS_SUCCESS) printf(" EG_deleteObject model = %d\n", status);
  EG_close(context);

  return 0;
}
