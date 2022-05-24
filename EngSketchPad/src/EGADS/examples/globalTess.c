/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Global Tessellation Tester
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"
#include <math.h>

//#define DISPLACE
#define GLOBALVERTS


int main(int argc, char *argv[])
{
  int          i, j, k, n, status, oclass, mtype, nbody, nvert, nface, nedge;
  int          pt, pi, plen, tlen, global, *sens;
  float        arg;
  double       params[3], box[6], coord[3], size;
#ifdef DISPLACE
  double       *displace = NULL;
#endif
  ego          context, model, geom, *bodies, tess, copy, *dum, *faces, *edges;
  const char   *OCCrev;
  const double *points, *uv;
  const int    *tris, *tric, *ptype, *pindex;

  if ((argc != 2) && (argc != 5)) {
    printf(" Usage: globalTess Model [angle relSide relSag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d %s\n\n", i, j, OCCrev);

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
                          &bodies, &sens);
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

  /* test each body */

  for (i = 0; i < nbody; i++) {

    status = EG_getTopology(bodies[i], &geom, &oclass, &mtype, NULL, &j, &dum,
                            &sens);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getTopology Body %d = %d\n", i+1, status);
      continue;
    }
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", i+1);
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", i+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", i+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", i+1);
    }

    status = EG_makeTessBody(bodies[i], params, &tess);
    if (status != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", i, status);
      continue;
    }
    status = EG_getBodyTopos(bodies[i], NULL, FACE, &nface, &faces);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos Face %d = %d\n", i+1, status);
      EG_deleteObject(tess);
      continue;
    }
    status = EG_getBodyTopos(bodies[i], NULL, EDGE, &nedge, &edges);
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos Edge %d = %d\n", i+1, status);
      EG_free(faces);
      EG_deleteObject(tess);
      continue;
    }

    status = EG_statusTessBody(tess, &geom, &j, &nvert);
    printf(" statusTessBody      = %d %d  npts = %d\n", status, j, nvert);
#ifdef DISPLACE
    displace = (double *) malloc(3*nvert*sizeof(double));
    if (displace == NULL) {
      printf(" Malloc of %d displacements\n", nvert);
      EG_free(edges);
      EG_free(faces);
      EG_deleteObject(tess);
      continue;
    }
    for (j = 0; j < nvert; j++) {
      displace[3*j  ] = 1.0;
      displace[3*j+1] = 0.0;
      displace[3*j+2] = 0.0;
    }
    status = EG_copyObject(tess, displace, &copy);
    free(displace);
#else
    status = EG_copyObject(tess, NULL, &copy);
#endif
    if (status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos Edge %d = %d\n", i+1, status);
      EG_free(edges);
      EG_free(faces);
      EG_deleteObject(tess);
      continue;
    }
    status = EG_statusTessBody(copy, &geom, &j, &nvert);
    printf(" statusTessBody copy = %d %d  npts = %d\n", status, j, nvert);

#ifdef GLOBALVERTS
    for (global = 1; global <= nvert; global++) {
        status = EG_getGlobal(copy, global, &pt, &pi, coord);
        if (status != EGADS_SUCCESS) {
          printf(" Body %d/Vert %d: EG_getGlobal = %d\n",
                 i+1, global, status);
        }
        printf("  GlobalF %6d: %lf %lf %lf\n", global, coord[0],
               coord[1], coord[2]);
    }
#endif

    for (j = 1; j <= nface; j++) {
      status = EG_getTessFace(tess, j, &plen, &points, &uv, &ptype, &pindex,
                              &tlen, &tris, &tric);
      if (status != EGADS_SUCCESS) {
        printf(" %d EG_getTessFace %d = %d\n", i+1, j, status);
        continue;
      }

      for (k = 0; k < plen; k++) {
        status = EG_localToGlobal(tess, j, k+1, &global);
        if (status != EGADS_SUCCESS) {
          printf(" Body %d/Face %d/Vert %d: EG_localToGlobal = %d\n",
                 i+1, j, k+1, status);
          continue;
        }
        status = EG_getGlobal(copy, global, &pt, &pi, coord);
        if (status != EGADS_SUCCESS) {
          printf(" Body %d/Face %d/Vert %d: EG_getGlobal = %d\n",
                 i+1, j, k+1, status);
          continue;
        }
        if (ptype[k] == -1) {
          if ((pt != -k-1) || (pi != j))
            printf("  GlobalF %6d: %d %d  %d %d\n", global, pt, pi, -k-1, j);
        } else {
          if ((pt != ptype[k]) || (pi != pindex[k]))
            printf("  GlobalF %6d: %d %d  %d %d\n", global, pt, pi,
                   ptype[k], pindex[k]);
        }
#ifdef DISPLACE
        coord[0] -= 1.0;
        if ((fabs(coord[0]-points[3*k  ]) < 1.e-8) &&
            (fabs(coord[1]-points[3*k+1]) < 1.e-8) &&
            (fabs(coord[2]-points[3*k+2]) < 1.e-8)) continue;
#else
        if ((coord[0] == points[3*k  ]) && (coord[1] == points[3*k+1]) &&
            (coord[2] == points[3*k+2])) continue;
#endif
        printf("  GlobalF %6d: %lf %lf %lf  %lf %lf %lf\n", global, coord[0],
               coord[1], coord[2], points[3*k  ], points[3*k+1], points[3*k+2]);
      }
    }

    for (j = 1; j <= nedge; j++) {
      status = EG_getTessEdge(tess, j, &plen, &points, &uv);
      if (status == EGADS_DEGEN) continue;
      if (status != EGADS_SUCCESS) {
        printf(" %d EG_getTessEdge %d = %d\n", i+1, j, status);
        continue;
      }

      for (k = 0; k < plen; k++) {
        status = EG_localToGlobal(tess, -j, k+1, &global);
        if (status != EGADS_SUCCESS) {
          printf(" Body %d/Edge %d/Vert %d: EG_localToGlobal = %d\n",
                 i+1, j, k+1, status);
          continue;
        }
        status = EG_getGlobal(copy, global, &pt, &pi, coord);
        if (status != EGADS_SUCCESS) {
          printf(" Body %d/Edge %d/Vert %d: EG_getGlobal = %d\n",
                 i+1, j, k+1, status);
          continue;
        }
        if ((pt == 0) && ((k != 0) && (k != plen-1)))
          printf("  GlobalE %6d:  Node @ %d/%d!\n", global, k, plen-1);
        if (pt != 0)
          if ((pt != k+1) || (pi != j))
            printf("  GlobalE %6d: %d %d  %d %d\n", global, pt, pi,
                   k+1, j);
#ifdef DISPLACE
        coord[0] -= 1.0;
        if ((fabs(coord[0]-points[3*k  ]) < 1.e-8) &&
            (fabs(coord[1]-points[3*k+1]) < 1.e-8) &&
            (fabs(coord[2]-points[3*k+2]) < 1.e-8)) continue;
#else
        if ((coord[0] == points[3*k  ]) && (coord[1] == points[3*k+1]) &&
            (coord[2] == points[3*k+2])) continue;
#endif
        printf("  GlobalE %6d: %lf %lf %lf  %lf %lf %lf\n", global, coord[0],
               coord[1], coord[2], points[3*k  ], points[3*k+1], points[3*k+2]);
      }
    }
    
    status = EG_getBodyTopos(bodies[i], NULL, NODE, &n, NULL);
    if (status == EGADS_SUCCESS)
      for (j = 1; j <= n; j++) {
        status = EG_localToGlobal(tess, 0, j, &global);
        if (status != EGADS_SUCCESS) {
          printf(" Body %d/Node %d: EG_localToGlobal = %d\n", i+1, j, status);
        } else {
          printf(" Body %d/Node %d: %d\n", i+1, j, global);
        }
      }

    EG_deleteObject(copy);
    EG_deleteObject(tess);
    EG_free(edges);
    EG_free(faces);
    printf("\n");
  }

  status = EG_deleteObject(model);
  if (status != EGADS_SUCCESS) printf(" EG_deleteObject = %d\n", status);
  EG_close(context);

  return 0;
}
