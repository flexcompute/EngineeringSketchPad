/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Cart3D Export Example
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

//#define BUILTIN


/*
 * 	calculates and returns a complete Body tessellation
 *
 * 	where:	body	- ego of a body tessellation
 *              nface   - number of faces in the body
 *              face    - the EGADS Face objects
 * 		nvert	- Number of vertices (returned)
 * 		verts	- coordinates (returned) 3*nverts in len -- freeable
 * 		ntriang	- number of triangles (returned)
 * 		triang	- triangle indices (returned) 3*ntriang in len
 *			  freeable
 *              comp    - Cart3D component ID per triangle -- freeable
 */

static int
bodyTessellation(ego tess, int nface, ego *faces, int *nvert, double **verts,
                 int *ntriang, int **triang, int **comp)
{
  int          status, i, j, k, base, npts, ntri, cID, *tri, *table, *compID;
  int          plen, tlen, atype, alen, nattr;
  const char   *string, *aname;
  const int    *tris, *tric, *ptype, *pindex, *ints;
  double       *xyzs;
  const double *points, *uv, *reals;

  *nvert  = *ntriang = 0;
  *verts  = NULL;
  *triang = NULL;
  *comp   = NULL;

  npts = ntri = 0;
  for (i = 1; i <= nface; i++) {
    status = EG_getTessFace(tess, i, &plen, &points, &uv, &ptype, &pindex,
                            &tlen, &tris, &tric);
    if (status != EGADS_SUCCESS) {
      printf(" Face %d: EG_getTessFace status = %d (bodyTessellation)!\n",
             i, status);
    } else {
      npts += plen;
      ntri += tlen;
    }
  }

  /* get the memory associated with the points */

  table = (int *) EG_alloc(2*npts*sizeof(int));
  if (table == NULL) {
    printf(" Error: Can not allocate node table (bodyTessellation)!\n");
    return EGADS_MALLOC;
  }
  xyzs = (double *) EG_alloc(3*npts*sizeof(double));
  if (xyzs == NULL) {
    printf(" Error: Can not allocate XYZs (bodyTessellation)!\n");
    EG_free(table);
    return EGADS_MALLOC;
  }

  /* zipper up the edges -- a Face at a time */

  npts = 0;
  for (j = 1; j <= nface; j++) {
    status = EG_getTessFace(tess, j, &plen, &points, &uv, &ptype, &pindex,
                            &tlen, &tris, &tric);
    if (status != EGADS_SUCCESS) continue;
    for (i = 0; i < plen; i++) {
      table[2*npts  ] = ptype[i];
      table[2*npts+1] = pindex[i];
      xyzs[3*npts  ]  = points[3*i  ];
      xyzs[3*npts+1]  = points[3*i+1];
      xyzs[3*npts+2]  = points[3*i+2];
      /* for non-interior pts -- try to match with others */
      if (ptype[i] != -1)
        for (k = 0; k < npts; k++)
          if ((table[2*k]==table[2*npts]) && (table[2*k+1]==table[2*npts+1])) {
            table[2*npts  ] = k;
            table[2*npts+1] = 0;
            break;
          }
      npts++;
    }
  }

  /* fill up the whole triangle list -- a Face at a time */

  tri = (int *) EG_alloc(3*ntri*sizeof(int));
  if (tri == NULL) {
    printf(" Error: Can not allocate triangles (bodyTessellation)!\n");
    EG_free(xyzs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  compID = (int *) EG_alloc(ntri*sizeof(int));
  if (compID == NULL) {
    printf(" Error: Can not allocate components (bodyTessellation)!\n");
    EG_free(tri);
    EG_free(xyzs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ntri; i++) compID[i] = 1;
  ntri = base = 0;
  for (j = 1; j <= nface; j++) {

    /* look for component ID */
    cID    = 1;
    status = EG_attributeNum(faces[j-1], &nattr);
    if (status == EGADS_SUCCESS) {
      for (k = 0; k < nattr; k++) {
        status = EG_attributeGet(faces[j-1], k+1, &aname, &atype, &alen, &ints,
                                 &reals, &string);
        if (status != EGADS_SUCCESS) continue;
        if (atype == ATTRINT) {
          if (strncmp(aname,"CartBC",6) == 0) {
            cID = ints[0];
            printf(" Face %d: Component ID = %d\n", j, cID);
          }
        } else if (atype == ATTRREAL) {
          if (strncmp(aname,"CartBC",6) == 0) {
            cID = reals[0] + 0.00001;
            printf(" Face %d: Component ID = %d\n", j, cID);
          }
        }
      }
    }

    /* get the face tessellation and store it away */
    status = EG_getTessFace(tess, j, &plen, &points, &uv, &ptype, &pindex,
                            &tlen, &tris, &tric);
    if (status != EGADS_SUCCESS) continue;
    for (i = 0; i < tlen; i++, ntri++) {
      k = tris[3*i  ] + base;
      if (table[2*k-1] == 0) {
        tri[3*ntri  ] = table[2*k-2] + 1;
      } else {
        tri[3*ntri  ] = k;
      }
      k = tris[3*i+1] + base;
      if (table[2*k-1] == 0) {
        tri[3*ntri+1] = table[2*k-2] + 1;
      } else {
        tri[3*ntri+1] = k;
      }
      k = tris[3*i+2] + base;
      if (table[2*k-1] == 0) {
        tri[3*ntri+2] = table[2*k-2] + 1;
      } else {
        tri[3*ntri+2] = k;
      }
      compID[ntri] = cID;
    }
    base += plen;
  }

  /* remove the unused points -- crunch the point list
   *    NOTE: the returned pointer verts has the full length (not realloc'ed)
   */

  for (i = 0; i <   npts; i++) table[i] = 0;
  for (i = 0; i < 3*ntri; i++) table[tri[i]-1]++;
  for (plen = i = 0; i < npts; i++) {
    if (table[i] == 0) continue;
    xyzs[3*plen  ] = xyzs[3*i  ];
    xyzs[3*plen+1] = xyzs[3*i+1];
    xyzs[3*plen+2] = xyzs[3*i+2];
    plen++;
    table[i] = plen;
  }
  /* reset the triangle indices */
  for (i = 0; i < 3*ntri; i++) {
    k      = tri[i]-1;
    tri[i] = table[k];
  }
  EG_free(table);

  *nvert   = plen;
  *verts   = xyzs;
  *ntriang = ntri;
  *triang  = tri;
  *comp    = compID;
  return EGADS_SUCCESS;
}


int main(int argc, char *argv[])
{
  int        i, j, status, oclass, mtype, nbody, nvert, ntriang, nface;
  int        *triang, *comp;
  char       filename[24];
  const char *OCCrev;
  float      arg;
  double     params[3], box[6], size, *verts;
  FILE       *fp;
  ego        context, model, geom, solid, *bodies, tess, *dum, *faces;

  if ((argc != 2) && (argc != 5)) {
    printf(" Usage: egads2cart Model [angle relSide relSag]\n\n");
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
                          &bodies, &triang);
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
    snprintf(filename, 24, "egads.%3.3d.a.tri", i+1);
    solid = bodies[i];

    mtype = 0;
    EG_getTopology(bodies[i], &geom, &oclass, &mtype, NULL, &j, &dum, &triang);
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

#ifdef BUILTIN
    status = EG_statusTessBody(tess, &geom, &j, &nvert);
    printf(" statusTessBody = %d %d  npts = %d\n", status, j, nvert);
#endif

    /* zip up the tessellation */

    status = bodyTessellation(tess, nface, faces, &nvert, &verts,
                              &ntriang, &triang, &comp);
#ifdef BUILTIN
    for (j = 0; j < nvert; j++) {
      int    stat, ptype, pindex;
      double coord[3];

      stat = EG_getGlobal(tess, j+1, &ptype, &pindex, coord);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getGlobal %d/%d = %d\n", j+1, nvert, stat);
        continue;
      }
      if ((coord[0] == verts[3*j  ]) && (coord[1] == verts[3*j+1]) &&
          (coord[2] == verts[3*j+2])) continue;
      printf("  Global %6d: %lf %lf %lf  %lf %lf %lf\n", j+1, coord[0], coord[1],
             coord[2], verts[3*j  ], verts[3*j+1], verts[3*j+2]);
    }
#endif
    EG_deleteObject(tess);
    EG_free(faces);
    if (status != EGADS_SUCCESS) continue;

    /* write it out */

    fp = fopen(filename, "w");
    if (fp == NULL) {
      printf(" Can not Open file %s! NO FILE WRITTEN\n", filename);
      EG_free(verts);
      EG_free(triang);
      if (solid != bodies[i]) EG_deleteObject(solid);
      continue;
    }
    printf("\nWriting Cart3D component file %s\n", filename);
    /* header */
    fprintf(fp, "%d  %d\n", nvert, ntriang);
    /* ...vertList     */
    for (j = 0; j < nvert; j++)
      fprintf(fp, " %20.13le %20.13le %20.13le\n",
              verts[3*j  ], verts[3*j+1], verts[3*j+2]);
    /* ...Connectivity */
    for (j = 0; j < ntriang; j++)
      fprintf(fp, "%6d %6d %6d\n",triang[3*j], triang[3*j+1], triang[3*j+2]);
    /* ...Component list*/
    for (j = 0; j < ntriang; j++)
      fprintf(fp, "%6d\n", comp[j]);
    fclose(fp);

    printf("      # verts = %d,  # tris = %d\n\n", nvert, ntriang);
    EG_free(verts);
    EG_free(triang);
    EG_free(comp);
    if (solid != bodies[i]) EG_deleteObject(solid);
  }

  status = EG_deleteObject(model);
  if (status != EGADS_SUCCESS) printf(" EG_deleteObject = %d\n", status);
  EG_close(context);

  return 0;
}
