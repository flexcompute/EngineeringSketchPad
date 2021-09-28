/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Display the Curvature using wv (the WebViewer)
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

#include "egads.h"
#include "wsserver.h"

#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif
#include <winsock2.h>
#endif

#define CROSS(a,b,c)  a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                      a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                      a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)     (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


/* structure to hold on to the EGADS triangulation per Body */
typedef struct {
  ego *faces;
  ego *edges;
  ego body;
  ego tess;
  int mtype;
  int nfaces;
  int nedges;
} bodyData;


/* globals used in these functions */
  static int        nbody, key = -1;
  static double     params[3];
  static float      focus[4], lims[2] = {-1.0, +1.0};
  static wvContext *cntxt;
  static bodyData  *bodydata;
  static int        sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static char       *keys[3] = {"GC", "RC1", "RC2"};


/* blue-white-red spectrum */
static float color_map[256*3] =
             { 0.0000, 0.0000, 1.0000,    0.0078, 0.0078, 1.0000,
               0.0156, 0.0156, 1.0000,    0.0234, 0.0234, 1.0000,
               0.0312, 0.0312, 1.0000,    0.0391, 0.0391, 1.0000,
               0.0469, 0.0469, 1.0000,    0.0547, 0.0547, 1.0000,
               0.0625, 0.0625, 1.0000,    0.0703, 0.0703, 1.0000,
               0.0781, 0.0781, 1.0000,    0.0859, 0.0859, 1.0000,
               0.0938, 0.0938, 1.0000,    0.1016, 0.1016, 1.0000,
               0.1094, 0.1094, 1.0000,    0.1172, 0.1172, 1.0000,
               0.1250, 0.1250, 1.0000,    0.1328, 0.1328, 1.0000,
               0.1406, 0.1406, 1.0000,    0.1484, 0.1484, 1.0000,
               0.1562, 0.1562, 1.0000,    0.1641, 0.1641, 1.0000,
               0.1719, 0.1719, 1.0000,    0.1797, 0.1797, 1.0000,
               0.1875, 0.1875, 1.0000,    0.1953, 0.1953, 1.0000,
               0.2031, 0.2031, 1.0000,    0.2109, 0.2109, 1.0000,
               0.2188, 0.2188, 1.0000,    0.2266, 0.2266, 1.0000,
               0.2344, 0.2344, 1.0000,    0.2422, 0.2422, 1.0000,
               0.2500, 0.2500, 1.0000,    0.2578, 0.2578, 1.0000,
               0.2656, 0.2656, 1.0000,    0.2734, 0.2734, 1.0000,
               0.2812, 0.2812, 1.0000,    0.2891, 0.2891, 1.0000,
               0.2969, 0.2969, 1.0000,    0.3047, 0.3047, 1.0000,
               0.3125, 0.3125, 1.0000,    0.3203, 0.3203, 1.0000,
               0.3281, 0.3281, 1.0000,    0.3359, 0.3359, 1.0000,
               0.3438, 0.3438, 1.0000,    0.3516, 0.3516, 1.0000,
               0.3594, 0.3594, 1.0000,    0.3672, 0.3672, 1.0000,
               0.3750, 0.3750, 1.0000,    0.3828, 0.3828, 1.0000,
               0.3906, 0.3906, 1.0000,    0.3984, 0.3984, 1.0000,
               0.4062, 0.4062, 1.0000,    0.4141, 0.4141, 1.0000,
               0.4219, 0.4219, 1.0000,    0.4297, 0.4297, 1.0000,
               0.4375, 0.4375, 1.0000,    0.4453, 0.4453, 1.0000,
               0.4531, 0.4531, 1.0000,    0.4609, 0.4609, 1.0000,
               0.4688, 0.4688, 1.0000,    0.4766, 0.4766, 1.0000,
               0.4844, 0.4844, 1.0000,    0.4922, 0.4922, 1.0000,
               0.5000, 0.5000, 1.0000,    0.5078, 0.5078, 1.0000,
               0.5156, 0.5156, 1.0000,    0.5234, 0.5234, 1.0000,
               0.5312, 0.5312, 1.0000,    0.5391, 0.5391, 1.0000,
               0.5469, 0.5469, 1.0000,    0.5547, 0.5547, 1.0000,
               0.5625, 0.5625, 1.0000,    0.5703, 0.5703, 1.0000,
               0.5781, 0.5781, 1.0000,    0.5859, 0.5859, 1.0000,
               0.5938, 0.5938, 1.0000,    0.6016, 0.6016, 1.0000,
               0.6094, 0.6094, 1.0000,    0.6172, 0.6172, 1.0000,
               0.6250, 0.6250, 1.0000,    0.6328, 0.6328, 1.0000,
               0.6406, 0.6406, 1.0000,    0.6484, 0.6484, 1.0000,
               0.6562, 0.6562, 1.0000,    0.6641, 0.6641, 1.0000,
               0.6719, 0.6719, 1.0000,    0.6797, 0.6797, 1.0000,
               0.6875, 0.6875, 1.0000,    0.6953, 0.6953, 1.0000,
               0.7031, 0.7031, 1.0000,    0.7109, 0.7109, 1.0000,
               0.7188, 0.7188, 1.0000,    0.7266, 0.7266, 1.0000,
               0.7344, 0.7344, 1.0000,    0.7422, 0.7422, 1.0000,
               0.7500, 0.7500, 1.0000,    0.7578, 0.7578, 1.0000,
               0.7656, 0.7656, 1.0000,    0.7734, 0.7734, 1.0000,
               0.7812, 0.7812, 1.0000,    0.7891, 0.7891, 1.0000,
               0.7969, 0.7969, 1.0000,    0.8047, 0.8047, 1.0000,
               0.8125, 0.8125, 1.0000,    0.8203, 0.8203, 1.0000,
               0.8281, 0.8281, 1.0000,    0.8359, 0.8359, 1.0000,
               0.8438, 0.8438, 1.0000,    0.8516, 0.8516, 1.0000,
               0.8594, 0.8594, 1.0000,    0.8672, 0.8672, 1.0000,
               0.8750, 0.8750, 1.0000,    0.8828, 0.8828, 1.0000,
               0.8906, 0.8906, 1.0000,    0.8984, 0.8984, 1.0000,
               0.9062, 0.9062, 1.0000,    0.9141, 0.9141, 1.0000,
               0.9219, 0.9219, 1.0000,    0.9297, 0.9297, 1.0000,
               0.9375, 0.9375, 1.0000,    0.9453, 0.9453, 1.0000,
               0.9531, 0.9531, 1.0000,    0.9609, 0.9609, 1.0000,
               0.9688, 0.9688, 1.0000,    0.9766, 0.9766, 1.0000,
               0.9844, 0.9844, 1.0000,    0.9922, 0.9922, 1.0000,
               1.0000, 1.0000, 1.0000,    1.0000, 0.9922, 0.9922,
               1.0000, 0.9844, 0.9844,    1.0000, 0.9766, 0.9766,
               1.0000, 0.9688, 0.9688,    1.0000, 0.9609, 0.9609,
               1.0000, 0.9531, 0.9531,    1.0000, 0.9453, 0.9453,
               1.0000, 0.9375, 0.9375,    1.0000, 0.9297, 0.9297,
               1.0000, 0.9219, 0.9219,    1.0000, 0.9141, 0.9141,
               1.0000, 0.9062, 0.9062,    1.0000, 0.8984, 0.8984,
               1.0000, 0.8906, 0.8906,    1.0000, 0.8828, 0.8828,
               1.0000, 0.8750, 0.8750,    1.0000, 0.8672, 0.8672,
               1.0000, 0.8594, 0.8594,    1.0000, 0.8516, 0.8516,
               1.0000, 0.8438, 0.8438,    1.0000, 0.8359, 0.8359,
               1.0000, 0.8281, 0.8281,    1.0000, 0.8203, 0.8203,
               1.0000, 0.8125, 0.8125,    1.0000, 0.8047, 0.8047,
               1.0000, 0.7969, 0.7969,    1.0000, 0.7891, 0.7891,
               1.0000, 0.7812, 0.7812,    1.0000, 0.7734, 0.7734,
               1.0000, 0.7656, 0.7656,    1.0000, 0.7578, 0.7578,
               1.0000, 0.7500, 0.7500,    1.0000, 0.7422, 0.7422,
               1.0000, 0.7344, 0.7344,    1.0000, 0.7266, 0.7266,
               1.0000, 0.7188, 0.7188,    1.0000, 0.7109, 0.7109,
               1.0000, 0.7031, 0.7031,    1.0000, 0.6953, 0.6953,
               1.0000, 0.6875, 0.6875,    1.0000, 0.6797, 0.6797,
               1.0000, 0.6719, 0.6719,    1.0000, 0.6641, 0.6641,
               1.0000, 0.6562, 0.6562,    1.0000, 0.6484, 0.6484,
               1.0000, 0.6406, 0.6406,    1.0000, 0.6328, 0.6328,
               1.0000, 0.6250, 0.6250,    1.0000, 0.6172, 0.6172,
               1.0000, 0.6094, 0.6094,    1.0000, 0.6016, 0.6016,
               1.0000, 0.5938, 0.5938,    1.0000, 0.5859, 0.5859,
               1.0000, 0.5781, 0.5781,    1.0000, 0.5703, 0.5703,
               1.0000, 0.5625, 0.5625,    1.0000, 0.5547, 0.5547,
               1.0000, 0.5469, 0.5469,    1.0000, 0.5391, 0.5391,
               1.0000, 0.5312, 0.5312,    1.0000, 0.5234, 0.5234,
               1.0000, 0.5156, 0.5156,    1.0000, 0.5078, 0.5078,
               1.0000, 0.5000, 0.5000,    1.0000, 0.4922, 0.4922,
               1.0000, 0.4844, 0.4844,    1.0000, 0.4766, 0.4766,
               1.0000, 0.4688, 0.4688,    1.0000, 0.4609, 0.4609,
               1.0000, 0.4531, 0.4531,    1.0000, 0.4453, 0.4453,
               1.0000, 0.4375, 0.4375,    1.0000, 0.4297, 0.4297,
               1.0000, 0.4219, 0.4219,    1.0000, 0.4141, 0.4141,
               1.0000, 0.4062, 0.4062,    1.0000, 0.3984, 0.3984,
               1.0000, 0.3906, 0.3906,    1.0000, 0.3828, 0.3828,
               1.0000, 0.3750, 0.3750,    1.0000, 0.3672, 0.3672,
               1.0000, 0.3594, 0.3594,    1.0000, 0.3516, 0.3516,
               1.0000, 0.3438, 0.3438,    1.0000, 0.3359, 0.3359,
               1.0000, 0.3281, 0.3281,    1.0000, 0.3203, 0.3203,
               1.0000, 0.3125, 0.3125,    1.0000, 0.3047, 0.3047,
               1.0000, 0.2969, 0.2969,    1.0000, 0.2891, 0.2891,
               1.0000, 0.2812, 0.2812,    1.0000, 0.2734, 0.2734,
               1.0000, 0.2656, 0.2656,    1.0000, 0.2578, 0.2578,
               1.0000, 0.2500, 0.2500,    1.0000, 0.2422, 0.2422,
               1.0000, 0.2344, 0.2344,    1.0000, 0.2266, 0.2266,
               1.0000, 0.2188, 0.2188,    1.0000, 0.2109, 0.2109,
               1.0000, 0.2031, 0.2031,    1.0000, 0.1953, 0.1953,
               1.0000, 0.1875, 0.1875,    1.0000, 0.1797, 0.1797,
               1.0000, 0.1719, 0.1719,    1.0000, 0.1641, 0.1641,
               1.0000, 0.1562, 0.1562,    1.0000, 0.1484, 0.1484,
               1.0000, 0.1406, 0.1406,    1.0000, 0.1328, 0.1328,
               1.0000, 0.1250, 0.1250,    1.0000, 0.1172, 0.1172,
               1.0000, 0.1094, 0.1094,    1.0000, 0.1016, 0.1016,
               1.0000, 0.0938, 0.0938,    1.0000, 0.0859, 0.0859,
               1.0000, 0.0781, 0.0781,    1.0000, 0.0703, 0.0703,
               1.0000, 0.0625, 0.0625,    1.0000, 0.0547, 0.0547,
               1.0000, 0.0469, 0.0469,    1.0000, 0.0391, 0.0391,
               1.0000, 0.0312, 0.0312,    1.0000, 0.0234, 0.0234,
               1.0000, 0.0156, 0.0156,    1.0000, 0.0078, 0.0078 };


static void spec_col(float scalar, float color[])
{
  int   indx;
  float frac;

  if (lims[0] == lims[1]) {
    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 0.0;
  } else if (scalar <= lims[0]) {
    color[0] = color_map[0];
    color[1] = color_map[1];
    color[2] = color_map[2];
  } else if (scalar >= lims[1]) {
    color[0] = color_map[3*255  ];
    color[1] = color_map[3*255+1];
    color[2] = color_map[3*255+2];
  } else {
    frac  = 255.0*(scalar - lims[0])/(lims[1] - lims[0]);
    if (frac < 0  ) frac = 0;
    if (frac > 255) frac = 255;
    indx  = frac;
    frac -= indx;
    if (indx == 255) {
      indx--;
      frac += 1.0;
    }
    color[0] = frac*color_map[3*(indx+1)  ] + (1.0-frac)*color_map[3*indx  ];
    color[1] = frac*color_map[3*(indx+1)+1] + (1.0-frac)*color_map[3*indx+1];
    color[2] = frac*color_map[3*(indx+1)+2] + (1.0-frac)*color_map[3*indx+2];
  }
}


int main(int argc, char *argv[])
{
  int          i, j, k, m, n, ibody, stat, oclass, mtype, len, ntri, sum;
  int          nloops, nseg, nledges, *segs, *lsenses, *senses, *esenses;
  int          nh, *heads;
  const int    *tris, *tric, *ptype, *pindex;
  float        arg, color[3], *colrs, *lsegs;
  double       box[6], result[18], size, dot, norm[3], *norms, *pu, *pv;
  const double *xyzs, *uvs, *ts;
  char         gpname[46], *startapp;
  const char   *OCCrev;
  ego          context, model, geom, *bodies, *dum, *ledges, *loops;
  wvData       items[5];
  float        eye[3]    = {0.0, 0.0, 7.0};
  float        center[3] = {0.0, 0.0, 0.0};
  float        up[3]     = {0.0, 1.0, 0.0};

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vTess filename [angle maxlen sag]\n\n");
    return 1;
  }

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d %s\n\n", i, j, OCCrev);

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

  focus[0] = 0.5*(box[0]+box[3]);
  focus[1] = 0.5*(box[1]+box[4]);
  focus[2] = 0.5*(box[2]+box[5]);
  focus[3] = size;

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyData *) malloc(nbody*sizeof(bodyData));
  if (bodydata == NULL) {
    printf(" MALLOC Error on Body storage!\n");
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
  printf(" NOTE: tParams = %lf %lf %lf\n\n", params[0], params[1], params[2]);

  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
    mtype = 0;
    EG_getTopology(bodies[ibody], &geom, &oclass,
                   &mtype, NULL, &j, &dum, &senses);
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", ibody+1);
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", ibody+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", ibody+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", ibody+1);
    }

    stat = EG_getBodyTopos(bodies[ibody], NULL, FACE,
                           &bodydata[ibody].nfaces, &bodydata[ibody].faces);
    i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE,
                           &bodydata[ibody].nedges, &bodydata[ibody].edges);
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      return 1;
    }

    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
  }
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) printf(" failed to create wvContext!\n");

  /* make the scene */
  for (sum = stat = ibody = 0; ibody < nbody; ibody++) {

    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      snprintf(gpname, 34, "Body %d Face %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      color[0]  = 1.0;
      color[1]  = ibody;
      color[1] /= nbody;
      color[2]  = 0.0;
      colrs = (float *) malloc(3*len*sizeof(float));
      if (colrs == NULL) continue;
      for (j = 0; j < len; j++) {
        colrs[3*j  ]  = 1.0;
        colrs[3*j+1]  = ibody;
        colrs[3*j+1] /= nbody;
        colrs[3*j+2]  = 0.0;
      }
      stat = wv_setData(WV_REAL32, len, (void *) colrs,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(colrs);
      for (nseg = j = 0; j < ntri; j++)
        for (k = 0; k < 3; k++)
          if (tric[3*j+k] < j+1) nseg++;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      for (nseg = j = 0; j < ntri; j++)
        for (k = 0; k < 3; k++)
          if (tric[3*j+k] < j+1) {
            segs[2*nseg  ] = tris[3*j+sides[k][0]];
            segs[2*nseg+1] = tris[3*j+sides[k][1]];
            nseg++;
          }
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", i, gpname);
      free(segs);
/*    color[0] = color[1] = color[2] = 0.8;  */
      color[0] = color[1] = color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", i, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
                         WV_ON|WV_ORIENTATION|WV_SHADING, 5, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      sum += ntri;
    }

    /* put normals of faces in "edge" slot */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
      norms = (double *) malloc(6*len*sizeof(double));
      if (norms == NULL) continue;
      snprintf(gpname, 34, "Body %d Edge %d", ibody+1, i+1);

      pu = &result[3];
      pv = &result[6];
      for (j = 0; j < len; j++) {
        norms[6*j  ] = norms[6*j+3] = xyzs[3*j  ];
        norms[6*j+1] = norms[6*j+4] = xyzs[3*j+1];
        norms[6*j+2] = norms[6*j+5] = xyzs[3*j+2];
        stat = EG_evaluate(bodydata[ibody].faces[i], &uvs[2*j], result);
        if (stat != EGADS_SUCCESS) continue;
        CROSS(norm, pu, pv);
        dot = sqrt(DOT(norm, norm))*bodydata[ibody].faces[i]->mtype;
        if (dot == 0.0) continue;
        norm[0]      *= 0.025*focus[3]/dot;
        norm[1]      *= 0.025*focus[3]/dot;
        norm[2]      *= 0.025*focus[3]/dot;
        norms[6*j+3] += norm[0];
        norms[6*j+4] += norm[1];
        norms[6*j+5] += norm[2];
      }
      stat = wv_setData(WV_REAL64, 2*len, (void *) norms,  WV_VERTICES,
                        &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      free(norms);
      color[0] = 0.0;
      color[1] = 0.0;
      color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, 0, 2, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
    }

    /* get loops */
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTopology(bodydata[ibody].faces[i], &geom, &oclass,
                            &mtype, NULL, &nloops, &loops, &lsenses);
      if (stat != EGADS_SUCCESS) continue;
      for (nh = j = 0; j < nloops; j++) {
        stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, NULL,
                              &nledges, &ledges, &esenses);
        if (stat != EGADS_SUCCESS) continue;

        /* count */
        for (ntri = nseg = k = 0; k < nledges; k++) {
          m = 0;
          while (ledges[k] != bodydata[ibody].edges[m]) {
            m++;
            if (m == bodydata[ibody].nedges) break;
          }
          /* assume that the edge is degenerate and removed */
          if (m == bodydata[ibody].nedges) continue;
          stat = EG_getTessEdge(bodydata[ibody].tess, m+1, &len,
                                &xyzs, &ts);
          if (stat != EGADS_SUCCESS) {
            printf(" EG_getTessEdge %d = %d!\n", m+1, stat);
            nseg = 0;
            break;
          }
          if (len == 2)
            if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
                (xyzs[2] == xyzs[5])) continue;
          nh++;
          nseg += len;
          ntri += 2*(len-1);
        }
        if (nseg == 0) continue;
        lsegs = (float *) malloc(3*nseg*sizeof(float));
        if (lsegs == NULL) {
          printf(" Can not allocate %d Segments!\n", nseg);
          continue;
        }
        segs = (int *) malloc(ntri*sizeof(int));
        if (segs == NULL) {
          printf(" Can not allocate %d Line Segments!\n", ntri);
          free(lsegs);
          continue;
        }
        heads = (int *) malloc(nh*sizeof(int));
        if (heads == NULL) {
          printf(" Can not allocate %d Heads!\n", nh);
          free(segs);
          free(lsegs);
          continue;
        }

        /* fill */
        for (nh = ntri = nseg = k = 0; k < nledges; k++) {
          m = 0;
          while (ledges[k] != bodydata[ibody].edges[m]) {
            m++;
            if (m == bodydata[ibody].nedges) break;
          }
          /* assume that the edge is degenerate and removed */
          if (m == bodydata[ibody].nedges) continue;
          EG_getTessEdge(bodydata[ibody].tess, m+1, &len, &xyzs, &ts);
          if (len == 2)
            if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
                (xyzs[2] == xyzs[5])) continue;
          if (esenses[k] == -1) heads[nh] = -ntri/2 - 1;
          for (n = 0; n < len-1; n++) {
            segs[ntri] = n+nseg+1;
            ntri++;
            segs[ntri] = n+nseg+2;
            ntri++;
          }
          if (esenses[k] ==  1) heads[nh] = ntri/2;
          for (n = 0; n < len; n++) {
            lsegs[3*nseg  ] = xyzs[3*n  ];
            lsegs[3*nseg+1] = xyzs[3*n+1];
            lsegs[3*nseg+2] = xyzs[3*n+2];
            nseg++;
          }
          nh++;
        }
        sprintf(gpname, "Body %d Loop %d/%d", ibody+1, i+1, j+1);
        stat = wv_setData(WV_REAL32, nseg, (void *) lsegs,  WV_VERTICES, &items[0]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
        wv_adjustVerts(&items[0], focus);
        free(lsegs);
        stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[1]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
        stat = wv_setData(WV_INT32, ntri, (void *) segs, WV_INDICES, &items[2]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
        free(segs);
        stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
        if (stat < 0) {
          printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
        } else {
          if (cntxt != NULL)
            if (cntxt->gPrims != NULL) {
              cntxt->gPrims[stat].lWidth = 1.0;
              n = wv_addArrowHeads(cntxt, stat, 0.05, nh, heads);
              if (n != 0) printf(" wv_addArrowHeads = %d\n", n);
            }
        }
        free(heads);
/*      wv_printGPrim(cntxt, stat);  */
      }
    }
  }
  printf(" ** %d gPrims with %d triangles **\n", stat+1, sum);

  /* start the server code */

  stat = 0;
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    /* we have a single valid server -- stay alive a long as we have a client */
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
/*    wv_broadcastText("I'm here!");  */
    }
  }
  wv_cleanupServers();

  /* finish up */
  for (ibody = 0; ibody < nbody; ibody++) {
    EG_deleteObject(bodydata[ibody].tess);
    EG_free(bodydata[ibody].edges);
    EG_free(bodydata[ibody].faces);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}


/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, char *text, /*@unused@*/ int lena)
{
  int          i, j, k, m, n, ibody, stat, sum, len, ntri, index, oclass, mtype;
  int          nloops, nseg, nledges, *segs, *lsenses, *esenses;
  int          nh, *heads;
  const int    *tris, *tric, *ptype, *pindex;
  char         gpname[43];
  float        *lsegs, *colrs, val, color[3];
  double       rc[8], rmin, rmax, result[18], dot, norm[3], *norms, *pu, *pv;
  const double *xyzs, *uvs, *ts;
  ego          geom, *ledges, *loops;
  wvData       items[4];

  printf(" RX: %s\n", text);

  if ((strcmp(text,"finer") != 0) && (strcmp(text,"coarser") != 0) &&
      (strcmp(text,"next")  != 0) && (strcmp(text,"limits")  != 0)) return;

  /* just change the color mapping */
  if ((strcmp(text,"next") == 0) || (strcmp(text,"limits") == 0)) {
    if (strcmp(text,"next") == 0) {
      key++;
      if (key > 2) key = 0;
    } else {
      printf(" Enter new limits [old = %e, %e]:", lims[0], lims[1]);
      scanf("%f %f", &lims[0], &lims[1]);
      printf(" new limits = %e %e\n", lims[0], lims[1]);
    }
    stat = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], keys[key]);
    if (stat < 0) printf(" wv_setKey = %d!\n", stat);

    for (ibody = 0; ibody < nbody; ibody++) {
      /* get faces */
      for (i = 0; i < bodydata[ibody].nfaces; i++) {
        stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                              &xyzs, &uvs, &ptype, &pindex, &ntri,
                              &tris, &tric);
        if (stat != EGADS_SUCCESS) continue;
        sprintf(gpname, "Body %d Face %d", ibody+1, i+1);
        index = wv_indexGPrim(cntxt, gpname);
        if (index < 0) {
          printf(" wv_indexGPrim = %d for %s (%d)!\n", i, gpname, index);
          continue;
        }
        if (len == 0) continue;
        colrs = (float *) malloc(3*len*sizeof(float));
        if (colrs == NULL) continue;
        for (j = 0; j < len; j++) {
          rc[0] = rc[4] = 0.0;
          stat = EG_curvature(bodydata[ibody].faces[i], &uvs[2*j], rc);
          if (stat != EGADS_SUCCESS)
            printf(" Face %d: %d EG_curvature = %d\n", i+1, j, stat);
          if (key == 0) {
            val = 0.0;
            if (rc[0]*rc[4] > 0.0) val =  1.0;
            if (rc[0]*rc[4] < 0.0) val = -1.0;
            val *= pow(fabs(rc[0]*rc[4]*focus[3]*focus[3]),0.25);
          } else if (key == 1) {
            val = rc[0];
            if (rc[4] < rc[0]) val = rc[4];
            val *= focus[3];
          } else {
            val = rc[0];
            if (rc[4] > rc[0]) val = rc[4];
            val *= focus[3];
          }
          spec_col(val, &colrs[3*j]);
        }
        stat = wv_setData(WV_REAL32, len, (void *) colrs, WV_COLORS, &items[0]);
        if (stat < 0) printf(" wv_setData = %d for %s/item color!\n", i, gpname);
        free(colrs);
        stat = wv_modGPrim(cntxt, index, 1, items);
        if (stat < 0)
          printf(" wv_modGPrim = %d for %s (%d)!\n", stat, gpname, index);
      }
    }
    return;
  }

  if (strcmp(text,"finer")   == 0) params[0] *= 0.5;
  if (strcmp(text,"coarser") == 0) params[0] *= 2.0;
  printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf,  key = %d\n",
         params[2], params[0], params[1], key);
  if (key == -1) {
    key  = 0;
    stat = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], keys[key]);
    if (stat < 0) printf(" wv_setKey = %d!\n", stat);
  }

  for (ibody = 0; ibody < nbody; ibody++) {
    EG_deleteObject(bodydata[ibody].tess);
    bodydata[ibody].tess = NULL;
    stat = EG_makeTessBody(bodydata[ibody].body, params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS)
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
  }

  /* make the scene */
  for (sum = stat = ibody = 0; ibody < nbody; ibody++) {

    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      sprintf(gpname, "Body %d Face %d", ibody+1, i+1);
      index = wv_indexGPrim(cntxt, gpname);
      if (index < 0) {
        printf(" wv_indexGPrim = %d for %s (%d)!\n", i, gpname, index);
        continue;
      }
      stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      for (nseg = j = 0; j < ntri; j++)
        for (k = 0; k < 3; k++)
          if (tric[3*j+k] < j+1) nseg++;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      for (nseg = j = 0; j < ntri; j++)
        for (k = 0; k < 3; k++)
          if (tric[3*j+k] < j+1) {
            segs[2*nseg  ] = tris[3*j+sides[k][0]];
            segs[2*nseg+1] = tris[3*j+sides[k][1]];
            nseg++;
          }
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(segs);
      colrs = (float *) malloc(3*len*sizeof(float));
      if (colrs == NULL) continue;
      for (j = 0; j < len; j++) {
        rc[0] = rc[4] = 0.0;
        stat = EG_curvature(bodydata[ibody].faces[i], &uvs[2*j], rc);
        if (stat != EGADS_SUCCESS)
          printf(" Face %d: %d EG_curvature = %d\n", i+1, j, stat);
        if (key == 0) {
          val  = 0.0;
          rmin = rmax = fabs(rc[0]);
          if (fabs(rc[4]) < rmin) rmin = fabs(rc[4]);
          if (fabs(rc[4]) > rmax) rmax = fabs(rc[4]);
          if (rmax != 0.0) {
            if (rc[0]*rc[4] > 0.0) val =  1.0;
            if (rc[0]*rc[4] < 0.0) val = -1.0;
            if (rmin/rmax < 1.e-5) val =  0.0;
          }
          val *= pow(fabs(rc[0]*rc[4]*focus[3]*focus[3]),0.25);
        } else if (key == 1) {
          val = rc[0];
          if (rc[4] < rc[0]) val = rc[4];
          val *= focus[3];
        } else {
          val = rc[0];
          if (rc[4] > rc[0]) val = rc[4];
          val *= focus[3];
        }
        spec_col(val, &colrs[3*j]);
      }
      stat = wv_setData(WV_REAL32, len, (void *) colrs, WV_COLORS, &items[3]);
      if (stat < 0) printf(" wv_setData = %d for %s/item color!\n", i, gpname);
      free(colrs);
      stat = wv_modGPrim(cntxt, index, 4, items);
      if (stat < 0)
        printf(" wv_modGPrim = %d for %s (%d)!\n", stat, gpname, index);
      sum += ntri;
    }

    /* put normals of faces in "edge" slot */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
      norms = (double *) malloc(6*len*sizeof(double));
      if (norms == NULL) continue;
      sprintf(gpname, "Body %d Edge %d", ibody+1, i+1);

      pu = &result[3];
      pv = &result[6];
      for (j = 0; j < len; j++) {
        norms[6*j  ] = norms[6*j+3] = xyzs[3*j  ];
        norms[6*j+1] = norms[6*j+4] = xyzs[3*j+1];
        norms[6*j+2] = norms[6*j+5] = xyzs[3*j+2];
        stat = EG_evaluate(bodydata[ibody].faces[i], &uvs[2*j], result);
        if (stat != EGADS_SUCCESS) continue;
        CROSS(norm, pu, pv);
        dot = sqrt(DOT(norm, norm))*bodydata[ibody].faces[i]->mtype;
        if (dot == 0.0) continue;
        norm[0]      *= 0.025*focus[3]/dot;
        norm[1]      *= 0.025*focus[3]/dot;
        norm[2]      *= 0.025*focus[3]/dot;
        norms[6*j+3] += norm[0];
        norms[6*j+4] += norm[1];
        norms[6*j+5] += norm[2];
      }
      stat = wv_setData(WV_REAL64, 2*len, (void *) norms,  WV_VERTICES,
                        &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      free(norms);
      color[0] = 0.0;
      color[1] = 0.0;
      color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, 0, 2, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
    }

    /* get loops */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTopology(bodydata[ibody].faces[i], &geom, &oclass,
                            &mtype, NULL, &nloops, &loops, &lsenses);
      if (stat != EGADS_SUCCESS) continue;
      for (nh = j = 0; j < nloops; j++) {
        stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, NULL,
                              &nledges, &ledges, &esenses);
        if (stat != EGADS_SUCCESS) continue;

        /* count */
        for (ntri = nseg = k = 0; k < nledges; k++) {
          m = 0;
          while (ledges[k] != bodydata[ibody].edges[m]) {
            m++;
            if (m == bodydata[ibody].nedges) break;
          }
          /* assume that the edge is degenerate and removed */
          if (m == bodydata[ibody].nedges) continue;
          stat = EG_getTessEdge(bodydata[ibody].tess, m+1, &len,
                                &xyzs, &ts);
          if (stat != EGADS_SUCCESS) {
            printf(" EG_getTessEdge %d = %d!\n", m+1, stat);
            nseg = 0;
            break;
          }
          if (len == 2)
            if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
                (xyzs[2] == xyzs[5])) continue;
          nh++;
          nseg += len;
          ntri += 2*(len-1);
        }
        if (nseg == 0) continue;
        lsegs = (float *) malloc(3*nseg*sizeof(float));
        if (lsegs == NULL) {
          printf(" Can not allocate %d Segments!\n", nseg);
          continue;
        }
        segs = (int *) malloc(ntri*sizeof(int));
        if (segs == NULL) {
          printf(" Can not allocate %d Line Segments!\n", ntri);
          free(lsegs);
          continue;
        }
        heads = (int *) malloc(nh*sizeof(int));
        if (heads == NULL) {
          printf(" Can not allocate %d Heads!\n", nh);
          free(segs);
          free(lsegs);
          continue;
        }

        /* fill */
        for (nh = ntri = nseg = k = 0; k < nledges; k++) {
          m = 0;
          while (ledges[k] != bodydata[ibody].edges[m]) {
            m++;
            if (m == bodydata[ibody].nedges) break;
          }
          /* assume that the edge is degenerate and removed */
          if (m == bodydata[ibody].nedges) continue;
          EG_getTessEdge(bodydata[ibody].tess, m+1, &len, &xyzs, &ts);
          if (len == 2)
            if ((xyzs[0] == xyzs[3]) && (xyzs[1] == xyzs[4]) &&
                (xyzs[2] == xyzs[5])) continue;
          if (esenses[k] == -1) heads[nh] = -ntri/2 - 1;
          for (n = 0; n < len-1; n++) {
            segs[ntri] = n+nseg+1;
            ntri++;
            segs[ntri] = n+nseg+2;
            ntri++;
          }
          if (esenses[k] ==  1) heads[nh] = ntri/2;
          for (n = 0; n < len; n++) {
            lsegs[3*nseg  ] = xyzs[3*n  ];
            lsegs[3*nseg+1] = xyzs[3*n+1];
            lsegs[3*nseg+2] = xyzs[3*n+2];
            nseg++;
          }
          nh++;
        }
        sprintf(gpname, "Body %d Loop %d/%d", ibody+1, i+1, j+1);
        index = wv_indexGPrim(cntxt, gpname);
        if (index < 0) {
          printf(" wv_indexGPrim = %d for %s (%d)!\n", i, gpname, index);
          continue;
        }
        stat = wv_setData(WV_REAL32, nseg, (void *) lsegs,  WV_VERTICES, &items[0]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
        wv_adjustVerts(&items[0], focus);
        free(lsegs);
        stat = wv_setData(WV_INT32, ntri, (void *) segs, WV_INDICES, &items[1]);
        if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
        free(segs);
        stat = wv_modGPrim(cntxt, index, 2, items);
        if (stat < 0) {
          printf(" wv_modGPrim = %d for %s!\n", stat, gpname);
        } else {
          n = wv_addArrowHeads(cntxt, index, 0.05, nh, heads);
          if (n != 0) printf(" wv_addArrowHeads = %d\n", n);
        }
        free(heads);
      }
    }
  }
  printf(" **  now with %d triangles **\n\n", sum);
}
