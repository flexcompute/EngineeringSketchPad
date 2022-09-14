/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Effective Topo Tessellation using wv
 *
 *      Copyright 2011-2022, Massachusetts Institute of Technology
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

//#define DISJOINTQUADS
//#define NOEFFECT
#define COLORING
#define SAVEMODEL


  static wvContext *cntxt;
  static ego       ebody, body, tess = NULL;
  static float     focus[4];
#ifdef COLORING
  static int       key = -1;
  static float     lims[2] = {-1.0, +1.0};
  static char      *keys[8] = {"U", "V", "dX/dU", "dY/dU", "dZ/dU",
                                         "dX/dV", "dY/dV", "dZ/dV"};


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
#endif


static void
parseOut(int level, ego object, int sense)
{
  int    i, j, stat, oclass, mtype, nobjs, periodic, nedge, nface;
  int    *senses, *ivec;
  ego    prev, next, geom, top, *objs, *edges, *faces;
  double limits[4], *rvec;
  static char *classType[36] = {"CONTEXT", "TRANSFORM", "TESSELLATION",
                                "NIL", "EMPTY", "REFERENCE", "", "",
                                "", "", "PCURVE", "CURVE", "SURFACE", "",
                                "", "", "", "", "", "", "NODE",
                                "EGDE", "LOOP", "FACE", "SHELL",
                                "BODY", "MODEL", "", "", "", "", "EEDGE",
                                "ELOOP", "EFACE", "ESHELL", "EBODY"};
  static char *curvType[9] = {"Line", "Circle", "Ellipse", "Parabola",
                              "Hyperbola", "Trimmed", "Bezier", "BSpline",
                              "Offset"};
  static char *surfType[11] = {"Plane", "Spherical", "Cylinder", "Revolution",
                               "Toroidal", "Trimmed" , "Bezier", "BSpline",
                               "Offset", "Conical", "Extrusion"};
 
#if defined(WIN32) && defined(_OCC64)
  long long pointer;
  
  pointer = (long long) object;
#else
  long pointer;
  
  pointer = (long) object;
#endif

  stat = EG_getInfo(object, &oclass, &mtype, &top, &prev, &next);
  if (stat != EGADS_SUCCESS) {
    printf(" parseOut: %d EG_getInfo return = %d\n", level, stat);
    return;
  }
  
  /* geometry */
  if ((oclass >= PCURVE) && (oclass <= SURFACE)) {
    stat = EG_getGeometry(object, &oclass, &mtype, &geom, &ivec, &rvec);
    if (stat != EGADS_SUCCESS) {
      printf(" parseOut: %d EG_getGeometry return = %d\n", level, stat);
      return;
    }
    stat = EG_getRange(object, limits, &periodic);

    /* output most data except the axes info */
    if (oclass != SURFACE) {

      for (i = 0; i < 2*level; i++) printf(" ");
#if defined(WIN32) && defined(_OCC64)
      printf("%s %llx  range = %le %le  per = %d\n",
             classType[oclass], pointer, limits[0], limits[1], periodic);
#else
      printf("%s %lx  range = %le %le  per = %d\n",
             classType[oclass], pointer, limits[0], limits[1], periodic);
#endif

      for (i = 0; i < 2*level+2; i++) printf(" ");
      if (oclass == PCURVE) {
        switch (mtype) {
          case CIRCLE:
            printf("%s  radius = %lf\n", curvType[mtype-1], rvec[6]);
            break;
          case ELLIPSE:
            printf("%s  major = %lf, minor = %lf\n",
                   curvType[mtype-1], rvec[6], rvec[7]);
            break;
          case PARABOLA:
            printf("%s  focus = %lf\n", curvType[mtype-1], rvec[6]);
            break;
          case HYPERBOLA:
            printf("%s  major = %lf, minor = %lf\n",
                   curvType[mtype-1], rvec[6], rvec[7]);
            break;
          case TRIMMED:
            printf("%s  first = %lf, last = %lf\n",
                   curvType[mtype-1], rvec[0], rvec[1]);
            break;
          case BEZIER:
            printf("%s  flags = %x, degree = %d, #CPs = %d\n",
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2]);
            break;
          case BSPLINE:
            printf("%s  flags = %x, degree = %d, #CPs = %d, #knots = %d\n",
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3]);
/*          for (i = 0; i < 2*level+2; i++) printf(" ");
            printf("knots =");
            for (i = 0; i < ivec[3]; i++) printf(" %lf", rvec[i]);
            printf("\n");  */
            break;
          case OFFSET:
            printf("%s  offset = %lf\n", curvType[mtype-1], rvec[0]);
            break;
          case 0:
            printf("unknown curve type!\n");
            break;
          default:
            printf("%s   %lf %lf   %lf %lf\n", curvType[mtype-1],
                   rvec[0], rvec[1], rvec[2], rvec[3]);
        }
      } else {
        switch (mtype) {
          case CIRCLE:
            printf("%s  radius = %lf\n", curvType[mtype-1], rvec[9]);
            break;
          case ELLIPSE:
            printf("%s  major = %lf, minor = %lf\n",
                   curvType[mtype-1], rvec[9], rvec[10]);
            break;
          case PARABOLA:
            printf("%s  focus = %lf\n", curvType[mtype-1], rvec[9]);
            break;
          case HYPERBOLA:
            printf("%s  major = %lf, minor = %lf\n",
                   curvType[mtype-1], rvec[9], rvec[10]);
            break;
          case TRIMMED:
            printf("%s  first = %lf, last = %lf\n",
                   curvType[mtype-1], rvec[0], rvec[1]);
            break;
          case BEZIER:
            printf("%s  flags = %x, degree = %d, #CPs = %d\n",
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2]);
            break;
          case BSPLINE:
            printf("%s  flags = %x, degree = %d, #CPs = %d, #knots = %d\n",
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3]);
/*          for (i = 0; i < 2*level+2; i++) printf(" ");
            printf("knots =");
            for (i = 0; i < ivec[3]; i++) printf(" %lf", rvec[i]);
            printf("\n");  */
            break;
          case OFFSET:
            printf("%s  offset = %lf\n", curvType[mtype-1], rvec[3]);
            break;
          case 0:
            printf("unknown curve type!\n");
            break;
          default:
            printf("%s\n", curvType[mtype-1]);
        }
      }

    } else {
      for (i = 0; i < 2*level; i++) printf(" ");
#if defined(WIN32) && defined(_OCC64)
      printf("%s %llx  Urange = %le %le  Vrange = %le %le  per = %d\n",
             classType[oclass], pointer, limits[0], limits[1],
                                         limits[2], limits[3], periodic);
#else
      printf("%s %lx  Urange = %le %le  Vrange = %le %le  per = %d\n",
             classType[oclass], pointer, limits[0], limits[1],
                                         limits[2], limits[3], periodic);
#endif

      for (i = 0; i < 2*level+2; i++) printf(" ");
      switch (mtype) {
        case SPHERICAL:
          printf("%s  radius = %lf\n", surfType[mtype-1], rvec[9]);
          break;
        case CONICAL:
          printf("%s  angle = %lf, radius = %lf\n",
                 surfType[mtype-1], rvec[12], rvec[13]);
          printf("    rvec = %lf %lf %lf   %lf %lf %lf  \n",
                 rvec[0], rvec[1], rvec[2], rvec[3], rvec[4],  rvec[5]);
          printf("           %lf %lf %lf   %lf %lf %lf  \n",
                 rvec[6], rvec[7], rvec[8], rvec[9], rvec[10], rvec[11]);
          break;
        case CYLINDRICAL:
          printf("%s  radius = %lf\n", surfType[mtype-1], rvec[12]);
          break;
        case TOROIDAL:
          printf("%s  major = %lf, minor = %lf\n",
                 surfType[mtype-1], rvec[12], rvec[13]);
          break;
        case BEZIER:
          printf("%s  flags = %x, U deg = %d #CPs = %d, V deg = %d #CPs = %d\n",
                 surfType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3], ivec[4]);
          break;
        case BSPLINE:
          printf("%s  flags = %x, U deg = %d #CPs = %d #knots = %d ",
                 surfType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3]);
          printf(" V deg = %d #CPs = %d #knots = %d\n",
                 ivec[4], ivec[5], ivec[6]);
/*        for (i = 0; i < 2*level+2; i++) printf(" ");
          printf("Uknots =");
          for (i = 0; i < ivec[3]; i++) printf(" %lf", rvec[i]);
          for (i = 0; i < 2*level+2; i++) printf(" ");
          printf("\nVknots =");
          for (i = 0; i < ivec[6]; i++) printf(" %lf", rvec[i+ivec[3]]);
          printf("\n"); */
          break;
        case TRIMMED:
          printf("%s  U trim = %lf %lf, V trim = %lf %lf\n",
                 surfType[mtype-1], rvec[0], rvec[1], rvec[2], rvec[3]);
          break;
        case OFFSET:
          printf("%s  offset = %lf\n", surfType[mtype-1], rvec[0]);
          break;
        case 0:
          printf("unknown surface type!\n");
          break;
        default:
          printf("%s\n", surfType[mtype-1]);
      }
    }

    if (ivec != NULL) EG_free(ivec);
    if (rvec != NULL) EG_free(rvec);
    if (geom != NULL) parseOut(level+1, geom, 0);
    return;
  }

  /* output class and pointer data */

  for (i = 0; i < 2*level; i++) printf(" ");
  if (sense == 0) {
#if defined(WIN32) && defined(_OCC64)
    printf("%s %llx  mtype = %d\n", classType[oclass], pointer, mtype);
#else
    printf("%s %lx  mtype = %d\n",  classType[oclass], pointer, mtype);
#endif
  } else {
#if defined(WIN32) && defined(_OCC64)
    printf("%s %llx  sense = %d\n", classType[oclass], pointer, sense);
#else
    printf("%s %lx  sense = %d\n",  classType[oclass], pointer, sense);
#endif
  }

  /* topology*/
  if ((oclass >= NODE) && (oclass <= MODEL)) {
    stat = EG_getTopology(object, &geom, &oclass, &mtype, limits, &nobjs,
                          &objs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" parseOut: %d EG_getTopology return = %d\n", level, stat);
      return;
    }
    if (oclass == NODE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      stat = EG_indexBodyTopo(body, object);
      printf("%3d  XYZ = %lf %lf %lf\n", stat, limits[0], limits[1], limits[2]);
    } else if (oclass == EDGE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      if (mtype == DEGENERATE) {
        printf("tRange = %lf %lf -- Degenerate!\n", limits[0], limits[1]);
      } else {
        printf("tRange = %lf %lf\n", limits[0], limits[1]);
      }
    } else if (oclass == FACE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      printf("uRange = %lf %lf, vRange = %lf %lf\n", limits[0], limits[1],
                                                     limits[2], limits[3]);
    }
    
    if ((geom != NULL) && (mtype != DEGENERATE)) parseOut(level+1, geom, 0);
    if (senses == NULL) {
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], 0);
    } else {
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], senses[i]);
    }
    if ((geom != NULL) && (oclass == LOOP))
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i+nobjs], 0);

  } else {
    /* effective topology */
    stat = EG_getTopology(object, &geom, &oclass, &mtype, limits, &nobjs,
                          &objs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" parseOut: %d EG_getTopology return = %d\n", level, stat);
      return;
    }
    if (oclass == EEDGE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      if (mtype == DEGENERATE) {
        printf("tRange = %lf %lf -- Degenerate!\n", limits[0], limits[1]);
      } else {
        printf("tRange = %lf %lf -- mtype = %d\n", limits[0], limits[1], mtype);
      }
      stat = EG_effectiveEdgeList(object, &nedge, &edges, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) {
        printf(" parseOut: %d EG_effectiveEdgeList return = %d\n", level, stat);
        return;
      }
      for (j = 0; j < nedge; j++) {
        for (i = 0; i < 2*level+2; i++) printf(" ");
        printf("tStart = %lf  sense = %d\n", rvec[j], ivec[j]);
        parseOut(level+1, edges[j], 0);
      }
      EG_free(edges);
      EG_free(ivec);
      EG_free(rvec);
    } else if (oclass == EFACE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      printf("uRange = %lf %lf, vRange = %lf %lf\n", limits[0], limits[1],
                                                     limits[2], limits[3]);
      stat = EG_getBodyTopos(ebody, object, FACE, &nface, &faces);
      if (stat != EGADS_SUCCESS) {
        printf(" parseOut: %d EG_getBodyTopos return = %d\n", level, stat);
        return;
      }
      for (j = 0; j < nface; j++) parseOut(level+1, faces[j], 0);
      EG_free(faces);
    }
    
    if (senses == NULL) {
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], 0);
    } else {
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], senses[i]);
    }
  }

}


/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *uPtr, /*@unused@*/ void *wsi, char *text,
                    /*@unused@*/ int lena)
{
  int          i, index, stat, ient, nnode;
  char         tag[5];
  double       d, dist, coord[3], data[3];
  ego          eface, *nodes;
#ifdef COLORING
  char         gpname[43];
  int          j, len, ntri, nface;
  float        val, *colrs;
  double       result[18];
  ego          *faces;
  const int    *tris, *tric, *ptype, *pindex;
  const double *xyzs, *uvs;
  wvData       items[1];
#endif
  
  if (strncmp(text,"Located: ", 9) == 0) {
    sscanf(&text[8], "%lf %lf %lf", &coord[0], &coord[1], &coord[2]);
    dist     = 1.e200;
    coord[0] = coord[0]*focus[3] + focus[0];
    coord[1] = coord[1]*focus[3] + focus[1];
    coord[2] = coord[2]*focus[3] + focus[2];
    printf(" Closest Node to %lf %lf %lf:\n", coord[0], coord[1], coord[2]);
    stat = EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos = %d\n", stat);
      return;
    }
    for (index = i = 0; i < nnode; i++) {
      stat = EG_evaluate(nodes[i], NULL, data);
      if (stat != EGADS_SUCCESS) continue;
      d    = sqrt((data[0]-coord[0])*(data[0]-coord[0]) +
                  (data[1]-coord[1])*(data[1]-coord[1]) +
                  (data[2]-coord[2])*(data[2]-coord[2]));
      if (d < dist) {
        index = i+1;
        dist  = d;
      }
    }
    EG_free(nodes);
    printf(" Nearest Node = %d  dist = %le\n", index, dist);
    return;
  }
  if (strncmp(text,"Picked: ", 8) == 0) {
    sscanf(&text[12], "%d %s %d", &i, tag, &ient);
    
    printf(" Picked: %s %d\n", tag, ient);
    if (strcmp(tag, "Face") == 0) {
      stat = EG_objectBodyTopo(ebody, EFACE, ient, &eface);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_objectBodyTopo = %d\n", stat);
        return;
      }
      parseOut(0, eface, 0);
      printf("\n");
    }
    return;
  }
  
#ifdef COLORING
  /* just change the color mapping */
  if ((strcmp(text,"next") == 0) || (strcmp(text,"limits") == 0)) {
    if (strcmp(text,"next") == 0) {
      key++;
      if (key > 7) key = 0;
    } else {
      printf(" Enter new limits [old = %e, %e]:", lims[0], lims[1]);
      scanf("%f %f", &lims[0], &lims[1]);
      printf(" new limits = %e %e\n", lims[0], lims[1]);
    }
    stat = wv_setKey(cntxt, 256, color_map, lims[0], lims[1], keys[key]);
    if (stat < 0) printf(" wv_setKey = %d!\n", stat);
    EG_getBodyTopos(ebody, NULL, EFACE, &nface, &faces);

    for (i = 0; i < nface; i++) {
      stat = EG_getTessFace(tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getTessFace %d/%d = %d\n", i+1, nface, stat);
        continue;
      }
      sprintf(gpname, "Body %d Face %d", 1, i+1);
      index = wv_indexGPrim(cntxt, gpname);
      if (index < 0) {
        printf(" wv_indexGPrim = %d for %s (%d)!\n", i, gpname, index);
        continue;
      }
      if (len == 0) continue;
      colrs = (float *) malloc(3*len*sizeof(float));
      if (colrs == NULL) continue;
      for (j = 0; j < len; j++) {
        EG_evaluate(faces[i], &uvs[2*j], result);
        if (key == 0) {
          val = uvs[2*j  ];
        } else if (key == 1) {
          val = uvs[2*j+1];
        } else {
          val = result[key+1];
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
    EG_free(faces);
    return;
  }
#endif

  printf(" browserMessage = %s\n", text);
}


int main(int argc, char *argv[])
{
  int          i, j, k, ibody, stat, oclass, mtype, len, ntri, ndum, sum;
  int          nenode, needge, neface, nface, nseg, *segs, *senses;
  int          nbody, nedge, nnode, ngp, atype, alen, quad;
  const int    *tris, *tric, *ptype, *pindex, *ints;
  float        arg, color[3];
  double       box[6], maxdev, size, tol, *realx, params[3], uv[2], result[18];
  const double *xyzs, *uvs, *ts, *reals;
  char         gpname[34], *startapp;
  const char   *OCCrev, *string;
  ego          context, model, geom, obj, *bodies, *dum, *faces, *edges, *nodes;
#ifdef DISJOINTQUADS
  ego          quads;
#endif
#ifdef SAVEMODEL
  ego          save[3];
#endif
  wvData       items[5];
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = {  0,     3,     4,     2};

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vEffect filename [angle maxlen sag]\n\n");
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
  stat = EG_copyObject(bodies[0], NULL, &body);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_copyObject = %d\n", stat);
    return 1;
  }
  for (i = nbody; i < mtype; i++) {
    if (bodies[i]->oclass != TESSELLATION) continue;
    stat = EG_statusTessBody(bodies[i], &geom, &j, &len);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_statusTessBody %d = %d\n", i, stat);
      continue;
    }
    if (geom != bodies[0]) continue;
    printf(" Found Tessellation for first Body @ %d!\n", i);
    stat = EG_copyObject(bodies[i], body, &tess);
    if (stat != EGADS_SUCCESS)
      printf(" EG_copyObject %d = %d\n", i, stat);
    break;
  }
  EG_deleteObject(model);

  params[0] =  0.025*size;
  params[1] =  0.005*size;
  params[2] = 15.0;
#ifdef DISJOINTQUADS
  params[0] =  0.100*size;
  params[1] =  0.010*size;
  params[2] = 24.0;
#endif
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
  printf(" Reference size = %le\n", size);
  
  mtype = 0;
  EG_getTopology(body, &geom, &oclass, &mtype, NULL, &j, &dum, &senses);
  stat = EG_tolerance(body, &tol);
  if (mtype == WIREBODY) {
    printf(" Body Type = WireBody   tol = %le\n", tol);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  } else if (mtype == FACEBODY) {
    printf(" Body Type = FaceBody   tol = %le\n", tol);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  } else if (mtype == SHEETBODY) {
    printf(" Body Type = SheetBody  tol = %le\n", tol);
  } else {
    printf(" Body Type = SolidBody  tol = %le\n", tol);
  }
  EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
  EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges);
  EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes);
  printf("          nNode = %d   nEdge = %d   nFace = %d\n\n",
         nnode, nedge, nface);
  
  /* mark Node/Edges to ".Keep" */
  do {
    printf(" Enter Keep Index [-Node/+Edge]: ");
    scanf("%d", &ndum);
    if (ndum == 0) break;
    if (ndum > 0) {
      if (ndum > nedge) {
        printf("   Error -- Edge index too big [1-%d]\n", nedge);
        continue;
      }
      obj = edges[ndum-1];
    } else {
      ndum = -ndum;
      if (ndum > nnode) {
        printf("   Error -- Node index too big [1-%d]\n", nnode);
        continue;
      }
      obj = nodes[ndum-1];
    }
    stat = EG_attributeAdd(obj, ".Keep", ATTRSTRING, 2, NULL, NULL, ".");
    if (stat != EGADS_SUCCESS)
      printf("   EG_attributeAdd = %d\n", stat);
  } while (ndum != 0);
  EG_free(nodes);
  EG_free(edges);
  printf(" \n");
  
  /* do our first tessellation */
  if (tess == NULL) {
    stat = EG_makeTessBody(body, params, &tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody  = %d\n", stat);
      EG_deleteObject(body);
      EG_close(context);
      exit(1);
    }
  }

#ifdef NOEFFECT
  EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes);
#else
  
  /* start the EBody */
  stat = EG_initEBody(tess, 3.0, &ebody);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_initEBody     = %d\n", stat);
    EG_deleteObject(tess);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  }
  
  /* make some composites! */
  stat = EG_getBodyTopos(ebody, NULL, FACE, &nface, &faces);
  do {
    EG_getBodyTopos(ebody, NULL, EFACE, &neface, NULL);
    EG_getBodyTopos(ebody, NULL, EEDGE, &needge, NULL);
    EG_getBodyTopos(ebody, NULL,  NODE, &nenode, NULL);
    printf("          nNode = %d  nEEdge = %d  nEFace = %d\n\n",
           nenode, needge, neface);
    printf(" Enter number of Faces: ");
    scanf("%d", &ndum);
    if (ndum <= 0) break;
    dum = (ego *) EG_alloc(ndum*sizeof(ego));
    if (dum == NULL) {
      printf("\n Error: Malloc of %d egos!\n", ndum);
      break;
    }
    if (ndum == nface) {
      for (i = 0; i < ndum; i++) dum[i] = faces[i];
    } else {
      printf(" Enter Faces Indices: ");
      for (i = 0; i < ndum; i++) {
        scanf("%d", &j);
        if ((j < 1) || (j > nface)) {
          printf("\n Error: Bad Index = %d [1-%d]\n", j, nface);
          break;
        }
        dum[i] = faces[j-1];
      }
    }
    stat = EG_makeEFace(ebody, ndum, dum, &geom);
    if (stat != EGADS_SUCCESS) printf(" EG_makeEFace = %d\n", stat);
    EG_free(dum);
    if (ndum == nface) break;
  } while (stat == EGADS_SUCCESS);
  EG_free(faces);

  /* close the EBody */
  stat = EG_finishEBody(ebody);
  EG_deleteObject(tess);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_finalize      = %d\n", stat);
    EG_deleteObject(ebody);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  }

  stat = EG_getBodyTopos(ebody, NULL, EFACE, &nface, &faces);
  i    = EG_getBodyTopos(ebody, NULL, EEDGE, &nedge, &edges);
  j    = EG_getBodyTopos(ebody, NULL, NODE,  &nnode, &nodes);
  if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
    printf(" EG_getBodyTopos Face = %d\n", stat);
    printf(" EG_getBodyTopos Edge = %d\n", i);
    EG_deleteObject(ebody);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  }
/*  printf("           nEFaces = %d   nEEdges = %d\n", nface, nedge);  */

  /* tessellate the EBody with the same parameters */
  stat = EG_makeTessBody(ebody, params, &tess);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTessBody = %d\n", stat);
    EG_free(nodes);
    EG_deleteObject(ebody);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  }
  maxdev = 0.0;
  for (i = 0; i < nedge; i++) {
    stat = EG_getTessEdge(tess, i+1, &len, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) continue;
    box[0] = box[1] = -1.0;
    stat = EG_evaluate(edges[i], &ts[0], result);
    if (stat == EGADS_DEGEN) continue;
    if (stat == EGADS_SUCCESS)
      box[0] = sqrt((xyzs[0]-result[0])*(xyzs[0]-result[0]) +
                    (xyzs[1]-result[1])*(xyzs[1]-result[1]) +
                    (xyzs[2]-result[2])*(xyzs[2]-result[2]));
    stat = EG_evaluate(edges[i], &ts[len-1], result);
    if (stat == EGADS_SUCCESS)
      box[1] = sqrt((xyzs[3*len-3]-result[0])*(xyzs[3*len-3]-result[0]) +
                    (xyzs[3*len-2]-result[1])*(xyzs[3*len-2]-result[1]) +
                    (xyzs[3*len-1]-result[2])*(xyzs[3*len-1]-result[2]));
    if ((box[0] < 0.0) || (box[1] < 0.0))
      printf(" EEdge %3d: %le %le\n", i+1, box[0], box[1]);
    if (box[0] > maxdev) maxdev = box[0];
    if (box[1] > maxdev) maxdev = box[1];
/*  if (box[0] > 1.e-8)
      printf(" EEdge = %3d  first  dist = %le\n", i+1, box[0]);
    if (box[1] > 1.e-8)
      printf(" EEdge = %3d  second dist = %le\n", i+1, box[1]);  */
  }
  printf("\n EEdge/Node  deviation = %le\n", maxdev);
  
  tol = maxdev = 0.0;
  for (i = 0; i < nface; i++) {
    stat = EG_getTessFace(tess, i+1, &len, &xyzs, &uvs, &ptype, &pindex, &ntri,
                          &tris, &tric);
    if (stat != EGADS_SUCCESS) continue;
    for (j = 0; j < len; j++) {
      if (ptype[j] <  0) continue;
      if (ptype[j] == 0) {
        stat = EG_evaluate(faces[i], &uvs[2*j], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EFace %3d: EG_evaluate stat = %d\n", i+1, stat);
          if (stat != EGADS_EXTRAPOL) continue;
        }
        box[0] = sqrt((xyzs[3*j  ]-result[0])*(xyzs[3*j  ]-result[0]) +
                      (xyzs[3*j+1]-result[1])*(xyzs[3*j+1]-result[1]) +
                      (xyzs[3*j+2]-result[2])*(xyzs[3*j+2]-result[2]));
        if (box[0] > maxdev) maxdev = box[0];
        continue;
      }
      
      uv[0] = uvs[2*j  ];
      uv[1] = uvs[2*j+1];
/*
      stat = EG_getTessEdge(tess, pindex[j], &alen, &reals, &ts);
      if (stat != EGADS_SUCCESS) {
        printf(" EFace %3d: EG_getTessEdge %d stat = %d\n",
               i+1, pindex[j], stat);
        continue;
      }
      stat = EG_getEdgeUV(faces[i], edges[pindex[j]-1], 0, ts[ptype[j]-1], uv);
      if (stat != EGADS_SUCCESS) {
        printf(" EFace %3d: EG_getEdgeUV %d %d stat = %d\n",
               i+1, pindex[j], ptype[j], stat);
        continue;
      }
 */
      stat = EG_evaluate(faces[i], uv, result);
      if (stat != EGADS_SUCCESS) {
        printf(" EFace %3d: EG_evaluate UV stat = %d\n", i+1, stat);
        if (stat != EGADS_EXTRAPOL) continue;
      }
      box[0] = sqrt((xyzs[3*j  ]-result[0])*(xyzs[3*j  ]-result[0]) +
                    (xyzs[3*j+1]-result[1])*(xyzs[3*j+1]-result[1]) +
                    (xyzs[3*j+2]-result[2])*(xyzs[3*j+2]-result[2]));
/*    if (box[0] > 1.e-5)
        printf(" EFace %3d: EEdge = %d  pnt = %d  dist = %le\n", i+1,
               pindex[j], ptype[j], box[0]);  */
      if (box[0] > tol) tol = box[0];
    }
  }
  printf(" EFace/Node  deviation = %le\n", maxdev);
  printf(" EFace/EEdge deviation = %le\n", tol);
  EG_free(edges);
  EG_free(faces);
                      
#ifdef DISJOINTQUADS
  stat = EG_quadTess(tess, &quads);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeQuads = %d\n", stat);
    EG_deleteObject(ebody);
    EG_deleteObject(body);
    EG_close(context);
    exit(1);
  }
  EG_deleteObject(tess);
  tess = quads;
#endif
#endif
  printf(" \n");

  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
    printf(" failed to create wvContext!\n");
    EG_deleteObject(tess);
#ifndef NOEFFECT
    EG_deleteObject(ebody);
#endif
    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
    printf(" EG_close          = %d\n", EG_close(context));
    return 1;
  }

  /* make the scene */
  for (ngp = sum = stat = ibody = 0; ibody < nbody; ibody++) {

    quad = 0;
    stat = EG_attributeRet(tess, ".tessType", &atype,
                           &alen, &ints, &reals, &string);
    if (stat == EGADS_SUCCESS)
      if (atype == ATTRSTRING)
        if (strcmp(string, "Quad") == 0) quad = 1;

    /* get faces */
    for (i = 0; i < nface; i++) {
      stat = EG_getTessFace(tess, i+1, &len, &xyzs, &uvs, &ptype, &pindex,
                            &ntri, &tris, &tric);
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
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) nseg++;
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) nseg++;
      }
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) {
              segs[2*nseg  ] = tris[3*j+sides[k][0]];
              segs[2*nseg+1] = tris[3*j+sides[k][1]];
              nseg++;
            }
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) {
              segs[2*nseg  ] = tris[6*j+sideq[k][0]];
              segs[2*nseg+1] = tris[6*j+sideq[k][1]];
              nseg++;
            }
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
      if (stat > 0) ngp = stat+1;
      sum += ntri;
    }

    /* get edges */
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < nedge; i++) {
      stat  = EG_getTessEdge(tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
      nseg = len-1;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Edge %d\n",
               nseg, ibody, i+1);
        continue;
      }
      for (j = 0; j < len-1; j++) {
        segs[2*j  ] = j + 1;
        segs[2*j+1] = j + 2;
      }

      snprintf(gpname, 34, "Body %d Edge %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs, WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS,   &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_INDICES,  &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(segs);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
      if (stat < 0) {
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      } else {
        if (cntxt != NULL)
          if (cntxt->gPrims != NULL) {
            cntxt->gPrims[stat].lWidth = 1.5;
            if (wv_addArrowHeads(cntxt, stat, 0.05, 1, &nseg) != 0)
              printf(" wv_addArrowHeads Error\n");
          }
        ngp = stat+1;
      }
    }
    
    /* get nodes */
    color[0] = color[1] = color[2] = 0.0;
    realx    = (double *) EG_alloc(3*nnode*sizeof(double));
    if (realx == NULL) {
      printf(" Node Malloc Error!\n");
    } else {
      for (i = 0; i < nnode; i++)
        EG_getTopology(nodes[i], &geom, &oclass, &mtype, &realx[3*i],
                       &ndum, &dum, &senses);
      snprintf(gpname, 34, "Body %d Loop %d", ibody+1, 0);
      stat = wv_setData(WV_REAL64, nnode, (void *) realx, WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS,   &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      free(realx);
      stat = wv_addGPrim(cntxt, gpname, WV_POINT, WV_ON, 2, items);
      if (stat < 0) {
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      } else {
        if (cntxt != NULL)
          if (cntxt->gPrims != NULL)
            cntxt->gPrims[stat].pSize = 8.0;
        ngp = stat+1;
      }
    }
  }
  printf(" ** %d gPrims with %d triangles **\n", ngp, sum);

  /* start the server code */

  stat = 0;
  wv_setCallBack(cntxt, browserMessage);
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    /* we have a single valid server -- stay alive a long as we have a client */
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
    }
  }
  wv_cleanupServers();

  /* finish up */
  EG_free(nodes);
  
#ifdef SAVEMODEL
  save[0] = body;
  save[1] = ebody;
  save[2] = tess;
  stat = EG_makeTopology(context, NULL, MODEL, 3, NULL, 1, save, NULL, &model);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTopology on Model = %d\n", stat);
    EG_deleteObject(ebody);
    EG_deleteObject(ebody);
    EG_close(context);
    return 1;
  }
  stat = EG_saveModel(model, "effect.egads");
  if (stat != EGADS_SUCCESS)
    printf(" EG_saveModel = %d\n", stat);
  EG_deleteObject(model);

#else

  EG_deleteObject(tess);
#ifndef NOEFFECT
  EG_deleteObject(ebody);
#endif
  EG_deleteObject(body);

#endif

  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
