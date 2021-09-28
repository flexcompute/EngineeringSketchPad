/*
 ************************************************************************
 *                                                                      *
 * Tessellate.h - header for Tessellate.c                               *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2018/2021  John F. Dannenhoffer, III (Syracuse University)
 *
 * This library is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 *    License along with this library; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *     MA  02110-1301  USA
 */

#ifndef _TESSELLATE_H_
#define _TESSELLATE_H_

/*
 ************************************************************************
 *                                                                      *
 * Structures                                                           *
 *                                                                      *
 ************************************************************************
 */
#define TESS_MAGIC 4431000

typedef struct oct_T {
    int           npnt;                 /* number of Points */
    double        *xyz;                 /* array  of Points */
    int           ntri;                 /* number of Triangles */
    int           *trip;                /* array  of Triangles */
    double        xcent;                /* x centroid */
    double        ycent;                /* y centroid */
    double        zcent;                /* z centroid */
    struct oct_T  *child;               /* array of 8 children (or NULL) */
} oct_T;

typedef struct {
    int           magic;                /* magic number for TESS */
    int           ntri;                 /* number of Triangles */
    int           mtri;                 /* maximum   Triangles (malloc size) */
    int           *trip;                /* Points associated with each Triangle */
                                        /*    trip[3*i  ]  first  Point for Triangle i (bias-0) */
                                        /*    trip[3*i+1]  second Point for Triangle i (bias-0) */
                                        /*    trip[3*i+2]  third  Point for Triangle i (bias-0) */
    int           *trit;                /* neighboring Triangles assoc with each Triangle */
                                        /*    trit[3*i  ]  neighbor opposite Point 0 */
                                        /*    trit[3*i+1]  neighbor opposite Point 1 */
                                        /*    trit[3*i+2]  neighbor opposite Point 2 */
    int           *ttyp;                /* flag associated with each Triangle (see constants below) */
    int           nhang;                /* total number of hanging sides */
    int           nlink;                /* total number of sides with links */
    int           ncolr;                /* maximum color */
    double        *bbox;                /* bounding boxes of Triangles */
                                        /*    bbod[6*i  ]  minimum x */
                                        /*    bbod[6*i+1]  maximum x */
                                        /*    bbod[6*i+2]  minimum y */
                                        /*    bbod[6*i+3]  maximum y */
                                        /*    bbod[6*i+4]  minimum z */
                                        /*    bbod[6*i+5]  maximum z */
    int           npnt;                 /* number of Points */
    int           mpnt;                 /* maximum   Points (malloc size) */
    double        *xyz;                 /* coordinates of Points */
                                        /*    xyz[3*i  ] x-coordinate of Point i */
                                        /*    xyz[3*i+1] y-coordinate of Point i */
                                        /*    xyz[3*i+2] z-coordinate of Point i */
    double        *uv;                  /* parametric coordinates of Points (or NULL) */
                                        /*    uv[2*i  ] u-coordinate of Point i */
                                        /*    uv[2*i+1] v-coordinate of Point i */
    int           *ptyp;                /* flag associated with each Point (see constants below) */
    oct_T         *octree;              /* pointer to root of octree (or NULL) */
} tess_T;

typedef struct {
    int           pnt;                  /* Point number (bias-0) */
    int           tri;                  /* Triangle just after Point */
} seg_T;

/*
 ************************************************************************
 *                                                                      *
 * Callable routines                                                    *
 *                                                                      *
 ************************************************************************
 */
/* add a Point */
int addPoint(tess_T  *tess,             /* (in)  pointer to TESS */
             double  x,                 /* (in)  x-coordinate */
             double  y,                 /* (in)  y-coordinate */
             double  z);                /* (inO  z-coordinate */

/* add a Triangle */
int addTriangle(tess_T  *tess,          /* (in)  pointer to TESS */
                int     ip0,            /* (in)  index of first  Point (bias-0) */
                int     ip1,            /* (in)  index of second Point (bias-0) */
                int     ip2,            /* (in)  index of third  Point (bias-0) */
                int     it0,            /* (in)  index of first  Triangle (or -1) */
                int     it1,            /* (in)  index of first  Triangle (or -1) */
                int     it2);           /* (in)  index of first  Triangle (or -1) */

/* create a Triangle that bridges gap between Point and Triangle */
int bridgeToPoint(tess_T  *tess,        /* (in)  pointer to TESS */
                  int     itri,         /* (in)  index of Triangle (bias-0) */
                  int     ipnt);        /* (in)  index of Point    (bias-0) */

/* create two Triangles that bridge gap between given Triangles */
int bridgeTriangles(tess_T  *tess,      /* (in)  pointer to TESS */
                    int     itri,       /* (in)  index of first  Triangle (bias-0) */
                    int     jtri);      /* (in)  index of secind Triangle (bias-0) */

/* check areas in UV */
int checkAreas(tess_T  *tess,           /* (in)  pointer to TESS */
               int     *nneg,           /* (out) number of negative areas */
               int     *npos);          /* (out) number of positive areas */

/* color a Triangle and its neighbors (up to links) */
int colorTriangles(tess_T  *tess,       /* (in)  pointer to TESS */
                   int     itri,        /* (in)  index of Triangle (bias-0) */
                   int     icolr);      /* (in)  color (0-255) */

/* copy a Tessellation */
int copyTess(tess_T  *src,              /* (in)  pointer to source TESS */
             tess_T  *tgt);             /* (in)  pointer to target TESS */

/* create a Link on one side of a Triangle */
int createLink(tess_T  *tess,           /* (in)  pointer to TESS */
               int     itri,            /* (in)  index of Triangle (bias-0) */
               int     isid);           /* (in)  side index (0-2) */

/* create Links between given Points */
int createLinks(tess_T  *tess,          /* (in)  pointer to TESS */
                int     isrc,           /* (in)  index of source Point (bias-0) */
                int     itgt);          /* (in)  index of target Point (bias-0) */

/* cut Triangles in a Tessellation through given Points */
int cutTriangles(tess_T  *tess,         /* (in)  pointer to TESS */
                 int     icolr,         /* (in)  color (or -1 for all) */
                 double  data[]);       /* (in)  cut is data[0]+x*data[1]+y*data[2]+z*data[3] */

/* delete a Triangle */
int deleteTriangle(tess_T  *tess,       /* (in)  pointer to TESS */
                   int     itri);       /* (in)  index of Triangle to delete (bias-0) */

/* detect and link creases */
int detectCreases(tess_T  *tess,        /* (in)  pointer to TESS */
                  double  angdeg);      /* (in)  maximum crease angle (deg) */

/* extend loop to given x/y/z */
int extendLoop(tess_T  *tess,           /* (in)  pointer to TESS */
               int     ipnt,            /* (in)  index of one Point on loop (bias-0) */
               int     itype,           /* (in)  type: 1=x, 2=y, 3=z */
               double  val);            /* (in)  constant x/y/z for extension */

/* create a new Tessellation from Triangles of a given color */
int extractColor(tess_T  *tess,         /* (in)  pointer to TESS */
                 int     icolr,         /* (in)  color to extract */
                 tess_T  *subTess);     /* (out) pointer to sub-TESS object */

/* fill a loop with Triangles */
int fillLoop(tess_T  *tess,             /* (in)  pointer to TESS */
             int     ipnt);             /* (in)  index of one Point on loop (bias-0) */

/* find the loops */
int findLoops(tess_T  *tess,            /* (in)  pointer to TESS */
              int     *nloop,           /* (int) maximum loops returned */
                                        /* (out) number of loops found */
              int     ibeg[],           /* (out) array of a point in each loop */
              double  alen[]);          /* (out) length (in XYZ)  of each loop */

/* flatten coordinates of a given color */
int flattenColor(tess_T  *tess,         /* (in)  pointer to TESS */
                 int     icolr,         /* (in)  color index */
                 double  tol);          /* (in)  tolerance */

/* apply floater to get UV at interior Points */
int floaterUV(tess_T  *tess);           /* (in)  pointer to TESS */

/* free a Tessellation */
int freeTess(tess_T  *tess);            /* (in)  pointer to TESS to be freed */

/* get the Points along a loop */
int getLoop(tess_T  *tess,              /* (in)  pointer to TESS */
            int     ipnt,               /* (in)  index of one Point on loop (bias-0) */
            int     *nseg,              /* (out) number of segments in loop */
            seg_T   *seg[]);            /* (out) array  of Segments (freeable) */

/* initialize a Tessellation */
int initialTess(tess_T  *tess);         /* (in)  pointer to TESS */

/* initialize UV by projecting to best-fit plane of boundary Points */
int initialUV(tess_T  *tess);           /* (in)  pointer to TESS */

/* join two points (at their average) */
int joinPoints(tess_T  *tess,           /* (in)  pointer to TESS */
               int     ipnt,            /* (in)  index of first  Point (bias-0) */
               int     jpnt);           /* (in)  index of second Point (bias-0) */

/* make links between colors */
int makeLinks(tess_T  *tess);           /* (in)  pointer to TESS */

/* find nearest approach betwen two line segments */
int nearApproach(double xyz_a[],        /* (in)  coordinates of first  point */
                 double xyz_b[],        /* (in)  coordinates of second point */
                 double xyz_c[],        /* (in)  coordinates of third  point */
                 double xyz_d[],        /* (in)  coordinates of fourth point */
                 double *s_ab,          /* (out) fraction of distance from a to b */
                 double *s_cd,          /* (out) fraction of distance from c to d */
                 double *dist);         /* (out) shortest distance */

/* find nearest point to Tessellation */
int nearestTo(tess_T  *tess,            /* (in)  pointer to TESS */
              double  dbest,            /* (in)  initial best distance */
              double  xyz_in[],         /* (in)  input point */
              int     *ibest,           /* (out) Triangle index (bias-0) */
              double  xyz_out[]);       /* (out) output point */

/* read an ASCII stl file */
int readStlAscii(tess_T  *tess,         /* (in)  pointer to TESS */
                 char    *filename);    /* (in)  name of file */

/* read a binary stl file */
int readStlBinary(tess_T  *tess,        /* (in)  pointer to TESS */
                  char    *filename);   /* (in)  name of file */

/* read an ASCII tri file */
int readTriAscii(tess_T  *tess,         /* (in)  pointer to TESS */
                 char    *filename);    /* (in)  name of file */

/* remove all links */
int removeLinks(tess_T  *tess);         /* (in)  pointer to TESS */

/* scribe between given Points */
int scribe(tess_T  *tess,               /* (in)  pointer to TESS */
           int     isrc,                /* (in)  index of source Point (bias-0) */
           int     itgt);               /* (in)  index of target Point (bias-0) */

/* set up Triangle neighbor information (generally used internally) */
int setupNeighbors(tess_T  *tess);      /* (in)  pointer to TESS */

/* sort Triangles by color */
int sortTriangles(tess_T  *tess);       /* (in)  pointer to TESS */

/* split Triangle and its neighbor */
int splitTriangle(tess_T  *tess,        /* (in) pointer to TESS */
                  int     itri,         /* (in) Triangle index (bias-0) */
                  int     ipnt,         /* (in) startig Point index (bias-0) */
                  double  frac);        /* (in) fractional distance on side to split */

/* swap Triangles */
int swapTriangles(tess_T  *tess,        /* (in)  pointer to TESS */
                  int     itri,         /* (in)  index of first  Triangle (bias-0) */
                  int     jtri);        /* (in)  index of second Triangle (bias-0) */

/* evaluate at a given parametric coordinate */
int UVtoXYZ(tess_T  *tess,              /* (in)  pointer ot TESS */
            int     icolr,              /* (in)  color of Triangles */
            double  uv_in[],            /* (in)  parametric coordinates */
            double  xyz_out[]);         /* (out) physical coordinates */

/* write an ASCII triangle file */
int writeTriAscii(tess_T  *tess,        /* (in)  pointer to TESS */
                  char    *filename);   /* (in)  name of file */

/* write an ASCII stl file */
int writeStlAscii(tess_T  *tess,        /* (in)  pointer to TESS */
                  char    *filename);   /* (in)  name of file */

/* write an ASCII stl file */
int writeStlBinary(tess_T  *tess,       /* (in)  pointer to TESS */
                   char    *filename);  /* (in)  name of file */

/* find point nearest to the Tessellation */
int XYZtoUVXYZ(tess_T  *tess,           /* (in)  pointer to TESS */
               int     icolr,           /* (in)  color of Triangles */
               double  xyz_in[],        /* (in)  physical coordinates */
               double  uv_out[],        /* (out) parametric coordinates (or UNDEF) */
               double  xyz_out[]);      /* (out) coordinates of output coordinates */


/*
 ************************************************************************
 *                                                                      *
 * Defined constants                                                    *
 *                                                                      *
 ************************************************************************
 */
#define TRI_COLOR        0x0000ffff
#define TRI_VISIBLE      0x00010000
#define TRI_ACTIVE       0x00020000
#define TRI_LINK         0x001C0000
#define TRI_T0_LINK      0x00040000
#define TRI_T1_LINK      0x00080000
#define TRI_T2_LINK      0x00100000
#define TRI_EDGE         0x00E00000
#define TRI_T0_EDGE      0x00200000
#define TRI_T1_EDGE      0x00400000
#define TRI_T2_EDGE      0x00800000

#define PNT_INDEX        0x00ffffff
#define PNT_NODE         0x01000000
#define PNT_EDGE         0x02000000
#define PNT_FACE         0x04000000

#define UNDEF            -12345.6789

/*
 ************************************************************************
 *                                                                      *
 * Return codes (errors are -601 to -699)                               *
 *                                                                      *
 ************************************************************************
 */
#define SUCCESS                                 0

#define TESS_NOT_A_TESS                      -601
#define TESS_BAD_POINT_INDEX                 -602
#define TESS_BAD_TRIANGLE_INDEX              -603
#define TESS_BAD_VALUE                       -604
#define TESS_BAD_FILE_NAME                   -605
#define TESS_NOT_AN_ASCII_FILE               -606
#define TESS_NO_PARAMETRIC_COORDINATES       -607
#define TESS_NOT_IMPLEMENTED                 -608
#define TESS_NOT_CONVERGED                   -609
#define TESS_INTERNAL_ERROR                  -699

#endif  /* _TESSELLATE_H_ */
