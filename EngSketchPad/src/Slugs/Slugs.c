/*
 ************************************************************************
 *                                                                      *
 * Slugs.c -- server for Static Legacy Unstructured Geometry System     *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2022  John F. Dannenhoffer, III (Syracuse University)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#ifdef WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include <winsock2.h>
#endif

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#ifdef GRAFIC
    #include "grafic.h"
#endif

#ifdef WIN32
    #define snprintf   _snprintf
#endif

#define STRNCPY(A, B, LEN) strncpy(A, B, LEN); A[LEN-1] = '\0';

#include "egads.h"
#include "common.h"
#include "Fitter.h"
#include "Tessellate.h"
#include "RedBlackTree.h"

#include "wsserver.h"
#include "emp.h"

/*
 ***********************************************************************
 *                                                                     *
 * macros (including those that go along with common.h)                *
 *                                                                     *
 ***********************************************************************
 */
#ifdef DEBUG
   #define DOPEN {if (dbg_fp == NULL) dbg_fp = fopen("Slugs.dbg", "w");}
   static  FILE *dbg_fp=NULL;
#endif

static void *realloc_temp=NULL;            /* used by RALLOC macro */

#define  RED(COLOR)      (float)(COLOR / 0x10000        ) / (float)(255)
#define  GREEN(COLOR)    (float)(COLOR / 0x00100 % 0x100) / (float)(255)
#define  BLUE(COLOR)     (float)(COLOR           % 0x100) / (float)(255)

#define MAX_EXPR_LEN     128
#define MAX_STR_LEN    32767

#define MIN3(A,B,C)   (MIN(MIN(A,B),C))
#define MAX3(A,B,C)   (MAX(MAX(A,B),C))

#define UINT32 unsigned int
#define UINT16 unsigned short int
#define REAL32 float

/*
 ***********************************************************************
 *                                                                     *
 * structures                                                          *
 *                                                                     *
 ***********************************************************************
 */
#ifdef FOO
static float color_map[256*3] =
{ 0.0000,    0.0000,    1.0000,     0.0000,    0.0157,    1.0000,
  0.0000,    0.0314,    1.0000,     0.0000,    0.0471,    1.0000,
  0.0000,    0.0627,    1.0000,     0.0000,    0.0784,    1.0000,
  0.0000,    0.0941,    1.0000,     0.0000,    0.1098,    1.0000,
  0.0000,    0.1255,    1.0000,     0.0000,    0.1412,    1.0000,
  0.0000,    0.1569,    1.0000,     0.0000,    0.1725,    1.0000,
  0.0000,    0.1882,    1.0000,     0.0000,    0.2039,    1.0000,
  0.0000,    0.2196,    1.0000,     0.0000,    0.2353,    1.0000,
  0.0000,    0.2510,    1.0000,     0.0000,    0.2667,    1.0000,
  0.0000,    0.2824,    1.0000,     0.0000,    0.2980,    1.0000,
  0.0000,    0.3137,    1.0000,     0.0000,    0.3294,    1.0000,
  0.0000,    0.3451,    1.0000,     0.0000,    0.3608,    1.0000,
  0.0000,    0.3765,    1.0000,     0.0000,    0.3922,    1.0000,
  0.0000,    0.4078,    1.0000,     0.0000,    0.4235,    1.0000,
  0.0000,    0.4392,    1.0000,     0.0000,    0.4549,    1.0000,
  0.0000,    0.4706,    1.0000,     0.0000,    0.4863,    1.0000,
  0.0000,    0.5020,    1.0000,     0.0000,    0.5176,    1.0000,
  0.0000,    0.5333,    1.0000,     0.0000,    0.5490,    1.0000,
  0.0000,    0.5647,    1.0000,     0.0000,    0.5804,    1.0000,
  0.0000,    0.5961,    1.0000,     0.0000,    0.6118,    1.0000,
  0.0000,    0.6275,    1.0000,     0.0000,    0.6431,    1.0000,
  0.0000,    0.6588,    1.0000,     0.0000,    0.6745,    1.0000,
  0.0000,    0.6902,    1.0000,     0.0000,    0.7059,    1.0000,
  0.0000,    0.7216,    1.0000,     0.0000,    0.7373,    1.0000,
  0.0000,    0.7529,    1.0000,     0.0000,    0.7686,    1.0000,
  0.0000,    0.7843,    1.0000,     0.0000,    0.8000,    1.0000,
  0.0000,    0.8157,    1.0000,     0.0000,    0.8314,    1.0000,
  0.0000,    0.8471,    1.0000,     0.0000,    0.8627,    1.0000,
  0.0000,    0.8784,    1.0000,     0.0000,    0.8941,    1.0000,
  0.0000,    0.9098,    1.0000,     0.0000,    0.9255,    1.0000,
  0.0000,    0.9412,    1.0000,     0.0000,    0.9569,    1.0000,
  0.0000,    0.9725,    1.0000,     0.0000,    0.9882,    1.0000,
  0.0000,    1.0000,    0.9961,     0.0000,    1.0000,    0.9804,
  0.0000,    1.0000,    0.9647,     0.0000,    1.0000,    0.9490,
  0.0000,    1.0000,    0.9333,     0.0000,    1.0000,    0.9176,
  0.0000,    1.0000,    0.9020,     0.0000,    1.0000,    0.8863,
  0.0000,    1.0000,    0.8706,     0.0000,    1.0000,    0.8549,
  0.0000,    1.0000,    0.8392,     0.0000,    1.0000,    0.8235,
  0.0000,    1.0000,    0.8078,     0.0000,    1.0000,    0.7922,
  0.0000,    1.0000,    0.7765,     0.0000,    1.0000,    0.7608,
  0.0000,    1.0000,    0.7451,     0.0000,    1.0000,    0.7294,
  0.0000,    1.0000,    0.7137,     0.0000,    1.0000,    0.6980,
  0.0000,    1.0000,    0.6824,     0.0000,    1.0000,    0.6667,
  0.0000,    1.0000,    0.6510,     0.0000,    1.0000,    0.6353,
  0.0000,    1.0000,    0.6196,     0.0000,    1.0000,    0.6039,
  0.0000,    1.0000,    0.5882,     0.0000,    1.0000,    0.5725,
  0.0000,    1.0000,    0.5569,     0.0000,    1.0000,    0.5412,
  0.0000,    1.0000,    0.5255,     0.0000,    1.0000,    0.5098,
  0.0000,    1.0000,    0.4941,     0.0000,    1.0000,    0.4784,
  0.0000,    1.0000,    0.4627,     0.0000,    1.0000,    0.4471,
  0.0000,    1.0000,    0.4314,     0.0000,    1.0000,    0.4157,
  0.0000,    1.0000,    0.4000,     0.0000,    1.0000,    0.3843,
  0.0000,    1.0000,    0.3686,     0.0000,    1.0000,    0.3529,
  0.0000,    1.0000,    0.3373,     0.0000,    1.0000,    0.3216,
  0.0000,    1.0000,    0.3059,     0.0000,    1.0000,    0.2902,
  0.0000,    1.0000,    0.2745,     0.0000,    1.0000,    0.2588,
  0.0000,    1.0000,    0.2431,     0.0000,    1.0000,    0.2275,
  0.0000,    1.0000,    0.2118,     0.0000,    1.0000,    0.1961,
  0.0000,    1.0000,    0.1804,     0.0000,    1.0000,    0.1647,
  0.0000,    1.0000,    0.1490,     0.0000,    1.0000,    0.1333,
  0.0000,    1.0000,    0.1176,     0.0000,    1.0000,    0.1020,
  0.0000,    1.0000,    0.0863,     0.0000,    1.0000,    0.0706,
  0.0000,    1.0000,    0.0549,     0.0000,    1.0000,    0.0392,
  0.0000,    1.0000,    0.0235,     0.0000,    1.0000,    0.0078,
  0.0078,    1.0000,    0.0000,     0.0235,    1.0000,    0.0000,
  0.0392,    1.0000,    0.0000,     0.0549,    1.0000,    0.0000,
  0.0706,    1.0000,    0.0000,     0.0863,    1.0000,    0.0000,
  0.1020,    1.0000,    0.0000,     0.1176,    1.0000,    0.0000,
  0.1333,    1.0000,    0.0000,     0.1490,    1.0000,    0.0000,
  0.1647,    1.0000,    0.0000,     0.1804,    1.0000,    0.0000,
  0.1961,    1.0000,    0.0000,     0.2118,    1.0000,    0.0000,
  0.2275,    1.0000,    0.0000,     0.2431,    1.0000,    0.0000,
  0.2588,    1.0000,    0.0000,     0.2745,    1.0000,    0.0000,
  0.2902,    1.0000,    0.0000,     0.3059,    1.0000,    0.0000,
  0.3216,    1.0000,    0.0000,     0.3373,    1.0000,    0.0000,
  0.3529,    1.0000,    0.0000,     0.3686,    1.0000,    0.0000,
  0.3843,    1.0000,    0.0000,     0.4000,    1.0000,    0.0000,
  0.4157,    1.0000,    0.0000,     0.4314,    1.0000,    0.0000,
  0.4471,    1.0000,    0.0000,     0.4627,    1.0000,    0.0000,
  0.4784,    1.0000,    0.0000,     0.4941,    1.0000,    0.0000,
  0.5098,    1.0000,    0.0000,     0.5255,    1.0000,    0.0000,
  0.5412,    1.0000,    0.0000,     0.5569,    1.0000,    0.0000,
  0.5725,    1.0000,    0.0000,     0.5882,    1.0000,    0.0000,
  0.6039,    1.0000,    0.0000,     0.6196,    1.0000,    0.0000,
  0.6353,    1.0000,    0.0000,     0.6510,    1.0000,    0.0000,
  0.6667,    1.0000,    0.0000,     0.6824,    1.0000,    0.0000,
  0.6980,    1.0000,    0.0000,     0.7137,    1.0000,    0.0000,
  0.7294,    1.0000,    0.0000,     0.7451,    1.0000,    0.0000,
  0.7608,    1.0000,    0.0000,     0.7765,    1.0000,    0.0000,
  0.7922,    1.0000,    0.0000,     0.8078,    1.0000,    0.0000,
  0.8235,    1.0000,    0.0000,     0.8392,    1.0000,    0.0000,
  0.8549,    1.0000,    0.0000,     0.8706,    1.0000,    0.0000,
  0.8863,    1.0000,    0.0000,     0.9020,    1.0000,    0.0000,
  0.9176,    1.0000,    0.0000,     0.9333,    1.0000,    0.0000,
  0.9490,    1.0000,    0.0000,     0.9647,    1.0000,    0.0000,
  0.9804,    1.0000,    0.0000,     0.9961,    1.0000,    0.0000,
  1.0000,    0.9882,    0.0000,     1.0000,    0.9725,    0.0000,
  1.0000,    0.9569,    0.0000,     1.0000,    0.9412,    0.0000,
  1.0000,    0.9255,    0.0000,     1.0000,    0.9098,    0.0000,
  1.0000,    0.8941,    0.0000,     1.0000,    0.8784,    0.0000,
  1.0000,    0.8627,    0.0000,     1.0000,    0.8471,    0.0000,
  1.0000,    0.8314,    0.0000,     1.0000,    0.8157,    0.0000,
  1.0000,    0.8000,    0.0000,     1.0000,    0.7843,    0.0000,
  1.0000,    0.7686,    0.0000,     1.0000,    0.7529,    0.0000,
  1.0000,    0.7373,    0.0000,     1.0000,    0.7216,    0.0000,
  1.0000,    0.7059,    0.0000,     1.0000,    0.6902,    0.0000,
  1.0000,    0.6745,    0.0000,     1.0000,    0.6588,    0.0000,
  1.0000,    0.6431,    0.0000,     1.0000,    0.6275,    0.0000,
  1.0000,    0.6118,    0.0000,     1.0000,    0.5961,    0.0000,
  1.0000,    0.5804,    0.0000,     1.0000,    0.5647,    0.0000,
  1.0000,    0.5490,    0.0000,     1.0000,    0.5333,    0.0000,
  1.0000,    0.5176,    0.0000,     1.0000,    0.5020,    0.0000,
  1.0000,    0.4863,    0.0000,     1.0000,    0.4706,    0.0000,
  1.0000,    0.4549,    0.0000,     1.0000,    0.4392,    0.0000,
  1.0000,    0.4235,    0.0000,     1.0000,    0.4078,    0.0000,
  1.0000,    0.3922,    0.0000,     1.0000,    0.3765,    0.0000,
  1.0000,    0.3608,    0.0000,     1.0000,    0.3451,    0.0000,
  1.0000,    0.3294,    0.0000,     1.0000,    0.3137,    0.0000,
  1.0000,    0.2980,    0.0000,     1.0000,    0.2824,    0.0000,
  1.0000,    0.2667,    0.0000,     1.0000,    0.2510,    0.0000,
  1.0000,    0.2353,    0.0000,     1.0000,    0.2196,    0.0000,
  1.0000,    0.2039,    0.0000,     1.0000,    0.1882,    0.0000,
  1.0000,    0.1725,    0.0000,     1.0000,    0.1569,    0.0000,
  1.0000,    0.1412,    0.0000,     1.0000,    0.1255,    0.0000,
  1.0000,    0.1098,    0.0000,     1.0000,    0.0941,    0.0000,
  1.0000,    0.0784,    0.0000,     1.0000,    0.0627,    0.0000,
  1.0000,    0.0471,    0.0000,     1.0000,    0.0314,    0.0000,
  1.0000,    0.0157,    0.0000,     1.0000,    0.0000,    0.0000 };
#endif

static float lims[2]  = {-2.0, +2.0};

/* data structures */
typedef struct {
    int                  ipnt;               /* Point index (bias-0) */
    int                  nedg;               /* number of incident Edges */
    double               x;                  /* x-coordinate */
    double               y;                  /* y-coordinate */
    double               z;                  /* z-coordinate */
    ego                  enode;              /* ego object */
} node_T;

typedef struct {
    int                  ibeg;               /* Node at beginning (bias-1) */
    int                  iend;               /* Node at end       (bias-1)*/
    int                  ileft;              /* Face on left      (bias-1) */
    int                  irite;              /* Face on rite      (bias-1)*/
    int                  npnt;               /* number of Points          along Edge */
    int                  *pnt;               /* array  of Point indices   along Edge */
    double               *xyz;               /* array  of Point locations along Edge */
    int                  ncp;                /* number of control points */
    double               *cp;                /* array  of control points */
    ego                  eedge;              /* ego object */
} edge_T;

typedef struct {
    int                  icol;               /* color  of associated Triangles */
    tess_T               tess;               /* TESS object */
    int                  nedg;               /* number of associated Edges */
    int                  *edg;               /* array  of associated Edges (signed, bias-1) */
    int                  nlup;               /* number of Loops */
    int                  *lup;               /* index (in *edg) of first Edge in Loop */
    int                  npnt;               /* number of interior Points in Face */
    double               *xyz;               /* array  of interior Points in Face */
    int                  ntrain;             /* number of training Points in Face */
    double               *xyztrain;          /* array  of training Points in Face */
    int                  ncp;                /* number of control points */
    double               *cp;                /* array  of control points */
    int                  done;               /* =1 when done */
    ego                  eface;              /* ego object */
} face_T;

typedef struct {
    int                  ibeg;               /* Point index at beg     (bias-0) */
    int                  iend;               /* Point index at end     (bias-0) */
    int                  prev;               /* previous Segment or -1 (bias-0) */
    int                  next;               /* next     Segment or -1 (bias-0) */
} sgmt_T;

typedef struct {
    void      *mutex;                        /* the mutex or NULL for single thread */
    long      master;                        /* master thread ID */
} emp_T;

/*
 ***********************************************************************
 *                                                                     *
 * global variables                                                    *
 *                                                                     *
 ***********************************************************************
 */
/* global variables associated with tessellation, ... */
static   tess_T          tess;               /* global TESS object */
static   tess_T          tess_undo;          /* undo copy of TESS object */

static   int             Mnode = 0;          /* maximum   Nodes */
static   int             Nnode = 0;          /* number of Nodes */
static   node_T          *node = NULL;       /* array  of Nodes */
static   int             Medge = 0;          /* maximum   Edges */
static   int             Nedge = 0;          /* number of Edges */
static   edge_T          *edge = NULL;       /* array  of Edges */
static   int             Mface = 0;          /* maximum   Faces */
static   int             Nface = 0;          /* number of Faces */
static   face_T          *face = NULL;       /* array  of Faces */

static int               outLevel = 1;       /* default output level */
static char              casename[257];      /* name of case */

/* global variables associated with graphical user interface (gui) */
static wvContext         *cntxt;             /* context for the WebViewer */
static int               port   = 7681;      /* port number */
static int               batch = 0;          /* =0 to enable visualization */
static float             focus[4];

/* global variables associated with scene graph meta-data */
#define MAX_METADATA_LEN 32000
static char              *sgFocusData = NULL;

/* global variables associated with scene graph */
static int               Tris_pend    =   0;  /* =1 if awaiting update */

/* global variables associated with CurPt */
static int               CurPt_index  = -1;  /* index (or -1 for none) */
static int               CurPt_gprim  = -1;  /* GPrim (or -1 for none) */
static int               CurPt_pend   =  0;  /* =1 if awaiting update */

/* global variables associated with Hangs */
static int               Hangs_gprim  = -1;  /* GPrim (or -1 for none) */
static int               Hangs_pend   = 0;   /* =1 if awaiting update */

/* global variables associated with Links */
static int               Links_gprim  = -1;  /* GPrim (or -1 for none) */
static int               Links_pend   = 0;   /* =1 if awaiting update */

/* global variables associated with the response buffer */
static int               max_resp_len = 0;
static int               response_len = 0;
static char              *response    = NULL;

/* global variables associated with journals */
static FILE              *jrnl_out = NULL;   /* output journal file */

/* global variables associated with subsampling */
static int               subsample = 1;
static int               nctrlpnt  = 0;

/*
 ***********************************************************************
 *                                                                     *
 * declarations for routines defined below                             *
 *                                                                     *
 ***********************************************************************
 */

static int        buildCurPt();
static int        buildHangs();
static int        buildLinks();
static int        buildTriangles();

static int        generateBrep(char message[]);
static int        generateEgads(char egadsname[], char message[]);
static int        generateFits(int ncp, char message[]);
static void       empFit2dCloud(void *struc);
static int        makeNodesAndEdges(int nsgmt, sgmt_T sgmt[],
                                    int ibeg, int iend, int nodnum[], int icolr, int jcolr);

static void       processMessage(char *text);
static int        getToken(char *text, int nskip, char *token);
static int        closestPoint(double xloc, double yloc, double zloc);
static int        closestTriangle(double xloc, double yloc, double zloc);
static int        storeUndo();
static void       setColor(int rgb, float color[]);
#ifdef FOO
static void       spec_col(float scalar, float lims[], float out[]);
#endif
       void       printEgo(ego obj);

//$$$static int        fitPlane(double xyz[], int n, double *a, double *b, double *c, double *d);

extern int        EG_getPlane(ego eloop, ego *eplane);

#ifdef GRAFIC
    static void   plotPoints_image(int*, void*, void*, void*, void*, void*,
                                   void*, void*, void*, void*, void*, float*, char*, int);
    extern int    plotCurve(int m, double XYZcloud[], /*@null@*/double Tcloud[],
                            int n, double cp[], double normf, double dotmin, int nmin);
    extern int    plotSurface(int m, double XYZcloud[], /*@null@*/double Ucloud[],
                              int n, double cp[], double normf, int nmin);
#endif


/*
 ***********************************************************************
 *                                                                     *
 *   Main program                                                      *
 *                                                                     *
 ***********************************************************************
 */
int
main(int       argc,                /* (in)  number of arguments */
     char      *argv[])             /* (in)  array of arguments */
{

    int       status, status2, i, bias, showUsage=0, ihand;
    int       iedge, iface;
    float     fov, zNear, zFar;
    float     eye[3]    = {0.0, 0.0, 7.0};
    float     center[3] = {0.0, 0.0, 0.0};
    float     up[3]     = {0.0, 1.0, 0.0};
    char      *filename=NULL, *jrnlname=NULL, *tempname=NULL;
    char      *text=NULL, *wv_start, test[80];
    void      *temp;
    FILE      *jrnl_in=NULL, *fp=NULL;

#ifdef GRAFIC
    int      io_kbd=5, io_scr=6;
    char     pltitl[255];
#endif

    ROUTINE(MAIN);

    /* --------------------------------------------------------------- */

    /* initialize the random number seed */
    srand(12345);

    /* dynamically allocated array so that everything is on heap (not stack) */
    MALLOC(sgFocusData, char, 256        );
    MALLOC(filename,    char, 257        );
    MALLOC(jrnlname,    char, 257        );
    MALLOC(tempname,    char, 257        );
    MALLOC(text,        char, MAX_STR_LEN);

    /* get the flags and casename(s) from the command line */
    casename[0] = '\0';
    jrnlname[0] = '\0';

    for (i = 1; i < argc; i++) {
        if        (strcmp(argv[i], "-port") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &port);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-jrnl") == 0) {
            if (i < argc-1) {
                STRNCPY(jrnlname, argv[++i], 257);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-outLevel") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &outLevel);
                if (outLevel < 0) outLevel = 0;
                if (outLevel > 3) outLevel = 3;
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-batch") == 0) {
            batch = 1;
       } else if (strcmp(argv[i], "-nctrlpnt") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &nctrlpnt);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-subsample") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &subsample);
                if (subsample < 1) subsample = 1;
            } else {
                showUsage = 1;
                break;
            }
        } else if (strlen(casename) == 0) {
            STRNCPY(casename, argv[i], 257);
        } else {
            SPRINT0(0, "two casenames given");
            showUsage = 1;
            break;
        }
    }

    if (showUsage) {
        SPRINT0(0, "proper usage: 'Slugs [-port X] [-jrnl jrnlname] [-outLevel X] [-batch] [casename[.stl]]'");
        SPRINT0(0, "STOPPING...\a");
        exit(0);
    }

    /* check size of various types used within .stl files */
    if        (sizeof(UINT16) != 2) {
        SPRINT0(0, "ERROR:: uint16 should have size 2");
        exit(0);
    } else if (sizeof(UINT32) != 4) {
        SPRINT0(0, "ERROR:: uint32 should have size 4");
        exit(0);
    } else if (sizeof(REAL32) != 4) {
        SPRINT0(0, "ERROR:: real32 should have size 4");
        exit(0);
    }

    /* welcome banner */
    SPRINT0(1, "**********************************************************");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*                    Program Slugs                       *");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*        written by John Dannenhoffer, 2013/2022         *");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "**********************************************************");

#ifdef GRAFIC
    /* initialize the grafics */
    sprintf(pltitl, "Program Fitter2D");
    grinit_(&io_kbd, &io_scr, pltitl, strlen(pltitl));
#endif

    /* allocate nominal response buffer */
    max_resp_len = 4096;
    MALLOC(response, char, max_resp_len);

    /* add .stl to filename if not present */
    if (strlen(casename) > 0) {
        STRNCPY(filename, casename, 257);
        if (strstr(casename, ".stl") == NULL && strstr(casename, ".tri") == NULL) {
            strcat(filename, ".stl");
        }
    } else {
        SPRINT0(0, "ERROR:: a casename must be given");
        exit(0);
    }

    /* read the .stl or .tri file */
    if (strstr(filename, ".tri") != NULL) {
        status = readTriAscii(&tess, filename);
        SPRINT1(3, "--> readTriAscii 0> status=%d", status);

    } else if (strstr(filename, ".stl") != NULL) {
        fp = fopen(filename, "rb");
        if (fp == NULL) {
            SPRINT1(0, "ERROR:: \"%s\" does not exist", filename);
            exit(0);
        }

        (void) fread(test, sizeof(char), 5, fp);
        fclose(fp);

        if (strncmp(test, "solid", 5) == 0) {
            SPRINT1(1, "--> \"%s\" is an ASCII file", filename);
            status = readStlAscii(&tess, filename);
            SPRINT1(3, "--> readStlAscii -> status=%d", status);
        } else {
            SPRINT1(1, "--> \"%s\" is a binary file", filename);
            status = readStlBinary(&tess, filename);
            SPRINT1(3, "--> readStlBinary -> status=%d", status);
        }
    } else {
        SPRINT1(0, "ERROR:: \"%s\" is not a .stl or .tri file", filename);
        exit(0);
    }

    /* make links between the colors */
    if (tess.ncolr > 0) {
        status = makeLinks(&tess);
        SPRINT1(3, "--> makeLinks -> status=%d", status);
    }

    /* report initial statistics */
    SPRINT0(2, "==> initialization complete");
    SPRINT1(2, "    npnt  = %d", tess.npnt );
    SPRINT1(2, "    ntri  = %d", tess.ntri );
    SPRINT1(2, "    ncolr = %d", tess.ncolr);
    SPRINT1(2, "    nhang = %d", tess.nhang);
    SPRINT1(2, "    nlink = %d", tess.nlink);

    /* initialize the undo copy */
    status = initialTess(&tess_undo);
    SPRINT1(3, "initialTess -> status=%d", status);

    /* open the output journal file */
    snprintf(tempname, 257, "port%d.jrnl", port);
    jrnl_out = fopen(tempname, "w");

    fprintf(jrnl_out, "# casename=%s\n\n", casename);
    fflush( jrnl_out);

    /* initialize the scene graph meta data */
    if (batch == 0) {
        SPLINT_CHECK_FOR_NULL(sgFocusData);

        sgFocusData[0] = '\0';
    }

    /* create the WebViewer context */
    if (batch == 0) {
        bias  =  0;
        fov   = 30.0;
        zNear =  1.0;
        zFar  = 10.0;
        cntxt = wv_createContext(bias, fov, zNear, zFar, eye, center, up);
        if (cntxt == NULL) {
            SPRINT0(0, "ERROR:: failed to create wvContext");
            exit(0);
        }
    }

    /* build the initial Scene Graph */
    if (batch == 0) {
        status = buildTriangles();
        SPRINT1(3, "--> buildTriangles -> status=%d", status);

        status = buildHangs();
        SPRINT1(3, "--> buildHangs -> status=%d", status);

        status = buildLinks();
        SPRINT1(3, "--> buildLinks -> status=%d", status);
    }

    /* process the input journal file if jrnlname exists */
    if (strlen(jrnlname) > 0) {
        SPRINT1(0, "==> Opening input journal file \"%s\"", jrnlname);

        jrnl_in = fopen(jrnlname, "r");
        if (jrnl_in == NULL) {
            SPRINT0(0, "ERROR:: Journal file cannot be opened");
            exit(0);
        } else {
            while (1) {
                temp = fgets(text, MAX_STR_LEN, jrnl_in);
                if (temp == NULL) break;

                if (feof(jrnl_in)) break;

                /* if there is a \n, convert it to a \0 */
                text[MAX_STR_LEN-1] = '\0';

                for (i = 0; i < strlen(text); i++) {
                    if (text[i] == '\n') {
                        text[i] =  '\0';
                        break;
                    }
                }

                if (strncmp(text, "##end##", 7) == 0) break;

                processMessage(text);
            }

            fclose(jrnl_in);

            SPRINT0(0, "==> Closing input journal file");
        }
    }

    /* get the command to start the client (if any) */
    if (batch == 0) {
        wv_start = getenv("SLUGS_START");
    }

    /* start the server */
    if (batch == 0) {
        status2 = SUCCESS;
        if (wv_startServer(port, NULL, NULL, NULL, 0, cntxt) == 0) {

            /* stay alive a long as we have a client */
            while (wv_statusServer(0)) {
                usleep(100000);

                /* start the browser if the first time through this loop */
                if (status2 == SUCCESS) {
                    if (wv_start != NULL) {
                        status2 += system(wv_start);
                    }

                    status2++;
                }

                /* start hand-shaking */
                if (Tris_pend  != 0 ||
                    Hangs_pend != 0 ||
                    Links_pend != 0 ||
                    CurPt_pend != 0   ) {
                    ihand = 1;
                    if (wv_handShake(cntxt) != 1) {
                        SPRINT0(0, "ERROR:: handShake out of Sync 1");
                    }
                } else {
                    ihand = 0;
                }

                /* update Triangles if there are changes pending */
                if (Tris_pend != 0) {
                    Tris_pend = 0;

                    status = buildTriangles();
                    SPRINT1(3, "buildTriangles -> status=%d", status);
                }

                /* update Hangs if there are changes pending */
                if (Hangs_pend != 0) {
                    Hangs_pend = 0;

                    status = buildHangs();
                    SPRINT1(3, "buildHangs -> status=%d", status);
                }

                /* update Links if there are changes pending */
                if (Links_pend != 0) {
                    Links_pend = 0;

                    status = buildLinks();
                    SPRINT1(3, "buildLinks -> status=%d", status);
                }

                /* update CurPt if there are changes pending */
                if (CurPt_pend != 0) {
                    CurPt_pend = 0;

                    status = buildCurPt();
                    SPRINT1(3, "buildCurPt -> status=%d", status);
                }

                /* complete hand-shaking since all changes have been made */
                if (ihand == 1) {
                    if (wv_handShake(cntxt) != 0) {
                        SPRINT0(0, "ERROR:: handShake out of Sync 0");
                    }
                }
            }
        }
    }

    /* cleanup and exit */
    if (jrnl_out != NULL) fclose(jrnl_out);
    jrnl_out = NULL;

    for (iface = 1; iface <= Nface; iface++) {
        status = freeTess(&(face[iface].tess));
        SPRINT2(2, "freeTess(iface=%d) -> status=%d", iface, status);

        FREE(face[iface].edg);
        FREE(face[iface].lup);
        FREE(face[iface].xyz);
        FREE(face[iface].xyztrain);
        FREE(face[iface].cp );
    }

    for (iedge = 1; iedge <= Nedge; iedge++) {
        FREE(edge[iedge].pnt);
        FREE(edge[iedge].xyz);
        FREE(edge[iedge].cp );
    }

    FREE(face);
    FREE(edge);
    FREE(node);

    status = freeTess(&tess);
    SPRINT1(2, "freeTess(tess) -> status=%d", status);

    status = freeTess(&tess_undo);
    SPRINT1(2, "freeTess(tess_undo) -> status=%d", status);

    wv_cleanupServers();
    status = SUCCESS;

    SPRINT0(1, "==> Slugs completed successfully");

cleanup:
    FREE(response   );
    FREE(text       );
    FREE(tempname   );
    FREE(jrnlname   );
    FREE(filename   );
    FREE(sgFocusData);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   buildCurPt - make/update CurPt in scene graph                     *
 *                                                                     *
 ***********************************************************************
 */
static int
buildCurPt()
{
    int       status = SUCCESS;

    int       attrs;
    float     *pnt=NULL;

    float     col[3];
    char      gpname[33];

    wvData    items[7];

    ROUTINE(buildCurPt);

    /* --------------------------------------------------------------- */

    if (CurPt_index < 0) {
        goto cleanup;
    }

    /* note: we need to update more than one Point because of
       an undocumented "feature" in wv */
    MALLOC(pnt, float, 6);

    pnt[0] = tess.xyz[3*CurPt_index  ];   pnt[3] = pnt[0];
    pnt[1] = tess.xyz[3*CurPt_index+1];   pnt[4] = pnt[1];
    pnt[2] = tess.xyz[3*CurPt_index+2];   pnt[5] = pnt[2];
    status = wv_setData(WV_REAL32, 2, (void*)pnt, WV_VERTICES, &(items[0]));
    SPRINT4(3, "wv_setData(VERTICES, %f, %f, %f) -> status=%d", pnt[0], pnt[1], pnt[2], status);

    wv_adjustVerts(&(items[0]), focus);

    if (CurPt_gprim < 0) {
        snprintf(gpname, 32, "CurPt");
        attrs = WV_ON;

        setColor(0x000000, col);
        status = wv_setData(WV_REAL32, 1, (void*)col, WV_COLORS, &(items[1]));
        SPRINT4(3, "wv_setData(COLORS, %f, %f, %f) -> status=%d", col[0], col[1], col[2], status);

        CurPt_gprim = wv_addGPrim(cntxt, gpname, WV_POINT, attrs, 2, items);
        SPRINT1(3, "wv_addGPrim(WV_POINT) -> CurPt_gprim=%d", CurPt_gprim);

        cntxt->gPrims[CurPt_gprim].pSize = 8.0;
    } else {
        SPRINT1(3, "CurPt_gprim=%d", CurPt_gprim);
        status = wv_modGPrim(cntxt, CurPt_gprim, 1, items);
        SPRINT2(3, "wv_modGPrim(CurPt_gprim=%d) -> status=%d", CurPt_gprim, status);
    }

cleanup:
    FREE(pnt);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   buildHangs - make/update Hangs in scene graph                     *
 *                                                                     *
 ***********************************************************************
 */
static int
buildHangs()
{
    int       status = SUCCESS;

    int       itri, jtri, nhang, ip0, ip1, ip2;
    int       attrs;
    float     *hang=NULL;

    float     color[3];
    char      gpname[33];

    wvData    items[7];

    ROUTINE(buildHangs);

    /* --------------------------------------------------------------- */

    /* determine the number of Hangs in the current tessellation */
    tess.nhang = 0;
    for (itri = 0; itri < tess.ntri; itri++) {
        if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;

        jtri = tess.trit[3*itri  ];
        if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
            tess.nhang++;
        }

        jtri = tess.trit[3*itri+1];
        if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
            tess.nhang++;
        }

        jtri = tess.trit[3*itri+2];
        if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
            tess.nhang++;
        }
    }

    /* plot Hangs (if there are any) */
    if (tess.nhang > 0) {
        snprintf(gpname, 32, "Hangs");
        attrs = WV_ON;

        MALLOC(hang, float, 6*tess.nhang);

        nhang = 0;
        for (itri = 0; itri < tess.ntri; itri++) {
            if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;

            jtri = tess.trit[3*itri  ];
            if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
                ip1 = tess.trip[3*itri+1];
                ip2 = tess.trip[3*itri+2];

                hang[6*nhang  ] = tess.xyz[3*ip1  ];
                hang[6*nhang+1] = tess.xyz[3*ip1+1];
                hang[6*nhang+2] = tess.xyz[3*ip1+2];
                hang[6*nhang+3] = tess.xyz[3*ip2  ];
                hang[6*nhang+4] = tess.xyz[3*ip2+1];
                hang[6*nhang+5] = tess.xyz[3*ip2+2];
                nhang++;
            }

            jtri = tess.trit[3*itri+1];
            if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
                ip2 = tess.trip[3*itri+2];
                ip0 = tess.trip[3*itri  ];

                hang[6*nhang  ] = tess.xyz[3*ip2  ];
                hang[6*nhang+1] = tess.xyz[3*ip2+1];
                hang[6*nhang+2] = tess.xyz[3*ip2+2];
                hang[6*nhang+3] = tess.xyz[3*ip0  ];
                hang[6*nhang+4] = tess.xyz[3*ip0+1];
                hang[6*nhang+5] = tess.xyz[3*ip0+2];
                nhang++;
            }

            jtri = tess.trit[3*itri+2];
            if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
                ip0 = tess.trip[3*itri  ];
                ip1 = tess.trip[3*itri+1];

                hang[6*nhang  ] = tess.xyz[3*ip0  ];
                hang[6*nhang+1] = tess.xyz[3*ip0+1];
                hang[6*nhang+2] = tess.xyz[3*ip0+2];
                hang[6*nhang+3] = tess.xyz[3*ip1  ];
                hang[6*nhang+4] = tess.xyz[3*ip1+1];
                hang[6*nhang+5] = tess.xyz[3*ip1+2];
                nhang++;
            }

            if (nhang > tess.nhang) {
                SPRINT2(0, "ERROR:: nhang=%d but tess.nhang=%d", nhang, tess.nhang);
                status = -999;
                goto cleanup;
            }
        }

        /* set up vertices for the Hangs */
        status = wv_setData(WV_REAL32, 2*nhang, (void*)hang, WV_VERTICES, &(items[0]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        wv_adjustVerts(&(items[0]), focus);

        /* hang color */
        setColor(0xff0000, color);
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[1]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        /* if a graphic primitive does not exist yet, make it */
        if (Hangs_gprim < 0) {
            Hangs_gprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, 2, items);
            SPRINT1(3, "wv_addGPrim(WV_LINE) -> Hangs_gprim=%d", Hangs_gprim);
            if (Hangs_gprim >= 0) {

                /* make line width 5 (does not work for ANGLE) */
                cntxt->gPrims[Hangs_gprim].lWidth = 5.0;
            }

        /* otherwise, simply modify the verticies in the previous one */
        } else {
            status = wv_modGPrim(cntxt, Hangs_gprim, 1, items);
            SPRINT2(3, "wv_modGPrim(Hangs_gprim=%d) -> status=%d", Hangs_gprim, status);
        }

    /* if no Hangs, remove the previous gprim */
    } else if (Hangs_gprim >= 0) {
        wv_removeGPrim(cntxt, Hangs_gprim);
        SPRINT0(3, "wv_removeGPrim");

        Hangs_gprim = -1;
    }

cleanup:
    FREE(hang);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   buildLinks - make/update Links in scene graph                     *
 *                                                                     *
 ***********************************************************************
 */
static int
buildLinks()
{
    int       status = SUCCESS;

    int       itri, nlink, ip0, ip1, ip2;
    int       attrs;
    float     *link=NULL;

    float     color[3];
    char      gpname[33];

    wvData    items[7];

    ROUTINE(buildLinks);

    /* --------------------------------------------------------------- */

    /* plot Links (if there are any) */
    if (tess.nlink > 0) {
        snprintf(gpname, 32, "Links");
        attrs = WV_ON;

        MALLOC(link, float, 6*tess.nlink);

        nlink = 0;
        for (itri = 0; itri < tess.ntri; itri++) {
            if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;

            if (tess.ttyp[itri] & TRI_T0_LINK) {
                if (itri > tess.trit[3*itri  ]) {
                    ip1 = tess.trip[3*itri+1];
                    ip2 = tess.trip[3*itri+2];

                    link[6*nlink  ] = tess.xyz[3*ip1  ];
                    link[6*nlink+1] = tess.xyz[3*ip1+1];
                    link[6*nlink+2] = tess.xyz[3*ip1+2];
                    link[6*nlink+3] = tess.xyz[3*ip2  ];
                    link[6*nlink+4] = tess.xyz[3*ip2+1];
                    link[6*nlink+5] = tess.xyz[3*ip2+2];
                    nlink++;
                }
            }
            if (tess.ttyp[itri] & TRI_T1_LINK) {
                if (itri > tess.trit[3*itri+1]) {
                    ip2 = tess.trip[3*itri+2];
                    ip0 = tess.trip[3*itri  ];

                    link[6*nlink  ] = tess.xyz[3*ip2  ];
                    link[6*nlink+1] = tess.xyz[3*ip2+1];
                    link[6*nlink+2] = tess.xyz[3*ip2+2];
                    link[6*nlink+3] = tess.xyz[3*ip0  ];
                    link[6*nlink+4] = tess.xyz[3*ip0+1];
                    link[6*nlink+5] = tess.xyz[3*ip0+2];
                    nlink++;
                }
            }
            if (tess.ttyp[itri] & TRI_T2_LINK) {
                if (itri > tess.trit[3*itri+2]) {
                    ip0 = tess.trip[3*itri  ];
                    ip1 = tess.trip[3*itri+1];

                    link[6*nlink  ] = tess.xyz[3*ip0  ];
                    link[6*nlink+1] = tess.xyz[3*ip0+1];
                    link[6*nlink+2] = tess.xyz[3*ip0+2];
                    link[6*nlink+3] = tess.xyz[3*ip1  ];
                    link[6*nlink+4] = tess.xyz[3*ip1+1];
                    link[6*nlink+5] = tess.xyz[3*ip1+2];
                    nlink++;
                }
            }

            if (nlink > tess.nlink) {
                SPRINT2(0, "ERROR:: nlink=%d but tess.nlink=%d", nlink, tess.nlink);
                status = -999;
                goto cleanup;
            }
        }

        /* set up vertices for the Links */
        status = wv_setData(WV_REAL32, 2*nlink, (void*)link, WV_VERTICES, &(items[0]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        wv_adjustVerts(&(items[0]), focus);

        /* link color */
        setColor(0xffffff, color);
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[1]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        /* if a graphic primitive does not exist yet, make it */
        if (Links_gprim < 0) {
            Links_gprim = wv_addGPrim(cntxt, gpname, WV_LINE, attrs, 2, items);
            SPRINT1(3, "wv_addGPrim(WV_LINE) -> Links_gprim=%d", Links_gprim);
            if (Links_gprim >= 0) {

                /* make line width 5 (does not work for ANGLE) */
                cntxt->gPrims[Links_gprim].lWidth = 5.0;
            }

        /* otherwise, simply modify the verticies in the previous one */
        } else {
            status = wv_modGPrim(cntxt, Links_gprim, 1, items);
            SPRINT2(3, "wv_modGPrim(Links_gprim=%d) -> status=%d", Links_gprim, status);
        }

    /* if no Links, remove the previous gprim */
    } else if (Links_gprim >= 0) {
        wv_removeGPrim(cntxt, Links_gprim);
        SPRINT0(3, "wv_removeGPrim");

        Links_gprim = -1;
    }

cleanup:
    FREE(link);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   buildTriangles - make/update Triangles in scene graph             *
 *                                                                     *
 ***********************************************************************
 */
static int
buildTriangles()
{
    int       status = SUCCESS;

    int       icolr, itri, jtri, ntri, ipnt, npnt, nseg;
    int       igprim, attrs;
    int       *tri=NULL, *seg=NULL;
    float     *xyz=NULL;

    float     color[3];
    double    bigbox[6], size;
    char      gpname[33];

    wvData    items[7];

    ROUTINE(buildTriangles);

    /* --------------------------------------------------------------- */

    /* remove any graphic primitives that already exist */
    wv_removeAll(cntxt);

    /* reset all the GPrim values since they do not exist anymore */
    CurPt_gprim = -1;
    Hangs_gprim = -1;
    Links_gprim = -1;

    /* find the values needed to adjust the vertices */
    bigbox[0] = bigbox[1] = bigbox[2] = +HUGEQ;
    bigbox[3] = bigbox[4] = bigbox[5] = -HUGEQ;

    for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
        if (tess.xyz[3*ipnt  ] < bigbox[0]) bigbox[0] = tess.xyz[3*ipnt  ];
        if (tess.xyz[3*ipnt+1] < bigbox[1]) bigbox[1] = tess.xyz[3*ipnt+1];
        if (tess.xyz[3*ipnt+2] < bigbox[2]) bigbox[2] = tess.xyz[3*ipnt+2];
        if (tess.xyz[3*ipnt  ] > bigbox[3]) bigbox[3] = tess.xyz[3*ipnt  ];
        if (tess.xyz[3*ipnt+1] > bigbox[4]) bigbox[4] = tess.xyz[3*ipnt+1];
        if (tess.xyz[3*ipnt+2] > bigbox[5]) bigbox[5] = tess.xyz[3*ipnt+2];
    }

                                    size = bigbox[3] - bigbox[0];
    if (size < bigbox[4]-bigbox[1]) size = bigbox[4] - bigbox[1];
    if (size < bigbox[5]-bigbox[2]) size = bigbox[5] - bigbox[2];

    focus[0] = (bigbox[0] + bigbox[3]) / 2;
    focus[1] = (bigbox[1] + bigbox[4]) / 2;
    focus[2] = (bigbox[2] + bigbox[5]) / 2;
    focus[3] = size;

    /* generate the scene graph focus data */
    snprintf(sgFocusData, 256, "sgFocus;[%20.12e,%20.12e,%20.12e,%20.12e]", focus[0], focus[1], focus[2], focus[3]);

    /* loop through the Colors */
    for (icolr = 0; icolr <= tess.ncolr; icolr++) {

        /* name and attributes */
        snprintf(gpname, 32, "Color %d", icolr);
        attrs = WV_ON | WV_ORIENTATION | WV_LINES;

        /* find number of Triangles with this color */
        ntri = 0;
        for (itri = 0; itri < tess.ntri; itri++) {
            if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;

            if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
                ntri++;
            }
        }

        /* if there are not Triangles with this color, skip the
           rest of the processing */
        if (ntri <= 0) continue;

        /* allocate storage for vertices and Triangles */
        MALLOC(xyz, float, 9*ntri);
        MALLOC(tri, int,   3*ntri);
        MALLOC(seg, int,   6*ntri);

        /* set up vertices and segments of this color */
        npnt = 0;
        ntri = 0;
        nseg = 0;
        for (itri = 0; itri < tess.ntri; itri++) {
            if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;

            if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {

                ipnt = tess.trip[3*itri  ];
                xyz[3*npnt  ] = tess.xyz[3*ipnt  ];
                xyz[3*npnt+1] = tess.xyz[3*ipnt+1];
                xyz[3*npnt+2] = tess.xyz[3*ipnt+2];
                npnt++;

                ipnt = tess.trip[3*itri+1];
                xyz[3*npnt  ] = tess.xyz[3*ipnt  ];
                xyz[3*npnt+1] = tess.xyz[3*ipnt+1];
                xyz[3*npnt+2] = tess.xyz[3*ipnt+2];
                npnt++;

                ipnt = tess.trip[3*itri+2];
                xyz[3*npnt  ] = tess.xyz[3*ipnt  ];
                xyz[3*npnt+1] = tess.xyz[3*ipnt+1];
                xyz[3*npnt+2] = tess.xyz[3*ipnt+2];
                npnt++;

                tri[3*ntri  ] = npnt - 3;
                tri[3*ntri+1] = npnt - 2;
                tri[3*ntri+2] = npnt - 1;
                ntri++;

                jtri = tess.trit[3*itri  ];
                if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
                } else  if ((tess.ttyp[itri] & TRI_T0_LINK) == 0) {
                    if (itri > tess.trit[3*itri  ]) {
                        seg[2*nseg  ] = npnt - 2;
                        seg[2*nseg+1] = npnt - 1;
                        nseg++;
                    }
                }

                jtri = tess.trit[3*itri+1];
                if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
                } else  if ((tess.ttyp[itri] & TRI_T1_LINK) == 0) {
                    if (itri > tess.trit[3*itri+1]) {
                        seg[2*nseg  ] = npnt - 1;
                        seg[2*nseg+1] = npnt - 3;
                        nseg++;
                    }
                }

                jtri = tess.trit[3*itri+2];
                if (jtri < 0 || (tess.ttyp[jtri] & TRI_ACTIVE) == 0) {
                } else if ((tess.ttyp[itri] & TRI_T2_LINK) == 0) {
                    if (itri > tess.trit[3*itri+2]) {
                        seg[2*nseg  ] = npnt - 3;
                        seg[2*nseg+1] = npnt - 2;
                        nseg++;
                    }
                }
            }
        }

        /* (indexed) verticies */
        status = wv_setData(WV_REAL32, npnt, (void*)xyz, WV_VERTICES, &(items[0]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        wv_adjustVerts(&(items[0]), focus);

        /* (indexed) Triangles */
        status = wv_setData(WV_INT32, 3*ntri, (void*)tri, WV_INDICES, &(items[1]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        /* constant Triangle colors */
        if         (icolr == 0) {
            setColor(0xcfcfcf, color);         /* light grey */
        } else if ((icolr-1)%12 ==  0) {
            setColor(0xffcfcf, color);         /* light red */
        } else if ((icolr-1)%12 ==  1) {
            setColor(0xcfffcf, color);         /* light green */
        } else if ((icolr-1)%12 ==  2) {
            setColor(0xcfcfff, color);         /* light blue */
        } else if ((icolr-1)%12 ==  3) {
            setColor(0xcfffff, color);         /* light cyan */
        } else if ((icolr-1)%12 ==  4) {
            setColor(0xffcfff, color);         /* light magenta */
        } else if ((icolr-1)%12 ==  5) {
            setColor(0xffffcf, color);         /* light yellow */
        } else if ((icolr-1)%12 ==  6) {
            setColor(0xff7f7f, color);         /* medium red */
        } else if ((icolr-1)%12 ==  7) {
            setColor(0x7fff7f, color);         /* medium green */
        } else if ((icolr-1)%12 ==  8) {
            setColor(0x7f7fff, color);         /* medium blue */
        } else if ((icolr-1)%12 ==  9) {
            setColor(0x7fffff, color);         /* medium cyan */
        } else if ((icolr-1)%12 == 10) {
            setColor(0xff7fff, color);         /* medium magenta */
        } else {
            setColor(0xffff7f, color);         /* medium yellow */
        }

        status = wv_setData(WV_REAL32, 1, (void*)color, WV_COLORS, &(items[2]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        /* (indexed) Triangle sides */
        status = wv_setData(WV_INT32, 2*nseg, (void*)seg, WV_LINDICES, &(items[3]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        /* segment colors */
        setColor(0x000000, color);
        status = wv_setData(WV_REAL32, 1, (void*)color, WV_LCOLOR, &(items[4]));
        if (status != SUCCESS) {
            SPRINT1(3, "wv_setData -> status=%d", status);
        }

        /* make graphic primitive */
        igprim = wv_addGPrim(cntxt, gpname, WV_TRIANGLE, attrs, 5, items);
        SPRINT1(3, "wv_addGPrim(WV_TRIANGLE) -> igprim=%d", igprim);
        if (igprim >= 0) {

            /* make line width 1 */
            cntxt->gPrims[igprim].lWidth = 1.0;
        }

        /* free storage associated with this color */
        FREE(xyz);
        FREE(tri);
        FREE(seg);
    }

cleanup:
    FREE(xyz);
    FREE(tri);
    FREE(seg);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   generateBrep - generate Brep based upon Triangle colors           *
 *                                                                     *
 *   sets up:                                                          *
 *      node[].ipnt                                                    *
 *      node[].nedg                                                    *
 *      node[].x                                                       *
 *      node[].y                                                       *
 *      node[].z                                                       *
 *      edge[].ibeg                                                    *
 *      edge[].iend                                                    *
 *      edge[].ileft                                                   *
 *      edge[].irite                                                   *
 *      edge[].npnt                                                    *
 *      edge[].xyz                                                     *
 *      face[].icol                                                    *
 *      face[].tess                                                    *
 *      face[].nedg                                                    *
 *      face[].edg                                                     *
 *      face[].nlup                                                    *
 *      face[].lup                                                     *
 *                                                                     *
 ***********************************************************************
 */
static int
generateBrep(char   message[])          /* (out) error message */
{
    int       status = SUCCESS;

    int       nchange, uncolored;
    int       icolr, jcolr, itri, jtri, ntri, ip0, ip1, ip2;
    int       isgmt, nsgmt, ibeg, iend, ilup, jlup, klup;
    int       ipnt, jpnt, count, done, i, j, swap;
    int       inode, iedge, jedge, iface;
    int       *nodnum=NULL, *edgtmp=NULL, *luptmp=NULL;
    double    areax, areay, areaz, amax, *area=NULL;
    char      filename[257];

    sgmt_T     *sgmt=NULL;

    ROUTINE(generateBrep);

    /* --------------------------------------------------------------- */

    strcpy(message, "okay");

    SPRINT0(1, "\nGenerating Brep...");

    /* determine if there are any coincident points */
    for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
        for (jpnt = ipnt+1; jpnt < tess.npnt; jpnt++) {
            if (fabs(tess.xyz[3*ipnt  ]-tess.xyz[3*jpnt  ]) < EPS06 &&
                fabs(tess.xyz[3*ipnt+1]-tess.xyz[3*jpnt+1]) < EPS06 &&
                fabs(tess.xyz[3*ipnt+2]-tess.xyz[3*jpnt+2]) < EPS06   ) {
                printf("duplicate point found\n");
                printf("ipnt=%6d %20.10f %20.10f %20.10f\n", ipnt, tess.xyz[3*ipnt], tess.xyz[3*ipnt+1], tess.xyz[3*ipnt+2]);
                printf("jpnt=%6d %20.10f %20.10f %20.10f\n", jpnt, tess.xyz[3*jpnt], tess.xyz[3*jpnt+1], tess.xyz[3*jpnt+2]);
            }
        }
    }

    /* clear the Node, Edge, and Face tables */
    Mnode = 0;
    Nnode = 0;
    FREE(node);

    for (iedge = 1; iedge <= Nedge; iedge++) {
        FREE(edge[iedge].pnt);
        FREE(edge[iedge].xyz);
        FREE(edge[iedge].cp );
    }
    Medge = 0;
    Nedge = 0;
    FREE(edge);

    for (iface = 1; iface <= Nface; iface++) {
        FREE(face[iface].edg);
        FREE(face[iface].lup);
        FREE(face[iface].xyz);
        FREE(face[iface].xyztrain);
        FREE(face[iface].cp );
    }

    Mface = 0;
    Nface = 0;
    FREE(face);

    /* make sure all Triangles are colored */
    uncolored = 0;
    for (itri = 0; itri < tess.ntri; itri++) {
        if ((tess.ttyp[itri] & TRI_COLOR) == 0) {
            uncolored++;
        }
    }

    if (uncolored > 0 && uncolored < tess.ntri) {
        SPRINT1(-1,   "ERROR:: there are %d uncolored Triangles", uncolored);
        if (uncolored < 20) {
            for (itri = 0; itri < tess.ntri; itri++) {
                if ((tess.ttyp[itri] & TRI_COLOR) == 0) {
                    ip0 = tess.trip[3*itri  ];
                    ip1 = tess.trip[3*itri+1];
                    ip2 = tess.trip[3*itri+2];
                    SPRINT4(0, "        itri=%5d at (%12.5f %12.5f %12.5f)", itri,
                            (tess.xyz[3*ip0  ]+tess.xyz[3*ip1  ]+tess.xyz[3*ip2  ])/3,
                            (tess.xyz[3*ip0+1]+tess.xyz[3*ip1+1]+tess.xyz[3*ip2+1])/3,
                            (tess.xyz[3*ip0+2]+tess.xyz[3*ip1+2]+tess.xyz[3*ip2+2])/3);
                }
            }
        }
        snprintf(message, 80, "there are %d uncolored Triangles", uncolored);
        status = -999;
        goto cleanup;
    }

    /* allocate Segment table (larger than needed) */
    MALLOC(sgmt, sgmt_T, tess.ntri);

    /* allocate a table which identifies the Node at any Point */
    MALLOC(nodnum, int, tess.npnt);
    for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
        nodnum[ipnt] = -1;
    }

    /* special processing if all Triangles are uncolored */
    if (uncolored == tess.ntri) {
        printf("    special processing for all uncolored Triangles\n");

        Mnode = 1;
        Medge = 1;
        Mface = 1;

        icolr = 1;

    /* determine the number of Faces */
    } else {
        Mface = 0;
        for (icolr = 1; icolr <= tess.ncolr; icolr++) {
            for (itri = 0; itri < tess.ntri; itri++) {
                if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
                    Mface++;
                    break;
                }
            }
        }
        SPRINT1(1, "   there are %d Faces", Mface);

        /* preallocate Node, Edge, and Face tables */
        Mnode = 100;
        Medge = 100;
//      Mface set above
    }

    MALLOC(node, node_T, Mnode+1);
    MALLOC(edge, edge_T, Medge+1);
    MALLOC(face, face_T, Mface+1);

    /* initialize the Faces */
    SPRINT0(1, "Initializing Faces...");
    if (tess.ncolr == 0) {
        Nface++;

        face[Nface].icol  = icolr;
//      face[Nface].tess  = NULL;
        face[Nface].nedg  = 0;
        face[Nface].edg   = NULL;
        face[Nface].nlup  = 0;
        face[Nface].lup   = NULL;
        face[Nface].npnt  = 0;
        face[Nface].xyz   = NULL;
        face[Nface].ntrain   = 0;
        face[Nface].xyztrain = NULL;
        face[Nface].ncp   = 0;
        face[Nface].cp    = NULL;
//      face[Nface].eface = NULL;

        if (Nface == 1) {
            status = copyTess(&tess, &(face[Nface].tess));
            CHECK_STATUS(copyTess);
        } else {
            status = extractColor(&tess, icolr, &(face[Nface].tess));
            CHECK_STATUS(extractColor);
        }

        SPRINT4(1, "   created Face %3d .icol=%6d, .npnt=%6d, .ntri=%6d", Nface,
                face[Nface].icol,
                face[Nface].tess.npnt,
                face[Nface].tess.ntri);

    } else {
        for (icolr = 1; icolr <= tess.ncolr; icolr++) {
            ntri = 0;
            for (itri = 0; itri < tess.ntri; itri++) {
                if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
                    ntri++;
                }
            }

            if (ntri > 0) {
                Nface++;

                face[Nface].icol  = icolr;
//              face[Nface].tess  = NULL;
                face[Nface].nedg  = 0;
                face[Nface].edg   = NULL;
                face[Nface].nlup  = 0;
                face[Nface].lup   = NULL;
                face[Nface].npnt  = 0;
                face[Nface].xyz   = NULL;
                face[Nface].ntrain   = 0;
                face[Nface].xyztrain = NULL;
                face[Nface].ncp   = 0;
                face[Nface].cp    = NULL;
//             face[Nface].eface = NULL;

                status = extractColor(&tess, icolr, &(face[Nface].tess));
                CHECK_STATUS(extractColor);

                SPRINT4(1, "   created Face %3d .icol=%6d, .npnt=%6d, .ntri=%6d", Nface,
                        face[Nface].icol,
                        face[Nface].tess.npnt,
                        face[Nface].tess.ntri);
            }
        }
    }

    /* loop through all color pairs to find possible Edges */
    SPRINT0(1, "Looking for possible Edges (and Nodes)...");
    if (tess.ncolr == 0) {
        int    icorn1, icorn2, icorn3, icorn4, kpnt, found, *interior=NULL, ninterior;
        int    nS=0, nN=0, nW=0, nE=0;
        double *xyzS=NULL, *xyzN=NULL, *xyzW=NULL, *xyzE=NULL;

        FILE   *fp;

        printf("Enter point numbers at 4 corners of the Face:\n");
        scanf("%d %d %d %d", &icorn1, &icorn2, &icorn3, &icorn4);

        MALLOC(interior, int, tess.npnt);

        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            interior[ipnt] = 1;
        }

        /* add .dat to filename */
        STRNCPY(filename, casename, 257);
        strcat(filename, ".dat");

        fp = fopen(filename, "w");

        /* find the segments between corners */
        for (iedge = 0; iedge < 4; iedge++) {
            if (iedge == 0) {
                jpnt = icorn1;
                kpnt = icorn2;
            } else if (iedge == 1) {
                jpnt = icorn2;
                kpnt = icorn3;
            } else if (iedge == 2) {
                jpnt = icorn3;
                kpnt = icorn4;
            } else {
                jpnt = icorn4;
                kpnt = icorn1;
            }

            nsgmt = 0;
            while (jpnt != kpnt) {
                found = 0;
                for (itri = 0; itri < tess.ntri; itri++) {
                    if        (tess.trip[3*itri  ] == jpnt && tess.trit[3*itri+1] < 0) {
                        sgmt[nsgmt].ibeg = tess.trip[3*itri  ];
                        sgmt[nsgmt].iend = tess.trip[3*itri+2];
                        sgmt[nsgmt].prev = nsgmt - 1;
                        sgmt[nsgmt].next = nsgmt + 1;
                        nsgmt++;

                        jpnt = tess.trip[3*itri+2];
                        found++;
                    } else if (tess.trip[3*itri+1] == jpnt && tess.trit[3*itri+2] < 0) {
                        sgmt[nsgmt].ibeg = tess.trip[3*itri+1];
                        sgmt[nsgmt].iend = tess.trip[3*itri  ];
                        sgmt[nsgmt].prev = nsgmt - 1;
                        sgmt[nsgmt].next = nsgmt + 1;
                        nsgmt++;

                        jpnt = tess.trip[3*itri  ];
                        found++;
                    } else if (tess.trip[3*itri+2] == jpnt && tess.trit[3*itri  ] < 0) {
                        sgmt[nsgmt].ibeg = tess.trip[3*itri+2];
                        sgmt[nsgmt].iend = tess.trip[3*itri+1];
                        sgmt[nsgmt].prev = nsgmt - 1;
                        sgmt[nsgmt].next = nsgmt + 1;
                        nsgmt++;

                        jpnt = tess.trip[3*itri+1];
                        found++;
                    }
                    if (jpnt == kpnt) break;
                }
                if (found == 0) break;
            }

            sgmt[      0].prev = -1;
            sgmt[nsgmt-1].next = -1;

            status = makeNodesAndEdges(nsgmt, sgmt, 0, nsgmt-1, nodnum, 1, 0);
            CHECK_STATUS(makeNodesAndEdges);

            if (iedge == 0) {
                nS = nsgmt + 1;
                MALLOC(xyzS, double, 3*nS);

                ipnt = sgmt[0].ibeg;
                xyzS[0] = tess.xyz[3*ipnt  ];
                xyzS[1] = tess.xyz[3*ipnt+1];
                xyzS[2] = tess.xyz[3*ipnt+2];
                interior[ipnt] = 0;

                for (isgmt = 0; isgmt < nsgmt; isgmt++) {
                    ipnt = sgmt[isgmt].iend;
                    xyzS[3*isgmt+3] = tess.xyz[3*ipnt  ];
                    xyzS[3*isgmt+4] = tess.xyz[3*ipnt+1];
                    xyzS[3*isgmt+5] = tess.xyz[3*ipnt+2];
                    interior[ipnt] = 0;
                }
            } else if (iedge == 1) {
                nE = nsgmt + 1;
                MALLOC(xyzE, double, 3*nE);

                ipnt = sgmt[0].ibeg;
                xyzE[0] = tess.xyz[3*ipnt  ];
                xyzE[1] = tess.xyz[3*ipnt+1];
                xyzE[2] = tess.xyz[3*ipnt+2];
                interior[ipnt] = 0;

                for (isgmt = 0; isgmt < nsgmt; isgmt++) {
                    ipnt = sgmt[isgmt].iend;
                    xyzE[3*isgmt+3] = tess.xyz[3*ipnt  ];
                    xyzE[3*isgmt+4] = tess.xyz[3*ipnt+1];
                    xyzE[3*isgmt+5] = tess.xyz[3*ipnt+2];
                    interior[ipnt] = 0;
                }
            } else if (iedge == 2) {
                nN = nsgmt + 1;
                MALLOC(xyzN, double, 3*nN);

                ipnt = sgmt[0].ibeg;
                xyzN[0] = tess.xyz[3*ipnt  ];
                xyzN[1] = tess.xyz[3*ipnt+1];
                xyzN[2] = tess.xyz[3*ipnt+2];
                interior[ipnt] = 0;

                for (isgmt = 0; isgmt < nsgmt; isgmt++) {
                    ipnt = sgmt[isgmt].iend;
                    xyzN[3*isgmt+3] = tess.xyz[3*ipnt  ];
                    xyzN[3*isgmt+4] = tess.xyz[3*ipnt+1];
                    xyzN[3*isgmt+5] = tess.xyz[3*ipnt+2];
                    interior[ipnt] = 0;
                }
            } else {
                nW = nsgmt + 1;
                MALLOC(xyzW, double, 3*nW);

                ipnt = sgmt[0].ibeg;
                xyzW[0] = tess.xyz[3*ipnt  ];
                xyzW[1] = tess.xyz[3*ipnt+1];
                xyzW[2] = tess.xyz[3*ipnt+2];
                interior[ipnt] = 0;

                for (isgmt = 0; isgmt < nsgmt; isgmt++) {
                    ipnt = sgmt[isgmt].iend;
                    xyzW[3*isgmt+3] = tess.xyz[3*ipnt  ];
                    xyzW[3*isgmt+4] = tess.xyz[3*ipnt+1];
                    xyzW[3*isgmt+5] = tess.xyz[3*ipnt+2];
                    interior[ipnt] = 0;
                }
            }
        }

        /* add the south boundary to .dat file */
        fprintf(fp, "%5d%5d  south\n", nS, 0);

        for (i = 0; i < nS; i++) {
            fprintf(fp, "%15.7f %15.7f %15.7f\n", xyzS[3*i], xyzS[3*i+1], xyzS[3*i+2]);
        }

        FREE(xyzS);

        /* add  the north boundary to .dat file */
        fprintf(fp, "%5d%5d  north\n", nN, 0);

        for (i = 0; i < nN; i++) {
            j = nN - 1 - i;
            fprintf(fp, "%15.7f %15.7f %15.7f\n", xyzN[3*j], xyzN[3*j+1], xyzN[3*j+2]);
        }

        FREE(xyzN);

        /* add the west boundary to the .dat file */
        fprintf(fp, "%5d%5d  west\n", nW, 0);

        for (i = 0; i < nW; i++) {
            j = nW - 1 - i;
            fprintf(fp, "%15.7f %15.7f %15.7f\n", xyzW[3*j], xyzW[3*j+1], xyzW[3*j+2]);
        }

        FREE(xyzW);

        /* add the east boundary to the .dat file */
        fprintf(fp, "%5d%5d  east\n", nE, 0);

        for (i = 0; i < nE; i++) {
            fprintf(fp, "%15.7f %15.7f %15.7f\n", xyzE[3*i], xyzE[3*i+1], xyzE[3*i+2]);
        }

        FREE(xyzE);

        /* add interior points to .dat file */
        ninterior = 0;
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (interior[ipnt] == 1) ninterior++;
        }
        fprintf(fp, "%5d%5d  interior\n", ninterior, 0);
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (interior[ipnt] == 1) {
                fprintf(fp, "%15.7f %15.7f %15.7f\n",
                        tess.xyz[3*ipnt], tess.xyz[3*ipnt+1], tess.xyz[3*ipnt+2]);
            }
        }

        fclose(fp);

        FREE(interior);

    } else {
        for (icolr = 1; icolr <= tess.ncolr; icolr++) {
            for (jcolr = icolr+1; jcolr <= tess.ncolr; jcolr++) {

                /* find all posible Segments */
                nsgmt = 0;
                for (itri = 0; itri < tess.ntri; itri++) {
                    if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
                        jtri = tess.trit[3*itri  ];
                        if (jtri >= 0) {
                            if ((tess.ttyp[jtri] & TRI_COLOR) == jcolr) {
                                sgmt[nsgmt].ibeg = tess.trip[3*itri+1];
                                sgmt[nsgmt].iend = tess.trip[3*itri+2];
                                sgmt[nsgmt].prev = -1;
                                sgmt[nsgmt].next = -1;
                                nsgmt++;
                            }
                        }

                        jtri = tess.trit[3*itri+1];
                        if (jtri >= 0) {
                            if ((tess.ttyp[jtri] & TRI_COLOR) == jcolr) {
                                sgmt[nsgmt].ibeg = tess.trip[3*itri+2];
                                sgmt[nsgmt].iend = tess.trip[3*itri  ];
                                sgmt[nsgmt].prev = -1;
                                sgmt[nsgmt].next = -1;
                                nsgmt++;
                            }
                        }

                        jtri = tess.trit[3*itri+2];
                        if (jtri >= 0) {
                            if ((tess.ttyp[jtri] & TRI_COLOR) == jcolr) {
                                sgmt[nsgmt].ibeg = tess.trip[3*itri  ];
                                sgmt[nsgmt].iend = tess.trip[3*itri+1];
                                sgmt[nsgmt].prev = -1;
                                sgmt[nsgmt].next = -1;
                                nsgmt++;
                            }
                        }
                    }
                }
                SPRINT3(2, "icolr=%2d  jcolr=%2d  nsgmt=%d", icolr, jcolr, nsgmt);

                /* if there are no Segments, there is nothing to do */
                if (nsgmt <= 0) continue;

                /* arrange the Segments head to tail.  the while loop is needed
                   because there may be more than one set of disjoint segments
                   between two colors */
                while (1) {

                    /* find the first segment that is not used */
                    ibeg = -1;
                    for (isgmt = 0; isgmt < nsgmt; isgmt++) {
                        if (sgmt[isgmt].prev == -1 && sgmt[isgmt].next == -1) {
                            ibeg = isgmt;         /* first Segment */
                            iend = isgmt;         /* last  Segment */
                            break;
                        }
                    }

                    /* if all the segments are used, we are done */
                    if (ibeg < 0) break;

                    /* keep adding Segments to end while possible */
                    nchange = 1;
                    while (nchange > 0) {
                        nchange = 0;

                        if (sgmt[ibeg].ibeg == sgmt[iend].iend) {
                            break;
                        }

                        for (isgmt = 1; isgmt < nsgmt; isgmt++) {
                            if (sgmt[isgmt].prev >= 0) {
                                continue;
                            } else if (sgmt[isgmt].ibeg == sgmt[iend].iend) {
                                sgmt[isgmt].prev = iend;
                                sgmt[iend ].next = isgmt;

                                iend = isgmt;
                                nchange++;
                            }
                        }
                    }

                    /* keep adding Segments to beginning while possible */
                    nchange = 1;
                    while (nchange > 0) {
                        nchange = 0;

                        if (sgmt[ibeg].ibeg == sgmt[iend].iend) {
                            break;
                        }

                        for (isgmt = 1; isgmt < nsgmt; isgmt++) {
                            if (sgmt[isgmt].next >= 0) {
                                continue;
                            } else if (sgmt[isgmt].iend == sgmt[ibeg].ibeg) {
                                sgmt[isgmt].next = ibeg;
                                sgmt[ibeg ].prev = isgmt;

                                ibeg = isgmt;
                                nchange++;
                            }
                        }
                    }

                    /* if only one Segment for this Edge, specially mark the Segment so that
                       it does not get used again */
                    if (ibeg == iend) {
                        sgmt[ibeg].prev = -2;
                        sgmt[ibeg].next = -2;
                    }

                    status = makeNodesAndEdges(nsgmt, sgmt, ibeg, iend, nodnum, icolr, jcolr);
                    CHECK_STATUS(makeNodesAndEdges);
                }
            }
        }
    }

    /* inform the Faces of their incident Edges */
    for (iface = 1; iface <= Nface; iface++) {

        MALLOC(face[iface].edg, int, face[iface].nedg  );
        MALLOC(face[iface].lup, int, face[iface].nedg+1);   // too big for safety

        for (jedge = 0; jedge < face[iface].nedg; jedge++) {
            face[iface].edg[jedge] = 0;           // needed to avoid clang warning
        }

        jedge = 0;
        for (iedge = 1; iedge <= Nedge; iedge++) {
            if (edge[iedge].ileft == iface) {
                face[iface].edg[jedge++] = +iedge;
            }
            if (edge[iedge].irite == iface) {
                face[iface].edg[jedge++] = -iedge;
            }
        }
    }

    /* reorder the Edges in each Face to form loops */
    for (iface = 1; iface <= Nface; iface++) {
        printf("Face %6d\n", iface);

        if (face[iface].nedg <= 0) continue;

        printf("...at beginning\n");
        for (i = 0; i < face[iface].nedg; i++) {
            iedge = face[iface].edg[i];
            if        (iedge > 0) {
                printf("     Edge %6d, npnt=%6d, ibeg=%3d, iend=%3d\n",
                       iedge, edge[+iedge].npnt, edge[+iedge].ibeg, edge[+iedge].iend);
            } else if (iedge < 0) {
                printf("     Edge %6d, npnt=%6d, iend=%3d, ibeg=%3d\n",
                       iedge, edge[-iedge].npnt, edge[-iedge].iend, edge[-iedge].ibeg);
            } else {
                printf("     Edge %6d, degenerate,   ibeg=%3d, iend=%3d\n",
                       iedge,                    edge[+iedge].ibeg, edge[+iedge].iend);
            }
        }

        done  = 0;
        face[iface].nlup = 0;

        /* while(1) to make multiple Loops */
        while (1) {

            face[iface].lup[face[iface].nlup] = done;

            /* start the next Loop at the first available Edge */
            iedge = face[iface].edg[done];
            if        (iedge > 0) {
                ibeg = edge[+iedge].ibeg;
                iend = edge[+iedge].iend;
            } else {
                ibeg = edge[-iedge].iend;
                iend = edge[-iedge].ibeg;
            }
            done++;

            /* while(1) to add Edges to current Loop until it closes */
            count=0;

            while (1) {
                for (i = done; i < face[iface].nedg; i++) {
                    iedge = face[iface].edg[i];
                    if        (iedge > 0) {
                        if (edge[+iedge].ibeg == iend) {
                            if (i > done) {
                                swap                  = face[iface].edg[done];
                                face[iface].edg[done] = face[iface].edg[i   ];
                                face[iface].edg[i   ] = swap;
                            }
                            iend = edge[+iedge].iend;
                            done++;
                            break;
                        }
                    } else {
                        if (edge[-iedge].iend == iend) {
                            if (i > done) {
                                swap                  = face[iface].edg[done];
                                face[iface].edg[done] = face[iface].edg[i   ];
                                face[iface].edg[i   ] = swap;
                            }
                            iend = edge[-iedge].ibeg;
                            done++;
                            break;
                        }
                    }
                }

                /* check if Loop is closed */
                if (iend == ibeg) {
                    break;
                }

                /* infinite loop safety */
                count++;
                if (count > 100) {
                    SPRINT0(-1,   "ERROR:: could not link Edges node to tail");
                    snprintf(message, 80, "could not link Edges node to tail");
                    status = -999;
                    goto cleanup;
                }
            }

            /* finish this Loop */
            face[iface].nlup++;

            /* if all Edges are used, we are done */
            if (done == face[iface].nedg) break;
        }

        /* mark end of last Loop */
        face[iface].lup[face[iface].nlup] = face[iface].nedg;

        printf("...after sorting into Loops\n");
        for (ilup = 0; ilup < face[iface].nlup; ilup++) {
            printf("   Loop %6d (%d:%d)\n", ilup, face[iface].lup[ilup], face[iface].lup[ilup+1]-1);

            for (i = face[iface].lup[ilup]; i < face[iface].lup[ilup+1]; i++) {
                iedge = face[iface].edg[i];
                if         (iedge > 0) {
                    printf("     Edge %6d, npnt=%6d, ibeg=%3d, iend=%3d\n",
                           iedge, edge[+iedge].npnt, edge[+iedge].ibeg, edge[+iedge].iend);
                } else if (iedge < 0) {
                    printf("     Edge %6d, npnt=%6d, iend=%3d, ibeg=%3d\n",
                           iedge, edge[-iedge].npnt, edge[-iedge].iend, edge[-iedge].ibeg);
                } else {
                    printf("     Edge %6d, degenerate,   ibeg=%3d, iend=%3d\n",
                           iedge,                    edge[+iedge].ibeg, edge[+iedge].iend);
                }
            }
        }

        /* sort the Loops based upon area */
        if (face[iface].nlup > 1) {
            MALLOC(area, double, face[iface].nlup);

            for (ilup = 0; ilup < face[iface].nlup; ilup++) {
                area[ilup] = 0;

                i     = face[iface].lup[ilup];
                iedge = face[iface].edg[i];

                if        (iedge > 0) {
                    ibeg = edge[+iedge].ibeg;
                    ip0  = edge[+iedge].pnt[0];
                } else {
                    ibeg = edge[-iedge].iend;
                    ip0  = edge[-iedge].pnt[edge[-iedge].npnt-1];
                }

                while (1) {
                    if        (iedge > 0) {
                        for (j = 0; j < edge[+iedge].npnt-1; j++) {
                            ip1 = edge[+iedge].pnt[j  ];
                            ip2 = edge[+iedge].pnt[j+1];

                            areax = (tess.xyz[3*ip1+1] - tess.xyz[3*ip0+1])
                                  * (tess.xyz[3*ip2+2] - tess.xyz[3*ip0+2])
                                  - (tess.xyz[3*ip2+1] - tess.xyz[3*ip0+1])
                                  * (tess.xyz[3*ip1+2] - tess.xyz[3*ip0+2]);
                            areay = (tess.xyz[3*ip1+2] - tess.xyz[3*ip0+2])
                                  * (tess.xyz[3*ip2  ] - tess.xyz[3*ip0  ])
                                  - (tess.xyz[3*ip2+2] - tess.xyz[3*ip0+2])
                                  * (tess.xyz[3*ip1  ] - tess.xyz[3*ip0  ]);
                            areaz = (tess.xyz[3*ip1  ] - tess.xyz[3*ip0  ])
                                  * (tess.xyz[3*ip2+1] - tess.xyz[3*ip0+1])
                                  - (tess.xyz[3*ip2  ] - tess.xyz[3*ip0  ])
                                  * (tess.xyz[3*ip1+1] - tess.xyz[3*ip0+1]);

                            area[ilup] += sqrt(areax*areax + areay*areay + areaz*areaz);
                        }

                        if (edge[+iedge].iend == ibeg) break;
                    } else {
                        for (j = edge[-iedge].npnt-1; j > 0; j--) {
                            ip1 = edge[-iedge].pnt[j  ];
                            ip2 = edge[-iedge].pnt[j-1];

                            areax = (tess.xyz[3*ip1+1] - tess.xyz[3*ip0+1])
                                  * (tess.xyz[3*ip2+2] - tess.xyz[3*ip0+2])
                                  - (tess.xyz[3*ip2+1] - tess.xyz[3*ip0+1])
                                  * (tess.xyz[3*ip1+2] - tess.xyz[3*ip0+2]);
                            areay = (tess.xyz[3*ip1+2] - tess.xyz[3*ip0+2])
                                  * (tess.xyz[3*ip2  ] - tess.xyz[3*ip0  ])
                                  - (tess.xyz[3*ip2+2] - tess.xyz[3*ip0+2])
                                  * (tess.xyz[3*ip1  ] - tess.xyz[3*ip0  ]);
                            areaz = (tess.xyz[3*ip1  ] - tess.xyz[3*ip0  ])
                                  * (tess.xyz[3*ip2+1] - tess.xyz[3*ip0+1])
                                  - (tess.xyz[3*ip2  ] - tess.xyz[3*ip0  ])
                                  * (tess.xyz[3*ip1+1] - tess.xyz[3*ip0+1]);

                            area[ilup] += sqrt(areax*areax + areay*areay + areaz*areaz);
                        }

                        if (edge[-iedge].ibeg == ibeg) break;
                    }

                    if (i == face[iface].nedg-1) break;

                    iedge = face[iface].edg[++i];
                }

                i     = face[iface].lup[ilup];
                iedge = face[iface].edg[i];
                printf("   Loop %2d starts at i=%2d (iedge=%4d) and has area %f\n", ilup, i, iedge, area[ilup]);
            }

            MALLOC(edgtmp, int, face[iface].nedg+1);
            MALLOC(luptmp, int, face[iface].nlup+1);

            for (i = 0; i < face[iface].nedg; i++) {
                edgtmp[i] = face[iface].edg[i];
            }
            for (i = 0; i <= face[iface].nlup; i++) {
                luptmp[i] = face[iface].lup[i];
            }

            j = 0;
            for (ilup = 0; ilup < face[iface].nlup; ilup++) {

                klup = -1;
                amax =  0;
                for (jlup = 0; jlup < face[iface].nlup; jlup++) {
                    if (area[jlup] > amax) {
                        klup = jlup;
                        amax = area[jlup];
                    }
                }

                area[klup] = -1;        /* so that it does not get picked again */

                for (i = luptmp[klup]; i < luptmp[klup+1]; i++) {
                    face[iface].edg[j++] = edgtmp[i];
                }

                face[iface].lup[ilup+1] = j;
            }

            FREE(edgtmp);
            FREE(luptmp);

            printf("...after sorting Loops so that largest area is first\n");
            for (ilup = 0; ilup < face[iface].nlup; ilup++) {
                printf("   Loop %6d (%d:%d)\n", ilup, face[iface].lup[ilup], face[iface].lup[ilup+1]-1);

                for (i = face[iface].lup[ilup]; i < face[iface].lup[ilup+1]; i++) {
                    iedge = face[iface].edg[i];
                    if        (iedge > 0) {
                        printf("     Edge %6d, npnt=%6d, ibeg=%3d, iend=%3d\n",
                               iedge, edge[+iedge].npnt, edge[+iedge].ibeg, edge[+iedge].iend);
                    } else if (iedge < 0) {
                        printf("     Edge %6d, npnt=%6d, iend=%3d, ibeg=%3d\n",
                               iedge, edge[-iedge].npnt, edge[-iedge].iend, edge[-iedge].ibeg);
                    } else {
                        printf("     Edge %6d, degenerate,   ibeg=%3d, iend=%3d\n",
                               iedge,                    edge[+iedge].ibeg, edge[+iedge].iend);
                    }
                }
            }

            FREE(area);

            printf("there are %d Loops, so no face* file created\n", face[iface].nlup);
            continue;
        }
    }

    /* set the coordinates for each Edge Point (including bounding Nodes) */
    for (iedge = 1; iedge <= Nedge; iedge++) {
        MALLOC(edge[iedge].xyz, double, 3*edge[iedge].npnt);

        for (i = 0; i < edge[iedge].npnt; i++) {
            ipnt = edge[iedge].pnt[i];
            edge[iedge].xyz[3*i  ] = tess.xyz[3*ipnt  ];
            edge[iedge].xyz[3*i+1] = tess.xyz[3*ipnt+1];
            edge[iedge].xyz[3*i+2] = tess.xyz[3*ipnt+2];
        }
    }

    /* set the ptyp for each Node, Edge, and Face Point */
    for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
        tess.ptyp[ipnt] = 0;
    }

    for (inode = 1; inode <= Nnode; inode++) {
        ipnt = node[inode].ipnt;
        tess.ptyp[ipnt] = PNT_NODE | inode;
    }

    for (iedge = 1; iedge <= Nedge; iedge++) {
        for (i = 1; i < edge[iedge].npnt-1; i++) {
            ipnt = edge[iedge].pnt[i];
            tess.ptyp[ipnt] = PNT_EDGE | iedge;
        }
    }

    for (itri = 0; itri < tess.ntri; itri++ ) {
        if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;
        icolr = tess.ttyp[itri] & TRI_COLOR;

        for (iface = 1; iface <= Nface; iface++) {
            if (Nface == 1 || face[iface].icol == icolr) {
                ipnt = tess.trip[3*itri  ];
                if (tess.ptyp[ipnt] == 0) tess.ptyp[ipnt] = PNT_FACE | iface;

                ipnt = tess.trip[3*itri+1];
                if (tess.ptyp[ipnt] == 0) tess.ptyp[ipnt] = PNT_FACE | iface;

                ipnt = tess.trip[3*itri+2];
                if (tess.ptyp[ipnt] == 0) tess.ptyp[ipnt] = PNT_FACE | iface;

                break;
            }
        }
    }

    /* finally print out the whole structure */
    SPRINT0(1, "\nSummary of Brep\n");
    SPRINT0(1, " inode   ipnt   nedg          x          y          z");
    for (inode = 1; inode <= Nnode; inode++) {
        SPRINT6(1, "%6d %6d %6d %10.4f %10.4f %10.4f", inode,
                node[inode].ipnt,
                node[inode].nedg,
                node[inode].x,
                node[inode].y,
                node[inode].z);
    }

    SPRINT0(1, " iedge   ibeg   iend  ileft  irite   npnt");
    for (iedge = 1; iedge <= Nedge; iedge++) {
        SPRINT6(1, "%6d %6d %6d %6d %6d %6d", iedge,
                edge[iedge].ibeg,
                edge[iedge].iend,
                edge[iedge].ileft,
                edge[iedge].irite,
                edge[iedge].npnt);
    }

    SPRINT0(1, " iface   icol   nlup   nedg    edg...");
    for (iface = 1; iface <= Nface; iface++) {
        SPRINT1x(1,  "%6d",      iface      );
        SPRINT1x(1, " %6d", face[iface].icol);
        SPRINT1x(1, " %6d", face[iface].nlup);
        SPRINT1x(1, " %6d", face[iface].nedg);
        for (iedge = 0; iedge < face[iface].nedg; iedge++) {
            SPRINT1x(1, " %6d", face[iface].edg[iedge]);
        }
        SPRINT0(1, " ");
    }

cleanup:
    FREE(sgmt  );
    FREE(nodnum);
    FREE(edgtmp);
    FREE(luptmp);
    FREE(area  );

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   generateEgads - generate EGADS Brep                               *
 *                                                                     *
 *   sets up:                                                          *
 *      node[].enode                                                   *
 *      edge[].eedge                                                   *
 *      face[].eface                                                   *
 *                                                                     *
 ***********************************************************************
 */
static int
generateEgads(char   egadsname[],       /* (in)  name of file to write */
              char   message[])         /* (out) error message */
{
    int       status = SUCCESS;

    int       inode, iedge, iloop, iface, closed=1, nfaces, ncp, ndata;
    int       periodic, nnode, nedge, nface, ipnt;

    int       i, j, k, ij, header[7], *senses=NULL, oclass, mtype, nchild, *senses2;
    double    xyz[3], *cpdata=NULL, tdata[4], data[18], data2[4];
    double    rms, rmstrain, uv_out[2], xyz_out[18];
    ego       context, eref, ecurv, esurf, *eloops=NULL, eshell, ebody, emodel;
    ego       enode, enodes[2], *eedges=NULL, *efaces=NULL, *echilds;
    FILE      *fpsum=NULL;

    ROUTINE(generateEgads);

    /* -------------------------------------------------2-------------- */

    SPRINT0(1, "Generating EGADS ...");

#ifdef DEBUG
    /* print whole data structure */
    for (inode = 1; inode <=Nnode; inode++) {
        printf("inode=%3d, nedg=%3d, xyz=%10.5f %10.5f %10.5f\n", inode,
               node[inode].nedg,
               node[inode].x, node[inode].y, node[inode].z);
    }
    for (iedge = 1; iedge <= Nedge; iedge++) {
        printf("iedge=%3d, ibeg=%3d, iend=%3d, ileft=%3d, irite=%3d\n", iedge,
               edge[iedge].ibeg,  edge[iedge].iend,
               edge[iedge].ileft, edge[iedge].irite);
    }
    for (iface = 1; iface <= Nface; iface++) {
        printf("iface=%3d, icol=%3d, nedg=%3d, nlup=%3d\n", iface,
               face[iface].icol, face[iface].nedg, face[iface].nlup);
        printf("           edg=");
        for (i = 0; i < face[iface].nedg; i++) {
            printf(" %3d", face[iface].edg[i]);
        }
        printf("\n           lup=");
        for (i = 0; i < face[iface].nlup; i++) {
            printf(" %3d", face[iface].lup[i]);
        }
        printf("\n");
    }
#endif

    strcpy(message, "okay");

    status = EG_open(&context);
    CHECK_STATUS(EG_open);

    status = EG_setOutLevel(context, outLevel);
    CHECK_STATUS(EG_setOutLevel);

    /* make each of the Nodes */
    for (inode = 1; inode <= Nnode; inode++) {
        SPRINT1(1, "\n*********\nworking on inode=%d\n*********", inode);

        /* create an array of the point associated with the Node */
        xyz[0] = node[inode].x;
        xyz[1] = node[inode].y;
        xyz[2] = node[inode].z;

        /* make the Node */
        status = EG_makeTopology(context, NULL, NODE, 0,
                                 xyz, 0, NULL, NULL,
                                 &(node[inode].enode));
        SPRINT2(1, "EG_makeTopology(NODE=%2d) -> status=%d", inode, status);
        CHECK_STATUS(EG_makeTopology);
    }

    /* make each of the Edges */
    for (iedge = 1; iedge <= Nedge; iedge++) {
        SPRINT1(1, "\n*********\nworking on iedge=%d\n*********", iedge);

        ncp = edge[iedge].ncp;

        /* create the Curve */
        header[0] = 0;            // bitflag
        header[1] = 3;            // degree
        header[2] = ncp;          // number of control points
        header[3] = ncp+4;        // number of knots

        ndata = (ncp+4) + 3 * ncp;
        MALLOC(cpdata, double, ndata);

        ndata = 0;

        /* knot vector */
        cpdata[ndata++] = 0;
        cpdata[ndata++] = 0;
        cpdata[ndata++] = 0;
        cpdata[ndata++] = 0;

        for (j = 1; j < ncp-3; j++) {
            cpdata[ndata++] = j;
        }

        cpdata[ndata++] = ncp-3;
        cpdata[ndata++] = ncp-3;
        cpdata[ndata++] = ncp-3;
        cpdata[ndata++] = ncp-3;

        /* control points */
        for (j = 0; j < ncp; j++) {
            cpdata[ndata++] = edge[iedge].cp[3*j  ];
            cpdata[ndata++] = edge[iedge].cp[3*j+1];
            cpdata[ndata++] = edge[iedge].cp[3*j+2];
        }

        status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                                 header, cpdata, &ecurv);
        printf("EG_makeGeometry(CURVE) -> status=%d\n", status);
        CHECK_STATUS(EG_makeGeometry);

        FREE(cpdata);

        /* create the Edge */
        status = EG_getRange(ecurv, tdata, &periodic);
        CHECK_STATUS(EG_getRange);

        enodes[0] = node[edge[iedge].ibeg].enode;
        enodes[1] = node[edge[iedge].iend].enode;

        status = EG_makeTopology(context, ecurv, EDGE, TWONODE,
                                 tdata, 2, enodes, NULL, &(edge[iedge].eedge));
        CHECK_STATUS(EG_makeTopology);

        SPRINT1(1, "EG_makeTopology(EDGE) -> status=%d", status);
    }

    /* make each of the Faces */
    MALLOC(efaces, ego, Nface);
    nfaces = 0;

    if (subsample > 1) {
        fpsum = fopen("subsample.summary", "a");
        fprintf(fpsum, "%5d", nctrlpnt);
    }

    for (iface = 1; iface <= Nface; iface++) {
        SPRINT1(1, "\n*********\nworking on iface=%d\n*********", iface);

        ncp = face[iface].ncp;

        /* see if Points are co-planar */
        if (face[iface].cp == NULL) {
            printf("planar\n");

            MALLOC(eedges, ego, face[iface].nedg);
            MALLOC(senses, int, face[iface].nedg);
            MALLOC(eloops, ego, face[iface].nlup);

            /* make the outer Loop associated with this Face.  start with
               the list of Edges associated with the out Loop */
            j = 0;
            for (i = face[iface].lup[0]; i < face[iface].lup[1]; i++) {
                if (face[iface].edg[i] > 0) {
                    eedges[j] = edge[+face[iface].edg[i]].eedge;
                    senses[j] = SFORWARD;
                } else {
                    eedges[j] = edge[-face[iface].edg[i]].eedge;
                    senses[j] = SREVERSE;
                }
                j++;
            }

            /* make the Loop (cannot use EG_makeLoop since sense of first Edge
               might be SREVERSE) */
            status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL,
                                     j, eedges, senses, &eloops[0]);
            CHECK_STATUS(EG_makeTopology);

            /* make the inner Loops associated with this Face */
            for (iloop = 1; iloop < face[iface].nlup; iloop++) {

                /* make the list of Edges associated with this Loop */
                j = 0;
                for (i = face[iface].lup[iloop]; i < face[iface].lup[iloop+1]; i++) {
                    if (face[iface].edg[i] > 0) {
                        eedges[j] = edge[+face[iface].edg[i]].eedge;
                        senses[j] = SFORWARD;
                    } else {
                        eedges[j] = edge[-face[iface].edg[i]].eedge;
                        senses[j] = SREVERSE;
                    }
                    j++;
                }

                /* make the Loop (cannot use EG_makeLoop since sense of first Edge
                   might be SREVERSE) */
                status = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL,
                                         j, eedges, senses, &eloops[iloop]);
                CHECK_STATUS(EG_makeTopology);
            }

            /* get the plane from the first Loop */
            status = EG_getPlane(eloops[0], &esurf);
            if (status < 0) {
                printEgo(eloops[0]);
            }
            CHECK_STATUS(EG_getPlane);

            /* make the Face */
            senses[0] = SFORWARD;       /* outer Loop */
            for (i = 1; i < face[iface].nlup; i++) {
                senses[i] = SREVERSE;   /* inner Loops */
            }

            status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL,
                                         face[iface].nlup, eloops, senses, &(face[iface].eface));
            CHECK_STATUS(EG_makeTopology);

            efaces[nfaces] = face[iface].eface;
            nfaces++;

            FREE(eedges);
            FREE(senses);
            FREE(eloops);

        /* for now, this only works with Faces bounded by 4 Edges */
        } else if (face[iface].lup[1] < 2 || face[iface].lup[1] > 4) {
            SPRINT3(-1,   "ERROR:: Face %d (color %d) has %d Edges (and is expecting 2, 3, or 4)\n", iface, face[iface].icol, face[iface].lup[1]);
            for (i = 0; i < face[iface].lup[1]; i++) {
                iedge = abs(face[iface].edg[i]);
                SPRINT5(0, "        iedge=%3d, ileft=%3d (color %3d), irite=%3d (color %3d)", iedge,
                        edge[iedge].ileft, face[edge[iedge].ileft].icol,
                        edge[iedge].irite, face[edge[iedge].irite].icol);
            }
            snprintf(message, 80, "Face %d (color %d) has %d Edges (and is expecting 2, 3, or 4)\n", iface, face[iface].icol, face[iface].lup[1]);

            closed = 0;

        /* Faces with 4 Edges in outer loop */
        } else {
            printf("non-planar\n");

            printf("iface=%d, nedg=%d, nlup=%d\n", iface, face[iface].nedg, face[iface].nlup);

            MALLOC(eedges, ego, 8+2*face[iface].nedg                 );
            MALLOC(senses, int, 4+  face[iface].nedg+face[iface].nlup);
            MALLOC(eloops, ego,                      face[iface].nlup);

            /* create the Surface */
            header[0] = 0;            // bitflag
            header[1] = 3;            // u degree
            header[2] = ncp;          // u number of control points
            header[3] = ncp + 4;      // u number of knots
            header[4] = 3;            // v degree
            header[5] = ncp;          // v number of control points
            header[6] = ncp + 4;      // v number of knots

            ndata = 2 * (ncp+4) + 3 * ncp * ncp;
            MALLOC(cpdata, double, ndata);

            ndata = 0;

            /* knot vectors for u and v */
            for (i = 0; i < 2; i++) {
                cpdata[ndata++] = 0;
                cpdata[ndata++] = 0;
                cpdata[ndata++] = 0;
                cpdata[ndata++] = 0;

                for (j = 1; j < ncp-3; j++) {
                    cpdata[ndata++] = j;
                }

                cpdata[ndata++] = ncp-3;
                cpdata[ndata++] = ncp-3;
                cpdata[ndata++] = ncp-3;
                cpdata[ndata++] = ncp-3;
            }

            /* control points */
            for (ij = 0; ij < ncp*ncp; ij++) {
                cpdata[ndata++] = face[iface].cp[3*ij  ];
                cpdata[ndata++] = face[iface].cp[3*ij+1];
                cpdata[ndata++] = face[iface].cp[3*ij+2];
            }

            status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                                     header, cpdata, &esurf);
            printf("EG_makeGeometry(SURFACE) -> status=%d\n", status);
            CHECK_STATUS(EG_makeGeometry);

            FREE(cpdata);

            /* create the PCurves associated with the Edges */
            iedge = face[iface].edg[0];
            if        (iedge > 0) {
                eedges[0] = edge[+iedge].eedge;
                senses[0] = SFORWARD;
                data[0]   = 0;
                data[1]   = 0;
                data[2]   = ncp-3;
                data[3]   = 0;
            } else if (iedge < 0) {
                eedges[0] = edge[-iedge].eedge;
                senses[0] = SREVERSE;
                data[0]   = ncp-3;
                data[1]   = 0;
                data[2]   = 3-ncp;
                data[3]   = 0;
            } else {
                assert(0);
            }
            status = EG_makeGeometry(context, PCURVE, LINE, esurf, NULL, data, &(eedges[4]));
            printf("EG_makeGeometry(PCURVE 0) -> status=%d\n", status);
            CHECK_STATUS(EG_makeGeometry);

            if (face[iface].lup[1] != 2) {
                iedge = face[iface].edg[1];
                if        (iedge > 0) {
                    eedges[1] = edge[+iedge].eedge;
                    senses[1] = SFORWARD;
                    data[0]   = ncp-3;
                    data[1]   = 0;
                    data[2]   = 0;
                    data[3]   = ncp-3;
                } else if (iedge < 0) {
                    eedges[1] = edge[-iedge].eedge;
                    senses[1] = SREVERSE;
                    data[0]   = ncp-3;
                    data[1]   = ncp-3;
                    data[2]   = 0;
                    data[3]   = 3-ncp;
                } else {
                    assert(0);
                }
            } else {
                enode = NULL;
                for (inode = 1; inode <= Nnode; inode++) {
                    if (fabs(node[inode].x-face[iface].cp[3*ncp-3]) < EPS06 &&
                        fabs(node[inode].y-face[iface].cp[3*ncp-2]) < EPS06 &&
                        fabs(node[inode].z-face[iface].cp[3*ncp-1]) < EPS06   ) {
                        enode = node[inode].enode;
                        break;
                    }
                }
                if (enode == NULL) {
                    printf("could not find degeneracy\n");
                    assert(0);
                }

                data2[0] = 0;
                data2[1] = ncp-3;
                status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                                         data2, 1, &enode, NULL, &(eedges[1]));
                printf("EG_makeTopology(DEGEN  1) -> status=%d\n", status);
                CHECK_STATUS(EG_makeTopology);

                senses[1] = SFORWARD;
                data[0]   = ncp-3;
                data[1]   = 0;
                data[2]   = 0;
                data[3]   = ncp-3;
            }
            status = EG_makeGeometry(context, PCURVE, LINE, esurf, NULL, data, &(eedges[5]));
            printf("EG_makeGeometry(PCURVE 1) -> status=%d\n", status);
            CHECK_STATUS(EG_makeGeometry);

            if (face[iface].lup[1] != 2) {
                iedge = face[iface].edg[2];
            } else {
                iedge = face[iface].edg[1];
            }
            if        (iedge > 0) {
                eedges[2] = edge[+iedge].eedge;
                senses[2] = SFORWARD;
                data[0]   = ncp-3;
                data[1]   = ncp-3;
                data[2]   = 3-ncp;
                data[3]   = 0;
            } else if (iedge < 0) {
                eedges[2] = edge[-iedge].eedge;
                senses[2] = SREVERSE;
                data[0]   = 0;
                data[1]   = ncp-3;
                data[2]   = ncp-3;
                data[3]   = 0;
            } else {
                assert(0);
            }
            status = EG_makeGeometry(context, PCURVE, LINE, esurf, NULL, data, &(eedges[6]));
            printf("EG_makeGeometry(PCURVE 2) -> status=%d\n", status);
            CHECK_STATUS(EG_makeGeometry);

            if (face[iface].lup[1] != 2 && face[iface].lup[1] != 3) {
                if (face[iface].lup[1] != 2) {
                    iedge = face[iface].edg[3];
                } else {
                    iedge = face[iface].edg[2];
                }
                if        (iedge > 0) {
                    eedges[3] = edge[+iedge].eedge;
                    senses[3] = SFORWARD;
                    data[0]   = 0;
                    data[1]   = ncp-3;
                    data[2]   = 0;
                    data[3]   = 3-ncp;
                } else if (iedge < 0) {
                    eedges[3] = edge[-iedge].eedge;
                    senses[3] = SREVERSE;
                    data[0]   = 0;
                    data[1]   = 0;
                    data[2]   = 0;
                    data[3]   = ncp-3;
                } else {
                    assert(0);
                }
            } else {
                enode = NULL;
                for (inode = 1; inode <= Nnode; inode++) {
                    if (fabs(node[inode].x-face[iface].cp[0]) < EPS06 &&
                        fabs(node[inode].y-face[iface].cp[1]) < EPS06 &&
                        fabs(node[inode].z-face[iface].cp[2]) < EPS06   ) {
                        enode = node[inode].enode;
                        break;
                    }
                }
                if (enode == NULL) {
                    printf("could not find degeneracy\n");
                    assert(0);
                }

                data2[0] = 0;
                data2[1] = ncp-3;
                status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                                         data2, 1, &enode, NULL, &(eedges[3]));
                printf("EG_makeTopology(DEGEN  3) -> status=%d\n", status);
                CHECK_STATUS(EG_makeTopology);

                senses[3] = SFORWARD;
                data[0]   = 0;
                data[1]   = ncp-3;
                data[2]   = 0;
                data[3]   = 3-ncp;
            }
            status = EG_makeGeometry(context, PCURVE, LINE, esurf, NULL, data, &(eedges[7]));
            printf("EG_makeGeometry(PCURVE 3) -> status=%d\n", status);
            CHECK_STATUS(EG_makeGeometry);

            /* create the outer Loop */
            status = EG_makeTopology(context, esurf, LOOP, CLOSED, NULL, 4, eedges, senses, &(eloops[0]));
            printf("EG_makeTopology(LOOP) -> status=%d\n", status);
            CHECK_STATUS(EG_makeTopology);

            /* add any inner Loops */
            for (iloop = 1; iloop < face[iface].nlup; iloop++) {

                j = 0;
                k = face[iface].lup[iloop+1] - face[iface].lup[iloop];
                for (i = face[iface].lup[iloop]; i < face[iface].lup[iloop+1]; i++) {
                    iedge = face[iface].edg[i];
                    if        (iedge > 0) {
                        eedges[j] = edge[+iedge].eedge;
                        senses[j] = SFORWARD;
                    } else {
                        eedges[j] = edge[-iedge].eedge;
                        senses[j] = SREVERSE;
                    }

                    /* pcurves */
                    status = EG_otherCurve(esurf, eedges[j], 1e-4, &eedges[k]);
                    printf("EG_otherCurve(%d) -> status=%d\n", j, status);
                    CHECK_STATUS(EG_otherCurve);

                    j++;
                    k++;
                }

                /* make the Loop (cannot use EG_makeLoop since sense of first Edge
                   might be SREVERSE) */
                status = EG_makeTopology(context, esurf, LOOP, CLOSED, NULL,
                                         j, eedges, senses, &eloops[iloop]);
                CHECK_STATUS(EG_makeTopology);
            }

            /* create the Face */
            senses[0] = SFORWARD;       /* outer Loop */
            for (i = 1; i < face[iface].nlup; i++) {
                senses[i] = SREVERSE;   /* inner Loops */
            }

            status = EG_makeTopology(context, esurf, FACE, SFORWARD, NULL,
                                     face[iface].nlup, eloops, senses, &(face[iface].eface));
            printf("EG_makeTopology(FACE) -> status=%d\n", status);
            CHECK_STATUS(EG_makeTopology);

            efaces[nfaces] = face[iface].eface;
            nfaces++;

            FREE(eedges);
            FREE(senses);
            FREE(eloops);
        }

        if (face[iface].npnt > 0) {

            /* measure accuracy for both training and testing points */
            rmstrain = 0;
            for (ipnt = 0; ipnt < face[iface].ntrain; ipnt++) {
                status = EG_invEvaluate(face[iface].eface,
                                        &(face[iface].xyztrain[3*ipnt]), uv_out, xyz_out);
                CHECK_STATUS(EG_invEvaluate);

                rmstrain += (face[iface].xyztrain[3*ipnt  ]-xyz_out[0]) * (face[iface].xyztrain[3*ipnt  ]-xyz_out[0])
                         +  (face[iface].xyztrain[3*ipnt+1]-xyz_out[1]) * (face[iface].xyztrain[3*ipnt+1]-xyz_out[1])
                         +  (face[iface].xyztrain[3*ipnt+2]-xyz_out[2]) * (face[iface].xyztrain[3*ipnt+2]-xyz_out[2]);
            }
            rmstrain = sqrt(rmstrain / face[iface].ntrain);

            rms = 0;
            for (ipnt = 0; ipnt < face[iface].npnt; ipnt++) {
                status = EG_invEvaluate(face[iface].eface,
                                        &(face[iface].xyz[3*ipnt]), uv_out, xyz_out);
                CHECK_STATUS(EG_invEvaluate);

                rms += (face[iface].xyz[3*ipnt  ]-xyz_out[0]) * (face[iface].xyz[3*ipnt  ]-xyz_out[0])
                    +  (face[iface].xyz[3*ipnt+1]-xyz_out[1]) * (face[iface].xyz[3*ipnt+1]-xyz_out[1])
                    +  (face[iface].xyz[3*ipnt+2]-xyz_out[2]) * (face[iface].xyz[3*ipnt+2]-xyz_out[2]);
            }
            rms = sqrt(rms / face[iface].npnt);

            printf("\niface=%3d  ntrain=%6d  rms=%11.3e\n", iface, face[iface].ntrain, rmstrain);
            printf(  "           npnt  =%6d  rms=%11.3e\n",        face[iface].npnt,   rms     );

            if (fpsum != NULL && face[iface].npnt > 0) {
                fprintf(fpsum, " %3d %11.3e %11.3e", iface, rmstrain, rms);
            }
        } else {
            printf("\niface=%3d  is planar, so accuracy is not computed\n", iface);
        }
    }

    if (fpsum != NULL) {
        fprintf(fpsum, "\n");
        fclose(fpsum);
    }

    if (closed < 0) {
        SPRINT0(-1, "ERROR:: not all Faces were made");
        goto cleanup;
    }

    if (Nface == 1) {
        SPRINT0(1, "\n*********\nspecial case to make SheetBody\n*********");

        status = EG_makeTopology(context, NULL, SHELL, OPEN,
                                 NULL, 1, efaces, NULL, &eshell);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                                 NULL, 1, &eshell, NULL, &ebody);
        CHECK_STATUS(EG_makeTopology);

    } else {
        /* assemble the Faces into a Shell */
        SPRINT0(1, "\n*********\nworking on shell and body\n*********");
        printf("before makeTopology(SHELL, nfaces=%d, closed=%d)\n", nfaces, closed);
        if (closed == 1) {
            status = EG_makeTopology(context, NULL, SHELL, CLOSED,
                                     NULL, nfaces, efaces, NULL, &eshell);
            SPRINT1(1, "EG_makeTopology(SHELL) -> status=%d", status);
            CHECK_STATUS(EG_makeTopology);
        } else {
            status = EG_makeTopology(context, NULL, SHELL, OPEN,
                                     NULL, nfaces, efaces, NULL, &eshell);
            SPRINT1(1, "EG_makeTopology(SHELL) -> status=%d", status);
            CHECK_STATUS(EG_makeTopology);
        }

        /* check is shell got properly closed */
        status = EG_getTopology(eshell, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses2);
        CHECK_STATUS(EG_getTopology);

        if (closed == 1 && mtype == OPEN) {
            printf("WARNING:: expecting shell to be closed but it is open (nchild=%d)\n", nchild);
            closed = 0;
        }

        /* if Shell was closed, make a SolidBody from the Shell */
        if (closed == 1) {
            status = EG_makeTopology(context, NULL, BODY, SOLIDBODY,
                                     NULL, 1, &eshell, NULL, &ebody);
            SPRINT1(1, "EG_makeTopology(BODY) -> status=%d", status);
            CHECK_STATUS(EG_makeTopology);
        } else {
            status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                                     NULL, 1, &eshell, NULL, &ebody);
            SPRINT1(1, "EG_makeTopology(BODY) -> status=%d", status);
            CHECK_STATUS(EG_makeTopology);
        }
    }

    status = EG_getBodyTopos(ebody, NULL, NODE, &nnode, NULL);
    CHECK_STATUS(EG_getBodyTopos);
    SPRINT1(1, "Body has %5d Nodes", nnode);

    status = EG_getBodyTopos(ebody, NULL, EDGE, &nedge, NULL);
    CHECK_STATUS(EG_getBodyTopos);
    SPRINT1(1, "Body has %5d Edges", nedge);

    status = EG_getBodyTopos(ebody, NULL, FACE, &nface, NULL);
    CHECK_STATUS(EG_getBodyTopos);
    SPRINT1(1, "Body has %5d Faces", nface);

    SPRINT0(2, "\nebody");
    if (outLevel > 1) printEgo(ebody);

    /* make a Model for the SolidBody */
    status = EG_makeTopology(context, NULL, MODEL, 0,
                             NULL, 1, &ebody, NULL, &emodel);
    SPRINT1(1, "EG_makeTopology(MODEL) -> status=%d", status);
    CHECK_STATUS(EG_makeTopology);

    SPRINT0(1, "\nemodel");
    if (outLevel > 0) printEgo(emodel);

    /* write out an egads file */
    status = remove(egadsname);
    if (status == 0) {
        SPRINT1(-1, "WARNING:: file \"%s\" is being overwritten", egadsname);
    } else {
        SPRINT1(-1, "File \"%s\" is being written", egadsname);
    }

    status = EG_saveModel(emodel, egadsname);
    CHECK_STATUS(EG_saveModel);

    /* clean up */
    status = EG_close(context);
    CHECK_STATUS(EG_close);

cleanup:
    FREE(eedges);
    FREE(senses);
    FREE(eloops);
    FREE(efaces);
    FREE(cpdata);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   generateFits - generate Bspline fits for each Edge and Node       *
 *                                                                     *
 *   updates:                                                          *
 *      edge[].cp                                                      *
 *   sets up:                                                          *
 *      face[].npnt  (non-planar Faces)                                *
 *      face[].xyz   (non-planar Faces)                                *
 *      face[].ncp   (non-planar Faces)                                *
 *      face[].cp    (non-planar Faces)                                *
 *      face[].done                                                    *
 *                                                                     *
 ***********************************************************************
 */
static int
generateFits(int    ncp,                /* (in)  number of control points */
             char   message[])          /* (out) error message */
{
    int       status = SUCCESS;

    int       npnt, ntrain, ipnt, iloop, iedge, iface;
    int       i, j, ii, jj, nmin, planar, numiter, bitflag, nsample;
    double    length, smooth, maxf;
    double    xx[4], yy[4], zz[4], prod, normf, dotmin, prodmax;
    double    dx1=0, dy1=0, dz1=0, dx2=0, dy2=0, dz2=0;
    double    areax, areay, areaz, area, xdegen, ydegen, zdegen;
    double    *Tcloud=NULL, *UVcloud=NULL;

    /* variables needed for multi-threading */
    int       nthread, ithread;
    long      start;
    void      **threads=NULL;
    emp_T     empFitter;
    clock_t   old_time, new_time;

    ROUTINE(generateFits);

    /* --------------------------------------------------------------- */

    SPRINT1(1, "\nGenerating Fits (ncp=%d) ...", ncp);

    strcpy(message, "okay");

    /* make each of the Edges */
    old_time = clock();

    for (iedge = 1; iedge <= Nedge; iedge++) {
        length = 0;
        for (ipnt = 1; ipnt < edge[iedge].npnt; ipnt++) {
            length += sqrt((edge[iedge].xyz[3*ipnt-3]-edge[iedge].xyz[3*ipnt  ])
                          *(edge[iedge].xyz[3*ipnt-3]-edge[iedge].xyz[3*ipnt  ])
                          +(edge[iedge].xyz[3*ipnt-2]-edge[iedge].xyz[3*ipnt+1])
                          *(edge[iedge].xyz[3*ipnt-2]-edge[iedge].xyz[3*ipnt+1])
                          +(edge[iedge].xyz[3*ipnt-1]-edge[iedge].xyz[3*ipnt+2])
                          *(edge[iedge].xyz[3*ipnt-1]-edge[iedge].xyz[3*ipnt+2]));
        }

        SPRINT3(1, "\n*********\nfitting iedge=%d (npnt=%3d, length=%10.5f)\n*********", iedge, edge[iedge].npnt, length);

        npnt = edge[iedge].npnt;

        /* allocate space for control points */
        edge[iedge].ncp = ncp;
        FREE(  edge[iedge].cp);
        MALLOC(edge[iedge].cp, double, 3*ncp);

        /* set the control points at its boundaries */
        edge[iedge].cp[0      ] = edge[iedge].xyz[0       ];
        edge[iedge].cp[1      ] = edge[iedge].xyz[1       ];
        edge[iedge].cp[2      ] = edge[iedge].xyz[2       ];

        edge[iedge].cp[3*ncp-3] = edge[iedge].xyz[3*npnt-3];
        edge[iedge].cp[3*ncp-2] = edge[iedge].xyz[3*npnt-2];
        edge[iedge].cp[3*ncp-1] = edge[iedge].xyz[3*npnt-1];

        /* allocate space for the T-parameters associated with the cloud */
        MALLOC(Tcloud, double, npnt);

        /* fit the cloud of points  with ncp control points */
        bitflag = 1;
        smooth  = 1;
        numiter = 1000;
        status = fit1dCloud(npnt, bitflag, edge[iedge].xyz, ncp, edge[iedge].cp, smooth,
                            Tcloud, &normf, &maxf, &dotmin, &nmin, &numiter, stdout);
        printf("fit1dCloud(npnt=%d, ncp=%d) -> status=%d,  numiter=%3d,  normf=%12.4e,  dotmin=%.4f,  nmin=%d\n",
               npnt, ncp, status, numiter, normf, dotmin, nmin);
        CHECK_STATUS(fit1dCloud);

#ifdef GRAFIC
        /* plot the fit */
        status = plotCurve(npnt, edge[iedge].xyz, Tcloud,
                           ncp,  edge[iedge].cp, normf, dotmin, nmin);
        printf("plotCurve -> status=%d\n", status);
        CHECK_STATUS(plotCurve);
#endif

        FREE(Tcloud);
    }

    new_time = clock();
    printf("generateFits(1D), CPU=%10.2f sec\n",
           (double)(new_time-old_time)/(double)(CLOCKS_PER_SEC));

    /* make each of the Faces */
    for (iface = 1; iface <= Nface; iface++) {
        SPRINT1(1, "\n*********\ninitializing iface=%d\n*********", iface);

        face[iface].done = 0;

        /* count the number of points associated with this Face */
        npnt = 0;
        for (iloop = 0; iloop < face[iface].nlup; iloop++) {
            for (i = face[iface].lup[iloop] ; i < face[iface].lup[iloop+1]; i++) {
                iedge = abs(face[iface].edg[i]);
                npnt += (edge[iedge].npnt - 1);
            }
        }
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (tess.ptyp[ipnt] == (PNT_FACE | iface)) {
                npnt++;
            }
        }

        assert(npnt > 0);

        /* determine size of training set (try to have at least 2 points per
           control net cell) */
        nsample = MIN(MAX(1, npnt / (2 * ncp * ncp)), subsample);

        /* create an array of the discrete points associated with this Face */
        face[iface].npnt = npnt;
        FREE(  face[iface].xyz);
        MALLOC(face[iface].xyz, double, 3*npnt);
        MALLOC(face[iface].xyztrain, double, 3*npnt);

        npnt = 0;
        ntrain = 0;
        for (iloop = 0; iloop < face[iface].nlup; iloop++) {
            for (i = face[iface].lup[iloop]; i < face[iface].lup[iloop+1]; i++) {
                iedge = face[iface].edg[i];
                if (iedge > 0) {
                    for (ipnt = 0; ipnt < edge[+iedge].npnt-1; ipnt++) {
                        face[iface].xyz[3*npnt  ] = edge[+iedge].xyz[3*ipnt  ];
                        face[iface].xyz[3*npnt+1] = edge[+iedge].xyz[3*ipnt+1];
                        face[iface].xyz[3*npnt+2] = edge[+iedge].xyz[3*ipnt+2];
                        npnt++;

                        if (rand()%nsample == 0) {
                            face[iface].xyztrain[3*ntrain  ] = edge[+iedge].xyz[3*ipnt  ];
                            face[iface].xyztrain[3*ntrain+1] = edge[+iedge].xyz[3*ipnt+1];
                            face[iface].xyztrain[3*ntrain+2] = edge[+iedge].xyz[3*ipnt+2];
                            ntrain++;
                        }
                    }
                } else if (iedge < 0) {
                    for (ipnt = edge[-iedge].npnt-1; ipnt > 0; ipnt--) {
                        face[iface].xyz[3*npnt  ] = edge[-iedge].xyz[3*ipnt  ];
                        face[iface].xyz[3*npnt+1] = edge[-iedge].xyz[3*ipnt+1];
                        face[iface].xyz[3*npnt+2] = edge[-iedge].xyz[3*ipnt+2];
                        npnt++;

                        if (rand()%nsample == 0) {
                            face[iface].xyztrain[3*ntrain  ] = edge[-iedge].xyz[3*ipnt  ];
                            face[iface].xyztrain[3*ntrain+1] = edge[-iedge].xyz[3*ipnt+1];
                            face[iface].xyztrain[3*ntrain+2] = edge[-iedge].xyz[3*ipnt+2];
                            ntrain++;
                        }
                    }
                }
            }
        }
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (tess.ptyp[ipnt] == (PNT_FACE | iface)) {
                face[iface].xyz[3*npnt  ] = tess.xyz[3*ipnt  ];
                face[iface].xyz[3*npnt+1] = tess.xyz[3*ipnt+1];
                face[iface].xyz[3*npnt+2] = tess.xyz[3*ipnt+2];
                npnt++;

                if (rand()%nsample == 0) {
                    face[iface].xyztrain[3*ntrain  ] = tess.xyz[3*ipnt  ];
                    face[iface].xyztrain[3*ntrain+1] = tess.xyz[3*ipnt+1];
                    face[iface].xyztrain[3*ntrain+2] = tess.xyz[3*ipnt+2];
                    ntrain++;
                }
            }
        }
        assert(npnt == face[iface].npnt);

        face[iface].ntrain = ntrain;

        printf("iface=%3d   npnt  =%6d\n",  iface, face[iface].npnt);
        printf("            ntrain=%6d (%3d%%)\n", face[iface].ntrain, 100*face[iface].ntrain/face[iface].npnt);

        /* create a plot of the "training" and "testing" points */
#ifdef GRAFIC
        {
            int    indgr=1+2+4+16+64+1024, itype=0;
            char   pltitl[80];

            sprintf(pltitl, "~u~v~Face %d", iface);

            grctrl_(plotPoints_image, &indgr, pltitl,
                    (void*)(&itype),
                    (void*)(&(face[iface].ntrain)),
                    (void*)(face[iface].xyztrain),
                    (void*)(&(face[iface].npnt)),
                    (void*)(face[iface].xyz),
                    (void*)NULL,
                    (void*)NULL,
                    (void*)NULL,
                    (void*)NULL,
                    (void*)NULL,
                    strlen(pltitl));
        }
#endif

        /* see if Points are co-planar by computing triple product
           of 4 Points at a time */
        planar = 1;

        ipnt = 0;
        xx[0] = face[iface].xyz[3*ipnt  ];
        yy[0] = face[iface].xyz[3*ipnt+1];
        zz[0] = face[iface].xyz[3*ipnt+2];

        ipnt = 1 * face[iface].npnt / 3;
        xx[1] = face[iface].xyz[3*ipnt  ];
        yy[1] = face[iface].xyz[3*ipnt+1];
        zz[1] = face[iface].xyz[3*ipnt+2];

        ipnt = 2 * face[iface].npnt / 3;
        xx[2] = face[iface].xyz[3*ipnt  ];
        yy[2] = face[iface].xyz[3*ipnt+1];
        zz[2] = face[iface].xyz[3*ipnt+2];

        dx1 = xx[1] - xx[0];
        dy1 = yy[1] - yy[0];
        dz1 = zz[1] - zz[0];
        dx2 = xx[2] - xx[0];
        dy2 = yy[2] - yy[0];
        dz2 = zz[2] - zz[0];

        areax = dy1 * dz2 - dz1 * dy2;
        areay = dz1 * dx2 - dx1 * dz2;
        areaz = dx1 * dy2 - dy1 * dx2;
        area  = sqrt(areax * areax + areay * areay + areaz * areaz);
        if (area < EPS06) {
            int iii;
            printf("points are colinear (%d, %d, %d)\n", 0, face[iface].npnt/3, 2*face[iface].npnt/3);
            for (iii = 0; iii < face[iface].npnt; iii++) {
                printf("%5d  %15.7f %15.7f %15.7f\n", iii, face[iface].xyz[3*iii], face[iface].xyz[3*iii+1], face[iface].xyz[3*iii+2]);
            }
            exit(0);
        } else {
            areax /= area;
            areay /= area;
            areaz /= area;
        }

        prodmax = 0;
        for (ipnt = 0; ipnt < face[iface].npnt; ipnt++) {
            xx[3] = face[iface].xyz[3*ipnt  ];
            yy[3] = face[iface].xyz[3*ipnt+1];
            zz[3] = face[iface].xyz[3*ipnt+2];

            /* if triple product is not zero, Points are not planar */
            prod = (xx[3] - xx[0]) * areax + (yy[3] - yy[0]) * areay + (zz[3] - zz[0]) * areaz;
            if (fabs(prod) > prodmax) prodmax = fabs(prod);
            if (fabs(prod) > EPS03) {
                printf("non-planar (fitting)  ipnt=%d, prod=%12.5f\n", ipnt, prod);
                planar = 0;
                break;
            }
        }

        if (planar == 1) {
            SPRINT1(1, "planar (skipping)  prodmax=%12.5f", prodmax);

            FREE(face[iface].xyz);
            face[iface].npnt = 0;

            continue;
        }

        /* for now, this only works with Faces bounded by 2, 3, or 4 Edges */
        if (face[iface].lup[1] < 2 || face[iface].lup[1] > 4) {
            printf("in generateFits\n");
            SPRINT3(-1,   "ERROR:: Face %d (color %d) has %d Edges (and is expecting 2, 3, or 4)",
                    iface, face[iface].icol, face[iface].lup[1]);
            for (i = 0; i < face[iface].lup[1]; i++) {
                iedge = abs(face[iface].edg[i]);
                SPRINT5(0, "        iedge=%3d, ileft=%3d (color %3d), irite=%3d (color %3d)", iedge,
                        edge[iedge].ileft, face[edge[iedge].ileft].icol,
                        edge[iedge].irite, face[edge[iedge].irite].icol);
            }
            snprintf(message, 80, "Face %d (color %d) has %d Edges (and is expecting 2, 3, or 4)",
                     iface, face[iface].icol, face[iface].lup[1]);

            FREE(face[iface].xyz);
            face[iface].npnt = 0;

            continue;
        }

        /* allocate space for control points */
        face[iface].ncp = ncp;
        FREE(  face[iface].cp);
        MALLOC(face[iface].cp, double, 3*ncp*ncp);

        /* set the control points at its boundaries */
        // south
        iedge = face[iface].edg[0];
        printf("extracting south control points from iedge=%5d\n", iedge);
        if (edge[abs(iedge)].ncp != ncp) {
            printf("mismatch 0\n");
        } else if (iedge > 0) {
            j  = 0;
            for (i = 0; i < ncp; i++) {
                ii = i;
                face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*ii  ];
                face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*ii+1];
                face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*ii+2];
            }
        } else {
            iedge = -iedge;
            j  = 0;
            for (i = 0; i < ncp; i++) {
                ii = ncp - 1 - i;
                face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*ii  ];
                face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*ii+1];
                face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*ii+2];
            }
        }

        // east
        if (face[iface].lup[1] == 3 || face[iface].lup[1] == 4) {
            iedge = face[iface].edg[1];
            printf("extracting east  control points from iedge=%5d\n", iedge);
            if (edge[abs(iedge)].ncp != ncp) {
                printf("mismatch 1\n");
            } else if (iedge > 0) {
                i = ncp - 1;
                for (j = 0; j < ncp; j++) {
                    jj = j;
                    face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*jj  ];
                    face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*jj+1];
                    face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*jj+2];
                }
            } else {
                iedge = -iedge;
                i = ncp - 1;
                for (j = 0; j < ncp; j++) {
                    jj = ncp - 1 - j;
                    face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*jj  ];
                    face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*jj+1];
                    face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*jj+2];
                }
            }
        } else {
            i = ncp - 1;
            xdegen = face[iface].cp[3*(i+ncp*0)  ];
            ydegen = face[iface].cp[3*(i+ncp*0)+1];
            zdegen = face[iface].cp[3*(i+ncp*0)+2];
            printf("copying    east  control points from degen     %10.4f %10.4f %10.4f\n",
                   xdegen, ydegen, zdegen);

            for (j = 1; j < ncp; j++) {
                face[iface].cp[3*(i+ncp*j)  ] = xdegen;
                face[iface].cp[3*(i+ncp*j)+1] = ydegen;
                face[iface].cp[3*(i+ncp*j)+2] = zdegen;
            }

            /* remove points that are at the degeneracy */
            for (ipnt = npnt-1; ipnt >= 0; ipnt--) {
                if (fabs(face[iface].xyz[3*ipnt  ]-xdegen) < EPS06 &&
                    fabs(face[iface].xyz[3*ipnt+1]-ydegen) < EPS06 &&
                    fabs(face[iface].xyz[3*ipnt+2]-zdegen) < EPS06   ) {
                    npnt--;
                    face[iface].xyz[3*ipnt  ] = face[iface].xyz[3*npnt  ];
                    face[iface].xyz[3*ipnt+1] = face[iface].xyz[3*npnt+1];
                    face[iface].xyz[3*ipnt+2] = face[iface].xyz[3*npnt+2];
                }
            }
            face[iface].npnt = npnt;
            if (face[iface].ntrain > face[iface].npnt) {
                face[iface].ntrain = face[iface].npnt;
            }
        }

        // north
        if (face[iface].lup[1] != 2) {
            iedge = face[iface].edg[2];
        } else {
            iedge = face[iface].edg[1];
        }
        printf("extracting north control points from iedge=%5d\n", iedge);
        if (edge[abs(iedge)].ncp != ncp) {
            printf("mismatch 2\n");
        } else if (iedge > 0) {
            j = ncp - 1;
            for (i = 0; i < ncp; i++) {
                ii = ncp - 1 - i;
                face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*ii  ];
                face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*ii+1];
                face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*ii+2];
            }
        } else {
            iedge = -iedge;
            j = ncp - 1;
            for (i = 0; i < ncp; i++) {
                ii = i;
                face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*ii  ];
                face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*ii+1];
                face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*ii+2];
            }
        }

        // west
        if (face[iface].lup[1] == 4) {
            iedge = face[iface].edg[3];
            printf("extracting west  control points from iedge=%5d\n", iedge);
            if (edge[abs(iedge)].ncp != ncp) {
                printf("mismatch 3\n");
            } else if (iedge > 0) {
                i = 0;
                for (j = 0; j < ncp; j++) {
                    jj = ncp - 1 - j;
                    face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*jj  ];
                    face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*jj+1];
                    face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*jj+2];
                }
            } else {
                iedge = -iedge;
                i = 0;
                for (j = 0; j < ncp; j++) {
                    jj = j;
                    face[iface].cp[3*(i+ncp*j)  ] = edge[iedge].cp[3*jj  ];
                    face[iface].cp[3*(i+ncp*j)+1] = edge[iedge].cp[3*jj+1];
                    face[iface].cp[3*(i+ncp*j)+2] = edge[iedge].cp[3*jj+2];
                }
            }
        } else {
#ifndef __clang_analyzer__
            i = 0;
            xdegen = face[iface].cp[3*(i+ncp*0)  ];
            ydegen = face[iface].cp[3*(i+ncp*0)+1];
            zdegen = face[iface].cp[3*(i+ncp*0)+2];
            printf("copying    west  control points from degen     %10.4f %10.4f %10.4f\n",
                   xdegen, ydegen, zdegen);

            for (j = 1; j < ncp; j++) {
                face[iface].cp[3*(i+ncp*j)  ] = xdegen;
                face[iface].cp[3*(i+ncp*j)+1] = ydegen;
                face[iface].cp[3*(i+ncp*j)+2] = zdegen;
            }

            /* remove points that are at the degeneracy */
            for (ipnt = npnt-1; ipnt >= 0; ipnt--) {
                if (fabs(face[iface].xyz[3*ipnt  ]-xdegen) < EPS06 &&
                    fabs(face[iface].xyz[3*ipnt+1]-ydegen) < EPS06 &&
                    fabs(face[iface].xyz[3*ipnt+2]-zdegen) < EPS06   ) {
                    npnt--;
                    face[iface].xyz[3*ipnt  ] = face[iface].xyz[3*npnt  ];
                    face[iface].xyz[3*ipnt+1] = face[iface].xyz[3*npnt+1];
                    face[iface].xyz[3*ipnt+2] = face[iface].xyz[3*npnt+2];
                }
            }
            face[iface].npnt = npnt;
            if (face[iface].ntrain > face[iface].npnt) {
                face[iface].ntrain = face[iface].npnt;
            }
#endif
        }
    }

    /* set up for multi-threading */
    empFitter.mutex  = NULL;
    empFitter.master = EMP_ThreadID();

    nthread = EMP_Init(&start);
    if (nthread > 4) nthread = 4;       // hyper-threading does not work on OSX

    SPRINT1(1, "\n*********\nstarting multi-threaded fits with %d threads\n*********", nthread);

    old_time = clock();

    if (nthread > 1) {

        /* create the mutex to handle list synchronization */
        empFitter.mutex = EMP_LockCreate();
        if (empFitter.mutex == NULL) {
            printf("EMPerror:: mutex creation = NULL\n");
            nthread = 1;
        } else {

            /* get storage for extra threads */
            threads = (void**) malloc((nthread-1)*sizeof(void*));
            if (threads == NULL) {
                EMP_LockDestroy(empFitter.mutex);
                nthread = 1;
            }
        }
    }

    /* create the threads and get going */
    if (threads != NULL) {
        for (ithread = 0; ithread < nthread-1; ithread++) {
            threads[ithread] = EMP_ThreadCreate(empFit2dCloud, &empFitter);
            if (threads[ithread] == NULL) {
                printf("EMPerror:: creating thread %d", ithread+1);
            }
        }
    }

    /* now run the fitter from the original thread */
    empFit2dCloud(&empFitter);

    /* wait for all others to return */
    if (threads != NULL) {
        for (ithread = 0; ithread < nthread-1; ithread++) {
            if (threads[ithread] != NULL) {
                EMP_ThreadWait(threads[ithread]);
            }
        }
    }

    /* cleanup the threads */
    if (threads != NULL) {
        for (ithread = 0; ithread < nthread-1; ithread++) {
            if (threads[ithread] != NULL) {
                EMP_ThreadDestroy(threads[ithread]);
            }
        }
    }

    if (empFitter.mutex != NULL) {
        EMP_LockDestroy(empFitter.mutex);
    }

    if (threads != NULL) {
        free(threads);
    }

    new_time = clock();
    printf("generateFits(2D), CPU=%10.2f sec\n",
           (double)(new_time-old_time)/(double)(CLOCKS_PER_SEC));

cleanup:
    FREE(Tcloud );
    FREE(UVcloud);

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   empFit2dCloud - does fitting in a thread                          *
 *                                                                     *
 ***********************************************************************
 */
static void
empFit2dCloud(void *struc)              /* (both) emp structure */
{
    int    status, iface, jface, numiter, nmin, nmax, bitflag;
    long   ID;
    double *UVcloud=NULL, smooth, maxf, normf;

    emp_T  *empFitter = (emp_T *)struc;

    ROUTINE(empFit2dCloud);

    /* --------------------------------------------------------------- */

    /* get our identifier */
    ID = EMP_ThreadID();
    if (ID == empFitter->master) {
        printf("ID %12lx: is master\n", ID);
    } else {
        printf("ID %12lx: start thread\n", ID);
    }

    /* look for next Face to process */
    while (1) {

        /* figure out which Face to do.  use a mutex to make sure
         that only one thread does one Face */
        if (empFitter->mutex != NULL) {
            EMP_LockSet(empFitter->mutex);
        }

        /* look for the Face with the most Points in the cloud */
        iface = -1;
        nmax  = -1;

        for (jface = 1; jface <= Nface; jface++) {
            if (face[jface].done == 0) {
                if (face[jface].npnt > nmax) {
                    iface = jface;
                    nmax  = face[jface].npnt;
                }
            }
        }

        /* release the mutex */
        if (empFitter->mutex != NULL) {
            EMP_LockRelease(empFitter->mutex);
        }

        /* if there are no Faces to process, break out of while(1) loop */
        if (iface < 1) break;

        /* mark this Face as being done (so that it does not get
           scheduled again) */
        face[iface].done = 1;
        printf("ID %12lx: iface %3d has %5d training points\n", ID, iface, face[iface].ntrain);

        /* do not process planar faces */
        if (face[iface].npnt == 0) {
            printf("ID %12lx: iface %3d skipped\n", ID, iface);
            continue;
        }

        /* fit the data */
        MALLOC(UVcloud, double, 2*face[iface].npnt);

        bitflag = 0;
        smooth  = 1;
        numiter = 100;

        if (outLevel > 1) {
            status = fit2dCloud(face[iface].ntrain,
                                bitflag,
                                face[iface].xyztrain,
                                face[iface].ncp,
                                face[iface].ncp,
                                face[iface].cp,
                                smooth,
                                UVcloud, &normf, &maxf, &nmin, &numiter, stdout);
        } else {
            status = fit2dCloud(face[iface].ntrain,
                                bitflag,
                                face[iface].xyztrain,
                                face[iface].ncp,
                                face[iface].ncp,
                                face[iface].cp,
                                smooth,
                                UVcloud, &normf, &maxf, &nmin, &numiter, NULL);
        }

        printf("ID %12lx: iface %3d complete with status=%d, numiter=%3d, normf=%12.4e, nmin=%d\n",
               ID, iface, status, numiter, normf, nmin);

#ifdef GRAFIC
        /* plot the fit */
        if (ID == empFitter->master) {
            status = plotSurface(face[iface].npnt, face[iface].xyz, UVcloud,
                                 face[iface].ncp,  face[iface].cp, normf, nmin);
            printf("plotSurface -> status=%d\n", status);
            CHECK_STATUS(plotSurface);
        }
#endif

        FREE(UVcloud);
    }


    /* all work is done, so exit (if not the master) */
    if (ID != empFitter->master) {
        printf("ID %12lx: stop  thread\n", ID);
        EMP_ThreadExit();
    }

cleanup:
    return;
}


/*
 ***********************************************************************
 *                                                                     *
 *   makeNodesAndEdges - make Nodes and Edges from sgements            *
 *                                                                     *
 ***********************************************************************
 */
static int
makeNodesAndEdges(int    nsgmt,
                  sgmt_T sgmt[],
                  int    ibeg,
                  int    iend,
                  int    nodnum[],
                  int    icolr,
                  int    jcolr)
{
    int status = SUCCESS;

    int   ipnt, isgmt, iface;

    ROUTINE(makeNodesAndEdges);

    /* --------------------------------------------------------------- */

    /* create a Node at the beginning if there is not one there already */
    ipnt= sgmt[ibeg].ibeg;
    if (nodnum[ipnt] < 0) {
        nodnum[ipnt] = Nnode + 1;

        if (Nnode >= Mnode) {
            Mnode += 100;
            RALLOC(node, node_T, Mnode+1);
        }

        Nnode++;

        node[Nnode].ipnt  = ipnt;
        node[Nnode].nedg  = 0;
        node[Nnode].x     = tess.xyz[3*ipnt  ];
        node[Nnode].y     = tess.xyz[3*ipnt+1];
        node[Nnode].z     = tess.xyz[3*ipnt+2];
        node[Nnode].enode = NULL;

        SPRINT6(1, "   created Node %3d .ipnt=%6d, .nedg=%6d, .x=%10.4f, .y=%10.4f, .z=%10.4f", Nnode,
                node[Nnode].ipnt,
                node[Nnode].nedg,
                node[Nnode].x,
                node[Nnode].y,
                node[Nnode].z);
    }

    /* create a Node at the end if there is not one there already */
    ipnt= sgmt[iend].iend;
    if (nodnum[ipnt] < 0) {
        nodnum[ipnt] = Nnode + 1;

        if (Nnode >= Mnode) {
            Mnode += 100;
            RALLOC(node, node_T, Mnode+1);
        }

        Nnode++;

        node[Nnode].ipnt  = ipnt;
        node[Nnode].nedg  = 0;
        node[Nnode].x     = tess.xyz[3*ipnt  ];
        node[Nnode].y     = tess.xyz[3*ipnt+1];
        node[Nnode].z     = tess.xyz[3*ipnt+2];
        node[Nnode].enode = NULL;

        SPRINT6(1, "   created Node %3d .ipnt=%6d, .nedg=%6d, .x=%10.4f, .y=%10.4f, .z=%10.4f", Nnode,
                node[Nnode].ipnt,
                node[Nnode].nedg,
                node[Nnode].x,
                node[Nnode].y,
                node[Nnode].z);
    }

    /* create the Edge */
    if (Nedge >= Medge) {
        Medge += 100;
        RALLOC(edge, edge_T, Medge+1);
    }

    Nedge++;

    edge[Nedge].ibeg  = nodnum[sgmt[ibeg].ibeg];
    edge[Nedge].iend  = nodnum[sgmt[iend].iend];
    edge[Nedge].ileft = 0;
    edge[Nedge].irite = 0;
    edge[Nedge].npnt  = nsgmt + 1;
    edge[Nedge].pnt   = NULL;
    edge[Nedge].xyz   = NULL;
    edge[Nedge].ncp   = 0;
    edge[Nedge].cp    = NULL;
    edge[Nedge].eedge = NULL;

    node[edge[Nedge].ibeg].nedg += 1;
    node[edge[Nedge].iend].nedg += 1;

    for (iface = 1; iface <= Nface; iface++) {
        if        (face[iface].icol == icolr) {
            edge[Nedge].ileft = iface;
            face[iface].nedg += 1;
        } else if (face[iface].icol == jcolr) {
            edge[Nedge].irite = iface;
            face[iface].nedg += 1;
        }
    }

    MALLOC(edge[Nedge].pnt, int, nsgmt+1);

    isgmt = ibeg;
    edge[Nedge].pnt[0] = sgmt[isgmt].ibeg;

    for (ipnt = 1; ipnt <= nsgmt; ipnt++) {
        edge[Nedge].pnt[ipnt] = sgmt[isgmt].iend;
        isgmt = sgmt[isgmt].next;
        if (isgmt < 0) {
            edge[Nedge].npnt = ipnt + 1;
            break;
        }
    }

    SPRINT6(1, "   created Edge %3d .ibeg=%6d, .iend=%6d, .ileft=%4d, .irite=%4d, .npnt=%6d", Nedge,
            edge[Nedge].ibeg,
            edge[Nedge].iend,
            edge[Nedge].ileft,
            edge[Nedge].irite,
            edge[Nedge].npnt);

cleanup:
return status;
}



/*
 ***********************************************************************
 *                                                                     *
 *   browserMessage - called when client sends a message to the server *
 *                                                                     *
 ***********************************************************************
 */
void
browserMessage(
  /*@unused@*/ void    *userPtr,
               void    *wsi,
               char    *text,
  /*@unused@*/ int     lena)
{
    int       status;

    ROUTINE(browserMessage);

    /* --------------------------------------------------------------- */

    /* process the Message */
    processMessage(text);

    /* send the response */
    SPRINT1(2, "response-> %s", response);
    wv_sendText(wsi, response);

    if (strlen(sgFocusData) > 0) {
        SPRINT1(3, "sgFocus-> %s", sgFocusData);
        wv_sendText(wsi, sgFocusData);

        /* nullify meta data so that it does not get sent again */
        sgFocusData[0] = '\0';
    }

    /* there is no key */
    status = wv_setKey(cntxt,   0, NULL,      lims[0], lims[1], NULL            );
    if (status != SUCCESS) {
        SPRINT1(3, "wv_setKet -> status=%d", status);
    }
}


/*
 ***********************************************************************
 *                                                                     *
 *   processMessage - process the message and create the response      *
 *                                                                     *
 ***********************************************************************
 */
static void
processMessage(char    *text)
{
    char      *arg1=NULL, *arg2=NULL, *arg3=NULL;
    char      *arg4=NULL, *arg5=NULL;
    char      filename[257];

    int       status, ipnt, itri;

    ROUTINE(processMessage);

#define CHECK_FOR_ERROR(COMMAND,SUBCOMMAND) \
    if (status < 0) {                                                                                \
        snprintf(response, max_resp_len, "ERROR:: %s :: %status=%d", #COMMAND, #SUBCOMMAND, status); \
        response_len = strlen(response);                                                             \
        goto cleanup;                                                                                \
    }
#define CHECK_FOR_ERROR2(COMMAND,MESSAGE)                                        \
    if (status < 0) {                                                            \
        snprintf(response, max_resp_len, "ERROR:: %s -> %s", #COMMAND, MESSAGE); \
        response_len = strlen(response);                                         \
        goto cleanup;                                                            \
    }


    /* --------------------------------------------------------------- */

    MALLOC(arg1, char, MAX_EXPR_LEN);
    MALLOC(arg2, char, MAX_EXPR_LEN);
    MALLOC(arg3, char, MAX_EXPR_LEN);
    MALLOC(arg4, char, MAX_EXPR_LEN);
    MALLOC(arg5, char, MAX_EXPR_LEN);

    SPRINT1(1, "==> processMessage(text=%s)", text);

    /* initialize the response */
    response_len = 0;
    response[0] = '\0';

    /* NO-OP */
    if (strlen(text) == 0) {

    /* "# casename" */
    } else if (strncmp(text, "# casename=", 11) == 0) {

        /* do nothing, since this is added when jrnl_out is opened */

    /* "# comment" */
    } else if (strncmp(text, "#", 1) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

    /* "identify;" */
    } else if (strncmp(text, "identify;", 9) == 0) {

        /* build the response */
        snprintf(response, max_resp_len, "identify;Slugs;");
        response_len = strlen(response);

    /* "automaticLinks;" */
    } else if (strncmp(text, "automaticLinks;", 16) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(automaticLinks, storeUndo);

        /* remove old Links and create new ones between colors */
        status = removeLinks(&tess);
        CHECK_FOR_ERROR(automaticLinks, removeLinks);

        status = makeLinks(&tess);
        CHECK_FOR_ERROR(automaticLinks, makeLinks);

        /* build the response */
        snprintf(response, max_resp_len, "automaticLinks;okay");
        response_len = strlen(response);

        /* set flags */
        Tris_pend  = 1;
        Links_pend = 1;

    /* "bridgeToPoint;x;y;z;" */
    } else if (strncmp(text, "bridgeToPoint;", 14) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(bridgeToPoint, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc  = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc  = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc  = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Triangle */
        itri = closestTriangle(xloc, yloc, zloc);
        SPRINT1(3, "closestTriangle -> itri=%d", itri);

        /* bridge between Triangle and current Point */
        status = bridgeToPoint(&tess, itri, CurPt_index);
        CHECK_FOR_ERROR(bridgeToPoint, bridgeToPoint);

        /* reset the current Point */
        CurPt_index = -1;

        /* build the response */
        snprintf(response, max_resp_len, "bridgeToPoint;okay");
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Hangs_pend = 1;
        Links_pend = 1;
        CurPt_pend = 1;

    /* "colorTriangles;x;y;z;icolor;" */
    } else if (strncmp(text, "colorTriangles;", 15) == 0) {
        int    icolr=0;
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(colorTriangles, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc  = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc  = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc  = strtod(arg3, NULL);
        if (getToken(text, 4, arg4)) icolr = strtol(arg4, NULL, 10);
        SPRINT4(3, "xloc=%f  yloc=%f  zloc=%f  icolr=%d", xloc, yloc, zloc, icolr);

        /* find the closest Triangle */
        itri = closestTriangle(xloc, yloc, zloc);
        SPRINT1(3, "closestTriangle -> itri=%d", itri);

        /* color the Triangle and its neighbors */
        status = colorTriangles(&tess, itri, icolr);
        CHECK_FOR_ERROR(colorTriangles, colorTriangles);

        /* update the latest color number */
        tess.ncolr = MAX(tess.ncolr, icolr);

        /* build the response */
        snprintf(response, max_resp_len, "colorTriangles;%d;okay", tess.ncolr);
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Links_pend = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "cutTriangles;icolr;itype;xloc;yloc;zloc;" */
    } else if (strncmp(text, "cutTriangles;", 13) == 0) {
        int    jpnt, icolr=-1, itype=0;
        double xloc=0, yloc=0, zloc=0, data[4];

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(cutTriangles, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) icolr = strtol(arg1, NULL, 10);
        if (getToken(text, 2, arg2)) itype = strtol(arg2, NULL, 10);
        if (getToken(text, 3, arg3)) xloc  = strtod(arg3, NULL);
        if (getToken(text, 4, arg4)) yloc  = strtod(arg4, NULL);
        if (getToken(text, 5, arg5)) zloc  = strtod(arg5, NULL);
        SPRINT5(3, "icolr=%d, type=%d, xloc=%f, yloc=%f, zloc=%f", icolr, itype, xloc, yloc, zloc);

        /* find the closest Point */
        ipnt = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> ipnt=%d", ipnt);
        jpnt = CurPt_index;
        SPRINT1(3, "currentPoint -> jpnt=%d", jpnt);

        /* set up data matrix that defines cutting plane
           independent of Y and Z (constant X) */
        if        (itype == 0) {
            data[0] = tess.xyz[3*ipnt  ];
            data[1] = -1;
            data[2] =  0;
            data[3] =  0;
        /* independent of X and Z (constant Y) */
        } else if (itype == 1) {
            data[0] = tess.xyz[3*ipnt+1];
            data[1] =  0;
            data[2] = -1;
            data[3] =  0;
        /* independent of X and Y (constant Z) */
        } else if (itype == 2) {
            data[0] = tess.xyz[3*ipnt+2];
            data[1] =  0;
            data[2] =  0;
            data[3] = -1;
        /* independent of X */
        } else if (itype == 3 && jpnt >= 0) {
            data[0] = 1;
            data[1] = 0;
            data[2] = (tess.xyz[3*ipnt+2] - tess.xyz[3*jpnt+2])
                    / (tess.xyz[3*ipnt+1] * tess.xyz[3*jpnt+2]
                     - tess.xyz[3*jpnt+1] * tess.xyz[3*ipnt+2]);
            data[3] = (tess.xyz[3*jpnt+1] - tess.xyz[3*ipnt+1])
                    / (tess.xyz[3*ipnt+1] * tess.xyz[3*jpnt+2]
                     - tess.xyz[3*jpnt+1] * tess.xyz[3*ipnt+2]);
        /* independent of Y */
        } else if (itype == 4 && jpnt >= 0) {
            data[0] = 1;
            data[1] = (tess.xyz[3*jpnt+2] - tess.xyz[3*ipnt+2])
                    / (tess.xyz[3*ipnt+2] * tess.xyz[3*jpnt  ]
                     - tess.xyz[3*jpnt+2] * tess.xyz[3*ipnt  ]);
            data[2] = 0;
            data[3] = (tess.xyz[3*ipnt  ] - tess.xyz[3*jpnt  ])
                    / (tess.xyz[3*ipnt+2] * tess.xyz[3*jpnt  ]
                     - tess.xyz[3*jpnt+2] * tess.xyz[3*ipnt  ]);
        /* independent of Z */
        } else if (itype == 5 && jpnt >= 0) {
            data[0] = 1;
            data[1] = (tess.xyz[3*ipnt+1] - tess.xyz[3*jpnt+1])
                    / (tess.xyz[3*ipnt  ] * tess.xyz[3*jpnt+1]
                     - tess.xyz[3*jpnt  ] * tess.xyz[3*ipnt+1]);
            data[2] = (tess.xyz[3*jpnt  ] - tess.xyz[3*ipnt  ])
                    / (tess.xyz[3*ipnt  ] * tess.xyz[3*jpnt+1]
                     - tess.xyz[3*jpnt  ] * tess.xyz[3*ipnt+1]);
            data[3] = 0;
        /* data vector given */
        } else {
            printf("Enter data: "); fflush(stdout);
            fscanf(stdin, "%lf %lf %lf %lf", &(data[0]), &(data[1]), &(data[2]), &(data[3]));
        }

        /* cut the Triangles */
        status = cutTriangles(&tess, icolr, data);
        CHECK_FOR_ERROR(cutTriangles, cutTriangles);

        /* build the response */
        snprintf(response, max_resp_len, "cutTriangles;%d;okay", tess.ncolr);
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Links_pend = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "deleteTriangle;x;y;z;" */
    } else if (strncmp(text, "deleteTriangle;", 15) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(deleteTriangle, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc  = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc  = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc  = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Triangle */
        itri = closestTriangle(xloc, yloc, zloc);
        SPRINT1(3, "closestTriangle -> itri=%d", itri);

        /* delete the Triangle */
        status = deleteTriangle(&tess, itri);
        CHECK_FOR_ERROR(deleteTriangle, deleteTriangle);

        /* reset the current Point */
        CurPt_index = -1;

        /* build the response */
        snprintf(response, max_resp_len, "deleteTriangle;okay");
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "fillHole;x;y;z;" */
    } else if (strncmp(text, "fillHole;", 9) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(fillHole, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc  = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc  = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc  = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Point */
        ipnt = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> ipnt=%d", ipnt);

        /* fill the hole adjacent to this Point */
        status = fillLoop(&tess, ipnt);
        CHECK_FOR_ERROR(fillHole, fillLoop);

        /* reset the current Point */
        CurPt_index = -1;

        /* build the response */
        snprintf(response, max_resp_len, "fillHole;okay");
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "flattenColor;x;y;z;" */
    } else if (strncmp(text, "flattenColor;", 13) == 0) {
        int    icolr;
        double xloc=0, yloc=0, zloc=0, tol=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(colorTriangles, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc  = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc  = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc  = strtod(arg3, NULL);
        if (getToken(text, 4, arg4)) tol   = strtod(arg4, NULL);
        SPRINT4(3, "xloc=%f  yloc=%f  zloc=%f  tol=%f", xloc, yloc, zloc, tol);

        /* find the closest Triangle */
        itri = closestTriangle(xloc, yloc, zloc);
        icolr = tess.ttyp[itri] & TRI_COLOR;
        SPRINT2(3, "closestTriangle -> itri=%d (color %d)", itri, icolr);

        /* color the Triangle and its neighbors */
        status = flattenColor(&tess, icolr, tol);
        CHECK_FOR_ERROR(flattenColor, flattenColor);

        /* build the response */
        snprintf(response, max_resp_len, "flattenColor;okay");
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Links_pend = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "generateEgads;filename;ncp;" */
    } else if (strncmp(text, "generateEgads;", 14) == 0) {
        int  ncp = 5;
        char message[80];

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        if (getToken(text, 1, arg1)) strncpy(filename, arg1, 80);
        if (getToken(text, 2, arg2)) ncp = strtol(arg2, NULL, 10);

        /* generate the Brep */
        status = generateBrep(message);
        CHECK_FOR_ERROR2(generateBrep, message);

        /* generate the fits */
        if (nctrlpnt > 0) {
            ncp = nctrlpnt;
            printf("WARNING:: overriding ncp=%d\n", ncp);
        }
        status = generateFits(ncp, message);
        CHECK_FOR_ERROR2(generateFits, message);

        /* generate the egads file */
        status = generateEgads(filename, message);
        CHECK_FOR_ERROR2(generateEgads, message);

        /* build the response */
        snprintf(response, max_resp_len, "generateEgads;%s", message);
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "joinPoints;x;y;z;" */
    } else if (strncmp(text, "joinPoints;", 11) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(joinPoints, storeUndo);

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Point */
        ipnt = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> ipnt=%d", ipnt);

        /* join the Points */
        status = joinPoints(&tess, ipnt, CurPt_index);
        CHECK_FOR_ERROR(joinPoints, joinPoints);

        /* build the response */
        snprintf(response, max_resp_len, "joinPoints;okay");
        response_len = strlen(response);

        /* set flags */
        CurPt_index = ipnt;

        Tris_pend  = 1;
        Hangs_pend = 1;
        Links_pend = 1;
        CurPt_pend = 1;

    /* "linkToPoint;x;y;z;" */
    } else if (strncmp(text, "linkToPoint;", 11) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        if (CurPt_index > 0) {
            status = storeUndo();
            CHECK_FOR_ERROR(linkToPoint, storeUndo);
        }

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Point */
        ipnt = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> ipnt=%d", ipnt);

        /* build the Links */
        if (CurPt_index >= 0) {
            status = createLinks(&tess, CurPt_index, ipnt);
            CHECK_FOR_ERROR(linkToPoint, createLinks);

            /* build the response */
            snprintf(response, max_resp_len, "linkToPoint;%d;okay", CurPt_index);
            response_len = strlen(response);
        } else {
            CurPt_index = ipnt;

            snprintf(response, max_resp_len, "pickPoint;%d;okay", CurPt_index);
            response_len = strlen(response);
        }

        /* set flags */
        CurPt_index = ipnt;

        Tris_pend  = 1;
        Links_pend = 1;
        CurPt_pend = 1;

    /* "markCreases;angdeg;" */
    } else if (strncmp(text, "markCreases;", 12) == 0) {
        double angdeg = 45;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        status = storeUndo();
        CHECK_FOR_ERROR(colorTriangles, storeUndo);

        /* extract argument */
        if (getToken(text, 1, arg1)) angdeg = strtod(arg1, NULL);

        /* mark the creases */
        status = detectCreases(&tess, angdeg);
        CHECK_FOR_ERROR(detectCreases, detectCreases);

        /* build the response */
        snprintf(response, max_resp_len, "markCreases;okay");
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Links_pend = 1;
        CurPt_pend = 1;

    /* "identifyPoint;x;y;z" */
    } else if (strncmp(text, "identifyPoint;", 14) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Point */
        ipnt = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> CurPt_index=%d", ipnt);

        /* print out information */
        printf("ipnt=%6d\n", ipnt);
        for (itri = 0; itri < tess.ntri; itri++) {
            if (tess.trip[3*itri  ] == ipnt ||
                tess.trip[3*itri+1] == ipnt ||
                tess.trip[3*itri+2] == ipnt   ) {
                printf("     itri=%6d: points= %6d %6d %6d, tris= %6d (%3d) %6d (%3d) %6d (%3d)\n", itri,
                       tess.trip[3*itri], tess.trip[3*itri+1], tess.trip[3*itri+2],
                       tess.trit[3*itri  ], tess.ttyp[tess.trit[3*itri  ]]&TRI_COLOR,
                       tess.trit[3*itri+1], tess.ttyp[tess.trit[3*itri+1]]&TRI_COLOR,
                       tess.trit[3*itri+2], tess.ttyp[tess.trit[3*itri+2]]&TRI_COLOR);
            }
        }

        /* build the response */
        snprintf(response, max_resp_len, "identifyPoint;%d;okay", CurPt_index);
        response_len = strlen(response);

        /* set flags */

    /* "pickPoint;x;y;z;" */
    } else if (strncmp(text, "pickPoint;", 10) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Point */
        CurPt_index = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> CurPt_index=%d", CurPt_index);

        /* build the response */
        snprintf(response, max_resp_len, "pickPoint;%d;okay", CurPt_index);
        response_len = strlen(response);

        /* set flags */
        CurPt_pend = 1;

    /* "scribeToPoint;x;y;z;" */
    } else if (strncmp(text, "scribeToPoint;", 14) == 0) {
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* make an undo copy */
        if (CurPt_index > 0) {
            status = storeUndo();
            CHECK_FOR_ERROR(linkToPoint, storeUndo);
        }

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Point */
        ipnt = closestPoint(xloc, yloc, zloc);
        SPRINT1(3, "closestPoint -> ipnt=%d", ipnt);

        /* build the Links */
        if (CurPt_index >= 0) {
            status = scribe(&tess, CurPt_index, ipnt);
            CHECK_FOR_ERROR(scribePoints, scribePoints);

            status = createLinks(&tess, CurPt_index, ipnt);
            CHECK_FOR_ERROR(linkToPoint, createLinks);

            /* build the response */
            snprintf(response, max_resp_len, "scribeToPoint;%d;okay", CurPt_index);
            response_len = strlen(response);
        } else {
            CurPt_index = ipnt;

            snprintf(response, max_resp_len, "scribeToPoint;%d;okay", CurPt_index);
            response_len = strlen(response);
        }

        /* set flags */
        CurPt_index = ipnt;

        Tris_pend  = 1;
        Links_pend = 1;
        CurPt_pend = 1;

    /* "identifyTriangle;x;y;z;" */
    } else if (strncmp(text, "identifyTriangle;", 17) == 0) {
        int    icolr;
        double xloc=0, yloc=0, zloc=0;

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract arguments */
        if (getToken(text, 1, arg1)) xloc  = strtod(arg1, NULL);
        if (getToken(text, 2, arg2)) yloc  = strtod(arg2, NULL);
        if (getToken(text, 3, arg3)) zloc  = strtod(arg3, NULL);
        SPRINT3(3, "xloc=%f  yloc=%f  zloc=%f", xloc, yloc, zloc);

        /* find the closest Triangle */
        itri = closestTriangle(xloc, yloc, zloc);
        SPRINT1(3, "closestTriangle -> itri=%d", itri);

        /* get its color */
        icolr = tess.ttyp[itri] & TRI_COLOR;

        /* build the response */
        snprintf(response, max_resp_len, "identifyTriangle;%d;%d;okay", itri, icolr);
        response_len = strlen(response);

        /* set flags */
        CurPt_index = -1;

        Tris_pend  = 1;
        Links_pend = 1;
        Hangs_pend = 1;
        CurPt_pend = 1;

    /* "undo;" */
    } else if (strncmp(text, "undo;", 5) == 0) {
        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* if an undo exists, use it now */
        SPRINT1(3, "tess_undo.ntri=%d", tess_undo.ntri);
        if (tess_undo.ntri > 0) {
            status = freeTess(&tess);
            CHECK_FOR_ERROR(undo, freeTess);

            status = copyTess(&tess_undo, &tess);
            CHECK_FOR_ERROR(undo, copyTess);

            status = freeTess(&tess_undo);
            CHECK_FOR_ERROR(undo, freeTess);

            /* build the response */
            snprintf(response, max_resp_len, "undo;okay");
            response_len = strlen(response);

            /* set flags */
            CurPt_index = -1;

            Tris_pend  = 1;
            Hangs_pend = 1;
            Links_pend = 1;
            CurPt_pend = 1;

        } else {
            /* build the response */
            snprintf(response, max_resp_len, "ERROR:: nothing to undo");
            response_len = strlen(response);
        }

    /* "writeStlFile;filename;" */
    } else if (strncmp(text, "writeStlFile;", 13) == 0) {

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "%s\n", text);
            fflush( jrnl_out);
        }

        /* extract argument */
        if (getToken(text, 1, arg1) == 0) arg1[0] = '\0';

        status = writeStlBinary(&tess, arg1);
        CHECK_FOR_ERROR(writeStlFile, writeStlBinary);

        /* build the response */
        snprintf(response, max_resp_len, "writeStlFile;okay");
        response_len = strlen(response);

    /* "addComment;comment;" */
    } else if (strncmp(text, "addComment;", 11) == 0) {

        /* extract argument */
        if (getToken(text, 1, arg1) == 0) arg1[0] = '\0';

        /* write journal entry */
        if (jrnl_out != NULL) {
            fprintf(jrnl_out, "# %s\n", arg1);
            fflush( jrnl_out);
        }
    }

cleanup:
    FREE(arg5);
    FREE(arg4);
    FREE(arg3);
    FREE(arg2);
    FREE(arg1);
}


/*
 ***********************************************************************
 *                                                                     *
 *   getToken - get a token from a string                              *
 *                                                                     *
 ***********************************************************************
 */
static int
getToken(char   *text,                  /* (in)  full text */
         int    nskip,                  /* (in)  tokens to skip */
         char   *token)                 /* (out) token */
{
    int    lentok, i, count, iskip;

    ROUTINE(getToken);

    /* --------------------------------------------------------------- */

    token[0] = '\0';
    lentok   = 0;

    /* count the number of semi-colons */
    count = 0;
    for (i = 0; i < strlen(text); i++) {
        if (text[i] == ';') {
            count++;
        }
    }

    if (count < nskip+1) return 0;

    /* skip over nskip tokens */
    i = 0;
    for (iskip = 0; iskip < nskip; iskip++) {
        while (text[i] != ';') {
            i++;
        }
        i++;
    }

    /* extract the token we are looking for */
    while (text[i] != ';') {
        token[lentok++] = text[i++];
        token[lentok  ] = '\0';
    }

    return strlen(token);
}


/*
 ***********************************************************************
 *                                                                     *
 *   closestPoint - find closest Point to given coordinates            *
 *                                                                     *
 ***********************************************************************
 */
static int
closestPoint(double xloc,
             double yloc,
             double zloc)
{
    int status = 0;

    int    ipnt;
    double dbest, dtest;

    ROUTINE(closestPoint);

    /* --------------------------------------------------------------- */

    dbest = pow(xloc - tess.xyz[0], 2)
          + pow(yloc - tess.xyz[1], 2)
          + pow(zloc - tess.xyz[2], 2);

    for (ipnt = 1; ipnt < tess.npnt; ipnt++) {
        dtest = pow(xloc - tess.xyz[3*ipnt  ], 2)
              + pow(yloc - tess.xyz[3*ipnt+1], 2)
              + pow(zloc - tess.xyz[3*ipnt+2], 2);

        if (dtest < dbest) {
            dbest  = dtest;
            status = ipnt;
        }
    }

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   closestTriangle - find closest Triangle to given coordinates      *
 *                                                                     *
 ***********************************************************************
 */
static int
closestTriangle(double xloc,
                double yloc,
                double zloc)
{
    int status = 0;

    int    itri, ip0, ip1, ip2;
    double xcent, ycent, zcent, dbest, dtest;

    ROUTINE(closestTriangle);

    /* --------------------------------------------------------------- */

    ip0 = tess.trip[0];
    ip1 = tess.trip[1];
    ip2 = tess.trip[2];

    xcent = (tess.xyz[3*ip0  ] + tess.xyz[3*ip1  ] + tess.xyz[3*ip2  ]) / 3;
    ycent = (tess.xyz[3*ip0+1] + tess.xyz[3*ip1+1] + tess.xyz[3*ip2+1]) / 3;
    zcent = (tess.xyz[3*ip0+2] + tess.xyz[3*ip1+2] + tess.xyz[3*ip2+2]) / 3;

    dbest = pow(xloc - xcent, 2)
          + pow(yloc - ycent, 2)
          + pow(zloc - zcent, 2);

    for (itri = 0; itri < tess.ntri; itri++) {
        if ((tess.ttyp[itri] & TRI_ACTIVE) == 0) continue;

        ip0 = tess.trip[3*itri  ];
        ip1 = tess.trip[3*itri+1];
        ip2 = tess.trip[3*itri+2];

        xcent = (tess.xyz[3*ip0  ] + tess.xyz[3*ip1  ] + tess.xyz[3*ip2  ]) / 3;
        ycent = (tess.xyz[3*ip0+1] + tess.xyz[3*ip1+1] + tess.xyz[3*ip2+1]) / 3;
        zcent = (tess.xyz[3*ip0+2] + tess.xyz[3*ip1+2] + tess.xyz[3*ip2+2]) / 3;

        dtest = pow(xloc - xcent, 2)
              + pow(yloc - ycent, 2)
              + pow(zloc - zcent, 2);

        if (dtest < dbest) {
            dbest  = dtest;
            status = itri;
        }
    }

    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   storeUndo - store an undo snapshot                                *
 *                                                                     *
 ***********************************************************************
 */
static int
storeUndo()
{
    int       status = SUCCESS;

    ROUTINE(storeUndo);

    /* --------------------------------------------------------------- */

    /* if there is previous undo information, free it now */
    if (tess_undo.ntri > 0) {
        status = freeTess(&tess_undo);
        if (status != SUCCESS) goto cleanup;
    }

    /* make a copy */
    status = copyTess(&tess, &tess_undo);

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 *   setColor - convert hex-code into a color                          *
 *                                                                     *
 ***********************************************************************
 */
static void
setColor(int    rgb,
         float  color[])
{
    ROUTINE(setColor);

    /* --------------------------------------------------------------- */

    color[0] = (float)((rgb & 0xff0000) / 0x10000) / 255.0;
    color[1] = (float)((rgb & 0x00ffff) / 0x00100) / 255.0;
    color[2] = (float)((rgb & 0x0000ff)          ) / 255.0;
}


#ifdef FOO
/*
 ***********************************************************************
 *                                                                     *
 *   spec_col - return color for a given scalar value                  *
 *                                                                     *
 ***********************************************************************
 */
static void
spec_col(float  scalar,
         float  lims[],
         float  color[])
{
    int   indx;
    float frac;

    ROUTINE(spec_col);

    /* --------------------------------------------------------------- */

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
        frac  = 255.0 * (scalar - lims[0]) / (lims[1] - lims[0]);
        if (frac < 0  ) frac = 0;
        if (frac > 255) frac = 255;
        indx  = frac;
        frac -= indx;
        if (indx == 255) {
            indx--;
            frac += 1.0;
        }

        color[0] = frac * color_map[3*(indx+1)  ] + (1.0-frac) * color_map[3*indx  ];
        color[1] = frac * color_map[3*(indx+1)+1] + (1.0-frac) * color_map[3*indx+1];
        color[2] = frac * color_map[3*(indx-1)+2] + (1.0-frac) * color_map[3*indx+2];
    }
}
#endif


/*
 ************************************************************************
 *                                                                      *
 *   printEgo - print EGADS topology for 5 levels                       *
 *                                                                      *
 ************************************************************************
 */
void
printEgo(ego    obj)                    /* (in)  ego to start */
{
    ego    eref0,    eref1,    eref2,    eref3,    eref4,   eref5;
    int    oclass0,  oclass1,  oclass2,  oclass3,  oclass4, oclass5;
    int    mtype0,   mtype1,   mtype2,   mtype3,   mtype4,  mtype5;
    int    nchild0,  nchild1,  nchild2,  nchild3,  nchild4, nchild5;
    int    ichild0,  ichild1,  ichild2,  ichild3,  ichild4, i;
    double data0[4], data1[4], data2[4], data3[4], data4[4], data5[4];
    ego    *ebodys0, *ebodys1, *ebodys2, *ebodys3, *ebodys4, *ebodys5;
    int    *senses0, *senses1, *senses2, *senses3, *senses4, *senses5;
    int    status, count=0;
    ego    topref, prev, next, refobj, context, prev1, next1;

    /* names take from egadsTypes.h */
    char *classname[27] = {"contxt",       "transform",  "tessellation",  "nil",
                           "empty",        "reference",  "ERROR 6",       "ERROR 7",
                           "ERROR 8",      "ERROR 9",    "pcurve",        "curve",
                           "surface",      "ERROR 13",   "ERROR 14",      "ERROR 15",
                           "ERROR 16",     "ERROR 17",   "ERROR 18",      "ERROR 19",
                           "node",         "edge",       "loop",          "face",
                           "shell",        "body",       "model"};

    /* pcurves and curves */
    char *mtypename1[10] = {"ERROR 0",     "line",       "circle",        "ellipse",
                            "parabola",    "hyperbola",  "trimmed",       "bezier",
                            "bspline",     "offset"};

    /* surfaces */
    char *mtypename2[12] = {"ERROR 0",     "plane",      "spherical",     "cylindrical",
                            "revolution",  "toroidal",   "trimmed",       "bezier",
                            "bspline",     "offset",    "conical",       "extrusion"};

    /* faces */
    char *mtypename3[ 3] = {"sreverse",    "nomtype",    "sforward"};

    /* other topology */
    char *mtypename4[10] = {"nomtype",     "onenode",    "twonode",       "open",
                            "closed",      "degenerate", "wirebody",      "facebody",
                            "sheetbody",   "solidbody"};

    char mtypename[32];

    /* --------------------------------------------------------------- */

    if (obj == NULL) {
        SPRINT0(0, "NULL");
        return;
    }

    status = EG_getContext(obj, &context);
    if (status < 0) {
        SPRINT1(0, "EG_getContext -> status=%d", status);
    }

    status = EG_getTopology(obj, &eref0, &oclass0, &mtype0,
                            data0, &nchild0, &ebodys0, &senses0);
    if (status != EGADS_SUCCESS) {
        EG_getInfo(obj, &oclass0, &mtype0, &topref, &prev, &next);
        nchild0 = -1;
        eref0   = topref;
    }

    if (oclass0 != NODE) {
        if (oclass0 == PCURVE || oclass0 == CURVE) {
            strcpy(mtypename, mtypename1[mtype0]);
        } else if (oclass0 == SURFACE) {
            strcpy(mtypename, mtypename2[mtype0]);
        } else if (oclass0 == FACE) {
            strcpy(mtypename, mtypename3[mtype0+1]);
        } else if (oclass0 == EDGE  || oclass0 == LOOP ||
                   oclass0 == SHELL || oclass0 == BODY   ) {
            strcpy(mtypename, mtypename4[mtype0]);
        } else {
            strcpy(mtypename, "");
        }
        SPRINT6(0, "oclass0=%3d (%s)  mtype0=%3d (%s)  obj=%llx,  eref0=%llx",
                oclass0, classname[oclass0], mtype0, mtypename,
                (long long)obj, (long long)eref0);
        if (oclass0 == LOOP || oclass0 == FACE) {
            SPRINT0x(0, "< senses=");
            for (i = 0; i < nchild0; i++) {
                SPRINT1x(0, "%2d ", senses0[i]);
            }
            SPRINT0(0, " ");
        }
    } else {
        SPRINT7(0,
                "oclass0=%3d (%s)  mtype0=%3d,  obj=%llx,  data0=%20.10e %20.10e %20.10e",
                oclass0, classname[oclass0], mtype0, (long long)obj, data0[0],
                data0[1], data0[2]);
    }

    next = obj->tref;
    prev = NULL;
    while (next != NULL) {
        refobj = next->attrs;
        if (refobj != context) {
            count++;
            EG_getInfo(refobj, &oclass1, &mtype1, &topref, &prev1, &next1);
            SPRINT5(0, "< refcount=%2d, refobj=%llx  (oclass=%2d (%s) mtype=%2d)",
                    count, (long long)refobj, oclass1, classname[oclass1], mtype1);
        }
        prev = next;
        next = prev->blind;
    }

    for (ichild0 = 0; ichild0 < nchild0; ichild0++) {
        status = EG_getTopology(ebodys0[ichild0], &eref1, &oclass1, &mtype1,
                                data1, &nchild1, &ebodys1, &senses1);
        if (status != EGADS_SUCCESS) {
            EG_getInfo(obj, &oclass1, &mtype1, &topref, &prev, &next);
            nchild1 = -1;
            eref1   = topref;
        }

        if (oclass1 != NODE) {
            if (oclass1 == PCURVE || oclass1 == CURVE) {
                strcpy(mtypename, mtypename1[mtype1]);
            } else if (oclass1 == SURFACE) {
                strcpy(mtypename, mtypename2[mtype1]);
            } else if (oclass1 == FACE) {
                strcpy(mtypename, mtypename3[mtype1+1]);
            } else if (oclass1 == EDGE  || oclass1 == LOOP ||
                       oclass1 == SHELL || oclass1 == BODY   ) {
                strcpy(mtypename, mtypename4[mtype1]);
            } else {
                strcpy(mtypename, "");
            }
            SPRINT6(0, ". oclass1=%3d (%s)  mtype1=%3d (%s)  obj=%llx,  eref1=%llx",
                    oclass1, classname[oclass1], mtype1, mtypename,
                    (long long)(ebodys0[ichild0]), (long long)eref1);
            if (oclass1 == LOOP || oclass1 == FACE) {
                SPRINT0x(0, ". < senses=");
                for (i = 0; i < nchild1; i++) {
                    SPRINT1x(0, "%2d ", senses1[i]);
                }
                SPRINT0(0, " ");
            }
        } else {
            SPRINT7(0, ". oclass1=%3d (%s)  mtype1=%3d,  obj=%llx,  data1=%20.10e %20.10e %20.10e",
                    oclass1, classname[oclass1], mtype1, (long long)(ebodys0[ichild0]),
                    data1[0], data1[1], data1[2]);
        }

        for (ichild1 = 0; ichild1 < nchild1; ichild1++) {
            status = EG_getTopology(ebodys1[ichild1], &eref2, &oclass2, &mtype2,
                                    data2, &nchild2, &ebodys2, &senses2);
            if (status != EGADS_SUCCESS) {
                EG_getInfo(obj, &oclass2, &mtype2, &topref, &prev, &next);
                nchild2 = -1;
                eref2   = topref;
            }

            if (oclass2 != NODE) {
                if (oclass2 == PCURVE || oclass2 == CURVE) {
                    strcpy(mtypename, mtypename1[mtype2]);
                } else if (oclass2 == SURFACE) {
                    strcpy(mtypename, mtypename2[mtype2]);
                } else if (oclass2 == FACE) {
                    strcpy(mtypename, mtypename3[mtype2+1]);
                } else if (oclass2 == EDGE  || oclass2 == LOOP ||
                           oclass2 == SHELL || oclass2 == BODY   ) {
                    strcpy(mtypename, mtypename4[mtype2]);
                } else {
                    strcpy(mtypename, "");
                }
                SPRINT6(0, ". . oclass2=%3d (%s)  mtype2=%3d (%s)  obj=%llx,  eref2=%llx",
                        oclass2, classname[oclass2], mtype2, mtypename,
                        (long long)(ebodys1[ichild1]), (long long)eref2);
                if (oclass2 == LOOP || oclass2 == FACE) {
                    SPRINT0x(0, ". . < senses=");
                    for (i=0; i < nchild2; i++) {
                        SPRINT1x(0, "%2d ", senses2[i]);
                    }
                    SPRINT0(0, " ");
                }
            } else {
                SPRINT7(0, ". . oclass2=%3d (%s)  mtype2=%3d,  obj=%llx,  data2=%20.10e %20.10e %20.10e",
                        oclass2, classname[oclass2], mtype2,
                        (long long)(ebodys1[ichild1]), data2[0], data2[1], data2[2]);
            }

            for (ichild2 = 0; ichild2 < nchild2; ichild2++) {
                status = EG_getTopology(ebodys2[ichild2], &eref3, &oclass3, &mtype3,
                               data3, &nchild3, &ebodys3, &senses3);
                if (status != EGADS_SUCCESS) {
                    EG_getInfo(obj, &oclass3, &mtype3, &topref, &prev, &next);
                    nchild3 = -1;
                    eref3   = topref;
                }

                if (oclass3 != NODE) {
                    if (oclass3 == PCURVE || oclass3 == CURVE) {
                        strcpy(mtypename, mtypename1[mtype3]);
                    } else if (oclass3 == SURFACE) {
                        strcpy(mtypename, mtypename2[mtype3]);
                    } else if (oclass3 == FACE) {
                        strcpy(mtypename, mtypename3[mtype3+1]);
                    } else if (oclass3 == EDGE  || oclass3 == LOOP ||
                               oclass3 == SHELL || oclass3 == BODY   ) {
                        strcpy(mtypename, mtypename4[mtype3]);
                    } else {
                        strcpy(mtypename, "");
                    }
                    SPRINT6(0, ". . . oclass3=%3d (%s)  mtype3=%3d (%s)  obj=%llx,  eref3=%llx",
                            oclass3, classname[oclass3], mtype3, mtypename,
                            (long long)(ebodys2[ichild2]), (long long)eref3);
                    if (oclass3 == LOOP || oclass3 == FACE) {
                        SPRINT0x(0, ". . . < senses=");
                        for (i=0; i < nchild3; i++) {
                            SPRINT1x(0, "%2d ", senses3[i]);
                        }
                        SPRINT0(0, " ");
                    }
                } else {
                    SPRINT7(0, ". . . oclass3=%3d (%s)  mtype3=%3d,  obj=%llx,  data3=%20.10e %20.10e %20.10e",
                            oclass3, classname[oclass3], mtype3,
                            (long long)(ebodys2[ichild2]), data3[0], data3[1], data3[2]);
                }

                for (ichild3 = 0; ichild3 < nchild3; ichild3++) {
                    status = EG_getTopology(ebodys3[ichild3], &eref4, &oclass4, &mtype4,
                                   data4, &nchild4, &ebodys4, &senses4);
                    if (status != EGADS_SUCCESS) {
                        EG_getInfo(obj, &oclass4, &mtype4, &topref, &prev, &next);
                        nchild4 = -1;
                        eref4   = topref;
                    }

                    if (oclass4 != NODE) {
                        if (oclass4 == PCURVE || oclass4 == CURVE) {
                            strcpy(mtypename, mtypename1[mtype4]);
                        } else if (oclass4 == SURFACE) {
                            strcpy(mtypename, mtypename2[mtype4]);
                        } else if (oclass4 == FACE) {
                            strcpy(mtypename, mtypename3[mtype4+1]);
                        } else if (oclass4 == EDGE  || oclass4 == LOOP ||
                                   oclass4 == SHELL || oclass4 == BODY   ) {
                            strcpy(mtypename, mtypename4[mtype4]);
                        } else {
                            strcpy(mtypename, "");
                        }
                        SPRINT6(0, ". . . . oclass4=%3d (%s)  mtype4=%3d (%s)  obj=%llx,  eref4=%llx",
                                oclass4, classname[oclass4], mtype4, mtypename,
                                (long long)(ebodys3[ichild3]), (long long)eref4);
                        if (oclass4 == LOOP || oclass4 == FACE) {
                            SPRINT0x(0, ". . . . < senses=");
                            for (i=0; i < nchild4; i++) {
                                SPRINT1x(0, "%2d ", senses4[i]);
                            }
                            SPRINT0(0, " ");
                        }
                    } else {
                        SPRINT7(0, ". . . . oclass4=%3d (%s)  mtype4=%3d,  obj=%llx,  data4=%20.10e %20.10e %20.10e",
                                oclass4, classname[oclass4], mtype4,
                                (long long)(ebodys3[ichild3]), data4[0], data4[1], data4[2]);
                    }

                    for (ichild4 = 0; ichild4 < nchild4; ichild4++) {
                        status = EG_getTopology(ebodys4[ichild4], &eref5, &oclass5, &mtype5,
                                                data5, &nchild5, &ebodys5, &senses5);
                        if (status != EGADS_SUCCESS) {
                            EG_getInfo(obj, &oclass5, &mtype5, &topref, &prev, &next);
                            nchild5 = -1;
                            eref5   = topref;
                        }

                        if (oclass5 != NODE) {
                            if (oclass5 == PCURVE || oclass5 == CURVE) {
                                strcpy(mtypename, mtypename1[mtype5]);
                            } else if (oclass5 == SURFACE) {
                                strcpy(mtypename, mtypename2[mtype5]);
                            } else if (oclass5 == FACE) {
                                strcpy(mtypename, mtypename3[mtype5+1]);
                            } else if (oclass5 == EDGE  || oclass5 == LOOP ||
                                       oclass5 == SHELL || oclass5 == BODY   ) {
                                strcpy(mtypename, mtypename4[mtype5]);
                            } else {
                                strcpy(mtypename, "");
                            }
                            SPRINT6(0, ". . . . . oclass5=%3d (%s)  mtype5=%3d (%s)  obj=%llx,  eref5=%llx",
                                    oclass5, classname[oclass5], mtype5, mtypename,
                                    (long long)(ebodys4[ichild4]), (long long)eref5);
                            if (oclass5 == LOOP || oclass5 == FACE) {
                                SPRINT0x(0, ". . . . . < senses=");
                                for (i=0; i < nchild5; i++) {
                                    SPRINT1x(0, "%2d ", senses5[i]);
                                }
                                SPRINT0(0, " ");
                            }
                        } else {
                            SPRINT7(0, ". . . . . oclass5=%3d (%s)  mtype5=%3d,  obj=%llx,  data5=%20.10e %20.10e %20.10e",
                                    oclass5, classname[oclass5], mtype5,
                                    (long long)(ebodys4[ichild4]), data5[0], data5[1], data5[2]);
                        }
                    }
                }
            }
        }
    }
}



//$$$/*
//$$$ ************************************************************************
//$$$ *                                                                      *
//$$$ *   fitPlane - least-square fit a*x + b*y + c*z + d to  points         *
//$$$ *                                                                      *
//$$$ ************************************************************************
//$$$ */
//$$$
//$$$static int
//$$$fitPlane(double xyz[],                  /* (in)  array  of points */
//$$$         int    n,                      /* (in)  number of points */
//$$$         double *a,                     /* (out) x-coefficient */
//$$$         double *b,                     /* (out) y-coefficient */
//$$$         double *c,                     /* (out) z-coefficient */
//$$$         double *d)                     /* (out) constant */
//$$${
//$$$    int       status = SUCCESS;         /* (out) return status */
//$$$
//$$$    int       i;
//$$$    double    xcent, ycent, zcent, xdet, ydet, zdet;
//$$$    double    xx, xy, xz, yy, yz, zz;
//$$$
//$$$    ROUTINE(fitPlane);
//$$$
//$$$    /* --------------------------------------------------------------- */
//$$$
//$$$    /* default returns */
//$$$    *a = 1;
//$$$    *b = 1;
//$$$    *c = 1;
//$$$    *d = 0;
//$$$
//$$$    /* make sure we have enough points */
//$$$    if (n < 3) {
//$$$        printf("not enough points\n");
//$$$        status = -999;
//$$$        goto cleanup;
//$$$    }
//$$$
//$$$    /* find the centroid of the points */
//$$$    xcent = 0;
//$$$    ycent = 0;
//$$$    zcent = 0;
//$$$
//$$$    for (i = 0; i < n; i++) {
//$$$        xcent += xyz[3*i  ];
//$$$        ycent += xyz[3*i+1];
//$$$        zcent += xyz[3*i+2];
//$$$    }
//$$$
//$$$    xcent /= n;
//$$$    ycent /= n;
//$$$    zcent /= n;
//$$$
//$$$    /* compute the covarience matrix (relative to the controid) */
//$$$    xx = 0;
//$$$    xy = 0;
//$$$    xz = 0;
//$$$    yy = 0;
//$$$    yz = 0;
//$$$    zz = 0;
//$$$
//$$$    for (i = 0; i < n; i++) {
//$$$        xx += (xyz[3*i  ] - xcent) * (xyz[3*i  ] - xcent);
//$$$        xy += (xyz[3*i  ] - xcent) * (xyz[3*i+1] - ycent);
//$$$        xz += (xyz[3*i  ] - xcent) * (xyz[3*i+2] - zcent);
//$$$        yy += (xyz[3*i+1] - ycent) * (xyz[3*i+1] - ycent);
//$$$        yz += (xyz[3*i+1] - ycent) * (xyz[3*i+2] - zcent);
//$$$        zz += (xyz[3*i+2] - zcent) * (xyz[3*i+2] - zcent);
//$$$    }
//$$$
//$$$    /* find the deteminants associated with assuming normx=1,
//$$$       normy=1, and normz=1 */
//$$$    xdet = yy * zz - yz * yz;
//$$$    ydet = zz * xx - xz * xz;
//$$$    zdet = xx * yy - xy * xy;
//$$$
//$$$    /* use the determinant with the best conditioning */
//$$$    if        (fabs(xdet) >= fabs(ydet) && fabs(xdet) >= fabs(zdet)) {
//$$$        *a = 1;
//$$$        *b = (xz * yz - xy * zz) / xdet;
//$$$        *c = (xy * yz - xz * yy) / xdet;
//$$$        *d = -xcent              / xdet;
//$$$    } else if (fabs(ydet) >= fabs(zdet) && fabs(ydet) >= fabs(xdet)) {
//$$$        *a = (yz * xz - xy * zz) / ydet;
//$$$        *b = 1;
//$$$        *c = (xy * yz - yz * xx) / ydet;
//$$$        *d = -ycent              / ydet;
//$$$    } else if (fabs(zdet) >= EPS12) {
//$$$        *a = (yz * xy - xz * yy) / zdet;
//$$$        *b = (xz * xy - yz * xx) / zdet;
//$$$        *c = 1;
//$$$        *d = -zcent              / zdet;
//$$$    } else {
//$$$        printf("cannot find big determinant (%f, %f, %f)\n", xdet, ydet, zdet);
//$$$        status = -999;
//$$$        goto cleanup;
//$$$    }
//$$$
//$$$cleanup:
//$$$    return status;
//$$$}
//$$$
//$$$//$$$ http://www.ilikebigbits.com/blog/2015/3/2/plane-from-points
//$$$//$$$// Constructs a plane from a collection of points
//$$$//$$$// so that the summed squared distance to all points is minimzized
//$$$//$$$fn plane_from_points(points: &[Vec3]) -> Plane {
//$$$//$$$    let n = points.len();
//$$$//$$$    assert!(n >= 3, "At least three points required");
//$$$//$$$
//$$$//$$$    let mut sum = Vec3{x:0.0, y:0.0, z:0.0};
//$$$//$$$    for p in points {
//$$$//$$$        sum = &sum + &p;
//$$$//$$$    }
//$$$//$$$    let centroid = &sum * (1.0 / (n as f64));
//$$$//$$$
//$$$//$$$    // Calc full 3x3 covariance matrix, excluding symmetries:
//$$$//$$$    let mut xx = 0.0; let mut xy = 0.0; let mut xz = 0.0;
//$$$//$$$    let mut yy = 0.0; let mut yz = 0.0; let mut zz = 0.0;
//$$$//$$$
//$$$//$$$    for p in points {
//$$$//$$$        let r = p - &centroid;
//$$$//$$$        xx += r.x * r.x;
//$$$//$$$        xy += r.x * r.y;
//$$$//$$$        xz += r.x * r.z;
//$$$//$$$        yy += r.y * r.y;
//$$$//$$$        yz += r.y * r.z;
//$$$//$$$        zz += r.z * r.z;
//$$$//$$$    }
//$$$//$$$
//$$$//$$$    let det_x = yy*zz - yz*yz;
//$$$//$$$    let det_y = xx*zz - xz*xz;
//$$$//$$$    let det_z = xx*yy - xy*xy;
//$$$//$$$
//$$$//$$$    let det_max = max3(det_x, det_y, det_z);
//$$$//$$$    assert!(det_max > 0.0, "The points don't span a plane");
//$$$//$$$
//$$$//$$$    // Pick path with best conditioning:
//$$$//$$$    let dir =
//$$$//$$$        if det_max == det_x {
//$$$//$$$            let a = (xz*yz - xy*zz) / det_x;
//$$$//$$$            let b = (xy*yz - xz*yy) / det_x;
//$$$//$$$            Vec3{x: 1.0, y: a, z: b}
//$$$//$$$        } else if det_max == det_y {
//$$$//$$$            let a = (yz*xz - xy*zz) / det_y;
//$$$//$$$            let b = (xy*xz - yz*xx) / det_y;
//$$$//$$$            Vec3{x: a, y: 1.0, z: b}
//$$$//$$$        } else {
//$$$//$$$            let a = (yz*xy - xz*yy) / det_z;
//$$$//$$$            let b = (xz*xy - yz*xx) / det_z;
//$$$//$$$            Vec3{x: a, y: b, z: 1.0}
//$$$//$$$        };
//$$$//$$$
//$$$//$$$    plane_from_point_and_normal(&centroid, &normalize(dir))
//$$$//$$$}


/*
 ************************************************************************
 *                                                                      *
 *   plotCurve_image - plot curve data (level 3)                        *
 *                                                                      *
 ************************************************************************
 */
#ifdef GRAFIC
static void
plotPoints_image(int   *ifunct,
                 void  *itypeP,
                 void  *ntrainP,
                 void  *xyztrainP,
                 void  *npntP,
                 void  *xyzP,
     /*@unused@*/void  *a5,
     /*@unused@*/void  *a6,
     /*@unused@*/void  *a7,
     /*@unused@*/void  *a8,
     /*@unused@*/void  *a9,
                float *scale,
                char  *text,
                int   textlen)
{
    int    *itype    = (int    *)itypeP;
    int    *ntrain   = (int    *)ntrainP;
    double *xyztrain = (double *)xyztrainP;
    int    *npnt     = (int    *)npntP;
    double *xyz      = (double *)xyzP;

    int    k;
    float  x4, y4, z4;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    int    icircle = GR_CIRCLE;
//$$$    int    istar   = GR_STAR;
    int    iplus   = GR_PLUS;
//$$$    int    isolid  = GR_SOLID;
//$$$    int    idotted = GR_DOTTED;
//$$$    int    igreen  = GR_GREEN;
    int    iblue   = GR_BLUE;
    int    ired    = GR_RED;
    int    iblack  = GR_BLACK;

    ROUTINE(plot1dBspine_image);

    /* --------------------------------------------------------------- */

    if (*ifunct == 0) {
        xmin = xyz[0];
        xmax = xyz[0];
        ymin = xyz[1];
        ymax = xyz[1];
        zmin = xyz[2];
        zmax = xyz[2];

        for (k = 1; k < *npnt; k++) {
            if (xyz[3*k  ] < xmin) xmin = xyz[3*k  ];
            if (xyz[3*k  ] > xmax) xmax = xyz[3*k  ];
            if (xyz[3*k+1] < ymin) ymin = xyz[3*k+1];
            if (xyz[3*k+1] > ymax) ymax = xyz[3*k+1];
            if (xyz[3*k+2] < zmin) zmin = xyz[3*k+2];
            if (xyz[3*k+2] > zmax) zmax = xyz[3*k+2];
        }

        if        (xmax-xmin >= zmax-zmin && ymax-ymin >= zmax-zmin) {
            *itype = 0;
            scale[0] = xmin - EPS06;
            scale[1] = xmax + EPS06;
            scale[2] = ymin - EPS06;
            scale[3] = ymax + EPS06;
        } else if (ymax-ymin >= xmax-xmin && zmax-zmin >= xmax-xmin) {
            *itype = 1;
            scale[0] = ymin - EPS06;
            scale[1] = ymax + EPS06;
            scale[2] = zmin - EPS06;
            scale[3] = zmax + EPS06;
        } else {
            *itype = 2;
            scale[0] = zmin - EPS06;
            scale[1] = zmax + EPS06;
            scale[2] = xmin - EPS06;
            scale[3] = xmax + EPS06;
        }

        strcpy(text, " ");

    } else if (*ifunct == 1) {

        /* training points */
        for (k = 0; k < *ntrain; k++) {
            x4 = xyztrain[3*k  ];
            y4 = xyztrain[3*k+1];
            z4 = xyztrain[3*k+2];

            grmov3_(&x4, &y4, &z4);
            grcolr_(&ired);
            grsymb_(&icircle);
        }



        /* cloud of points */
        for (k = 0; k < *npnt; k++) {
            x4 = xyz[3*k  ];
            y4 = xyz[3*k+1];
            z4 = xyz[3*k+2];

            grmov3_(&x4, &y4, &z4);

            grcolr_(&iblue);
            grsymb_(&iplus);
        }

        grcolr_(&iblack);

    } else {
        printf("illegal option\n");
    }
}
#endif
