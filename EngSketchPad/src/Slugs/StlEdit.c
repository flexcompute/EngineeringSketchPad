/************************************************************************/
/*                                                                      */
/*   StlEdit - stl file editor                                          */
/*                                                                      */
/*   written by John Dannenhoffer @ Syracuse University                 */
/*                                                                      */
/************************************************************************/

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
#include <math.h>
#include <string.h>
#include <assert.h>

#include "gv.h"
#include "Graphics.h"
#include "common.h"
#include "Tessellate.h"
#include "RedBlackTree.h"

#ifdef GRAFIC
    #include "grafic.h"
#endif

#define MIN3(A,B,C)   (MIN(MIN(A,B),C))
#define MAX3(A,B,C)   (MAX(MAX(A,B),C))

#define UINT32 unsigned int
#define UINT16 unsigned short int
#define REAL32 float

/* window defines */
#define  DataBase        1
#define  TwoD            2
#define  ThreeD          3
#define  Dials           4
#define  Key             5

/* external routines and global variables (in gv) */
extern   GWinAtt         gv_wAux;
extern   GWinAtt         gv_wDial;
extern   float           gv_black[3];
extern   float           gv_white[3];
extern   float           gv_lincor;
extern   double          gv_located[4];

extern   GWinAtt         gv_w3d;
extern   float           gv_xform[4][4];
extern   GvGraphic**     gv_list;
extern   GvGraphic*      gv_picked;
extern   int             gv_pickmask;
extern   void            PickGraphic(float x, float y, int flag);
extern   void            LocateGraphic(float x, float y, int flag);

static   void            *realloc_temp = NULL;       /* used by RALLOC macro */

/* data structures */
typedef struct {
    int                  ipnt;               /* Point index (bias-0) */
    int                  nedg;               /* number of incident Edges */
} nod_T;

typedef struct {
    int                  ibeg;               /* Node at beginning */
    int                  iend;               /* Node at end */
    int                  ileft;              /* Face (or color) on left */
    int                  irite;              /* Face (or color) on rite */
    int                  npnt;               /* number of Points along Edge */
    int                  *pnt;               /* array  of Points along Edge */
    int                  mark;               /* 0=vmin, 1=umax, 2=vmax, 3=umin, -1=unmarked */
} edg_T;

typedef struct {
    int                  icol;               /* color  of associated Triangles */
    tess_T               tess;               /* TESS object */
    int                  nedg;               /* number of associated Edges */
    int                  *edg;               /* array  of associated Edges */

    int                  imax;               /* first  surface dimension */
    int                  jmax;               /* second surface dimension */
    double               *xsrf;              /* x coordinates for surface */
    double               *ysrf;              /* y coordinates for surface */
    double               *zsrf;              /* z coordinates for surface */
} fac_T;

/* global variables */
static   int             new_data   =  1;    /* = 0 when GvGraphic objects are up to date
                                                = 1 when data is updated */
static   int             gridOn     =  0;    /* =1 to show grids */
static   int             fly_mode   =  1;    /* =0 arrow keys rotate, =1 arrow keys translate */
static   int             ngrobj     =  0;    /* number of GvGraphic objects */
static   GvGraphic       **grobjs;           /* array  of GvGraphic objects */
static   float           *saveit    = NULL;  /* cached scalar */
static   int             markedPnt  = -1;    /* index of marked Point    (or -1) */
static   int             markedTri  = -1;    /* index of marked Triangle (or -1) */
static   int             numarg     = -1;    /* numeric argument (-1 if not set) */
static   FILE            *script    = NULL;  /* pointer to script file (or NULL) */
static   FILE            *fpdump    = NULL;  /* pointer to dump file */

static   tess_T          tess;               /* global TESS object */

static   int             Nnod  = 0;          /* number of Nodes */
static   nod_T           *nod  = NULL;       /* array  of Nodes */
static   int             Nedg  = 0;          /* number of Edges */
static   edg_T           *edg  = NULL;       /* array  of Edges */
static   int             Nfac  = 0;          /* number of Faces */
static   fac_T           *fac  = NULL;       /* array  of faces */

/* forward declaration of routines defined below */
static int                      addEdge(int ileft, int irite, int ntemp, int itemp[], int ibeg, int iend);
static int                      addFace(int icolr);
static int                      addNode(int ipnt, int pntNod[]);
static int                      getModelSize(double box[]);
static int                      findEdge(int xcsr, int yscr, int *ipnt);
static int                      findFace(int xscr, int yscr, int *itri);
static int                      findPoint(int xscr, int yscr);
static int                      findTriangle(int xscr, int yscr);
static int                      makeCurve(int iedg);
static int                      makeCut(double zcut, int ifac, double xform[3][3], int *ncut,
                                        double *xcut[], double *ycut[]);
static int                      makeSurface(int ifac, int itype, int imax, int jmax);
static int                      makeSurface1(int ifac, int imax, int jmax, double xyzs[], int senw[]);
static int                      makeSurface2(int ifac, int imax, int jmax, double xyzs[], int senw[]);
static int                      makeTopology();
static int                      projectToFace(int ifac, double xyz_in[], double xyz_out[]);
static int                      splitEdge(int iedg, int ipnt);
static int                      trim(char *s);
static int                      writeEgads(char *filename);
//     void                     gvdraw(int phase);
//     int                      gvdupdate();
//     void                     gvdata(int ngraphics, GvGraphic *graphic[]);
//     int                      gvscalar(int key, GvGraphic *graphic, int lenm float *scalar);
//     void                     gvevent(int *win, int *type, int *xscr, int *yscr, int *state);
static int                      getInt(char *prompt);
static double                   getDbl(char *prompt);
static void                     getStr(char *prompt, char *answer);

#ifdef GRAFIC
static void                     plot2D(int*, void*, void*, void*, void*, void*, void*,
                                             void*, void*, void*, void*, float*, char*, int);
//$$$static void                     plot3D(int*, void*, void*, void*, void*, void*, void*,
//$$$                                             void*, void*, void*, void*, float*, char*, int);
#endif



/**********************************************************************/
int
main(int    argc,                       /* (in)  number of arguments */
     char   *argv[])                    /* (in)  array  of arguments */
{
    int    status;
    int    mtflag = -1;             /* multi-threading flag */
    float  focus[4];
    double box[6];
    char   filename[255], test[80];
    FILE   *fp = NULL;

    int     keys[4] =   {'w',            'x',            'y',            'z'             };
    int    types[4] =   {GV_SURFFACET,   GV_SURF,        GV_SURF,        GV_SURF         };
    float   lims[8] =   {0.0, 1.0,       0.0, 1.0,       0.0, 1.0,       0.0, 1.0        };
    char   titles[65] = "Facet colors    X coordinate    Y coordinate    Z coordinate    ";

//    ROUTINE(main);

    /*----------------------------------------------------------------*/

    if        (sizeof(UINT16) != 2) {
        printf("ERR: uint16 should have size 2\a\n");
        exit(0);
    } else if (sizeof(UINT32) != 4) {
        printf("ERR: uint32 should have size 4\a\n");
        exit(0);
    } else if (sizeof(REAL32) != 4) {
        printf("ERR: real32 should have size 4\a\n");
        exit(0);
    }

    /* welcome banner */
    printf("\n");
    printf("*****************************************************\n");
    printf("*                                                   *\n");
    printf("*                  Program StlEdit                  *\n");
    printf("*                                                   *\n");
    printf("*        written by John Dannenhoffer, 2013/2022    *\n");
    printf("*                                                   *\n");
    printf("*****************************************************\n");
    printf("\n");

    /* read stl file(s) */
    if (argc < 2) {
        printf("Proper call is: StlEdit filename [journal]\a\n");
        exit(0);
    } else {
        strncpy(filename, argv[1], 254);

        if        (strstr(filename, ".tri") != NULL) {
            readTriAscii(&tess, filename);

        } else if (strstr(filename, ".stl") != NULL) {
            fp = fopen(filename, "rb");
            if (fp == NULL) {
                printf("ERR: \"%s\" does not exist\a\n", filename);
                exit(0);
            }

            (void) fread(test, sizeof(char), 5, fp);
            fclose(fp);

            if (strncmp(test, "solid", 5) == 0) {
                printf("    \"%s\" is an ASCII file\n", filename);
                readStlAscii(&tess, filename);
            } else {
                printf("    \"%s\" is a binary file\n", filename);
                readStlBinary(&tess, filename);
            }
        } else {
            printf("ERR: \"%s\" is not a .stl or .tri file\n", filename);
            exit(0);
        }
    }

    /* get size of model for sizing gv viewer */
    getModelSize(box);

    focus[0] = (float)(0.5*(box[0]+box[3]));
    focus[1] = (float)(0.5*(box[1]+box[4]));
    focus[2] = (float)(0.5*(box[2]+box[5]));
    focus[3] = (float)(sqrt( (box[0]-box[3])*(box[0]-box[3])
                            +(box[1]-box[4])*(box[1]-box[4])
                            +(box[2]-box[5])*(box[2]-box[5]) ));

    /* open the dump file */
    fpdump = fopen("StlEdit.dump", "w");

    /* automatically fire journal (if specified) */
    if (argc == 3) {
        printf("    Opening journal %s...\n", argv[2]);
        script = fopen(argv[2], "r");
    }

    /* make background white and foreground black */
    gv_black[0] = gv_black[1] = gv_black[2] = 1.0;
    gv_white[0] = gv_white[1] = gv_white[2] = 0.0;

    /* start the viewer */
    status = gv_init("StlEditor", mtflag, 4, keys, types, lims,
                     titles, focus);
    printf("gv_init -> status=%d\n", status);

    /* clean up */
    fprintf(fpdump, "$\n");
    fclose(fpdump);

//cleanup:
    return 0;
}


/******************************************************************************/

static int
addEdge(int    ileft,                   /* (in)  color on left */
        int    irite,                   /* (in)  color on rite */
        int    ntemp,                   /* (in)  number of Points along Edge */
        int    itemp[],                 /* (in)  array  of Points along Edge */
        int    ibeg,                    /* (in)  Node at beginning */
        int    iend)                    /* (in)  Node at end */
{
    int    status = SUCCESS;            /* (out) return status */

    int    i;

    ROUTINE(addEdge);

    /*----------------------------------------------------------------*/

    assert(ntemp > 0);                            // needed to avoid clang warning

    /* make room for the new Edge */
    if (Nedg == 0) {
        MALLOC(edg, edg_T,      1);
    } else {
        RALLOC(edg, edg_T, Nedg+1);
    }

    /* store Edge info */
    edg[Nedg].ibeg  = ibeg;
    edg[Nedg].iend  = iend;
    edg[Nedg].ileft = ileft;
    edg[Nedg].irite = irite;
    edg[Nedg].npnt  = ntemp;
    edg[Nedg].pnt   = NULL;
    edg[Nedg].mark  = -1;

    MALLOC(edg[Nedg].pnt, int, ntemp);

    for (i = 0; i < ntemp; i++) {
        edg[Nedg].pnt[i] = itemp[i];
    }

    /* increase valence of .ibeg and .iend */
    (nod[ibeg].nedg)++;
    (nod[iend].nedg)++;

    /* increment number of Edges */
    Nedg++;

cleanup:
    return status;
}


/******************************************************************************/

static int
addFace(int    icolr)                   /* (in)  color of Triangles in Face */
{
    int    status = SUCCESS;            /* (out) return status */

    int    itri, isid, ipnt, jpnt, iedg, ip0, ip1, ip2;
    LONG   key1, key2, key3;
    rbt_T  *ntree = NULL;

    ROUTINE(addFace);

    /*----------------------------------------------------------------*/

    /* make room for the new Face */
    if (Nfac == 0) {
        MALLOC(fac, fac_T,      1);
    } else {
        RALLOC(fac, fac_T, Nfac+1);
    }

    /* store Face info */
    fac[Nfac].icol = icolr;

    fac[Nfac].nedg = 0;
    fac[Nfac].edg  = NULL;

    fac[Nfac].imax = 0;
    fac[Nfac].jmax = 0;
    fac[Nfac].xsrf = NULL;
    fac[Nfac].ysrf = NULL;
    fac[Nfac].zsrf = NULL;

    status = initialTess(&(fac[Nfac].tess));
    CHECK_STATUS(initialTess);

    /* store the Edges associated with this Face */
    for (iedg = 0; iedg < Nedg; iedg++) {
        if (edg[iedg].ileft == icolr) {
            (fac[Nfac].nedg)++;
        }

        if (edg[iedg].irite == icolr) {
            (fac[Nfac].nedg)++;
        }
    }

    MALLOC(fac[Nfac].edg, int, fac[Nfac].nedg);

    fac[Nfac].nedg = 0;
    for (iedg = 0; iedg < Nedg; iedg++) {
        if (edg[iedg].ileft == icolr) {
            fac[Nfac].edg[fac[Nfac].nedg] = iedg;
            (fac[Nfac].nedg)++;
        }

        if (edg[iedg].irite == icolr) {
            fac[Nfac].edg[fac[Nfac].nedg] = iedg;
            (fac[Nfac].nedg)++;
        }
    }

    /* count the number of Triangles of this color */
    fac[Nfac].tess.mtri = 0;
    for (itri = 0; itri < tess.ntri; itri++) {
        if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
            (fac[Nfac].tess.mtri)++;
        }
    }

    assert(fac[Nfac].tess.mtri > 0);              // needed to avoid clang warning

    fac[Nfac].tess.mpnt = 3 * fac[Nfac].tess.mtri;  /* larger than needed */

    RALLOC(fac[Nfac].tess.trip, int,    3*fac[Nfac].tess.mtri);
    RALLOC(fac[Nfac].tess.trit, int,    3*fac[Nfac].tess.mtri);
    RALLOC(fac[Nfac].tess.ttyp, int,      fac[Nfac].tess.mtri);
    RALLOC(fac[Nfac].tess.bbox, double, 6*fac[Nfac].tess.mtri);
    RALLOC(fac[Nfac].tess.xyz,  double, 3*fac[Nfac].tess.mpnt);
    RALLOC(fac[Nfac].tess.uv,   double, 2*fac[Nfac].tess.mpnt);

    /* get a red-black tree in which the Points will be stored */
    status = rbtCreate(3*fac[Nfac].tess.mtri, &ntree);
    if (ntree == NULL) {
        printf("ERR: ntree could not allocated in makeCut\a\n");
        exit(0);
    }

    /* create the Triangles (and Points) */
    fac[Nfac].tess.ntri = 0;
    for (itri = 0; itri < tess.ntri; itri++) {
        if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
            for (isid = 0; isid < 3; isid++) {
                ipnt = tess.trip[3*itri+isid];

                /* see if the Point already exists */
                key1 = (LONG)(tess.xyz[3*ipnt  ] * 100000);
                key2 = (LONG)(tess.xyz[3*ipnt+1] * 100000);
                key3 = (LONG)(tess.xyz[3*ipnt+2] * 100000);
                jpnt = rbtSearch(ntree, key1, key2, key3);

                /* create a new Point if not found */
                if (jpnt < 0) {
                    fac[Nfac].tess.xyz[3*fac[Nfac].tess.npnt  ] = tess.xyz[3*ipnt  ];
                    fac[Nfac].tess.xyz[3*fac[Nfac].tess.npnt+1] = tess.xyz[3*ipnt+1];
                    fac[Nfac].tess.xyz[3*fac[Nfac].tess.npnt+2] = tess.xyz[3*ipnt+2];
                    (fac[Nfac].tess.npnt)++;

                    jpnt = rbtInsert(ntree, key1, key2, key3);
                    if (jpnt != (fac[Nfac].tess.npnt-1)) {
                        printf("ERR: Trouble inserting into tree, jpnt=%d, .pnt=%d\a\n", jpnt, fac[Nfac].tess.npnt);
                        exit(0);
                    }
                }

                /* remember the Point's id */
                fac[Nfac].tess.trip[3*fac[Nfac].tess.ntri+isid] = jpnt;
            }

            /* create the Triangle */
            fac[Nfac].tess.trit[3*fac[Nfac].tess.ntri  ] = -1;
            fac[Nfac].tess.trit[3*fac[Nfac].tess.ntri+1] = -1;
            fac[Nfac].tess.trit[3*fac[Nfac].tess.ntri+2] = -1;
            fac[Nfac].tess.ttyp[  fac[Nfac].tess.ntri  ] = TRI_ACTIVE | TRI_VISIBLE;

            (fac[Nfac].tess.ntri)++;
        }
    }

    /* free up the red-black tree */
    rbtDelete(ntree);

    /* reallocate arrays with smaller arrays */
    fac[Nfac].tess.mpnt = fac[Nfac].tess.npnt;

    RALLOC(fac[Nfac].tess.xyz, double, 3*fac[Nfac].tess.mpnt);
    RALLOC(fac[Nfac].tess.uv,  double, 2*fac[Nfac].tess.mpnt);

    /* set up the bounding boxes of the Triangles */
    for (itri = 0; itri < fac[Nfac].tess.ntri; itri++) {
        ip0 = fac[Nfac].tess.trip[3*itri  ];
        ip1 = fac[Nfac].tess.trip[3*itri+1];
        ip2 = fac[Nfac].tess.trip[3*itri+2];

        fac[Nfac].tess.bbox[6*itri  ] = MIN3(fac[Nfac].tess.xyz[3*ip0  ],
                                             fac[Nfac].tess.xyz[3*ip1  ],
                                             fac[Nfac].tess.xyz[3*ip2  ]);
        fac[Nfac].tess.bbox[6*itri+1] = MAX3(fac[Nfac].tess.xyz[3*ip0  ],
                                             fac[Nfac].tess.xyz[3*ip1  ],
                                             fac[Nfac].tess.xyz[3*ip2  ]);
        fac[Nfac].tess.bbox[6*itri+2] = MIN3(fac[Nfac].tess.xyz[3*ip0+1],
                                             fac[Nfac].tess.xyz[3*ip1+1],
                                             fac[Nfac].tess.xyz[3*ip2+1]);
        fac[Nfac].tess.bbox[6*itri+3] = MAX3(fac[Nfac].tess.xyz[3*ip0+1],
                                             fac[Nfac].tess.xyz[3*ip1+1],
                                             fac[Nfac].tess.xyz[3*ip2+1]);
        fac[Nfac].tess.bbox[6*itri+4] = MIN3(fac[Nfac].tess.xyz[3*ip0+2],
                                             fac[Nfac].tess.xyz[3*ip1+2],
                                             fac[Nfac].tess.xyz[3*ip2+2]);
        fac[Nfac].tess.bbox[6*itri+5] = MAX3(fac[Nfac].tess.xyz[3*ip0+2],
                                             fac[Nfac].tess.xyz[3*ip1+2],
                                             fac[Nfac].tess.xyz[3*ip2+2]);
    }

    /* set up the neighbors */
    setupNeighbors(&(fac[Nfac].tess));

    /* increment number of Faces */
    Nfac++;

cleanup:
    FREE(ntree);
    return status;
}


/******************************************************************************/

static int
addNode(int    ipnt,                    /* (in)  Point index (bias-0) */
        int    pntNod[])                /* (in)  Nodes associated with each Point */
{
    int    status = SUCCESS;            /* (out) return status */

    ROUTINE(addNode);

    /*----------------------------------------------------------------*/

    /* if a Node already exists at this Point, simply return */
    if (pntNod != NULL) {
        if (pntNod[ipnt] >= 0) {
            goto cleanup;
        }
    }

    /* otherwise, make room for new Node */
    if (Nnod == 0) {
        MALLOC(nod, nod_T,      1);
    } else {
        RALLOC(nod, nod_T, Nnod+1);
    }

    /* add the Node (with an incidence of 1) */
    nod[Nnod].ipnt = ipnt;
    nod[Nnod].nedg = 0;

    /* remember the Node at this Point */
    if (pntNod != NULL) {
        pntNod[ipnt] = Nnod;
    }

    /* increment number of Nodes */
    Nnod++;

cleanup:
    return status;
}


/******************************************************************************/

static int
getModelSize(double box[])              /* (out) xyz min, xyz max */
{
    int    status = SUCCESS;            /* (out) return status */

    int    ipnt;

//    ROUTINE(getModelSize);

    /*----------------------------------------------------------------*/

    box[0] = tess.xyz[0];
    box[1] = tess.xyz[1];
    box[2] = tess.xyz[2];
    box[3] = tess.xyz[0];
    box[4] = tess.xyz[1];
    box[5] = tess.xyz[2];

    for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
        if (tess.xyz[3*ipnt  ] < box[0]) box[0] = tess.xyz[3*ipnt  ];
        if (tess.xyz[3*ipnt+1] < box[1]) box[1] = tess.xyz[3*ipnt+1];
        if (tess.xyz[3*ipnt+2] < box[2]) box[2] = tess.xyz[3*ipnt+2];
        if (tess.xyz[3*ipnt  ] > box[3]) box[3] = tess.xyz[3*ipnt  ];
        if (tess.xyz[3*ipnt+1] > box[4]) box[4] = tess.xyz[3*ipnt+1];
        if (tess.xyz[3*ipnt+2] > box[5]) box[5] = tess.xyz[3*ipnt+2];
    }

//cleanup:
    return status;
}


/******************************************************************************/

static int
findEdge(int    xscr,                   /* (in)  x-screen location */
         int    yscr,                   /* (in)  y-screen location */
         int    *ipnt)                  /* (out) index along Edge */
{
    int    iedg;                        /* (out) Point index (bias-0) */

    int    saved_pickmask, jpnt, jedg;
    float  xc, yc;
    double dbest, dtest;

//    ROUTINE(findPoint);

    /*----------------------------------------------------------------*/

    xc   = (float) ((2.0 * xscr) / (gv_w3d.xsize - 1.0) - 1.0);
    yc   = (float) ((2.0 * yscr) / (gv_w3d.ysize - 1.0) - 1.0);

    saved_pickmask = gv_pickmask;
    gv_pickmask = -1;
    LocateGraphic(xc, -yc, 0);
    gv_pickmask = saved_pickmask;

    /* find the closest Point to the located position */
    iedg  = -1;
    *ipnt = -1;
    dbest = HUGEQ;

    for (jedg = 0; jedg < Nedg; jedg++) {
        for (jpnt = 0; jpnt < edg[jedg].npnt; jpnt++) {
            dtest = SQR(tess.xyz[3*edg[jedg].pnt[jpnt]  ] - gv_located[0])
                  + SQR(tess.xyz[3*edg[jedg].pnt[jpnt]+1] - gv_located[1])
                  + SQR(tess.xyz[3*edg[jedg].pnt[jpnt]+2] - gv_located[2]);

            if (dtest < dbest) {
                iedg  = jedg;
                *ipnt = jpnt;
                dbest = dtest;
            }
        }
    }

//cleanup:
    return iedg;
}


/******************************************************************************/

static int
findFace(int    xscr,                   /* (in)  x-screen location */
         int    yscr,                   /* (in)  y-screen location */
         int    *itri)                  /* (out) Triangle index */
{
    int    ifac;                        /* (out) Face index (bias-1) */
                                        /*       -1 if nothing picked */

    int    saved_pickmask;
    float  xc, yc;

//    ROUTINE(findTriangle);

    /*----------------------------------------------------------------*/

    xc   = (float) ((2.0 * xscr) / (gv_w3d.xsize - 1.0) - 1.0);
    yc   = (float) ((2.0 * yscr) / (gv_w3d.ysize - 1.0) - 1.0);

    saved_pickmask = gv_pickmask;
    gv_pickmask = -1;
    PickGraphic(xc, -yc, 0);
    gv_pickmask = saved_pickmask;

    /* nothing is picked */
    if (gv_picked == NULL) {
        ifac  = -1;
        *itri = -1;

    /* picked entity is not the desired type */
    } else if (gv_picked->utype != 4) {
        ifac  = -1;
        *itri = -1;

    /* return the Triangle */
    } else {
        ifac  = gv_picked->uindex;
//$$$        *itri = gv_list[0]->object[0].type.distris.pick;
        *itri = 0;
    }

//cleanup:
    return ifac;
}


/******************************************************************************/

static int
findPoint(int    xscr,                  /* (in)  x-screen location */
          int    yscr)                  /* (in)  y-screen location */
{
    int    ipnt;                        /* (out) Point index (bias-0) */

    int    saved_pickmask, jpnt;
    float  xc, yc;
    double dbest, dtest;

//    ROUTINE(findPoint);

    /*----------------------------------------------------------------*/

    xc   = (float) ((2.0 * xscr) / (gv_w3d.xsize - 1.0) - 1.0);
    yc   = (float) ((2.0 * yscr) / (gv_w3d.ysize - 1.0) - 1.0);

    saved_pickmask = gv_pickmask;
    gv_pickmask = -1;
    LocateGraphic(xc, -yc, 0);
    gv_pickmask = saved_pickmask;

    /* find the closest Point to the located position */
    ipnt = -1;
    dbest = HUGEQ;

    for (jpnt = 0; jpnt < tess.npnt; jpnt++) {
        dtest = SQR(tess.xyz[3*jpnt  ] - gv_located[0])
              + SQR(tess.xyz[3*jpnt+1] - gv_located[1])
              + SQR(tess.xyz[3*jpnt+2] - gv_located[2]);

        if (dtest < dbest) {
            ipnt  = jpnt;
            dbest = dtest;
        }
    }

//cleanup:
    return ipnt;
}


/******************************************************************************/

static int
findTriangle(int    xscr,               /* (in)  x-screen location */
             int    yscr)               /* (in)  y-screen location */
{
    int    itri;                        /* (out) Triangle index (bias-1) */
                                        /*       -1 if nothing picked */

    int    saved_pickmask;
    float  xc, yc;

//    ROUTINE(findTriangle);

    /*----------------------------------------------------------------*/

    xc   = (float) ((2.0 * xscr) / (gv_w3d.xsize - 1.0) - 1.0);
    yc   = (float) ((2.0 * yscr) / (gv_w3d.ysize - 1.0) - 1.0);

    saved_pickmask = gv_pickmask;
    gv_pickmask = -1;
    PickGraphic(xc, -yc, 0);
    gv_pickmask = saved_pickmask;

    /* nothing is picked */
    if (gv_picked == NULL) {
        itri = -1;

    /* picked entity is not the desired type */
    } else if (gv_picked->utype != 1) {
        itri = -1;

    /* return the Triangle */
    } else {
        itri = gv_list[0]->object[0].type.distris.pick;
    }

//cleanup:
    return itri;
}


/******************************************************************************/

static int
makeCurve(int    iedg)                  /* (in)  Edge index (bias-0) */
{
    int    status = SUCCESS;            /* (out) return status */
//    ROUTINE(makeCurve);

    /*----------------------------------------------------------------*/

    printf("makeCurve)iedg=%d\n", iedg);

//cleanup:
    return status;
}


/******************************************************************************/

static int
makeCut(double zcut,                    /* (in)  transformed z for cut */
        int    ifac,                    /* (in)  Face index (bias-0) color */
        double xform[3][3],             /* (in)  transformation matrix */
        int    *ncut,                   /* (out) number of cut points */
        double *xcut[],                 /* (out) transformed x coordinates */
        double *ycut[])                 /* (out) transformed y coordinates */
{
    int    status = SUCCESS;            /* (out) return status */

    typedef struct {
        double    xbeg;
        double    ybeg;
        double    xend;
        double    yend;
        int       prev;
        int       next;
        int       part;
    } seg_TT;

    int    itri, ip0, ip1, ip2, iseg, jseg, mseg, nseg, nlnk, ibest, jbest, part;
    double frac, dtest, dbest, x0, x1, x2, y0, y1, y2, z0, z1, z2;
    seg_TT *seg = NULL;

    ROUTINE(makeCut);

    /*----------------------------------------------------------------*/

    /* default return values */
    *ncut = 0;
    *xcut = NULL;
    *ycut = NULL;

    /* initial Segment table */
    nseg = 0;
    mseg = 1000;
    MALLOC(seg, seg_TT, mseg);

    /* find Segments by looking at Triangles for ifac */
    for (itri = 0; itri < fac[ifac].tess.ntri; itri++) {

        if (nseg > mseg-2) {
            mseg += 1000;
            RALLOC(seg, seg_TT, mseg);
        }

        ip0 = fac[ifac].tess.trip[3*itri  ];
        ip1 = fac[ifac].tess.trip[3*itri+1];
        ip2 = fac[ifac].tess.trip[3*itri+2];

        x0 = xform[0][0] * fac[ifac].tess.xyz[3*ip0  ]
           + xform[0][1] * fac[ifac].tess.xyz[3*ip0+1]
           + xform[0][2] * fac[ifac].tess.xyz[3*ip0+2];
        y0 = xform[1][0] * fac[ifac].tess.xyz[3*ip0  ]
           + xform[1][1] * fac[ifac].tess.xyz[3*ip0+1]
           + xform[1][2] * fac[ifac].tess.xyz[3*ip0+2];
        z0 = xform[2][0] * fac[ifac].tess.xyz[3*ip0  ]
           + xform[2][1] * fac[ifac].tess.xyz[3*ip0+1]
           + xform[2][2] * fac[ifac].tess.xyz[3*ip0+2];

        x1 = xform[0][0] * fac[ifac].tess.xyz[3*ip1  ]
           + xform[0][1] * fac[ifac].tess.xyz[3*ip1+1]
           + xform[0][2] * fac[ifac].tess.xyz[3*ip1+2];
        y1 = xform[1][0] * fac[ifac].tess.xyz[3*ip1  ]
           + xform[1][1] * fac[ifac].tess.xyz[3*ip1+1]
           + xform[1][2] * fac[ifac].tess.xyz[3*ip1+2];
        z1 = xform[2][0] * fac[ifac].tess.xyz[3*ip1  ]
           + xform[2][1] * fac[ifac].tess.xyz[3*ip1+1]
           + xform[2][2] * fac[ifac].tess.xyz[3*ip1+2];

        x2 = xform[0][0] * fac[ifac].tess.xyz[3*ip2  ]
           + xform[0][1] * fac[ifac].tess.xyz[3*ip2+1]
           + xform[0][2] * fac[ifac].tess.xyz[3*ip2+2];
        y2 = xform[1][0] * fac[ifac].tess.xyz[3*ip2  ]
           + xform[1][1] * fac[ifac].tess.xyz[3*ip2+1]
           + xform[1][2] * fac[ifac].tess.xyz[3*ip2+2];
        z2 = xform[2][0] * fac[ifac].tess.xyz[3*ip2  ]
           + xform[2][1] * fac[ifac].tess.xyz[3*ip2+1]
           + xform[2][2] * fac[ifac].tess.xyz[3*ip2+2];

        if        (z0 > zcut && z1 < zcut && z2 < zcut) {
            frac = (zcut - z0) / (z1 - z0);
            seg[nseg].xbeg = (1-frac) * x0 + frac * x1;
            seg[nseg].ybeg = (1-frac) * y0 + frac * y1;

            frac = (zcut - z2) / (z0 - z2);
            seg[nseg].xend = (1-frac) * x2 + frac * x0;
            seg[nseg].yend = (1-frac) * y2 + frac * y0;

            seg[nseg].prev = -1;
            seg[nseg].next = -1;
            nseg++;
        } else if (z0 < zcut && z1 > zcut && z2 > zcut) {
            frac = (zcut - z2) / (z0 - z2);
            seg[nseg].xbeg = (1-frac) * x2 + frac * x0;
            seg[nseg].ybeg = (1-frac) * y2 + frac * y0;

            frac = (zcut - z0) / (z1 - z0);
            seg[nseg].xend = (1-frac) * x0 + frac * x1;
            seg[nseg].yend = (1-frac) * y0 + frac * y1;

            seg[nseg].prev = -1;
            seg[nseg].next = -1;
            nseg++;
        } else if (z0 < zcut && z1 > zcut && z2 < zcut) {
            frac = (zcut - z1) / (z2 - z1);
            seg[nseg].xbeg = (1-frac) * x1 + frac * x2;
            seg[nseg].ybeg = (1-frac) * y1 + frac * y2;

            frac = (zcut - z0) / (z1 - z0);
            seg[nseg].xend = (1-frac) * x0 + frac * x1;
            seg[nseg].yend = (1-frac) * y0 + frac * y1;

            seg[nseg].prev = -1;
            seg[nseg].next = -1;
            nseg++;
        } else if (z0 > zcut && z1 < zcut && z2 > zcut) {
            frac = (zcut - z0) / (z1 - z0);
            seg[nseg].xbeg = (1-frac) * x0 + frac * x1;
            seg[nseg].ybeg = (1-frac) * y0 + frac * y1;

            frac = (zcut - z1) / (z2 - z1);
            seg[nseg].xend = (1-frac) * x1 + frac * x2;
            seg[nseg].yend = (1-frac) * y1 + frac * y2;

            seg[nseg].prev = -1;
            seg[nseg].next = -1;
            nseg++;
        } else if (z0 < zcut && z1 < zcut && z2 > zcut) {
            frac = (zcut - z2) / (z0 - z2);
            seg[nseg].xbeg = (1-frac) * x2 + frac * x0;
            seg[nseg].ybeg = (1-frac) * y2 + frac * y0;

            frac = (zcut - z1) / (z2 - z1);
            seg[nseg].xend = (1-frac) * x1 + frac * x2;
            seg[nseg].yend = (1-frac) * y1 + frac * y2;

            seg[nseg].prev = -1;
            seg[nseg].next = -1;
            nseg++;
        } else if (z0 > zcut && z1 > zcut && z2 < zcut) {
            frac = (zcut - z1) / (z2 - z1);
            seg[nseg].xbeg = (1-frac) * x1 + frac * x2;
            seg[nseg].ybeg = (1-frac) * y1 + frac * y2;

            frac = (zcut - z2) / (z0 - z2);
            seg[nseg].xend = (1-frac) * x2 + frac * x0;
            seg[nseg].yend = (1-frac) * y2 + frac * y0;

            seg[nseg].prev = -1;
            seg[nseg].next = -1;
            nseg++;
        }
    }

    /* link Segments with coincident endpoints */
    nlnk = 0;
    for (iseg = 0; iseg < nseg; iseg++) {
        if (seg[iseg].next >= 0) continue;

        for (jseg = 0; jseg < nseg; jseg++) {
            if (iseg == jseg || seg[jseg].prev >= 0) continue;

            if (fabs(seg[iseg].xend - seg[jseg].xbeg) < EPS06 &&
                fabs(seg[iseg].yend - seg[jseg].ybeg) < EPS06   ) {
                seg[iseg].next = jseg;
                seg[jseg].prev = iseg;
                nlnk++;
                break;
            }
        }
    }

    /* if we are going to need pseudo-Segments, assign an index
       to each part of the intersection curve */
    if (nlnk < nseg-1) {
        part = 0;
        for (iseg = 0; iseg < nseg; iseg++) {
            seg[iseg].part = -1;
        }

        for (jseg = 0; jseg < nseg; jseg++) {
            if (seg[jseg].prev < 0 && seg[jseg].part < 0) {

                iseg = jseg;
                while (1) {
                    seg[iseg].part = part;
                    iseg = seg[iseg].next;
                    if (iseg < 0) break;
                }
                part++;
            }
        }
    }

    /* if nlnk is less than nseg-1, we need to create
       pseudo-Segments and pairs of links */
    while (nlnk < nseg-1) {
        dbest = HUGEQ;
        ibest = -1;
        jbest = -1;

        for (iseg = 0; iseg < nseg; iseg++) {
            if (seg[iseg].next >= 0) continue;

            for (jseg = 0; jseg < nseg; jseg++) {
                if (seg[jseg].prev >= 0) continue;
                if (seg[iseg].part == seg[jseg].part) continue;

                dtest = MAX(fabs(seg[iseg].xend - seg[jseg].xbeg),
                            fabs(seg[iseg].yend - seg[jseg].ybeg));
                if (dtest < dbest) {
                    dbest = dtest;
                    ibest = iseg;
                    jbest = jseg;
                }
            }
        }

        if (ibest < 0 || jbest < 0) {
            printf("ERR: could not find place for pseudo-Segment\a\n");
            goto cleanup;
        }

        if (nseg > mseg-2) {
            mseg += 1000;
            RALLOC(seg, seg_TT, mseg);
        }

        seg[nseg].xbeg = seg[ibest].xend;
        seg[nseg].ybeg = seg[ibest].yend;
        seg[nseg].xend = seg[jbest].xbeg;
        seg[nseg].yend = seg[jbest].ybeg;

        seg[nseg ].prev = ibest;
        seg[ibest].next = nseg;

        seg[nseg ].next = jbest;
        seg[jbest].prev = nseg;

        nseg++;
        nlnk += 2;

        /* combined the two parts */
        ibest = seg[ibest].part;
        jbest = seg[jbest].part;

        for (iseg = 0; iseg < nseg; iseg++) {
            if (seg[iseg].part == jbest) {
                seg[iseg].part =  ibest;
            }
        }
    }

    /* allocate memory woth room to add possibly two additional points */
    MALLOC(*xcut, double, nseg+3);
    MALLOC(*ycut, double, nseg+3);

    *ncut = 0;

    /* find the Segment that has an unfilled .prev */
    ibest = -1;
    for (iseg = 0; iseg < nseg; iseg++) {
        if (seg[iseg].prev < 0) {
            ibest = iseg;
            break;
        }
    }

    if (ibest < 0) {
        printf("ERR: no Segment with unfilled .prev\a\n");
        goto cleanup;
    }

    /* build the output arrays */
    iseg = ibest;
    (*xcut)[*ncut] = seg[iseg].xbeg;
    (*ycut)[*ncut] = seg[iseg].ybeg;
    (*ncut)++;

    while (1) {
        (*xcut)[*ncut] = seg[iseg].xend;
        (*ycut)[*ncut] = seg[iseg].yend;
        (*ncut)++;

        iseg = seg[iseg].next;
        if (iseg < 0) break;
    }

    if (*ncut != (nseg+1)) {
        printf("ERR: *ncut != (nseg+1)\a\n");
        goto cleanup;
    }

cleanup:
    FREE(seg);
    return status;
}


/******************************************************************************/

static int
makeSurface(int    ifac,                /* (in)  Face index (bias-0) */
            int    itype,               /* (in)  parameterization type */
            int    imax,                /* (in)  points between west and east */
            int    jmax)                /* (in)  points between south and north */
{
    int    status = SUCCESS;            /* (out) return status */

    int    iedg, jedg, ipnt, i, j, k, swap, last;
    int    senw[4], nnn, nlist, list[100];
    double *xyzs = NULL;

    ROUTINE(makeSurface);

    /*----------------------------------------------------------------*/

    printf("\n\nmakeSurface(ifac=%d, itype=%d)\n", ifac, itype);

    senw[0] = -1;                                 // needed to avoid clang warning
    senw[1] = -1;
    senw[2] = -1;
    senw[3] = -1;

    /* it is an error to have no Edges with .mark=0 */
    last = 0;
    for (iedg = 0; iedg < Nedg; iedg++) {
        if (edg[iedg].mark == 0) last = 1;
    }

    if (last == 0) {
        printf("ERR: no Edges with ,mark=0\a\n");
        goto cleanup;
    }

    /* make a list of all marked Edges */
    nlist = 0;
    for (iedg = 0; iedg < Nedg; iedg++) {
        if (edg[iedg].mark >= 0 && edg[iedg].mark <= 3) {
            list[2*nlist] = iedg;
            if        (edg[iedg].ileft == ifac) {
                list[2*nlist+1] = +1;
            } else if (edg[iedg].irite == ifac) {
                list[2*nlist+1] = -1;
            } else {
                list[2*nlist+1] =  0;
            }
            nlist++;
        }
    }

    /* make sure that the list starts with an Edge that is associated
       with ifac */
    for (i = 0; i < nlist; i++) {
        if (list[2*i+1] != 0) {
            swap        = list[0    ];
            list[0    ] = list[2*i  ];
            list[2*i  ] = swap;

            swap        = list[1];
            list[1    ] = list[2*i+1];
            list[2*i+1] = swap;

            break;
        }
    }

    /* order the list head to tail */
    for (i = 0; i < nlist; i++) {
        iedg = list[2*i];
        for (j = i+1; j < nlist; j++) {
            jedg = list[2*j];

            if (list[2*i+1] >= 0 && list[2*j+1] >= 0 && edg[iedg].iend == edg[jedg].ibeg) {
                list[2*i+1] = +1;
                list[2*j+1] = +1;

                swap        = list[2*i+2];
                list[2*i+2] = list[2*j  ];
                list[2*j  ] = swap;

                swap        = list[2*i+3];
                list[2*i+3] = list[2*j+1];
                list[2*j+1] = swap;

                break;
            }
            if (list[2*i+1] >= 0 && list[2*j+1] <= 0 && edg[iedg].iend == edg[jedg].iend) {
                list[2*i+1] = +1;
                list[2*j+1] = -1;

                swap        = list[2*i+2];
                list[2*i+2] = list[2*j  ];
                list[2*j  ] = swap;

                swap        = list[2*i+3];
                list[2*i+3] = list[2*j+1];
                list[2*j+1] = swap;

                break;
            }
            if (list[2*i+1] <= 0 && list[2*j+1] >= 0 && edg[iedg].ibeg == edg[jedg].ibeg) {
                list[2*i+1] = -1;
                list[2*j+1] = +1;

                swap        = list[2*i+2];
                list[2*i+2] = list[2*j  ];
                list[2*j  ] = swap;

                swap        = list[2*i+3];
                list[2*i+3] = list[2*j+1];
                list[2*j+1] = swap;

                break;
            }
            if (list[2*i+1] <= 0 && list[2*j+1] <= 0 && edg[iedg].ibeg == edg[jedg].iend) {
                list[2*i+1] = -1;
                list[2*j+1] = -1;

                swap        = list[2*i+2];
                list[2*i+2] = list[2*j  ];
                list[2*j  ] = swap;

                swap        = list[2*i+3];
                list[2*i+3] = list[2*j+1];
                list[2*j+1] = swap;

                break;
            }
        }
    }

    /* count the total number of Points in all marked Edges */
    nnn = 0;
    for (iedg = 0; iedg < Nedg; iedg++) {
        if (edg[iedg].mark >= 0 && edg[iedg].mark <= 3) {
            nnn += edg[iedg].npnt;
        }
    }

    /* allocate an array to hold the Points in the marked Edges */
    MALLOC(xyzs, double, 4*nnn);

    /* find first list entity with mark==0 (south) following
       another list entity with mark!=0 (not south) */
    last = -1;
    for (i = 0; i < nlist; i++) {
        iedg = list[2*i];

        j = (i+nlist-1) % nlist;
        jedg = list[2*j];
        if (edg[jedg].mark != 0 && edg[iedg].mark == 0) {
            nnn     =  0;
            senw[0] = -1;
            senw[1] = -1;
            senw[2] = -1;
            senw[3] = -1;

            for (j = 0; j < nlist; j++) {
                jedg = list[2*((i+j)%nlist)];

                if (edg[jedg].mark != last) {
                    last = edg[jedg].mark;

                    if (list[2*((i+j)%nlist)+1] > 0) {
                        ipnt = 0;
                    } else {
                        ipnt = edg[jedg].npnt - 1;
                    }

                    xyzs[4*nnn  ] = tess.xyz[3*edg[jedg].pnt[ipnt]  ];
                    xyzs[4*nnn+1] = tess.xyz[3*edg[jedg].pnt[ipnt]+1];
                    xyzs[4*nnn+2] = tess.xyz[3*edg[jedg].pnt[ipnt]+2];
                    xyzs[4*nnn+3] = 0;
                    nnn++;
                }

                for (k = 1; k < edg[jedg].npnt; k++) {
                    if (list[2*((i+j)%nlist)+1] > 0) {
                        ipnt = edg[jedg].pnt[k];
                    } else {
                        ipnt = edg[jedg].pnt[edg[jedg].npnt-k-1];
                    }

                    xyzs[4*nnn  ] = tess.xyz[3*ipnt  ];
                    xyzs[4*nnn+1] = tess.xyz[3*ipnt+1];
                    xyzs[4*nnn+2] = tess.xyz[3*ipnt+2];
                    xyzs[4*nnn+3] = xyzs[4*nnn-1] + sqrt(SQR(xyzs[4*nnn  ]-xyzs[4*nnn-4])
                                                        +SQR(xyzs[4*nnn+1]-xyzs[4*nnn-3])
                                                        +SQR(xyzs[4*nnn+2]-xyzs[4*nnn-2]));
                    nnn++;

                    senw[last] = nnn;
                }
            }
            break;
        }
    }

    if (senw[0] < 0) senw[0] = 0;
    if (senw[1] < 0) senw[1] = senw[0];
    if (senw[2] < 0) senw[2] = senw[1];
    if (senw[3] < 0) senw[3] = senw[2];

    /* allocate space for the surface */
    fac[ifac].imax = imax;
    fac[ifac].jmax = jmax;

    MALLOC(fac[ifac].xsrf, double, imax*jmax);
    MALLOC(fac[ifac].ysrf, double, imax*jmax);
    MALLOC(fac[ifac].zsrf, double, imax*jmax);

    /* fill the Grid based upon the type */
    if (itype == 1) {
        makeSurface1(ifac, imax, jmax, xyzs, senw);
    } else if (itype == 2) {
        makeSurface2(ifac, imax, jmax, xyzs, senw);
    }

    /* unmark all Edges */
    for (iedg = 0; iedg < Nedg; iedg++) {
        edg[iedg].mark = -1;
    }

 cleanup:
    FREE(xyzs);

    return status;
}


/******************************************************************************/

static int
makeSurface1(int    ifac,               /* (in)  Face index (bias-0) */
             int    imax,               /* (in)  points between west and east */
             int    jmax,               /* (in)  points between south and north */
             double xyzs[],             /* (in)  points around Edges */
             int    senw[])             /* (in)  number of points on each side */
{
    int    status = SUCCESS;            /* (out) return status */

    int    i, j, k, ij, ibeg, iend, jbeg, jend, iedg, ipnt, jpnt, ntri_old, itri;
    double stgt, fraci, fracj, xyz_in[3], xyz_out[3], xedg, yedg, zedg;

    ROUTINE(makeSurface1);

    /*----------------------------------------------------------------*/

    printf("\n\nmakeSurface1(ifac=%d)\n", ifac);

    /* find Edges that are associated with ifac but which are not marked */
    ipnt = -1;
    for (iedg = 0; iedg < Nedg; iedg++) {
        if (edg[iedg].ileft == ifac || edg[iedg].irite == ifac) {
            if (edg[iedg].mark < 0) {
                xedg = tess.xyz[3*edg[iedg].pnt[1]  ];
                yedg = tess.xyz[3*edg[iedg].pnt[1]+1];
                zedg = tess.xyz[3*edg[iedg].pnt[1]+2];
                for (jpnt = 0; jpnt < fac[ifac].tess.npnt; jpnt++) {
                    if (fabs(fac[ifac].tess.xyz[3*jpnt  ]-xedg) < EPS06 &&
                        fabs(fac[ifac].tess.xyz[3*jpnt+1]-yedg) < EPS06 &&
                        fabs(fac[ifac].tess.xyz[3*jpnt+2]-zedg) < EPS06   ) {
                        ipnt = jpnt;
                        printf("Edge %d is probably part of an inner loop, ipnt=%d\n", iedg, ipnt);
                        break;
                    }
                }
            }
        }
    }

    /* fill the Loop that contains ipnt */
    ntri_old = fac[ifac].tess.ntri;
    if (ipnt >= 0) {
        status = fillLoop(&(fac[ifac].tess), ipnt);
        CHECK_STATUS(fillLoop);
    }

    /* equidistribute points along nondegenerate south Edges */
    if (senw[0] > 0) {
        k = 1;
        j = 0;
        for (i = 0; i < imax; i++) {
            stgt = xyzs[4*senw[0]-1] * (double)(i) / (double)(imax-1);

            while (stgt > xyzs[4*k+3] && k < senw[0]-1) {
                k++;
            }

            fraci = (stgt - xyzs[4*k-1]) / (xyzs[4*k+3] - xyzs[4*k-1]);
            fac[ifac].xsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-4] + fraci * xyzs[4*k  ];
            fac[ifac].ysrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-3] + fraci * xyzs[4*k+1];
            fac[ifac].zsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-2] + fraci * xyzs[4*k+2];
        }
    }

    /* equidistribute points along nondegenerate east Edges */
    if (senw[1] > senw[0]) {
        k = senw[0] + 1;
        i = imax - 1;
        for (j = 0; j < jmax; j++) {
            stgt = xyzs[4*senw[1]-1] * (double)(j) / (double)(jmax-1);

            while (stgt > xyzs[4*k+3] && k < senw[1]-1) {
                k++;
            }

            fraci = (stgt - xyzs[4*k-1]) / (xyzs[4*k+3] - xyzs[4*k-1]);
            fac[ifac].xsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-4] + fraci * xyzs[4*k  ];
            fac[ifac].ysrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-3] + fraci * xyzs[4*k+1];
            fac[ifac].zsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-2] + fraci * xyzs[4*k+2];
        }
    }

    /* equidistribute points along nondegenerate north Edges */
    if (senw[2] > senw[1]) {
        k = senw[1] + 1;
        j = jmax - 1;
        for (i = imax-1; i >= 0; i--) {
            stgt = xyzs[4*senw[2]-1] * (double)(imax-1-i) / (double)(imax-1);

            while (stgt > xyzs[4*k+3] && k < senw[2]-1) {
                k++;
            }

            fraci = (stgt - xyzs[4*k-1]) / (xyzs[4*k+3] - xyzs[4*k-1]);
            fac[ifac].xsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-4] + fraci * xyzs[4*k  ];
            fac[ifac].ysrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-3] + fraci * xyzs[4*k+1];
            fac[ifac].zsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-2] + fraci * xyzs[4*k+2];
        }
    }

    /* equidistribute points along nondegenerate west Edges */
    if (senw[3] > senw[2]) {
        k = senw[2] + 1;
        i = 0;
        for (j = jmax-1; j >= 0; j--) {
            stgt = xyzs[4*senw[3]-1] * (double)(jmax-1-j) / (double)(jmax-1);

            while (stgt > xyzs[4*k+3] && k < senw[3]-1) {
                k++;
            }

            fraci = (stgt - xyzs[4*k-1]) / (xyzs[4*k+3] - xyzs[4*k-1]);
            fac[ifac].xsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-4] + fraci * xyzs[4*k  ];
            fac[ifac].ysrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-3] + fraci * xyzs[4*k+1];
            fac[ifac].zsrf[(i)+imax*(j)] = (1-fraci) * xyzs[4*k-2] + fraci * xyzs[4*k+2];
        }
    }

    /* if any side have no Edges, set up points by linear interpolation */
    if (senw[0] <= 0) {
        j = 0;
        for (i = 1; i < imax-1; i++) {
            fraci = (double)(i) / (double)(imax-1);

            ij = i + imax * j;
            fac[ifac].xsrf[ij] = (1-fraci) * fac[ifac].xsrf[(0     )+imax*(j)]
                               +    fraci  * fac[ifac].xsrf[(imax-1)+imax*(j)];
            fac[ifac].ysrf[ij] = (1-fraci) * fac[ifac].ysrf[(0     )+imax*(j)]
                               +    fraci  * fac[ifac].ysrf[(imax-1)+imax*(j)];
            fac[ifac].zsrf[ij] = (1-fraci) * fac[ifac].zsrf[(0     )+imax*(j)]
                               +    fraci  * fac[ifac].zsrf[(imax-1)+imax*(j)];
        }
    }

    if (senw[1] <= senw[0]) {
        i = imax - 1;
        for (j = 1; j < jmax-1; j++) {
            fracj = (double)(j) / (double)(jmax-1);

            ij = i + imax * j;
            fac[ifac].xsrf[ij] = (1-fracj) * fac[ifac].xsrf[(i)+imax*(0     )]
                               +    fracj  * fac[ifac].xsrf[(i)+imax*(jmax-1)];
            fac[ifac].ysrf[ij] = (1-fracj) * fac[ifac].ysrf[(i)+imax*(0     )]
                               +    fracj  * fac[ifac].ysrf[(i)+imax*(jmax-1)];
            fac[ifac].zsrf[ij] = (1-fracj) * fac[ifac].zsrf[(i)+imax*(0     )]
                               +    fracj  * fac[ifac].zsrf[(i)+imax*(jmax-1)];
        }
    }

    if (senw[2] <= senw[1]) {
        j = jmax - 1;
        for (i = 1; i < imax-1; i++) {
            fraci = (double)(i) / (double)(imax-1);

            ij = i + imax * j;
            fac[ifac].xsrf[ij] = (1-fraci) * fac[ifac].xsrf[(0     )+imax*(jmax-1)]
                               +    fraci  * fac[ifac].xsrf[(imax-1)+imax*(jmax-1)];
            fac[ifac].ysrf[ij] = (1-fraci) * fac[ifac].ysrf[(0     )+imax*(jmax-1)]
                               +    fraci  * fac[ifac].ysrf[(imax-1)+imax*(jmax-1)];
            fac[ifac].zsrf[ij] = (1-fraci) * fac[ifac].zsrf[(0     )+imax*(jmax-1)]
                               +    fraci  * fac[ifac].zsrf[(imax-1)+imax*(jmax-1)];
        }
    }

    if (senw[3] <= senw[2]) {
        i = 0;
        for (j = 1; j < jmax-1; j++) {
            fracj = (double)(j) / (double)(jmax-1);

            ij = i + imax * j;
            fac[ifac].xsrf[ij] = (1-fracj) * fac[ifac].xsrf[(i)+imax*(0     )]
                               +    fracj  * fac[ifac].xsrf[(i)+imax*(jmax-1)];
            fac[ifac].ysrf[ij] = (1-fracj) * fac[ifac].ysrf[(i)+imax*(0     )]
                               +    fracj  * fac[ifac].ysrf[(i)+imax*(jmax-1)];
            fac[ifac].zsrf[ij] = (1-fracj) * fac[ifac].zsrf[(i)+imax*(0     )]
                               +    fracj  * fac[ifac].zsrf[(i)+imax*(jmax-1)];
        }
    }

    /* interior points with projection (in rings) */
    ibeg = 0;
    iend = imax - 1;
    jbeg = 0;
    jend = jmax - 1;

    while (iend > ibeg || jend > jbeg) {
        for (j = jbeg+1; j < jend; j++) {
            for (i = ibeg+1; i < iend; i++) {

                if (i > ibeg+1 && i < iend-1 && j > jbeg+1 && j < jend-1) continue;

                fraci = (double)(i-ibeg) / (double)(iend-ibeg);
                fracj = (double)(j-jbeg) / (double)(jend-jbeg);

                xyz_in[0] = (1-fraci)             * fac[ifac].xsrf[(ibeg)+imax*(j   )]
                          +    fraci              * fac[ifac].xsrf[(iend)+imax*(j   )]
                          +             (1-fracj) * fac[ifac].xsrf[(i   )+imax*(jbeg)]
                          +                fracj  * fac[ifac].xsrf[(i   )+imax*(jend)]
                          - (1-fraci) * (1-fracj) * fac[ifac].xsrf[(ibeg)+imax*(jbeg)]
                          -    fraci  * (1-fracj) * fac[ifac].xsrf[(iend)+imax*(jbeg)]
                          - (1-fraci) *    fracj  * fac[ifac].xsrf[(ibeg)+imax*(jend)]
                          -    fraci  *    fracj  * fac[ifac].xsrf[(iend)+imax*(jend)];
                xyz_in[1] = (1-fraci)             * fac[ifac].ysrf[(ibeg)+imax*(j   )]
                          +    fraci              * fac[ifac].ysrf[(iend)+imax*(j   )]
                          +             (1-fracj) * fac[ifac].ysrf[(i   )+imax*(jbeg)]
                          +                fracj  * fac[ifac].ysrf[(i   )+imax*(jend)]
                          - (1-fraci) * (1-fracj) * fac[ifac].ysrf[(ibeg)+imax*(jbeg)]
                          -    fraci  * (1-fracj) * fac[ifac].ysrf[(iend)+imax*(jbeg)]
                          - (1-fraci) *    fracj  * fac[ifac].ysrf[(ibeg)+imax*(jend)]
                          -    fraci  *    fracj  * fac[ifac].ysrf[(iend)+imax*(jend)];
                xyz_in[2] = (1-fraci)             * fac[ifac].zsrf[(ibeg)+imax*(j   )]
                          +    fraci              * fac[ifac].zsrf[(iend)+imax*(j   )]
                          +             (1-fracj) * fac[ifac].zsrf[(i   )+imax*(jbeg)]
                          +                fracj  * fac[ifac].zsrf[(i   )+imax*(jend)]
                          - (1-fraci) * (1-fracj) * fac[ifac].zsrf[(ibeg)+imax*(jbeg)]
                          -    fraci  * (1-fracj) * fac[ifac].zsrf[(iend)+imax*(jbeg)]
                          - (1-fraci) *    fracj  * fac[ifac].zsrf[(ibeg)+imax*(jend)]
                          -    fraci  *    fracj  * fac[ifac].zsrf[(iend)+imax*(jend)];

                projectToFace(ifac, xyz_in, xyz_out);

                ij = i + imax * j;
                fac[ifac].xsrf[ij] = xyz_out[0];
                fac[ifac].ysrf[ij] = xyz_out[1];
                fac[ifac].zsrf[ij] = xyz_out[2];
            }
        }

        ibeg++;
        iend--;
        jbeg++;
        jend--;
    }

    /* remove Triangles added during the fill above */
    for (itri = fac[ifac].tess.ntri-1; itri >= ntri_old; itri--) {
        status = deleteTriangle(&(fac[ifac].tess), itri);
        CHECK_STATUS(deleteTriangle);
    }

cleanup:
    return status;
}


/******************************************************************************/

static int
makeSurface2(int    ifac,               /* (in)  Face index (bias-0) */
             int    imax,               /* (in)  points between west and east */
             int    jmax,               /* (in)  points between south and north */
             double xyzs[],             /* (in)  points around Edges */
             int    senw[])             /* (in)  number of points on each side */
{
    int    status = SUCCESS;            /* (out) return status */

    int    ipnt, i, j, k, ncut;
    double xwest, ywest, zwest, xeast, yeast, zeast;
    double stot, stgt, sbeg, send, fraci, fracj, xcut, ycut, zcut, zcut1, zcut2;
    double xform[3][3], dx, dy, dz, zcmin, zcmax;
    double *x = NULL, *y = NULL;

//    ROUTINE(makeSurface2);

    /*----------------------------------------------------------------*/

    printf("\n\nmakeSurface2(ifac=%d)\n", ifac);

    /* find the midpoint along the west and east Edges */
    xwest = (xyzs[4*senw[3]-4] + xyzs[4*senw[2]-4]) / 2;
    ywest = (xyzs[4*senw[3]-3] + xyzs[4*senw[2]-3]) / 2;
    zwest = (xyzs[4*senw[3]-2] + xyzs[4*senw[2]-2]) / 2;

    xeast = (xyzs[4*senw[0]-4] + xyzs[4*senw[1]-4]) / 2;
    yeast = (xyzs[4*senw[0]-3] + xyzs[4*senw[1]-3]) / 2;
    zeast = (xyzs[4*senw[0]-2] + xyzs[4*senw[1]-2]) / 2;

    /* take cuts perpendicular to the line between the west
       and east midpoints */
    dx = xeast - xwest;
    dy = yeast - ywest;
    dz = zeast - zwest;

    /* constant x cuts */
    if (fabs(dx) > fabs(dy) && fabs(dx) > fabs(dz)) {
        xform[0][0] = 0;
        xform[0][1] = 1;
        xform[0][2] = 0;
        xform[1][0] = 0;
        xform[1][1] = 0;
        xform[1][2] = 1;
        xform[2][0] = 1;
        xform[2][1] = 0;
        xform[2][2] = 0;

    /* constant y cuts */
    } else if (fabs(dy) > fabs(dx) && fabs(dy) > fabs(dz)) {
        xform[0][0] = 0;
        xform[0][1] = 0;
        xform[0][2] = 1;
        xform[1][0] = 1;
        xform[1][1] = 0;
        xform[1][2] = 0;
        xform[2][0] = 0;
        xform[2][1] = 1;
        xform[2][2] = 0;

    /* constant z cuts */
    } else {
        xform[0][0] = 1;
        xform[0][1] = 0;
        xform[0][2] = 0;
        xform[1][0] = 0;
        xform[1][1] = 1;
        xform[1][2] = 0;
        xform[2][0] = 0;
        xform[2][1] = 0;
        xform[2][2] = 1;
    }

    /* find the extrema in the cut coordinates */
    zcmin = +HUGEQ;
    zcmax = -HUGEQ;

    for (ipnt = 0; ipnt < fac[ifac].tess.npnt; ipnt++) {
        zcut = xform[2][0] * fac[ifac].tess.xyz[3*ipnt  ]
             + xform[2][1] * fac[ifac].tess.xyz[3*ipnt+1]
             + xform[2][2] * fac[ifac].tess.xyz[3*ipnt+2];

        if (zcut < zcmin) zcmin = zcut;
        if (zcut > zcmax) zcmax = zcut;
    }

    /* the i=0 Grid plane is equispaced along the west Edges */
    i = 0;
    for (j = 0; j < jmax; j++) {
        fac[ifac].xsrf[i+imax*j] = xwest;
        fac[ifac].ysrf[i+imax*j] = ywest;
        fac[ifac].zsrf[i+imax*j] = zwest;
    }

    /* make imax cuts perpendicular to the line between points W and E
       and equi-space between W and E */
    for (i = 1; i < imax-1; i++) {
        fraci = (double)(i) / (double)(imax-1);
        zcut  = (1-fraci) * zcmin + fraci * zcmax;

        makeCut(zcut, ifac, xform, &ncut, &x, &y);

#ifdef GRAFIC
        if (0) {
            int     io_kbd = 5;
            int     io_scr = 6;
            int     indgr  = 1 + 2 + 4 + 16 + 64;
            int     ilin[10], isym[10], nper[10];
            int     nline, m;
            float   xplot[10000], yplot[10000];
            char    pltitl[80];

            sprintf(pltitl, "~xbar~ybar~ zcut=%f", zcut);

            for (m = 0; m < ncut; m++) {
                xplot[m] = (float)x[m];
                yplot[m] = (float)y[m];
            }
            ilin[0] = +1;
            isym[0] = -1;
            nper[0] = ncut;

            xplot[ncut] = xplot[ncut-1];
            yplot[ncut] = yplot[ncut-1];
            ilin[1] = -2;
            isym[1] = +2;
            nper[1] =  1;

            nline = 2;

            grinit_(&io_kbd, &io_scr, "makeSurface2", strlen("makeSurface2"));
            grline_(ilin, isym, &nline, pltitl, &indgr, xplot, yplot, nper, strlen(pltitl));
        }
#endif

        /* find the intersection of the south Edges with the cut plane */
        k = 0;
        xcut  = xform[0][0] * xyzs[4*k] + xform[0][1] * xyzs[4*k+1] + xform[0][2] * xyzs[4*k+2];
        ycut  = xform[1][0] * xyzs[4*k] + xform[1][1] * xyzs[4*k+1] + xform[1][2] * xyzs[4*k+2];
        zcut1 = xform[2][0] * xyzs[4*k] + xform[2][1] * xyzs[4*k+1] + xform[2][2] * xyzs[4*k+2];
        for (k = 1; k < senw[0]; k++) {
            zcut2 = xform[2][0] * xyzs[4*k] + xform[2][1] * xyzs[4*k+1] + xform[2][2] * xyzs[4*k+2];
            fracj = (zcut - zcut1) / (zcut2 - zcut1);

            if (fracj >= -EPS06 && fracj <= 1+EPS06) {
                xcut  = (1-fracj) * (xform[0][0] * xyzs[4*k-4] + xform[0][1] * xyzs[4*k-3] + xform[0][2] * xyzs[4*k-2])
                      +    fracj  * (xform[0][0] * xyzs[4*k  ] + xform[0][1] * xyzs[4*k+1] + xform[0][2] * xyzs[4*k+2]);
                ycut  = (1-fracj) * (xform[1][0] * xyzs[4*k-4] + xform[1][1] * xyzs[4*k-3] + xform[1][2] * xyzs[4*k-2])
                      +    fracj  * (xform[1][0] * xyzs[4*k  ] + xform[1][1] * xyzs[4*k+1] + xform[1][2] * xyzs[4*k+2]);

                break;
            }
            zcut1 = zcut2;
        }

        /* if the south intersection is not part of the cut, add it to the
           beginning of the end (whichever is closest) */
        if        (fabs(x[0]-xcut) < EPS06 && fabs(y[0]-ycut) < EPS06) {

        } else if (fabs(x[ncut-1]-xcut) < EPS06 && fabs(y[ncut-1]-ycut) < EPS06) {

        } else if (MAX(fabs(x[0]-xcut), fabs(y[0]-ycut)) <
                   MAX(fabs(x[ncut-1]-xcut), fabs(y[ncut-1]-ycut))) {
            for (k = ncut; k > 0; k--) {
                x[k] = x[k-1];
                y[k] = y[k-1];
            }
            x[0] = xcut;
            y[0] = ycut;
            ncut++;
        } else {
            x[ncut] = xcut;
            y[ncut] = ycut;
            ncut++;
        }

        /* find the intersection of the north Edges with the cut plane */
        k = senw[1];
        xcut  = xform[0][0] * xyzs[4*k] + xform[0][1] * xyzs[4*k+1] + xform[0][2] * xyzs[4*k+2];
        ycut  = xform[1][0] * xyzs[4*k] + xform[1][1] * xyzs[4*k+1] + xform[1][2] * xyzs[4*k+2];
        zcut1 = xform[2][0] * xyzs[4*k] + xform[2][1] * xyzs[4*k+1] + xform[2][2] * xyzs[4*k+2];
        for (k = senw[1]+1; k < senw[2]; k++) {
            zcut2 = xform[2][0] * xyzs[4*k] + xform[2][1] * xyzs[4*k+1] + xform[2][2] * xyzs[4*k+2];
            fracj = (zcut - zcut1) / (zcut2 - zcut1);

            if (fracj >= -EPS06 && fracj <= 1+EPS06) {
                xcut  = (1-fracj) * (xform[0][0] * xyzs[4*k-4] + xform[0][1] * xyzs[4*k-3] + xform[0][2] * xyzs[4*k-2])
                      +    fracj  * (xform[0][0] * xyzs[4*k  ] + xform[0][1] * xyzs[4*k+1] + xform[0][2] * xyzs[4*k+2]);
                ycut  = (1-fracj) * (xform[1][0] * xyzs[4*k-4] + xform[1][1] * xyzs[4*k-3] + xform[1][2] * xyzs[4*k-2])
                      +    fracj  * (xform[1][0] * xyzs[4*k  ] + xform[1][1] * xyzs[4*k+1] + xform[1][2] * xyzs[4*k+2]);

                break;
            }
            zcut1 = zcut2;
        }

        /* if the north intersection is not part of the cut, add it to the
           beginning of the end (whichever is closest) */
        if        (fabs(x[0]-xcut) < EPS06 && fabs(y[0]-ycut) < EPS06) {

        } else if (fabs(x[ncut-1]-xcut) < EPS06 && fabs(y[ncut-1]-ycut) < EPS06) {

        } else if (MAX(fabs(x[0]-xcut), fabs(y[0]-ycut)) <
                   MAX(fabs(x[ncut-1]-xcut), fabs(y[ncut-1]-ycut))) {
            for (k = ncut; k > 0; k--) {
                x[k] = x[k-1];
                y[k] = y[k-1];
            }
            x[0] = xcut;
            y[0] = ycut;
            ncut++;
        } else {
            x[ncut] = xcut;
            y[ncut] = ycut;
            ncut++;
        }

        /* compute the total length of the cut */
        stot = 0;
        for (k = 1; k < ncut; k++) {
            stot += sqrt(SQR(x[k]-x[k-1]) + SQR(y[k]-y[k-1]));
        }

        /* equidistribute points along this cut */
        k    = 1;
        sbeg = 0;
        send = sqrt(SQR(x[k]-x[k-1]) + SQR(y[k]-y[k-1]));

        for (j = 0; j < jmax; j++) {
            stgt  = stot * (double)(j) / (double)(jmax-1);

            while (send < stgt && k < ncut-1) {
                k++;
                sbeg = send;
                send = sbeg + sqrt(SQR(x[k]-x[k-1]) + SQR(y[k]-y[k-1]));
            }

            fracj = (stgt - sbeg) / (send - sbeg);
            xcut  = (1-fracj) * x[k-1] + fracj * x[k];
            ycut  = (1-fracj) * y[k-1] + fracj * y[k];

            fac[ifac].xsrf[i+imax*j] = xform[0][0] * xcut + xform[1][0] * ycut + xform[2][0] * zcut;
            fac[ifac].ysrf[i+imax*j] = xform[0][1] * xcut + xform[1][1] * ycut + xform[2][1] * zcut;
            fac[ifac].zsrf[i+imax*j] = xform[0][2] * xcut + xform[1][2] * ycut + xform[2][2] * zcut;
        }

        FREE(y);
        FREE(x);
    }

    /* the i=imax-1 Grid plane is equispaced along the east Edges */
    i = imax - 1;
    for (j = 0; j < jmax; j++) {
        fac[ifac].xsrf[i+imax*j] = xeast;
        fac[ifac].ysrf[i+imax*j] = yeast;
        fac[ifac].zsrf[i+imax*j] = zeast;
    }

//cleanup:
    FREE(y);
    FREE(x);

    return status;
}

//$$$
//$$$/******************************************************************************/
//$$$
//$$$static void
//$$$makeSurface(int    ifac)                /* (in)  Face index (bias-0) */
//$$${
//$$$    typedef struct {
//$$$        int     beg;                         /* first Point in Loop */
//$$$        int     npnt;                        /* number of Vertex in Loop */
//$$$        double  alen;                        /* length (in XYZ) of Loop */
//$$$    } lup_t;
//$$$
//$$$    typedef struct {
//$$$        int     next;                        /* next boundary Vertex */
//$$$        int     lup;                         /* Loop number */
//$$$        double  ang;                         /* angle at Loop Vertices */
//$$$    } bnd_t;
//$$$
//$$$    int    ipnt, itri;
//$$$
//$$$    double *amat   = NULL;
//$$$    double *asmf   = NULL;
//$$$    int    *ismf   = NULL;
//$$$    double *urhs   = NULL;
//$$$    double *vrhs   = NULL;
//$$$    tri_T  *newTri = NULL;
//$$$    pnt_T  *newPnt = NULL;
//$$$    int    *ihole  = NULL;
//$$$    double *ahole  = NULL;
//$$$
//$$$    lup_t  *lups = NULL;                /* array  of Loops */
//$$$    int    nlup;                        /* number of Loops */
//$$$    bnd_t  *bnds = NULL;                /* array of boundary Vertices */
//$$$
//$$$    double ang0, ang1, ang2, ang3, ang4, ang5;
//$$$    double area, amin;
//$$$    double d1sq, d2sq, d3sq;
//$$$    double d01sq, d12sq, d20sq, d23sq;
//$$$    double d30sq, d04sq, d41sq, d15sq, d52sq;
//$$$    double du, dv;
//$$$    double dist;
//$$$    double err=1, erru, errv;
//$$$    int    found;
//$$$    int    i, j, k, ii, ij, nn, im1;
//$$$    int    imin;
//$$$    int    iter, itmax=1000;
//$$$    int    oldNtri, tmp;
//$$$    int    ip0, ip1, ip2, ip3, ip4, ip5;
//$$$    int    lup;
//$$$    int    oldIpnt;
//$$$    int    newNtri, count;
//$$$    int    nneg;
//$$$    int    outer;
//$$$    double sum;
//$$$    double xold, yold, zold, told, gold, hold;
//$$$    double xnew, ynew, znew, tnew, unew, vnew;
//$$$    int    nhole;
//$$$    int    MAXLOOPS = 500;
//$$$
//$$$    double errtol = EPS06;
//$$$    double frac = 0.25;
//$$$    double omega = 1.6;
//$$$
//$$$    ROUTINE(makeSurface);
//$$$
//$$$    /*----------------------------------------------------------------*/
//$$$
//$$$    printf("\n\nmakeSurface(ifac=%d)\n", ifac);
//$$$
//$$$    /* allocate storage that will be used to keep track of the Loop
//$$$       number and the next Vertex in the Loop */
//$$$    MALLOC(bnds, bnd_t, fac[ifac].npnt);
//$$$    MALLOC(lups, lup_t, MAXLOOPS      );
//$$$
//$$$    /* for each Point that is along a boundary of the Quilt, keep
//$$$       track of the next Point on the boundary */
//$$$    printf("finding next Point\n");
//$$$    for (ipnt = 0; ipnt < fac[ifac].npnt; ipnt++) {
//$$$        bnds[ipnt].next = -1;
//$$$        bnds[ipnt].lup  = -1;
//$$$    }
//$$$
//$$$    nn = 0;
//$$$    for (itri = 0; itri < fac[ifac].ntri; itri++) {
//$$$        ip0 = fac[ifac].tri[itri].p[0];
//$$$        ip1 = fac[ifac].tri[itri].p[1];
//$$$        ip2 = fac[ifac].tri[itri].p[2];
//$$$
//$$$        if (fac[ifac].tri[itri].t[0] < 0) {
//$$$            nn++;
//$$$            if (bnds[ip1].next == -1) {
//$$$                bnds[ip1].next = ip2;
//$$$            }
//$$$        }
//$$$        if (fac[ifac].tri[itri].t[1] < 0) {
//$$$            nn++;
//$$$            if (bnds[ip2].next == -1) {
//$$$                bnds[ip2].next = ip0;
//$$$            }
//$$$        }
//$$$        if (fac[ifac].tri[itri].t[2] < 0) {
//$$$            nn++;
//$$$            if (bnds[ip0].next == -1) {
//$$$                bnds[ip0].next = ip1;
//$$$            }
//$$$        }
//$$$    }
//$$$    printf("nn=%d\n", nn);
//$$$
//$$$    /* if there are no boundary Triangles, then return without
//$$$       setting UVtype */
//$$$    if (nn == 0) {
//$$$        printf("ERR: no boundary Tiangles were found\a\n");
//$$$        goto cleanup;
//$$$    }
//$$$
//$$$    /* associate each boundary Vertex with a Loop */
//$$$    printf("associating boundary Vertices with Loops\n");
//$$$    nlup = 0;
//$$$    while (1) {
//$$$
//$$$        /* find a Point that has a next but not a Loop */
//$$$        found = 0;
//$$$        for (ipnt = 0; ipnt < fac[ifac].npnt; ipnt++) {
//$$$            if ((bnds[ipnt].next >= 0) && (bnds[ipnt].lup < 0)) {
//$$$                bnds[ipnt].lup = nlup;
//$$$                bnds[ipnt].ang = 0;
//$$$
//$$$                lups[nlup].beg  = ipnt;
//$$$                lups[nlup].npnt = 1;
//$$$                lups[nlup].alen = 0;
//$$$
//$$$                found++;
//$$$                break;
//$$$            }
//$$$        }
//$$$
//$$$        /* if no unassigned boundary Point is found, then exit */
//$$$        if (found == 0) {
//$$$            break;
//$$$        }
//$$$
//$$$        /* assign all the Points in the Loop */
//$$$        count   = 0;
//$$$        oldIpnt = ipnt;
//$$$        ipnt    = bnds[ipnt].next;
//$$$        while (ipnt != lups[nlup].beg) {
//$$$            bnds[ipnt].lup = nlup;
//$$$            bnds[ipnt].ang = 0;
//$$$
//$$$            lups[nlup].npnt++;
//$$$            lups[nlup].alen += sqrt(  SQR(fac[ifac].xyz[3*ipnt  ] - fac[ifac].xyz[3*oldIpnt  ])
//$$$                                    + SQR(fac[ifac].xyz[3*ipnt+1] - fac[ifac].xyz[3*oldIpnt+1])
//$$$                                    + SQR(fac[ifac].xyz[3*ipnt+2] - fac[ifac].xyz[3*oldIpnt+2]));
//$$$
//$$$            oldIpnt = ipnt;
//$$$            ipnt    = bnds[ipnt].next;
//$$$            count++;
//$$$            if (count > fac[ifac].npnt+1) {
//$$$                printf("ERR: cannot form loop\a\n");
//$$$                goto cleanup;
//$$$            }
//$$$        }
//$$$
//$$$        printf("lups[%d].beg=%6d   .npnt=%6d   .alen=%f\n",
//$$$                nlup, lups[nlup].beg, lups[nlup].npnt, lups[nlup].alen);
//$$$
//$$$        nlup++;
//$$$    }
//$$$
//$$$    /* the outer Loop is the Loop with the greatest length (in XYZ) */
//$$$    printf("selecting outer Loop\n");
//$$$    outer = 0;
//$$$    for (lup = 1; lup < nlup; lup++) {
//$$$        if (lups[lup].alen > lups[outer].alen) {
//$$$            outer = lup;
//$$$        }
//$$$    }
//$$$    printf("outer Loop is %d\n", outer);
//$$$
//$$$    /* count the number of new Triangles we will have when we
//$$$       fill the inner Loops */
//$$$    newNtri = fac[ifac].ntri;
//$$$    for (lup = 0; lup < nlup; lup++) {
//$$$        if (lup != outer) {
//$$$            newNtri += lups[lup].npnt - 2;
//$$$        }
//$$$    }
//$$$    printf("newNtri=%d\n", newNtri);
//$$$
//$$$    /* get new tri array that is big enough for the new Triangles */
//$$$    MALLOC(newTri, tri_T, newNtri);
//$$$    memcpy(newTri, fac[ifac].tri, fac[ifac].ntri*sizeof(tri_T));
//$$$
//$$$    /* fill the inner Loops with Triangles */
//$$$    printf("filling inner Loops with Triangles\n");
//$$$    oldNtri = fac[ifac].ntri;
//$$$    if (nlup > 1) {
//$$$
//$$$        /* allocate arrays for closing holes */
//$$$        MALLOC(ihole, int,    fac[ifac].npnt);
//$$$        MALLOC(ahole, double, fac[ifac].npnt);
//$$$
//$$$        /* fill the inner Loops */
//$$$        for (lup = 0; lup < nlup; lup++) {
//$$$            if (lup == outer) continue;
//$$$
//$$$            /* make a list of the hole Vertices (ie, the Vertices
//$$$               that make up the Loop */
//$$$            nhole    = lups[lup].npnt;
//$$$            ihole[0] = lups[lup].beg;
//$$$
//$$$            for (i = 1; i < nhole; i++) {
//$$$                ihole[i] = bnds[ihole[i-1]].next;
//$$$            }
//$$$
//$$$            /* one-by-one, cut off the hole Vertex with the smallest angle */
//$$$            while (nhole > 3) {
//$$$
//$$$                /* compute the angle at each hole Vertex */
//$$$                for (i = 0; i < nhole; i++) {
//$$$                    im1 = i - 1;
//$$$                    ip1 = i + 1;
//$$$                    if (im1 <  0    ) im1 += nhole;
//$$$                    if (ip1 == nhole) ip1 -= nhole;
//$$$
//$$$                    d1sq = SQR(fac[ifac].xyz[3*ihole[i  ]  ] - fac[ifac].xyz[3*ihole[im1]  ])
//$$$                         + SQR(fac[ifac].xyz[3*ihole[i  ]+1] - fac[ifac].xyz[3*ihole[im1]+1])
//$$$                         + SQR(fac[ifac].xyz[3*ihole[i  ]+2] - fac[ifac].xyz[3*ihole[im1]+2]);
//$$$                    d2sq = SQR(fac[ifac].xyz[3*ihole[ip1]  ] - fac[ifac].xyz[3*ihole[i  ]  ])
//$$$                         + SQR(fac[ifac].xyz[3*ihole[ip1]+1] - fac[ifac].xyz[3*ihole[i  ]+1])
//$$$                         + SQR(fac[ifac].xyz[3*ihole[ip1]+2] - fac[ifac].xyz[3*ihole[i  ]+2]);
//$$$                    d3sq = SQR(fac[ifac].xyz[3*ihole[im1]  ] - fac[ifac].xyz[3*ihole[ip1]  ])
//$$$                         + SQR(fac[ifac].xyz[3*ihole[im1]+1] - fac[ifac].xyz[3*ihole[ip1]+1])
//$$$                         + SQR(fac[ifac].xyz[3*ihole[im1]+2] - fac[ifac].xyz[3*ihole[ip1]+1]);
//$$$
//$$$                    ahole[i] = ACOS((d1sq + d2sq - d3sq) / 2 / sqrt(d1sq * d2sq));
//$$$                }
//$$$
//$$$                /* minimum angle (which will be cut off) */
//$$$                amin = ahole[0];
//$$$                imin = 0;
//$$$
//$$$                for (i = 2; i < nhole; i++) {
//$$$                    if (ahole[i] < amin) {
//$$$                        amin = ahole[i];
//$$$                        imin = i;
//$$$                    }
//$$$                }
//$$$
//$$$                /* make new Triangle */
//$$$                im1 = imin - 1;
//$$$                ip1 = imin + 1;
//$$$                if (im1 <  0    ) im1 += nhole;
//$$$                if (ip1 == nhole) ip1 -= nhole;
//$$$
//$$$                newTri[newNtri].p[0] = ihole[imin];
//$$$                newTri[newNtri].p[1] = ihole[im1 ];
//$$$                newTri[newNtri].p[2] = ihole[ip1 ];
//$$$                newTri[newNtri].t[0] = -1;
//$$$                newTri[newNtri].t[1] = 1;
//$$$                newTri[newNtri].t[2] = 1;
//$$$                newNtri++;
//$$$
//$$$                /* remove imin from hole list */
//$$$                nhole--;
//$$$                for (j = imin; j < nhole; j++) {
//$$$                    ihole[j] = ihole[j+1];
//$$$                    ahole[j] = ahole[j+1];
//$$$                }
//$$$            }
//$$$
//$$$            /* make the final Triangle with the last 3 entries in
//$$$               the hole list */
//$$$            newTri[newNtri].p[0] = ihole[3];
//$$$            newTri[newNtri].p[1] = ihole[2];
//$$$            newTri[newNtri].p[2] = ihole[1];
//$$$            newTri[newNtri].t[0] = -1;
//$$$            newTri[newNtri].t[1] = -1;
//$$$            newTri[newNtri].t[2] = -1;
//$$$            newNtri++;
//$$$        }
//$$$
//$$$        /* set up the Triangle neighbor information for the new Triangles
//$$$           (note the backwards loop since it is most likely that the neighnor
//$$$           will be near the end of the list) */
//$$$        for (i = oldNtri; i < fac[ifac].ntri; i++) {
//$$$            ip0 = newTri[i].p[0];
//$$$            ip1 = newTri[i].p[1];
//$$$            ip2 = newTri[i].p[2];
//$$$
//$$$            for (j = fac[ifac].ntri-1; j >= 0; j--) {
//$$$                if        (newTri[j].p[2] == ip1 && newTri[j].p[1] == ip2) {
//$$$                    newTri[i].t[0] = j;
//$$$                    newTri[j].t[0] = i;
//$$$                } else if (newTri[j].p[0] == ip1 && newTri[j].p[2] == ip2) {
//$$$                    newTri[i].t[0] = j;
//$$$                    newTri[j].t[1] = i;
//$$$                } else if (newTri[j].p[1] == ip1 && newTri[j].p[0] == ip2) {
//$$$                    newTri[i].t[0] = j;
//$$$                    newTri[j].t[2] = i;
//$$$                }
//$$$            }
//$$$
//$$$            for (j = fac[ifac].ntri-1; j >= 0; j--) {
//$$$                if        (newTri[j].p[2] == ip2 && newTri[j].p[1] == ip0) {
//$$$                    newTri[i].t[1] = j;
//$$$                    newTri[j].t[0] = i;
//$$$                } else if (newTri[j].p[0] == ip2 && newTri[j].p[2] == ip0) {
//$$$                    newTri[i].t[1] = j;
//$$$                    newTri[j].t[1] = i;
//$$$                } else if (newTri[j].p[1] == ip2 && newTri[j].p[0] == ip0) {
//$$$                    newTri[i].t[1] = j;
//$$$                    newTri[j].t[2] = i;
//$$$                }
//$$$            }
//$$$
//$$$            for (j = fac[ifac].ntri-1; j >= 0; j--) {
//$$$                if        (newTri[j].p[2] == ip0 && newTri[j].p[1] == ip1) {
//$$$                    newTri[i].t[2] = j;
//$$$                    newTri[j].t[0] = i;
//$$$                } else if (newTri[j].p[0] == ip0 && newTri[j].p[2] == ip1) {
//$$$                    newTri[i].t[2] = j;
//$$$                    newTri[j].t[1] = i;
//$$$                } else if (newTri[j].p[1] == ip0 && newTri[j].p[0] == ip1) {
//$$$                    newTri[i].t[2] = j;
//$$$                    newTri[j].t[2] = i;
//$$$                }
//$$$            }
//$$$        }
//$$$    }
//$$$
//$$$    /* compute the angle at each boundary Vertex associated
//$$$       with the outer Loop */
//$$$    printf("computing boundary angles\n");
//$$$    for (itri = 0; itri < fac[ifac].ntri; itri++) {
//$$$        ip0 = newTri[itri].p[0];
//$$$        ip1 = newTri[itri].p[1];
//$$$        ip2 = newTri[itri].p[2];
//$$$
//$$$        if (bnds[ip0].lup >= 0 || bnds[ip1].lup >= 0 || bnds[ip2].lup >= 0) {
//$$$            d01sq = SQR(fac[ifac].xyz[3*ip0  ] - fac[ifac].xyz[3*ip1  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip0+1] - fac[ifac].xyz[3*ip1+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip0+2] - fac[ifac].xyz[3*ip1+2]);
//$$$            d12sq = SQR(fac[ifac].xyz[3*ip1  ] - fac[ifac].xyz[3*ip2  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip1+1] - fac[ifac].xyz[3*ip2+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip1+2] - fac[ifac].xyz[3*ip2+2]);
//$$$            d20sq = SQR(fac[ifac].xyz[3*ip2  ] - fac[ifac].xyz[3*ip0  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip2+1] - fac[ifac].xyz[3*ip0+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip2+2] - fac[ifac].xyz[3*ip0+2]);
//$$$
//$$$            ang0 = ACOS((d20sq + d01sq - d12sq) / 2 / sqrt(d20sq * d01sq));
//$$$            ang1 = ACOS((d01sq + d12sq - d20sq) / 2 / sqrt(d01sq * d12sq));
//$$$            ang2 = ACOS((d12sq + d20sq - d01sq) / 2 / sqrt(d12sq * d20sq));
//$$$
//$$$            if (bnds[ip0].lup >= 0) {bnds[ip0].ang += ang0;}
//$$$            if (bnds[ip1].lup >= 0) {bnds[ip1].ang += ang1;}
//$$$            if (bnds[ip2].lup >= 0) {bnds[ip2].ang += ang2;}
//$$$        }
//$$$    }
//$$$
//$$$    /* lay out the outer Loop (in UV) by marching along the Loop */
//$$$    printf("laying out outer Loop\n");
//$$$    ipnt                   = lups[outer].beg;
//$$$    fac[ifac].uv[2*ipnt  ] = 0;
//$$$    fac[ifac].uv[2*ipnt+1] = 0;
//$$$
//$$$    xold = fac[ifac].xyz[3*ipnt  ];
//$$$    yold = fac[ifac].xyz[3*ipnt+1];
//$$$    zold = fac[ifac].xyz[3*ipnt+2];
//$$$    gold = 0;
//$$$    hold = 0;
//$$$    told = 0;
//$$$
//$$$    do {
//$$$        ipnt = bnds[ipnt].next;
//$$$
//$$$        xnew = fac[ifac].xyz[3*ipnt  ];
//$$$        ynew = fac[ifac].xyz[3*ipnt+1];
//$$$        znew = fac[ifac].xyz[3*ipnt+2];
//$$$
//$$$        dist = sqrt(SQR(xnew - xold) + SQR(ynew - yold) + SQR(znew - zold));
//$$$
//$$$        unew = gold + dist * cos(told);
//$$$        vnew = hold + dist * sin(told);
//$$$        tnew = told + PI - bnds[ipnt].ang;
//$$$
//$$$        fac[ifac].uv[2*ipnt  ] = unew;
//$$$        fac[ifac].uv[2*ipnt+1] = vnew;
//$$$
//$$$        xold = xnew;
//$$$        yold = ynew;
//$$$        zold = znew;
//$$$        gold = unew;
//$$$        hold = vnew;
//$$$        told = tnew;
//$$$    } while (ipnt != lups[outer].beg);
//$$$
//$$$    {
//$$$        int    io_kbd = 5;
//$$$        int    io_scr = 6;
//$$$        int    indgr  = 1 + 2 + 4 + 16 + 64;
//$$$        int    nline  = 1;
//$$$        int    ilin[0], isym[1], nper[1];
//$$$        float  uplot[10000], vplot[10000];
//$$$        char   pltitl[80];
//$$$
//$$$        strcpy(pltitl, "~u~v~initial outer boundary");
//$$$
//$$$        ipnt = lups[outer].beg;
//$$$
//$$$        nper[0] =  0;
//$$$        ilin[0] = +1;
//$$$        isym[0] = -1;
//$$$        nline   =  1;
//$$$
//$$$        do {
//$$$            uplot[nper[0]] = (float)(fac[ifac].uv[2*ipnt  ]);
//$$$            vplot[nper[0]] = (float)(fac[ifac].uv[2*ipnt+1]);
//$$$            nper[0]++;
//$$$
//$$$            ipnt = bnds[ipnt].next;
//$$$        } while (ipnt != lups[outer].beg);
//$$$
//$$$        grinit_(&io_kbd, &io_scr, "makeSurface", strlen("makeSurface"));
//$$$        grline_(ilin, isym, &nline, pltitl, &indgr, uplot, vplot, nper, strlen(pltitl));
//$$$    }
//$$$
//$$$    /* determine how close we came to closing the Loop (in UV) */
//$$$    ipnt = lups[outer].beg;
//$$$    dist = sqrt(SQR(fac[ifac].uv[2*ipnt]) + SQR(fac[ifac].uv[2*ipnt+1]));
//$$$    printf("closing dist=%f\n", dist);
//$$$
//$$$    /* if the Loop is almost closed, close it by linearly
//$$$       adjusting the UV */
//$$$    printf("closing outer Loop\n");
//$$$    if (dist < frac*lups[outer].alen) {
//$$$        du = fac[ifac].uv[2*ipnt  ] / lups[outer].npnt;
//$$$        dv = fac[ifac].uv[2*ipnt+1] / lups[outer].npnt;
//$$$        printf("adjusting by du=%f,  dv=%f\n", du, dv);
//$$$
//$$$        ii = 0;
//$$$        do {
//$$$            ipnt = bnds[ipnt].next;
//$$$            ii++;
//$$$
//$$$            fac[ifac].uv[2*ipnt  ] -= (ii * du);
//$$$            fac[ifac].uv[2*ipnt+1] -= (ii * dv);
//$$$        } while (ipnt != lups[outer].beg);
//$$$
//$$$    /* otherwise, remap the outer Loop to the unit circle */
//$$$    } else {
//$$$        printf("mapping to circle\n");
//$$$        ipnt = lups[outer].beg;
//$$$
//$$$        xold = fac[ifac].xyz[3*ipnt  ];
//$$$        yold = fac[ifac].xyz[3*ipnt+1];
//$$$        zold = fac[ifac].xyz[3*ipnt+2];
//$$$        told = 0.0;
//$$$
//$$$        fac[ifac].uv[2*ipnt  ] = cos(told);
//$$$        fac[ifac].uv[2*ipnt+1] = sin(told);
//$$$
//$$$        ipnt = bnds[ipnt].next;
//$$$        while (ipnt != lups[outer].beg) {
//$$$            xnew = fac[ifac].xyz[3*ipnt  ];
//$$$            ynew = fac[ifac].xyz[3*ipnt+1];
//$$$            znew = fac[ifac].xyz[3*ipnt+2];
//$$$
//$$$            dist = sqrt(SQR(xnew - xold) + SQR(ynew - yold)
//$$$                                         + SQR(znew - zold));
//$$$            told = told + TWOPI * dist / lups[outer].alen;
//$$$
//$$$            fac[ifac].uv[2*ipnt  ] = cos(told);
//$$$            fac[ifac].uv[2*ipnt+1] = sin(told);
//$$$
//$$$            ipnt = bnds[ipnt].next;
//$$$            xold = xnew;
//$$$            yold = ynew;
//$$$            zold = znew;
//$$$        }
//$$$    }
//$$$
//$$$    {
//$$$        int    io_kbd = 5;
//$$$        int    io_scr = 6;
//$$$        int    indgr  = 1 + 2 + 4 + 16 + 64;
//$$$        int    nline  = 1;
//$$$        int    ilin[0], isym[1], nper[1];
//$$$        float  uplot[10000], vplot[10000];
//$$$        char   pltitl[80];
//$$$
//$$$        strcpy(pltitl, "~u~v~outer boundary");
//$$$
//$$$        ipnt = lups[outer].beg;
//$$$
//$$$        nper[0] =  0;
//$$$        ilin[0] = +1;
//$$$        isym[0] = -1;
//$$$        nline   =  1;
//$$$
//$$$        do {
//$$$            uplot[nper[0]] = (float)(fac[ifac].uv[2*ipnt  ]);
//$$$            vplot[nper[0]] = (float)(fac[ifac].uv[2*ipnt+1]);
//$$$            nper[0]++;
//$$$
//$$$            ipnt = bnds[ipnt].next;
//$$$        } while (ipnt != lups[outer].beg);
//$$$
//$$$        grinit_(&io_kbd, &io_scr, "makeSurface", strlen("makeSurface"));
//$$$        grline_(ilin, isym, &nline, pltitl, &indgr, uplot, vplot, nper, strlen(pltitl));
//$$$    }
//$$$return;
//$$$
//$$$    /* initialize the amat matrix (which will be used to solve for the UV
//$$$       at all interior Vertices) */
//$$$    MALLOC(amat, double, fac[ifac].npnt*fac[ifac].npnt);
//$$$
//$$$    for (ij = 0; ij < fac[ifac].npnt*fac[ifac].npnt; ij++) {
//$$$        amat[ij] = 0;
//$$$    }
//$$$
//$$$    /* find the mean value weights at each Vertex due to its neighbors.  the
//$$$       weigt function is the "mean value coordinates" function as described
//$$$       by M. Floater
//$$$
//$$$                        ip3--------ip2---------ip5
//$$$                          \         /\a5       /
//$$$                           \       /a2\       /
//$$$                            \     /    \     /
//$$$                             \a3 / itri \   /
//$$$                              \ /a0    a1\ /
//$$$                              ip0--------ip1
//$$$                                \     a4 /
//$$$                                 \      /
//$$$                                  \    /
//$$$                                   \  /
//$$$                                    \/
//$$$                                   ip4
//$$$     */
//$$$    printf("finding mean value weights\n");
//$$$    for (itri = 0; itri < fac[ifac].ntri; itri++) {
//$$$        ip0 = newTri[itri].p[0];
//$$$        ip1 = newTri[itri].p[1];
//$$$        ip2 = newTri[itri].p[2];
//$$$
//$$$        d01sq = SQR(fac[ifac].xyz[3*ip0  ] - fac[ifac].xyz[3*ip1  ])
//$$$              + SQR(fac[ifac].xyz[3*ip0+1] - fac[ifac].xyz[3*ip1+1])
//$$$              + SQR(fac[ifac].xyz[3*ip0+2] - fac[ifac].xyz[3*ip1+2]);
//$$$        d12sq = SQR(fac[ifac].xyz[3*ip1  ] - fac[ifac].xyz[3*ip2  ])
//$$$              + SQR(fac[ifac].xyz[3*ip1+1] - fac[ifac].xyz[3*ip2+1])
//$$$              + SQR(fac[ifac].xyz[3*ip1+2] - fac[ifac].xyz[3*ip2+2]);
//$$$        d20sq = SQR(fac[ifac].xyz[3*ip2  ] - fac[ifac].xyz[3*ip0  ])
//$$$              + SQR(fac[ifac].xyz[3*ip2+1] - fac[ifac].xyz[3*ip0+1])
//$$$              + SQR(fac[ifac].xyz[3*ip2+2] - fac[ifac].xyz[3*ip0+2]);
//$$$
//$$$        ang0 = ACOS((d20sq + d01sq - d12sq) / 2 / sqrt(d20sq * d01sq));
//$$$        ang1 = ACOS((d01sq + d12sq - d20sq) / 2 / sqrt(d01sq * d12sq));
//$$$        ang2 = ACOS((d12sq + d20sq - d01sq) / 2 / sqrt(d12sq * d20sq));
//$$$
//$$$        if (newTri[itri].t[1] > 0) {
//$$$            tmp = newTri[itri].t[1];
//$$$            if (       newTri[tmp].t[0] == itri) {
//$$$                ip3 =  newTri[tmp].p[0];
//$$$            } else if (newTri[tmp].t[1] == itri) {
//$$$                ip3 =  newTri[tmp].p[1];
//$$$            } else {
//$$$                ip3 =  newTri[tmp].p[2];
//$$$            }
//$$$
//$$$            d23sq = SQR(fac[ifac].xyz[3*ip2  ] - fac[ifac].xyz[3*ip3  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip2+1] - fac[ifac].xyz[3*ip3+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip2+2] - fac[ifac].xyz[3*ip3+2]);
//$$$            d30sq = SQR(fac[ifac].xyz[3*ip3  ] - fac[ifac].xyz[3*ip0  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip3+1] - fac[ifac].xyz[3*ip0+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip3+2] - fac[ifac].xyz[3*ip0+2]);
//$$$
//$$$            ang3 = ACOS((d30sq + d20sq - d23sq) / 2 / sqrt(d30sq * d20sq));
//$$$
//$$$            amat[ip0+ip2*fac[ifac].npnt] = (tan(ang0/2) + tan(ang3/2)) / sqrt(d20sq);
//$$$        }
//$$$
//$$$        if (newTri[itri].t[2] > 0) {
//$$$            tmp = newTri[itri].t[2];
//$$$            if (       newTri[tmp].t[0] == itri) {
//$$$                ip4 =  newTri[tmp].p[0];
//$$$            } else if (newTri[tmp].t[1] == itri) {
//$$$                ip4 =  newTri[tmp].p[1];
//$$$            } else {
//$$$                ip4 =  newTri[tmp].p[2];
//$$$            }
//$$$
//$$$            d04sq = SQR(fac[ifac].xyz[3*ip0  ] - fac[ifac].xyz[3*ip4  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip0+1] - fac[ifac].xyz[3*ip4+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip0+2] - fac[ifac].xyz[3*ip4+2]);
//$$$            d41sq = SQR(fac[ifac].xyz[3*ip4  ] - fac[ifac].xyz[3*ip1  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip4+1] - fac[ifac].xyz[3*ip1+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip4+2] - fac[ifac].xyz[3*ip1+2]);
//$$$
//$$$            ang4 = ACOS((d41sq + d01sq - d04sq) / 2 / sqrt(d41sq * d01sq));
//$$$
//$$$            amat[ip1+ip0*fac[ifac].npnt] = (tan(ang1/2) + tan(ang4/2)) / sqrt(d01sq);
//$$$        }
//$$$
//$$$        if (newTri[itri].t[0] > 0) {
//$$$            tmp = newTri[itri].t[0];
//$$$            if (       newTri[tmp].t[0] == itri) {
//$$$                ip5 =  newTri[tmp].p[0];
//$$$            } else if (newTri[tmp].t[1] == itri) {
//$$$                ip5 =  newTri[tmp].p[1];
//$$$            } else {
//$$$                ip5 =  newTri[tmp].p[2];
//$$$            }
//$$$
//$$$            d15sq = SQR(fac[ifac].xyz[3*ip1  ] - fac[ifac].xyz[3*ip5  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip1+1] - fac[ifac].xyz[3*ip5+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip1+2] - fac[ifac].xyz[3*ip5+2]);
//$$$            d52sq = SQR(fac[ifac].xyz[3*ip5  ] - fac[ifac].xyz[3*ip2  ])
//$$$                  + SQR(fac[ifac].xyz[3*ip5+1] - fac[ifac].xyz[3*ip2+1])
//$$$                  + SQR(fac[ifac].xyz[3*ip5+2] - fac[ifac].xyz[3*ip2+2]);
//$$$
//$$$            ang5 = ACOS((d52sq + d12sq - d15sq) / 2 / sqrt(d52sq * d12sq));
//$$$
//$$$            amat[ip2+ip1*fac[ifac].npnt] = (tan(ang2/2) + tan(ang5/2)) / sqrt(d12sq);
//$$$        }
//$$$    }
//$$$
//$$$    /* set up the final amat and the right-hand sides */
//$$$    MALLOC(urhs, double, fac[ifac].npnt);
//$$$    MALLOC(vrhs, double, fac[ifac].npnt);
//$$$
//$$$    for (i = 0; i < fac[ifac].npnt; i++) {
//$$$
//$$$        /* for interior Vertices, normalize the weights, set the
//$$$           diagonal to 1, and zero-out the right-hand sides */
//$$$        if (bnds[i].lup != outer) {
//$$$            sum = 0;
//$$$            for (j = 0; j < fac[ifac].npnt; j++) {
//$$$                sum -= amat[i+j*fac[ifac].npnt];
//$$$            }
//$$$            for (j = 0; j < fac[ifac].npnt; j++) {
//$$$                amat[i+j*fac[ifac].npnt] /= sum;
//$$$            }
//$$$
//$$$            amat[i+i*fac[ifac].npnt] = 1;
//$$$
//$$$            urhs[i] = 0;
//$$$            vrhs[i] = 0;
//$$$
//$$$            fac[ifac].uv[2*i  ] = 0;
//$$$            fac[ifac].uv[2*i+1] = 0;
//$$$
//$$$        /* for boundary Vertices, zero-out all amat elements, set the diagonal
//$$$           to 1, and store the boundary values in the RHS */
//$$$        } else {
//$$$            for (j = 0; j < fac[ifac].npnt; j++) {
//$$$                amat[i+j*fac[ifac].npnt] = 0;
//$$$            }
//$$$
//$$$            amat[i+i*fac[ifac].npnt] = 1;
//$$$
//$$$            urhs[i] = fac[ifac].uv[2*i  ];
//$$$            vrhs[i] = fac[ifac].uv[2*i+1];
//$$$        }
//$$$    }
//$$$
//$$$    /* count the number of non-zero entries in amat */
//$$$    nn = 0;
//$$$    for (i = 0; i < fac[ifac].npnt; i++) {
//$$$        for (j = 0; j < fac[ifac].npnt; j++) {
//$$$            if (fabs(amat[i+j*fac[ifac].npnt]) < EPS20) {
//$$$                nn++;
//$$$            }
//$$$        }
//$$$    }
//$$$
//$$$    /* allocate matrices for sparse-matrix form.  this form is as
//$$$       described in 'Numerical Recipes' */
//$$$    printf("setting up sparse matrix\n");
//$$$    MALLOC(asmf, double, nn);
//$$$    MALLOC(ismf, int,    nn);
//$$$
//$$$    /* store the matrix in sparse-matrix form (which will make
//$$$       the iterations below much more efficient) */
//$$$    for (j = 0; j < fac[ifac].npnt; j++) {
//$$$        asmf[j] = amat[j+j*fac[ifac].npnt];
//$$$    }
//$$$
//$$$    ismf[0] = fac[ifac].npnt + 2;
//$$$
//$$$    k = fac[ifac].npnt + 1;
//$$$    for (i = 0; i < fac[ifac].npnt; i++) {
//$$$        for (j = 0; j < fac[ifac].npnt; j++) {
//$$$            if (i!=j && fabs(amat[i+j*fac[ifac].npnt]) > EPS20) {
//$$$                k++;
//$$$                asmf[k] = amat[i+j*fac[ifac].npnt];
//$$$                ismf[k] = j;
//$$$            }
//$$$        }
//$$$
//$$$        ismf[i+1] = k + 1;
//$$$    }
//$$$
//$$$    /* solve for the Gs and Hs using successive-over-relaxation */
//$$$    printf("solving sparse matrix\n");
//$$$    for (iter = 0; iter < itmax; iter++) {
//$$$
//$$$        /* apply successive-over-relaxation */
//$$$        for (i = 0; i < fac[ifac].npnt; i++) {
//$$$            du = urhs[i] - asmf[i] * fac[ifac].uv[2*i  ];
//$$$            dv = vrhs[i] - asmf[i] * fac[ifac].uv[2*i+1];
//$$$
//$$$            for (k = ismf[i]; k < ismf[i+1]; k++) {
//$$$                j    = ismf[k];
//$$$                du  -= asmf[k] * fac[ifac].uv[2*j  ];
//$$$                dv  -= asmf[k] * fac[ifac].uv[2*j+1];
//$$$            }
//$$$
//$$$            fac[ifac].uv[2*i  ] += omega * du / asmf[i];
//$$$            fac[ifac].uv[2*i+1] += omega * dv / asmf[i];
//$$$        }
//$$$
//$$$        /* compute the norm of the residual */
//$$$        err = 0;
//$$$        for (i = 0; i < fac[ifac].npnt; i++) {
//$$$            erru = urhs[i] - asmf[i] * fac[ifac].uv[2*i  ];
//$$$            errv = vrhs[i] - asmf[i] * fac[ifac].uv[2*i+1];
//$$$
//$$$            for (k = ismf[i]; k < ismf[i+1]; k++) {
//$$$                j     = ismf[k];
//$$$                erru -= asmf[k] * fac[ifac].uv[2*j  ];
//$$$                errv -= asmf[k] * fac[ifac].uv[2*j+1];
//$$$            }
//$$$
//$$$            err += SQR(erru) + SQR(errv);
//$$$        }
//$$$        err = sqrt(err);
//$$$
//$$$        /* exit if converged */
//$$$        printf("iter=%5d   err=%15.8e\n", iter, err);
//$$$        if (err < errtol) break;
//$$$    }
//$$$
//$$$    /* find the number of Triangles with negative areas (in UV) */
//$$$    printf("checking final Triangles\n");
//$$$    nneg = 0;
//$$$    for (itri = 0; itri < oldNtri; itri++) {
//$$$        ip0 = newTri[itri].p[0];
//$$$        ip1 = newTri[itri].p[1];
//$$$        ip2 = newTri[itri].p[2];
//$$$
//$$$        area = (fac[ifac].uv[2*ip1  ] - fac[ifac].uv[2*ip0  ]) * (fac[ifac].uv[2*ip2+1] - fac[ifac].uv[2*ip0+1])
//$$$             - (fac[ifac].uv[2*ip1+1] - fac[ifac].uv[2*ip0+1]) * (fac[ifac].uv[2*ip2  ] - fac[ifac].uv[2*ip0  ]);
//$$$
//$$$        if (area < 0) nneg++;
//$$$    }
//$$$    printf("floater parameterization has nneg=%d of %d\n", nneg, fac[ifac].ntri);
//$$$
//$$$    /* if we converged and none of the Triangles have negative areas,
//$$$       then set the UVtype */
//$$$    if (err >= errtol || nneg != 0) {
//$$$        printf("ERR: floater did not succeed\a\n");
//$$$    }
//$$$
//$$$ cleanup:
//$$$    FREE(bnds  );
//$$$    FREE(lups  );
//$$$    FREE(urhs  );
//$$$    FREE(vrhs  );
//$$$    FREE(amat  );
//$$$    FREE(asmf  );
//$$$    FREE(ismf  );
//$$$    FREE(newTri);
//$$$    FREE(newPnt);
//$$$    FREE(ihole );
//$$$    FREE(ahole );
//$$$    return;
//$$$}


/******************************************************************************/

static int
makeTopology()
{
    int    status = SUCCESS;            /* (out) return status */

    int    ipnt, itri, jtri, itri1, icol1, icol2, isid, i, ilup, imid;
    int    inod, iedg, ifac, ibeg, iend, *pntNod = NULL;
    int    itemp[10000], ntemp;

    ROUTINE(makeTopology);

    /*----------------------------------------------------------------*/

    /* remove prior topologies */
    FREE(nod);
    Nnod = 0;

    for (iedg = 0; iedg < Nedg; iedg++) {
        FREE(edg[iedg].pnt);
    }
    FREE(edg);
    Nedg = 0;

    for (ifac = 0; ifac < Nfac; ifac++) {
        status = initialTess(&(fac[ifac].tess));
        CHECK_STATUS(initialTess);

        FREE(fac[ifac].edg );

        FREE(fac[ifac].xsrf);
        FREE(fac[ifac].ysrf);
        FREE(fac[ifac].zsrf);
    }
    FREE(fac);
    Nfac = 0;

    ntemp = 0;

    /* allocate an array in which we will keep track of the Nodes */
    MALLOC(pntNod, int, tess.npnt);

    for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
        pntNod[ipnt] = -1;
    }

    /* loop until there are no Sides that have Triangles with different colors */
    while (1) {

        /* look for a Side that has Triangles of different colors */
        itri1 = -1;
        icol1 = -1;
        icol2 = -1;

        for (itri = 0; itri < tess.ntri; itri++) {
            for (isid = 0; isid < 3; isid++) {
                if (isid == 0 && (tess.ttyp[itri] & TRI_T0_EDGE)) continue;
                if (isid == 1 && (tess.ttyp[itri] & TRI_T1_EDGE)) continue;
                if (isid == 2 && (tess.ttyp[itri] & TRI_T2_EDGE)) continue;

                jtri = tess.trit[3*itri+isid];
                if (jtri >= 0) {
                    if ((tess.ttyp[itri] & TRI_COLOR) != (tess.ttyp[jtri] & TRI_COLOR)) {
                        itri1 = itri;
                        icol1 = tess.ttyp[itri] & TRI_COLOR;
                        icol2 = tess.ttyp[jtri] & TRI_COLOR;

                        itemp[0] = tess.trip[3*itri1+(isid+1)%3];
                        itemp[1] = tess.trip[3*itri1+(isid+2)%3];
                        ntemp = 2;

                        if (isid == 0) tess.ttyp[itri] |= TRI_T0_EDGE;
                        if (isid == 1) tess.ttyp[itri] |= TRI_T1_EDGE;
                        if (isid == 2) tess.ttyp[itri] |= TRI_T2_EDGE;

                        if (tess.trit[3*jtri  ] == itri) tess.ttyp[jtri] |= TRI_T0_EDGE;
                        if (tess.trit[3*jtri+1] == itri) tess.ttyp[jtri] |= TRI_T1_EDGE;
                        if (tess.trit[3*jtri+2] == itri) tess.ttyp[jtri] |= TRI_T2_EDGE;

                        break;
                    }
                }
            }
            if (icol1 >= 0) break;
        }

        /* if no Triangles with different colors were found, we are done */
        if (icol1 < 0) break;

        /* propagate at the beginning of the current Edge until we have a new Node */
        ilup =  0;
        inod = -1;
        jtri = itri1;
        while (inod < 0) {

            /* check for a loop */
            if (itemp[0] == itemp[ntemp-1]) {
                addNode(itemp[0], pntNod);
                ilup = 1;
                break;
            }

            /* look through all Sides of jtri */
            for (isid = 0; isid < 3; isid++) {
                if (tess.trip[3*jtri+isid] != itemp[0]) continue;

                /* if the same as icol1, keep spinning around the Point */
                if ((tess.ttyp[jtri] & TRI_COLOR) == icol1) {
                    jtri = tess.trit[3*jtri+(isid+1)%3];

                    /* we have reached a non-manifold boundary */
                    if (jtri < 0) {
                        addNode(itemp[0], pntNod);
                        inod = Nnod;
                    }

                /* if we match icol2, then extend the Edge */
                } else if ((tess.ttyp[jtri] & TRI_COLOR) == icol2) {
                    for (i = ntemp; i > 0; i--) {
                        itemp[i] = itemp[i-1];
                    }
                    itemp[0] = tess.trip[3*jtri+(isid+1)%3];
                    ntemp++;
                    if (ntemp >= 9999) {
                        printf("ERR: recompile with bigger itemp\a\n");
                        goto cleanup;
                    }

                    if (isid == 1) tess.ttyp[jtri] |= TRI_T0_EDGE;
                    if (isid == 2) tess.ttyp[jtri] |= TRI_T1_EDGE;
                    if (isid == 0) tess.ttyp[jtri] |= TRI_T2_EDGE;

                    jtri = tess.trit[3*jtri+(isid+2)%3];

                    if (tess.trip[3*jtri+1] == itemp[0] &&
                        tess.trip[3*jtri+2] == itemp[1]   ) tess.ttyp[jtri] |= TRI_T0_EDGE;
                    if (tess.trip[3*jtri+2] == itemp[0] &&
                        tess.trip[3*jtri  ] == itemp[1]   ) tess.ttyp[jtri] |= TRI_T1_EDGE;
                    if (tess.trip[3*jtri  ] == itemp[0] &&
                        tess.trip[3*jtri+1] == itemp[1]   ) tess.ttyp[jtri] |= TRI_T2_EDGE;

                /* some other color means we need to add a Node */
                } else {
                    addNode(itemp[0], pntNod);
                    inod = Nnod;
                }

                break;
            }
        }

        /* propagate at the end of the current Edge until we have a new Node */
        inod = -1;
        jtri = itri1;
        while (inod < 0 && ilup == 0) {

            /* look through all Sides of jtri */
            for (isid = 0; isid < 3; isid++) {
                if (tess.trip[3*jtri+isid] != itemp[ntemp-1]) continue;

                /* if the same as icol1, keep spinning around the Point */
                if ((tess.ttyp[jtri]& TRI_COLOR) == icol1) {
                    jtri = tess.trit[3*jtri+(isid+2)%3];

                    /* we have reached a non-manifold boundary */
                    if (jtri < 0) {
                        addNode(itemp[ntemp-1], pntNod);
                        inod = Nnod;
                    }

                /* if we match icol2, then extend the Edge */
                } else if ((tess.ttyp[jtri] & TRI_COLOR) == icol2) {
                    itemp[ntemp++] = tess.trip[3*jtri+(isid+2)%3];
                    if (ntemp >= 9999) goto cleanup;

                    if (isid == 2) tess.ttyp[jtri] |= TRI_T0_EDGE;
                    if (isid == 0) tess.ttyp[jtri] |= TRI_T1_EDGE;
                    if (isid == 1) tess.ttyp[jtri] |= TRI_T2_EDGE;

                    jtri = tess.trit[3*jtri+(isid+1)%3];

                    if (tess.trip[3*jtri+1] == itemp[ntemp-2] &&
                        tess.trip[3*jtri+2] == itemp[ntemp-1]   ) tess.ttyp[jtri] |= TRI_T0_EDGE;
                    if (tess.trip[3*jtri+2] == itemp[ntemp-2] &&
                        tess.trip[3*jtri  ] == itemp[ntemp-1]   ) tess.ttyp[jtri] |= TRI_T1_EDGE;
                    if (tess.trip[3*jtri  ] == itemp[ntemp-2] &&
                        tess.trip[3*jtri+1] == itemp[ntemp-1]   ) tess.ttyp[jtri] |= TRI_T2_EDGE;

                    /* some other color means we need to add a Node */
                } else {
                    addNode(itemp[ntemp-1], pntNod);
                    inod = Nnod;
                }

                break;
            }
        }

        /* if not periodic, make a new Edge */
        if (itemp[0] != itemp[ntemp-1]) {
            ibeg = pntNod[itemp[0      ]];
            iend = pntNod[itemp[ntemp-1]];

            status = addEdge(icol1, icol2, ntemp, itemp, ibeg, iend);
            CHECK_STATUS(addEdge);

        /* otherwise, add a Node at the beginning and make two Edges */
        } else {
            imid = ntemp / 2;

            addNode(itemp[imid], pntNod);
            ibeg = pntNod[itemp[0   ]];
            iend = pntNod[itemp[imid]];

            status = addEdge(icol1, icol2, imid+1,       itemp       , ibeg, iend);
            CHECK_STATUS(addEdge);

            status = addEdge(icol1, icol2, ntemp-imid, &(itemp[imid]), iend, ibeg);
            CHECK_STATUS(addEdge);
        }
    }

    /* count the number of Triangles associated with each color */
    for (i = 0; i <= tess.ncolr; i++) {
        itemp[i] = 0;
    }

    for (itri = 0; itri < tess.ntri; itri++) {
        (itemp[tess.ttyp[itri] & TRI_COLOR])++;
    }

    /* make a Face for each used color */
    for (i = 0; i <= tess.ncolr; i++) {
        if (itemp[i] > 0) {
            status = addFace(i);
            CHECK_STATUS(addFace);
        }
    }

    /* remember the Face associated with each color */
    for (ifac = 0; ifac < Nfac; ifac++) {
        itemp[fac[ifac].icol] = ifac;
    }

    /* adjust the Edge-to-Face pointers */
    for (iedg = 0; iedg < Nedg; iedg++) {
        edg[iedg].ileft = itemp[edg[iedg].ileft];
        edg[iedg].irite = itemp[edg[iedg].irite];
    }

cleanup:
    FREE(pntNod);

    return status;
}


/******************************************************************************/

static int
projectToFace(int    ifac,              /* (in)  Face index (bias-0) */
              double xyz_in[],          /* (in)  input coordinates */
              double xyz_out[])         /* (out) output coordinates */
{
    int    status = SUCCESS;            /* (out) return status */
    int       jtri, ip0, ip1, ip2;
    double    s0, s1, s01, dbest, dbest2, dtest2, xtest, ytest, ztest;
    double    x02, y02, z02, x12, y12, z12, xx2, yy2, zz2, A, B, C, D, E, F, G;

//    ROUTINE(projectToFace);

    /*----------------------------------------------------------------*/

    /* default output */
    xyz_out[0] = xyz_in[0];
    xyz_out[1] = xyz_in[1];
    xyz_out[2] = xyz_in[2];

    dbest  = 1000;
    dbest2 = SQR(dbest);

    /* loop through all triangles in iface */
    for (jtri = 0; jtri < fac[ifac].tess.ntri; jtri++) {

        /* bounding box check */
        if (xyz_in[0] < fac[ifac].tess.bbox[6*jtri  ]-dbest ||
            xyz_in[0] > fac[ifac].tess.bbox[6*jtri+1]+dbest ||
            xyz_in[1] < fac[ifac].tess.bbox[6*jtri+2]-dbest ||
            xyz_in[1] > fac[ifac].tess.bbox[6*jtri+3]+dbest ||
            xyz_in[2] < fac[ifac].tess.bbox[6*jtri+4]-dbest ||
            xyz_in[2] > fac[ifac].tess.bbox[6*jtri+5]+dbest   ) continue;

        /* determine barycentric coordinates of point closest to xyz_in */
        ip0 = fac[ifac].tess.trip[3*jtri  ];
        ip1 = fac[ifac].tess.trip[3*jtri+1];
        ip2 = fac[ifac].tess.trip[3*jtri+2];

        x02 = fac[ifac].tess.xyz[3*ip0  ] - fac[ifac].tess.xyz[3*ip2  ];
        y02 = fac[ifac].tess.xyz[3*ip0+1] - fac[ifac].tess.xyz[3*ip2+1];
        z02 = fac[ifac].tess.xyz[3*ip0+2] - fac[ifac].tess.xyz[3*ip2+2];
        x12 = fac[ifac].tess.xyz[3*ip1  ] - fac[ifac].tess.xyz[3*ip2  ];
        y12 = fac[ifac].tess.xyz[3*ip1+1] - fac[ifac].tess.xyz[3*ip2+1];
        z12 = fac[ifac].tess.xyz[3*ip1+2] - fac[ifac].tess.xyz[3*ip2+2];
        xx2 = xyz_in[0]                   - fac[ifac].tess.xyz[3*ip2  ];
        yy2 = xyz_in[1]                   - fac[ifac].tess.xyz[3*ip2+1];
        zz2 = xyz_in[2]                   - fac[ifac].tess.xyz[3*ip2+2];

        A = x02 * x02 + y02 * y02 + z02 * z02;
        B = x12 * x02 + y12 * y02 + z12 * z02;
        C = B;
        D = x12 * x12 + y12 * y12 + z12 * z12;
        E = xx2 * x02 + yy2 * y02 + zz2 * z02;
        F = xx2 * x12 + yy2 * y12 + zz2 * z12;
        G = A * D - B * C;

        if (fabs(G) < EPS20) continue;

        s0 = (E * D - B * F) / G;
        s1 = (A * F - E * C) / G;

        /* clip the barycentric coordinates and evaluate */
        s0 = MINMAX(0, s0, 1);
        s1 = MINMAX(0, s1, 1);

        s01 = s0 + s1;
        if (s01 > 1) {
            s0 /= s01;
            s1 /= s01;
        }

        xtest = fac[ifac].tess.xyz[3*ip2  ] + s0 * x02 + s1 * x12;
        ytest = fac[ifac].tess.xyz[3*ip2+1] + s0 * y02 + s1 * y12;
        ztest = fac[ifac].tess.xyz[3*ip2+2] + s0 * z02 + s1 * z12;

        /* if point is closer than anything we have seen so far, remember this */
        dtest2 = SQR(xtest-xyz_in[0]) + SQR(ytest-xyz_in[1]) + SQR(ztest-xyz_in[2]);

        if (dtest2 < dbest2) {
            xyz_out[0] = xtest;
            xyz_out[1] = ytest;
            xyz_out[2] = ztest;
            dbest2     = dtest2;
            dbest      = sqrt(dbest2);
        }
    }

//cleanup:
    return status;
}


/******************************************************************************/

static int
splitEdge(int    iedg,                  /* (in)  Edge index (bias-0) */
          int    ipnt)                  /* (in)  Point on Edge (bias-0) */
{
    int    status = SUCCESS;            /* (out) return status */

    int    imid, npnt, iend, ifac;

    ROUTINE(splitEdge);

    /*----------------------------------------------------------------*/

    iend = edg[iedg].iend;

    npnt = edg[iedg].npnt - ipnt;

    /* create a Node at ipnt */
    addNode(edg[iedg].pnt[ipnt], NULL);
    imid = Nnod - 1;

    /* modify iend of iedg to point to new Node */
    edg[iedg].iend = imid;

    /* modify the number of Points in iedg */
    edg[iedg].npnt = ipnt + 1;

    /* update the valence of the new Node and old iend */
    (nod[iend].nedg)--;
    (nod[imid].nedg)++;

    /* create a new Edge with the last Point of iedg */
    status = addEdge(edg[iedg].ileft, edg[iedg].irite, npnt, &(edg[iedg].pnt[ipnt]), imid, iend);
    CHECK_STATUS(addEdge);

    /* add the new Edge to the Faces that are pointed to by iedg */
#ifndef __clang_analyzer__
    ifac = edg[iedg].ileft;
    (fac[ifac].nedg)++;
    RALLOC(fac[ifac].edg, int, fac[ifac].nedg);
    fac[ifac].edg[fac[ifac].nedg-1] = Nedg - 1;

    ifac = edg[iedg].irite;
    (fac[ifac].nedg)++;
    RALLOC(fac[ifac].edg, int, fac[ifac].nedg);
    fac[ifac].edg[fac[ifac].nedg-1] = Nedg - 1;
#endif

cleanup:
    return status;
}


/******************************************************************************/

static int
trim(char   *s)                         /* (in)  string to trim */
{
    int n;                              /* (out) length of trimmed string */

    int i;

//    ROUTINE(trim);

    /*----------------------------------------------------------------*/

    /* remove leading white space */
    for (n = 0; n < strlen(s); n++) {
        if (s[n] != ' ' && s[n] != '\t' && s[n] != '\n') {
            break;
        }
    }

    for (i = 0; (i+n) < strlen(s); i++) {
        s[i] = s[i+n];
    }

    /* remove trailing white space */
    for (n = strlen(s)-1; n >= 1; n--) {
        if (s[n] != ' ' && s[n] != '\t' && s[n] != '\n') {
            break;
        }
    }

    s[n+1] = '\0';

//cleanup:
    return n;
}


/******************************************************************************/

static int
writeEgads(char   *filename)            /* (in)  file name */
{
    int    status = SUCCESS;            /* (out) return status */

    int    inod, iedg, ifac;
    FILE   *fp = NULL;

//    ROUTINE(writeEgads);

    /*----------------------------------------------------------------*/

    /* return if a topology was not created */
    if (Nnod <= 0 || Nedg <= 0 || Nfac <= 0) {
        printf("ERR: topology does not exist\a\n");
        goto cleanup;
    }

    /* return if all Faces do not have an associated surface */
    for (ifac = 0; ifac < Nfac; ifac++) {
        if (fac[ifac].imax <= 0 || fac[ifac].jmax <= 0) {
            printf("ERR: face %d does not have a surface\a\n", ifac);
            goto cleanup;
        }
    }

    /* create a curve for each Edge */
    for (iedg = 0; iedg < Nedg; iedg++) {
        makeCurve(iedg);
    }

    /* create a new EGADS contect */

    /* create the EGADS Nodes */
    for (inod = 0; inod < Nnod; inod++) {
    }

    /* create the EGADS Edges */
    for (iedg = 0; iedg < Nedg; iedg++) {
    }

    /* create the EGADS Faces */
    for (ifac = 0; ifac < Nfac; ifac++) {
    }

    /* create a shell */

    /* create a SOLIDBODY */

    /* create a MODEL */

    /* write the EGADS file */
    fp = fopen(filename, "w");
    if (fp != NULL) {
        fprintf(fp, "this is not an EGADS file\n");
        fclose(fp);
    } else {
        printf("ERR: file could not be opened\a\n");
        goto cleanup;
    }

    /* remove the context and everything in it */

cleanup:
    return status;
}


/******************************************************************************/

void
gvdraw(int    phase)                    /* (in)  =1 for DataBase drawing
                                                 =3 for 2D & 3D window drawing
                                                 =4 for Dials drawing
                                                 =5 for Key window drawing */

    /* called when the appropriate drawing phase is initialized
       used for programmer augmented rendering */
{
    if (phase == 99999) exit(0);
}


/******************************************************************************/

int
gvupdate()
{
    /* used for single process operation to allow the changing of data */

    int        answer;                  /* (out) =0 if no data has changed
                                                 >0 new number of graphic objects */

    int        i, itri, ifac;
    static int first = 1;

    /* set background to white and foreground to black */
    if (first == 1) {
        (void) gv_allocfamily("Nodes");
        (void) gv_allocfamily("Edges");
        (void) gv_allocfamily("Faces");
        (void) gv_allocfamily("Grids");
        (void) gv_allocfamily("Lines");

        GraphicGCSetFB(gv_wAux.GCs, gv_white, gv_black);
        GraphicBGColor(gv_wAux.wid, gv_black);

        GraphicGCSetFB(gv_wDial.GCs, gv_white, gv_black);
        GraphicBGColor(gv_wDial.wid, gv_black);

        first = 0;
    }

    /* simply return if no new data */
    if (new_data == 0) {
        answer = 0;
        goto cleanup;
    } else {
        new_data = 0;
    }

    /* remove any prevous graphic objects */
    for (i = 0; i < ngrobj; i++) {
        gv_free(grobjs[i], 2);
    }
    ngrobj = 0;

    /* remove the cached scalar if it exists */
    if (saveit != NULL) FREE(saveit);
    saveit = NULL;

    /* count number of hanging Sides */
    tess.nhang = 0;
    for (itri = 0; itri < tess.ntri; itri++) {
        if (tess.ttyp[itri] & TRI_VISIBLE) {
            if (tess.trit[3*itri  ] < 0) (tess.nhang)++;
            if (tess.trit[3*itri+1] < 0) (tess.nhang)++;
            if (tess.trit[3*itri+2] < 0) (tess.nhang)++;
        }
    }

    /* count number of linked Sides */
    tess.nlink = 0;
    for (itri = 0; itri < tess.ntri; itri++) {
        if (tess.ttyp[itri] & TRI_VISIBLE) {
            if (tess.ttyp[itri] & TRI_T0_LINK) (tess.nlink)++;
            if (tess.ttyp[itri] & TRI_T1_LINK) (tess.nlink)++;
            if (tess.ttyp[itri] & TRI_T2_LINK) (tess.nlink)++;
        }
    }
    printf("    npnt=%d   ntri=%d   nhang=%d   nlink=%d   ncolr=%d\n",
           tess.npnt, tess.ntri, tess.nhang, tess.nlink, tess.ncolr);

    /* return the number of graphic objects */
    answer = 0;

    /* Triangles */
    answer++;

    /* Nodes */
    answer += Nnod;

    /* Edges */
    answer += Nedg;

    /* Faces */
    answer += Nfac;

    /* Grids */
    for (ifac = 0; ifac < Nfac; ifac++) {
        if (fac[ifac].imax > 0 && fac[ifac].jmax > 0) {
            answer++;
        }
    }

    /* hanging Sides */
    if (tess.nhang > 0) answer++;

    /* Links */
    if (tess.nlink > 0) answer++;

    /* marked Point */
    if (markedPnt >= 0) answer++;

    /* marked Triangle */
    if (markedTri >= 0) answer++;

    /* return the number of graphic objects */
cleanup:
    return answer;
}


/******************************************************************************/

void
gvdata(int       ngraphics,             /* (in)  number of graphic objects to fill */
       GvGraphic *graphic[])            /* (both) graphic objects */
{
    /* used to (re)set the graphics objects to be used in plotting */

    int      ng, k, itri, ip0, ip1, ip2, ipnt, inod, iedg, ifac, imax, jmax;
    int      attr, utype, uindex;
    char     title[80];
    GvColor  color;

    ngrobj = ngraphics;
    grobjs =  graphic;
    ng     = 0;

    /* Triangles */
    if (tess.ntri > 0) {
        color.red   = 1.0;
        color.green = 0.0;
        color.blue  = 0.0;

        if (Nfac > 0) {
            attr    =             GV_ORIENTATION | GV_FACETLIGHT | GV_FORWARD;
        } else if (gridOn == 1) {
            attr    = GV_SCALAR | GV_ORIENTATION | GV_FACETLIGHT | GV_FORWARD | GV_MESH;
        } else {
            attr    = GV_SCALAR | GV_ORIENTATION | GV_FACETLIGHT | GV_FORWARD;
        }
        utype       = 1;
        uindex      = 0;
        graphic[ng] = gv_alloc(GV_INDEXED, GV_DISJOINTTRIANGLES,
                               attr, color, "Triangles", utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->back.red   = 0.5;
            graphic[ng]->back.green = 0.5;
            graphic[ng]->back.blue  = 0.5;
            graphic[ng]->mesh.red   = 0.0;
            graphic[ng]->mesh.green = 0.0;
            graphic[ng]->mesh.blue  = 0.0;
            graphic[ng]->number     = 1;
            graphic[ng]->fdata      = (float *)malloc(3*(tess.npnt+1)*sizeof(float));

            k = 0;
            for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
                graphic[ng]->fdata[k++] = (float)tess.xyz[3*ipnt  ];
                graphic[ng]->fdata[k++] = (float)tess.xyz[3*ipnt+1];
                graphic[ng]->fdata[k++] = (float)tess.xyz[3*ipnt+2];
            }
            graphic[ng]->fdata[k++] = (float)0;
            graphic[ng]->fdata[k++] = (float)0;
            graphic[ng]->fdata[k++] = (float)0;

            graphic[ng]->object->length = tess.ntri;
            graphic[ng]->object->type.distris.index = (int *)malloc(3*tess.ntri*sizeof(int));

            k = 0;
            for (itri = 0; itri < tess.ntri; itri++) {
                if (tess.ttyp[itri] & TRI_VISIBLE) {
                    graphic[ng]->object->type.distris.index[k++] = tess.trip[3*itri  ];
                    graphic[ng]->object->type.distris.index[k++] = tess.trip[3*itri+1];
                    graphic[ng]->object->type.distris.index[k++] = tess.trip[3*itri+2];
                } else {
                    graphic[ng]->object->type.distris.index[k++] = tess.npnt;
                    graphic[ng]->object->type.distris.index[k++] = tess.npnt;
                    graphic[ng]->object->type.distris.index[k++] = tess.npnt;
                }
            }
        }
        ng++;
    }

    /* Nodes */
    for (inod = 0; inod < Nnod; inod++) {
        color.red   = 1.0;
        color.green = 0.0;
        color.blue  = 1.0;

        sprintf(title, "Node %4d", inod);

        attr        = GV_FOREGROUND | GV_FORWARD;
        utype       = 2;
        uindex      = inod;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_POINTS,
                               attr, color, title, utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->pointSize = 10;

            graphic[ng]->object->length = 1;
            graphic[ng]->fdata = (float *)malloc(3*sizeof(float));

            graphic[ng]->fdata[0] = (float)tess.xyz[3*nod[inod].ipnt  ];
            graphic[ng]->fdata[1] = (float)tess.xyz[3*nod[inod].ipnt+1];
            graphic[ng]->fdata[2] = (float)tess.xyz[3*nod[inod].ipnt+2];

            gv_adopt("Nodes", graphic[ng]);
        }
        ng++;
    }

    /* Edges */
    for (iedg = 0; iedg < Nedg; iedg++) {
        color.red   = 1.0;
        color.green = 0.0;
        color.blue  = 1.0;

        sprintf(title, "Edge %4d", iedg);

        if (edg[iedg].mark >= 0) {
            attr    = GV_FOREGROUND | GV_FORWARD | GV_MESH;
        } else {
            attr    = GV_FOREGROUND | GV_FORWARD;
        }
        utype       = 3;
        uindex      = iedg;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_POLYLINES,
                               attr, color, title, utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->mesh.red   = 0.0;
            graphic[ng]->mesh.green = 0.0;
            graphic[ng]->mesh.blue  = 0.0;
            graphic[ng]->number     = 1;
            graphic[ng]->lineWidth  = 3;
            graphic[ng]->pointSize  = 5;

            graphic[ng]->object->length = 1;
            graphic[ng]->object->type.plines.len = (int *)malloc(sizeof(int));
            graphic[ng]->object->type.plines.len[0] = edg[iedg].npnt;

            graphic[ng]->fdata = (float *)malloc(3*edg[iedg].npnt*sizeof(float));
            for (ipnt = 0; ipnt < edg[iedg].npnt; ipnt++) {
                graphic[ng]->fdata[3*ipnt  ] = (float)tess.xyz[3*edg[iedg].pnt[ipnt]  ];
                graphic[ng]->fdata[3*ipnt+1] = (float)tess.xyz[3*edg[iedg].pnt[ipnt]+1];
                graphic[ng]->fdata[3*ipnt+2] = (float)tess.xyz[3*edg[iedg].pnt[ipnt]+2];
            }

            gv_adopt("Edges", graphic[ng]);
        }
        ng++;
    }

    /* Faces */
    for (ifac = 0; ifac < Nfac; ifac++) {
        color.red   = 1.0;
        color.green = 1.0;
        color.blue  = 0.5;

        sprintf(title, "Face %4d", ifac);

        attr   = GV_FOREGROUND | GV_ORIENTATION | GV_FACETLIGHT | GV_FORWARD;
        utype  = 4;
        uindex = ifac;
        graphic[ng] = gv_alloc(GV_INDEXED, GV_DISJOINTTRIANGLES,
                               attr, color, title, utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv+alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->back.red   = 0.5;
            graphic[ng]->back.green = 0.5;
            graphic[ng]->back.blue  = 0.5;
            graphic[ng]->mesh.red   = 0.0;
            graphic[ng]->mesh.green = 0.0;
            graphic[ng]->mesh.blue  = 0.0;
            graphic[ng]->number     = 1;
            graphic[ng]->fdata      = (float *)malloc(3*fac[ifac].tess.npnt*sizeof(float));

            k = 0;
            for (ipnt = 0; ipnt < fac[ifac].tess.npnt; ipnt++) {
                graphic[ng]->fdata[k++] = (float)(fac[ifac].tess.xyz[3*ipnt  ]);
                graphic[ng]->fdata[k++] = (float)(fac[ifac].tess.xyz[3*ipnt+1]);
                graphic[ng]->fdata[k++] = (float)(fac[ifac].tess.xyz[3*ipnt+2]);
            }

            graphic[ng]->object->length = fac[ifac].tess.ntri;
            graphic[ng]->object->type.distris.index = (int *)malloc(3*fac[ifac].tess.ntri*sizeof(int));

            k = 0;
            for (itri = 0; itri < fac[ifac].tess.ntri; itri++) {
                graphic[ng]->object->type.distris.index[k++] = fac[ifac].tess.trip[3*itri  ];
                graphic[ng]->object->type.distris.index[k++] = fac[ifac].tess.trip[3*itri+1];
                graphic[ng]->object->type.distris.index[k++] = fac[ifac].tess.trip[3*itri+2];
            }

            gv_adopt("Faces", graphic[ng]);
        }
        ng++;
    }

    /* Grids */
    for (ifac = 0; ifac < Nfac; ifac++) {
        imax = fac[ifac].imax;
        jmax = fac[ifac].jmax;

        if (imax <= 0 || jmax <= 0) continue;

        color.red   = 1.0;
        color.green = 1.0;
        color.blue  = 0.0;

        sprintf(title, "Grid %4d", ifac);

        attr        = GV_ORIENTATION | GV_FACETLIGHT | GV_FORWARD | GV_MESH;
        utype       = 5;
        uindex      = ifac;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_QUADMESHS,
                               attr, color, title, utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->back.red   = 0.5;
            graphic[ng]->back.green = 0.5;
            graphic[ng]->back.blue  = 0.5;
            graphic[ng]->mesh.red   = 0.0;
            graphic[ng]->mesh.green = 0.0;
            graphic[ng]->mesh.blue  = 0.0;
            graphic[ng]->number     = 1;
            graphic[ng]->fdata      = (float *)malloc(3*imax*jmax*sizeof(float));

            for (k = 0; k < imax*jmax; k++) {
                graphic[ng]->fdata[3*k  ] = (float)fac[ifac].xsrf[k];
                graphic[ng]->fdata[3*k+1] = (float)fac[ifac].ysrf[k];
                graphic[ng]->fdata[3*k+2] = (float)fac[ifac].zsrf[k];
            }

            graphic[ng]->object->length = 1;
            graphic[ng]->object->type.qmeshes.size = (int *)malloc(2*sizeof(int));
            graphic[ng]->object->type.qmeshes.size[0] = imax;
            graphic[ng]->object->type.qmeshes.size[1] = jmax;

            gv_adopt("Grids", graphic[ng]);
        }
        ng++;
    }

    /* hanging Sides */
    if (tess.nhang > 0) {
        color.red   = 0.0;
        color.green = 1.0;
        color.blue  = 1.0;

        attr        = GV_FOREGROUND | GV_FORWARD;
        utype       = 11;
        uindex      =  0;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_DISJOINTLINES,
                               attr, color, "Hanging", utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->mesh.red   = 0.0;
            graphic[ng]->mesh.green = 0.0;
            graphic[ng]->mesh.blue  = 0.0;
            graphic[ng]->number     = 1;
            graphic[ng]->lineWidth  = 3;
            graphic[ng]->pointSize  = 10;

            graphic[ng]->object->length = tess.nhang;
            graphic[ng]->fdata = (float *)malloc(6*tess.nhang*sizeof(float));

            k = 0;
            for (itri = 0; itri < tess.ntri; itri++) {
                if (tess.ttyp[itri] & TRI_VISIBLE) {
                    ip0 = tess.trip[3*itri  ];
                    ip1 = tess.trip[3*itri+1];
                    ip2 = tess.trip[3*itri+2];

                    if (tess.trit[3*itri  ] < 0) {
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+2];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+2];
                    }
                    if (tess.trit[3*itri+1] < 0) {
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+2];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+2];
                    }
                    if (tess.trit[3*itri+2] < 0) {
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+2];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+2];
                    }
                }
            }
        }
        ng++;
    }

    /* linked Sides */
    if (tess.nlink > 0) {
        color.red   = 1.0;
        color.green = 1.0;
        color.blue  = 1.0;

        attr        = GV_FOREGROUND | GV_FORWARD;
        utype       = 12;
        uindex      =  0;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_DISJOINTLINES,
                               attr, color, "Linked", utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->mesh.red   = 0.0;
            graphic[ng]->mesh.green = 0.0;
            graphic[ng]->mesh.blue  = 0.0;
            graphic[ng]->number     = 1;
            graphic[ng]->lineWidth  = 3;
            graphic[ng]->pointSize  = 5;

            graphic[ng]->object->length = tess.nlink;
            graphic[ng]->fdata = (float *)malloc(6*tess.nlink*sizeof(float));

            k = 0;
            for (itri = 0; itri < tess.ntri; itri++) {
                if (tess.ttyp[itri] & TRI_VISIBLE) {
                    ip0 = tess.trip[3*itri  ];
                    ip1 = tess.trip[3*itri+1];
                    ip2 = tess.trip[3*itri+2];

                    if (tess.ttyp[itri] & TRI_T0_LINK) {
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+2];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+2];
                    }
                    if (tess.ttyp[itri] & TRI_T1_LINK) {
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip2+2];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+2];
                    }
                    if (tess.ttyp[itri] & TRI_T2_LINK) {
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip0+2];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1  ];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+1];
                        graphic[ng]->fdata[k++] = (float)tess.xyz[3*ip1+2];
                    }
                }
            }
        }
        ng++;
    }

    /* marked Point */
    if (markedPnt >= 0) {
        color.red   = 0.0;
        color.green = 0.0;
        color.blue  = 0.0;

        attr        = GV_FOREGROUND | GV_FORWARD;
        utype       = 13;
        uindex      = markedPnt;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_POINTS,
                               attr, color, "marked Pnt", utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR:  gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->pointSize  = 10;

            graphic[ng]->object->length = 1;
            graphic[ng]->fdata = (float *)malloc(3*sizeof(float));

            graphic[ng]->fdata[0] = (float)tess.xyz[3*markedPnt  ];
            graphic[ng]->fdata[1] = (float)tess.xyz[3*markedPnt+1];
            graphic[ng]->fdata[2] = (float)tess.xyz[3*markedPnt+2];
        }
        ng++;
    }

    /* marked Triangle */
    if (markedTri >= 0) {
        color.red   = 0.0;
        color.green = 1.0;
        color.blue  = 1.0;

        attr        = GV_FOREGROUND | GV_FORWARD;
        utype       = 14;
        uindex      = markedTri;
        graphic[ng] = gv_alloc(GV_NONINDEXED, GV_POINTS,
                               attr, color, "marked Tri", utype, uindex);
        if (graphic[ng] == NULL) {
            printf("ERR: gv_alloc error on graphic[%d]\a\n", ng);
        } else {
            graphic[ng]->pointSize  = 10;

            graphic[ng]->object->length = 1;
            graphic[ng]->fdata = (float *)malloc(3*sizeof(float));

            ip0 = tess.trip[3*markedTri  ];
            ip1 = tess.trip[3*markedTri+1];
            ip2 = tess.trip[3*markedTri+2];

            graphic[ng]->fdata[0] = (float)(tess.xyz[3*ip0  ]
                                           +tess.xyz[3*ip1  ]
                                           +tess.xyz[3*ip2  ])/3;
            graphic[ng]->fdata[1] = (float)(tess.xyz[3*ip0+1]
                                           +tess.xyz[3*ip1+1]
                                           +tess.xyz[3*ip2+1])/3;
            graphic[ng]->fdata[2] = (float)(tess.xyz[3*ip0+2]
                                           +tess.xyz[3*ip1+2]
                                           +tess.xyz[3*ip2+2])/3;
        }
        ng++;
    }

    /* set update flag to done */
    new_data = 0;
}


/******************************************************************************/

int
gvscalar(int        key,                /* (in)  scalar index (from gv_init) */
/*@unused@*/GvGraphic  *graphic,        /* (in)  the GvGraphic structure for scalar fill */
         int        len,                /* (in)  length of sclar to be filled */
         float      *scalar )           /* (both) scalar to be filled */
{
    int    status = SUCCESS;            /* (out) return status */

    int    ipnt, itri, answer=0;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    ROUTINE(gvscalar);

    /*----------------------------------------------------------------*/

    /* use cached scalar if it exists */
    if (saveit != NULL) {
        memcpy(scalar, saveit, len*sizeof(float));
        answer = 1;
        goto cleanup;
    }

    /* build a new scalar */
    if        (key == 0) {
        if (len != tess.ntri) goto cleanup;
        if (tess.ncolr == 0) {
            for (itri = 0; itri < tess.ntri; itri++) {
                scalar[itri] = 0;
            }
        } else {
            for (itri = 0; itri < tess.ntri; itri++) {
                scalar[itri] = (float)(tess.ttyp[itri] & TRI_COLOR) / (float)(tess.ncolr);
            }
        }
    } else if (key == 1) {
        if (len != tess.npnt) goto cleanup;
        xmin = tess.xyz[3*0  ];
        xmax = tess.xyz[3*0  ];
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (tess.xyz[3*ipnt  ] < xmin) xmin = tess.xyz[3*ipnt  ];
            if (tess.xyz[3*ipnt  ] > xmax) xmax = tess.xyz[3*ipnt  ];
        }
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            scalar[ipnt] = (float)(tess.xyz[3*ipnt  ] - xmin) / (float)(xmax - xmin);
        }
    } else if (key == 2) {
        if (len != tess.npnt) goto cleanup;
        ymin = tess.xyz[3*0+1];
        ymax = tess.xyz[3*0+1];
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (tess.xyz[3*ipnt+1] < ymin) ymin = tess.xyz[3*ipnt+1];
            if (tess.xyz[3*ipnt+1] > ymax) ymax = tess.xyz[3*ipnt+1];
        }
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            scalar[ipnt] = (float)(tess.xyz[3*ipnt+1] - ymin) / (float)(ymax - ymin);
        }
    } else if (key == 3) {
        if (len != tess.npnt) goto cleanup;
        zmin = tess.xyz[3*0+2];
        zmax = tess.xyz[3*0+2];
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            if (tess.xyz[3*ipnt+2] < zmin) zmin = tess.xyz[3*ipnt+2];
            if (tess.xyz[3*ipnt+2] > zmax) zmax = tess.xyz[3*ipnt+2];
        }
        for (ipnt = 0; ipnt < tess.npnt; ipnt++) {
            scalar[ipnt] = (float)(tess.xyz[3*ipnt+2] - zmin) / (float)(zmax - zmin);
        }
    } else {
        for (ipnt = 0; ipnt < len; ipnt++) {
            scalar[ipnt] = 0;
        }
    }

    /* cache result for later use */
    MALLOC(saveit, float, len);
    memcpy(saveit, scalar, len*sizeof(float));

    answer = 1;

cleanup:
    if (status < 0) {
        return status;
    } else {
        return answer;
    }
}


/******************************************************************************/

void
gvevent(int    *win,                    /* (in)  window that fielded event */
        int    *type,                   /* (in)  type of event */
        int    *xscr,                   /* (in)  x-screen location */
        int    *yscr,                   /* (in)  y-screen location */
        int    *state)                  /* (in)  state of event */
{
    /* allows modification of events */

    int    status = 0;

    int    ipnt, itri, jtri, nchange, ndelete, swap, i, j, itype, icolr;
    int    inod, iedg, ifac, imax, jmax;
    double dum, temp0, temp1, temp2, size, cosrot, sinrot, box[6], val;
    char   scriptName[81], comment[81], filename[81];
    FILE   *fp = NULL;

    do {
        /* get the next script line if we are reading a script (and insert a '$' if
              we have detected an EOF) */
        if (script != NULL) {
            if (fscanf(script, "%1s", (char *)state) != 1) {
                *state = '$';
            }
            *win  = ThreeD;
            *type = KeyPress;
        }

        if ((*win == ThreeD) && (*type == KeyPress)) {

            if (*state == '\0') {
                dum = getDbl("Dummy should not be called");
                printf("dum=%f\n", dum);

            /* 'b' - bridge to marked Triangle */
            } else if (*state == 'b') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (markedTri >= 0 && itri >= 0) {
                    printf("==> Option 'b' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d b\n", itri);

                    bridgeTriangles(&tess, itri, markedTri);

                    markedTri = -1;
                    new_data  = 1;
                } else if (markedTri < 0) {
                    printf("ERR: No marked Triangle\a\n");
                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'c' - color Triangle and neighbors */
            } else if (*state == 'c') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0) {
                    printf("==> Option 'c' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d c\n", itri);

                    printf("    current color = %d\n", tess.ttyp[itri] & TRI_COLOR);

                    icolr = getInt("Enter color (or -1): ");

                    if (icolr >= 0) {
                        tess.ncolr = MAX(tess.ncolr, icolr);
                        colorTriangles(&tess, itri, icolr);
                        new_data = 1;
                    } else {
                        printf("ERR: Color not applied\a\n");
                    }

                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'C' - sort Triangles by color */
            } else if (*state == 'C') {
                printf("==> Option 'C'\n");
                fprintf(fpdump, " C\n");

                sortTriangles(&tess);

                new_data = 1;

            /* 'd' - delete a Triangle */
            } else if (*state == 'd') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0) {
                    printf("==> Option 'd' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d d\n", itri);

                    deleteTriangle(&tess, itri);
                    new_data = 1;
                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'e' - classify Edge */
            } else if (*state == 'e') {
                if (numarg >= 0) {
                    iedg   = numarg;
                    numarg = -1;
                } else {
                    iedg   = findEdge(*xscr, *yscr, &ipnt);
                }

                if (iedg >= 0) {
                    printf("==> Option 'e' (iedg=%6d)\n", iedg);
                    fprintf(fpdump, "%6d e\n", iedg);

                    edg[iedg].mark = getInt("Enter 0=vmin, 1=umax, 2=vmax, 3=umin: ");
                    new_data = 1;
                } else {
                    printf("ERR: No Edge found\a\n");
                }

            /* 'E' - write an EGADS file */
            } else if (*state == 'E') {
                printf("==> Option 'E' \n");
                fprintf(fpdump, "E\n");

                getStr("Enter filename: ", filename);
                trim(filename);

                writeEgads(filename);

            /* 'f' - fill loop adjacent to Point */
            } else if (*state == 'f') {
                if (numarg >= 0) {
                    ipnt   = numarg;
                    numarg = -1;
                } else {
                    ipnt   = findPoint(*xscr, *yscr);
                }

                if (ipnt >= 0) {
                    printf("==> Option 'f' (ipnt=%6d)\n", ipnt);
                    fprintf(fpdump, "%6d f\n", ipnt);

                    fillLoop(&tess, ipnt);
                    new_data = 1;
                } else {
                    printf("ERR: No Point found\a\n");
                }

            /* 'F' - create Surface for Face */
            } else if (*state == 'F') {
                if (numarg >= 0) {
                    ifac   = numarg;
                    numarg = -1;
                    itri   = getInt("Enter Triangle index: ");
                } else {
                    ifac   = findFace(*xscr, *yscr, &itri);
                }

                if (ifac >= 0 && itri >= 0) {
                    printf("==> Option 'F' (ifac=%6d, otri=%6d)\n", ifac, itri);
                    fprintf(fpdump, "%6d F %6d\n", ifac, itri);

                    itype = getInt("Enter 1=TFI, 2=cuts: ");
                    imax  = getInt("Enter imax: ");
                    jmax  = getInt("Enter jmax: ");
                    if (itype == 1 || itype == 2) {
                        makeSurface(ifac, itype, imax, jmax);
                    } else {
                        printf("ERR: Bad syrface type\a\n");
                    }

                    new_data = 1;
                } else {
                    printf("ERR: No Face and/or riangle found\a\n");
                }

            /* 'G' - toggle grid visibility */
            } else if (*state == 'G') {
                printf("==> Option 'G' \n");

                gridOn   = 1 - gridOn;
                new_data = 1;

            /* 'h' - hide Triangle and neighbors */
            } else if (*state == 'h') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0) {
                    printf("==> Option 'h' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d h\n", itri);

                    tess.ttyp[itri] &= ~TRI_VISIBLE;

                    /* flood-fill to the neighbors */
                    nchange = 1;
                    while (nchange > 0) {
                        nchange = 0;
                        for (itri = 0; itri < tess.ntri; itri++) {
                            if ((tess.ttyp[itri] & TRI_VISIBLE) == 0) {
                                jtri = tess.trit[3*itri  ];
                                if (jtri >= 0) {
                                    if (tess.ttyp[jtri] &   TRI_VISIBLE) {
                                        tess.ttyp[jtri] &= ~TRI_VISIBLE;
                                        nchange++;
                                    }
                                }

                                jtri = tess.trit[3*itri+1];
                                if (jtri >= 0) {
                                    if (tess.ttyp[jtri] &   TRI_VISIBLE) {
                                        tess.ttyp[jtri] &= ~TRI_VISIBLE;
                                        nchange++;
                                    }
                                }

                                jtri = tess.trit[3*itri+2];
                                if (jtri >= 0) {
                                    if (tess.ttyp[jtri] &   TRI_VISIBLE) {
                                        tess.ttyp[jtri] &= ~TRI_VISIBLE;
                                        nchange++;
                                    }
                                }
                            }
                        }
                    }

                    new_data = 1;
                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'H' - hide Triangle with same color */
            } else if (*state == 'H') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0) {
                    printf("==> Option 'H' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d h\n", itri);

                    tess.ttyp[itri] &= ~TRI_VISIBLE;

                    icolr = tess.ttyp[itri] & TRI_COLOR;

                    for (itri = 0; itri < tess.ntri; itri++) {
                        if ((tess.ttyp[itri] & TRI_COLOR) == icolr) {
                            tess.ttyp[itri] &= ~TRI_VISIBLE;
                        }
                    }

                    new_data = 1;
                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'i' - invert all visible Triangles */
            } else if (*state == 'i') {
                printf("==> Option 'i' \n");
                fprintf(fpdump, "i\n");

                for (itri = 0; itri < tess.ntri; itri++) {
                    if (tess.ttyp[itri] & TRI_VISIBLE) {
                        swap                = tess.trip[3*itri+1];
                        tess.trip[3*itri+1] = tess.trip[3*itri+2];
                        tess.trip[3*itri+2] = swap;

                        swap                = tess.trip[3*itri+1];
                        tess.trip[3*itri+1] = tess.trip[3*itri+2];
                        tess.trip[3*itri+2] = swap;
                    }
                }

                new_data = 1;

            /* 'j' - join Points */
            } else if (*state == 'j' || *state == 'J') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0 && markedTri >= 0) {
                    printf("==> Option 'j' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d j\n", itri);

                    joinPoints(&tess, itri, markedTri);

                    markedTri = -1;
                    new_data = 1;
                } else {
                    printf("ERR: No Triangle found\a\n");
                }


            /* 'k'  - kill Triangles with same color */
            } else if (*state == 'k') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri   = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0) {
                    printf("==> Option 'k' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d k\n", itri);

                    icolr = tess.ttyp[itri] & TRI_COLOR;

                    ndelete = 0;
                    for (jtri = 0; jtri < tess.ntri; jtri++) {
                        if ((tess.ttyp[jtri] & TRI_COLOR) == icolr) {
                            deleteTriangle(&tess, jtri);
                            ndelete++;
                        }
                    }
                    printf("    %d Triangles deleted\n", ndelete);

                    new_data = 1;
                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'l' - link Points */
            } else if (*state == 'l') {
                if (numarg >= 0) {
                    ipnt   = numarg;
                    numarg = -1;
                } else {
                    ipnt    = findPoint(*xscr, *yscr);
                }

                if (ipnt >= 0 && markedPnt >= 0) {
                    fprintf(fpdump, "%6d l\n", ipnt);
                    printf("==> Option 'l' (ipnt=%6d)\n", ipnt);

                    createLinks(&tess, markedPnt, ipnt);

                    markedPnt = ipnt;
                } else if (ipnt < 0) {
                    printf("ERR: No Point found\a\n");
                } else {
                    printf("ERR: No marked Point\a\n");
                }

                new_data = 1;

            /* 'L' - make Links between colors */
            } else if (*state == 'L') {
                makeLinks(&tess);

                new_data = 1;

            /* 'm' - make topology */
            } else if (*state == 'm') {
                fprintf(fpdump, "m\n");
                printf("==> Option 'm' \n");

                status = makeTopology();
                if (status < 0) printf("makeTopology -> status=%d\n", status);

                new_data = 1;

            /* 'p' - toggle Point mark */
            } else if (*state == 'p') {
                if (numarg >= 0) {
                    ipnt   = numarg;
                    numarg = -1;
                } else {
                    ipnt   = findPoint(*xscr, *yscr);
                }

                if (ipnt >= 0) {
                    printf("==> Option 'p' (ipnt=%6d)\n", ipnt);
                    fprintf(fpdump, "%6d p\n", ipnt);

                    if (ipnt != markedPnt) {
                        markedPnt = ipnt;
                    } else {
                        printf("    Unmarking Point\a\n");
                        markedPnt = -1;
                    }

                    new_data = 1;
                } else {
                    printf("ERR: No Point found\a\n");
                }

            /* 'q' - query topology */
            } else if (*state == 'q') {
                printf("==> Option 'q' \n");
                fprintf(fpdump, "q\n");

                for (inod = 0; inod < Nnod; inod++) {
                    ipnt = nod[inod].ipnt;
                    printf("    Node[%4d].ipnt=%6d  .nedg=%2d  .x=%10.4f  .y=%10.4f  .z=%10.4f\n",
                           inod, ipnt, nod[inod].nedg, tess.xyz[3*ipnt  ], tess.xyz[3*ipnt+1], tess.xyz[3*ipnt+2]);
                }
                for (iedg = 0; iedg < Nedg; iedg++) {
                    printf("    Edge[%4d].ibeg=%6d  .iend=%6d  .ileft=%6d  .irite=%6d  .npnt=%6d   .mark=%2d\n",
                           iedg, edg[iedg].ibeg, edg[iedg].iend,
                           edg[iedg].ileft, edg[iedg].irite, edg[iedg].npnt, edg[iedg].mark);
                }
                for (ifac = 0; ifac < Nfac; ifac++) {
                    printf("    Face[%4d].icol=%2d  .ntri=%6d  .npnt=%6d  .nedg=%2d:",
                           ifac, fac[ifac].icol, fac[ifac].tess.ntri, fac[ifac].tess.npnt, fac[ifac].nedg);
                    for (iedg = 0; iedg < fac[ifac].nedg; iedg++) {
                        printf("  %6d", fac[ifac].edg[iedg]);
                    }
                    printf("\n");
                }

            /* 's' - show all Triangles */
            } else if (*state == 's') {
                printf("==> Option 's' \n");
                fprintf(fpdump, "s\n");

                for (itri = 0; itri < tess.ntri; itri++) {
                    if (tess.ttyp[itri] &  TRI_ACTIVE) {
                        tess.ttyp[itri] |= TRI_VISIBLE;
                    }
                }

                new_data = 1;

            /* 'S' - split Edge */
            } else if (*state == 'S') {
                if (numarg >= 0) {
                    iedg   = numarg;
                    numarg = -1;
                    ipnt   = getInt("Enter Point index: ");
                } else {
                    iedg   = findEdge(*xscr, *yscr, &ipnt);
                }

                if (iedg >= 0 && ipnt >= 0) {
                    printf("==> Option 'S' (iedg=%6d, ipnt=%6d)\n", iedg, ipnt);
                    fprintf(fpdump, "%6d S %6d\n", iedg, ipnt);

                    status = splitEdge(iedg, ipnt);
                    if (status < SUCCESS) printf("splitEdge -> status=%d\n", status);

                    new_data = 1;
                } else {
                    printf("ERR: No Edge and/or Point found\a\n");
                }


            /* 't' - toggle Triangle mark */
            } else if (*state == 't') {
                if (numarg >= 0) {
                    itri   = numarg;
                    numarg = -1;
                } else {
                    itri = findTriangle(*xscr, *yscr);
                }

                if (itri >= 0) {
                    printf("==> Option 't' (itri=%6d)\n", itri);
                    fprintf(fpdump, "%6d t\n", itri);

                    if (itri != markedTri) {
                        markedTri = itri;
                    } else {
                        printf("    Unmarking Triangle\a\n");
                        markedTri = -1;
                    }

                    new_data = 1;
                } else {
                    printf("ERR: No Triangle found\a\n");
                }

            /* 'T' - test face creation */
            } else if (*state == 'T') {
                if (numarg >= 0) {
                    icolr  = numarg;
                    numarg = -1;
                } else {
                    icolr = getInt("Enter icolr: ");
                }

                printf("==> Option 'T' (icolr=%d)\n", icolr);
                fprintf(fpdump, "%6d T\n", icolr);

                {
                    tess_T subTess1, subTess2;
//$$$                    int    iview=1;
                    int    nneg, npos, iloop, nloop, ibeg[10];
                    double alen[10];
#ifdef GRAFIC
                    int    io_kbd=5, io_scr=6, indgr=1+2+4+16+64, igrid=1, isymb=1;
                    char   pltitl[255];
#endif

                    /* make a TESS from the selected color */
                    status = initialTess(&subTess1);
                    printf("initialTess -> status=%d\n", status);

                    status = extractColor(&tess, icolr, &subTess1);
                    printf("extractColor(icolr=%d) -> status=%d\n", icolr, status);

                    printf("    subTess1.npnt = %d\n", subTess1.npnt);
                    printf("    subTess1.ntri = %d\n", subTess1.ntri);

//$$$                    sprintf(pltitl, "~x~y~ 3D plot of color %d", icolr);
//$$$                    grinit_(&io_kbd, &io_scr, "test", strlen("test"));
//$$$                    grctrl_(plot3D, &indgr, pltitl,
//$$$                            (void*)(&subTess1),
//$$$                            (void*)(&iview),
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            strlen(pltitl));

                    /* make a copy (which might have filled-in loops) */
                    status = copyTess(&subTess1, &subTess2);
                    printf("copyTess -> status=%d\n", status);

                    /* if there are inner loops, fill them */
                    nloop = 10;
                    status = findLoops(&subTess2, &nloop, ibeg, alen);
                    printf("findLoops -> status=%d, nloop=%d\n", status, nloop);

                    for (iloop = 1; iloop < nloop; iloop++) {
                        status = fillLoop(&subTess2, ibeg[iloop]);

                        printf("fillLoop(ibeg=%d) -> status=%d\n", ibeg[iloop], status);
                    }

                    /* generate the initial UV by projection */
                    status = initialUV(&subTess2);
                    printf("initialUV -> status=%d\n", status);

                    /* flip this over in U if there are more negartive areas than positive areas */
                    status = checkAreas(&subTess2, &nneg, &npos);
                    printf("checkAreas -> status=%d, nneg=%d, npos=%d, ntri=%d\n",
                           status, nneg, npos, subTess1.ntri);

                    if (nneg > npos) {
                        printf("flipping...\n");
                        for (ipnt = 0; ipnt < subTess1.npnt; ipnt++) {
                            subTess2.uv[2*ipnt] *= -1;
                        }

                        status = checkAreas(&subTess2, &nneg, &npos);
                        printf("checkAreas -> status=%d, nneg=%d, npos=%d, ntri=%d\n",
                               status, nneg, npos, subTess1.ntri);
                    }

//$$$                    sprintf(pltitl, "~u~v~ 2D plot of color %d", icolr);
//$$$                    grinit_(&io_kbd, &io_scr, "test", strlen("test"));
//$$$                    grctrl_(plot2D, &indgr, pltitl,
//$$$                            (void*)(&subTess2),
//$$$                            (void*)(&igrid),
//$$$                            (void*)(&isymb),
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            (void*)NULL,
//$$$                            strlen(pltitl));

                    /* perform floater smoothing */
                    status = floaterUV(&subTess2);
                    printf("floaterUV -> status=%d\n", status);

                    status = checkAreas(&subTess2, &nneg, &npos);
                    printf("checkAreas -> status=%d, nneg=%d, npos=%d, ntri=%d\n",
                           status, nneg, npos, subTess1.ntri);

                    /* copy the UV in the TESS with the filled-in loops back
                       into the origian TESS */
                    for (ipnt = 0; ipnt < subTess1.npnt; ipnt++) {
                        subTess1.uv[2*ipnt  ] = subTess2.uv[2*ipnt  ];
                        subTess1.uv[2*ipnt+1] = subTess2.uv[2*ipnt+1];
                    }

#ifdef GRAFIC
                    sprintf(pltitl, "~u~v~ 2D plot of color %d", icolr);
                    grinit_(&io_kbd, &io_scr, "test", strlen("test"));
                    grctrl_(plot2D, &indgr, pltitl,
                            (void*)(&subTess1),
                            (void*)(&igrid),
                            (void*)(&isymb),
                            (void*)NULL,
                            (void*)NULL,
                            (void*)NULL,
                            (void*)NULL,
                            (void*)NULL,
                            (void*)NULL,
                            (void*)NULL,
                            strlen(pltitl));
#endif

                    /* create the bicubic splines if all the Traingle areas are positive */
                    if (nneg == 0) {

                    }

                    /* clean up */
                    status = freeTess(&subTess2);
                    printf("freeTess -> status=%d\n", status);

                    status = freeTess(&subTess1);
                    printf("freeTess -> status=%d\n", status);
                }

                new_data = 0;

            /* 'v' - toggle Triangle visibility */
            } else if (*state == 'v') {
                printf("==> Option 'v' \n");
                fprintf(fpdump, "v\n");

                for (itri = 0; itri < tess.ntri; itri++) {
                    if (tess.ttyp[itri] & TRI_ACTIVE) {
                        if (tess.ttyp[itri] &   TRI_VISIBLE) {
                            tess.ttyp[itri] &= ~TRI_VISIBLE;
                        } else {
                            tess.ttyp[itri] |=  TRI_VISIBLE;
                        }
                    }
                }

                new_data = 1;

            /* 'w' - plot facet colors */
            } else if (*state == 'w') {
                if (saveit != NULL) FREE(saveit);
                saveit = NULL;

            /* 'W' - write stl file */
            } else if (*state == 'W') {
                printf("==> Option 'W' \n");
                fprintf(fpdump, "W\n");

                itype = getInt("Enter 0=asc, 1-bin, 2=tris: ");

                getStr("Enter filename: ", filename);
                trim(filename);

                if        (itype == 0) {
                    writeStlAscii(&tess, filename);
                } else if (itype == 1) {
                    writeStlBinary(&tess, filename);
                } else if (itype == 2) {
                    writeTriAscii(&tess, filename);
                }

                new_data = 1;

            /* 'x' - plot x coordinates */
            } else if (*state == 'x') {
                if (saveit != NULL) FREE(saveit);
                saveit = NULL;

            /* 'X' - extend loop adjacent to Point */
            } else if (*state == 'X') {
                if (numarg >= 0) {
                    ipnt   = numarg;
                    numarg = -1;
                } else {
                    ipnt   = findPoint(*xscr, *yscr);
                }

                if (ipnt >= 0) {
                    printf("==> Option 'X' (ipnt=%6d)\n", ipnt);
                    fprintf(fpdump, "%6d X\n", ipnt);

                    itype = getInt("Enter 1=x, 2=y, 3=z: ");
                    val   = getDbl("Enter value: ");

                    extendLoop(&tess, ipnt, itype, val);
                    new_data = 1;
                } else {
                    printf("ERR: No Point found\a\n");
                }

            /* 'y' - plot y coordinates */
            } else if (*state == 'y') {
                if (saveit != NULL) FREE(saveit);
                saveit = NULL;

            /* 'z' - plot z coordinates */
            } else if (*state == 'z') {
                if (saveit != NULL) FREE(saveit);
                saveit = NULL;

            /* '0' - append "0" to numarg */
            } else if (*state == '0') {
                if (numarg >= 0) {
                    numarg = 0 + numarg * 10;
                } else {
                    numarg = 0;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '1' - append "1" to numarg */
            } else if (*state == '1') {
                if (numarg >= 0) {
                    numarg = 1 + numarg * 10;
                } else {
                    numarg = 1;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '2' - append "2" to numarg */
            } else if (*state == '2') {
                if (numarg >= 0) {
                    numarg = 2 + numarg * 10;
                } else {
                    numarg = 2;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '3' - append "3" to numarg */
            } else if (*state == '3') {
                if (numarg >= 0) {
                    numarg = 3 + numarg * 10;
                } else {
                    numarg = 3;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '4' - append "4" to numarg */
            } else if (*state == '4') {
                if (numarg >= 0) {
                    numarg = 4 + numarg * 10;
                } else {
                    numarg = 4;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '5' - append "5" to numarg */
            } else if (*state == '5') {
                if (numarg >= 0) {
                    numarg = 5 + numarg * 10;
                } else {
                    numarg = 5;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '6' - append "6" to numarg */
            } else if (*state == '6') {
                if (numarg >= 0) {
                    numarg = 6 + numarg * 10;
                } else {
                    numarg = 6;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '7' - append "7" to numarg */
            } else if (*state == '7') {
                if (numarg >= 0) {
                    numarg = 7 + numarg * 10;
                } else {
                    numarg = 7;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '8' - append "8" to numarg */
            } else if (*state == '8') {
                if (numarg >= 0) {
                    numarg = 8 + numarg * 10;
                } else {
                    numarg = 8;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '9' - append "9" to numarg */
            } else if (*state == '9') {
                if (numarg >= 0) {
                    numarg = 9 + numarg * 10;
                } else {
                    numarg = 9;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

           /* <bksp> - erase last digit of numarg */
            } else if (*state == 65288) {
                if (numarg > 0) {
                    numarg = numarg / 10;
                } else {
                    numarg = -1;
                }
                if (script == NULL) {
                    printf("    numarg = %d\n", numarg);
                }

            /* '>' - write current viewpoint */
            } else if (*state == '>') {
                printf("==> Option '>' \n");
                if (numarg >= 0) {
                    sprintf(filename, "ViewMatrix%d.dat", numarg);
                    numarg = -1;
                } else {
                    sprintf(filename, "ViewMatrix.dat");
                }

                fp = fopen(filename, "w");
                if (fp != NULL) {
                    fprintf(fp, "%f %f %f %f\n", gv_xform[0][0], gv_xform[1][0],
                                                 gv_xform[2][0], gv_xform[3][0]);
                    fprintf(fp, "%f %f %f %f\n", gv_xform[0][1], gv_xform[1][1],
                                                 gv_xform[2][1], gv_xform[3][1]);
                    fprintf(fp, "%f %f %f %f\n", gv_xform[0][2], gv_xform[1][2],
                                                 gv_xform[2][2], gv_xform[3][2]);
                    fprintf(fp, "%f %f %f %f\n", gv_xform[0][3], gv_xform[1][3],
                                                 gv_xform[2][3], gv_xform[3][3]);
                    fclose(fp);

                    printf("    Current view transformation saved\n");
                } else {
                    printf("ERR: Could not open file\a\n");
                }

            /* '<' - read saved viewpoint */
            } else if (*state == '<') {
                printf("==> Option '<' \n");
                if (numarg >= 0) {
                    sprintf(filename, "ViewMatrix%d.dat", numarg);
                    numarg = -1;
                } else {
                    sprintf(filename, "ViewMatrix.dat");
                }

                fp = fopen(filename, "r");
                if (fp != NULL) {
                    fscanf(fp, "%f%f%f%f", &(gv_xform[0][0]), &(gv_xform[1][0]),
                                           &(gv_xform[2][0]), &(gv_xform[3][0]));
                    fscanf(fp, "%f%f%f%f", &(gv_xform[0][1]), &(gv_xform[1][1]),
                                           &(gv_xform[2][1]), &(gv_xform[3][1]));
                    fscanf(fp, "%f%f%f%f", &(gv_xform[0][2]), &(gv_xform[1][2]),
                                           &(gv_xform[2][2]), &(gv_xform[3][2]));
                    fscanf(fp, "%f%f%f%f", &(gv_xform[0][3]), &(gv_xform[1][3]),
                                           &(gv_xform[2][3]), &(gv_xform[3][3]));
                    fclose(fp);
                } else {
                    printf("ERR: Could not open file\a\n");
                }

                new_data = 1;

            /* <home> - original viewpoint */
            } else if (*state == 65360) {
                getModelSize(box);

                size = 0.5 * sqrt(  pow(box[3] - box[0], 2)
                                  + pow(box[4] - box[1], 2)
                                  + pow(box[5] - box[2], 2));

                gv_xform[0][0] = +1 / size;
                gv_xform[1][0] =  0;
                gv_xform[2][0] =  0;
                gv_xform[3][0] = -(box[0] + box[3]) / 2 / size;
                gv_xform[0][1] =  0;
                gv_xform[1][1] = +1 / size;
                gv_xform[2][1] =  0;
                gv_xform[3][1] = -(box[1] + box[4]) / 2 / size;
                gv_xform[0][2] =  0;
                gv_xform[1][2] =  0;
                gv_xform[2][2] = +1 / size;
                gv_xform[3][2] = -(box[2] + box[5]) / 2 / size;
                gv_xform[0][3] =  0;
                gv_xform[1][3] =  0;
                gv_xform[2][3] =  0;
                gv_xform[3][3] =  1;

            /* & - toggle flying mode */
            } else if (*state == '&') {
                if (fly_mode == 0) {
                    printf("==> Option '&' (turning fly mode on)\n");
                    fly_mode = 1;
                } else {
                    printf("==> Option '&' (turning fly mode off)\n");
                    fly_mode = 0;
                }

            /* <left> - rotate viewpoint left 30 deg or translate left */
            } else if (*state == 65361) {
                if (fly_mode == 0) {
                    cosrot = cos( PI/6);
                    sinrot = sin( PI/6);

                    for (i = 0; i < 4; i++) {
                        temp0 = gv_xform[i][0];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][0] = cosrot * temp0 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp0 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][0] -= 0.5;
                }

            /* <up> - rotate viewpoint up 30 deg or translate up */
            } else if (*state == 65362) {
                if (fly_mode == 0) {
                    cosrot = cos(-PI/6);
                    sinrot = sin(-PI/6);

                    for (i = 0; i < 4; i++) {
                        temp1 = gv_xform[i][1];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][1] = cosrot * temp1 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp1 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][1] += 0.5;
                }

            /* <rite> - rotate viewpoint right 30 deg or translate right */
            } else if (*state == 65363) {
                if (fly_mode == 0) {
                    cosrot = cos(-PI/6);
                    sinrot = sin(-PI/6);

                    for (i = 0; i < 4; i++) {
                        temp0 = gv_xform[i][0];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][0] = cosrot * temp0 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp0 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][0] += 0.5;
                }

            /* <down> - rotate viewpoint down 30 deg */
            } else if (*state == 65364) {
                if (fly_mode == 0) {
                    cosrot = cos(PI/6);
                    sinrot = sin(PI/6);

                    for (i = 0; i < 4; i++) {
                        temp1 = gv_xform[i][1];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][1] = cosrot * temp1 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp1 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][1] -= 0.5;
                }

            /* <PgUp> - zoom in */
            } else if (*state == 65365) {
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 3; j++) {
                        gv_xform[i][j] *= 2;
                    }
                }

            /* <PgDn- zoom out */
            } else if (*state == 65366) {
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 3; j++) {
                        gv_xform[i][j] /= 2;
                    }
                }

            /* '%' - comment */
            } else if (*state == '%') {
                printf("==> Option '%%'\n");

                getStr("Enter comment: ", comment);
                trim(comment);

                fprintf(fpdump, "%% %s\n", comment);

            /* '$' - start or stop reading script file */
            } else if (*state == '$') {
                printf("==> Option '$' \n");
                if (script == NULL) {
                    printf("==> Enter script name: ");
                    scanf("%s", scriptName);

                    printf("    Opening script file \"%s\" ...\n", scriptName);

                    script = fopen(scriptName, "r");
                    if (script == NULL) {
                        printf("ERR: unsuccessful\a\n");
                    }
                } else {
                    fclose(script);
                    printf("    Closing script file\n");

                    script = NULL;
                }

            /* '?' - print out help */
            } else if (*state == '?') {
                printf("==> Option '?' \n");
                printf("                                             \n");
                printf("t  # toggle Triangle mark                    \n");
                printf("c  # color Triangle and neighbors            \n");
                printf("m    make topology                           \n");
                printf("q    query topology                          \n");
                printf("S  # split Edge                              \n");
                printf("e  # mark Edge                               \n");
                printf("                                             \n");
                printf("i  - invert all visible Triangles            \n");
                printf("k  # kill Triangle and neighbors             \n");
                printf("d  # delete a Triangle                       \n");
                printf("u  - undelete last Triangle                  \n");
                printf("b  # bridge to marked Triangle               \n");
                printf("f  # fill loop adjacent to Triangle          \n");
                printf("                                             \n");
                printf("G  - toggle grid visibility                  \n");
                printf("h  # hide Triangle and neighbors             \n");
                printf("s  - show all Triangles                      \n");
                printf("v  - toggle Triangle visibility              \n");
                printf("                                             \n");
                printf("p  # toggle Point mark                       \n");
                printf("l  # link Points                             \n");
                printf("L  - create Links between colors             \n");
                printf("j  # join Points                             \n");
                printf("                                             \n");
                printf("g  # generate Grid for Triangle and neighbors\n");
                printf("W  - write stl file                          \n");
                printf("                                             \n");
                printf("0-9  add digit to numeric argument           \n");
                printf("Bksp erase last digit from numeric argument  \n");
                printf("                                             \n");
                printf("w  - plot facet colors                       \n");
                printf("x  - plot x coordinates                      \n");
                printf("y  - plot y coordinates                      \n");
                printf("z  - plot z coordinates                      \n");
                printf("                                             \n");
                printf("Home original viewpoint                      \n");
                printf("Left rotate viewpoint or transleft left      \n");
                printf("Rite rotate viewpoint or translate rite      \n");
                printf("Up   rotate viewpoint or translate up        \n");
                printf("Down rotate viewpoint or translate down      \n");
                printf("PgUP zoom in                                 \n");
                printf("PgDn zoom out                                \n");
                printf(">  # write current viewpoint                 \n");
                printf("<  # read  saved   viewpoint                 \n");
                printf("&  - toggle flying mode                      \n");
                printf("                                             \n");
                printf("%%  - comment                                 \n");
                printf("$  - start or stop reading script file       \n");
                printf("?  - print help                              \n");
            }
        }

        /* repeat as long as we are in a script */
    } while ((script != NULL) && (*type == KeyPress));
}


/******************************************************************************/

static int
getInt(char   *prompt)                  /* (in)  user prompt */
{
    int    answer;                      /* (out) value entered by user */

//    ROUTINE(getInt);

    /*----------------------------------------------------------------*/

    if (script == NULL) {
        printf("%s ", prompt);
        scanf("%d", &answer);
    } else {
        fscanf(script, "%d", &answer);

        printf("==> %s %d\n", prompt, answer);
    }

    fprintf(fpdump, "%d\n", answer);

//cleanup:
    return answer;
}


/******************************************************************************/

static double
getDbl(char   *prompt)                  /* (in)  user prompt */
{
    double answer;                      /* (out) value entered by user */

//    ROUTINE(getDbl);

    /*----------------------------------------------------------------*/

    if (script == NULL) {
        printf("%s ", prompt);
        scanf("%lf", &answer);
    } else {
        fscanf(script, "%lf", &answer);

        printf("==> %s %f\n", prompt, answer);
    }

    fprintf(fpdump, "%f\n", answer);

//cleanup:
    return answer;
}


/******************************************************************************/

static void
getStr(char   *prompt,                  /* (in)  use prompt */
       char   *answer)                  /* (out) value enetred by user */
{
//    ROUTINE(getStr);

    /*----------------------------------------------------------------*/

    if (script == NULL) {
        printf("%s ", prompt); fflush(0);
        scanf("%s", answer);
    } else {
        fscanf(script, "%s", answer);

        printf("==> %s %s\n", prompt, answer);
    }

    fprintf(fpdump, "%s", answer);

//cleanup:
    return;
}


/******************************************************************************/

//$$$static void
//$$$plot3D(int   *ifunct,
//$$$       void  *myTessP,
//$$$       void  *iviewP,
//$$$/*@unused@*/void  *a2,
//$$$/*@unused@*/void  *a3,
//$$$/*@unused@*/void  *a4,
//$$$/*@unused@*/void  *a5,
//$$$/*@unused@*/void  *a6,
//$$$/*@unused@*/void  *a7,
//$$$/*@unused@*/void  *a8,
//$$$/*@unused@*/void  *a9,
//$$$       float *scale,
//$$$       char  *text,
//$$$       int   textlen)
//$$${
//$$$    tess_T    *myTess = (tess_T *) myTessP;
//$$$    int       *iview  = (int    *) iviewP;
//$$$
//$$$    int      itri, ipnt, i;
//$$$    float    x4, y4, z4, fdum;
//$$$    double   xmin, xmax, ymin, ymax, zmin, zmax;
//$$$
//$$$    int       icirc  = GR_CIRCLE;
//$$$    int       iblack = GR_BLACK;
//$$$    int       ired   = GR_RED;
//$$$    int       ione   = 1;
//$$$
//$$$    /* --------------------------------------------------------------- */
//$$$
//$$$    /* ---------- return scales ----------*/
//$$$    if (*ifunct == 0) {
//$$$        xmin = myTess->xyz[0];
//$$$        xmax = myTess->xyz[0];
//$$$        ymin = myTess->xyz[1];
//$$$        ymax = myTess->xyz[1];
//$$$        zmin = myTess->xyz[2];
//$$$        zmax = myTess->xyz[2];
//$$$
//$$$        for (ipnt = 0; ipnt < myTess->npnt; ipnt++) {
//$$$            if (myTess->xyz[3*ipnt  ] < xmin) xmin = myTess->xyz[3*ipnt  ];
//$$$            if (myTess->xyz[3*ipnt  ] > xmax) xmax = myTess->xyz[3*ipnt  ];
//$$$            if (myTess->xyz[3*ipnt+1] < ymin) ymin = myTess->xyz[3*ipnt+1];
//$$$            if (myTess->xyz[3*ipnt+1] > ymax) ymax = myTess->xyz[3*ipnt+1];
//$$$            if (myTess->xyz[3*ipnt+2] < zmin) zmin = myTess->xyz[3*ipnt+2];
//$$$            if (myTess->xyz[3*ipnt+2] > zmax) zmax = myTess->xyz[3*ipnt+2];
//$$$        }
//$$$
//$$$        if (*iview == 1) {
//$$$            scale[0] = xmin - 0.05 * (xmax - xmin);
//$$$            scale[1] = xmax + 0.05 * (xmax - xmin);
//$$$            scale[2] = ymin - 0.05 * (ymax - ymin);
//$$$            scale[3] = ymax + 0.05 * (ymax - ymin);
//$$$        } else if (*iview == 2) {
//$$$            scale[0] = ymin - 0.05 * (ymax - ymin);
//$$$            scale[1] = ymax + 0.05 * (ymax - ymin);
//$$$            scale[2] = zmin - 0.05 * (zmax - zmin);
//$$$            scale[3] = zmax + 0.05 * (zmax - zmin);
//$$$        } else {
//$$$            scale[0] = zmin - 0.05 * (zmax - zmin);
//$$$            scale[1] = zmax + 0.05 * (zmax - zmin);
//$$$            scale[2] = xmin - 0.05 * (xmax - xmin);
//$$$            scale[3] = xmax + 0.05 * (xmax - xmin);
//$$$        }
//$$$
//$$$        strncpy(text, "nextView                        ", textlen-1);
//$$$
//$$$    /* ---------- plot image ---------- */
//$$$    } else if (*ifunct == 1) {
//$$$
//$$$        /* Triangles */
//$$$        for (itri = 0; itri < myTess->ntri; itri++) {
//$$$            ipnt = myTess->trip[3*itri+2];
//$$$
//$$$            x4 = myTess->xyz[3*ipnt  ];
//$$$            y4 = myTess->xyz[3*ipnt+1];
//$$$            z4 = myTess->xyz[3*ipnt+2];
//$$$
//$$$            if        (*iview == 1) {
//$$$                grmov2_(&x4, &y4);
//$$$            } else if (*iview == 2) {
//$$$                grmov2_(&y4, &z4);
//$$$            } else {
//$$$                grmov2_(&z4, &x4);
//$$$            }
//$$$
//$$$            for (i = 0; i < 3; i++) {
//$$$                ipnt = myTess->trip[3*itri+i];
//$$$
//$$$                x4 = myTess->xyz[3*ipnt  ];
//$$$                y4 = myTess->xyz[3*ipnt+1];
//$$$                z4 = myTess->xyz[3*ipnt+2];
//$$$
//$$$                if        (*iview == 1) {
//$$$                    grdrw2_(&x4, &y4);
//$$$                } else if (*iview == 2) {
//$$$                    grdrw2_(&y4, &z4);
//$$$                } else {
//$$$                    grdrw2_(&z4, &x4);
//$$$                }
//$$$            }
//$$$        }
//$$$
//$$$        /* hanging Points (as red circles) */
//$$$        grcolr_(&ired);
//$$$        for (itri = 0; itri < myTess->ntri; itri++) {
//$$$            for (i = 0; i < 3; i++) {
//$$$                if (myTess->trit[3*itri+i] < 0) {
//$$$                    ipnt = myTess->trip[3*itri+(i+1)%3];
//$$$
//$$$                    x4 = myTess->xyz[3*ipnt  ];
//$$$                    y4 = myTess->xyz[3*ipnt+1];
//$$$                    z4 = myTess->xyz[3*ipnt+2];
//$$$
//$$$                    if        (*iview == 1) {
//$$$                        grmov2_(&x4, &y4);
//$$$                    } else if (*iview == 2) {
//$$$                        grmov2_(&y4, &z4);
//$$$                    } else {
//$$$                        grmov2_(&z4, &x4);
//$$$                    }
//$$$                    grsymb_(&icirc);
//$$$
//$$$                    ipnt = myTess->trip[3*itri+(i+2)%3];
//$$$
//$$$                    x4 = myTess->xyz[3*ipnt  ];
//$$$                    y4 = myTess->xyz[3*ipnt+1];
//$$$                    z4 = myTess->xyz[3*ipnt+2];
//$$$
//$$$                    if        (*iview == 1) {
//$$$                        grmov2_(&x4, &y4);
//$$$                    } else if (*iview == 2) {
//$$$                        grmov2_(&y4, &z4);
//$$$                    } else {
//$$$                        grmov2_(&z4, &x4);
//$$$                    }
//$$$                    grsymb_(&icirc);
//$$$                }
//$$$            }
//$$$        }
//$$$        grcolr_(&iblack);
//$$$
//$$$    /* ---------- options ---------- */
//$$$    } else if (*ifunct == -14) {
//$$$        *iview = *iview + 1;
//$$$        if (*iview > 3) *iview = 1;
//$$$
//$$$        if        (*iview == 1) {
//$$$            grvalu_("LABLGR", &ione, &fdum, "~x~y~ ", strlen("LABLGR"), strlen("~x~y~ "));
//$$$        } else if (*iview == 2) {
//$$$            grvalu_("LABLGR", &ione, &fdum, "~y~z~ ", strlen("LABLGR"), strlen("~y~z~ "));
//$$$        } else {
//$$$            grvalu_("LABLGR", &ione, &fdum, "~z~x~ ", strlen("LABLGR"), strlen("~z~x~ "));
//$$$        }
//$$$
//$$$        grscpt_(&ione, "O", strlen("O"));
//$$$
//$$$    } else {
//$$$        printf("Illegal option selected\a\n");
//$$$    }
//$$$    return;
//$$$}


/******************************************************************************/

#ifdef GRAFIC
static void
plot2D(int   *ifunct,
       void  *myTessP,
       void  *igridP,
       void  *isymbP,
/*@unused@*/void  *a3,
/*@unused@*/void  *a4,
/*@unused@*/void  *a5,
/*@unused@*/void  *a6,
/*@unused@*/void  *a7,
/*@unused@*/void  *a8,
/*@unused@*/void  *a9,
       float *scale,
       char  *text,
       int   textlen)
{
    tess_T    *myTess = (tess_T *) myTessP;
    int       *igrid  = (int    *) igridP;
    int       *isymb  = (int    *) isymbP;

    int      itri, ipnt, i;
    float    u4, v4, upoly[3], vpoly[3], area;
    double   umin, umax, vmin, vmax;

    int       icirc   = GR_CIRCLE;
    int       isquare = GR_SQUARE;
    int       iblack  = GR_BLACK;
    int       ired    = GR_RED;
    int       iblue   = GR_BLUE;
    int       ione    = 1;
    int       ithree  = 3;

    /* --------------------------------------------------------------- */

    /* ---------- return scales ----------*/
    if (*ifunct == 0) {
        umin = myTess->uv[0];
        umax = myTess->uv[0];
        vmin = myTess->uv[1];
        vmax = myTess->uv[1];

        for (ipnt = 0; ipnt < myTess->npnt; ipnt++) {
            if (myTess->uv[2*ipnt  ] < umin) umin = myTess->uv[2*ipnt  ];
            if (myTess->uv[2*ipnt  ] > umax) umax = myTess->uv[2*ipnt  ];
            if (myTess->uv[2*ipnt+1] < vmin) vmin = myTess->uv[2*ipnt+1];
            if (myTess->uv[2*ipnt+1] > vmax) vmax = myTess->uv[2*ipnt+1];
        }

        scale[0] = umin - 0.05 * (umax - umin);
        scale[1] = umax + 0.05 * (umax - umin);
        scale[2] = vmin - 0.05 * (vmax - vmin);
        scale[3] = vmax + 0.05 * (vmax - vmin);

        strncpy(text, "Grid Symb                       ", textlen-1);

    /* ---------- plot image ---------- */
    } else if (*ifunct == 1) {

        /* negative areas filled in red */
        grcolr_(&ired);
        for (itri = 0; itri < myTess->ntri; itri++) {
            upoly[0] = myTess->uv[2*myTess->trip[3*itri  ]  ];
            vpoly[0] = myTess->uv[2*myTess->trip[3*itri  ]+1];
            upoly[1] = myTess->uv[2*myTess->trip[3*itri+1]  ];
            vpoly[1] = myTess->uv[2*myTess->trip[3*itri+1]+1];
            upoly[2] = myTess->uv[2*myTess->trip[3*itri+2]  ];
            vpoly[2] = myTess->uv[2*myTess->trip[3*itri+2]+1];

            area = (upoly[1] - upoly[0]) * (vpoly[2] - vpoly[0])
                 - (vpoly[1] - vpoly[0]) * (upoly[2] - upoly[0]);
            if (area < 0) {
                grfil2_(upoly, vpoly, &ithree, &ired);

                u4 = (upoly[0] + upoly[1] + upoly[2]) / 3;
                v4 = (vpoly[0] + vpoly[1] + vpoly[2]) / 3;
                grmov2_(&u4, &v4);
                grsymb_(&isquare);
            }
        }
        grcolr_(&iblack);

        /* Triangles */
        if (*igrid > 0) {
            for (itri = 0; itri < myTess->ntri; itri++) {
                ipnt = myTess->trip[3*itri+2];

                u4 = myTess->uv[2*ipnt  ];
                v4 = myTess->uv[2*ipnt+1];
                grmov2_(&u4, &v4);

                for (i = 0; i < 3; i++) {
                    ipnt = myTess->trip[3*itri+i];

                    u4 = myTess->uv[2*ipnt  ];
                    v4 = myTess->uv[2*ipnt+1];
                    grdrw2_(&u4, &v4);
                }
            }
        }

        /* hanging Points (as blue circles) */
        if (*isymb > 0) {
            grcolr_(&iblue);
            for (itri = 0; itri < myTess->ntri; itri++) {
                for (i = 0; i < 3; i++) {
                    if (myTess->trit[3*itri+i] < 0) {
                        ipnt = myTess->trip[3*itri+(i+1)%3];

                        u4 = myTess->uv[2*ipnt  ];
                        v4 = myTess->uv[2*ipnt+1];
                        grmov2_(&u4, &v4);
                        grsymb_(&icirc);

                        ipnt = myTess->trip[3*itri+(i+2)%3];

                        u4 = myTess->uv[2*ipnt  ];
                        v4 = myTess->uv[2*ipnt+1];
                        grmov2_(&u4, &v4);
                        grsymb_(&icirc);
                    }
                }
            }
            grcolr_(&iblack);
        }

    /* ---------- options ---------- */
    } else if (*ifunct == -7) {
        *igrid = 1 - *igrid;

        grscpt_(&ione, "R", strlen("R"));
    } else if (*ifunct == -19) {
        *isymb = 1 - *isymb;

        grscpt_(&ione, "R", strlen("R"));
    } else {
        printf("Illegal option selected\a\n");
    }
    return;
}
#endif


/******************************************************************************/

