/*
 ************************************************************************
 *                                                                      *
 * eggPDT2 --- parametric Delaunay triangulator                         *
 *                                                                      *
 *           Written by John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2018/2022  John F. Dannenhoffer, III (Syracuse University)
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

/*
 * Note: This code is for demonstration purposes.  In order to make
 *       it as readable as possible, it has NOT been optimized/organized
 *       for performance.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include "egg.h"

/* enable this to show details during eggGenerate */
//#define DEBUG2 1

/*
 ***********************************************************************
 *                                                                     *
 * Definitions                                                         *
 *                                                                     *
 ***********************************************************************
 */

/* "Point" the locations in physical space at which the Grid is defined

   the Points are ordered so that the Boundary Points come first,
   followed by the ConvexHull Points (if they exist), followed by the
   inserted field Points
*/

typedef struct {
    double        u;                    /* u-coordinate */
    double        v;                    /* v-coordinate */
    double        s;                    /* minimum local spacing */
    int           p[3];                 /* parent Points (bias-0) */
                                        /*           -1 if Boundary, -2 if ConvexHull */
} pnt_T;

/* "Triangle" the 3-sided cells that fill the space that is bounded by the Boundary

                 n[2]
                 / \
                /   \
       t[1]    /     \     t[0]
              /       \
             /         \
            /           \
           /             \
         n[0]------------n[1]

                t[2]
*/

typedef struct {
    int           p[3];                 /* Point    indices (bias-0) */
    int           t[3];                 /* Triangle indices (bias-0) */
    double        uc;                   /* u-coordinate of circumcircle */
    double        vc;                   /* v-coordinate of circumcircle */
    double        rr;                   /* radius^2     of circumcircle */
                                        /*     if rr<0, then Triangle is deleted */
} tri_T;

/* "Grid" a structure containing Points, Triangles, and other data */

typedef struct {
    int           mpnt;                 /* maximum number of Points */
    int           npnt;                 /* number of Points */
    int           nbnd;                 /* number of Boundary Points */
    pnt_T         *pnt;                 /* array  of Points */
    double        *uv;
    int           *p;

    int           mtri;                 /* maximum number of Triangles */
    int           ntri;                 /* number of Triangles */
    tri_T         *tri;                 /* array  of Triangles */
    int           *tris;
} grid_T;

/* "Boundary" a list of Points that enclose the domain to be gridded.  the
       Boundary should traverse the domain in a counterclockwise direction
       (ie, the domain should be to the left of the Boundary */

/* "Loop" a sequence of Triangle edges that bound a hole that is left when
       Triangles are deleted */

/* "ConvexHull" an initial set of Triangles that fully enclose all Boundary Points */

#define  EPS06               1.0e-6
#define  MIN(A,B)            (((A) < (B)) ? (A) : (B))
#define  MAX(A,B)            (((A) < (B)) ? (B) : (A))
#define  SQR(A)              ((A) * (A))

/*
 ***********************************************************************
 *                                                                     *
 * declarations                                                        *
 *                                                                     *
 ***********************************************************************
 */

/* declarations for support routines defined below */
static int    addToLoop(int ip0, int ip1, int iloop[], int *nloop);
static double computeArea(grid_T *grid, int itri);
static int    computeParameters(grid_T *grid, int itri);
static int    createTriangle(grid_T *grid, int ip0, int ip1, int ip2);
static double distance(grid_T *grid, int ip0, int ip1);
static int    insertPoint(grid_T *grid, int ipnt);
static int    intersect(grid_T *grid, int ip0, int ip1, int ip2, int ip3);
static int    newPoint(grid_T *grid);
static int    newTriangle(grid_T *grid);
static int    swapDiagonals(grid_T *grid, int i1, int ib);

#ifdef  GRAFIC
#include "grafic.h"
static void plotGrid(int*, void*, void*, void*, void*, void*,
                           void*, void*, void*, void*, void*, float*, char*, int);
#endif


/*
 ***********************************************************************
 *                                                                     *
 * eggGenerate - generate a Grid                                       *
 *                                                                     *
 ***********************************************************************
 */
int
eggGenerate(double uv[],                /* (in)  array of coordinates for Boundary Points */
                                        /* (out) array of coordinates for Grid     Points */
            int    lup[],               /* (in)  number of Boundary Points in each loop */
            void   **gridP)             /* (out) pointer to Grid (which is a grid_T*) */
{
    int    status = SUCCESS;            /* (out) default return */

    int    nlup, ilast, jpnt, ipnt, ibeg, iend, ilup, itri, ip0, ip1, ip2, it0, it1;
    int    ipass, nsave, ii, isid, jtri, jsid, nchange;
    int    ntrinew, *itrinew=NULL, *valence=NULL, nextpnt;
    double umin, umax, vmin, vmax, dist2, domsiz;
    double alfa=0.80;                   /* tolerance on how close a new Point can be
                                           to the vertices of a current Triangle */
    double beta=5.00;                   /* tolerance on how close a new Point can be
                                           from any other Point added in the current pass */
    grid_T *grid;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("eggGenerate()\n");
#endif

    /* create the Grid structure */
    grid = (grid_T *) malloc(sizeof(grid_T));
    if (grid == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }

#ifdef GRAFIC
    {
        int io_kbd=5, io_scr=6;

        grinit_(&io_kbd, &io_scr, "eggGenerate", strlen("eggGenerate"));
    }
#endif

    /* return a pointer to the Grid structure */
    *gridP = (void *) grid;

    /* create initial Point and Triangle tables */
    grid->npnt = 0;
    grid->mpnt = 1000;
    grid->pnt  = (pnt_T *) malloc(grid->mpnt*sizeof(pnt_T));
    if (grid->pnt == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }
    grid->uv = NULL;
    grid->p  = NULL;

    grid->ntri = 0;
    grid->mtri = 1000;
    grid->tri  = (tri_T *) malloc(grid->mtri*sizeof(tri_T));
    if (grid->tri == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }
    grid->tris = NULL;

    /* determine (and store) the number of Boundary Points */
    nlup        = 0;
    grid->nbnd  = 0;

#ifndef __clang_analyzer__
    ilast       = 1;
    while(ilast > 0) {
        ilast       = lup[nlup];
        grid->nbnd  = grid->nbnd + ilast;
        nlup++;
    }
#endif

    if (grid->nbnd < 3) {
        status = NUMBER_OF_POINT_MISMATCH;
        goto cleanup;
    }

    /* store Boundary Points */
    for (jpnt = 0; jpnt < grid->nbnd; jpnt++) {
        ipnt = newPoint(grid);
        grid->pnt[ipnt].u    = uv[2*jpnt  ];
        grid->pnt[ipnt].v    = uv[2*jpnt+1];
        grid->pnt[ipnt].s    = 0;
        grid->pnt[ipnt].p[0] = -1;   // Boundary
        grid->pnt[ipnt].p[1] = -1;
        grid->pnt[ipnt].p[2] = -1;
    }
#ifdef DEBUG2
    printf("Initialization                  npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

    /* compute average spacing to each Point's two neighbors, taking
       care to wrap around approriately at the beginning and end of
       the Boundary */
    iend = -1;
    for (ilup = 0; ilup < nlup; ilup++) {
        ibeg = iend + 1;
        iend = ibeg + lup[ilup] - 1;

        ip0 = iend;
        for (ip1 = ibeg; ip1 <= iend; ip1++) {
            dist2 = distance(grid, ip0, ip1) / 2;

            grid->pnt[ip0].s += dist2;
            grid->pnt[ip1].s += dist2;

            ip0 = ip1;
        }
    }
#ifdef DEBUG2
    printf("Spacings computed:              npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

    /* find the extrema of the configuration Points */
    umin = grid->pnt[0].u;
    umax = grid->pnt[0].u;
    vmin = grid->pnt[0].v;
    vmax = grid->pnt[0].v;

    for (ipnt = 1; ipnt < grid->nbnd; ipnt++) {
        umin = MIN(umin, grid->pnt[ipnt].u);
        umax = MAX(umax, grid->pnt[ipnt].u);
        vmin = MIN(vmin, grid->pnt[ipnt].v);
        vmax = MAX(vmax, grid->pnt[ipnt].v);
    }
    domsiz = MAX(umax-umin, vmax-vmin);

    /* set up a ConvexHull and associated initial Triangles that will
       surround the entire configuration */

    /* create the ConvexHull Points (notice that a right trapezoid is
       used so that the initial triangulation will be uniquely delaunay) */

    /* note that the spacing (%s) associated with the ConvexHull Points
       is set very large so that the Triangles that are "outside" the
       configuration will not be refined by the Point insertion
       algorithm */
    ipnt                 = newPoint(grid);
    grid->pnt[ipnt].u    =  1.05 * umin - 0.05 * umax;
    grid->pnt[ipnt].v    =  1.05 * vmin - 0.05 * vmax;
    grid->pnt[ipnt].s    =  10.0 * domsiz;
    grid->pnt[ipnt].p[0] = -2;    // ConvexHull

    ipnt                 = newPoint(grid);
    grid->pnt[ipnt].u    = -0.10 * umin + 1.10 * umax;
    grid->pnt[ipnt].v    =  1.05 * vmin - 0.05 * vmax;
    grid->pnt[ipnt].s    =  10.0 * domsiz;
    grid->pnt[ipnt].p[0] = -2;    // ConvexHull

    ipnt                 = newPoint(grid);
    grid->pnt[ipnt].u    = -0.05 * umin + 1.05 * umax;
    grid->pnt[ipnt].v    = -0.05 * vmin + 1.05 * vmax;
    grid->pnt[ipnt].s    =  10.0 * domsiz;
    grid->pnt[ipnt].p[0] = -2;    // ConvexHull

    ipnt                 = newPoint(grid);
    grid->pnt[ipnt].u    =  1.05 * umin - 0.05 * umax;
    grid->pnt[ipnt].v    = -0.05 * vmin + 1.05 * vmax;
    grid->pnt[ipnt].s    =  10.0 * domsiz;
    grid->pnt[ipnt].p[0] = -2;    // ConvexHull

    /* create the initial delaunay triangulation */
    status = createTriangle(grid, grid->npnt-4, grid->npnt-3, grid->npnt-2);
    if (status < SUCCESS) goto cleanup;

    status = createTriangle(grid, grid->npnt-2, grid->npnt-1, grid->npnt-4);
    if (status < SUCCESS) goto cleanup;

#ifdef DEBUG2
    printf("ConvexHull and base Triangles:  npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

#ifdef GRAFIC
    {
        int  indgr    = 1+2+4+16+64;
        char pltitl[] = "~u~v~ConvexHull and base Triangles";

        grctrl_(plotGrid, &indgr, pltitl,
                (void*)(grid),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                strlen(pltitl));
    }
#endif

    /* add all Points associated with the Boundary to the triangulation */
    for (ipnt = 0; ipnt < grid->nbnd; ipnt++) {
        status = insertPoint(grid, ipnt);
        if (status < SUCCESS) goto cleanup;
    }
#ifdef DEBUG2
    printf("Boundary Points added:          npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

#ifdef GRAFIC
    {
        int  indgr    = 1+2+4+16+64;
        char pltitl[] = "~u~v~Boundary Points added";

        grctrl_(plotGrid, &indgr, pltitl,
                (void*)(grid),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                strlen(pltitl));
    }
#endif

    /* add Points in the field so that the resulting Triangles will be
       shaped well and will be reasonably regular.  this will be done
       in multiple passes until no Points are added in the last pass */
    nsave = grid->npnt - 1;             // number of Points at start of this pass
    for (ipass = 0; ipass < grid->nbnd*10; ipass++) {
        if (grid->npnt <= nsave) break;

        nsave = grid->npnt;

        /* visit each Triangle and tentatively place a Point at its centroid */
        for (itri = 0; itri < grid->ntri; itri++) {
            if (grid->tri[itri].rr < 0) continue;

            ip0 = grid->tri[itri].p[0];
            ip1 = grid->tri[itri].p[1];
            ip2 = grid->tri[itri].p[2];

            ipnt = newPoint(grid);
            grid->pnt[ipnt].u    = (grid->pnt[ip0].u + grid->pnt[ip1].u + grid->pnt[ip2].u) / 3.0;
            grid->pnt[ipnt].v    = (grid->pnt[ip0].v + grid->pnt[ip1].v + grid->pnt[ip2].v) / 3.0;
            grid->pnt[ipnt].s    = (grid->pnt[ip0].s + grid->pnt[ip1].s + grid->pnt[ip2].s) / 3.0;
            grid->pnt[ipnt].p[0] = ip0;
            grid->pnt[ipnt].p[1] = ip1;
            grid->pnt[ipnt].p[2] = ip2;

            /* if the new Point is closer to any of the Triangle's current
               Points than allowed by the spacing (%s), reject the Point */
            if (distance(grid, ip0, ipnt) < alfa*grid->pnt[ip0].s ||
                distance(grid, ip1, ipnt) < alfa*grid->pnt[ip1].s ||
                distance(grid, ip2, ipnt) < alfa*grid->pnt[ip2].s   ) {
                grid->npnt--;
                continue;
            }

            /* if the new Point is close to any of the Points being added on
               this pass, reject the Point */
            for (ii = nsave; ii < grid->npnt-1; ii++) {
                if (distance(grid, ii, ipnt) < beta*grid->pnt[ipnt].s) {
                    grid->npnt--;
                    break; // out of do ii loop
                }
            }
        }

        /* now that we know all the Points that should be added in the
           current pass, insert them one by one */
        for (ipnt = nsave; ipnt < grid->npnt; ipnt++) {
            status = insertPoint(grid, ipnt);
            if (status < SUCCESS) goto cleanup;
        }

        /* next pass*/
#ifdef DEBUG2
        printf(".....pass %3d                   npnt=%5d\n", ipass, grid->npnt);
        printf("                                ntri=%5d\n",        grid->ntri);
#endif
    }
#ifdef GRAFIC
    {
        int  indgr    = 1+2+4+16+64;
        char pltitl[80];

        sprintf(pltitl, "~u~v~after pass %3d", ipass);

        grctrl_(plotGrid, &indgr, pltitl,
                (void*)(grid),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                strlen(pltitl));
    }
#endif

    /* loop through all Triangles and find another Triangle that
       shares it side */
    for (itri = 0; itri < grid->ntri; itri++) {
        if (grid->tri[itri].rr < 0) continue;

        for (isid = 0; isid < 3; isid++) {
            if (grid->tri[itri].t[isid] >= 0) continue;

            for (jtri = itri+1; jtri < grid->ntri; jtri++) {
                if (grid->tri[jtri].rr < 0) continue;

                for (jsid = 0; jsid < 3; jsid++) {
                    if (grid->tri[jtri].t[jsid] >= 0) continue;

                    if (grid->tri[itri].p[(isid+1)%3] == grid->tri[jtri].p[(jsid+2)%3] &&
                        grid->tri[itri].p[(isid+2)%3] == grid->tri[jtri].p[(jsid+1)%3]   ) {
                        grid->tri[itri].t[isid] = jtri;
                        grid->tri[jtri].t[jsid] = itri;
                    }
                }
            }
        }
    }

#ifdef DEBUG2
    for (itri = 0; itri < grid->ntri; itri++) {
        if (grid->tri[itri].rr < 0) continue;

        for (isid = 0; isid < 3; isid++) {
            if (grid->tri[itri].t[isid] < 0) {
                printf("unmatched side for itri=%5d, isid=%d between %5d and %5d\n",
                       itri, isid, grid->tri[itri].p[(isid+1)%3], grid->tri[itri].p[(isid+2)%3]);
            }
        }
    }
    printf("Neighbor table generated:       npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

    /* recover the Boundary by performing edge swaps if necessary
       for each segment of the Boundary (including the one that
       connects the first and last point of the Boundary)... */
    iend = -1;
    for (ilup = 0; ilup < nlup; ilup++) {
        ibeg = iend + 1;
        iend = ibeg +  lup[ilup] - 1;

        ip0 = iend;
        for (ip1 = ibeg; ip1 <= iend; ip1++) {

            /* look for a Triangle that has this Boundary segment
               on one of its sides */
            for (itri = 0; itri < grid->ntri; itri++) {
                if (grid->tri[itri].rr < 0) continue;

                if (grid->tri[itri].p[0] == ip0  && grid->tri[itri].p[1] == ip1) goto next_seg;
                if (grid->tri[itri].p[1] == ip0  && grid->tri[itri].p[2] == ip1) goto next_seg;
                if (grid->tri[itri].p[2] == ip0  && grid->tri[itri].p[0] == ip1) goto next_seg;
            }

            /* if no Triangles were found, then we need to swap edges.  in
               this case, look for the two Triangles that are configured as:

                              ip3
                             / | \
                            /  |  \
                           /   |   \
                          /    |    \
                         / it0 | it1 \
                        /      |      \
                      ip0 ' ' '|' ' ' ip1   <-----Boundary segment that we want
                        \      |      /
                         \     |     /
                          \    |    /
                           \   |   /
                            \  |  /
                             \ | /
                              ip2
            */
            for (it0 = 0; it0 < grid->ntri; it0++) {
                if (grid->tri[it0].rr < 0) continue;

                for (isid = 0; isid < 3; isid++) {
                    if (grid->tri[it0].p[isid] == ip0) {

                        for (it1 = 0; it1 < grid->ntri; it1++) {
                            if (grid->tri[it1].rr <  0  ) continue;
                            if (    it1     == it0) continue;

                            for (jsid = 0; jsid < 3; jsid++) {
                                if (grid->tri[it1].p[jsid] == ip1) {

                                    if (grid->tri[it0].p[(isid+1)%3] == grid->tri[it1].p[(jsid+2)%3] &&
                                        grid->tri[it0].p[(isid+2)%3] == grid->tri[it1].p[(jsid+1)%3]   ) {
                                        status = swapDiagonals(grid, it0, it1);
                                        if (status < SUCCESS) goto cleanup;

                                        /* undo swap if negative area resulted */
                                        if (computeArea(grid, it0) <= 0 || computeArea(grid, it1) <= 0.0) {
                                            status = swapDiagonals(grid, it0, it1);
                                            if (status < SUCCESS) goto cleanup;
                                        } else {
                                            goto next_seg;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /* the above techniques did not work, so try flipping diagonals
               until we have succeeded.  Note that this algorithm does not
               check for folds */
            nchange = 0;
            for (ipass = 0; ipass < grid->ntri; ipass++) {
                nchange = 0;

                /* look for the Triangle which is attached to ip0 and whose side
                   is pierced by a line between ip0 and ip1 */
                jtri = -1;
                for (itri = 0; itri < grid->ntri; itri++) {
                    if (grid->tri[itri].rr < 0) continue;

                    if (grid->tri[itri].p[0] == ip0) {
                        if (intersect(grid, ip0, ip1, grid->tri[itri].p[1], grid->tri[itri].p[2]) == 1) {
                            jtri = grid->tri[itri].t[0];
                            status = swapDiagonals(grid, itri, jtri);
                            nchange++;
                        }
                    } else if (grid->tri[itri].p[1] == ip0) {
                        if (intersect(grid, ip0, ip1, grid->tri[itri].p[2], grid->tri[itri].p[0]) == 1) {
                            jtri = grid->tri[itri].t[1];
                            status = swapDiagonals(grid, itri, jtri);
                            nchange++;
                        }
                    } else if (grid->tri[itri].p[2] == ip0) {
                        if (intersect(grid, ip0, ip1, grid->tri[itri].p[0], grid->tri[itri].p[1]) == 1) {
                            jtri = grid->tri[itri].t[2];
                            status = swapDiagonals(grid, itri, jtri);
                            nchange++;
                        }
                    } else if (grid->tri[itri].p[0] == ip1) {
                        if (intersect(grid, ip0, ip1, grid->tri[itri].p[1], grid->tri[itri].p[2]) == 1) {
                            jtri = grid->tri[itri].t[0];
                            status = swapDiagonals(grid, itri, jtri);
                            nchange++;
                        }
                    } else if (grid->tri[itri].p[1] == ip1) {
                        if (intersect(grid, ip0, ip1, grid->tri[itri].p[2], grid->tri[itri].p[0]) == 1) {
                            jtri = grid->tri[itri].t[1];
                            status = swapDiagonals(grid, itri, jtri);
                            nchange++;
                        }
                    } else if (grid->tri[itri].p[2] == ip1) {
                        if (intersect(grid, ip0, ip1, grid->tri[itri].p[0], grid->tri[itri].p[1]) == 1) {
                            jtri = grid->tri[itri].t[2];
                            status = swapDiagonals(grid, itri, jtri);
                            nchange++;
                        }
                    }
                }

                /* check to see if we have recovered the boundary segment yet */
                if (jtri >= 0) {
                    if        (grid->tri[jtri].p[0] == ip0 && grid->tri[jtri].p[1] == ip1) {
                        goto next_seg;
                    } else if (grid->tri[jtri].p[1] == ip0 && grid->tri[jtri].p[2] == ip1) {
                        goto next_seg;
                    } else if (grid->tri[jtri].p[2] == ip0 && grid->tri[jtri].p[0] == ip1) {
                        goto next_seg;
                    } else if (grid->tri[jtri].p[0] == ip1 && grid->tri[jtri].p[1] == ip0) {
                        goto next_seg;
                    } else if (grid->tri[jtri].p[1] == ip1 && grid->tri[jtri].p[2] == ip0) {
                        goto next_seg;
                    } else if (grid->tri[jtri].p[2] == ip1 && grid->tri[jtri].p[0] == ip0) {
                        goto next_seg;
                    }
                }
            }

            /* if we have gotten here, then we have found no set of
               Triangles for which an edge swap will recover the Boundary,
               so inform the user */

            if (nchange == 0) {
                printf("ERROR:: problem between %5d and %5d\n", ip0, ip1);
                status = COULD_NOT_RECOVER_BND;
                goto cleanup;
            }

            /* next Boundary segment */
        next_seg:
            ip0 = ip1;
        }
    }
#ifdef DEBUG2
    printf("Boundaries recovered:           npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

#ifdef GRAFIC
    {
        int  indgr    = 1+2+4+16+64;
        char pltitl[] = "~u~v~Boundaries recovered";

        grctrl_(plotGrid, &indgr, pltitl,
                (void*)(grid),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                strlen(pltitl));
    }
#endif

    /* break the neighbor information for all Triangle edges that
       correspond to Boundary segments */
    iend = -1;
    for (ilup = 0; ilup < nlup; ilup++) {
        ibeg = iend + 1;
        iend = ibeg + lup[ilup] - 1;

        ip0 = iend;
        for (ip1 = ibeg; ip1 <= iend; ip1++) {

            for (itri = 0; itri < grid->ntri; itri++) {
                if (grid->tri[itri].rr < 0) continue;

                /* these are adjacent to boundary on interior */
                if (       grid->tri[itri].p[0] == ip0  && grid->tri[itri].p[1] == ip1) {
                    grid->tri[itri].t[2] = 0;
                } else if (grid->tri[itri].p[1] == ip0  && grid->tri[itri].p[2] == ip1) {
                    grid->tri[itri].t[0] = 0;
                } else if (grid->tri[itri].p[2] == ip0  && grid->tri[itri].p[0] == ip1) {
                    grid->tri[itri].t[1] = 0;
                }

                /* these are adjacent to boundary on exterior */
                if (       grid->tri[itri].p[0] == ip1  && grid->tri[itri].p[1] == ip0) {
                    grid->tri[itri].t[2] = 0;
                    grid->tri[itri].rr   = -1;
                } else if (grid->tri[itri].p[1] == ip1  && grid->tri[itri].p[2] == ip0) {
                    grid->tri[itri].t[0] = 0;
                    grid->tri[itri].rr   = -1;
                } else if (grid->tri[itri].p[2] == ip1  && grid->tri[itri].p[0] == ip0) {
                    grid->tri[itri].t[1] = 0;
                    grid->tri[itri].rr   = -1;
                }
            }
            ip0 = ip1;
        }
    }

    /* any Triangle adjacent to an exterior Triangle is also exterior */
    for (ipass = 0; ipass < grid->ntri; ipass++) {
        nchange = 0;

        /* loop through all undeleted Triangles */
        for (itri = 0; itri < grid->ntri; itri++) {
            if (grid->tri[itri].rr < 0) continue;

            /* if a neighbor is deleted, then mark this Triangle for deletion also */
            for (isid = 0; isid < 3; isid++) {
                if (grid->tri[itri].t[isid] > 0) {
                    if (grid->tri[grid->tri[itri].t[isid]].rr < 0) {
                        grid->tri[itri].rr = -1;
                        nchange++;
                    }
                }
            }
        }
#ifdef DEBUG2
        printf(".....pass %5d   nchange=%5d\n", ipass, nchange);
#endif

        if (nchange == 0) break;
    }

#ifdef DEBUG2
    printf("External Triangles marked:      npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

    /* remove the unused Triangles from the data strucure.  start by
       determine the "new" location of each Triangle */
    itrinew = (int *) malloc(grid->ntri*sizeof(int));
    if (itrinew == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }
    ntrinew = 0;

    /* move non-deleted Triangles up i nthe list */
    for (itri = 0; itri < grid->ntri; itri++) {
        if (grid->tri[itri].rr >= 0) {
            itrinew[itri] = ntrinew;

            if (itri != ntrinew) {
                grid->tri[ntrinew].p[0] = grid->tri[itri].p[0];
                grid->tri[ntrinew].p[1] = grid->tri[itri].p[1];
                grid->tri[ntrinew].p[2] = grid->tri[itri].p[2];
                grid->tri[ntrinew].t[0] = grid->tri[itri].t[0];
                grid->tri[ntrinew].t[1] = grid->tri[itri].t[1];
                grid->tri[ntrinew].t[2] = grid->tri[itri].t[2];
                grid->tri[ntrinew].uc   = grid->tri[itri].uc;
                grid->tri[ntrinew].vc   = grid->tri[itri].vc;
                grid->tri[ntrinew].rr   = grid->tri[itri].rr;

                ntrinew++;
            }
        } else {
            itrinew[itri] = -1;
        }
    }

    /* adjust the neighbor info based upon the new locations */
    for (itri = 0; itri < ntrinew; itri++) {
        if (grid->tri[itri].rr < 0) continue;

        if (grid->tri[itri].t[0] >= 0) {
            grid->tri[itri].t[0] = itrinew[grid->tri[itri].t[0]];
        }
        if (grid->tri[itri].t[1] >= 0) {
            grid->tri[itri].t[1] = itrinew[grid->tri[itri].t[1]];
        }
        if (grid->tri[itri].t[2] >= 0) {
            grid->tri[itri].t[2] = itrinew[grid->tri[itri].t[2]];
        }
    }
    grid->ntri = ntrinew;

    free(itrinew);

#ifdef DEBUG2
    printf("External Triangles deleted:     npnt=%5d\n", grid->npnt);
    printf("                                ntri=%5d\n", grid->ntri);
#endif

#ifdef GRAFIC
    {
        int  indgr    = 1+2+4+16+64;
        char pltitl[] = "~u~v~External Triangles deleted";

        grctrl_(plotGrid, &indgr, pltitl,
                (void*)(grid),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                (void*)(NULL),
                strlen(pltitl));
    }
#endif

    /* count the valence of each Point */
    valence = (int *) malloc(grid->npnt*sizeof(int));
    if (valence == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }

    for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
        valence[ipnt] = 0;
    }

    for (itri = 0; itri < grid->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            ipnt = grid->tri[itri].p[isid];
            if (ipnt >= 0) {
                valence[ipnt]++;
            }
        }
    }

    /* determine the new location for each Point */
    nextpnt = 0;
    for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
        if (valence[ipnt] > 0) {
            valence[ipnt] = nextpnt;

            grid->pnt[nextpnt].u    = grid->pnt[ipnt].u;
            grid->pnt[nextpnt].v    = grid->pnt[ipnt].v;
            grid->pnt[nextpnt].s    = grid->pnt[ipnt].s;
            grid->pnt[nextpnt].p[0] = grid->pnt[ipnt].p[0];
            grid->pnt[nextpnt].p[1] = grid->pnt[ipnt].p[1];
            grid->pnt[nextpnt].p[2] = grid->pnt[ipnt].p[2];

            nextpnt++;
        } else {
            valence[ipnt] = -1;
        }
    }

    grid->npnt = nextpnt;

    /* adjust the pointers in the Triangles */
    for (itri = 0; itri < grid->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            ipnt = grid->tri[itri].p[isid];
            if (ipnt >= 0) {
                grid->tri[itri].p[isid] = valence[ipnt];
            }
        }
    }

    /* adjust the pointers in the Points */
    for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
        for (isid = 0; isid < 3; isid++) {
            jpnt = grid->pnt[ipnt].p[isid];
            if (jpnt >= 0) {
                grid->pnt[ipnt].p[isid] = valence[jpnt];
            }
        }
    }

    free(valence);

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * eggMorph - morph a Grid                                             *
 *                                                                     *
 ***********************************************************************
 */
int
eggMorph(void  *gridP,                  /* (in)  pointer to grid */
         double uvnew[],                /* (in)  array of new coordinates for Boundary Points */
                                        /* (out) array of new coordinates */
         void   **newGridP)             /* (in)  pointer to new grid (which is a grid_T*) */
{
    int    status = SUCCESS;            /* (out) Point index (bias-0) */

    int    ipnt, ip0, ip1, ip2, itri;

    grid_T *grid = (grid_T *) gridP;

    grid_T *newGrid;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("eggMorph(gridP=%llx)\n", (long long)gridP);
#endif

    /* create the grid structure and copy in everything except
       Point locations (which are updated below) */
    newGrid = (grid_T *) malloc(sizeof(grid_T));
    if (newGrid == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }

    *newGridP = (void *) newGrid;

    newGrid->mpnt = grid->npnt;
    newGrid->npnt = grid->npnt;
    newGrid->nbnd = grid->nbnd;
    newGrid->pnt  = (pnt_T *) malloc(newGrid->mpnt*sizeof(pnt_T));
    if (newGrid->pnt == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }
    newGrid->uv = NULL;
    newGrid->p  = NULL;

    for (ipnt = 0; ipnt < newGrid->npnt; ipnt++) {
        newGrid->pnt[ipnt].s    = grid->pnt[ipnt].s;
        newGrid->pnt[ipnt].p[0] = grid->pnt[ipnt].p[0];
        newGrid->pnt[ipnt].p[1] = grid->pnt[ipnt].p[1];
        newGrid->pnt[ipnt].p[2] = grid->pnt[ipnt].p[2];
    }

    newGrid->mtri = grid->ntri;
    newGrid->ntri = grid->ntri;
    newGrid->tri  = (tri_T *) malloc(grid->mtri*sizeof(tri_T));
    if (newGrid->tri == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }
    newGrid->tris = NULL;

    for (itri = 0; itri < newGrid->ntri; itri++) {
        newGrid->tri[itri].p[0] = grid->tri[itri].p[0];
        newGrid->tri[itri].p[1] = grid->tri[itri].p[1];
        newGrid->tri[itri].p[2] = grid->tri[itri].p[2];
        newGrid->tri[itri].t[0] = grid->tri[itri].t[0];
        newGrid->tri[itri].t[1] = grid->tri[itri].t[1];
        newGrid->tri[itri].t[2] = grid->tri[itri].t[2];
        newGrid->tri[itri].uc   = grid->tri[itri].uc;
        newGrid->tri[itri].vc   = grid->tri[itri].vc;
        newGrid->tri[itri].rr   = grid->tri[itri].rr;
    }

    /* copy in the Boundary Points */
    for (ipnt = 0; ipnt < newGrid->nbnd; ipnt++) {
        newGrid->pnt[ipnt].u = uvnew[2*ipnt  ];
        newGrid->pnt[ipnt].v = uvnew[2*ipnt+1];
    }

    /* update the interior Points */
#ifndef __clang_analyzer__
    for (ipnt = newGrid->nbnd; ipnt < newGrid->npnt; ipnt++) {
        ip0 = newGrid->pnt[ipnt].p[0];
        ip1 = newGrid->pnt[ipnt].p[1];
        ip2 = newGrid->pnt[ipnt].p[2];

        newGrid->pnt[ipnt].u = (newGrid->pnt[ip0].u + newGrid->pnt[ip1].u + newGrid->pnt[ip2].u) / 3;
        newGrid->pnt[ipnt].v = (newGrid->pnt[ip0].v + newGrid->pnt[ip1].v + newGrid->pnt[ip2].v) / 3;

        uvnew[2*ipnt  ] = newGrid->pnt[ipnt].u;
        uvnew[2*ipnt+1] = newGrid->pnt[ipnt].v;
    }
#endif

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * eggInfo - get info about a Grid                                     *
 *                                                                     *
 ***********************************************************************
 */
int
eggInfo(void   *gridP,                  /* (in)  pointer to grid */
        int    *npnt,                   /* (out) number of Points */
        int    *nbnd,                   /* (out) number of boundary Points */
  const double *uv[],                   /* (out) parametric coordinates */
  const int    *p[],                    /* (out) parent Points */
        int    *ntri,                   /* (out) number of Triangles */
  const int    *tris[])                 /* (out) indices of Points in Triangles (bias-0) */
{
    int    status = SUCCESS;

    int    ipnt, itri;

    grid_T *grid = (grid_T *) gridP;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("eggInfo(gridP=%llx)\n", (long long)gridP);
#endif

    if (grid->uv == NULL) {
        grid->uv   = (double *) malloc(2*grid->npnt*sizeof(double));
        grid->p    = (int    *) malloc(3*grid->npnt*sizeof(int   ));
        grid->tris = (int    *) malloc(3*grid->ntri*sizeof(int   ));

        if (grid->uv == NULL || grid->p == NULL || grid->tris == NULL) {
            status = MALLOC_ERROR;
            goto cleanup;
        }

        for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
            grid->uv[2*ipnt  ] = grid->pnt[ipnt].u;
            grid->uv[2*ipnt+1] = grid->pnt[ipnt].v;
            grid->p [3*ipnt  ] = grid->pnt[ipnt].p[0];
            grid->p [3*ipnt+1] = grid->pnt[ipnt].p[1];
            grid->p [3*ipnt+2] = grid->pnt[ipnt].p[2];
        }

        for (itri = 0; itri < grid->ntri; itri++) {
            grid->tris[3*itri  ] = grid->tri[itri].p[0];
            grid->tris[3*itri+1] = grid->tri[itri].p[1];
            grid->tris[3*itri+2] = grid->tri[itri].p[2];
        }
    }

    *npnt = grid->npnt;
    *nbnd = grid->nbnd;
    *uv   = grid->uv;
    *p    = grid->p;
    *ntri = grid->ntri;
    *tris = grid->tris;

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * eggDump - dump a Grid to a file                                     *
 *                                                                     *
 ***********************************************************************
 */
int
eggDump(void   *gridP,                  /* (in)  pointer to grid */
        FILE   *fp)                     /* (in)  pointer to ASCII file */
{
    int    status = SUCCESS;

    int    ipnt, itri;

    grid_T *grid = (grid_T *) gridP;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("eggDump(gridP=%llx, fp=%llx)\n", (long long)gridP, (long long)fp);
#endif

    /* if grid is NULL, then write an empty header and return */
    if (grid == NULL) {
        fprintf(fp, "%7d %7d %7d\n", 0, 0, 0);
        goto cleanup;
    }

    /* write the header */
    fprintf(fp, "%7d %7d %7d\n",
            grid->npnt, grid->nbnd, grid->ntri);

    /* write the Point table */
    for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
        fprintf(fp, "%21.15e %21.15e %21.15e %7d %7d %7d\n",
                grid->pnt[ipnt].u,
                grid->pnt[ipnt].v,
                grid->pnt[ipnt].s,
                grid->pnt[ipnt].p[0],
                grid->pnt[ipnt].p[1],
                grid->pnt[ipnt].p[2]);
    }

    /* write the Triangle table */
    for (itri = 0; itri < grid->ntri; itri++) {
        fprintf(fp, "%7d %7d %7d %7d %7d %7d %21.15e %21.15e %21.15e\n",
                grid->tri[itri].p[0],
                grid->tri[itri].p[1],
                grid->tri[itri].p[2],
                grid->tri[itri].t[0],
                grid->tri[itri].t[1],
                grid->tri[itri].t[2],
                grid->tri[itri].uc,
                grid->tri[itri].vc,
                grid->tri[itri].rr);
    }

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * eggLoad - load Grid from a file                                     *
 *                                                                     *
 ***********************************************************************
 */
int
eggLoad(FILE   *fp,                     /* (in)  pointer to ASCII file */
        void   **gridP)                 /* (out) pointer to *grid */
{
    int    status = SUCCESS;

    int    npnt, nbnd, ntri, ipnt, itri;

    grid_T *grid;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("eggLoad(fp=%llx)\n", (long long)fp);
#endif

    /* read the header */
    fscanf(fp, "%d %d %d\n", &npnt, &nbnd, &ntri);

    /* if there are no Points or Triangles, return a NULL */
    if (npnt <= 0 || ntri <= 0) {
        *gridP = NULL;
        goto cleanup;
    }

    /* create the Grid structure */
    grid = (grid_T *) malloc(sizeof(grid_T));
    if (grid == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }

    /* return a pointer to grid */
    *gridP = (void *) grid;

    grid->npnt = npnt;
    grid->nbnd = nbnd;
    grid->ntri = ntri;

    /* create initial Point and Tringle tables */
    grid->mpnt = grid->npnt;
    grid->pnt  = (pnt_T *) malloc(grid->mpnt*sizeof(pnt_T));
    if (grid->pnt == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }

    grid->mtri = grid->ntri;
    grid->tri  = (tri_T *) malloc(grid->mtri*sizeof(tri_T));
    if (grid->tri == NULL) {
        status = MALLOC_ERROR;
        goto cleanup;
    }

    /* read the Point table */
    for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
        fscanf(fp, "%lf %lf %lf %d %d %d\n",
               &(grid->pnt[ipnt].u),
               &(grid->pnt[ipnt].v),
               &(grid->pnt[ipnt].s),
               &(grid->pnt[ipnt].p[0]),
               &(grid->pnt[ipnt].p[1]),
               &(grid->pnt[ipnt].p[2]));
    }

    /* read the Triangle table */
    for (itri = 0; itri < grid->ntri; itri++) {
        fscanf(fp, "%d %d %d %d %d %d %lf %lf %lf\n",
               &(grid->tri[itri].p[0]),
               &(grid->tri[itri].p[1]),
               &(grid->tri[itri].p[2]),
               &(grid->tri[itri].t[0]),
               &(grid->tri[itri].t[1]),
               &(grid->tri[itri].t[2]),
               &(grid->tri[itri].uc),
               &(grid->tri[itri].vc),
               &(grid->tri[itri].rr));
    }

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * eggFree - free a Grid structure                                     *
 *                                                                     *
 ***********************************************************************
 */
int
eggFree(void   *gridP)                  /* (in)  pointer to grid */
{
    int    status = SUCCESS;

    grid_T *grid = (grid_T *) gridP;

    /* -------------------------------------------------------------- */

#ifdef DEBUG
    printf("eggFree(gridP=%llx)\n", (long long)gridP);
#endif

    if (grid != NULL) {
        if (grid->pnt  != NULL) free(grid->pnt);
        if (grid->tri  != NULL) free(grid->tri);

        if (grid->uv   != NULL) free(grid->uv);
        if (grid->p    != NULL) free(grid->p );
        if (grid->tris != NULL) free(grid->tris);

        free(grid);
    }

    grid = NULL;

//cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * addToLoop - add an edge defined by (ip0,ip1) to the loop table,     *
 *             taking care to delete duplicate (but reversed) edges    *
 *                                                                     *
 ***********************************************************************
 */
static int
addToLoop(int    ip0,                   /* (in)  Point at beginning of new Edge (bias-0) */
          int    ip1,                   /* (in)  Point at end       of new Edge (bias-0) */
          int    iloop[],               /* (both)table of loop of Point pairs   */
          int    *nloop)                /* (both)number of Point pairs in table */
{
    int    status = SUCCESS;            /* (out) default return */

    int    i, j;

    /* -------------------------------------------------------------- */

    /* if edge is already in loop (but in reverse order), delete it from
       the loop list (and return immediately) */
    for (i = 0; i < *nloop; i++) {
        if (iloop[2*i] == ip1 && iloop[2*i+1] == ip0) {
            (*nloop)--;

            for (j = i; j < *nloop; j++) {
                iloop[2*j  ] = iloop[2*j+2];
                iloop[2*j+1] = iloop[2*j+3];
            }
            goto cleanup;
        }
    }

    /* if we got here, the edge was not in the list, so add it to the
       end of the list (and return) */
    iloop[2*(*nloop)  ] = ip0;
    iloop[2*(*nloop)+1] = ip1;
    (*nloop)++;

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * computeArea --- compute Triangle area                               *
 *                                                                     *
 ***********************************************************************
 */
static double
computeArea(grid_T *grid,               /* (in)  pointer to grid */
            int    itri)                /* (in)  Triangle index (bias-0) */
{
    double result = 0;                  /* Triangle area */

    int    ip0, ip1, ip2;

    /* -------------------------------------------------------------- */

    /* get a local copy of the coordinates of the Points */
    ip0 = grid->tri[itri].p[0];
    ip1 = grid->tri[itri].p[1];
    ip2 = grid->tri[itri].p[2];

    /* compute area */
    result = ( (grid->pnt[ip1].u-grid->pnt[ip0].u) * (grid->pnt[ip2].v-grid->pnt[ip0].v)
             - (grid->pnt[ip1].v-grid->pnt[ip0].v) * (grid->pnt[ip2].u-grid->pnt[ip0].u)) / 2;

//cleanup:
    return result;
}


/*
 ***********************************************************************
 *                                                                     *
 * computeParameters - compute Triangle parameters                     *
 *                                                                     *
 ***********************************************************************
 */
static int
computeParameters(grid_T *grid,         /* (in)  pointer to grid */
                  int    itri)          /* (in)  Triangle index (bias-0) */
{
    int    status = SUCCESS;            /* (out) default return */

    int    ip0, ip1, ip2;
    double u0, v0, u1, v1, u2, v2, s;

    /* -------------------------------------------------------------- */

    if (itri < 0 || itri >= grid->ntri) {
        status = BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    /* get coordinates */
    ip0 = grid->tri[itri].p[0];
    u0  = grid->pnt[ip0].u;
    v0  = grid->pnt[ip0].v;

    ip1 = grid->tri[itri].p[1];
    u1  = grid->pnt[ip1].u;
    v1  = grid->pnt[ip1].v;

    ip2 = grid->tri[itri].p[2];
    u2  = grid->pnt[ip2].u;
    v2  = grid->pnt[ip2].v;

    /* find the intersections of the perpendicular bisectors of
       edges (0-1) and (1-2) to get the circumcenter */
    s = ((u2-u0) * (u1-u2) - (v2-v0) * (v2-v1))
      / ((v0-v1) * (u1-u2) - (v2-v1) * (u1-u0));

    grid->tri[itri].uc = (u0 + u1 + s * (v0 - v1)) / 2.0;
    grid->tri[itri].vc = (v0 + v1 + s * (u1 - u0)) / 2.0;

    /* compute the radius^2 of the circumcircle */
    grid->tri[itri].rr = SQR(grid->tri[itri].uc-u0) + SQR(grid->tri[itri].vc-v0);

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * createTriangle - create a Triangle that joins p01, ip1, and ip2     *
 *                                                                     *
 ***********************************************************************
 */
static int
createTriangle(grid_T *grid,            /* (in)  pointer to grid */
               int    ip0,              /* (in)  first  Point index (bias-0) */
               int    ip1,              /* (in)  second Point index (bias-0) */
               int    ip2)              /* (in)  third  Point index (bias-0) */
{
    int    status = SUCCESS;            /* (out) default return */

    int    i, inew;

    /* -------------------------------------------------------------- */

    /* if there is a Triangle marked for deletion, overwrite it with
       the new Triangle */
    inew = -1;
    for (i = 0; i < grid->ntri; i++) {
        if (inew<0 && grid->tri[i].rr<0) {
            inew = i;
        }
    }

    /* if no deleted Triangles were found, add the new Triangle to
       the end of the Triangle array */
    if (inew < 0) {
        inew = newTriangle(grid);
    }

    /* store the Point numbers in the new Triangle */
    grid->tri[inew].p[0] = ip0;
    grid->tri[inew].p[1] = ip1;
    grid->tri[inew].p[2] = ip2;
    grid->tri[inew].t[0] = -1;
    grid->tri[inew].t[1] = -1;
    grid->tri[inew].t[2] = -1;

    /* compute Triangle parameters */
    status = computeParameters(grid, inew);

//cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * distance - compute distance between two Points                      *
 *                                                                     *
 ***********************************************************************
 */
static double
distance(grid_T *grid,                  /* (in)  pointer to grid */
         int    ip0,                    /* (in)  first  Point index (bias-0) */
         int    ip1)                    /* (in)  second Point index (bias-0) */
{
    double result = 0;                  /* (out)  distance */

    /* -------------------------------------------------------------- */

#ifndef __clang_analyzer__
    result = sqrt( SQR(grid->pnt[ip0].u - grid->pnt[ip1].u)
                 + SQR(grid->pnt[ip0].v - grid->pnt[ip1].v));
#endif

//cleanup:
    return result;
}


/*
 ***********************************************************************
 *                                                                     *
 * insertPoint - insert Point into current traingulation               *
 *                                                                     *
 ***********************************************************************
 */
static int
insertPoint(grid_T *grid,               /* (in)  pointer to grid */
           int    ipnt)                 /* (in)  Point index (bias-0) */
{
    int    status = SUCCESS;            /* (out) default return */

    /*
      local data structure to store a list of pairs of Points that
      correspond to the boundary of the hole formed when the
      Triangles whose circumcircles contain ipnt are deleted
    */

#define MLOOP 10000
    int    nloop, iloop[2*MLOOP];

    int    itri, jloop;
    double rtest;

    /* -------------------------------------------------------------- */

    /*  clear the loop table */
    nloop = 0;

    /* find all Triangles whose circumcircles contain the new point */
    for (itri = 0; itri < grid->ntri; itri++) {
        if (grid->tri[itri].rr < 0) continue;

        rtest = SQR(grid->pnt[ipnt].u-grid->tri[itri].uc) + SQR(grid->pnt[ipnt].v-grid->tri[itri].vc);

        if (rtest < grid->tri[itri].rr) {

            /* put the appropriate edges in the loop table (note that
               addToLoop will take care of deleting duplicate edges) */
            status = addToLoop(grid->tri[itri].p[0],
                               grid->tri[itri].p[1], iloop, &nloop);
            if (status < SUCCESS) goto cleanup;

            status = addToLoop(grid->tri[itri].p[1],
                               grid->tri[itri].p[2], iloop, &nloop);
            if (status < SUCCESS) goto cleanup;

            status = addToLoop(grid->tri[itri].p[2],
                               grid->tri[itri].p[0], iloop, &nloop);
            if (status < SUCCESS) goto cleanup;

            grid->tri[itri].rr = -1.0;        // mark Triangle as deleted
         }
      }

    /* make new Triangles using the new Point and the Point pairs in the loop */
    for (jloop = 0; jloop < nloop; jloop++) {
        status = createTriangle(grid, iloop[2*jloop], iloop[2*jloop+1], ipnt);
        if (status < SUCCESS) goto cleanup;
    }

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * intersect - check if two line segments intersect                    *
 *                                                                     *
 ***********************************************************************
 */
static int
intersect(grid_T *grid,                 /* (in)  pointer to grid */
          int    ip0,                   /* (in)  first  Point on first  line */
          int    ip1,                   /* (in)  second Point on first  line */
          int    ip2,                   /* (in)  first  Point on second line */
          int    ip3)                   /* (in)  second Point on second line */
{
    int    status = SUCCESS;            /* (out) =0 no intersection */
                                        /*       =1    intersection */

    double u0, v0, u1, v1, u2, v2, u3, v3, d, s, t;

    /* -------------------------------------------------------------- */

    if (ip0 >= 0 && ip0 < grid->ntri) {
        u0 = grid->pnt[ip0].u;
        v0 = grid->pnt[ip0].v;
    } else {
        status = BAD_POINT_INDEX;
        goto cleanup;
    }

    if (ip1 >= 0 && ip1 < grid->ntri) {
        u1 = grid->pnt[ip1].u;
        v1 = grid->pnt[ip1].v;
    } else {
        status = BAD_POINT_INDEX;
        goto cleanup;
    }

    if (ip2 >= 0 && ip2 < grid->ntri) {
        u2 = grid->pnt[ip2].u;
        v2 = grid->pnt[ip2].v;
    } else {
        status = BAD_POINT_INDEX;
        goto cleanup;
    }

    if (ip3 >= 0 && ip3 < grid->ntri) {
        u3 = grid->pnt[ip3].u;
        v3 = grid->pnt[ip3].v;
    } else {
        status = BAD_POINT_INDEX;
        goto cleanup;
    }

    /* compute diagonal.  if nearly zero, then lines are parallel */
    d = ((u1-u0)*(v2-v3) - (v1-v0)*(u2-u3));
    if (fabs(d) < EPS06) {
        status = 0;
        goto cleanup;
    }

    /* compute fractional distance from ip0 to ip1 */
    s = ((u2-u0)*(v2-v3) - (v2-v0)*(u2-u3)) / d;
    if (s < 0  ||  s > 1) {
        status = 0;
        goto cleanup;
    }

    /* compute fractional distance from ip2 to ip3 */
    t = ((u1-u0)*(v2-v0) - (v1-v0)*(u2-u0)) / d;
    if (t < 0.0  ||  t > 1.0) {
        status = 0;
        goto cleanup;
    }

    /* we have passed all checks, so the lines intersect */
    status = 1;

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * newPoint - get a new Point                                          *
 *                                                                     *
 ***********************************************************************
 */
static int
newPoint(grid_T *grid)                   /* (in)  pointer to grid */
{
    int    status = SUCCESS;            /* (out) Point index (bias-0) */

    /* -------------------------------------------------------------- */

    /* increase size of "pnt" array if required */
    if (grid->npnt >= grid->mpnt-1) {
        grid->mpnt += 1000;
        grid->pnt = (pnt_T*) realloc(grid->pnt, grid->mpnt*sizeof(pnt_T));
        if (grid->pnt == NULL) {
            status = MALLOC_ERROR;
            goto cleanup;
        }
    }

    /* add new Point */
    status = grid->npnt;
    grid->npnt++;

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * newTriangle - get a new Triangle                                    *
 *                                                                     *
 ***********************************************************************
 */
static int
newTriangle(grid_T *grid)               /* (in)  pointer to grid */
{
    int    status = SUCCESS;            /* (out) Triangle index (bias-0) */

    /* -------------------------------------------------------------- */

    /* increase size of "tri" array if required */
    if (grid->ntri >= grid->mtri-1) {
        grid->mtri += 1000;
        grid->tri = (tri_T*) realloc(grid->tri, grid->mtri*sizeof(tri_T));
        if (grid->tri == NULL) {
            status = MALLOC_ERROR;
            goto cleanup;
        }
    }

    /* add new Triangle */
    grid->tri[grid->ntri].p[0] = -1;
    grid->tri[grid->ntri].p[1] = -1;
    grid->tri[grid->ntri].p[2] = -1;
    grid->tri[grid->ntri].t[0] = -1;
    grid->tri[grid->ntri].t[1] = -1;
    grid->tri[grid->ntri].t[2] = -1;
    grid->ntri++;

    status = grid->ntri-1;

cleanup:
    return status;
}


/*
 ***********************************************************************
 *                                                                     *
 * swapDiagonals - swap diagonals of two adjacent Triangles            *
 *                                                                     *
 ***********************************************************************
 */
static int
swapDiagonals(grid_T *grid,             /* (in)  pointer to grid */
              int    ia,                /* (in)  first  Triangle index (bias-0) */
              int    ib)                /* (in)  second Triangle index (bias-0) */
{
    int    status = SUCCESS;            /* (out) default return */

    int    ia0, ia1, ia2, ja0, ja1, ja2, ib0, ib1, ib2, jb0, jb1, jb2;

    /* -------------------------------------------------------------- */

    if (ia < 0 || ia >= grid->ntri) {
        status = BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    ia0 = grid->tri[ia].p[0];
    ia1 = grid->tri[ia].p[1];
    ia2 = grid->tri[ia].p[2];
    ib0 = grid->tri[ib].p[0];
    ib1 = grid->tri[ib].p[1];
    ib2 = grid->tri[ib].p[2];

    if (ib < 0 || ib >= grid->ntri) {
        status = BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    ja0 = grid->tri[ia].t[0];
    ja1 = grid->tri[ia].t[1];
    ja2 = grid->tri[ia].t[2];
    jb0 = grid->tri[ib].t[0];
    jb1 = grid->tri[ib].t[1];
    jb2 = grid->tri[ib].t[2];

    /*
        ja1  ia2    ib1  jb2             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia0    |    |    ib0   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja2  ia1    ib2  jb1             n[0]
    */

    if (    ia1 == ib2  &&  ia2 == ib1) {
        grid->tri[ia].p[0] = ia1;
        grid->tri[ia].p[1] = ib0;
        grid->tri[ia].p[2] = ia0;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja2;
        grid->tri[ia].t[2] = jb1;

        grid->tri[ib].p[0] = ia2;
        grid->tri[ib].p[1] = ia0;
        grid->tri[ib].p[2] = ib0;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb2;
        grid->tri[ib].t[2] = ja1;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb1 >= 0) {
            if (grid->tri[jb1].t[0] == ib) grid->tri[jb1].t[0] = ia;
            if (grid->tri[jb1].t[1] == ib) grid->tri[jb1].t[1] = ia;
            if (grid->tri[jb1].t[2] == ib) grid->tri[jb1].t[2] = ia;
        }
        if (ja1 >= 0) {
            if (grid->tri[ja1].t[0] == ia) grid->tri[ja1].t[0] = ib;
            if (grid->tri[ja1].t[1] == ia) grid->tri[ja1].t[1] = ib;
            if (grid->tri[ja1].t[2] == ia) grid->tri[ja1].t[2] = ib;
        }

    /*
        ja1  ia2    ib0  jb1             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia0    |    |    ib2   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja2  ia1    ib1  jb0             n[0]
    */
    } else if (ia1 == ib1  &&  ia2 == ib0) {
        grid->tri[ia].p[0] = ia1;
        grid->tri[ia].p[1] = ib2;
        grid->tri[ia].p[2] = ia0;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja2;
        grid->tri[ia].t[2] = jb0;

        grid->tri[ib].p[0] = ia2;
        grid->tri[ib].p[1] = ia0;
        grid->tri[ib].p[2] = ib2;
        grid->tri[ib].t[0] = ib;
        grid->tri[ib].t[1] = jb1;
        grid->tri[ib].t[2] = ja1;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb0 >= 0) {
            if (grid->tri[jb0].t[0] == ib) grid->tri[jb0].t[0] = ia;
            if (grid->tri[jb0].t[1] == ib) grid->tri[jb0].t[1] = ia;
            if (grid->tri[jb0].t[2] == ib) grid->tri[jb0].t[2] = ia;
        }
        if (ja1 >= 0) {
            if (grid->tri[ja1].t[0] == ia) grid->tri[ja1].t[0] = ib;
            if (grid->tri[ja1].t[1] == ia) grid->tri[ja1].t[1] = ib;
            if (grid->tri[ja1].t[2] == ia) grid->tri[ja1].t[2] = ib;
        }

    /*
        ja1  ia2    ib2  jb0             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia0    |    |    ib1   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja2  ia1    ib0  jb2             n[0]
    */
    } else if (ia1 == ib0  &&  ia2 == ib2) {
        grid->tri[ia].p[0] = ia1;
        grid->tri[ia].p[1] = ib1;
        grid->tri[ia].p[2] = ia0;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja2;
        grid->tri[ia].t[2] = jb2;

        grid->tri[ib].p[0] = ia2;
        grid->tri[ib].p[1] = ia0;
        grid->tri[ib].p[2] = ib1;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb0;
        grid->tri[ib].t[2] = ja1;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb2 >= 0) {
            if (grid->tri[jb2].t[0] == ib) grid->tri[jb2].t[0] = ia;
            if (grid->tri[jb2].t[1] == ib) grid->tri[jb2].t[1] = ia;
            if (grid->tri[jb2].t[2] == ib) grid->tri[jb2].t[2] = ia;
        }
        if (ja1 >= 0) {
            if (grid->tri[ja1].t[0] == ia) grid->tri[ja1].t[0] = ib;
            if (grid->tri[ja1].t[1] == ia) grid->tri[ja1].t[1] = ib;
            if (grid->tri[ja1].t[2] == ia) grid->tri[ja1].t[2] = ib;
        }

    /*
        ja0  ia1    ib1  jb2             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia2    |    |    ib0   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja1  ia0    ib2  jb1             n[0]
    */
    } else if (ia0 == ib2  &&  ia1 == ib1) {
        grid->tri[ia].p[0] = ia0;
        grid->tri[ia].p[1] = ib0;
        grid->tri[ia].p[2] = ia2;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja1;
        grid->tri[ia].t[2] = jb1;

        grid->tri[ib].p[0] = ia1;
        grid->tri[ib].p[1] = ia2;
        grid->tri[ib].p[2] = ib0;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb2;
        grid->tri[ib].t[2] = ja0;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb1 >= 0) {
            if (grid->tri[jb1].t[0] == ib) grid->tri[jb1].t[0] = ia;
            if (grid->tri[jb1].t[1] == ib) grid->tri[jb1].t[1] = ia;
            if (grid->tri[jb1].t[2] == ib) grid->tri[jb1].t[2] = ia;
        }
        if (ja0 >= 0) {
            if (grid->tri[ja0].t[0] == ia) grid->tri[ja0].t[0] = ib;
            if (grid->tri[ja0].t[1] == ia) grid->tri[ja0].t[1] = ib;
            if (grid->tri[ja0].t[2] == ia) grid->tri[ja0].t[2] = ib;
        }

    /*
        ja0  ia1    ib0  jb1             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia2    |    |    ib2   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja1  ia0    ib1  jb0             n[0]
    */
    } else if (ia0 == ib1  &&  ia1 == ib0) {
        grid->tri[ia].p[0] = ia0;
        grid->tri[ia].p[1] = ib2;
        grid->tri[ia].p[2] = ia2;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja1;
        grid->tri[ia].t[2] = jb0;

        grid->tri[ib].p[0] = ia1;
        grid->tri[ib].p[1] = ia2;
        grid->tri[ib].p[2] = ib2;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb1;
        grid->tri[ib].t[2] = ja0;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb0 >= 0) {
            if (grid->tri[jb0].t[0] == ib) grid->tri[jb0].t[0] = ia;
            if (grid->tri[jb0].t[1] == ib) grid->tri[jb0].t[1] = ia;
            if (grid->tri[jb0].t[2] == ib) grid->tri[jb0].t[2] = ia;
        }
        if (ja0 >= 0) {
            if (grid->tri[ja0].t[0] == ia) grid->tri[ja0].t[0] = ib;
            if (grid->tri[ja0].t[1] == ia) grid->tri[ja0].t[1] = ib;
            if (grid->tri[ja0].t[2] == ia) grid->tri[ja0].t[2] = ib;
        }

    /*
        ja0  ia1    ib2  jb0             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia2    |    |    ib1   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja1  ia0    ib0  jb2             n[0]
    */
    } else if (ia0 == ib0  &&  ia1 == ib2) {
        grid->tri[ia].p[0] = ia0;
        grid->tri[ia].p[1] = ib1;
        grid->tri[ia].p[2] = ia2;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja1;
        grid->tri[ia].t[2] = jb2;

        grid->tri[ib].p[0] = ia1;
        grid->tri[ib].p[1] = ia2;
        grid->tri[ib].p[2] = ib1;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb0;
        grid->tri[ib].t[2] = ja0;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb2 >= 0) {
            if (grid->tri[jb2].t[0] == ib) grid->tri[jb2].t[0] = ia;
            if (grid->tri[jb2].t[1] == ib) grid->tri[jb2].t[1] = ia;
            if (grid->tri[jb2].t[2] == ib) grid->tri[jb2].t[2] = ia;
        }
        if (ja0 >= 0) {
            if (grid->tri[ja0].t[0] == ia) grid->tri[ja0].t[0] = ib;
            if (grid->tri[ja0].t[1] == ia) grid->tri[ja0].t[1] = ib;
            if (grid->tri[ja0].t[2] == ia) grid->tri[ja0].t[2] = ib;
        }

    /*
        ja2  ia0    ib1  jb2             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia1    |    |    ib0   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja0  ia2    ib2  jb1             n[0]
    */
    } else if (ia2 == ib2  &&  ia0 == ib1) {
        grid->tri[ia].p[0] = ia2;
        grid->tri[ia].p[1] = ib0;
        grid->tri[ia].p[2] = ia1;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja0;
        grid->tri[ia].t[2] = jb1;

        grid->tri[ib].p[0] = ia0;
        grid->tri[ib].p[1] = ia1;
        grid->tri[ib].p[2] = ib0;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb2;
        grid->tri[ib].t[2] = ja2;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb1 >= 0) {
            if (grid->tri[jb1].t[0] == ib) grid->tri[jb1].t[0] = ia;
            if (grid->tri[jb1].t[1] == ib) grid->tri[jb1].t[1] = ia;
            if (grid->tri[jb1].t[2] == ib) grid->tri[jb1].t[2] = ia;
        }
        if (ja2 >= 0) {
            if (grid->tri[ja2].t[0] == ia) grid->tri[ja2].t[0] = ib;
            if (grid->tri[ja2].t[1] == ia) grid->tri[ja2].t[1] = ib;
            if (grid->tri[ja2].t[2] == ia) grid->tri[ja2].t[2] = ib;
        }

    /*
        ja2  ia0    ib0  jb1             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia1    |    |    ib2   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja0  ia2    ib1  jb0             n[0]
    */
    } else if (ia2 == ib1  &&  ia0 == ib0) {
        grid->tri[ia].p[0] = ia2;
        grid->tri[ia].p[1] = ib2;
        grid->tri[ia].p[2] = ia1;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja0;
        grid->tri[ia].t[2] = jb0;

        grid->tri[ib].p[0] = ia0;
        grid->tri[ib].p[1] = ia1;
        grid->tri[ib].p[2] = ib2;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb1;
        grid->tri[ib].t[2] = ja2;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb0 >= 0) {
            if (grid->tri[jb0].t[0] == ib) grid->tri[jb0].t[0] = ia;
            if (grid->tri[jb0].t[1] == ib) grid->tri[jb0].t[1] = ia;
            if (grid->tri[jb0].t[2] == ib) grid->tri[jb0].t[2] = ia;
        }
        if (ja2 >= 0) {
            if (grid->tri[ja2].t[0] == ia) grid->tri[ja2].t[0] = ib;
            if (grid->tri[ja2].t[1] == ia) grid->tri[ja2].t[1] = ib;
            if (grid->tri[ja2].t[2] == ia) grid->tri[ja2].t[2] = ib;
        }

    /*
        ja2  ia0    ib2  jb0             n[0]
             / |    | \                  / \
            /  |    |  \                /ib\
           /   |    |   \              /     \
          /    |    |    \          n[1]-----n[2]
        ia1    |    |    ib1   ==>
          \    |    |    /          n[2]-----n[1]
           \   |    |   /              \     /
            \  |    |  /                \ia/
             \ |    | /                  \ /
        ja0  ia2    ib0  jb2             n[0]
    */
    } else if (ia2 == ib0  &&  ia0 == ib2) {
        grid->tri[ia].p[0] = ia2;
        grid->tri[ia].p[1] = ib1;
        grid->tri[ia].p[2] = ia1;
        grid->tri[ia].t[0] = ib;
        grid->tri[ia].t[1] = ja0;
        grid->tri[ia].t[2] = jb2;

        grid->tri[ib].p[0] = ia0;
        grid->tri[ib].p[1] = ia1;
        grid->tri[ib].p[2] = ib1;
        grid->tri[ib].t[0] = ia;
        grid->tri[ib].t[1] = jb0;
        grid->tri[ib].t[2] = ja2;

        status = computeParameters(grid, ia);
        if (status != SUCCESS) goto cleanup;
        status = computeParameters(grid, ib);
        if (status != SUCCESS) goto cleanup;

        if (jb2 >= 0) {
            if (grid->tri[jb2].t[0] == ib) grid->tri[jb2].t[0] = ia;
            if (grid->tri[jb2].t[1] == ib) grid->tri[jb2].t[1] = ia;
            if (grid->tri[jb2].t[2] == ib) grid->tri[jb2].t[2] = ia;
        }
        if (ja2 >= 0) {
            if (grid->tri[ja2].t[0] == ia) grid->tri[ja2].t[0] = ib;
            if (grid->tri[ja2].t[1] == ia) grid->tri[ja2].t[1] = ib;
            if (grid->tri[ja2].t[2] == ia) grid->tri[ja2].t[2] = ib;
        }

    } else {
        status = CANNOT_SWAP;
    }

cleanup:
    return status;
}


/***********************************************************************/
/*                                                                     */
/*   plotGrid - level 3 GRAFIC routine                                 */
/*                                                                     */
/***********************************************************************/

#ifdef GRAFIC
static void
plotGrid(int    *ifunct,                /* (in)  GRAFIC function indicator */
         void   *gridP,                 /* (in)  pointer to grid */
/*@unused@*/void   *a1,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a2,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a3,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a4,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a5,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a6,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a7,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a8,                    /* (in)  dummy GRAFIC argument */
/*@unused@*/void   *a9,                    /* (in)  dummy GRAFIC argument */
         float  *scale,                 /* (out) array of scales */
         char   *text,                  /* (out) help text */
         int    textlen)                /* (in)  length of text */
{
    grid_T  *grid   = (grid_T *) gridP;

    int     ipnt, itri, ip0, ip1, ip2;
    float   u4[3], v4[3];
    double  umin, umax, vmin, vmax;

    int     iblack  = GR_BLACK;
    int     igreen  = GR_GREEN;
    int     icircle = GR_CIRCLE;
    int     isquare = GR_SQUARE;
    int     ithree  = 3;

    /* --------------------------------------------------------------- */

    /* return scales */
    if (*ifunct == 0) {
        umin = grid->pnt[0].u;
        umax = grid->pnt[0].u;
        vmin = grid->pnt[0].v;
        vmax = grid->pnt[0].v;

        for (ipnt = 0; ipnt < grid->npnt; ipnt++) {
            if (grid->pnt[ipnt].u < umin) umin = grid->pnt[ipnt].u;
            if (grid->pnt[ipnt].u > umax) umax = grid->pnt[ipnt].u;
            if (grid->pnt[ipnt].v < vmin) vmin = grid->pnt[ipnt].v;
            if (grid->pnt[ipnt].v > vmax) vmax = grid->pnt[ipnt].v;
        }

        scale[0] = umin;
        scale[1] = umax;
        scale[2] = vmin;
        scale[3] = vmax;

        strcpy(text, "Grid Nearest");

    /* plot image */
    } else if (*ifunct == 1) {

        grcolr_(&igreen);
        for (ipnt = 0; ipnt < grid->nbnd; ipnt++) {
            u4[0] = grid->pnt[ipnt].u;
            v4[0] = grid->pnt[ipnt].v;

            grmov2_(u4, v4);
            if (ipnt == 0) {
                grsymb_(&isquare);
            } else {
                grsymb_(&icircle);
            }
        }

        grcolr_(&iblack);
        for (itri = 0; itri < grid->ntri; itri++) {
            ip0 = grid->tri[itri].p[0];
            ip1 = grid->tri[itri].p[1];
            ip2 = grid->tri[itri].p[2];

            u4[0] = grid->pnt[ip0].u;
            v4[0] = grid->pnt[ip0].v;
            u4[1] = grid->pnt[ip1].u;
            v4[1] = grid->pnt[ip1].v;
            u4[2] = grid->pnt[ip2].u;
            v4[2] = grid->pnt[ip2].v;

            grply2_(u4, v4, &ithree);
        }

        grcolr_(&iblack);

    } else {
        printf("Illegal option selected\n");
    }

//cleanup:
    return;
}
#endif
