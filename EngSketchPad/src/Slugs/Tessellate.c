/*
 ************************************************************************
 *                                                                      *
 * create and manage tessellations                                      *
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
#include <math.h>
#include <string.h>
#include <assert.h>

#include "egads.h"        // needed for EG_alloc, ...
#include "common.h"
#include "Tessellate.h"
#include "RedBlackTree.h"

/* definitions needed to read/write binary stl files */
#define UINT32 unsigned int
#define UINT16 unsigned short int
#define REAL32 float

#define MIN3(A,B,C)     (MIN(MIN(A,B),C))
#define MAX3(A,B,C)     (MAX(MAX(A,B),C))
#define ACOS(A)         acos(MINMAX(-1, (A), +1))

static  void            *realloc_temp = NULL;       /* used by RALLOC macro */

typedef struct {
    int     nrow;
    int     nent;
    int     ment;
    double  *a;
    int     *icol;
    int     *next;
} smf_T;

/* forward declarations of static routines defined below */
static int    connectNeighbors(tess_T *tess, int itri);
static int    dijkstra(tess_T *tess, int isrc, int itgt, int prev[], int link[]);
static double distance(tess_T *tess, int ipnt, int jpnt);
static int    eigen(double a[], int n, double eval[], double evec[]);
static int    smfAdd( smf_T *smf, int irow, int icol);
static int    smfFree(smf_T *smf);
static int    smfInit(smf_T *smf, int nrow);
static void   triNormal(tess_T *tess, int ip0, int ip1, int ip2, double *area, double norm[]);
static double turn(tess_T *tess, int ipnt, int jpnt, int kpnt, int itri);

static int    buildOctree(tess_T *tess, int nmax, oct_T *tree);
static int    refineOctree(oct_T *tree, int nmax, int depth);
static int    removeOctree(oct_T *tree);


/*
 ******************************************************************************
 *                                                                            *
 * addPoint - add a Point                                                     *
 *                                                                            *
 ******************************************************************************
 */
int
addPoint(tess_T  *tess,                 /* (in)  pointer to TESS */
         double  x,                     /* (in)  x-coordinate */
         double  y,                     /* (in)  y-coordinate */
         double  z)                     /* (inO  z-coordinate */
{
    int    status = 0;                  /* (out) return status */

    int    ipnt;

    ROUTINE(addPoint);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* make room for the new Point (if needed) */
    if (tess->npnt >= tess->mpnt-1) {
        (tess->mpnt) += 1000;
        RALLOC(tess->xyz,  double, 3*tess->mpnt);
        RALLOC(tess->uv,   double, 2*tess->mpnt);
        RALLOC(tess->ptyp, int,      tess->mpnt);
    }

    /* create the new Point */
    ipnt = tess->npnt;

    tess->xyz[3*ipnt  ] = x;
    tess->xyz[3*ipnt+1] = y;
    tess->xyz[3*ipnt+2] = z;

    tess->uv[2*ipnt  ] = 0;
    tess->uv[2*ipnt+1] = 0;

    tess->ptyp[ ipnt  ] = 0;

    (tess->npnt)++;

    /* return the new Point's index */
    status = ipnt;

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * addTriangle - add a Triangle                                               *
 *                                                                            *
 ******************************************************************************
 */
int
addTriangle(tess_T  *tess,              /* (in)  pointer to TESS */
            int     ip0,                /* (in)  index of first  Point (bias-0) */
            int     ip1,                /* (in)  index of second Point (bias-0) */
            int     ip2,                /* (in)  index of third  Point (bias-0) */
            int     it0,                /* (in)  index of first  Triangle (or -1) */
            int     it1,                /* (in)  index of first  Triangle (or -1) */
            int     it2)                /* (in)  index of first  Triangle (or -1) */
{
    int    status = 0;                  /* (out) return status */

    int    itri;

    ROUTINE(addTriangle);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (ip0 < 0 || ip0 >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    } else if (ip1 < 0 || ip1 >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    } else if (ip2 < 0 || ip2 >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    } else if (it0 < -1 || it0 >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (it1 < -1 || it1 >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (it2 < -1 || it2 >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    /* make room for the new Triangle (if needed) */
    if (tess->ntri >= tess->mtri-1) {
        (tess->mtri) += 1000;

        RALLOC(tess->trip, int,    3*tess->mtri);
        RALLOC(tess->trit, int,    3*tess->mtri);
        RALLOC(tess->ttyp, int,      tess->mtri);
        RALLOC(tess->bbox, double, 6*tess->mtri);
    }

    /* create the new Triangle */
    itri = tess->ntri;

    tess->trip[3*itri  ] = ip0;
    tess->trip[3*itri+1] = ip1;
    tess->trip[3*itri+2] = ip2;
    tess->trit[3*itri  ] = it0;
    tess->trit[3*itri+1] = it1;
    tess->trit[3*itri+2] = it2;
    tess->ttyp[  itri  ] = TRI_ACTIVE | TRI_VISIBLE;

    tess->bbox[6*itri  ] = MIN3(tess->xyz[3*ip0  ], tess->xyz[3*ip1  ], tess->xyz[3*ip2  ]);
    tess->bbox[6*itri+1] = MAX3(tess->xyz[3*ip0  ], tess->xyz[3*ip1  ], tess->xyz[3*ip2  ]);
    tess->bbox[6*itri+2] = MIN3(tess->xyz[3*ip0+1], tess->xyz[3*ip1+1], tess->xyz[3*ip2+1]);
    tess->bbox[6*itri+3] = MAX3(tess->xyz[3*ip0+1], tess->xyz[3*ip1+1], tess->xyz[3*ip2+1]);
    tess->bbox[6*itri+4] = MIN3(tess->xyz[3*ip0+2], tess->xyz[3*ip1+2], tess->xyz[3*ip2+2]);
    tess->bbox[6*itri+5] = MAX3(tess->xyz[3*ip0+2], tess->xyz[3*ip1+2], tess->xyz[3*ip2+2]);

    (tess->ntri)++;

    /* connect this Triangle with its neighbors */
    status = connectNeighbors(tess, itri);
    CHECK_STATUS(connectNeighbors);

    /* return the new Triangle's index */
    status = itri;

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * bridgeToPoint - create one Triangle that bridges gap between Triangle and Point *
 *                                                                            *
 ******************************************************************************
 */
int
bridgeToPoint(tess_T  *tess,            /* (in)  pointer to TESS */
              int     itri,             /* (in)  index of Triangle (bias-0) */
              int     ipnt)             /* (in)  index of Point    (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    ipa, ipb, jpa, jpb, itria, itrib, jtri, iside;
    double xbar, ybar, zbar, dbest, dtest;

    ROUTINE(bridgeTriangles);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (ipnt < 0 || ipnt >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    }

    /* find the two Points for the bridging in itri */
    ipa = -1;
    ipb = -1;
    dbest = HUGEQ;

    for (iside = 0; iside < 3; iside++) {
        if (tess->trit[3*itri+iside] < 0) {
            jpa   = tess->trip[3*itri+(iside+1)%3];
            jpb   = tess->trip[3*itri+(iside+2)%3];
            xbar  = (tess->xyz[3*jpa  ] + tess->xyz[3*jpb  ]) / 2;
            ybar  = (tess->xyz[3*jpa+1] + tess->xyz[3*jpb+1]) / 2;
            zbar  = (tess->xyz[3*jpa+2] + tess->xyz[3*jpb+2]) / 2;
            dtest = pow(tess->xyz[3*ipnt  ]-xbar, 2)
                  + pow(tess->xyz[3*ipnt+1]-ybar, 2)
                  + pow(tess->xyz[3*ipnt+2]-zbar, 2);
            if (dtest < dbest) {
                ipa   = jpa;
                ipb   = jpb;
                dbest = dtest;
            }
        }
    }

    if (ipa < 0 || ipb < 0) {
        printf("ERROR:: itri=%d does not adjoin a Loop\a\n", itri);
        goto cleanup;
    }

    /* identify itria and itrib if they exist */
    itria = -1;
    itrib = -1;

    for (jtri = 0; jtri < tess->ntri; jtri++) {
        if ((tess->ttyp[jtri] & TRI_ACTIVE) == 0) continue;

        for (iside = 0; iside < 3; iside++) {
            jpa = tess->trip[3*jtri+(iside+1)%3];
            jpb = tess->trip[3*jtri+(iside+2)%3];

            if        (jpa == ipb  && jpb == ipnt) {
                itria = jtri;
            } else if (jpa == ipnt && jpb == ipa ) {
                itrib = jtri;
            }
        }
    }

    /* create the Triangle connecting ipnt, ipb, and ipa */
    status = addTriangle(tess, ipnt, ipb,   ipa,
                               itri, itrib, itria);
    CHECK_STATUS(addTriangle);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * bridgeTriangles - create two Triangles that bridge gap between given Triangles *
 *                                                                            *
 ******************************************************************************
 */
int
bridgeTriangles(tess_T  *tess,          /* (in)  pointer to TESS */
                int     itri,           /* (in)  index of first  Triangle (bias-0) */
                int     jtri)           /* (in)  index of second Triangle (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    ipa, ipb, jpa, jpb;

    ROUTINE(bridgeTriangles);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (jtri < 0 || jtri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    /* find the two Points for the bridging in itri */
    if        (tess->trit[3*itri  ] < 0) {
        ipa = tess->trip[3*itri+1];
        ipb = tess->trip[3*itri+2];
    } else if (tess->trit[3*itri+1] < 0) {
        ipa = tess->trip[3*itri+2];
        ipb = tess->trip[3*itri  ];
    } else if (tess->trit[3*itri+2] < 0) {
        ipa = tess->trip[3*itri  ];
        ipb = tess->trip[3*itri+1];
    } else {
        printf("ERROR:: itri=%d does not adjoin a Loop\a\n", itri);
        goto cleanup;
    }

    /* find the two Points for the bridging in jtri */
    if        (tess->trit[3*jtri  ] < 0) {
        jpa = tess->trip[3*jtri+1];
        jpb = tess->trip[3*jtri+2];
    } else if (tess->trit[3*jtri+1] < 0) {
        jpa = tess->trip[3*jtri+2];
        jpb = tess->trip[3*jtri  ];
    } else if (tess->trit[3*jtri+2] < 0) {
        jpa = tess->trip[3*jtri  ];
        jpb = tess->trip[3*jtri+1];
    } else {
        printf("ERROR:: jtri=%d does not adjoin a Loop\a\n", jtri);
        goto cleanup;
    }

    /* create the Triangle connecting ipa, jpb, and jpa */
    status = addTriangle(tess, ipa,  jpb,          jpa,
                               jtri, tess->ntri+1, -1 );
    CHECK_STATUS(addTriangle);

    /* create the Triangle connecting jpa, ipb, and ipa */
    status = addTriangle(tess, jpa,  ipb,          ipa,
                               itri, tess->ntri-1, -1 );
    CHECK_STATUS(addTriangle);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * checkAreas - check areas in UV                                             *
 *                                                                            *
 ******************************************************************************
 */
int
checkAreas(tess_T  *tess,               /* (in)  pointer to TESS */
           int     *nneg,               /* (out) number of negative areas */
           int     *npos)               /* (out) number of positive area */
{
    int    status = 0;                  /* (out) return status */

    int      itri;
    double   u0, v0, u1, v1, u2, v2, area;

//    ROUTINE(checkAreas);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* find the number of Triangles with negative and positive areas */
    *nneg = 0;
    *npos = 0;

    for (itri = 0; itri < tess->ntri; itri++) {
        u0 = tess->uv[2*tess->trip[3*itri  ]  ];
        v0 = tess->uv[2*tess->trip[3*itri  ]+1];
        u1 = tess->uv[2*tess->trip[3*itri+1]  ];
        v1 = tess->uv[2*tess->trip[3*itri+1]+1];
        u2 = tess->uv[2*tess->trip[3*itri+2]  ];
        v2 = tess->uv[2*tess->trip[3*itri+2]+1];

        area = (u1 - u0) * (v2 - v0) - (v1 - v0) * (u2 - u0);

        if (area < 0) (*nneg)++;
        if (area > 0) (*npos)++;
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * colorTriangles - color a Triangle and its neighbors (up to links)          *
 *                                                                            *
 ******************************************************************************
 */
int
colorTriangles(tess_T  *tess,           /* (in)  pointer to TESS */
               int     itri,            /* (in)  index of Triangle (bias-0) */
               int     icolr)           /* (in)  color (0-255) */
{
    int    status = 0;                  /* (out) return status */

    int    jcolr, jtri, ktri, nchange;

//    ROUTINE(colorTriangle);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (icolr < 0 || icolr > 255) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* do nothing if old and new colors match */
    jcolr = tess->ttyp[itri] & TRI_COLOR;
    if (icolr == jcolr) {
        goto cleanup;
    }

    /* color the first Triangle */
    tess->ttyp[itri] = (tess->ttyp[itri] & ~TRI_COLOR) | icolr;

    /* recursively flood-fill up to Sides that are links */
    for (ktri = 0; ktri < tess->ntri+10; ktri++) {
        nchange = 0;

        for (jtri = 0; jtri < tess->ntri; jtri++) {
            if ((tess->ttyp[jtri] & TRI_ACTIVE) == 0    ) continue;
            if ((tess->ttyp[jtri] & TRI_COLOR)  != icolr) continue;

            if ((tess->ttyp[jtri] & TRI_T0_LINK) == 0) {
                ktri = tess->trit[3*jtri  ];
                if (ktri >= 0) {
                    if ((tess->ttyp[ktri] & TRI_COLOR) == jcolr) {
                        tess->ttyp[ktri] = (tess->ttyp[ktri] & ~TRI_COLOR) | icolr;
                        nchange++;
                    }
                }
            }

            if ((tess->ttyp[jtri] & TRI_T1_LINK) == 0) {
                ktri = tess->trit[3*jtri+1];
                if (ktri >= 0) {
                    if ((tess->ttyp[ktri] & TRI_COLOR) == jcolr) {
                        tess->ttyp[ktri] = (tess->ttyp[ktri] & ~TRI_COLOR) | icolr;
                        nchange++;
                    }
                }
            }

            if ((tess->ttyp[jtri] & TRI_T2_LINK) == 0) {
                ktri = tess->trit[3*jtri+2];
                if (ktri >= 0) {
                    if ((tess->ttyp[ktri] & TRI_COLOR) == jcolr) {
                        tess->ttyp[ktri] = (tess->ttyp[ktri] & ~TRI_COLOR) | icolr;
                        nchange++;
                    }
                }
            }
        }

        if (nchange == 0) goto cleanup;
    }

    /* gettting here means that we never broke out of above loop */
    printf("ERROR:: reached maximum iterations while coloring\n");
    status = TESS_NOT_CONVERGED;

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * connectNeighbors - set up the neighbor info for the neighbors of the given Triangle *
 *                                                                            *
 ******************************************************************************
 */
static int
connectNeighbors(tess_T  *tess,         /* (in)  pointer to TESS */
                 int     itri)          /* (in)  Triangle index (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    ipnt, jtri;

//    ROUTINE(connectNeighbors);

    /*----------------------------------------------------------------*/

    ipnt = tess->trip[3*itri  ];
    jtri = tess->trit[3*itri+1];
    if (jtri >= 0) {
        if        (tess->trip[3*jtri  ] == ipnt) {
            tess->trit[3*jtri+2] = itri;
        } else if (tess->trip[3*jtri+1] == ipnt) {
            tess->trit[3*jtri  ] = itri;
        } else if (tess->trip[3*jtri+2] == ipnt) {
            tess->trit[3*jtri+1] = itri;
        } else {
            printf("ERROR:: Trouble stitching things up (A)\a\n");
            status = TESS_INTERNAL_ERROR;
            goto cleanup;
        }
    }

    ipnt = tess->trip[3*itri+1];
    jtri = tess->trit[3*itri+2];
    if (jtri >= 0) {
        if        (tess->trip[3*jtri  ] == ipnt) {
            tess->trit[3*jtri+2] = itri;
        } else if (tess->trip[3*jtri+1] == ipnt) {
            tess->trit[3*jtri  ] = itri;
        } else if (tess->trip[3*jtri+2] == ipnt) {
            tess->trit[3*jtri+1] = itri;
        } else {
            printf("ERROR:: Trouble stitching things up (B)\a\n");
            status = TESS_INTERNAL_ERROR;
            goto cleanup;
        }
    }

    ipnt = tess->trip[3*itri+2];
    jtri = tess->trit[3*itri  ];
    if (jtri >= 0) {
        if        (tess->trip[3*jtri  ] == ipnt) {
            tess->trit[3*jtri+2] = itri;
        } else if (tess->trip[3*jtri+1] == ipnt) {
            tess->trit[3*jtri  ] = itri;
        } else if (tess->trip[3*jtri+2] == ipnt) {
            tess->trit[3*jtri+1] = itri;
        } else {
            printf("ERROR:: Trouble stitching things up (C)\a\n");
            status = TESS_INTERNAL_ERROR;
            goto cleanup;
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * copyTess - copy a Tessellation                                             *
 *                                                                            *
 ******************************************************************************
 */
int
copyTess(tess_T  *src,                  /* (in)  pointer to source TESS */
         tess_T  *tgt)                  /* (in)  pointer to target TESS */
{
    int    status = 0;                  /* (out) return status */

    ROUTINE(copyTess);

    /* --------------------------------------------------------------- */

    if (src == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (src->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* if the target is already a TESS, free up its arrays */
    if (tgt->magic == TESS_MAGIC) {
        FREE(tgt->trip);
        FREE(tgt->trit);
        FREE(tgt->ttyp);
        FREE(tgt->bbox);
        FREE(tgt->xyz );
        FREE(tgt->uv  );
        FREE(tgt->ptyp);

    /* otherwise nullify th epointers */
    } else {
        tgt->trip = NULL;
        tgt->trit = NULL;
        tgt->ttyp = NULL;
        tgt->bbox = NULL;
        tgt->xyz  = NULL;
        tgt->uv   = NULL;
        tgt->ptyp = NULL;
    }

    /* initialize */
    tgt->magic  = src->magic;
    tgt->ntri   = src->ntri;
    tgt->mtri   = src->mtri;
    tgt->nhang  = src->nhang;
    tgt->nlink  = src->nlink;
    tgt->ncolr  = src->ncolr;
    tgt->npnt   = src->npnt;
    tgt->mpnt   = src->mpnt;
    tgt->octree = NULL;

    MALLOC(tgt->trip, int,    3*tgt->mtri);
    MALLOC(tgt->trit, int,    3*tgt->mtri);
    MALLOC(tgt->ttyp, int,      tgt->mtri);
    MALLOC(tgt->bbox, double, 6*tgt->mtri);
    MALLOC(tgt->xyz,  double, 3*tgt->mpnt);
    MALLOC(tgt->uv,   double, 2*tgt->mpnt);
    MALLOC(tgt->ptyp, int,      tgt->mpnt);

    /* copy in the data */
    memcpy(tgt->trip, src->trip, 3*tgt->mtri*sizeof(int   ));
    memcpy(tgt->trit, src->trit, 3*tgt->mtri*sizeof(int   ));
    memcpy(tgt->ttyp, src->ttyp,   tgt->mtri*sizeof(int   ));
    memcpy(tgt->bbox, src->bbox, 6*tgt->mtri*sizeof(double));
    memcpy(tgt->xyz,  src->xyz,  3*tgt->mpnt*sizeof(double));
    memcpy(tgt->uv,   src->uv,   2*tgt->mpnt*sizeof(double));
    memcpy(tgt->ptyp, src->ptyp,   tgt->mpnt*sizeof(int   ));

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * createLink - create a Link on one side of a Triangle                       *
 *                                                                            *
 ******************************************************************************
 */
int
createLink(tess_T  *tess,               /* (in)  pointer to TESS */
           int     itri,                /* (in)  index of Triangle (bias-0) */
           int     isid)                /* (in)  side index (0-2) */
{
    int    status = 0;                  /* (out) return status */

    int    jtri;

//    ROUTINE(createLink);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (isid < 0 || isid >= 3) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* return immediately if Link already exists */
    if (isid == 0 && (tess->ttyp[itri] & TRI_T0_LINK) != 0) goto cleanup;
    if (isid == 1 && (tess->ttyp[itri] & TRI_T1_LINK) != 0) goto cleanup;
    if (isid == 2 && (tess->ttyp[itri] & TRI_T2_LINK) != 0) goto cleanup;

    /* create the Links if the companion Triangle exists */
    jtri = tess->trit[3*itri+isid];

    if (jtri >= 0) {
        if (isid == 0) tess->ttyp[itri] |= TRI_T0_LINK;
        if (isid == 1) tess->ttyp[itri] |= TRI_T1_LINK;
        if (isid == 2) tess->ttyp[itri] |= TRI_T2_LINK;

        if (tess->trit[3*jtri  ] == itri) tess->ttyp[jtri] |= TRI_T0_LINK;
        if (tess->trit[3*jtri+1] == itri) tess->ttyp[jtri] |= TRI_T1_LINK;
        if (tess->trit[3*jtri+2] == itri) tess->ttyp[jtri] |= TRI_T2_LINK;

        (tess->nlink)++;
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * createLinks - create Links between given Points                            *
 *                                                                            *
 ******************************************************************************
 */
int
createLinks(tess_T  *tess,              /* (in)  pointer to TESS */
            int     isrc,               /* (in)  index of source Point (bias-0) */
            int     itgt)               /* (in)  index of target Point (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    ipnt, ipm1, itri;
    int    *prev = NULL, *link = NULL;

    ROUTINE(createLinks);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* allocate storage for Dijkstra's algorithm */
    MALLOC(prev, int,    tess->npnt);
    MALLOC(link, int,    tess->npnt);

    /* find the path via dijkstra */
    status = dijkstra(tess, isrc, itgt, prev, link);
    CHECK_STATUS(dijkstra);

    /* create the Links by traversing from the target to the source */
    ipnt = itgt;
    while (prev[ipnt] >= 0) {
        itri = link[ipnt];
        ipm1 = prev[ipnt];

        if        ((tess->trip[3*itri+1] == ipm1 && tess->trip[3*itri+2] == ipnt) ||
                   (tess->trip[3*itri+2] == ipm1 && tess->trip[3*itri+1] == ipnt)   ) {
            if ((tess->ttyp[itri] & TRI_T0_LINK) == 0) {
                tess->ttyp[itri] |= TRI_T0_LINK;
                (tess->nlink)++;
            }
                itri = tess->trit[3*itri  ];
        } else if ((tess->trip[3*itri+2] == ipm1 && tess->trip[3*itri  ] == ipnt) ||
                   (tess->trip[3*itri  ] == ipm1 && tess->trip[3*itri+2] == ipnt)   ) {
            if ((tess->ttyp[itri] & TRI_T1_LINK) == 0) {
                tess->ttyp[itri] |= TRI_T1_LINK;
                (tess->nlink)++;
            }
            itri = tess->trit[3*itri+1];
        } else if ((tess->trip[3*itri  ] == ipm1 && tess->trip[3*itri+1] == ipnt) ||
                   (tess->trip[3*itri+1] == ipm1 && tess->trip[3*itri  ] == ipnt)   ) {
            if ((tess->ttyp[itri] & TRI_T2_LINK) == 0) {
                tess->ttyp[itri] |= TRI_T2_LINK;
                (tess->nlink)++;
            }
            itri = tess->trit[3*itri+2];
        }

        if (itri >= 0) {
            if        ((tess->trip[3*itri+1] == ipm1 && tess->trip[3*itri+2] == ipnt) ||
                       (tess->trip[3*itri+2] == ipm1 && tess->trip[3*itri+1] == ipnt)   ) {
                tess->ttyp[itri] |= TRI_T0_LINK;
            } else if ((tess->trip[3*itri+2] == ipm1 && tess->trip[3*itri  ] == ipnt) ||
                       (tess->trip[3*itri  ] == ipm1 && tess->trip[3*itri+2] == ipnt)   ) {
                tess->ttyp[itri] |= TRI_T1_LINK;
            } else if ((tess->trip[3*itri  ] == ipm1 && tess->trip[3*itri+1] == ipnt) ||
                       (tess->trip[3*itri+1] == ipm1 && tess->trip[3*itri  ] == ipnt)   ) {
                tess->ttyp[itri] |= TRI_T2_LINK;
            }
        }

        ipnt = ipm1;
    }

cleanup:
    FREE(link);
    FREE(prev);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * cutTriangles - cut Triangles through given Points                          *
 *                                                                            *
 ******************************************************************************
 */
int
cutTriangles(tess_T  *tess,             /* (in)  pointer to TESS */
             int     icolr,             /* (in)  color (or -1 for all) */
             double  data[])            /* (in)  cut is data[0]+x*data[1]+y*data[2]+z*data[3] */
{
    int    status = 0;                  /* (out) return status */

    int    ipnt, itri, jtri, isid, jsid, mcut;
    int    ip0, ip1, ip2, ip3, ipnew, it2, it3;
    double frac, xx, yy, zz, *cut = NULL;

    ROUTINE(cutTriangles);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* create a rotated copy of the Points, such that cut=0
       is where Triangle should be cut */
    mcut = tess->npnt + 100;
    MALLOC(cut, double, mcut);

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        cut[ipnt] = data[0] + tess->xyz[3*ipnt  ] * data[1]
                            + tess->xyz[3*ipnt+1] * data[2]
                            + tess->xyz[3*ipnt+2] * data[3];
    }

    /* loop through all Triangles or the current color, looking for a
       side that straddles cut=0 */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;
        if (icolr >= 0 && (tess->ttyp[itri] & TRI_COLOR) != icolr) continue;

        /*
                 ip3                      ip3
                 / \                      /:\
                /   \                    / : \
           it0 /     \ it3          it0 /  :  \ it3
              /  jtri \                /   :   \
             /         \              /    :    \
            /           \            / jtri: n-1 \
          ip2-----------ip1   ==>  ip2---ipnew---ip1
            \           /            \ itri: n-2 /
             \         /              \    :    /
              \  itri /                \   :   /
           it1 \     / it2          it1 \  :  / it2
                \   /                    \ : /
                 \ /                      \:/
                 ip0                      ip0
        */

        for (isid = 0; isid < 3; isid++) {
            ip0 = tess->trip[3*itri+(isid  )%3];
            ip1 = tess->trip[3*itri+(isid+1)%3];
            ip2 = tess->trip[3*itri+(isid+2)%3];
            ip3 = -1;

            if ((cut[ip1] < 0 && cut[ip2] > 0) ||
                (cut[ip1] > 0 && cut[ip2] < 0)   ) {
                jtri = tess->trit[3*itri+(isid  )%3];
//                it0  = -1;
//                it1  = tess->trit[3*itri+(isid+1)%3];
                it2  = tess->trit[3*itri+(isid+2)%3];
                it3  = -1;

                for (jsid = 0; jsid < 3; jsid++) {
                    if (tess->trit[3*jtri+jsid] == itri) {
                        ip3 = tess->trip[3*jtri+(jsid  )%3];
                        it3 = tess->trit[3*jtri+(jsid+1)%3];
//                        it0 = tess->trit[3*jtri+(jsid+2)%3];
                        break;
                    }
                }

                if (it3 < 0 || ip3 < 0) {
                    break;
                }

                /* add a Point at the crossing (if not too close to
                   an existing Point */
                frac = cut[ip1] / (cut[ip1] - cut[ip2]);
                if (frac < EPS06 || frac > 1-EPS06) continue;

                xx   = (1-frac) * tess->xyz[3*ip1  ] + frac * tess->xyz[3*ip2  ];
                yy   = (1-frac) * tess->xyz[3*ip1+1] + frac * tess->xyz[3*ip2+1];
                zz   = (1-frac) * tess->xyz[3*ip1+2] + frac * tess->xyz[3*ip2+2];

                status = addPoint(tess, xx, yy, zz);
                CHECK_STATUS(addPoint);

                ipnew = tess->npnt - 1;

                if (ipnew >= mcut) {
                    mcut += 100;
                    RALLOC(cut, double, mcut);
                }

                cut[ipnew] = 0;

                /* modify itri and jtri */
                tess->trip[3*itri+(isid+1)%3] = ipnew;
                tess->trit[3*itri+(isid+2)%3] = -1;

                tess->trip[3*jtri+(jsid+2)%3] = ipnew;
                tess->trit[3*jtri+(jsid+1)%3] = -1;

                /* create the new Triangles (and hold off neighbor
                   information amongst them) */
                status = addTriangle(tess, ipnew, ip0, ip1, it2, -1, -1);
                CHECK_STATUS(addTriangle);

                status = addTriangle(tess, ipnew, ip1, ip3, it3, -1, -1);
                CHECK_STATUS(addTriangle);

                /* now that all Triangles are made, set up neighbor information */
                tess->trit[3*itri+(isid+2)%3] = tess->ntri - 2;
                tess->trit[3*jtri+(jsid+1)%3] = tess->ntri - 1;

                tess->trit[3*(tess->ntri-2)+1] = tess->ntri - 1;
                tess->trit[3*(tess->ntri-2)+2] = itri;

                tess->trit[3*(tess->ntri-1)+1] = jtri;
                tess->trit[3*(tess->ntri-1)+2] = tess->ntri - 2;

                /* color the new Triangles */
                tess->ttyp[tess->ntri-2] = (tess->ttyp[tess->ntri-2] & ~ TRI_COLOR) | (tess->ttyp[itri] & TRI_COLOR);
                tess->ttyp[tess->ntri-1] = (tess->ttyp[tess->ntri-1] & ~ TRI_COLOR) | (tess->ttyp[jtri] & TRI_COLOR);

                /* decreast itri so that it gets looked at again */
                itri--;
            }
        }
    }

cleanup:
    FREE(cut);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * deleteTriangle - delete a Triangle                                         *
 *                                                                            *
 ******************************************************************************
 */
int
deleteTriangle(tess_T  *tess,           /* (in)  pointer to TESS */
               int     itri)            /* (in)  index of Triangle to delete (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    jtri;

//    ROUTINE(deleteTriangle);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    /* mark the Triangle as deleted */
    tess->ttyp[itri] &= ~(TRI_ACTIVE | TRI_VISIBLE);

    /* remove the neighbor pointers from the neighboring Triangles */
    jtri = tess->trit[3*itri  ];
    if (jtri >= 0) {
        if (tess->trit[3*jtri  ] == itri) tess->trit[3*jtri  ] = -1;
        if (tess->trit[3*jtri+1] == itri) tess->trit[3*jtri+1] = -1;
        if (tess->trit[3*jtri+2] == itri) tess->trit[3*jtri+2] = -1;
    }

    jtri = tess->trit[3*itri+1];
    if (jtri >= 0) {
        if (tess->trit[3*jtri  ] == itri) tess->trit[3*jtri  ] = -1;
        if (tess->trit[3*jtri+1] == itri) tess->trit[3*jtri+1] = -1;
        if (tess->trit[3*jtri+2] == itri) tess->trit[3*jtri+2] = -1;
    }

    jtri = tess->trit[3*itri+2];
    if (jtri >= 0) {
        if (tess->trit[3*jtri  ] == itri) tess->trit[3*jtri  ] = -1;
        if (tess->trit[3*jtri+1] == itri) tess->trit[3*jtri+1] = -1;
        if (tess->trit[3*jtri+2] == itri) tess->trit[3*jtri+2] = -1;
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * dijkstra - find shorest path between isrc and itgt                         *
 *                                                                            *
 ******************************************************************************
 */
static int
dijkstra(tess_T  *tess,                 /* (in)  pointer to TESS */
         int     isrc,                  /* (in)  index of source Point (bias-0) */
         int     itgt,                  /* (in)  index of target Point (bias-0) */
         int     prev[],                /* (out) previous Points */
         int     link[])                /* (out) Triangle back to prev */
{
    int    status = 0;                  /* (out) return status */

    int    ipnt, itri, ip0, ip1, ip2, ipass, nchange;
    double d01, d12, d20, dmin, dmax;
    double *dist = NULL;

    ROUTINE(dijkstra);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* allocate storage for distances from isrc */
    MALLOC(dist, double, tess->npnt);

    dmax = 2 * distance(tess, isrc, itgt);

    /* initialize distance to all Points */
    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        prev[ipnt] = -1;
        link[ipnt] = -1;
        dist[ipnt] = HUGEQ;
    }

    dist[isrc] = 0;

    /* make passes through Triangles until no distances are updated */
    for (ipass = 0; ipass < tess->ntri; ipass++) {
        nchange = 0;

        for (itri = 0; itri < tess->ntri; itri++) {
            if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

            ip0 = tess->trip[3*itri  ];
            ip1 = tess->trip[3*itri+1];
            ip2 = tess->trip[3*itri+2];

            if (dist[ip0] < dmax || dist[ip1] < dmax) {
                d01 = distance(tess, ip0, ip1);

                dmin = dist[ip1] + d01;
                if (dmin < dist[ip0]) {
                    dist[ip0] = dmin;
                    link[ip0] = itri;
                    prev[ip0] = ip1;
                    nchange++;
                }

                dmin = dist[ip0] + d01;
                if (dmin < dist[ip1]) {
                    dist[ip1] = dmin;
                    link[ip1] = itri;
                    prev[ip1] = ip0;
                }
            }

            if (dist[ip1] < dmax || dist[ip2] < dmax) {
                d12 = distance(tess, ip1, ip2);

                dmin = dist[ip2] + d12;
                if (dmin < dist[ip1]) {
                    dist[ip1] = dmin;
                    link[ip1] = itri;
                    prev[ip1] = ip2;
                    nchange++;
                }

                dmin = dist[ip2] + d12;
                if (dmin < dist[ip2]) {
                    dist[ip2] = dmin;
                    link[ip2] = itri;
                    prev[ip2] = ip1;
                    nchange++;
                }
            }

            if (dist[ip2] < dmax || dist[ip0] < dmax) {
                d20 = distance(tess, ip2, ip0);

                dmin = dist[ip0] + d20;
                if (dmin < dist[ip2]) {
                    dist[ip2] = dmin;
                    link[ip2] = itri;
                    prev[ip2] = ip0;
                    nchange++;
                }

                dmin = dist[ip2] + d20;
                if (dmin < dist[ip0]) {
                    dist[ip0] = dmin;
                    link[ip0] = itri;
                    prev[ip0] = ip2;
                    nchange++;
                }
            }
        }

        if (nchange == 0) {
            break;
        }
    }

cleanup:
    FREE(dist);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * distance - find the distance between two Points                            *
 *                                                                            *
 ******************************************************************************
 */
static double
distance(tess_T  *tess,                 /* (in)  pointer to TESS */
         int     ipnt,                  /* (in)  first  Point index (bias-0) */
         int     jpnt)                  /* (in)  second Point index (bias-0) */
{
    double answer;                      /* (out) distance */

//    ROUTINE(distance);

    /*----------------------------------------------------------------*/

    answer = sqrt( SQR(tess->xyz[3*ipnt  ] - tess->xyz[3*jpnt  ])
                 + SQR(tess->xyz[3*ipnt+1] - tess->xyz[3*jpnt+1])
                 + SQR(tess->xyz[3*ipnt+2] - tess->xyz[3*jpnt+2]));

//cleanup:
    return answer;
}


/*
 ******************************************************************************
 *                                                                            *
 * eigen - find eigenvalues and eigenvectors of real symmetric matrix         *
 *                                                                            *
 ******************************************************************************
 */
static int
eigen(double a[],                       /* (in)  (assumed) symmetric matrix, stored by row */
                                        /* (out) mangled form of a */
      int    n,                         /* (in)  order of a */
      double eval[],                    /* (out) sorted eigenvalues (large to small) */
      double evec[])                    /* (out) sorted eigenvectors (columns of evec) */
{
    int status = 0;                     /* (out) default return value */

    #define A(   I,J) a[   I*n+J]
    #define EVEC(I,J) evec[I*n+J]

    int i, p, q, r, count, isweep, qmax;
    double theta, sintht, costht, tantht, tau, temp1, temp2, sum, dmax, swap;

    /* --------------------------------------------------------------- */

    /* note that we assume that a is symmetric and only use the
       diaginal and super-diagonal elements */

//    ROUTINE(eigen);

    /*----------------------------------------------------------------*/

    /* initialize eval (the eigenvalues) to the diagonal of a */
    for (p = 0; p < n; p++) {
        eval[p] = A(p,p);
    }

    /* initialize the evec matrix to the identity matrix */
    for (p = 0; p < n; p++) {
        for (q = 0; q < n; q++) {
            if (p == q) {
                EVEC(p,q) = 1;
            } else {
                EVEC(p,q) = 0;
            }
        }
    }

    /* take up to 50 sweeps through the matirx */
    for (isweep = 0; isweep < 50; isweep++) {

        /* determine if we are done by looking at the super-diagonal elements */
        count = 0;
        for (p = 0; p < n-1; p++) {
            for (q = p+1; q < n; q++) {
                if (A(p,q) != 0) {
                    count++;
                }
            }
        }

        /* if they are all zero, we are done */
        if (count == 0) {

            /* normalize the eigenvectors (stored in the columns of evec) */
            for (q = 0; q < n; q++) {
                sum = 0;
                for (p = 0; p < n; p++) {
                    sum += EVEC(p,q) * EVEC(p,q);
                }
                for (p = 0; p < n; p++) {
                    EVEC(p,q) /= sqrt(sum);
                }
            }

            /* order the eigenvalues (and associated eigenvectors) from largest to
               smallest eigenvalue with a simple insertion sort */
            for (q = 0; q < n-1; q++) {
                qmax = q;
                dmax = fabs(eval[q]);
                for (i = q+1; i < n; i++) {
                    if (fabs(eval[i]) > dmax) {
                        qmax = i;
                        dmax = fabs(eval[i]);
                    }
                }

                if (qmax != q) {
                    swap       = eval[q   ];
                    eval[q   ] = eval[qmax];
                    eval[qmax] = swap;

                    for (p = 0; p < n; p++) {
                        swap         = EVEC(p,q   );
                        EVEC(p,q   ) = EVEC(p,qmax);
                        EVEC(p,qmax) = swap;
                    }
                }
            }

            goto cleanup;
        }

        /* perform jacobi rotations in super-diagonal part of a (across rows, starting at top) */
        for (p = 0; p < n-1; p++) {
            for (q = p+1; q < n; q++) {

                if (A(p,q) == 0) continue;

                /* find theta that annihilates A(p,q) */
                theta = (eval[q] - eval[p]) / 2 / A(p,q);

                if (fabs(theta) > 1e10) {
                    tantht = +1 / (2 * theta);
                } else if (theta > 0) {
                    tantht = +1 / (+theta + sqrt(1 + theta * theta));
                } else {
                    tantht = -1 / (-theta + sqrt(1 + theta * theta));
                }

                costht = 1 / sqrt(1 + tantht * tantht);
                sintht = tantht * costht;
                tau    = sintht / (1 + costht);

                /* update diagonal (and eigenvalues) */
                A(p,p)  -= tantht * A(p,q);
                A(q,q)  += tantht * A(p,q);
                eval[p] -= tantht * A(p,q);
                eval[q] += tantht * A(p,q);

                /* perform jacobi rotation to annihilate A(p,q)*/
                A(p,q) = 0;

                /* columns p and q above row p */
                for (r = 0; r < p; r++) {
                    temp1  = A(r,p);
                    temp2  = A(r,q);
                    A(r,p) = temp1 - sintht * (temp2 + temp1 * tau);
                    A(r,q) = temp2 + sintht * (temp1 - temp2 * tau);
                }

                /* row p between columns p and q and column q between ros p and q */
                for (r = p+1; r < q; r++) {
                    temp1  = A(p,r);
                    temp2  = A(r,q);
                    A(p,r) = temp1 - sintht * (temp2 + temp1 * tau);
                    A(r,q) = temp2 + sintht * (temp1 - temp2 * tau);
                }

                /* rows p and q after column q */
                for (r = q+1; r < n; r++) {
                    temp1  = A(p,r);
                    temp2  = A(q,r);
                    A(p,r) = temp1 - sintht * (temp2 + temp1 * tau);
                    A(q,r) = temp2 + sintht * (temp1 - temp2 * tau);
                }

                /* accumualating this rotation into the eigenvectors */
                for (r = 0; r < n; r++) {
                    temp1     = EVEC(r,p);
                    temp2     = EVEC(r,q);
                    EVEC(r,p) = temp1 - sintht * (temp2 + temp1 * tau);
                    EVEC(r,q) = temp2 + sintht * (temp1 - temp2 * tau);
                }
            }
        }
    }

    /* we did not converge */
    status = TESS_NOT_CONVERGED;

cleanup:
    #undef EVEC
    #undef A

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * detectCreases - create links at creases                                    *
 *                                                                            *
 ******************************************************************************
 */
int
detectCreases(tess_T  *tess,            /* (in)  pointer to TESS */
              double  angdeg)           /* (in)  maximum crease angle (deg) */
{
    int    status = 0;                  /* (out) return status */

    int     itri, isid, jtri;
    double  cosang, area, normi[3], normj[3], dot;

    ROUTINE(detectCreases);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    cosang = cos(angdeg * PIo180);

    /* loop through all Triangle pairs */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

        for (isid = 0; isid < 3; isid++) {
            jtri = tess->trit[3*itri+isid];
            if (jtri < itri                         ) continue;
            if ((tess->ttyp[jtri] & TRI_ACTIVE) == 0) continue;

            /* find the two normals and their dot product */
            triNormal(tess, tess->trip[3*itri  ],
                            tess->trip[3*itri+1],
                            tess->trip[3*itri+2], &area, normi);

            triNormal(tess, tess->trip[3*jtri  ],
                            tess->trip[3*jtri+1],
                            tess->trip[3*jtri+2], &area, normj);

            dot = normi[0] * normj[0] + normi[1] * normj[1] + normi[2] * normj[2];

            /* if the dot product is less than the tolerace, create a Link
               on this side */
            if (dot < cosang) {
                status = createLink(tess, itri, isid);
                CHECK_STATUS(createLink);
            }
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * extendLoop - extend loop to given x/y/z                                    *
 *                                                                            *
 ******************************************************************************
 */
int
extendLoop(tess_T  *tess,               /* (in)  pointer to TESS */
           int     ipnt,                /* (in)  index of one Point on loop (bias-0) */
           int     itype,               /* (in)  type: 1=x, 2=y, 3=z */
           double  val)                 /* (in)  constant x/y/z for extension */
{
    int    status = 0;                  /* (out) return status */

    int     nseg, iseg, ip1, jp1, itri, jpnt, jtri, npnt_save, ntri_save;
    seg_T   *seg = NULL;

    ROUTINE(extendLoop);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (ipnt < 0 || ipnt >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    } else if (itype < 0 || itype > 3) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* get the loop following the hanging sides from ipnt */
    status = getLoop(tess, ipnt, &nseg, &seg);
    CHECK_STATUS(getLoop);

    /* add a copy of each of the Points in the loop */
    npnt_save = tess->npnt;
    ntri_save = tess->ntri;

    for (iseg = 0; iseg < nseg; iseg++) {
        ipnt = seg[iseg].pnt;

        if        (itype == 1) {
            status = addPoint(tess, val, tess->xyz[3*ipnt+1], tess->xyz[3*ipnt+2]);
            CHECK_STATUS(addPoint);
        } else if (itype == 2) {
            status = addPoint(tess, tess->xyz[3*ipnt  ], val, tess->xyz[3*ipnt+2]);
            CHECK_STATUS(addPoint);
        } else {
            status = addPoint(tess, tess->xyz[3*ipnt  ], tess->xyz[3*ipnt+1], val);
            CHECK_STATUS(addPoint);
        }
    }

    /* create Triangles that connect the old loop Points with the
       newly created Points (above)

                      jpnt       jp1
                        *---------*
                        | \       |
                   itri |   \  2  | jtri
                        |  1  \   |
                        |       \ |
                   =====*=========*=====
                      ipnt       ip1
                      seg[iseg].tri
    */

    for (iseg = 0; iseg < nseg; iseg++) {
        ipnt = seg[iseg].pnt;
        jpnt = npnt_save + iseg;

        if (iseg == 0) {
            ip1  = seg[iseg+1].pnt;
            jp1  = jpnt + 1;
            itri = -1;
            jtri = -1;
        } else if (iseg < nseg-1) {
            ip1  = seg[iseg+1].pnt;
            jp1  = jpnt + 1;
            itri = tess->ntri - 1;
            jtri = -1;
        } else {
            ip1  = seg[0     ].pnt;
            jp1  = npnt_save;
            itri = tess->ntri - 1;
            jtri = ntri_save;
        }

        status = addTriangle(tess, ip1,  ipnt, jpnt,
                                   itri, -1,   seg[iseg].tri);
        CHECK_STATUS(addTriangle);

        status = addTriangle(tess, jpnt, jp1,          ip1,
                                   jtri, tess->ntri-1, -1 );
        CHECK_STATUS(addTriangle);
    }

cleanup:
    FREE(seg);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * extractColor - create a new Tessellation from Triangles of a given color   *
 *                                                                            *
 ******************************************************************************
 */
int
extractColor(tess_T  *tess,             /* (in)  pointer to TESS */
             int     icolr,             /* (in)  color to select */
             tess_T  *subTess)          /* (out) pointer to sub-TESS object */
{
    int    status = 0;                  /* (out) return status */

    int    itri, ipnt, isid, *mapPnt=NULL;

    ROUTINE(extractColor);

    /* --------------------------------------------------------------- */

    if        (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (subTess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (icolr < 0 || icolr > 255) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* initailize the subTess */
    status = initialTess(subTess);
    CHECK_STATUS(initialTess);

    /* make and initialize array to keep track of new Points */
    MALLOC(mapPnt, int, tess->npnt);

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        mapPnt[ipnt] = -1;
    }

    /* loop through tess and add the Points associated with icolr */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_COLOR) == icolr) {
            for (isid = 0; isid < 3; isid++) {
                ipnt = tess->trip[3*itri+isid];

                if (mapPnt[ipnt] < 0) {
                    status = addPoint(subTess, tess->xyz[3*ipnt  ],
                                               tess->xyz[3*ipnt+1],
                                               tess->xyz[3*ipnt+2]);
                    CHECK_STATUS(addPoint);

                    mapPnt[ipnt] = status;
                }
            }
        }
    }

    /* create the necessary Triangles (without neighbor info) */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_COLOR) == icolr) {
            status = addTriangle(subTess, mapPnt[tess->trip[3*itri  ]],
                                          mapPnt[tess->trip[3*itri+1]],
                                          mapPnt[tess->trip[3*itri+2]], -1, -1, -1);
            CHECK_STATUS(addTriangle);
        }
    }

    /* set up the neighbor info */
    status = setupNeighbors(subTess);
    CHECK_STATUS(setupNeighbors);

cleanup:
    FREE(mapPnt);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * fillLoop - fill a loop with Triangles                                      *
 *                                                                            *
 ******************************************************************************
 */
int
fillLoop(tess_T  *tess,                 /* (in)  pointer to TESS */
         int     ipnt)                  /* (in)  index of one Point on loop (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int     nseg, imin, iseg, jseg, im1, ip1;
    double  amin, atst;
    seg_T   *seg = NULL;

    ROUTINE(fillLoop);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (ipnt < 0 || ipnt >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    }

    /* get the loop following the hanging sides from ipnt */
    status = getLoop(tess, ipnt, &nseg, &seg);
    CHECK_STATUS(getLoop);

    /* see if there are any repeated Points in the segment loop */
    for (iseg = 0; iseg < nseg; iseg++) {
        for (jseg = iseg+1; jseg < nseg; jseg++) {
            if (seg[iseg].pnt == seg[jseg].pnt) {
                printf("    segments iseg=%5d and jseg=%5d share Point %d\n",
                       iseg, jseg, seg[iseg].pnt);
            }
        }
    }

    /* as long as there are more than 3 Segments in the loop, cut off
          the Point with the smallest turn */
    while (nseg > 3) {

        /* find the Point with the smallest turn */
        imin = -1;
        amin = HUGEQ;

        for (iseg = 0; iseg < nseg; iseg++) {
            im1 = iseg - 1;
            ip1 = iseg + 1;
            if (im1 <  0   ) im1 += nseg;
            if (ip1 == nseg) ip1 -= nseg;

            atst = turn(tess, seg[im1].pnt, seg[iseg].pnt, seg[ip1].pnt, seg[iseg].tri);
            if (atst < amin) {
                imin = iseg;
                amin = atst;
            }
        }

        /* create a Triangle that cuts off the smallest angle Point */
        im1 = imin - 1;
        ip1 = imin + 1;
        if (im1 <  0   ) im1 += nseg;
        if (ip1 == nseg) ip1 -= nseg;

        status = addTriangle(tess, seg[ip1].pnt, seg[imin].pnt, seg[im1 ].pnt,
                                   seg[im1].tri, -1,            seg[imin].tri);
        CHECK_STATUS(addTriangle);

        /* remove the Point that was just cut off from the loop */
        seg[im1].tri = tess->ntri - 1;
        for (iseg = imin; iseg < nseg-1; iseg++) {
            seg[iseg] = seg[iseg+1];
        }
        nseg--;
    }

    /* make a Triangle with the final 3 Segments */
    status = addTriangle(tess, seg[0].pnt, seg[2].pnt, seg[1].pnt,
                               seg[1].tri, seg[0].tri, seg[2].tri);
    CHECK_STATUS(addTriangle);

cleanup:
    FREE(seg);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * findLoops - find the loops                                                 *
 *                                                                            *
 ******************************************************************************
 */
int
findLoops(tess_T  *tess,                /* (in)  pointer to TESS */
          int     *nloop,               /* (in)  maximum loops returned */
                                        /* (out) number of loops returned */
          int     ibeg[],               /* (out) array of a point in each loop */
          double  alen[])               /* (out) length (in XYZ)  of each loop */
{
    int    status = 0;                  /* (out) return status */

    int     maxloop, again, ipnt, jpnt, itri, ilup, jlup, isid, nseg, iseg, iswap;
    int     *stat=NULL;
    double  swap;
    seg_T   *seg=NULL;

    ROUTINE(findLoops);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (*nloop <= 0) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    maxloop = *nloop;
    *nloop  = 0;

    /* create an array that has -1 for all boundary Points and
       -2 for all interior Points */
    MALLOC(stat, int, tess->npnt);

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        stat[ipnt] = -2;
    }

    for (itri = 0; itri < tess->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            if (tess->trit[3*itri+isid] < 0) {
                stat[tess->trip[3*itri+(isid+1)%3]] = -1;
                stat[tess->trip[3*itri+(isid+2)%3]] = -1;
            }
        }
    }

    /* iteratively look for Points with stat[] == -1 */
    *nloop = 0;
    again  = 1;
    while (*nloop < maxloop-1 && again > 0) {
        again = 0;

        for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
            if (stat[ipnt] == -1) {
                ibeg[*nloop] = ipnt;
                alen[*nloop] = 0;

                status = getLoop(tess, ipnt, &nseg, &seg);
                CHECK_STATUS(getLoop);

                for (iseg = 0; iseg < nseg; iseg++) {
                    stat[seg[iseg].pnt] = *nloop;
                }

                FREE(seg);

                (*nloop)++;
                again++;
                break;
            }
        }
    }

    /* compute the lengths of each loop */
    for (itri = 0; itri < tess->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            if (tess->trit[3*itri+isid] < 0) {
                ipnt = tess->trip[3*itri+(isid+1)%3];
                jpnt = tess->trip[3*itri+(isid+2)%3];
                ilup = stat[ipnt];

                alen[ilup] += sqrt(SQR(tess->xyz[3*ipnt  ] - tess->xyz[3*jpnt  ])
                                  +SQR(tess->xyz[3*ipnt+1] - tess->xyz[3*jpnt+1])
                                  +SQR(tess->xyz[3*ipnt+2] - tess->xyz[3*jpnt+2]));
            }
        }
    }

    /* sort the loops from longest to shortest */
    for (ilup = 0; ilup < *nloop-1; ilup++) {
        for (jlup = ilup; jlup < *nloop; jlup++) {
            if (alen[jlup] > alen[ilup]) {
                iswap      = ibeg[ilup];
                ibeg[ilup] = ibeg[jlup];
                ibeg[jlup] = iswap;

                swap       = alen[ilup];
                alen[ilup] = alen[jlup];
                alen[jlup] = swap;
            }
        }
    }

cleanup:
    FREE(stat);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * flatten - flatten coordinates of a given color                             *
 *                                                                            *
 ******************************************************************************
 */
int
flattenColor(tess_T  *tess,                  /* (in)  pointer to TESS */
             int     icolr,                  /* (in)  color index */
             double  tol)                    /* (in)  tolerance */
{
    int         status = 0;                  /* (out)  return status */

    int         itri, isid, ipnt, navg;
    double      xmin, xmax, xavg, ymin, ymax, yavg, zmin, zmax, zavg;

    ROUTINE(flattenColor);

    /* ----------------------------------------------------------------------- */

    /* find the bounding box around the color */
    xmin = +HUGEQ;
    xmax = -HUGEQ;
    xavg = 0;
    ymin = +HUGEQ;
    ymax = -HUGEQ;
    yavg = 0;
    zmin = +HUGEQ;
    zmax = -HUGEQ;
    zavg = 0;
    navg = 0;

    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;
        if ((tess->ttyp[itri] & TRI_COLOR) != icolr) continue;

        for (isid = 0; isid < 3; isid++) {
            ipnt = tess->trip[3*itri+isid];

            if (tess->xyz[3*ipnt  ] < xmin) xmin = tess->xyz[3*ipnt  ];
            if (tess->xyz[3*ipnt  ] > xmax) xmax = tess->xyz[3*ipnt  ];
            if (tess->xyz[3*ipnt+1] < ymin) ymin = tess->xyz[3*ipnt+1];
            if (tess->xyz[3*ipnt+1] > ymax) ymax = tess->xyz[3*ipnt+1];
            if (tess->xyz[3*ipnt+2] < zmin) zmin = tess->xyz[3*ipnt+2];
            if (tess->xyz[3*ipnt+2] > zmax) zmax = tess->xyz[3*ipnt+2];

            xavg += tess->xyz[3*ipnt  ];
            yavg += tess->xyz[3*ipnt+1];
            zavg += tess->xyz[3*ipnt+2];
            navg += 1;
        }
    }
    xavg /= navg;
    yavg /= navg;
    zavg /= navg;

    printf("xmin=%12.5f  xmax=%12.5f  xavg=%12.5f\n", xmin, xmax, xavg);
    printf("ymin=%12.5f  ymax=%12.5f  yavg=%12.5f\n", ymin, ymax, yavg);
    printf("zmin=%12.5f  zmax=%12.5f  zavg=%12.5f\n", zmin, zmax, zavg);

    /* if the bounding box in one direction is much smaller than the others,
       flatten in that direction */
    if        ((xmax-xmin) < tol*(ymax-ymin) && (xmax-xmin) < tol*(zmax-zmin)) {
        printf("flattening to X=%12.5f\n", xavg);
        for (itri = 0; itri < tess->ntri; itri++) {
            if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;
            if ((tess->ttyp[itri] & TRI_COLOR) != icolr) continue;

            for (isid = 0; isid < 3; isid++) {
                ipnt = tess->trip[3*itri+isid];

                tess->xyz[3*ipnt  ] = xavg;
            }
        }
    } else if ((ymax-ymin) < tol*(zmax-zmin) && (ymax-ymin) < tol*(xmax-xmin)) {
        printf("flattening to Y=%12.5f\n", yavg);
        for (itri = 0; itri < tess->ntri; itri++) {
            if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;
            if ((tess->ttyp[itri] & TRI_COLOR) != icolr) continue;

            for (isid = 0; isid < 3; isid++) {
                ipnt = tess->trip[3*itri+isid];

                tess->xyz[3*ipnt+1] = yavg;
            }
        }
    } else if ((zmax-zmin) < tol*(xmax-xmin) && (zmax-zmin) < tol*(ymax-ymin)) {
        printf("flattening to Z=%12.5f\n", zavg);
        for (itri = 0; itri < tess->ntri; itri++) {
            if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;
            if ((tess->ttyp[itri] & TRI_COLOR) != icolr) continue;

            for (isid = 0; isid < 3; isid++) {
                ipnt = tess->trip[3*itri+isid];

                tess->xyz[3*ipnt+2] = zavg;
            }
        }
    } else {
        status = TESS_BAD_VALUE;
    }

//cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * floaterUV - apply floater to get UV at interior Points                     *
 *                                                                            *
 ******************************************************************************
 */
int
floaterUV(tess_T  *tess)                /* (in)  pointer to TESS */
{
    int         status = 0;                  /* (out)  return status */

    int         *onbound = NULL;
    double      *urhs    = NULL;
    double      *vrhs    = NULL;

    smf_T       smf;

    double      ang0, ang1, ang2, ang3, ang4, ang5;
    double      d01sq, d12sq, d20sq, d23sq;
    double      d30sq, d04sq, d41sq, d15sq, d52sq;
    double      umin, umax, vmin, vmax;
    double      sum, du, dv, erru, errv, err=1, errmax, errtol, omega;
    int         status2, isid, ient, iter, itmax;
    int         itri, jtri, ipnt, jpnt, ip0, ip1, ip2, ip3, ip4, ip5;

    ROUTINE(floaterUV);

    /* ----------------------------------------------------------------------- */

    smf.a    = NULL;
    smf.icol = NULL;
    smf.next = NULL;

    /* determine if Points are interior or on a boundary */
    MALLOC(onbound, int, tess->npnt);

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        onbound[ipnt] = 0;
    }

    for (itri = 0; itri < tess->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            if (tess->trit[3*itri+isid] < 0) {
                onbound[tess->trip[3*itri+(isid+1)%3]] = 1;
                onbound[tess->trip[3*itri+(isid+2)%3]] = 1;
            }
        }
    }

    /* find the extrema in u and v of the boundary Points */
    umin = +1e+10;
    umax = -1e+10;
    vmin = +1e+10;
    vmax = -1e+10;

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        if (onbound[ipnt] > 0) {
            if (tess->uv[2*ipnt  ] < umin) umin = tess->uv[2*ipnt  ];
            if (tess->uv[2*ipnt  ] > umax) umax = tess->uv[2*ipnt  ];
            if (tess->uv[2*ipnt+1] < vmin) vmin = tess->uv[2*ipnt+1];
            if (tess->uv[2*ipnt+1] > vmax) vmax = tess->uv[2*ipnt+1];
        }
    }

    /* initialize the sparse-matrix form */
    status = smfInit(&smf, tess->npnt);
    CHECK_STATUS(smfInit);

    /* find the mean value weights at each Point due to its neighbors.  the
         weight function is the "mean value coordinates" function as described
         by M. Floater

                        ip3--------ip2---------ip5
                          \         /\a5       /
                           \       /a2\       /
                            \     /    \     /
                             \a3 / itri \   /
                              \ /a0    a1\ /
                              ip0--------ip1
                                \     a4 /
                                 \      /
                                  \    /
                                   \  /
                                    \/
                                   ip4
     */
    for (itri = 0; itri < tess->ntri; itri++) {
        ip0 = tess->trip[3*itri  ];
        ip1 = tess->trip[3*itri+1];
        ip2 = tess->trip[3*itri+2];

        d01sq = SQR(tess->xyz[3*ip0  ] - tess->xyz[3*ip1  ])
              + SQR(tess->xyz[3*ip0+1] - tess->xyz[3*ip1+1])
              + SQR(tess->xyz[3*ip0+2] - tess->xyz[3*ip1+2]);
        d12sq = SQR(tess->xyz[3*ip1  ] - tess->xyz[3*ip2  ])
              + SQR(tess->xyz[3*ip1+1] - tess->xyz[3*ip2+1])
              + SQR(tess->xyz[3*ip1+2] - tess->xyz[3*ip2+2]);
        d20sq = SQR(tess->xyz[3*ip2  ] - tess->xyz[3*ip0  ])
              + SQR(tess->xyz[3*ip2+1] - tess->xyz[3*ip0+1])
              + SQR(tess->xyz[3*ip2+2] - tess->xyz[3*ip0+2]);

        ang0 = ACOS((d20sq + d01sq - d12sq) / 2 / sqrt(d20sq * d01sq));
        ang1 = ACOS((d01sq + d12sq - d20sq) / 2 / sqrt(d01sq * d12sq));
        ang2 = ACOS((d12sq + d20sq - d01sq) / 2 / sqrt(d12sq * d20sq));

        jtri = tess->trit[3*itri  ];
        if (jtri >= 0) {
            if (       tess->trit[3*jtri  ] == itri) {
                ip5 =  tess->trip[3*jtri  ];
            } else if (tess->trit[3*jtri+1] == itri) {
                ip5 =  tess->trip[3*jtri+1];
            } else {
                ip5 =  tess->trip[3*jtri+2];
            }

            d15sq = SQR(tess->xyz[3*ip1  ] - tess->xyz[3*ip5  ])
                  + SQR(tess->xyz[3*ip1+1] - tess->xyz[3*ip5+1])
                  + SQR(tess->xyz[3*ip1+2] - tess->xyz[3*ip5+2]);
            d52sq = SQR(tess->xyz[3*ip5  ] - tess->xyz[3*ip2  ])
                  + SQR(tess->xyz[3*ip5+1] - tess->xyz[3*ip2+1])
                  + SQR(tess->xyz[3*ip5+2] - tess->xyz[3*ip2+2]);

            ang5 = ACOS((d52sq + d12sq - d15sq) / 2 / sqrt(d52sq * d12sq));

            status = smfAdd(&smf, ip2, ip1);
            CHECK_STATUS(smfAdd);

            smf.a[status] = (tan(ang2/2) + tan(ang5/2)) / sqrt(d12sq);
        }

        jtri = tess->trit[3*itri+1];
        if (jtri >= 0) {
            if (       tess->trit[3*jtri  ] == itri) {
                ip3 =  tess->trip[3*jtri  ];
            } else if (tess->trit[3*jtri+1] == itri) {
                ip3 =  tess->trip[3*jtri+1];
            } else {
                ip3 =  tess->trip[3*jtri+2];
            }

            d23sq = SQR(tess->xyz[3*ip2  ] - tess->xyz[3*ip3  ])
                  + SQR(tess->xyz[3*ip2+1] - tess->xyz[3*ip3+1])
                  + SQR(tess->xyz[3*ip2+2] - tess->xyz[3*ip3+2]);
            d30sq = SQR(tess->xyz[3*ip3  ] - tess->xyz[3*ip0  ])
                  + SQR(tess->xyz[3*ip3+1] - tess->xyz[3*ip0+1])
                  + SQR(tess->xyz[3*ip3+2] - tess->xyz[3*ip0+2]);

            ang3 = ACOS((d30sq + d20sq - d23sq) / 2 / sqrt(d30sq * d20sq));

            status = smfAdd(&smf, ip0, ip2);
            CHECK_STATUS(smfAdd);

            smf.a[status] = (tan(ang0/2) + tan(ang3/2)) / sqrt(d20sq);
        }

        jtri = tess->trit[3*itri+2];
        if (jtri >= 0) {
            if (       tess->trit[3*jtri  ] == itri) {
                ip4 =  tess->trip[3*jtri  ];
            } else if (tess->trit[3*jtri+1] == itri) {
                ip4 =  tess->trip[3*jtri+1];
            } else {
                ip4 =  tess->trip[3*jtri+2];
            }

            d04sq = SQR(tess->xyz[3*ip0  ] - tess->xyz[3*ip4  ])
                  + SQR(tess->xyz[3*ip0+1] - tess->xyz[3*ip4+1])
                  + SQR(tess->xyz[3*ip0+2] - tess->xyz[3*ip4+2]);
            d41sq = SQR(tess->xyz[3*ip4  ] - tess->xyz[3*ip1  ])
                  + SQR(tess->xyz[3*ip4+1] - tess->xyz[3*ip1+1])
                  + SQR(tess->xyz[3*ip4+2] - tess->xyz[3*ip1+2]);

            ang4 = ACOS((d41sq + d01sq - d04sq) / 2 / sqrt(d41sq * d01sq));

            status = smfAdd(&smf, ip1, ip0);
            CHECK_STATUS(smfAdd);

            smf.a[status] = (tan(ang1/2) + tan(ang4/2)) / sqrt(d01sq);
        }
    }

    /* set up the final matrix and the right-hand sides */
    assert(tess->npnt > 0);                       // needed to avoid clang warning

    MALLOC(urhs, double, tess->npnt);
    MALLOC(vrhs, double, tess->npnt);

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {

        /* for interior Points, normalize the weights, set the diagonal to 1,
           and zero-out the right-hand sides */
        if (onbound[ipnt] == 0) {
            sum = 0;
            ient = ipnt;
            while (ient >= 0) {
                sum -= smf.a[ient];

                ient = smf.next[ient];
            }

            ient = ipnt;
            while (ient >= 0) {
                smf.a[ient] /= sum;

                ient = smf.next[ient];
            }

            smf.a[ipnt] = 1;

            urhs[ipnt] = 0;
            vrhs[ipnt] = 0;

        /* for boundary Points, zero-out all matrix elements, set the diagonal
           to 1, and store the boundary values in the RHS */
        } else {
            ient = ipnt;
            while (ient >= 0) {
                smf.a[ient] = 0;

                ient = smf.next[ient];
            }

            smf.a[ipnt] = 1;

            urhs[ipnt] = tess->uv[2*ipnt  ];
            vrhs[ipnt] = tess->uv[2*ipnt+1];
        }
    }

    /* solve for the Gs and Hs using successive-over-relaxation */
    itmax  = MAX(1000, tess->npnt);
    errtol = 1e-6 * MAX(umax-umin, vmax-vmin);
    omega  = 0.80;
    errmax = 0;
    for (iter = 0; iter < itmax; iter++) {

        /* apply successive-over-relaxation */
        for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
            du = urhs[ipnt];
            dv = vrhs[ipnt];

            ient = ipnt;
            while (ient >= 0) {
                jpnt = smf.icol[ient];

                du -= smf.a[ient] * tess->uv[2*jpnt  ];
                dv -= smf.a[ient] * tess->uv[2*jpnt+1];

                ient = smf.next[ient];
            }

            tess->uv[2*ipnt  ] += omega * du / smf.a[ipnt];
            tess->uv[2*ipnt+1] += omega * dv / smf.a[ipnt];
        }

        /* compute the norm of the residual */
        err = 0;
        for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
            erru = urhs[ipnt];
            errv = vrhs[ipnt];

            ient = ipnt;
            while (ient >= 0) {
                jpnt = smf.icol[ient];

                erru -= smf.a[ient] * tess->uv[2*jpnt  ];
                errv -= smf.a[ient] * tess->uv[2*jpnt+1];

                ient = smf.next[ient];
            }
            err += SQR(erru) + SQR(errv);
        }
        err = sqrt(err);

        if (err < errtol || iter%100 == 0) {
            printf(    "iter=%5d   omega=%8.3f   err=%11.4e\n", iter, omega, err);
        }

        /* exit if converged */
        if (err < errtol) {
            printf("    converged\n");
            goto cleanup;

        /* reset omega if we are diverging */
        } else if (err >= errmax) {
            omega  = 0.80;
            errmax = MAX(errmax, err);

        /* otherwise increase omega a little bit (but not more than 1.4) */
        } else {
            omega = MIN(1.001*omega, 1.40);
        }
    }

    /* we did not converge */
    status = TESS_NOT_CONVERGED;

 cleanup:
    status2 = smfFree(&smf);
    if (status2 < SUCCESS) {
        status = status2;
    }

    FREE(onbound);
    FREE(urhs);
    FREE(vrhs);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * freeTess - free a Tessellation                                             *
 *                                                                            *
 ******************************************************************************
 */
int
freeTess(tess_T  *tess)                 /* (in)  pointer to TESS to be freed */
{
    int    status = 0;                  /* (out) return status */

    ROUTINE(freeTess);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* remove all internal arrays */
    FREE(tess->trip);
    FREE(tess->trit);
    FREE(tess->ttyp);
    FREE(tess->bbox);
    FREE(tess->xyz );
    FREE(tess->uv  );
    FREE(tess->ptyp);

    /* initialize all the counters */
    tess->ntri  = 0;
    tess->npnt  = 0;
    tess->nlink = 0;
    tess->ncolr = 0;

    /* remove the octree if it exists */
    status = removeOctree(tess->octree);
    CHECK_STATUS(removeOctree);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * getLoop - get the Points along a loop                                      *
 *                                                                            *
 ******************************************************************************
 */
int
getLoop(tess_T  *tess,                  /* (in)  pointer to TESS */
        int     ipnt,                   /* (in)  index of one Point on loop (bias-0) */
        int     *nseg,                  /* (out) number of segments in loop */
        seg_T   *seg[])                 /* (out) arrau  of Segments (freeable) */
{
    int    status = 0;                  /* (out) return status */

    int    mseg, itri, jtri, jpnt, jseg, kseg;

    ROUTINE(getLoop);

    /* --------------------------------------------------------------- */

    /* default return */
    *seg = NULL;

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (ipnt < 0 || ipnt >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    }

    /* initial allocation for Segment table */
    *nseg = 0;
    mseg  = 1000;
    MALLOC(*seg, seg_T, mseg);

    /* find a Triangle that uses ipnt.  set up first Segment
       (and Point for second Segment) */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

        if        (tess->trip[3*itri  ] == ipnt && tess->trit[3*itri+2] < 0) {
            (*seg)[*nseg].pnt = ipnt;
            (*seg)[*nseg].tri = itri;
            (*nseg)++;

            (*seg)[*nseg].pnt = tess->trip[3*itri+1];
            break;
        } else if (tess->trip[3*itri+1] == ipnt && tess->trit[3*itri  ] < 0) {
            (*seg)[*nseg].pnt = ipnt;
            (*seg)[*nseg].tri = itri;
            (*nseg)++;

            (*seg)[*nseg].pnt = tess->trip[3*itri+2];
            break;
        } else if (tess->trip[3*itri+2] == ipnt && tess->trit[3*itri+1] < 0) {
            (*seg)[*nseg].pnt = ipnt;
            (*seg)[*nseg].tri = itri;
            (*nseg)++;

            (*seg)[*nseg].pnt = tess->trip[3*itri  ];
            break;
        }
    }

    if (*nseg == 0) {
        printf("ERROR:: could not find beginning Triangle\a\n");
        status = TESS_INTERNAL_ERROR;
        goto cleanup;
    }

    /* add Segments to the loop by spinning around "next" Point */
    while ((*seg)[*nseg].pnt != (*seg)[0].pnt) {
        ipnt = (*seg)[*nseg  ].pnt;
        itri = (*seg)[*nseg-1].tri;

        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) {
            printf("ERROR:: itri=%d is not active\a\n", itri);
            status = TESS_INTERNAL_ERROR;
            goto cleanup;
        }

        while (1) {
            if        (tess->trip[3*itri  ] == ipnt) {
                jpnt = tess->trip[3*itri+1];
                jtri = tess->trit[3*itri+2];
            } else if (tess->trip[3*itri+1] == ipnt) {
                jpnt = tess->trip[3*itri+2];
                jtri = tess->trit[3*itri  ];
            } else if (tess->trip[3*itri+2] == ipnt) {
                jpnt = tess->trip[3*itri  ];
                jtri = tess->trit[3*itri+1];
            } else {
                printf("ERROR:: Loop could not be constructed\a\n");
                status = TESS_INTERNAL_ERROR;
                goto cleanup;
            }

            if (jtri >= 0) {
                itri = jtri;
            } else {
                /* see if Point already exists in segments */
                jseg = -1;
                for (kseg = *nseg-1; kseg > 0; kseg--) {
                    if ((*seg)[kseg].pnt == ipnt) {
                        jseg = kseg;
                        break;
                    }
                }

                /* if it exists, cut out the part of the loop
                   back to jseg */
                if (jseg >= 0) {
                    (*seg)[jseg].tri = itri;
                    *nseg            = jseg + 1;

                /* otherwise add a new segment */
                } else {
                    (*seg)[*nseg].pnt = ipnt;
                    (*seg)[*nseg].tri = itri;
                    (*nseg)++;
                }

                if (*nseg > mseg-10) {
                    mseg += 1000;
                    RALLOC(*seg, seg_T, mseg);
                }

                (*seg)[*nseg].pnt = jpnt;
                break;
            }
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * initialTess - initialize a Tessellation                                    *
 *                                                                            *
 ******************************************************************************
 */
int
initialTess(tess_T  *tess)              /* (in)  pointer to TESS */
{
    int    status = 0;                  /* (out) return status */

    ROUTINE(initialTess);

    /* --------------------------------------------------------------- */

    /* initialize */
    tess->magic = TESS_MAGIC;
    tess->ntri  = 0;
    tess->mtri  = 1;
    tess->nhang = 0;          /* not computed within Tessellate.c */
    tess->nlink = 0;
    tess->ncolr = 0;
    tess->npnt  = 0;
    tess->mpnt  = 1;

    tess->trip   = NULL;
    tess->trit   = NULL;
    tess->ttyp   = NULL;
    tess->bbox   = NULL;
    tess->xyz    = NULL;
    tess->uv     = NULL;
    tess->ptyp   = NULL;
    tess->octree = NULL;

    MALLOC(tess->trip, int,    3*tess->mtri);
    MALLOC(tess->trit, int,    3*tess->mtri);
    MALLOC(tess->ttyp, int,      tess->mtri);
    MALLOC(tess->bbox, double, 6*tess->mtri);
    MALLOC(tess->xyz,  double, 3*tess->mpnt);
    MALLOC(tess->uv,   double, 2*tess->mpnt);
    MALLOC(tess->ptyp, int,      tess->mpnt);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * initialUV - initialize UV by projecting to best-fit plane of boundary Points *
 *                                                                            *
 ******************************************************************************
 */
int
initialUV(tess_T *tess)                 /* (in)  pointer to TESS */
{
    int    status = 0;                  /* (out) return status */

    int    nbound, itri, ipnt, isid;
    double xcent, ycent, zcent, Mat[9], evalue[3], evector[9];

    ROUTINE(initialUV);

    /* --------------------------------------------------------------- */

    /* find the centroid of the boundary Points */
    nbound = 0;
    xcent  = 0;
    ycent  = 0;
    zcent  = 0;

    for (itri = 0; itri < tess->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            if (tess->trit[3*itri+(isid  )%3] < 0 ||
                tess->trit[3*itri+(isid+1)%3] < 0   ) {
                ipnt = tess->trip[3*itri+(isid+2)%3];
                xcent += tess->xyz[3*ipnt  ];
                ycent += tess->xyz[3*ipnt+1];
                zcent += tess->xyz[3*ipnt+2];
                nbound++;
            }
        }
    }

    xcent /= nbound;
    ycent /= nbound;
    zcent /= nbound;

    /* create the sum-squares matrix of the boundary Points */
    Mat[0] = 0;          // Sxx
    Mat[1] = 0;          // Sxy
    Mat[2] = 0;          // Sxz
    Mat[4] = 0;          // Syy
    Mat[5] = 0;          // Syz
    Mat[8] = 0;          // Szz

    for (itri = 0; itri < tess->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            if (tess->trit[3*itri+(isid  )%3] < 0 ||
                tess->trit[3*itri+(isid+1)%3] < 0   ) {
                ipnt = tess->trip[3*itri+(isid+2)%3];

                Mat[0] += (tess->xyz[3*ipnt  ] - xcent) * (tess->xyz[3*ipnt  ] - xcent);
                Mat[1] += (tess->xyz[3*ipnt  ] - xcent) * (tess->xyz[3*ipnt+1] - ycent);
                Mat[2] += (tess->xyz[3*ipnt  ] - xcent) * (tess->xyz[3*ipnt+2] - zcent);
                Mat[4] += (tess->xyz[3*ipnt+1] - ycent) * (tess->xyz[3*ipnt+1] - ycent);
                Mat[5] += (tess->xyz[3*ipnt+1] - ycent) * (tess->xyz[3*ipnt+2] - zcent);
                Mat[8] += (tess->xyz[3*ipnt+2] - zcent) * (tess->xyz[3*ipnt+2] - zcent);
            }
        }
    }
    Mat[3] = Mat[1];     // Syx = Sxy
    Mat[6] = Mat[2];     // Szx = Sxz
    Mat[7] = Mat[5];     // Szy = Syz

    /* find the eigenvalues and eigenvectors of the sum-squares matrix */
    status = eigen(Mat, 3, evalue, evector);
    CHECK_STATUS(eigen);

    /* project all the Points to the plane defined by the first two eigenvectors */
    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        tess->uv[2*ipnt  ] = evector[0] * (tess->xyz[3*ipnt  ] - xcent)
                           + evector[3] * (tess->xyz[3*ipnt+1] - ycent)
                           + evector[6] * (tess->xyz[3*ipnt+2] - zcent);
        tess->uv[2*ipnt+1] = evector[1] * (tess->xyz[3*ipnt  ] - xcent)
                           + evector[4] * (tess->xyz[3*ipnt+1] - ycent)
                           + evector[7] * (tess->xyz[3*ipnt+2] - zcent);
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * joinPoints - join two points (at their average)                            *
 *                                                                            *
 ******************************************************************************
 */
int
joinPoints(tess_T  *tess,               /* (in)  pointer to TESS */
           int     ipnt,                /* (in)  index of first  Point (bias-0) */
           int     jpnt)                /* (in)  index of second Point (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    itri;

//    ROUTINE(joinPoints);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (ipnt < 0 || ipnt >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    } else if (jpnt < 0 || jpnt >= tess->npnt) {
        status = TESS_BAD_POINT_INDEX;
        goto cleanup;
    }

    /* join the points (at the average) */
    tess->xyz[3*ipnt  ] = (tess->xyz[3*ipnt  ] + tess->xyz[3*jpnt  ]) / 2;
    tess->xyz[3*ipnt+1] = (tess->xyz[3*ipnt+1] + tess->xyz[3*jpnt+1]) / 2;
    tess->xyz[3*ipnt+2] = (tess->xyz[3*ipnt+2] + tess->xyz[3*jpnt+2]) / 2;

    /* change the Point index from jpnt to ipnt */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

        if (tess->trip[3*itri  ] == jpnt) tess->trip[3*itri  ] = ipnt;
        if (tess->trip[3*itri+1] == jpnt) tess->trip[3*itri+1] = ipnt;
        if (tess->trip[3*itri+2] == jpnt) tess->trip[3*itri+2] = ipnt;
    }

    /* update the neighbors (to eliminate any degenerate loops that
       might be created */
    status = setupNeighbors(tess);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * makeLinks - make links between colors                                      *
 *                                                                            *
 ******************************************************************************
 */
int
makeLinks(tess_T  *tess)                /* (in)  pointer to TESS */
{
    int    status = 0;                  /* (out) return status */

    int    itri, jtri, isid;

    ROUTINE(makeLinks);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* make a link between any two Triangles with different colors
     note that we do not check for a link, since recreating one is
     not a problem */
    for (itri = 0; itri < tess->ntri; itri++) {
        for (isid = 0; isid < 3; isid++) {
            jtri = tess->trit[3*itri+isid];
            if ((tess->ttyp[itri] & TRI_COLOR) != (tess->ttyp[jtri] & TRI_COLOR)) {
                status = createLink(tess, itri, isid);
                CHECK_STATUS(createLink);
            }
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 *   nearApproach - find nearest approach between two line segments           *
 *                                                                            *
 ******************************************************************************
 */
int
nearApproach(double xyz_a[],            /* (in)  coordinates of first  point */
             double xyz_b[],            /* (in)  coordinates of second point */
             double xyz_c[],            /* (in)  coordinates of third  point */
             double xyz_d[],            /* (in)  coordinates of fourth point */
             double *s_ab,              /* (out) fraction of distance from a to b */
             double *s_cd,              /* (out) fraction of distance from c to d */
             double *dist)              /* (out) shortest distance */
{
    int       status = SUCCESS;         /* (out) return status */

    double    xba, yba, zba, xcd, ycd, zcd, xca, yca, zca;
    double    A, B, C, D, E, F;

    ROUTINE(nearApproach);

    /* --------------------------------------------------------------- */

    /* default output */
    *s_ab = -1;
    *s_cd = -1;
    *dist = -1;

    /* necessary differences */
    xba = xyz_b[0] - xyz_a[0];
    yba = xyz_b[1] - xyz_a[1];
    zba = xyz_b[2] - xyz_a[2];

    xcd = xyz_c[0] - xyz_d[0];
    ycd = xyz_c[1] - xyz_d[1];
    zcd = xyz_c[2] - xyz_d[2];

    xca = xyz_c[0] - xyz_a[0];
    yca = xyz_c[1] - xyz_a[1];
    zca = xyz_c[2] - xyz_a[2];

    /* set up matrix */
    A = xba * xba + yba * yba + zba * zba;
    B = xcd * xba + ycd * yba + zcd * zba;
    C = B;
    D = xcd * xcd + ycd * ycd + zcd * zcd;
    E = xba * xca + yba * yca + zba * zca;
    F = xcd * xca + ycd * yca + zcd * zca;

    if (fabs(A*D-B*C) > EPS20) {
        *s_ab = (E * D - B * F) / (A * D - B * C);
        *s_cd = (A * F - E * C) / (A * D - B * C);

        xba = (1-(*s_ab)) * xyz_a[0] + (*s_ab)*xyz_b[0];
        yba = (1-(*s_ab)) * xyz_a[1] + (*s_ab)*xyz_b[1];
        zba = (1-(*s_ab)) * xyz_a[2] + (*s_ab)*xyz_b[2];

        xcd = (1-(*s_cd)) * xyz_c[0] + (*s_cd)*xyz_d[0];
        ycd = (1-(*s_cd)) * xyz_c[1] + (*s_cd)*xyz_d[1];
        zcd = (1-(*s_cd)) * xyz_c[2] + (*s_cd)*xyz_d[2];

        *dist = sqrt( (xba - xcd) * (xba - xcd)
                    + (yba - ycd) * (yba - ycd)
                    + (zba - zcd) * (zba - zcd));
    }

//cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 *   nearestTo - find nearest point to Tessellation                           *
 *                                                                            *
 ******************************************************************************
 */
int
nearestTo(tess_T  *tess,                /* (in)  pointer to TESS */
          double  dbest,                /* (in)  initial best distance */
          double  xyz_in[],             /* (in)  input point */
          int    *ibest,                /* (out) Triangle index (bias-0) */
          double  xyz_out[])            /* (out) output point */
{
    int       status = SUCCESS;         /* (out) return status */

    int       itri, ip0, ip1, ip2;
    double    s0, s1, s01, dbest2, dtest2, xtest, ytest, ztest;
    double    x02, y02, z02, x12, y12, z12, xx2, yy2, zz2, A, B, C, D, E, F;
    oct_T     *octree;

    ROUTINE(nearestTo);

    /* --------------------------------------------------------------- */

    /* default output */
    *ibest     = -1;
    xyz_out[0] = xyz_in[0];
    xyz_out[1] = xyz_in[1];
    xyz_out[2] = xyz_in[2];

    dbest2 = dbest * dbest;

    /* if an octree exists, only use the Triangles in this octant */
    if (tess->octree != NULL) {
        octree = tess->octree;

        while (octree->child != NULL) {
            if (xyz_in[2] < octree->zcent) {
                if (xyz_in[1] < octree->ycent) {
                    if (xyz_in[0] < octree->xcent) {
                        octree = &(octree->child[0]);
                    } else {
                        octree = &(octree->child[1]);
                    }
                } else {
                    if (xyz_in[0] < octree->xcent) {
                        octree = &(octree->child[2]);
                    } else {
                        octree = &(octree->child[3]);
                    }
                }
            } else {
                if (xyz_in[1] < octree->ycent) {
                    if (xyz_in[0] < octree->xcent) {
                        octree = &(octree->child[4]);
                    } else {
                        octree = &(octree->child[5]);
                    }
                } else {
                    if (xyz_in[0] < octree->xcent) {
                        octree = &(octree->child[6]);
                    } else {
                        octree = &(octree->child[7]);
                    }
                }
            }
        }

        /* found the octant, so loop through its Triangles */
        for (itri = 0; itri < octree->ntri; itri++) {

            /* determine barycentric coordinates of point closest to xyz_in */
            ip0 = octree->trip[3*itri  ];
            ip1 = octree->trip[3*itri+1];
            ip2 = octree->trip[3*itri+2];

            x02 = octree->xyz[3*ip0  ] - octree->xyz[3*ip2  ];
            y02 = octree->xyz[3*ip0+1] - octree->xyz[3*ip2+1];
            z02 = octree->xyz[3*ip0+2] - octree->xyz[3*ip2+2];
            x12 = octree->xyz[3*ip1  ] - octree->xyz[3*ip2  ];
            y12 = octree->xyz[3*ip1+1] - octree->xyz[3*ip2+1];
            z12 = octree->xyz[3*ip1+2] - octree->xyz[3*ip2+2];
            xx2 = xyz_in[0]            - octree->xyz[3*ip2  ];
            yy2 = xyz_in[1]            - octree->xyz[3*ip2+1];
            zz2 = xyz_in[2]            - octree->xyz[3*ip2+2];

            A = x02 * x02 + y02 * y02 + z02 * z02;
            B = x12 * x02 + y12 * y02 + z12 * z02;
            C = B;
            D = x12 * x12 + y12 * y12 + z12 * z12;
            E = xx2 * x02 + yy2 * y02 + zz2 * z02;
            F = xx2 * x12 + yy2 * y12 + zz2 * z12;

            if (fabs(A*D-B*C) < EPS20) continue;

            s0 = (E * D - B * F) / (A * D - B * C);
            s1 = (A * F - E * C) / (A * D - B * C);

            /* clip the barycentric coordinates and evaluate */
            s0 = MINMAX(0, s0, 1);
            s1 = MINMAX(0, s1, 1);

            s01 = s0 + s1;
            if (s01 > 1) {
                s0 /= s01;
                s1 /= s01;
            }

            xtest = octree->xyz[3*ip2  ] + s0 * x02 + s1 * x12;
            ytest = octree->xyz[3*ip2+1] + s0 * y02 + s1 * y12;
            ztest = octree->xyz[3*ip2+2] + s0 * z02 + s1 * z12;

            /* if point is closer than anything we have seen so far, remember this */
            dtest2 = (xtest - xyz_in[0]) * (xtest - xyz_in[0])
                   + (ytest - xyz_in[1]) * (ytest - xyz_in[1])
                   + (ztest - xyz_in[2]) * (ztest - xyz_in[2]);

            if (dtest2 < dbest2) {
                *ibest     = itri;
                xyz_out[0] = xtest;
                xyz_out[1] = ytest;
                xyz_out[2] = ztest;
                dbest2     = dtest2;
            }
        }

    /* no octree exists, so use all the Triangles */
    } else {

        /* loop through all Triangles */
        for (itri = 0; itri < tess->ntri; itri++) {
            if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

            /* determine barycentric coordinates of point closest to xyz_in */
            ip0 = tess->trip[3*itri  ];
            ip1 = tess->trip[3*itri+1];
            ip2 = tess->trip[3*itri+2];

            x02 = tess->xyz[3*ip0  ] - tess->xyz[3*ip2  ];
            y02 = tess->xyz[3*ip0+1] - tess->xyz[3*ip2+1];
            z02 = tess->xyz[3*ip0+2] - tess->xyz[3*ip2+2];
            x12 = tess->xyz[3*ip1  ] - tess->xyz[3*ip2  ];
            y12 = tess->xyz[3*ip1+1] - tess->xyz[3*ip2+1];
            z12 = tess->xyz[3*ip1+2] - tess->xyz[3*ip2+2];
            xx2 = xyz_in[0]          - tess->xyz[3*ip2  ];
            yy2 = xyz_in[1]          - tess->xyz[3*ip2+1];
            zz2 = xyz_in[2]          - tess->xyz[3*ip2+2];

            A = x02 * x02 + y02 * y02 + z02 * z02;
            B = x12 * x02 + y12 * y02 + z12 * z02;
            C = B;
            D = x12 * x12 + y12 * y12 + z12 * z12;
            E = xx2 * x02 + yy2 * y02 + zz2 * z02;
            F = xx2 * x12 + yy2 * y12 + zz2 * z12;

            if (fabs(A*D-B*C) < EPS20) continue;

            s0 = (E * D - B * F) / (A * D - B * C);
            s1 = (A * F - E * C) / (A * D - B * C);

            /* clip the barycentric coordinates and evaluate */
            s0 = MINMAX(0, s0, 1);
            s1 = MINMAX(0, s1, 1);

            s01 = s0 + s1;
            if (s01 > 1) {
                s0 /= s01;
                s1 /= s01;
            }

            xtest = tess->xyz[3*ip2  ] + s0 * x02 + s1 * x12;
            ytest = tess->xyz[3*ip2+1] + s0 * y02 + s1 * y12;
            ztest = tess->xyz[3*ip2+2] + s0 * z02 + s1 * z12;

            /* if point is closer than anything we have seen so far, remember this */
            dtest2 = (xtest - xyz_in[0]) * (xtest - xyz_in[0])
                   + (ytest - xyz_in[1]) * (ytest - xyz_in[1])
                   + (ztest - xyz_in[2]) * (ztest - xyz_in[2]);

            if (dtest2 < dbest2) {
                *ibest     = itri;
                xyz_out[0] = xtest;
                xyz_out[1] = ytest;
                xyz_out[2] = ztest;
                dbest2     = dtest2;
            }
        }
    }

//cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * readStlAscii - read an ASCII stl file                                      *
 *                                                                            *
 ******************************************************************************
 */
int
readStlAscii(tess_T  *tess,             /* (in)  pointer to TESS */
             char    *filename)         /* (in)  name of file */
{
    int    status = 0;                  /* (out) return status */

    int    itri, isid, ipnt;
    LONG   key1, key2, key3;
    double xin, yin, zin;
    char   string[255];
    rbt_T  *ntree = NULL;
    FILE   *fp = NULL;

    ROUTINE(readStlAscii);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* initialize the TESS */
    status = initialTess(tess);
    CHECK_STATUS(initialTess);

    /* get a red-black tree in which the Points will be stored */
    status = rbtCreate(1000, &ntree);
    if (status < SUCCESS || ntree == NULL) {
        printf("ERROR:: ntree could not be allocated in routine readStlAscii\a\n");
        exit(0);
    }

    /* count the number of triangles by reading file and looking for "facet" */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        status = TESS_BAD_FILE_NAME;
        goto cleanup;
    }

    while (feof(fp) != 1) {
        (void) fgets(string, 255, fp);
        if (strncmp(string,    "facet", 5) == 0 ||
            strncmp(string,  "  facet", 7) == 0 ||
            strncmp(string, "   facet", 8) == 0   ) {
            (tess->ntri)++;
        }
    }
    fclose(fp);

    /* make room for the Triangles */
    tess->mtri = tess->ntri;

    assert(tess->mtri > 0);                       // needed to avoid clang warning

    RALLOC(tess->trip, int,    3*tess->mtri);
    RALLOC(tess->trit, int,    3*tess->mtri);
    RALLOC(tess->ttyp, int,      tess->mtri);
    RALLOC(tess->bbox, double, 6*tess->mtri);

    /* open the file and read its header */
    fp = fopen(filename, "r");

    (void) fgets(string, 255, fp);                  /* solid */
    (void) fgets(string, 255, fp);                  /* facet -or- endsolid */

    /* read the data associated with a facet (Triangle) */
    itri = 0;
    while (strncmp(string,    "facet", 5) == 0 ||
           strncmp(string,  "  facet", 7) == 0 ||
           strncmp(string, "   facet", 8) == 0   ) {
        (void) fgets(string, 255, fp);              /* outer */

        /* read each of the Triangle's Points */
        for (isid = 0; isid < 3; isid++) {
            fscanf(fp, "%s %lf %lf %lf", string, &xin, &yin, &zin);

            /* see if the point already exists */
            key1 = (LONG)(xin * 10000000);
            key2 = (LONG)(yin * 10000000);
            key3 = (LONG)(zin * 10000000);
            ipnt = rbtSearch(ntree, key1, key2, key3);

            /* create a new Point if needed */
            if (ipnt < 0) {
                status = addPoint(tess, xin, yin, zin);
                CHECK_STATUS(addPoint);

                ipnt = rbtInsert(ntree, key1, key2, key3);
                if (ipnt != (tess->npnt-1)) {
                    printf("ERROR:: Trouble with inserting in tree, ipnt=%d, npnt=%d\a\n", ipnt, tess->npnt);
                    status = TESS_INTERNAL_ERROR;
                    goto cleanup;
                }
            }

            /* remember the Point's id */
            if        (isid == 0) {
                tess->trip[3*itri  ] = ipnt;
            } else if (isid == 1) {
                tess->trip[3*itri+1] = ipnt;
            } else {
                tess->trip[3*itri+2] = ipnt;
            }
        }

        /* create the Triangle and get ready for next "read" */
        tess->trit[3*itri  ] = -1;
        tess->trit[3*itri+1] = -1;
        tess->trit[3*itri+2] = -1;
        tess->ttyp[  itri  ] = TRI_ACTIVE | TRI_VISIBLE;

        (void) fgets(string, 255, fp);              /* throw away -eol- */
        (void) fgets(string, 255, fp);              /* endloop */
        (void) fgets(string, 255, fp);              /* endfacet */
        (void) fgets(string, 255, fp);              /* facet -or- endsolid */

        if (itri % 10000 == 0) {
            printf(".");
            fflush(stdout);
        }

        itri++;
    }
    printf("\n");
    fclose(fp);

    rbtDelete(ntree);

    printf("    After reading: npnt = %8d\n", tess->npnt);
    printf("                   ntri = %8d\n", tess->ntri);

    /* set up Triangle neighbors */
    status = setupNeighbors(tess);
    CHECK_STATUS(setupNeighbors);

cleanup:
    FREE(ntree);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * readStlBinary - read a binary stl file                                     *
 *                                                                            *
 ******************************************************************************
 */
int
readStlBinary(tess_T  *tess,            /* (in)  pointer to TESS */
              char    *filename)        /* (in)  name of file */
{
    int    status = 0;                  /* (out) return status */

    int    isid, ipnt, itri;
    UINT16 nattr;
    UINT32 ntri32;
    LONG   key1, key2, key3;
    double xin, yin, zin;
    REAL32 normal[3], vertex[3];
    char   header[80];
    rbt_T  *ntree = NULL;
    FILE   *fp = NULL;

    ROUTINE(readStlBinary);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* initialize the TESS */
    status = initialTess(tess);
    CHECK_STATUS(initialTess);

    /* get a red-black tree in which the Points will be stored */
    status = rbtCreate(1000, &ntree);
    if (status < SUCCESS || ntree == NULL) {
        printf("ERROR:: ntree could not be allocated in routine readStlBinary\a\n");
        exit(0);
    }

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        status = TESS_BAD_FILE_NAME;
        goto cleanup;
    }

    /* read the header */
    (void) fread(header, sizeof(char), 80, fp);

    /* get the number of triangles */
    (void) fread(&ntri32, sizeof(UINT32), 1, fp);
    tess->ntri = ntri32;

    /* make room for the Triangles */
    tess->mtri = tess->ntri;

    RALLOC(tess->trip, int,    3*tess->mtri);
    RALLOC(tess->trit, int,    3*tess->mtri);
    RALLOC(tess->ttyp, int,      tess->mtri);
    RALLOC(tess->bbox, double, 6*tess->mtri);

    /* read the Triangles */
    for (itri = 0; itri < tess->ntri; itri++) {
        (void) fread(normal, sizeof(REAL32), 3, fp);

        for (isid = 0; isid < 3; isid++) {
            (void) fread(vertex, sizeof(REAL32), 3, fp);

            xin = vertex[0];
            yin = vertex[1];
            zin = vertex[2];

            /* see if the point already exists */
            key1 = (LONG)(xin * 10000000);
            key2 = (LONG)(yin * 10000000);
            key3 = (LONG)(zin * 10000000);
            ipnt = rbtSearch(ntree, key1, key2, key3);

            /* create a new Point if needed */
            if (ipnt < 0) {
                status = addPoint(tess, xin, yin, zin);
                CHECK_STATUS(addPoint);

                ipnt = rbtInsert(ntree, key1, key2, key3);
                if (ipnt != (tess->npnt-1)) {
                    printf("ERROR:: Trouble with inserting in tree, ipnt=%d, npnt=%d\a\n", ipnt, tess->npnt);
                    status = TESS_INTERNAL_ERROR;
                    goto cleanup;
                }
            }

            /* remember the Point's id */
            if        (isid == 0) {
                tess->trip[3*itri  ] = ipnt;
            } else if (isid == 1) {
                tess->trip[3*itri+1] = ipnt;
            } else {
                tess->trip[3*itri+2] = ipnt;
            }
        }

        (void) fread(&nattr, sizeof(UINT16), 1, fp);

        /* create the Triangle and get ready for next "read" */
        tess->trit[3*itri  ] = -1;
        tess->trit[3*itri+1] = -1;
        tess->trit[3*itri+2] = -1;
        tess->ttyp[  itri  ] = TRI_ACTIVE | TRI_VISIBLE | nattr;

        tess->ncolr = MAX(tess->ncolr, nattr);

        if (itri % 10000 == 0) {
            printf(".");
            fflush(stdout);
        }
    }
    printf("\n");
    fclose(fp);

    rbtDelete(ntree);

    printf("    After reading: npnt = %8d\n", tess->npnt );
    printf("                   ntri = %8d\n", tess->ntri );
    printf("                  ncolr = %8d\n", tess->ncolr);

    /* set up Triangle neighbors */
    status = setupNeighbors(tess);
    CHECK_STATUS(setupNeighbors);

cleanup:
    FREE(ntree);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * readTriAscii - read an ASCII tri file                                      *
 *                                                                            *
 ******************************************************************************
 */
int
readTriAscii(tess_T  *tess,             /* (in)  pointer to TESS */
             char    *filename)         /* (in)  name of file */
{
    int    status = 0;                  /* (out) return status */

    int    jtri, ntri, jpnt, npnt, ibody, jbody;
    int    jface, ip0, ip1, ip2, it0, it1, it2;
    double uin, vin, xin, yin, zin;
    rbt_T  *ntree = NULL;
    FILE   *fp = NULL;

    ROUTINE(readTriAscii);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* initialize the TESS */
    status = initialTess(tess);
    CHECK_STATUS(initialTess);

    /* get a red-black tree in which the Points will be stored */
    status = rbtCreate(1000, &ntree);
    if (ntree == NULL) {
        printf("ERROR:: ntree could not be allocated in routine readStlAscii\a\n");
        exit(0);
    }

    printf("Enter ibody: "); scanf("%d", &ibody);

    /* read the file until the requested Body is found */
    fp = fopen(filename, "r");

    while (feof(fp) != 1) {
        fscanf(fp, "%d %d %d", &jbody, &npnt, &ntri);

        /* skip if ibody does not match jbody */
        if (jbody != ibody) {
            printf("skipping Body %d\n", jbody);

            for (jpnt = 0; jpnt < npnt; jpnt++) {
                fscanf(fp, "%lf %lf %lf %lf %lf",
                       &xin, &yin, &zin, &uin, &vin);
            }
            for (jtri = 0; jtri < ntri; jtri++) {
                fscanf(fp, "%d %d %d %d %d %d %d",
                       &jface, &ip0, &ip1, &ip2, &it0, &it1, &it2);
            }

        /* jbody matches ibody */
        } else {
            printf("reading  Body %d (npnt=%d, ntri=%d)\n", jbody, npnt, ntri);

            /* make room for the Triangles */
            tess->ntri = ntri;
            tess->mtri = tess->ntri;

            RALLOC(tess->trip, int,    3*tess->mtri);
            RALLOC(tess->trit, int,    3*tess->mtri);
            RALLOC(tess->ttyp, int,      tess->mtri);
            RALLOC(tess->bbox, double, 6*tess->mtri);

            /* read the Points */
            printf("Points "); fflush(stdout);
            for (jpnt = 0; jpnt < npnt; jpnt++) {
                if (jpnt%100000 == 0) {printf("."); fflush(stdout);}

                fscanf(fp, "%lf %lf %lf %lf %lf",
                       &xin, &yin, &zin, &uin, &vin);

                /* create a new Point if needed */
                status = addPoint(tess, xin, yin, zin);
                CHECK_STATUS(addPoint);
            }
            printf(" done\n");

            /* read the Triangles */
            printf("Triangles "); fflush(stdout);
            for (jtri = 0; jtri < ntri; jtri++) {
                if (jtri%100000 == 0) {printf("."); fflush(stdout);}

                fscanf(fp, "%d %d %d %d %d %d %d",
                       &jface, &ip0, &ip1, &ip2, &it0, &it1, &it2);

                /* create the Triangle and get ready for next "read" */
                tess->trip[3*jtri  ] = ip0;
                tess->trip[3*jtri+1] = ip1;
                tess->trip[3*jtri+2] = ip2;
                tess->trit[3*jtri  ] = it0;
                tess->trit[3*jtri+1] = it1;
                tess->trit[3*jtri+2] = it2;
                tess->ttyp[  jtri  ] = TRI_ACTIVE | TRI_VISIBLE | jface;

                tess->ncolr = MAX(tess->ncolr, jface);
            }
            printf(" done\n");

            fclose(fp);
            goto cleanup;
        }
    }

    rbtDelete(ntree);

    printf("    After reading: npnt = %8d\n", tess->npnt);
    printf("                   ntri = %8d\n", tess->ntri);

cleanup:
    FREE(ntree);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * removeLinks - remove all Links                                             *
 *                                                                            *
 ******************************************************************************
 */
int
removeLinks(tess_T  *tess)              /* (in)  pointer to TESS */
{
    int    status = 0;                  /* (out) return status */

    int    itri;

//    ROUTINE(removeLinks);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    for (itri = 0; itri < tess->ntri; itri++) {
        tess->ttyp[itri] &= ~TRI_LINK;
    }

    tess->nlink = 0;

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * scribe - scribe between given Points                                       *
 *                                                                            *
 ******************************************************************************
 */
int
scribe(tess_T  *tess,                   /* (in)  pointer to TESS */
       int     isrc,                    /* (in)  index of source Point (bias-0) */
       int     itgt)                    /* (in)  index of target Point (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int    ipath, ipnt, itri, *prev=NULL, *link=NULL, ismth, nsmth=1001;
    int    ip0, ip1, ip2, npath, ibest, imax=0, jmax=-1;
    double alen, *path=NULL, dtest, dbest, sbest, dxyztol, dxyzmax;
    double xyz_in[3], xyz_out[3], xyz_a[3], xyz_b[3], xyz_c[3], xyz_d[3], s_ab, s_cd;

    ROUTINE(scribe);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* create an octree to make inverse evaluations faster */
    MALLOC(tess->octree, oct_T, 1);

    status = buildOctree(tess, 1000, tess->octree);
    CHECK_STATUS(buildOctree);

    /* allocate storage for Dijkstra's algorithm */
    MALLOC(prev, int, tess->npnt);
    MALLOC(link, int, tess->npnt);

    /* find the path via dijkstra (done backwards to make extraction
       of the path easier) */
    status = dijkstra(tess, itgt, isrc, prev, link);
    CHECK_STATUS(dijkstra);

    /* determine the number of Points in the path */
    npath = 2;
    ipnt  = isrc;
    while (prev[ipnt] != itgt) {
        ipnt = prev[ipnt];
        npath++;
    }

    /* create the path */
    MALLOC(path, double, 3*npath);

    npath = 0;
    ipnt  = isrc;
    path[3*npath  ] = tess->xyz[3*ipnt  ];
    path[3*npath+1] = tess->xyz[3*ipnt+1];
    path[3*npath+2] = tess->xyz[3*ipnt +2];
    npath++;
    while (prev[ipnt] != itgt) {
        ipnt = prev[ipnt];
        path[3*npath  ] = tess->xyz[3*ipnt  ];
        path[3*npath+1] = tess->xyz[3*ipnt+1];
        path[3*npath+2] = tess->xyz[3*ipnt+2];
        npath++;
    }
    ipnt  = itgt;
    path[3*npath  ] = tess->xyz[3*ipnt  ];
    path[3*npath+1] = tess->xyz[3*ipnt+1];
    path[3*npath+2] = tess->xyz[3*ipnt +2];
    npath++;

    /* find the average length of segment in the path */
    alen = 0;
    for (ipath = 1; ipath < npath; ipath++) {
        alen += sqrt( (path[3*ipath  ] - path[3*ipath-3]) * (path[3*ipath  ] - path[3*ipath-3])
                    + (path[3*ipath+1] - path[3*ipath-2]) * (path[3*ipath+1] - path[3*ipath-2])
                    + (path[3*ipath+2] - path[3*ipath-1]) * (path[3*ipath+2] - path[3*ipath-1]));
    }
    alen /= (npath-1);

    dxyztol = EPS06 * alen;
    nsmth   = 10000;

    /* smooth the path while making sure it stays on surface */
    for (ismth = 0; ismth < nsmth; ismth++) {

        dxyzmax = 0;
        for (ipath = 1; ipath < npath-1; ipath++) {
            xyz_in[0] = (path[3*ipath-3] + 2*path[3*ipath  ] + path[3*ipath+3]) / 4;
            xyz_in[1] = (path[3*ipath-2] + 2*path[3*ipath+1] + path[3*ipath+4]) / 4;
            xyz_in[2] = (path[3*ipath-1] + 2*path[3*ipath+2] + path[3*ipath+5]) / 4;

            status = nearestTo(tess, alen, xyz_in, &itri, xyz_out);
            if (status < 0 || itri < 0) {printf("problem\n"); exit(0);}

            if (fabs(path[3*ipath  ]-xyz_out[0]) > dxyzmax) {
                imax = ipath;
                jmax = 0;
                dxyzmax = fabs(path[3*ipath  ]-xyz_out[0]);
            }
            if (fabs(path[3*ipath+1]-xyz_out[1]) > dxyzmax) {
                imax = ipath;
                jmax = 1;
                dxyzmax = fabs(path[3*ipath+1]-xyz_out[1]);
            }
            if (fabs(path[3*ipath+2]-xyz_out[2]) > dxyzmax) {
                imax = ipath;
                jmax = 2;
                dxyzmax = fabs(path[3*ipath+2]-xyz_out[2]);
            }

            path[3*ipath  ] = xyz_out[0];
            path[3*ipath+1] = xyz_out[1];
            path[3*ipath+2] = xyz_out[2];
        }

        if (ismth%100 == 0) {
            printf("%5d  dxyzmax=%12.5e (%3d,%d) %10.4f %10.4f %10.4f\n", ismth, dxyzmax, imax, jmax,
                   path[3*imax], path[3*imax+1], path[3*imax+2]);
        }

        /* check (and print) convergence */
        if (dxyzmax < dxyztol) {
            printf("converged\n");
            printf("%5d  dxyzmax=%12.5e (%3d,%d) %10.4f %10.4f %10.4f\n", ismth, dxyzmax, imax, jmax,
                   path[3*imax], path[3*imax+1], path[3*imax+2]);
            break;
        }
    }

    /* start at src */
    ipnt = isrc;

    for (ipath = 1; ipath < npath; ipath++) {
        xyz_a[0] = path[3*ipath-3];
        xyz_a[1] = path[3*ipath-2];
        xyz_a[2] = path[3*ipath-1];
        xyz_b[0] = path[3*ipath  ];
        xyz_b[1] = path[3*ipath+1];
        xyz_b[2] = path[3*ipath+2];

        /* move to the next point if a and b are very close to each other */
        if (fabs(xyz_a[0]-xyz_b[0]) < EPS06 &&
            fabs(xyz_a[1]-xyz_b[1]) < EPS06 &&
            fabs(xyz_a[2]-xyz_b[2]) < EPS06   ) continue;

        /* find a Triangle that contains ipnt, which intersects
           the segment from a to b, and which has the minimum
           distance */
        ibest = -1;
        dbest = alen;
        sbest = -1;

        for (itri = 0; itri < tess->ntri; itri++) {
            if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

            ip0 = tess->trip[3*itri  ];
            ip1 = tess->trip[3*itri+1];
            ip2 = tess->trip[3*itri+2];

            if        (ip0 == ipnt) {
                xyz_c[0] = tess->xyz[3*ip1  ];
                xyz_c[1] = tess->xyz[3*ip1+1];
                xyz_c[2] = tess->xyz[3*ip1+2];

                xyz_d[0] = tess->xyz[3*ip2  ];
                xyz_d[1] = tess->xyz[3*ip2+1];
                xyz_d[2] = tess->xyz[3*ip2+2];

                status =  nearApproach(xyz_a, xyz_b, xyz_c, xyz_d, &s_ab, &s_cd, &dtest);
                CHECK_STATUS(nearApproach);

                if (s_ab > 0 && s_ab <= 1+EPS06 && s_cd >= 0 && s_cd <= 1 && dtest < dbest) {
                    dbest = dtest;
                    sbest = s_cd;
                    ibest = itri;
                }
            } else if (ip1 == ipnt) {
                xyz_c[0] = tess->xyz[3*ip2  ];
                xyz_c[1] = tess->xyz[3*ip2+1];
                xyz_c[2] = tess->xyz[3*ip2+2];

                xyz_d[0] = tess->xyz[3*ip0  ];
                xyz_d[1] = tess->xyz[3*ip0+1];
                xyz_d[2] = tess->xyz[3*ip0+2];

                status =  nearApproach(xyz_a, xyz_b, xyz_c, xyz_d, &s_ab, &s_cd, &dtest);
                CHECK_STATUS(nearApproach);

                if (s_ab > 0 && s_ab <= 1+EPS06 && s_cd >= 0 && s_cd <= 1 && dtest < dbest) {
                    dbest = dtest;
                    sbest = s_cd;
                    ibest = itri;
                }
            } else if (ip2 == ipnt) {
                xyz_c[0] = tess->xyz[3*ip0  ];
                xyz_c[1] = tess->xyz[3*ip0+1];
                xyz_c[2] = tess->xyz[3*ip0+2];

                xyz_d[0] = tess->xyz[3*ip1  ];
                xyz_d[1] = tess->xyz[3*ip1+1];
                xyz_d[2] = tess->xyz[3*ip1+2];

                status =  nearApproach(xyz_a, xyz_b, xyz_c, xyz_d, &s_ab, &s_cd, &dtest);
                CHECK_STATUS(nearApproach);

                if (s_ab > 0 && s_ab <= 1+EPS06 && s_cd >= 0 && s_cd <= 1 && dtest < dbest) {
                    dbest = dtest;
                    sbest = s_cd;
                    ibest = itri;
                }
            }
        }

        /* if we found a good Triangle, split it, advance ipnt,
           and decrement ipath so that this segment of path
           gets used again */
        if (ibest >= 0) {

            /* if we are at the end of the path, make sure that
               we use target point */
            if (ipath == npath+1) {
                



            /* otherwise, if sbest is either very small or big,
               just use the closest point in ibest */
            } else if (sbest < 0.0002) {
                if        (ipnt == tess->trip[3*ibest  ]) {
                    ipnt = tess->trip[3*ibest+1];
                } else if (ipnt == tess->trip[3*ibest+1]) {
                    ipnt = tess->trip[3*ibest+2];
                } else {
                    ipnt = tess->trip[3*ibest  ];
                }
            } else if (sbest > 0.9998) {
                if        (ipnt == tess->trip[3*ibest  ]) {
                    ipnt = tess->trip[3*ibest+2];
                } else if (ipnt == tess->trip[3*ibest+1]) {
                    ipnt = tess->trip[3*ibest  ];
                } else {
                    ipnt = tess->trip[3*ibest+1];
                }
            } else {
                status = splitTriangle(tess, ibest, ipnt, sbest);
                CHECK_STATUS(splitTriangle);

                ipnt = tess->npnt - 1;
            }

            /* update path so that we do not use (again) the part
               we have already used */
            ipath--;
            path[3*ipath  ] = tess->xyz[3*ipnt  ];
            path[3*ipath+1] = tess->xyz[3*ipnt+1];
            path[3*ipath+2] = tess->xyz[3*ipnt+2];

        /* otherwise, try the next segment of the path */
        } else {
            continue;
        }
    }

    /* remove the octree, since future operations might change
       the tessellation (and hence make the octree incorrect */
    status = removeOctree(tess->octree);
    CHECK_STATUS(removeOctree);

    FREE(tess->octree);

cleanup:
    FREE(path);
    FREE(link);
    FREE(prev);

    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * setupNeighbors - set up Triangle neighbor information (generally used internally) *
 *                                                                            *
 ******************************************************************************
 */
int
setupNeighbors(tess_T  *tess)           /* (in)  pointer to TESS */
{
    int    status = 0;                  /* (out) return status */

    typedef struct {
        int         ltri;
        int         lsid;
        int         rtri;
        int         rsid;
        int         bpnt;
        int         epnt;
    } Side;

    int    ip0, ip1, ip2, nsid = 0;
    int    isid, itri;
    rbt_T  *stree = NULL;
    Side   *sid = NULL;

    ROUTINE(setupNeighbors);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    MALLOC(sid, Side, 3*tess->ntri);

    /* get a red-black tree for the Sides */
    status = rbtCreate(1000, &stree);
    if (stree == NULL) {
        printf("ERROR:: stree could not be allocated in routine setupNeighbors\a\n");
        exit(0);
    }

    /* loop through all Sides of all Triangles and check to see if the Side
          is already in the Side list.  if not, add a new Side. */
    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

        ip0 = tess->trip[3*itri  ];
        ip1 = tess->trip[3*itri+1];
        ip2 = tess->trip[3*itri+2];

        isid = rbtSearch(stree, ip1, ip0, 0);
        if (isid >= 0) {
            sid[isid].rtri = itri;
            sid[isid].rsid = 2;
        } else {
            isid = rbtInsert(stree, ip0, ip1, 0);

            sid[isid].ltri = itri;
            sid[isid].lsid =  2;
            sid[isid].rtri = -1;
            sid[isid].rsid = -1;
            sid[isid].bpnt = ip0;
            sid[isid].epnt = ip1;

            nsid = MAX(nsid, isid+1);
        }

        isid = rbtSearch(stree, ip2, ip1, 0);
        if (isid >= 0) {
            sid[isid].rtri = itri;
            sid[isid].rsid = 0;
        } else {
            isid = rbtInsert(stree, ip1, ip2, 0);

            sid[isid].ltri = itri;
            sid[isid].lsid =  0;
            sid[isid].rtri = -1;
            sid[isid].rsid = -1;
            sid[isid].bpnt = ip1;
            sid[isid].epnt = ip2;

            nsid = MAX(nsid, isid+1);
        }

        isid = rbtSearch(stree, ip0, ip2, 0);
        if (isid >= 0) {
            sid[isid].rtri = itri;
            sid[isid].rsid = 1;
        } else {
            isid = rbtInsert(stree, ip2, ip0, 0);

            sid[isid].ltri = itri;
            sid[isid].lsid =  1;
            sid[isid].rtri = -1;
            sid[isid].rsid = -1;
            sid[isid].bpnt = ip2;
            sid[isid].epnt = ip0;

            nsid = MAX(nsid, isid+1);
        }
    }

    /* initialize the neighbors */
    for (itri = 0; itri < tess->ntri; itri++) {
        tess->trit[3*itri  ] = -1;
        tess->trit[3*itri+1] = -1;
        tess->trit[3*itri+2] = -1;
    }

    /* apply the neighbor information to the Triangles */
    for (isid = 0; isid < nsid; isid++) {
        if (       sid[isid].lsid == 0) {
            tess->trit[3*sid[isid].ltri  ] = sid[isid].rtri;
        } else if (sid[isid].lsid == 1) {
            tess->trit[3*sid[isid].ltri+1] = sid[isid].rtri;
        } else if (sid[isid].lsid == 2) {
            tess->trit[3*sid[isid].ltri+2] = sid[isid].rtri;
        }

        if (       sid[isid].rsid == 0) {
            tess->trit[3*sid[isid].rtri  ] = sid[isid].ltri;
        } else if (sid[isid].rsid == 1) {
            tess->trit[3*sid[isid].rtri+1] = sid[isid].ltri;
        } else if (sid[isid].rsid == 2) {
            tess->trit[3*sid[isid].rtri+2] = sid[isid].ltri;
        }
    }

    /* release the Side tree */
    rbtDelete(stree);

cleanup:
    FREE(stree);
    FREE(sid  );

    return status;
}



/*
 ******************************************************************************
 *                                                                            *
 * smfAdd - add an element to a sparse matrix                                 *
 *                                                                            *
 ******************************************************************************
 */
static int
smfAdd(smf_T   *smf,                    /* (in)  pointer to SMF */
       int     irow,                    /* (in)  row index (bias-0) */
       int     icol)                    /* (in)  column index (bias-0) */
{
    int     status = -1;                /* (out) entity index (bias-0) */

    int     ient, ilast;

    ROUTINE(smfSet);

    /* --------------------------------------------------------------- */

    ient  = irow;
    ilast = -1;

    while (ient >= 0) {
        if (smf->icol[ient] == icol) {
            status = ient;
            goto cleanup;
        }

        ilast = ient;
        ient  = smf->next[ient];
    }

    if (smf->nent >= smf->ment-2) {
        (smf->ment) += 1000;

        RALLOC(smf->a,    double, smf->ment);
        RALLOC(smf->icol, int,    smf->ment);
        RALLOC(smf->next, int,    smf->ment);
    }

    status = smf->nent;

    smf->a[   smf->nent] = 0;
    smf->icol[smf->nent] = icol;
    smf->next[smf->nent] = -1;

    smf->next[ilast] = smf->nent;

    smf->nent = smf->nent + 1;

cleanup:
    return status;
}



/*
 ******************************************************************************
 *                                                                            *
 * smfFree - free up sparse matrix storage                                    *
 *                                                                            *
 ******************************************************************************
 */
static int
smfFree(smf_T   *smf)                   /* (in)  pointer to SMF */
{
    int     status = 0;                 /* (out) return status */

    ROUTINE(smfFree);

    /* --------------------------------------------------------------- */

    FREE(smf->a   );
    FREE(smf->icol);
    FREE(smf->next);

    smf->nrow = 0;
    smf->nent = 0;
    smf->ment = 0;

//cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * smfInit - initialize sparse matrix storage                                 *
 *                                                                            *
 ******************************************************************************
 */
static int
smfInit(smf_T      *smf,                /* (in)  pointer to SMF */
        int        nrow)                /* (in)  number of rows */
{
    int     status = 0;                 /* (out) return status */

    int     ient;

    ROUTINE(smfInit);

    /* --------------------------------------------------------------- */

    smf->nrow = nrow;
    smf->nent = nrow;
    smf->ment = nrow + 1000;
    smf->a    = NULL;
    smf->icol = NULL;
    smf->next = NULL;

    MALLOC(smf->a,    double, smf->ment);
    MALLOC(smf->icol, int,    smf->ment);
    MALLOC(smf->next, int,    smf->ment);

    for (ient = 0; ient < nrow; ient++) {
        smf->a[   ient] = 0;
        smf->icol[ient] = ient;
        smf->next[ient] = -1;
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * sortTriangles - sortTriangles by color                                     *
 *                                                                            *
 ******************************************************************************
 */
int
sortTriangles(tess_T  *tess)            /* (in)  pointer ot TESS */
{
    int    status = 0;                  /* (out) return status */

    int    icolr, itri, jtri;

    ROUTINE(sortTriangles);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* loop through all colors */
    itri = 0;
    for (icolr = 0; icolr <= tess->ncolr; icolr++) {
        jtri = tess->ntri - 1;

        /* look for pairs of Triangles to swap */
        while (itri < jtri) {

            /* move itri to point to first Triangle with color != icolr */
            while (itri <= jtri && (tess->ttyp[itri] & TRI_COLOR) == icolr) {
                itri++;
            }

            /* move jtri to point to last Triangle with color == icolr */
            while (itri <= jtri && (tess->ttyp[jtri] & TRI_COLOR) != icolr) {
                jtri--;
            }

            /* swap Triangles and advance itri and jtri */
            if (itri < jtri) {
                status = swapTriangles(tess, itri, jtri);
                CHECK_STATUS(swapTriangles);

                itri++;
                jtri--;
            } else {
                break;
            }
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * splitTriangle - split Triangle and its neighbor                            *
 *                                                                            *
 ******************************************************************************
 */
int
splitTriangle(tess_T  *tess,            /* (in)  pointer to TESS */
              int     itri,             /* (in)  Triangle index (bias-0) */
              int     ipnt,             /* (in)  starting Point index (bias-0) */
              double  frac)             /* (in)  fractional distance on side to split */
{
    int    status = 0;                  /* (out) return status */

    int    ip0, ip1, ip2, ip3, il0=0, il1=0, il2=0, il3=0, il4=0;
    int    it0, it1, it2, it3, it4, it5, it6;
    double xnew, ynew, znew;

    ROUTINE(splitTriangle);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (frac < EPS06 || frac > 1-EPS06) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* find the identity of the various Triangles and Points
              ip1                           ip1
              /|\                           /|\
       it1  /  |  \  it4             it1  /  |  \  it4
          /    |    \                   /    |    \
        / itri |      \               / itri | it6  \
    ipnt       |       ip2   ==>  ipnt------ip3-------ip2
        \      | it0  /               \  it5 | it0  /
          \    |    /                   \    |    /
       it2  \  |  /  it3             it2  \  |  /  it3
              \|/                           \|/
              ip0                           ip0
    */

    if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (ipnt == tess->trip[3*itri  ]) {
        ip0 = tess->trip[3*itri+1];
        ip1 = tess->trip[3*itri+2];
        it0 = tess->trit[3*itri  ];
        it1 = tess->trit[3*itri+1];
        it2 = tess->trit[3*itri+2];

        if ((tess->ttyp[itri] & TRI_T0_LINK) != 0) il0 = 1;
        if ((tess->ttyp[itri] & TRI_T1_LINK) != 0) il1 = 1;
        if ((tess->ttyp[itri] & TRI_T2_LINK) != 0) il2 = 1;
    } else if (ipnt == tess->trip[3*itri+1]) {
        ip0 = tess->trip[3*itri+2];
        ip1 = tess->trip[3*itri  ];
        it0 = tess->trit[3*itri+1];
        it1 = tess->trit[3*itri+2];
        it2 = tess->trit[3*itri  ];

        if ((tess->ttyp[itri] & TRI_T0_LINK) != 0) il2 = 1;
        if ((tess->ttyp[itri] & TRI_T1_LINK) != 0) il0 = 1;
        if ((tess->ttyp[itri] & TRI_T2_LINK) != 0) il1 = 1;
    } else if (ipnt == tess->trip[3*itri+2]) {
        ip0 = tess->trip[3*itri  ];
        ip1 = tess->trip[3*itri+1];
        it0 = tess->trit[3*itri+2];
        it1 = tess->trit[3*itri  ];
        it2 = tess->trit[3*itri+1];

        if ((tess->ttyp[itri] & TRI_T0_LINK) != 0) il1 = 1;
        if ((tess->ttyp[itri] & TRI_T1_LINK) != 0) il2 = 1;
        if ((tess->ttyp[itri] & TRI_T2_LINK) != 0) il0 = 1;
    } else {
        printf("ipnt=%d not associated with itri=%d\n", ipnt, itri);
        status = -999;
        goto cleanup;
    }

    if (it0 < 0 || it0 >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if ((tess->ttyp[it0] & TRI_ACTIVE) == 0) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (itri == tess->trit[3*it0  ]) {
        ip2 = tess->trip[3*it0  ];
        it3 = tess->trit[3*it0+1];
        it4 = tess->trit[3*it0+2];

        if ((tess->ttyp[it0] & TRI_T1_LINK) != 0) il3 = 1;
        if ((tess->ttyp[it0] & TRI_T2_LINK) != 0) il4 = 1;
    } else if (itri == tess->trit[3*it0+1]) {
        ip2 = tess->trip[3*it0+1];
        it3 = tess->trit[3*it0+2];
        it4 = tess->trit[3*it0  ];

        if ((tess->ttyp[it0] & TRI_T0_LINK) != 0) il4 = 1;
        if ((tess->ttyp[it0] & TRI_T2_LINK) != 0) il3 = 1;
    } else if (itri == tess->trit[3*it0+2]) {
        ip2 = tess->trip[3*it0+2];
        it3 = tess->trit[3*it0  ];
        it4 = tess->trit[3*it0+1];

        if ((tess->ttyp[it0] & TRI_T0_LINK) != 0) il3 = 1;
        if ((tess->ttyp[it0] & TRI_T1_LINK) != 0) il4 = 1;
    } else {
        printf("itri=%d not associated with it0=%d\n", itri, it0);
        status = -999;
        goto cleanup;
    }

    /* create the new Point */
    xnew = (1-frac) * tess->xyz[3*ip0  ] + frac * tess->xyz[3*ip1  ];
    ynew = (1-frac) * tess->xyz[3*ip0+1] + frac * tess->xyz[3*ip1+1];
    znew = (1-frac) * tess->xyz[3*ip0+2] + frac * tess->xyz[3*ip1+2];

    status = ip3 = addPoint(tess, xnew, ynew, znew);
    CHECK_STATUS(addPoint);

    /* make itri and it0 use the new Point */
    tess->trip[3*itri  ] = ip3;
    tess->trip[3*itri+1] = ip1;
    tess->trip[3*itri+2] = ipnt;

    tess->trip[3*it0  ] = ip3;
    tess->trip[3*it0+1] = ip0;
    tess->trip[3*it0+2] = ip2;

    /* create the new Triangles */
    status = it5 = addTriangle(tess, ip3, ipnt, ip0, -1, -1, -1);
    CHECK_STATUS(addTriangle);

    status = it6 = addTriangle(tess, ip3, ip2,  ip1, -1, -1, -1);
    CHECK_STATUS(addTriangle);

    tess->ttyp[it5] = (tess->ttyp[it5] & ~TRI_COLOR) | (tess->ttyp[itri] & TRI_COLOR);
    tess->ttyp[it6] = (tess->ttyp[it6] & ~TRI_COLOR) | (tess->ttyp[it0 ] & TRI_COLOR);

    /* set up the neighbor information for the interior Triangles */
    tess->trit[3*itri  ] = it1;
    tess->trit[3*itri+1] = it5;
    tess->trit[3*itri+2] = it6;

    tess->trit[3*it0  ] = it3;
    tess->trit[3*it0+1] = it6;
    tess->trit[3*it0+2] = it5;

    tess->trit[3*it5  ] = it2;
    tess->trit[3*it5+1] = it0;
    tess->trit[3*it5+2] = itri;

    tess->trit[3*it6  ] = it4;
    tess->trit[3*it6+1] = itri;
    tess->trit[3*it6+2] = it0;

    /* set up the link information for the interior Triangles */
    tess->ttyp[itri] &= ~TRI_LINK;
    tess->ttyp[it0 ] &= ~TRI_LINK;
    tess->ttyp[it5 ] &= ~TRI_LINK;
    tess->ttyp[it6 ] &= ~TRI_LINK;

    if (il0 == 1) {
        status = createLink(tess, itri, 2);
        CHECK_STATUS(createLink);
        status = createLink(tess, it0,  2);
        CHECK_STATUS(createLink);
        status = createLink(tess, it5,  1);
        CHECK_STATUS(createLink);
        status = createLink(tess, it6,  1);
        CHECK_STATUS(createLink);
    }
    if (il1 == 1) {
        status = createLink(tess, itri, 0);
        CHECK_STATUS(createLink);
    }
    if (il2 == 1) {
        status = createLink(tess, it5,  0);
        CHECK_STATUS(createLink);
    }
    if (il3 == 1) {
        status = createLink(tess, it0,  0);
        CHECK_STATUS(createLink);
    }
    if (il4 == 1) {
        status = createLink(tess, it6,  0);
        CHECK_STATUS(createLink);
    }

    /* make sure interior Triangles are connected properly */
    status = connectNeighbors(tess, itri);
    CHECK_STATUS(connectNeighbors);

    status = connectNeighbors(tess, it0 );
    CHECK_STATUS(connectNeighbors);

    status = connectNeighbors(tess, it5 );
    CHECK_STATUS(connectNeighbors);

    status = connectNeighbors(tess, it6 );
    CHECK_STATUS(connectNeighbors);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * swapTriangles - swap Triangles                                             *
 *                                                                            *
 ******************************************************************************
 */
int
swapTriangles(tess_T  *tess,            /* (in)  pointer ot TESS */
              int     itri,             /* (in)  index of first  Triangle (bias-0) */
              int     jtri)             /* (in)  index of second Triangle (bias-0) */
{
    int    status = 0;                  /* (out) return status */

    int   swap, ktri;

//    ROUTINE(swapTriangles);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (itri < 0 || itri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    } else if (jtri < 0 || jtri >= tess->ntri) {
        status = TESS_BAD_TRIANGLE_INDEX;
        goto cleanup;
    }

    /* swap the data */
    swap                 = tess->trip[3*itri  ];
    tess->trip[3*itri  ] = tess->trip[3*jtri  ];
    tess->trip[3*jtri  ] = swap;

    swap                 = tess->trip[3*itri+1];
    tess->trip[3*itri+1] = tess->trip[3*jtri+1];
    tess->trip[3*jtri+1] = swap;

    swap                 = tess->trip[3*itri+2];
    tess->trip[3*itri+2] = tess->trip[3*jtri+2];
    tess->trip[3*jtri+2] = swap;

    swap                 = tess->trit[3*itri  ];
    tess->trit[3*itri  ] = tess->trit[3*jtri  ];
    tess->trit[3*jtri  ] = swap;

    swap                 = tess->trit[3*itri+1];
    tess->trit[3*itri+1] = tess->trit[3*jtri+1];
    tess->trit[3*jtri+1] = swap;

    swap                 = tess->trit[3*itri+2];
    tess->trit[3*itri+2] = tess->trit[3*jtri+2];
    tess->trit[3*jtri+2] = swap;

    swap                 = tess->ttyp[  itri  ];
    tess->ttyp[  itri  ] = tess->ttyp[  jtri  ];
    tess->ttyp[  jtri  ] = swap;

    /* neighbors that used to point to itri should now point to jtri */
    ktri = tess->trit[3*jtri  ];
    if (ktri >= 0) {
        if (tess->trit[3*ktri  ] == itri) tess->trit[3*ktri  ] = jtri;
        if (tess->trit[3*ktri+1] == itri) tess->trit[3*ktri+1] = jtri;
        if (tess->trit[3*ktri+2] == itri) tess->trit[3*ktri+2] = jtri;
    }

    ktri = tess->trit[3*jtri+1];
    if (ktri >= 0) {
        if (tess->trit[3*ktri  ] == itri) tess->trit[3*ktri  ] = jtri;
        if (tess->trit[3*ktri+1] == itri) tess->trit[3*ktri+1] = jtri;
        if (tess->trit[3*ktri+2] == itri) tess->trit[3*ktri+2] = jtri;
    }

    ktri = tess->trit[3*jtri+2];
    if (ktri >= 0) {
        if (tess->trit[3*ktri  ] == itri) tess->trit[3*ktri  ] = jtri;
        if (tess->trit[3*ktri+1] == itri) tess->trit[3*ktri+1] = jtri;
        if (tess->trit[3*ktri+2] == itri) tess->trit[3*ktri+2] = jtri;
    }

    /* neighbors that used to point to jtri should now point to itri */
    ktri = tess->trit[3*itri  ];
    if (ktri >= 0) {
        if (tess->trit[3*ktri  ] == jtri) tess->trit[3*ktri  ] = itri;
        if (tess->trit[3*ktri+1] == jtri) tess->trit[3*ktri+1] = itri;
        if (tess->trit[3*ktri+2] == jtri) tess->trit[3*ktri+2] = itri;
    }

    ktri = tess->trit[3*itri+1];
    if (ktri >= 0) {
        if (tess->trit[3*ktri  ] == jtri) tess->trit[3*ktri  ] = itri;
        if (tess->trit[3*ktri+1] == jtri) tess->trit[3*ktri+1] = itri;
        if (tess->trit[3*ktri+2] == jtri) tess->trit[3*ktri+2] = itri;
    }

    ktri = tess->trit[3*itri+2];
    if (ktri >= 0) {
        if (tess->trit[3*ktri  ] == jtri) tess->trit[3*ktri  ] = itri;
        if (tess->trit[3*ktri+1] == jtri) tess->trit[3*ktri+1] = itri;
        if (tess->trit[3*ktri+2] == jtri) tess->trit[3*ktri+2] = itri;
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * triNormal - compute the area and normal of a Triangle                      *
 *                                                                            *
 ******************************************************************************
 */
static void
triNormal(tess_T  *tess,                /* (in)  pointer to TESS */
          int     ip0,                  /* (in)  first  Point index (bias-0) */
          int     ip1,                  /* (in)  second Point index (bias-0) */
          int     ip2,                  /* (in)  third  Point index (bias-0) */
          double  *area,                /* (out) area of Triangle */
          double  norm[])               /* (out) unit normal */
{
//    ROUTINE(triNormal);

    /*----------------------------------------------------------------*/

    if ((ip0 < 0) || (ip1 < 0) || (ip2 < 0)) {
        *area   = 0.0;
        norm[0] = 0.0;
        norm[1] = 0.0;
        norm[2] = 0.0;

    } else {
        norm[0] = (tess->xyz[3*ip1+1] - tess->xyz[3*ip0+1])
                * (tess->xyz[3*ip2+2] - tess->xyz[3*ip0+2])
                - (tess->xyz[3*ip2+1] - tess->xyz[3*ip0+1])
                * (tess->xyz[3*ip1+2] - tess->xyz[3*ip0+2]);
        norm[1] = (tess->xyz[3*ip1+2] - tess->xyz[3*ip0+2])
                * (tess->xyz[3*ip2  ] - tess->xyz[3*ip0  ])
                - (tess->xyz[3*ip2+2] - tess->xyz[3*ip0+2])
                * (tess->xyz[3*ip1  ] - tess->xyz[3*ip0  ]);
        norm[2] = (tess->xyz[3*ip1  ] - tess->xyz[3*ip0  ])
                * (tess->xyz[3*ip2+1] - tess->xyz[3*ip0+1])
                - (tess->xyz[3*ip2  ] - tess->xyz[3*ip0  ])
                * (tess->xyz[3*ip1+1] - tess->xyz[3*ip0+1]);

        *area = sqrt(SQR(norm[0]) + SQR(norm[1]) + SQR(norm[2]));

        if ((*area) > 0) {
            norm[0] /= (*area);
            norm[1] /= (*area);
            norm[2] /= (*area);
        }
    }

//cleanup:
    return;
}


/*
 ******************************************************************************
 *                                                                            *
 * turn - compute turning angle                                               *
 *                                                                            *
 ******************************************************************************
 */
static double
turn(tess_T  *tess,                     /* (in)  pointer to TESS */
     int     ipnt,                      /* (in)  first  Point index (bias-0) */
     int     jpnt,                      /* (in)  second Point index (bias-0) */
     int     kpnt,                      /* (in)  third  Point index (bias-0) */
     int     itri)                      /* (in)  index of Tringle for normal */
{
//    ROUTINE(turn);

    int    ip0, ip1, ip2;
    double veca[3], vecb[3], vecc[3], triple, dot, answer;

    /*----------------------------------------------------------------*/

    /* find normalized vectors a=i->j abd b=j->k */
    veca[0]  = tess->xyz[3*jpnt  ] - tess->xyz[3*ipnt  ];
    veca[1]  = tess->xyz[3*jpnt+1] - tess->xyz[3*ipnt+1];
    veca[2]  = tess->xyz[3*jpnt+2] - tess->xyz[3*ipnt+2];

    vecb[0] = tess->xyz[3*kpnt  ] - tess->xyz[3*jpnt  ];
    vecb[1] = tess->xyz[3*kpnt+1] - tess->xyz[3*jpnt+1];
    vecb[2] = tess->xyz[3*kpnt+2] - tess->xyz[3*jpnt+2];

    /* find normalized Triangle normal */
    ip0 = tess->trip[3*itri  ];
    ip1 = tess->trip[3*itri+1];
    ip2 = tess->trip[3*itri+2];

    vecc[0]  = (tess->xyz[3*ip1+1] - tess->xyz[3*ip0+1])
             * (tess->xyz[3*ip2+2] - tess->xyz[3*ip0+2])
             - (tess->xyz[3*ip2+1] - tess->xyz[3*ip0+1])
             * (tess->xyz[3*ip1+2] - tess->xyz[3*ip0+2]);
    vecc[1]  = (tess->xyz[3*ip1+2] - tess->xyz[3*ip0+2])
             * (tess->xyz[3*ip2  ] - tess->xyz[3*ip0  ])
             - (tess->xyz[3*ip2+2] - tess->xyz[3*ip0+2])
             * (tess->xyz[3*ip1  ] - tess->xyz[3*ip0  ]);
    vecc[2]  = (tess->xyz[3*ip1  ] - tess->xyz[3*ip0  ])
             * (tess->xyz[3*ip2+1] - tess->xyz[3*ip0+1])
             - (tess->xyz[3*ip2  ] - tess->xyz[3*ip0  ])
             * (tess->xyz[3*ip1+1] - tess->xyz[3*ip0+1]);

    /* compute normalized triple product */
    triple = ( veca[0]*vecb[1]*vecc[2] + veca[1]*vecb[2]*vecc[0] + veca[2]*vecb[0]*vecc[1]
             - veca[0]*vecb[2]*vecc[1] - veca[1]*vecb[0]*vecc[2] - veca[2]*vecb[1]*vecc[0])
           / sqrt(veca[0]*veca[0] + veca[1]*veca[1] + veca[2]*veca[2])
           / sqrt(vecb[0]*vecb[0] + vecb[1]*vecb[1] + vecb[2]*vecb[2])
           / sqrt(vecc[0]*vecc[0] + vecc[1]*vecc[1] + vecc[2]*vecc[2]);

    /* compute dot product */
    dot = veca[0]*vecb[0] + veca[1]*vecb[1] + veca[2]*vecb[2];

    if (dot > 0) {
        answer = triple;
    } else if (triple > 0) {
        answer = +PI - triple;
    } else {
        answer = -PI - triple;
    }

//cleanup:
    return answer;
}


/*
 ******************************************************************************
 *                                                                            *
 * UVtoXYZ - evaluate at a given parametric coordinate                        *
 *                                                                            *
 ******************************************************************************
 */
int
UVtoXYZ(tess_T  *tess,                  /* (in)  pointer ot TESS */
        int     icolr,                  /* (in)  color of Triangles */
        double  uv_in[],                /* (in)  parametric coordinates */
        double  xyz_out[])              /* (out) physical coordinates */
{
    int    status = 0;                  /* (out) return status */

    int    jtri, ip0, ip1, ip2;
    double u02, u12, uu2, v02, v12, vv2, D, s0, s1, s2;

//    ROUTINE(UVtoXYZ);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (icolr < 0 || icolr > 255) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* default output */
    xyz_out[0] = HUGEQ;
    xyz_out[1] = HUGEQ;
    xyz_out[2] = HUGEQ;

    /* loop through all Triangles with icolr */
    for (jtri = 0; jtri < tess->ntri; jtri++) {
        if ((tess->ttyp[jtri] & TRI_COLOR) != icolr) continue;

        /* determine barycentric coordinates of point closest to xyz_in */
        ip0 = tess->trip[3*jtri  ];
        ip1 = tess->trip[3*jtri+1];
        ip2 = tess->trip[3*jtri+2];

        u02 = tess->uv[2*ip0  ] - tess->uv[2*ip2  ];
        v02 = tess->uv[2*ip0+1] - tess->uv[2*ip2+1];
        u12 = tess->uv[2*ip1  ] - tess->uv[2*ip2  ];
        v12 = tess->uv[2*ip1+1] - tess->uv[2*ip2+1];
        uu2 = uv_in[0]          - tess->uv[2*ip2  ];
        vv2 = uv_in[1]          - tess->uv[2*ip2+1];

        D  =  u02 * v12 - v02 * u12;
        if (fabs(D) > EPS20) {
            s0 = (uu2 * v12 - vv2 * u12) / D;
            s1 = (u02 * vv2 - v02 * uu2) / D;
            s2 = 1 - s0 - s1;

            if (s0 > -EPS06 && s1 > -EPS06 && s2 > -EPS06) {
                xyz_out[0] = s0 * tess->xyz[3*ip0  ]
                           + s1 * tess->xyz[3*ip1  ]
                           + s2 * tess->xyz[3*ip2  ];
                xyz_out[1] = s0 * tess->xyz[3*ip0+1]
                           + s1 * tess->xyz[3*ip1+1]
                           + s2 * tess->xyz[3*ip2+1];
                xyz_out[2] = s0 * tess->xyz[3*ip0+2]
                           + s1 * tess->xyz[3*ip1+2]
                           + s2 * tess->xyz[3*ip2+2];
                goto cleanup;
            }
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * writeTriAscii - write an ASCII triangle file                               *
 *                                                                            *
 ******************************************************************************
 */
int
writeTriAscii(tess_T  *tess,            /* (in)  pointer to TESS */
              char    *filename)        /* (in)  name of file */
{
    int    status = 0;                  /* (out) return status */

    int    ipnt, itri;
    FILE   *fp = NULL;

//    ROUTINE(writeTriAscii);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* open the output file */
    fp = fopen(filename, "w");
    if (fp == NULL) {
        status = TESS_BAD_FILE_NAME;
        goto cleanup;
    }

    /* write the Points */
    fprintf(fp, "%10d\n", tess->npnt);

    for (ipnt = 0; ipnt < tess->npnt; ipnt++) {
        fprintf(fp, "%20.10e%20.10e%20.10e\n",
                tess->xyz[3*ipnt], tess->xyz[3*ipnt+1], tess->xyz[3*ipnt+2]);
    }

    /* write the Triangles */
    fprintf(fp, "%10d\n", tess->ntri);

    for (itri = 0; itri < tess->ntri; itri++) {
        fprintf(fp, "%10d%10d%10d%10d%10d%10d%5d\n",
                tess->trip[3*itri], tess->trip[3*itri+1], tess->trip[3*itri+2],
                tess->trit[3*itri], tess->trit[3*itri+1], tess->trit[3*itri+2],
                tess->ttyp[itri] & TRI_COLOR);
    }

    /* close the file */
    fclose(fp);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * writeStlAscii - write an ASCII stl file                                    *
 *                                                                            *
 ******************************************************************************
 */
int
writeStlAscii(tess_T  *tess,            /* (in)  pointer to TESS */
              char    *filename)        /* (in)  name of file */
{
    int    status = 0;                  /* (out) return status */

    int    itri, ip0, ip1, ip2;
    double area, norm[3];
    FILE   *fp = NULL;

//    ROUTINE(writeStlAscii);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* open the output file */
    fp = fopen(filename, "w");
    if (fp == NULL) {
        status = TESS_BAD_FILE_NAME;
        goto cleanup;
    }

    /* create header */
    fprintf(fp, "solid OBJECT\n");

    /* write the Triangles */
    for (itri = 0; itri < tess->ntri; itri++) {
        if (tess->ttyp[itri] & TRI_VISIBLE) {
            ip0 = tess->trip[3*itri  ];
            ip1 = tess->trip[3*itri+1];
            ip2 = tess->trip[3*itri+2];

            triNormal(tess, ip0, ip1, ip2, &area, norm);

            fprintf(fp, "  facet normal %14.6e %14.6e %14.6e\n",
                    norm[0], norm[1], norm[2]);
            fprintf(fp, "    outer loop\n");
            fprintf(fp, "      vertex   %14.6e %14.6e %14.6e\n",
                    tess->xyz[3*ip0  ], tess->xyz[3*ip0+1], tess->xyz[3*ip0+2]);
            fprintf(fp, "      vertex   %14.6e %14.6e %14.6e\n",
                    tess->xyz[3*ip1  ], tess->xyz[3*ip1+1], tess->xyz[3*ip1+2]);
            fprintf(fp, "      vertex   %14.6e %14.6e %14.6e\n",
                    tess->xyz[3*ip2  ], tess->xyz[3*ip2+1], tess->xyz[3*ip2+2]);
            fprintf(fp, "    endloop\n");
            fprintf(fp, "  endfacet\n");
        }
    }

    fprintf(fp, "endsolid OBJECT\n");
    fclose(fp);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * writeStlBinary - write a binary stl file                                   *
 *                                                                            *
 ******************************************************************************
 */
int
writeStlBinary(tess_T  *tess,           /* (in)  pointer to TESS */
               char    *filename)       /* (in)  name of file */
{
    int    status = 0;                  /* (out) return status */

    int    itri, ip0, ip1, ip2;
    UINT16 icolr;
    UINT32 ntri_visible;
    double area, normal[3];
    REAL32 norm[3], vert[3];
    char   header[81];
    FILE   *fp = NULL;

//    ROUTINE(writeStlBinary);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    }

    /* open the output file */
    fp = fopen(filename, "wb");
    if (fp == NULL) {
        status = TESS_BAD_FILE_NAME;
        goto cleanup;
    }

    /* create header */
    sprintf(header, "written by StlEdit, ncolr=%d", tess->ncolr);
    (void) fwrite(header, sizeof(char), 80, fp);

    /* number of Triangles */
    ntri_visible = 0;
    for (itri = 0; itri < tess->ntri; itri++) {
        if (tess->ttyp[itri] & TRI_VISIBLE) {
            ntri_visible++;
        }
    }
    (void) fwrite(&ntri_visible, sizeof(UINT32), 1, fp);

    /* write the Triangles */
    for (itri = 0; itri < tess->ntri; itri++) {
        if (tess->ttyp[itri] & TRI_VISIBLE) {
            ip0 = tess->trip[3*itri  ];
            ip1 = tess->trip[3*itri+1];
            ip2 = tess->trip[3*itri+2];

            triNormal(tess, ip0, ip1, ip2, &area, normal);
            norm[0] = normal[0];
            norm[1] = normal[1];
            norm[2] = normal[2];
            (void) fwrite(norm, sizeof(REAL32), 3, fp);

            vert[0] = tess->xyz[3*ip0  ];
            vert[1] = tess->xyz[3*ip0+1];
            vert[2] = tess->xyz[3*ip0+2];
            (void) fwrite(vert, sizeof(REAL32), 3, fp);

            vert[0] = tess->xyz[3*ip1  ];
            vert[1] = tess->xyz[3*ip1+1];
            vert[2] = tess->xyz[3*ip1+2];
            (void) fwrite(vert, sizeof(REAL32), 3, fp);

            vert[0] = tess->xyz[3*ip2  ];
            vert[1] = tess->xyz[3*ip2+1];
            vert[2] = tess->xyz[3*ip2+2];
            (void) fwrite(vert, sizeof(REAL32), 3, fp);

            icolr = tess->ttyp[itri] & TRI_COLOR;
            (void) fwrite(&icolr, sizeof(UINT16), 1, fp);
        }
    }

    fclose(fp);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * XYZtoUVXYZ - find point nearest to the Tessellation                        *
 *                                                                            *
 ******************************************************************************
 */
int
XYZtoUVXYZ(tess_T  *tess,               /* (in)  pointer to TESS */
           int     icolr,               /* (in)  color of Triangles */
           double  xyz_in[],            /* (in)  physical coordinates */
           double  uv_out[],            /* (out) parametric coordinates (or UNDEF) */
           double  xyz_out[])           /* (out) coordinates of output coordinates */
{
    int    status = 0;                  /* (out) return status */

    int       jtri, ip0, ip1, ip2;
    double    s0, s1, s01, dbest, dbest2, dtest2, xtest, ytest, ztest;
    double    x02, y02, z02, x12, y12, z12, xx2, yy2, zz2, A, B, C, D, E, F, G;

//    ROUTINE(XYZtoUVXYZ);

    /* --------------------------------------------------------------- */

    if (tess == NULL) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (tess->magic != TESS_MAGIC) {
        status = TESS_NOT_A_TESS;
        goto cleanup;
    } else if (icolr < 0 || icolr > 255) {
        status = TESS_BAD_VALUE;
        goto cleanup;
    }

    /* default output */
    xyz_out[0] = xyz_in[0];
    xyz_out[1] = xyz_in[1];
    xyz_out[2] = xyz_in[2];

    dbest  = 1000;
    dbest2 = SQR(dbest);

    /* loop through all Triangles with icolr */
    for (jtri = 0; jtri < tess->ntri; jtri++) {
        if ((tess->ttyp[jtri] & TRI_COLOR) != icolr) continue;

        /* bounding box check */
        if (xyz_in[0] < tess->bbox[6*jtri  ]-dbest ||
            xyz_in[0] > tess->bbox[6*jtri+1]+dbest ||
            xyz_in[1] < tess->bbox[6*jtri+2]-dbest ||
            xyz_in[1] > tess->bbox[6*jtri+3]+dbest ||
            xyz_in[2] < tess->bbox[6*jtri+4]-dbest ||
            xyz_in[2] > tess->bbox[6*jtri+5]+dbest   ) continue;

        /* determine barycentric coordinates of point closest to xyz_in */
        ip0 = tess->trip[3*jtri  ];
        ip1 = tess->trip[3*jtri+1];
        ip2 = tess->trip[3*jtri+2];

        x02 = tess->xyz[3*ip0  ] - tess->xyz[3*ip2  ];
        y02 = tess->xyz[3*ip0+1] - tess->xyz[3*ip2+1];
        z02 = tess->xyz[3*ip0+2] - tess->xyz[3*ip2+2];
        x12 = tess->xyz[3*ip1  ] - tess->xyz[3*ip2  ];
        y12 = tess->xyz[3*ip1+1] - tess->xyz[3*ip2+1];
        z12 = tess->xyz[3*ip1+2] - tess->xyz[3*ip2+2];
        xx2 = xyz_in[0]          - tess->xyz[3*ip2  ];
        yy2 = xyz_in[1]          - tess->xyz[3*ip2+1];
        zz2 = xyz_in[2]          - tess->xyz[3*ip2+2];

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

        xtest = tess->xyz[3*ip2  ] + s0 * x02 + s1 * x12;
        ytest = tess->xyz[3*ip2+1] + s0 * y02 + s1 * y12;
        ztest = tess->xyz[3*ip2+2] + s0 * z02 + s1 * z12;

        /* if point is closer than anything we have seen so far, remember this */
        dtest2 = SQR(xtest-xyz_in[0]) + SQR(ytest-xyz_in[1]) + SQR(ztest-xyz_in[2]);

        if (dtest2 < dbest2) {
            xyz_out[0] = xtest;
            xyz_out[1] = ytest;
            xyz_out[2] = ztest;

            uv_out[0] =    s0     * tess->uv[2*ip0  ]
                      +       s1  * tess->uv[2*ip1  ]
                      + (1-s0-s1) * tess->uv[2*ip2  ];
            uv_out[1] =    s0     * tess->uv[2*ip0+1]
                      +       s1  * tess->uv[2*ip1+1]
                      + (1-s0-s1) * tess->uv[2*ip2+1];

            dbest2     = dtest2;
            dbest      = sqrt(dbest2);
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * buildOctree - build octree with no more than given number of triangles     *
 *                                                                            *
 ******************************************************************************
 */
static int
buildOctree(tess_T *tess,
            int    nmax,
            oct_T  *tree)
{
    int    status = 0;                  /* (out) return status */

    int    itri, ipnt;

    ROUTINE(buildOctree);

    /* --------------------------------------------------------------- */

    /* initialize this tree */
    tree->npnt = tess->npnt;
    tree->xyz  = tess->xyz;
    tree->ntri = 0;
    tree->trip = NULL;
    tree->xcent = 0;
    tree->ycent = 0;
    tree->zcent = 0;
    tree->child = NULL;

    /* copy the Triangles into the base tree */
    MALLOC(tree->trip, int, 3*tess->ntri);

    for (itri = 0; itri < tess->ntri; itri++) {
        if ((tess->ttyp[itri] & TRI_ACTIVE) == 0) continue;

        tree->trip[3*tree->ntri  ] = tess->trip[3*itri  ];
        tree->trip[3*tree->ntri+1] = tess->trip[3*itri+1];
        tree->trip[3*tree->ntri+2] = tess->trip[3*itri+2];
        tree->ntri++;
    }

    /* compute the centroid */
    for (ipnt = 1; ipnt < tess->npnt; ipnt++) {
        tree->xcent += tess->xyz[3*ipnt  ];
        tree->ycent += tess->xyz[3*ipnt+1];
        tree->zcent += tess->xyz[3*ipnt+2];
    }

    tree->xcent /= tess->npnt;
    tree->ycent /= tess->npnt;
    tree->zcent /= tess->npnt;

    /* now refine the tree until all children have fewer than nmax Triangles */
    status = refineOctree(tree, nmax, 0);
    CHECK_STATUS(refineOctree);

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * refineOctree - refine octree until no leaf has more than given number of triangles *
 *                                                                            *
 ******************************************************************************
 */
static int
refineOctree(oct_T  *tree,
             int    nmax,
             int    depth)
{
    int    status = 0;                  /* (out) return status */

    int    ichild, itri, jtri, ip0, ip1, ip2, ntri_parent;
    double xmin, xmax, ymin, ymax, zmin, zmax;

    ROUTINE(refineOctree);

    /* --------------------------------------------------------------- */

    ntri_parent = tree->ntri;

    /* octant is a leaf if it contains nmax or fewer Triangles */
    if (tree->ntri <= nmax) {
        goto cleanup;
    }

    /* we need to refine, so set up the 8 children */
    MALLOC(tree->child, oct_T, 8);

    /* initialize the children, with more Triangle storage than needed */
    for (ichild = 0; ichild < 8; ichild++) {
        tree->child[ichild].npnt  = tree->npnt;
        tree->child[ichild].xyz   = tree->xyz;
        tree->child[ichild].ntri  = 0;
        tree->child[ichild].trip  = NULL;
        tree->child[ichild].xcent = 0;
        tree->child[ichild].ycent = 0;
        tree->child[ichild].zcent = 0;
        tree->child[ichild].child = NULL;

        MALLOC(tree->child[ichild].trip, int, 3*tree->ntri);
    }

    /* loop through all the Triangles and assign it to one or more children */
    for (itri = 0; itri < tree->ntri; itri++) {
        ip0 = tree->trip[3*itri  ];
        ip1 = tree->trip[3*itri+1];
        ip2 = tree->trip[3*itri+2];

        xmin = MIN(MIN(tree->xyz[3*ip0  ], tree->xyz[3*ip1  ]), tree->xyz[3*ip2  ]);
        xmax = MAX(MAX(tree->xyz[3*ip0  ], tree->xyz[3*ip1  ]), tree->xyz[3*ip2  ]);
        ymin = MIN(MIN(tree->xyz[3*ip0+1], tree->xyz[3*ip1+1]), tree->xyz[3*ip2+1]);
        ymax = MAX(MAX(tree->xyz[3*ip0+1], tree->xyz[3*ip1+1]), tree->xyz[3*ip2+1]);
        zmin = MIN(MIN(tree->xyz[3*ip0+2], tree->xyz[3*ip1+2]), tree->xyz[3*ip2+2]);
        zmax = MAX(MAX(tree->xyz[3*ip0+2], tree->xyz[3*ip1+2]), tree->xyz[3*ip2+2]);

        /* octant 0 */
        if (xmin < tree->xcent && ymin < tree->ycent && zmin < tree->zcent) {
            jtri = tree->child[0].ntri;

            tree->child[0].xcent += (xmin + xmax) / 2;
            tree->child[0].ycent += (ymin + ymax) / 2;
            tree->child[0].zcent += (zmin + zmax) / 2;

            tree->child[0].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[0].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[0].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[0].ntri++;
        }

        /* octant 1 */
        if (xmax > tree->xcent && ymin < tree->ycent && zmin < tree->zcent) {
            jtri = tree->child[1].ntri;

            tree->child[1].xcent += (xmin + xmax) / 2;
            tree->child[1].ycent += (ymin + ymax) / 2;
            tree->child[1].zcent += (zmin + zmax) / 2;

            tree->child[1].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[1].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[1].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[1].ntri++;
        }

        /* octant 2 */
        if (xmin < tree->xcent && ymax > tree->ycent && zmin < tree->zcent) {
            jtri = tree->child[2].ntri;

            tree->child[2].xcent += (xmin + xmax) / 2;
            tree->child[2].ycent += (ymin + ymax) / 2;
            tree->child[2].zcent += (zmin + zmax) / 2;

            tree->child[2].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[2].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[2].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[2].ntri++;
        }

        /* octant 3 */
        if (xmax > tree->xcent && ymax > tree->ycent && zmin < tree->zcent) {
            jtri = tree->child[3].ntri;

            tree->child[3].xcent += (xmin + xmax) / 2;
            tree->child[3].ycent += (ymin + ymax) / 2;
            tree->child[3].zcent += (zmin + zmax) / 2;

            tree->child[3].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[3].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[3].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[3].ntri++;
        }

        /* octant 4 */
        if (xmin < tree->xcent && ymin < tree->ycent && zmax > tree->zcent) {
            jtri = tree->child[4].ntri;

            tree->child[4].xcent += (xmin + xmax) / 2;
            tree->child[4].ycent += (ymin + ymax) / 2;
            tree->child[4].zcent += (zmin + zmax) / 2;

            tree->child[4].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[4].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[4].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[4].ntri++;
        }

        /* octant 5 */
        if (xmax > tree->xcent && ymin < tree->ycent && zmax > tree->zcent) {
            jtri = tree->child[5].ntri;

            tree->child[5].xcent += (xmin + xmax) / 2;
            tree->child[5].ycent += (ymin + ymax) / 2;
            tree->child[5].zcent += (zmin + zmax) / 2;

            tree->child[5].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[5].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[5].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[5].ntri++;
        }

        /* octant 6 */
        if (xmin < tree->xcent && ymax > tree->ycent && zmax > tree->zcent) {
            jtri = tree->child[6].ntri;

            tree->child[6].xcent += (xmin + xmax) / 2;
            tree->child[6].ycent += (ymin + ymax) / 2;
            tree->child[6].zcent += (zmin + zmax) / 2;

            tree->child[6].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[6].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[6].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[6].ntri++;
        }

        /* octant 7 */
        if (xmax > tree->xcent && ymax > tree->ycent && zmax > tree->zcent) {
            jtri = tree->child[7].ntri;

            tree->child[7].xcent += (xmin + xmax) / 2;
            tree->child[7].ycent += (ymin + ymax) / 2;
            tree->child[7].zcent += (zmin + zmax) / 2;

            tree->child[7].trip[3*jtri  ] = tree->trip[3*itri  ];
            tree->child[7].trip[3*jtri+1] = tree->trip[3*itri+1];
            tree->child[7].trip[3*jtri+2] = tree->trip[3*itri+2];
            tree->child[7].ntri++;
        }
    }

    /* normalize the centroid in each child */
    for (ichild = 0; ichild < 8; ichild++) {
        if (tree->child[ichild].ntri > 0) {
            tree->child[ichild].xcent /= tree->child[ichild].ntri;
            tree->child[ichild].ycent /= tree->child[ichild].ntri;
            tree->child[ichild].zcent /= tree->child[ichild].ntri;
        } else {
            tree->child[ichild].xcent  = 0;
            tree->child[ichild].ycent  = 0;
            tree->child[ichild].zcent  = 0;
        }
    }

    /* free up Triangle storage in tree */
    tree->ntri = 0;
    FREE(tree->trip);

    /* try to refine each child */
    for (ichild = 0; ichild < 8; ichild++) {
        if        (tree->child[ichild].ntri == 0) {
        } else if (tree->child[ichild].ntri >= ntri_parent) {
            printf("WARNING:: recursion stopping because tree->child[%d].ntri=%7d and ntri_parent=%7d  (depth=%d)\n",
                   ichild, tree->child[ichild].ntri, ntri_parent, depth);
        } else {
            status = refineOctree(&(tree->child[ichild]), nmax, depth+1);
            CHECK_STATUS(refineOctree);
        }
    }

cleanup:
    return status;
}


/*
 ******************************************************************************
 *                                                                            *
 * removeOctree - remove the contents of an octree                            *
 *                                                                            *
 ******************************************************************************
 */
static int
removeOctree(oct_T  *tree)
{
    int    status = 0;                  /* (out) return status */

    int    ichild;

    ROUTINE(removeOctree);

    /* --------------------------------------------------------------- */

    /* if the octree does not exist, return immediately */
    if (tree == NULL) {
        goto cleanup;
    }

    /* if the octree has children, remove them */
    if (tree->child != NULL) {
        for (ichild = 0; ichild < 8; ichild++) {
            status = removeOctree(&(tree->child[ichild]));
            CHECK_STATUS(removeOctree);
        }

        FREE(tree->child);
    }

    /* remove the triangle storage */
    tree->ntri = 0;
    FREE(tree->trip);

cleanup:
    return status;
}
