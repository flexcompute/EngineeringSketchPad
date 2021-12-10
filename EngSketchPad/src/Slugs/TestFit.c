/*
 ************************************************************************
 *                                                                      *
 * TestFit.c -- test Fitter on Edges (and maybe Faces)                  *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
 *                                                                      *
 ************************************************************************
*/

/*
 * Copyright (C) 2013/2021  John F. Dannenhoffer, III (Syracuse University)
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

#define  DEBUG  0

#ifdef GRAFIC
    #include "grafic.h"
#endif

#ifdef WIN32
    #define snprintf   _snprintf
#endif

#define STRNCPY(A, B, LEN) strncpy(A, B, LEN); A[LEN-1] = '\0';

#define  EPS06      1.0e-6
#define  MIN(A,B)   (((A) < (B)) ? (A) : (B))
#define  MAX(A,B)   (((A) < (B)) ? (B) : (A))
#define  SQR(A)     ((A) * (A))

#include "egads.h"
#include "Fitter.h"

/*
 ***********************************************************************
 *                                                                     *
 * macros                                                              *
 *                                                                     *
 ***********************************************************************
 */

/* macros for error checking */
#define ROUTINE(NAME) char routine[] = #NAME ;\
    if (routine[0] == '\0') printf("bad routine(%s)\n", routine);

#define CHECK_STATUS(X)                                                 \
    if (status < FIT_SUCCESS) {                                         \
        printf( "ERROR:: BAD STATUS = %d from %s (called from %s:%d)\n", status, #X, routine, __LINE__); \
        goto cleanup;                                                   \
    }

#define MALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR != NULL) {                                                  \
        printf("ERROR:: MALLOC overwrites for %s=%llx (called from %s:%d)\n", #PTR, (long long)PTR, routine, __LINE__); \
        status = FIT_MALLOC;                                            \
        goto cleanup;                                                   \
    }                                                                   \
    PTR = (TYPE *) malloc((SIZE) * sizeof(TYPE));                       \
    if (PTR == NULL) {                                                  \
        printf("ERROR:: MALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
        status = FIT_MALLOC;                                            \
        goto cleanup;                                                   \
    }

#define RALLOC(PTR,TYPE,SIZE)                                           \
    if (PTR == NULL) {                                                  \
        MALLOC(PTR,TYPE,SIZE);                                          \
    } else {                                                            \
       realloc_temp = realloc(PTR, (SIZE) * sizeof(TYPE));              \
       if (PTR == NULL) {                                               \
           printf("ERROR:: RALLOC PROBLEM for %s (called from %s:%d)\n", #PTR, routine, __LINE__); \
           status = FIT_MALLOC;                                         \
           goto cleanup;                                                \
       } else {                                                         \
           PTR = (TYPE *)realloc_temp;                                  \
       }                                                                \
    }

#define FREE(PTR)                                               \
    if (PTR != NULL) {                                          \
        free(PTR);                                              \
    }                                                           \
    PTR = NULL;

#define SPLINT_CHECK_FOR_NULL(X)                                        \
    if ((X) == NULL) {                                                  \
        printf("ERROR:: SPLINT found %s is NULL (called from %s:%d)\n", #X, routine, __LINE__); \
        status = -999;                                                  \
        goto cleanup;                                                   \
    }

//static void *realloc_temp=NULL;            /* used by RALLOC macro */

/*
 ***********************************************************************
 *                                                                     *
 * structures                                                          *
 *                                                                     *
 ***********************************************************************
 */

/*
 ***********************************************************************
 *                                                                     *
 * global variables                                                    *
 *                                                                     *
 ***********************************************************************
 */

static int               outLevel = 1;       /* default output level */

ego    context;

/*
 ***********************************************************************
 *                                                                     *
 * declarations for routines defined below                             *
 *                                                                     *
 ***********************************************************************
 */

/* debug GRAFIC routies */
#ifdef GRAFIC
    int           plotCurve2(int iedge, int m, double XYZcloud[], /*@null@*/double Tcloud[],
                             int n, double cp[], double normf, double dotmin, int nmin);
    static void   plotCurve_image(int*, void*, void*, void*, void*, void*,
                                  void*, void*, void*, void*, void*, float*, char*, int);
    int           plotSurface2(int iface, int m, double XYZcloud[], /*@null@*/double UVcloud[],
                               int nu, int nv, double cp[], double normf, int nmin);
    static void   plotSurface_image(int*, void*, void*, void*, void*, void*,
                                    void*, void*, void*, void*, void*, float*, char*, int);

    static int    eval1dBspline(double T, int n, double cp[], double XYZ[],
                                /*@null@*/double dXYZdT[], /*@null@*/double dXYZdP[]);
    static int    eval2dBspline(double U, double V, int nu, int nv, double P[], double XYZ[],
                                /*@null@*/double dXYZdU[], /*@null@*/double dXYZdV[], /*@null@*/double dXYZdP[]);
    static int    cubicBsplineBases(int ncp, double t, double N[], double dN[]);
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

    int       status;
    int       nnode=0, inode, nedge=0, iedge, nface=0, iface, npnt=0, ipnt, ntri=0, itri;
    int       ncp_in=7, nu, nv, nmin, numiter, bitflag=1;
    int       ip0, ip1, ip2, it0, it1, it2, ptype, pindx, i, j, ii, jj, next;
    int       *ncp=NULL;
    double    *xyz=NULL, **edgecp=NULL, **edget=NULL, *cp=NULL, *uv=NULL;
    double    smooth=1, xnode, ynode, znode, tedge, normf, maxf, dotmin;
    double    uface, vface;
    char      *casename=NULL, *filename=NULL;
    clock_t   old_time, new_time, fit1d_time=0, fit2d_time=0;

    int       header[8];
    double    data[2000], xmin, xmax, ymin, ymax, zmin, zmax, size;
    ego       *esurfs=NULL, *efaces=NULL, ebody, emodel;

    FILE      *fp, *fp_ncp, *fp_points;

#ifdef GRAFIC
    int      io_kbd=5, io_scr=6;
    char     pltitl[255];
#endif

    ROUTINE(MAIN);

    /* --------------------------------------------------------------- */

    MALLOC(casename, char, 257);
    MALLOC(filename, char, 257);

    /* get the flags and casename(s) from the command line */
    if (argc < 2) {
        printf("Proper usage: TestFit casename [ncp=%d [smooth=%f [outLevel=%d]]]\n",
               ncp_in, smooth, outLevel);
        exit(0);
    } else {
        STRNCPY(casename, argv[1], 256);
    }

    if (argc >= 3) {
        sscanf(argv[2], "%d", &ncp_in);
    }

    if (argc >= 4) {
        sscanf(argv[3], "%lf", &smooth);
    }

    if (argc >= 5) {
        sscanf(argv[4], "%d", &outLevel);
    }

    /* welcome banner */
    printf("**********************************************************\n");
    printf("*                                                        *\n");
    printf("*                   Program TestFit                      *\n");
    printf("*                                                        *\n");
    printf("*           written by John Dannenhoffer, 2021           *\n");
    printf("*                                                        *\n");
    printf("**********************************************************\n");

#ifdef GRAFIC
    /* initialize the grafics */
    sprintf(pltitl, "Program TestFit.  casename=%s", casename);
    grinit_(&io_kbd, &io_scr, pltitl, strlen(pltitl));
#endif

    /* set up an EGADS context */
    status = EG_open(&context);
    CHECK_STATUS(EG_open);

    status = EG_setOutLevel(context, outLevel);
    CHECK_STATUS(EG_setOutLevel);

    /* read the input file */
    STRNCPY(filename, casename, 257);
    strcat(filename, ".tess");

    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("File \"%s\" does not exist\n", filename);
        exit(0);
    }

    /* open the plotdata file */
    STRNCPY(filename, casename, 257);
    strcat(filename, ".points");

    fp_points = fopen(filename, "w");
    if (fp == NULL) {
        printf("File \"%s\" could not be opened\n", filename);
        exit(0);
    }

    /* read the header */
    fscanf(fp, "%d %d %d\n", &nnode, &nedge, &nface);

    /* set up storgae for egos */
    MALLOC(esurfs, ego, nface);
    MALLOC(efaces, ego, nface);

    /* skip past the Nodes */
    for (inode = 0; inode < nnode; inode++) {
        fscanf(fp, "%lf %lf %lf\n", &xnode, &ynode, &znode);
    }

    /* create structure for holding Edge control points (and enough
       room for two degenerate Edges) */
    MALLOC(ncp,    int,     nedge+2);
    MALLOC(edgecp, double*, nedge+2);
    MALLOC(edget,  double*, nedge+2);

    for (iedge = 0; iedge < nedge+2; iedge++) {
        edgecp[iedge] = NULL;
        edget[ iedge] = NULL;
    }

    /* open .ncp file if ncp_in was negative */
    if (ncp_in < 0) {
        STRNCPY(filename, casename, 257);
        strcat(filename, ".ncp");

        fp_ncp = fopen(filename, "r");
        if (fp_ncp == NULL) {
            printf("File \"%s\" could not be opened\n", filename);
            exit(0);
        }
    } else {
        fp_ncp = NULL;
    }

    /* read the Edges and process them */
    for (iedge = 0; iedge < nedge; iedge++) {
        printf("Processing Edge %d (of %d)\n\n", iedge+1, nedge);

        if (ncp_in > 0) {
            ncp[iedge] = ncp_in;
        } else if (ncp_in < 0) {
            fscanf(fp_ncp, "%d", &(ncp[iedge]));
        } else {
            printf("Enter ncp for Edge %d:\n", iedge+1);
            scanf("%d", &(ncp[iedge]));
        }

        fscanf(fp, "%d", &npnt);

#ifndef __clang_analyzer__
        MALLOC(edgecp[iedge], double, 3*ncp[iedge]);
        MALLOC(edget[ iedge], double, npnt);

        MALLOC(xyz, double, 3*npnt);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            fscanf(fp, "%lf %lf %lf %lf", &(xyz[3*ipnt]), &(xyz[3*ipnt+1]), &(xyz[3*ipnt+2]), &tedge);
        }

        edgecp[iedge][             0] = xyz[       0];
        edgecp[iedge][             1] = xyz[       1];
        edgecp[iedge][             2] = xyz[       2];
        edgecp[iedge][3*ncp[iedge]-3] = xyz[3*npnt-3];
        edgecp[iedge][3*ncp[iedge]-2] = xyz[3*npnt-2];
        edgecp[iedge][3*ncp[iedge]-1] = xyz[3*npnt-1];

        old_time = clock();
        bitflag = 1;
        numiter = 0;
        status = fit1dCloud(npnt-2, bitflag, &(xyz[3]), ncp[iedge], edgecp[iedge], smooth, edget[iedge],
                            &normf, &maxf, &dotmin, &nmin, &numiter, stdout);
        new_time = clock();
        fit1d_time += (new_time - old_time);
        printf("fit1dCloud -> status=%d, normf=%10.3e, maxf=%10.3e, dotmin=%.3f, nmin=%d, CPU=%.4f\n\n",
               status, normf, maxf, dotmin, nmin, (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

#ifdef GRAFIC
        if (status == FIT_SUCCESS) {
            status = plotCurve2(iedge+1, npnt-2, &(xyz[3]), edget[iedge], ncp[iedge], edgecp[iedge], normf, dotmin, nmin);
            if (outLevel > 0) printf("plotCurve(iedge=%d) -> status=%d\n", iedge, status);
            CHECK_STATUS(plotCurve);
        }
#endif
#endif

        FREE(xyz);
    }

    /* number of control points in degenerate Edges */
    ncp[nedge  ] = 7;
    fscanf(fp_ncp, "%d", &(ncp[nedge]));
    ncp[nedge+1] = ncp[nedge];

    if (fp_ncp != NULL) fclose(fp_ncp);

    /* read the Faces and process them */
    for (iface = 0; iface < nface; iface++) {
        int iedge0 = 0;
        int iedge1 = 0;
        int iedge2 = 0;
        int iedge3 = 0;
        int iedge4 = 0;
        int ncloud = 0;

        fscanf(fp, "%d %d", &npnt, &ntri);

        fprintf(fp_points, "%5d %5d Points_%d\n", npnt, 0, iface+1);

        MALLOC(xyz, double, 3*npnt);
        MALLOC(uv,  double, 2*npnt);

        /* read cloud of points and determine if there are 4 surrounding Edges */
        assert (npnt > 0);

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            fscanf(fp, "%lf %lf %lf %lf %lf %d %d",
                   &(xyz[3*ncloud]), &(xyz[3*ncloud+1]), &(xyz[3*ncloud+2]), &uface, &vface, &ptype, &pindx);

            fprintf(fp_points, "%15.7f %15.7f %15.7f\n", xyz[3*ncloud], xyz[3*ncloud+1], xyz[3*ncloud+2]);

            if (ptype == -1) {
                ncloud++;
            } else if (ptype != abs(iedge4)) {
                if (pindx == 2) {
                    iedge0 = iedge1;
                    iedge1 = iedge2;
                    iedge2 = iedge3;
                    iedge3 = iedge4;
                    iedge4 = ptype;
                } else if (pindx > 2) {
                    iedge0 = iedge1;
                    iedge1 = iedge2;
                    iedge2 = iedge3;
                    iedge3 = iedge4;
                    iedge4 = -ptype;
                }
            }
        }

        /* Faces with more the 4 Edges cannot be processed */
        if (iedge0 != 0) {
            printf("Face %d CANNOT be processed\n", iface);
            printf("iedge0=%2d, iedge1=%2d, iedge2=%2d, iedge3=%2d, iedge4=%2d\n",
                   iedge0, iedge1, iedge2, iedge3, iedge4);

            for (itri = 0; itri < ntri; itri++) {
                fscanf(fp, "%d %d %d %d %d %d", &ip0, &ip1, &ip2, &it0, &it1, &it2);
            }

            continue;
        }

        /* if the Face only has two Edges, create degenerate Edges now */
#ifndef __clang_analyzer__
        if (iedge1 == 0 && iedge2 == 0) {
            iedge1 = nedge + 1;
            iedge2 = iedge3;
            iedge3 = nedge + 2;

            for (iedge = nedge; iedge <= nedge+1; iedge++) {
                FREE(edgecp[iedge]);

                MALLOC(edgecp[iedge], double, 3*ncp[iedge]);

                if ((iedge4 > 0 && iedge == nedge+1) ||
                    (iedge4 < 0 && iedge == nedge  )   ) {
                    for (i = 0;  i < ncp[iedge]; i++) {
                        edgecp[iedge][3*i  ] = edgecp[abs(iedge4)-1][                     0];
                        edgecp[iedge][3*i+1] = edgecp[abs(iedge4)-1][                     1];
                        edgecp[iedge][3*i+2] = edgecp[abs(iedge4)-1][                     2];
                    }
                } else {
                    for (i = 0;  i < ncp[iedge]; i++) {
                        edgecp[iedge][3*i  ] = edgecp[abs(iedge4)-1][3*ncp[abs(iedge4)-1]-3];
                        edgecp[iedge][3*i+1] = edgecp[abs(iedge4)-1][3*ncp[abs(iedge4)-1]-2];
                        edgecp[iedge][3*i+2] = edgecp[abs(iedge4)-1][3*ncp[abs(iedge4)-1]-1];
                    }
                }
            }
        }
#endif

        if (DEBUG) {
            if (iedge1 > 0) {
                for (i = 0; i < ncp[+iedge1-1]; i++) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge1, i, edgecp[+iedge1-1][3*i], edgecp[+iedge1-1][3*i+1], edgecp[+iedge1-1][3*i+2]);
                }
            } else {
                for (i = ncp[-iedge1-1]-1; i >= 0; i--) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge1, i, edgecp[-iedge1-1][3*i], edgecp[-iedge1-1][3*i+1], edgecp[-iedge1-1][3*i+2]);
                }
            }
            printf("\n");

            if (iedge2 > 0) {
                for (i = 0; i < ncp[+iedge2-1]; i++) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge2-1, i, edgecp[+iedge2-1][3*i], edgecp[+iedge2-1][3*i+1], edgecp[+iedge2-1][3*i+2]);
                }
            } else {
                for (i = ncp[-iedge2-1]-1; i >= 0; i--) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge2-1, i, edgecp[-iedge2-1][3*i], edgecp[-iedge2-1][3*i+1], edgecp[-iedge2-1][3*i+2]);
                }
            }
            printf("\n");

            if (iedge3 > 0) {
                for (i = 0; i < ncp[+iedge3-1]; i++) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge3-1, i, edgecp[+iedge3-1][3*i], edgecp[+iedge3-1][3*i+1], edgecp[+iedge3-1][3*i+2]);
                }
            } else {
                for (i = ncp[-iedge3-1]-1; i >= 0; i--) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge3-1, i, edgecp[-iedge3-1][3*i], edgecp[-iedge3-1][3*i+1], edgecp[-iedge3-1][3*i+2]);
                }
            }
            printf("\n");

            if (iedge4 > 0) {
                for (i = 0; i < ncp[+iedge4-1]; i++) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge4-1, i, edgecp[+iedge4-1][3*i], edgecp[+iedge4-1][3*i+1], edgecp[+iedge4-1][3*i+2]);
                }
            } else {
                for (i = ncp[-iedge4-1]-1; i >= 0; i--) {
                    printf("%3d %3d %10.5f %10.5f %10.5f\n", iedge4-1, i, edgecp[-iedge4-1][3*i], edgecp[-iedge4-1][3*i+1], edgecp[-iedge4-1][3*i+2]);
                }
            }
            printf("\n");
        }

        /* process the Face */
        printf("Processing Face %d (of %d)\n\n", iface+1, nface);
        printf("iedge0=%2d, iedge1=%2d, iedge2=%2d, iedge3=%2d, iedge4=%2d\n",
               iedge0, iedge1, iedge2, iedge3, iedge4);

        xmin = xyz[0];
        xmax = xyz[0];
        ymin = xyz[1];
        ymax = xyz[1];
        zmin = xyz[2];
        zmax = xyz[2];

        for (i = 0; i < ncloud; i++) {
            if (xyz[3*i  ] < xmin) xmin = xyz[3*i  ];
            if (xyz[3*i  ] > xmax) xmax = xyz[3*i  ];
            if (xyz[3*i+1] < ymin) ymin = xyz[3*i+1];
            if (xyz[3*i+1] > ymax) ymax = xyz[3*i+1];
            if (xyz[3*i+2] < zmin) zmin = xyz[3*i+2];
            if (xyz[3*i+2] > zmax) zmax = xyz[3*i+2];
        }

        size = sqrt(  (xmax - xmin) * (xmax - xmin)
                    + (ymax - ymin) * (ymax - ymin)
                    + (zmax - zmin) * (zmax - zmin));
        printf("size=%10.4e\n", size);

        if (ncp[abs(iedge1)-1] == ncp[abs(iedge3)-1]) {
            nu = ncp[abs(iedge1)-1];
        } else {
            printf("iedge1=%d has %d points but iedge3=%d has %d points\n",
                   iedge1, ncp[abs(iedge1)-1], iedge3, ncp[abs(iedge3)-1]);
            exit(0);
        }

        if (ncp[abs(iedge2)-1] == ncp[abs(iedge4)-1]) {
            nv = ncp[abs(iedge2)-1];
        } else {
            printf("iedge2=%d has %d points but iedge4=%d has %d points\n",
                   iedge2, ncp[abs(iedge2)-1], iedge4, ncp[abs(iedge4)-1]);
            exit(0);
        }

        printf("nu=%d,  nv=%d\n", nu, nv);
        assert (nu > 0 && nv > 0);

        /* set up outline of control net */
        MALLOC(cp, double, 3*nu*nv);

#ifndef __clang_analyzer__
        if (iedge1 > 0) {
            j  = 0;
            ii = 0;
            for (i = 0; i < nu; i++) {
                cp[3*(i+j*nu)  ] = edgecp[+iedge1-1][3*ii  ];
                cp[3*(i+j*nu)+1] = edgecp[+iedge1-1][3*ii+1];
                cp[3*(i+j*nu)+2] = edgecp[+iedge1-1][3*ii+2];
                ii++;
            }
        } else {
            j  = 0;
            ii = nu - 1;
            for (i = 0; i < nu; i++) {
                cp[3*(i+j*nu)  ] = edgecp[-iedge1-1][3*ii  ];
                cp[3*(i+j*nu)+1] = edgecp[-iedge1-1][3*ii+1];
                cp[3*(i+j*nu)+2] = edgecp[-iedge1-1][3*ii+2];
                ii--;
                }
        }

        if (iedge2 > 0) {
            i  = nu - 1;
            jj = 0;
            for (j = 0; j < nv; j++) {
                cp[3*(i+j*nu)  ] = edgecp[+iedge2-1][3*jj  ];
                cp[3*(i+j*nu)+1] = edgecp[+iedge2-1][3*jj+1];
                cp[3*(i+j*nu)+2] = edgecp[+iedge2-1][3*jj+2];
                jj++;
            }
        } else {
            i  = nu - 1;
            jj = nv - 1;
            for (j = 0; j < nv; j++) {
                cp[3*(i+j*nu)  ] = edgecp[-iedge2-1][3*jj  ];
                cp[3*(i+j*nu)+1] = edgecp[-iedge2-1][3*jj+1];
                cp[3*(i+j*nu)+2] = edgecp[-iedge2-1][3*jj+2];
                jj--;
            }
        }

        if (iedge3 > 0) {
            j  = nv - 1;
            ii = nu - 1;
            for (i = 0; i < nu; i++) {
                cp[3*(i+j*nu)  ] = edgecp[+iedge3-1][3*ii  ];
                cp[3*(i+j*nu)+1] = edgecp[+iedge3-1][3*ii+1];
                cp[3*(i+j*nu)+2] = edgecp[+iedge3-1][3*ii+2];
                ii--;
            }
        } else {
            j  = nv - 1;
            ii = 0;
            for (i = 0; i < nu; i++) {
                cp[3*(i+j*nu)  ] = edgecp[-iedge3-1][3*ii  ];
                cp[3*(i+j*nu)+1] = edgecp[-iedge3-1][3*ii+1];
                cp[3*(i+j*nu)+2] = edgecp[-iedge3-1][3*ii+2];
                ii++;
            }
        }

        if (iedge4 > 0) {
            i  = 0;
            jj = nv - 1;
            for (j = 0; j < nv; j++) {
                cp[3*(i+j*nu)  ] = edgecp[+iedge4-1][3*jj  ];
                cp[3*(i+j*nu)+1] = edgecp[+iedge4-1][3*jj+1];
                cp[3*(i+j*nu)+2] = edgecp[+iedge4-1][3*jj+2];
                jj--;
            }
        } else {
            i  = 0;
            jj = 0;
            for (j = 0; j < nv; j++) {
                cp[3*(i+j*nu)  ] = edgecp[-iedge4-1][3*jj  ];
                cp[3*(i+j*nu)+1] = edgecp[-iedge4-1][3*jj+1];
                cp[3*(i+j*nu)+2] = edgecp[-iedge4-1][3*jj+2];
                jj++;
            }
        }

        /* perform the fit and plot the results */
        old_time = clock();
        bitflag = 0;
        numiter = 0;
        status = fit2dCloud(ncloud, bitflag, xyz, nu, nv, cp, smooth, uv, &normf, &maxf, &nmin, &numiter, stdout);
        new_time = clock();
        fit2d_time += (new_time - old_time);
        printf("fit2dcloud -> status=%d, normf=%10.3e, maxf=%10.3e, nmin=%d, numiter=%d, CPU=%.4f\n\n",
               status, normf, maxf, nmin, numiter, (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

#ifdef GRAFIC
        status = plotSurface2(iface+1, ncloud, xyz, uv, nu, nv, cp, normf, nmin);
        CHECK_STATUS(plotSurface2);
#endif
#endif

        for (itri = 0; itri < ntri; itri++) {
            fscanf(fp, "%d %d %d %d %d %d", &ip0, &ip1, &ip2, &it0, &it1, &it2);
        }

        /* make the surface and Face */
        header[0] = 0;         // bitflag
        header[1] = 3;         // uDegree
        header[2] = nu;        // uCp
        header[3] = nu+4;      // uKnots
        header[4] = 3;         // vDegree
        header[5] = nv;        // vCp
        header[6] = nv+4;      // vKnots

        next = 0;
        data[next++] = 0;      // uKnots
        data[next++] = 0;
        data[next++] = 0;
        for (i = 0; i <= nu-3; i++) {
            data[next++] = i;
        }
        data[next++] = nu-3;
        data[next++] = nu-3;
        data[next++] = nu-3;

        data[next++] = 0;      // vKnots
        data[next++] = 0;
        data[next++] = 0;
        for (i = 0; i <= nv-3; i++) {
            data[next++] = i;
        }
        data[next++] = nv-3;
        data[next++] = nv-3;
        data[next++] = nv-3;

#ifndef __clang_analyzer__
        for (i = 0; i < 3*nu*nv; i++) {
            data[next++] = cp[i];
        }
#endif

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, data, &esurfs[iface]);
        CHECK_STATUS(EG_makeGeometry);

        data[0] = 0;
        data[1] = nu-3;
        data[2] = 0;
        data[3] = nv-3;

        status = EG_makeFace(esurfs[iface], SFORWARD, data, &efaces[iface]);
        CHECK_STATUS(eg_makeface);

        FREE(xyz);
        FREE(uv );
        FREE(cp );
    }

    /* close file and exit */
    fclose(fp);
    fclose(fp_points);

    /* sew Faces into a model */
    if (nface == 1) {
        header[0] =  SFORWARD;
        status = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL,
                                 1, efaces, header, &ebody);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, MODEL, 0, NULL,
                                 1, &ebody, NULL, &emodel);
        CHECK_STATUS(EG_makeTopology);
    } else {
        status = EG_sewFaces(nface, efaces, 0, 0, &emodel);
        CHECK_STATUS(EG_sewFaces);
    }

    STRNCPY(filename, casename, 257);
    strcat(filename, ".egads");

    remove(filename);

    status = EG_saveModel(emodel, filename);
    CHECK_STATUS(EG_saveModel);

    status = EG_close(context);
    CHECK_STATUS(EG_close);

    printf("Totl 1D fitting time = %10.4f sec\n", (double)(fit1d_time) / (double)(CLOCKS_PER_SEC));
    printf("Totl 2D fitting time = %10.4f sec\n", (double)(fit2d_time) / (double)(CLOCKS_PER_SEC));
    printf("==> TestFit completed successfully (for \"%s\")\n", casename);

cleanup:
    FREE(casename);
    FREE(filename);

    SPLINT_CHECK_FOR_NULL(edgecp);
    SPLINT_CHECK_FOR_NULL(edget );
    
    for (iedge = 0; iedge < nedge; iedge++) {
        FREE(edgecp[iedge]);
        FREE(edget[ iedge]);
    }
    FREE(edgecp);
    FREE(edget );
    FREE(ncp   );
    FREE(cp    );
    FREE(uv    );

    FREE(xyz);

    FREE(esurfs);
    FREE(efaces);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotCurve2 - plot curve data                                       *
 *                                                                      *
 ************************************************************************
 */
#ifdef GRAFIC
int
plotCurve2(int    iedge,                /* (in)  Edge number */
           int    m,                    /* (in)  number of points in cloud */
           double XYZcloud[],           /* (in)  array  of points in cloud */
 /*@null@*/double Tcloud[],             /* (in)  T-parameters of points in cloud (or NULL) */
           int    n,                    /* (in)  number of control points */
           double cp[],                 /* (in)  array  of control points */
           double normf,                /* (in)  RMS of distances between cloud and fit */
           double dotmin,               /* (in)  minimum normalized dot product of control polygon */
           int    nmin)                 /* (in)  minimum number of control points in any interval */
{
    int    status = FIT_SUCCESS;         /* (out) return status */

    int    indgr=1+2+4+16+64+1024, itype=0;
    char   pltitl[255];

    ROUTINE(plotCurve2);

    /* --------------------------------------------------------------- */

    sprintf(pltitl, "~x~y~ m=%d,  n=%d,  normf=%10.3e,  dotmin=%.4f,  nmin=%d",
            m, n, normf, dotmin, nmin);

    grctrl_(plotCurve_image, &indgr, pltitl,
            (void*)(&iedge),
            (void*)(&itype),
            (void*)(&m),
            (void*)(XYZcloud),
            (void*)(Tcloud),
            (void*)(&n),
            (void*)(cp),
            (void*)NULL,
            (void*)NULL,
            (void*)NULL,
            strlen(pltitl));

//cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotCurve_image - plot curve data (level 3)                        *
 *                                                                      *
 ************************************************************************
 */
static void
plotCurve_image(int   *ifunct,
                void  *iedgeP,
                void  *itypeP,
                void  *mP,
                void  *XYZcloudP,
                void  *TcloudP,
                void  *nP,
                void  *cpP,
    /*@unused@*/void  *a7,
    /*@unused@*/void  *a8,
    /*@unused@*/void  *a9,
                float *scale,
                char  *text,
                int   textlen)
{
    int    *iedge    = (int    *)iedgeP;
    int    *itype    = (int    *)itypeP;
    int    *m        = (int    *)mP;
    double *XYZcloud = (double *)XYZcloudP;
    double *Tcloud   = (double *)TcloudP;
    int    *n        = (int    *)nP;
    double *cp       = (double *)cpP;

    int    j, k, status;
    float  x4, y4, z4;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double TT, XYZ[3];
    char   lablgr[80];

    int    icircle = GR_CIRCLE;
    int    istar   = GR_STAR;
    int    isolid  = GR_SOLID;
    int    idotted = GR_DOTTED;
    int    igreen  = GR_GREEN;
    int    iblue   = GR_BLUE;
    int    ired    = GR_RED;
    int    iorange = GR_ORANGE;
    int    iblack  = GR_BLACK;

    ROUTINE(plot1dBspine_image);

    /* --------------------------------------------------------------- */

    if (*ifunct == 0) {
        xmin = XYZcloud[0];
        xmax = XYZcloud[0];
        ymin = XYZcloud[1];
        ymax = XYZcloud[1];
        zmin = XYZcloud[2];
        zmax = XYZcloud[2];

        for (k = 1; k < *m; k++) {
            if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
            if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
            if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
            if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
            if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
            if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
        }

        for (j = 0; j < *n; j++) {
            if (cp[3*j  ] < xmin) xmin = cp[3*j  ];
            if (cp[3*j  ] > xmax) xmax = cp[3*j  ];
            if (cp[3*j+1] < ymin) ymin = cp[3*j+1];
            if (cp[3*j+1] > ymax) ymax = cp[3*j+1];
            if (cp[3*j+2] < zmin) zmin = cp[3*j+2];
            if (cp[3*j+2] > zmax) zmax = cp[3*j+2];
        }

        if        (xmax-xmin >= zmax-zmin && ymax-ymin >= zmax-zmin) {
            *itype = 0;
            scale[0] = xmin - EPS06;
            scale[1] = xmax + EPS06;
            scale[2] = ymin - EPS06;
            scale[3] = ymax + EPS06;

            int iset=1;
            float zero=0;
            sprintf(lablgr, "~x~y~Edge %d", *iedge);
            grvalu_("LABLGR", &iset, &zero, lablgr, strlen("LABLGR"), strlen(lablgr));
        } else if (ymax-ymin >= xmax-xmin && zmax-zmin >= xmax-xmin) {
            *itype = 1;
            scale[0] = ymin - EPS06;
            scale[1] = ymax + EPS06;
            scale[2] = zmin - EPS06;
            scale[3] = zmax + EPS06;

            int iset=1;
            float zero=0;
            sprintf(lablgr, "~y~z~Edge %d", *iedge);
            grvalu_("LABLGR", &iset, &zero, lablgr, strlen("LABLGR"), strlen(lablgr));
        } else {
            *itype = 2;
            scale[0] = zmin - EPS06;
            scale[1] = zmax + EPS06;
            scale[2] = xmin - EPS06;
            scale[3] = xmax + EPS06;

            int iset=1;
            float zero=0;
            sprintf(lablgr, "~z~x~Edge %d", *iedge);
            grvalu_("LABLGR", &iset, &zero, lablgr, strlen("LABLGR"), strlen(lablgr));
        }

        strcpy(text, " ");

    } else if (*ifunct == 1) {

        /* cloud of points */
        grcolr_(&igreen);

        for (k = 0; k < *m; k++) {
            x4 = XYZcloud[3*k  ];
            y4 = XYZcloud[3*k+1];
            z4 = XYZcloud[3*k+2];
            if (*itype == 0) {
                grmov3_(&x4, &y4, &z4);
            } else if (*itype == 1) {
                grmov3_(&y4, &z4, &x4);
            } else {
                grmov3_(&z4, &x4, &y4);
            }
            grsymb_(&icircle);
        }

        /* control points */
        grcolr_(&iblue);
        grdash_(&idotted);

        x4 = cp[0];
        y4 = cp[1];
        z4 = cp[2];
        if (*itype == 0) {
            grmov3_(&x4, &y4, &z4);
        } else if (*itype == 1) {
            grmov3_(&y4, &z4, &x4);
        } else {
            grmov3_(&z4, &x4, &y4);
        }
        grsymb_(&istar);

        for (j = 1; j < *n; j++) {
            x4 = cp[3*j  ];
            y4 = cp[3*j+1];
            z4 = cp[3*j+2];
            if (*itype == 0) {
                grdrw3_(&x4, &y4, &z4);
            } else if (*itype == 1) {
                grdrw3_(&y4, &z4, &x4);
            } else {
                grdrw3_(&z4, &x4, &y4);
            }
            grsymb_(&istar);
        }

        /* Bspline curve */
        grcolr_(&iblack);
        grdash_(&isolid);

        x4 = cp[0];
        y4 = cp[1];
        z4 = cp[2];
        if (*itype == 0) {
            grmov3_(&x4, &y4, &z4);
        } else if (*itype == 1) {
            grmov3_(&y4, &z4, &x4);
        } else {
            grmov3_(&z4, &x4, &y4);
        }

        for (j = 1; j < 201; j++) {
            TT = (*n-3) * (double)(j) / (double)(201-1);

            status = eval1dBspline(TT, *n, cp, XYZ, NULL, NULL);
            CHECK_STATUS(eval1dBspline);

            x4 = XYZ[0];
            y4 = XYZ[1];
            z4 = XYZ[2];
            if (*itype == 0) {
                grdrw3_(&x4, &y4, &z4);
            } else if (*itype == 1) {
                grdrw3_(&y4, &z4, &x4);
            } else {
                grdrw3_(&z4, &x4, &y4);
            }
        }

        /* distance from point to curve */
        if (Tcloud != NULL) {
            grcolr_(&ired);

            for (k = 0; k < *m; k++) {
                x4 = XYZcloud[3*k  ];
                y4 = XYZcloud[3*k+1];
                z4 = XYZcloud[3*k+2];
                if (*itype == 0) {
                    grmov3_(&x4, &y4, &z4);
                } else if (*itype == 1) {
                    grmov3_(&y4, &z4, &x4);
                } else {
                    grmov3_(&z4, &x4, &y4);
                }

                status = eval1dBspline(Tcloud[k], *n, cp, XYZ, NULL, NULL);
                CHECK_STATUS(eval1dBspline);

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                if (*itype == 0) {
                    grdrw3_(&x4, &y4, &z4);
                } else if (*itype == 1) {
                    grdrw3_(&y4, &z4, &x4);
                } else {
                    grdrw3_(&z4, &x4, &y4);
                }
            }
        }

        /* distance of internal control points from average */
        grcolr_(&iorange);

        for (j = 1; j < *n-1; j++) {
            x4 = (cp[3*(j-1)  ] + cp[3*(j+1)  ]) / 2;
            y4 = (cp[3*(j-1)+1] + cp[3*(j+1)+1]) / 2;
            z4 = (cp[3*(j-1)+2] + cp[3*(j+1)+2]) / 2;
            grmov3_(&x4, &y4, &z4);

            x4 = cp[3*j  ];
            y4 = cp[3*j+1];
            z4 = cp[3*j+2];
            grdrw3_(&x4, &y4, &z4);
        }

        grcolr_(&iblack);

    } else {
        printf("ERROR:: illegal option\n");
    }

cleanup:
    return;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotSurface2 - plot surface data                                   *
 *                                                                      *
 ************************************************************************
 */
int
plotSurface2(int    iface,              /* (in)  Face nummber */
             int    m,                  /* (in)  number of points in cloud */
             double XYZcloud[],         /* (in)  array  of points in cloud */
   /*@null@*/double UVcloud[],          /* (in)  UV-parameters of points in cloud (or NULL) */
             int    nu,                 /* (in)  number of control points in U direction */
             int    nv,                 /* (in)  number of control points in V direction */
             double cp[],               /* (in)  array  of control points */
             double normf,              /* (in)  RMS of distances between cloud and fit */
             int    nmin)               /* (in)  minimum number of cloud points in any interval */
{
    int    status = FIT_SUCCESS;        /* (out) return status */

    int    indgr=1+2+4+16+64+1024;
    char   pltitl[255];

    ROUTINE(plotSurface2);

    /* --------------------------------------------------------------- */

    sprintf(pltitl, "~x~y~Face %d: m=%d,  nu=%d,  nv=%d,  normf=%10.3e,  nmin=%d",
            iface, m, nu, nv, normf, nmin);

    grctrl_(plotSurface_image, &indgr, pltitl,
            (void*)(&iface),
            (void*)(&m),
            (void*)(XYZcloud),
            (void*)(UVcloud),
            (void*)(&nu),
            (void*)(&nv),
            (void*)(cp),
            (void*)NULL,
            (void*)NULL,
            (void*)NULL,
            strlen(pltitl));

//cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotSurface_image - plot surface data (level 3)                    *
 *                                                                      *
 ************************************************************************
 */
static void
plotSurface_image(int    *ifunct,
                  void   *ifaceP,
                  void   *mP,
                  void   *XYZcloudP,
                  void   *UVcloudP,
                  void   *nuP,
                  void   *nvP,
                  void   *cpP,
        /*@null@*/void   *a7,
        /*@null@*/void   *a8,
        /*@null@*/void   *a9,
                  float  *scale,
                  char   *text,
                  int    textlen)
{
    int      *iface    = (int    *) ifaceP;
    int      *m        = (int    *) mP;
    double   *XYZcloud = (double *) XYZcloudP;
    double   *UVcloud  = (double *) UVcloudP;
    int      *nu       = (int    *) nuP;
    int      *nv       = (int    *) nvP;
    double   *cp       = (double *) cpP;

    int      status = FIT_SUCCESS;

    int      i, j, k;
    float    x4, y4, z4;
    double   xmin, xmax, ymin, ymax, UV[2], XYZ[3];

    int      icircle  = GR_CIRCLE;

    int      ired     = GR_RED;
    int      igreen   = GR_GREEN;
    int      iblue    = GR_BLUE;
    int      iyellow  = GR_YELLOW;
    int      iorange  = GR_ORANGE;
    int      iblack   = GR_BLACK;

    int      isolid   = GR_SOLID;
    int      idotted  = GR_DOTTED;

    ROUTINE(plotSurface_image);

    /* --------------------------------------------------------------- */

    /* return scales */
    if (*ifunct == 0) {
        xmin = XYZcloud[0];
        xmax = XYZcloud[0];
        ymin = XYZcloud[1];
        ymax = XYZcloud[1];

        for (k = 0; k < *m; k++) {
            if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
            if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
            if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
            if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
        }

        k = 0;
        for (j = 0; j < *nv; j++) {
            for (i = 0; i < *nu; i++) {
                if (cp[3*k  ] < xmin) xmin = cp[3*k  ];
                if (cp[3*k  ] > xmax) xmax = cp[3*k  ];
                if (cp[3*k+1] < ymin) ymin = cp[3*k+1];
                if (cp[3*k+1] > ymax) ymax = cp[3*k+1];

                k++;
            }
        }

        scale[0] = xmin;
        scale[1] = xmax;
        scale[2] = ymin;
        scale[3] = ymax;

        text[0] = '\0';


    /* plot image */
    } else if (*ifunct == 1) {

        /* points in cloud */
        grcolr_(&igreen);
        for (k = 0; k < *m; k++) {
            x4 = XYZcloud[3*k  ];
            y4 = XYZcloud[3*k+1];
            z4 = XYZcloud[3*k+2];

            grmov3_(&x4, &y4, &z4);
            grsymb_(&icircle);
        }

        /* Bspline surface */
        grcolr_(&iyellow);

        for (j = 0; j < 21; j++) {
            UV[1] = (*nv-3) * (double)(j) / (double)(21-1);

            status= eval2dBspline(0, UV[1], *nu, *nv, cp, XYZ, NULL, NULL, NULL);
            CHECK_STATUS(eval2dBspline);

            x4 = XYZ[0];
            y4 = XYZ[1];
            z4 = XYZ[2];
            grmov3_(&x4, &y4, &z4);

            for (i = 1; i < 21; i++) {
                UV[0] = (*nu-3) * (double)(i) / (double)(21-1);

                status= eval2dBspline(UV[0], UV[1], *nu, *nv, cp, XYZ, NULL, NULL, NULL);
                CHECK_STATUS(eval2dBspline);

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                grdrw3_(&x4, &y4, &z4);
            }
        }

        for (i = 0; i < 21; i++) {
            UV[0] = (*nu-3) * (double)(i) / (double)(21-1);

            status= eval2dBspline(UV[0], 0, *nu, *nv, cp, XYZ, NULL, NULL, NULL);
            CHECK_STATUS(eval2dBspline);

            x4 = XYZ[0];
            y4 = XYZ[1];
            z4 = XYZ[2];
            grmov3_(&x4, &y4, &z4);

            for (j = 1; j < 21; j++) {
                UV[1] = (*nv-3) * (double)(j) / (double)(21-1);

                status= eval2dBspline(UV[0], UV[1], *nu, *nv, cp, XYZ, NULL, NULL, NULL);
                CHECK_STATUS(eval2dBspline);

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                grdrw3_(&x4, &y4, &z4);
            }
        }

        /* draw the control net */
        grcolr_(&iblue);
        grdash_(&idotted);

        for (j = 0; j < *nv; j++) {

            x4 = cp[3*((0)+(*nu)*(j))  ];
            y4 = cp[3*((0)+(*nu)*(j))+1];
            z4 = cp[3*((0)+(*nu)*(j))+2];
            grmov3_(&x4, &y4, &z4);

            for (i = 1; i < *nu; i++) {
                x4 = cp[3*((i)+(*nu)*(j))  ];
                y4 = cp[3*((i)+(*nu)*(j))+1];
                z4 = cp[3*((i)+(*nu)*(j))+2];
                grdrw3_(&x4, &y4, &z4);
            }
        }

        for (i = 0; i < *nu; i++) {
            grcolr_(&iblack);

            x4 = cp[3*((i)+(*nu)*(0))  ];
            y4 = cp[3*((i)+(*nu)*(0))+1];
            z4 = cp[3*((i)+(*nu)*(0))+2];
            grmov3_(&x4, &y4, &z4);

            for (j = 1; j < *nv; j++) {
                x4 = cp[3*((i)+(*nu)*(j))  ];
                y4 = cp[3*((i)+(*nu)*(j))+1];
                z4 = cp[3*((i)+(*nu)*(j))+2];
                grdrw3_(&x4, &y4, &z4);
            }
        }
        grdash_(&isolid);

        /* distance from point to surface */
        if (UVcloud != NULL) {
            grcolr_(&ired);

            for (k = 0; k < *m; k++) {
                x4 = XYZcloud[3*k  ];
                y4 = XYZcloud[3*k+1];
                z4 = XYZcloud[3*k+2];
                grmov3_(&x4, &y4, &z4);

                status = eval2dBspline(UVcloud[2*k], UVcloud[2*k+1], *nu, *nv, cp, XYZ, NULL, NULL, NULL);
                CHECK_STATUS(eval2dBspline);

                x4 = XYZ[0];
                y4 = XYZ[1];
                z4 = XYZ[2];
                grdrw3_(&x4, &y4, &z4);
            }
        }

        /* distance from interior control points to their average */
        grcolr_(&iorange);

        for (j = 1; j < *nv-1; j++) {
            for (i = 1; i < *nu-1; i++) {
                x4 = (cp[3*((i-1)+(*nu)*(j  ))  ] + cp[3*((i+1)+(*nu)*(j  ))  ]
                     +cp[3*((i  )+(*nu)*(j-1))  ] + cp[3*((i  )+(*nu)*(j+1))  ]) / 2
                    -(cp[3*((i-1)+(*nu)*(j-1))  ] + cp[3*((i+1)+(*nu)*(j-1))  ]
                     +cp[3*((i-1)+(*nu)*(j+1))  ] + cp[3*((i+1)+(*nu)*(j+1))  ]) / 4;
                y4 = (cp[3*((i-1)+(*nu)*(j  ))+1] + cp[3*((i+1)+(*nu)*(j  ))+1]
                     +cp[3*((i  )+(*nu)*(j-1))+1] + cp[3*((i  )+(*nu)*(j+1))+1]) / 2
                    -(cp[3*((i-1)+(*nu)*(j-1))+1] + cp[3*((i+1)+(*nu)*(j-1))+1]
                     +cp[3*((i-1)+(*nu)*(j+1))+1] + cp[3*((i+1)+(*nu)*(j+1))+1])/4;
                z4 = (cp[3*((i-1)+(*nu)*(j  ))+2] + cp[3*((i+1)+(*nu)*(j  ))+2]
                     +cp[3*((i  )+(*nu)*(j-1))+2] + cp[3*((i  )+(*nu)*(j+1))+2]) / 2
                    -(cp[3*((i-1)+(*nu)*(j-1))+2] + cp[3*((i+1)+(*nu)*(j-1))+2]
                     +cp[3*((i-1)+(*nu)*(j+1))+2] + cp[3*((i+1)+(*nu)*(j+1))+2]) / 4;
                grmov3_(&x4, &y4, &z4);

                x4 = cp[3*(i+(*nu)*j)  ];
                y4 = cp[3*(i+(*nu)*j)+1];
                z4 = cp[3*(i+(*nu)*j)+2];
                grdrw3_(&x4, &y4, &z4);
            }
        }

        grcolr_(&iblack);
    } else {
        printf("ERROR:: illegal option (iface=%d)\n", *iface);
    }

cleanup:
    return;
}


/*
 ************************************************************************
 *                                                                      *
 *   eval1dBspline - evaluate cubic Bspline and its derivatives         *
 *                                                                      *
 ************************************************************************
 */
static int
eval1dBspline(double T,                 /* (in)  independent variable */
              int    n,                 /* (in)  number of control points */
              double cp[],              /* (in)  array  of control points */
              double XYZ[],             /* (out) dependent variables */
    /*@null@*/double dXYZdT[],          /* (out) derivative wrt T (or NULL) */
    /*@null@*/double dXYZdP[])          /* (out) derivative wrt P (or NULL) */
{
    int    status = FIT_SUCCESS;        /* (out) return status */

    int    i, span;
    double N[4], dN[4];

    ROUTINE(eval1dBspline);

    /* --------------------------------------------------------------- */

    assert (n > 3);

    XYZ[0] = 0;
    XYZ[1] = 0;
    XYZ[2] = 0;

    /* set up the Bspline bases */
    status = cubicBsplineBases(n, T, N, dN);
    CHECK_STATUS(cubicBsplineBases);

    span = MIN(floor(T), n-4);

    /* find the dependent variable */
    for (i = 0; i < 4; i++) {
        XYZ[0] += N[i] * cp[3*(i+span)  ];
        XYZ[1] += N[i] * cp[3*(i+span)+1];
        XYZ[2] += N[i] * cp[3*(i+span)+2];
    }

    /* find the deriviative wrt T */
    if (dXYZdT != NULL) {
        dXYZdT[0] = 0;
        dXYZdT[1] = 0;
        dXYZdT[2] = 0;

        for (i = 0; i < 4; i++) {
            dXYZdT[0] += dN[i] * cp[3*(i+span)  ];
            dXYZdT[1] += dN[i] * cp[3*(i+span)+1];
            dXYZdT[2] += dN[i] * cp[3*(i+span)+2];
        }
    }

    /* find the derivative wrt P */
    if (dXYZdP != NULL) {
        for (i = 0; i < n; i++) {
            dXYZdP[i] = 0;
        }
        for (i = 0; i < 4; i++) {
            dXYZdP[i+span] = N[i];
        }
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   eval2dBspline - evaluate Spline and its derivative                 *
 *                                                                      *
 ************************************************************************
 */
static int
eval2dBspline(double U,                 /* (in)  first  independent variable */
              double V,                 /* (in)  second independent variable */
              int    nu,                /* (in)  number of U control points */
              int    nv,                /* (in)  number of V control points */
              double cp[],              /* (in)  array  of control points */
              double XYZ[],             /* (out) dependent variable */
    /*@null@*/double dXYZdU[],          /* (out) derivative wrt U (or NULL) */
    /*@null@*/double dXYZdV[],          /* (out) derivative wrt V (or NULL) */
    /*@null@*/double dXYZdP[])          /* (out) derivative wrt P (or NULL) */
{
    int    status = FIT_SUCCESS;         /* (out) return status */

    int    i, j, spanu, spanv;
    double Nu[4], dNu[4], Nv[4], dNv[4];

    ROUTINE(eval2dBspline);

    /* --------------------------------------------------------------- */

    assert(nu > 3);
    assert(nv > 3);

    XYZ[0] = 0;
    XYZ[1] = 0;
    XYZ[2] = 0;

    /* set up the Bspline bases */
    status = cubicBsplineBases(nu, U, Nu, dNu);
    CHECK_STATUS(cubicBsplineBases);

    status = cubicBsplineBases(nv, V, Nv, dNv);
    CHECK_STATUS(cubicBsplineBases);

    spanu = MIN(floor(U), nu-4);
    spanv = MIN(floor(V), nv-4);

    /* find the dependent variable */
    for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) {
            XYZ[0] += Nu[i] * Nv[j] * cp[3*((i+spanu)+nu*(j+spanv))  ];
            XYZ[1] += Nu[i] * Nv[j] * cp[3*((i+spanu)+nu*(j+spanv))+1];
            XYZ[2] += Nu[i] * Nv[j] * cp[3*((i+spanu)+nu*(j+spanv))+2];
        }
    }

    /* find the derivative wrt U */
    if (dXYZdU != NULL) {
        dXYZdU[0] = 0;
        dXYZdU[1] = 0;
        dXYZdU[2] = 0;

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                dXYZdU[0] += dNu[i] * Nv[j] * cp[3*((i+spanu)+nu*(j+spanv))  ];
                dXYZdU[1] += dNu[i] * Nv[j] * cp[3*((i+spanu)+nu*(j+spanv))+1];
                dXYZdU[2] += dNu[i] * Nv[j] * cp[3*((i+spanu)+nu*(j+spanv))+2];
            }
        }
    }

    /* find the derivative wrt V */
    if (dXYZdV != NULL) {
        dXYZdV[0] = 0;
        dXYZdV[1] = 0;
        dXYZdV[2] = 0;

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                dXYZdV[0] += Nu[i] * dNv[j] * cp[3*((i+spanu)+nu*(j+spanv))  ];
                dXYZdV[1] += Nu[i] * dNv[j] * cp[3*((i+spanu)+nu*(j+spanv))+1];
                dXYZdV[2] += Nu[i] * dNv[j] * cp[3*((i+spanu)+nu*(j+spanv))+2];
            }
        }
    }

    /* find the derivative wrt P */
    if (dXYZdP != NULL) {
        for (j = 0; j < nv; j++) {
            for (i = 0; i < nu; i++) {
                dXYZdP[i+nu*j] = 0;
            }
        }

        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++) {
                dXYZdP[(i+spanu)+nu*(j+spanv)] = Nu[i] * Nv[j];
            }
        }
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   cubicBsplineBases - basis function values for cubic Bspline        *
 *                                                                      *
 ************************************************************************
 */
static int
cubicBsplineBases(int    ncp,           /* (in)  number of control points */
                  double T,             /* (in)  independent variable (0<=T<=(ncp-3) */
                  double N[],           /* (out) bases */
                  double dN[])          /* (out) d(bases)/d(T) */
{
    int       status = FIT_SUCCESS;     /* (out) return status */

    int      i, r, span;
    double   saved, dsaved, num, dnum, den, dden, temp, dtemp;
    double   left[4], dleft[4], rite[4], drite[4];

    ROUTINE(cubicBsplineBases);

    /* --------------------------------------------------------------- */

    span = MIN(floor(T)+3, ncp-1);

    N[ 0] = 1.0;
    dN[0] = 0;

    for (i = 1; i <= 3; i++) {
        left[ i] = T - MAX(0, span-2-i);
        dleft[i] = 1;

        rite[ i] = MIN(ncp-3,span-3+i) - T;
        drite[i] =                     - 1;

        saved  = 0;
        dsaved = 0;

        for (r = 0; r < i; r++) {
            num   = N[ r];
            dnum  = dN[r];

            den   = rite[ r+1] + left[ i-r];
            dden  = drite[r+1] + dleft[i-r];

            temp  = num / den;
            dtemp = (dnum * den - dden * num) / den / den;

            N[ r] = saved  + rite[ r+1] * temp;
            dN[r] = dsaved + drite[r+1] * temp + rite[r+1] * dtemp;

            saved  = left[ i-r] * temp;
            dsaved = dleft[i-r] * temp + left[i-r] * dtemp;
        }

        N[ i] = saved;
        dN[i] = dsaved;
    }

//cleanup:
    return status;
}

#endif
