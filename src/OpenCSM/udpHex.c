/*
 ************************************************************************
 *                                                                      *
 * udpHex -- general hexahedron                                         *
 *                                                                      *
 *            Written by John Dannenhoffer@ Syracuse University         *
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

#define NUMUDPINPUTBODYS 0
#define NUMUDPARGS       6

/* set up the necessary structures (uses NUMUDPARGS) */
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define CORNERS(IUDP,I,J)  ((double *) (udps[IUDP].arg[0].val))[I+3*J]
#define UKNOTS( IUDP,I  )  ((double *) (udps[IUDP].arg[1].val))[I]
#define VKNOTS( IUDP,I  )  ((double *) (udps[IUDP].arg[2].val))[I]
#define WKNOTS( IUDP,I  )  ((double *) (udps[IUDP].arg[3].val))[I]
#define AREA(   IUDP    )  ((double *) (udps[IUDP].arg[4].val))[0]
#define VOLUME( IUDP    )  ((double *) (udps[IUDP].arg[5].val))[0]

static char  *argNames[NUMUDPARGS] = {"corners", "uknots", "vknots", "wknots", "area",    "volume", };
static int    argTypes[NUMUDPARGS] = {ATTRREAL,  ATTRREAL, ATTRREAL, ATTRREAL, -ATTRREAL, -ATTRREAL,};
static int    argIdefs[NUMUDPARGS] = {0,         0,        0,        0,        0,         0,        };
static double argDdefs[NUMUDPARGS] = {0.,        0.,       0.,       0.,       0.,        0.,       };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define MAX(A,B) (((A) < (B)) ? (B) : (A))


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     header[7], iford, oclass, mtype, nchild, *senses;
    int     nuknots, nvknots, nwknots, i, j, k, ndata;
    double  *uknots=NULL, *vknots=NULL, *wknots=NULL, *data=NULL;
    double  bounds[4], fraci, fracj, frack;
    ego     esurf, efaces[6], emodel, eref, *echilds;

    ROUTINE(udpExecute);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    for (i = 0; i < udps[0].arg[0].size/3; i++) {
        printf("corners[%3d]= %10.5f %10.5f %10.5f\n", i,
               CORNERS(0,0,i), CORNERS(0,1,i), CORNERS(0,2,i));
    }
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf("uknots[ %3d]=%10.5f\n", i, UKNOTS(0,i));
    }
    for (j = 0; j < udps[0].arg[2].size; j++) {
        printf("vknots[ %3d]=%10.5f\n", j, VKNOTS(0,j));
    }
    for (k = 0; k < udps[0].arg[3].size; k++) {
        printf("wknots[ %3d]=%10.5f\n", k, WKNOTS(0,k));
    }
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    /* check arguments */
    if (udps[0].arg[0].size != 24) {
        printf(" udpExecute: corners should contain 24 values\n");
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp();
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    for (i = 0; i < udps[numUdp].arg[0].size/3; i++) {
        printf("corners[%3d]= %10.5f %10.5f %10.5f\n", i,
               CORNERS(numUdp,0,i), CORNERS(numUdp,1,i), CORNERS(numUdp,2,i));
    }
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf("uknots[ %3d]=%10.5f\n", i, UKNOTS(numUdp,i));
    }
    for (j = 0; j < udps[numUdp].arg[2].size; j++) {
        printf("vknots[ %3d]=%10.5f\n", j, VKNOTS(numUdp,j));
    }
    for (k = 0; k < udps[numUdp].arg[3].size; k++) {
        printf("wknots[ %3d]=%10.5f\n", k, WKNOTS(numUdp,k));
    }
#endif

    /* set up uknots */
    if (udps[numUdp].arg[1].size < 2) {
        nuknots = 4;
    } else {
        nuknots = udps[numUdp].arg[1].size + 2;
    }

    MALLOC(uknots, double, nuknots);

    if (udps[numUdp].arg[1].size < 2) {
        uknots[1] = 0;
        uknots[2] = 1;
    } else {
        for (i = 0; i < nuknots-2; i++) {
            uknots[i+1] = UKNOTS(0,i);
        }
    }

    uknots[0        ] = uknots[1        ];
    uknots[nuknots-1] = uknots[nuknots-2];

    /* set up vknots */
    if (udps[numUdp].arg[2].size < 2) {
        nvknots = 4;
    } else {
        nvknots = udps[numUdp].arg[2].size + 2;
    }

    MALLOC(vknots, double, nvknots);

    if (udps[numUdp].arg[2].size < 2) {
        vknots[1] = 0;
        vknots[2] = 1;
    } else {
        for (j = 0; j < nvknots-2; j++) {
            vknots[j+1] = VKNOTS(0,j);
        }
    }

    vknots[0        ] = vknots[1        ];
    vknots[nvknots-1] = vknots[nvknots-2];

    /* set up wknots */
    if (udps[numUdp].arg[3].size < 2) {
        nwknots = 4;
    } else {
        nwknots = udps[numUdp].arg[3].size + 2;
    }

    MALLOC(wknots, double, nwknots);

    if (udps[numUdp].arg[3].size < 2) {
        wknots[1] = 0;
        wknots[2] = 1;
    } else {
        for (k = 0; k < nwknots-2; k++) {
            wknots[k+1] = WKNOTS(0,k);
        }
    }

    wknots[0        ] = wknots[1        ];
    wknots[nwknots-1] = wknots[nwknots-2];

#ifdef DEBUG
    for (i = 0; i < nuknots; i++) {
        printf("uknots[%2d]=%10.5f\n", i, uknots[i]);
    }
    for (j = 0; j < nvknots; j++) {
        printf("vknots[%2d]=%10.5f\n", j, vknots[j]);
    }
    for (k = 0; k < nwknots; k++) {
        printf("wknots[%2d]=%10.5f\n", k, wknots[k]);
    }
#endif

    ndata = MAX(MAX(nuknots, nvknots), nwknots);
    ndata = ndata * (2 + 3 * ndata);

    MALLOC(data, double, ndata);

    /*
                ^ Y,V
                |
                2----------3
               /:         /|
              / :        / |
             /  :       /  |
            6----------7   |
            |   0------|---1  --> X,U
            |  '       |  /
            | '        | /
            |'         |/
            4----------5
           /
          Z,W

    */

    /* create the 6 Faces */
    bounds[0] = 0;
    bounds[1] = 1;
    bounds[2] = 0;
    bounds[3] = 1;

    /* face 0: xmin (0,2,4,6) */
    header[0] = 0;            /* bitflag */
    header[1] = 1;            /* uDegree */
    header[2] = nvknots-2;    /* uCP */
    header[3] = nvknots;      /* uKnots */
    header[4] = 1;            /* vDegree */
    header[5] = nwknots-2;    /* vCP */
    header[6] = nwknots;      /* vKnots */

    ndata = 0;

    for (j = 0; j < nvknots; j++) {
        data[ndata++] = vknots[j];
    }
    for (k = 0; k < nwknots; k++) {
        data[ndata++] = wknots[k];
    }

    for (k = 0; k < nwknots-2; k++) {
        for (j = 0; j < nvknots-2; j++) {
            fracj = vknots[j+1] / vknots[nvknots-1];
            frack = wknots[k+1] / wknots[nwknots-1];
            data[ndata++] = (1-fracj) * (1-frack) * CORNERS(numUdp,0,0)
                          + (  fracj) * (1-frack) * CORNERS(numUdp,0,2)
                          + (1-fracj) * (  frack) * CORNERS(numUdp,0,4)
                          + (  fracj) * (  frack) * CORNERS(numUdp,0,6);
            data[ndata++] = (1-fracj) * (1-frack) * CORNERS(numUdp,1,0)
                          + (  fracj) * (1-frack) * CORNERS(numUdp,1,2)
                          + (1-fracj) * (  frack) * CORNERS(numUdp,1,4)
                          + (  fracj) * (  frack) * CORNERS(numUdp,1,6);
            data[ndata++] = (1-fracj) * (1-frack) * CORNERS(numUdp,2,0)
                          + (  fracj) * (1-frack) * CORNERS(numUdp,2,2)
                          + (1-fracj) * (  frack) * CORNERS(numUdp,2,4)
                          + (  fracj) * (  frack) * CORNERS(numUdp,2,6);
        }
    }

    status = EG_makeGeometry(context, SURFACE, BSPLINE,NULL, header, data, &esurf);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_makeFace(esurf, SREVERSE, bounds, &(efaces[0]));
    CHECK_STATUS(EG_makeFace);

    iford = 1;
    status = EG_attributeAdd(efaces[0], "_iford", ATTRINT, 1, &iford, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* face 1: xmax (1,3,5,7) */
    header[0] = 0;            /* bitflag */
    header[1] = 1;            /* uDegree */
    header[2] = nvknots-2;    /* uCP */
    header[3] = nvknots;      /* uKnots */
    header[4] = 1;            /* vDegree */
    header[5] = nwknots-2;    /* vCP */
    header[6] = nwknots;      /* vKnots */

    ndata = 0;

    for (j = 0; j < nvknots; j++) {
        data[ndata++] = vknots[j];
    }
    for (k = 0; k < nwknots; k++) {
        data[ndata++] = wknots[k];
    }

    for (k = 0; k < nwknots-2; k++) {
        for (j = 0; j < nvknots-2; j++) {
            fracj = vknots[j+1] / vknots[nvknots-1];
            frack = wknots[k+1] / wknots[nwknots-1];
            data[ndata++] = (1-fracj) * (1-frack) * CORNERS(numUdp,0,1)
                          + (  fracj) * (1-frack) * CORNERS(numUdp,0,3)
                          + (1-fracj) * (  frack) * CORNERS(numUdp,0,5)
                          + (  fracj) * (  frack) * CORNERS(numUdp,0,7);
            data[ndata++] = (1-fracj) * (1-frack) * CORNERS(numUdp,1,1)
                          + (  fracj) * (1-frack) * CORNERS(numUdp,1,3)
                          + (1-fracj) * (  frack) * CORNERS(numUdp,1,5)
                          + (  fracj) * (  frack) * CORNERS(numUdp,1,7);
            data[ndata++] = (1-fracj) * (1-frack) * CORNERS(numUdp,2,1)
                          + (  fracj) * (1-frack) * CORNERS(numUdp,2,3)
                          + (1-fracj) * (  frack) * CORNERS(numUdp,2,5)
                          + (  fracj) * (  frack) * CORNERS(numUdp,2,7);
        }
    }

    status = EG_makeGeometry(context, SURFACE, BSPLINE,NULL, header, data, &esurf);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_makeFace(esurf, SFORWARD, bounds, &(efaces[1]));
    CHECK_STATUS(EG_makeFace);

    iford = 2;
    status = EG_attributeAdd(efaces[1], "_iford", ATTRINT, 1, &iford, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* face 2: ymin (0,4,1,5) */
    header[0] = 0;            /* bitflag */
    header[1] = 1;            /* uDegree */
    header[2] = nwknots-2;    /* uCP */
    header[3] = nwknots;      /* uKnots */
    header[4] = 1;            /* vDegree */
    header[5] = nuknots-2;    /* vCP */
    header[6] = nuknots;      /* vKnots */

    ndata = 0;

    for (k = 0; k < nwknots; k++) {
        data[ndata++] = wknots[k];
    }
    for (i = 0; i < nuknots; i++) {
        data[ndata++] = uknots[i];
    }

    for (i = 0; i < nuknots-2; i++) {
        for (k = 0; k < nwknots-2; k++) {
            frack = wknots[k+1] / wknots[nwknots-1];
            fraci = uknots[i+1] / uknots[nuknots-1];
            data[ndata++] = (1-frack) * (1-fraci) * CORNERS(numUdp,0,0)
                          + (  frack) * (1-fraci) * CORNERS(numUdp,0,4)
                          + (1-frack) * (  fraci) * CORNERS(numUdp,0,1)
                          + (  frack) * (  fraci) * CORNERS(numUdp,0,5);
            data[ndata++] = (1-frack) * (1-fraci) * CORNERS(numUdp,1,0)
                          + (  frack) * (1-fraci) * CORNERS(numUdp,1,4)
                          + (1-frack) * (  fraci) * CORNERS(numUdp,1,1)
                          + (  frack) * (  fraci) * CORNERS(numUdp,1,5);
            data[ndata++] = (1-frack) * (1-fraci) * CORNERS(numUdp,2,0)
                          + (  frack) * (1-fraci) * CORNERS(numUdp,2,4)
                          + (1-frack) * (  fraci) * CORNERS(numUdp,2,1)
                          + (  frack) * (  fraci) * CORNERS(numUdp,2,5);
        }
    }

    status = EG_makeGeometry(context, SURFACE, BSPLINE,NULL, header, data, &esurf);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_makeFace(esurf, SFORWARD, bounds, &(efaces[2]));
    CHECK_STATUS(EG_makeFace);

    iford = 3;
    status = EG_attributeAdd(efaces[2], "_iford", ATTRINT, 1, &iford, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* face 3: ymax (2,6,3,7) */
    header[0] = 0;            /* bitflag */
    header[1] = 1;            /* uDegree */
    header[2] = nwknots-2;    /* uCP */
    header[3] = nwknots;      /* uKnots */
    header[4] = 1;            /* vDegree */
    header[5] = nuknots-2;    /* vCP */
    header[6] = nuknots;      /* vKnots */

    ndata = 0;

    for (k = 0; k < nwknots; k++) {
        data[ndata++] = wknots[k];
    }
    for (i = 0; i < nuknots; i++) {
        data[ndata++] = uknots[i];
    }

    for (i = 0; i < nuknots-2; i++) {
        for (k = 0; k < nwknots-2; k++) {
            frack = wknots[k+1] / wknots[nwknots-1];
            fraci = uknots[i+1] / uknots[nuknots-1];
            data[ndata++] = (1-frack) * (1-fraci) * CORNERS(numUdp,0,2)
                          + (  frack) * (1-fraci) * CORNERS(numUdp,0,6)
                          + (1-frack) * (  fraci) * CORNERS(numUdp,0,3)
                          + (  frack) * (  fraci) * CORNERS(numUdp,0,7);
            data[ndata++] = (1-frack) * (1-fraci) * CORNERS(numUdp,1,2)
                          + (  frack) * (1-fraci) * CORNERS(numUdp,1,6)
                          + (1-frack) * (  fraci) * CORNERS(numUdp,1,3)
                          + (  frack) * (  fraci) * CORNERS(numUdp,1,7);
            data[ndata++] = (1-frack) * (1-fraci) * CORNERS(numUdp,2,2)
                          + (  frack) * (1-fraci) * CORNERS(numUdp,2,6)
                          + (1-frack) * (  fraci) * CORNERS(numUdp,2,3)
                          + (  frack) * (  fraci) * CORNERS(numUdp,2,7);
        }
    }

    status = EG_makeGeometry(context, SURFACE, BSPLINE,NULL, header, data, &esurf);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_makeFace(esurf, SREVERSE, bounds, &(efaces[3]));
    CHECK_STATUS(EG_makeFace);

    iford = 4;
    status = EG_attributeAdd(efaces[3], "_iford", ATTRINT, 1, &iford, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* face 4: zmin (0,1,2,3) */
    header[0] = 0;            /* bitflag */
    header[1] = 1;            /* uDegree */
    header[2] = nuknots-2;    /* uCP */
    header[3] = nuknots;      /* uKnots */
    header[4] = 1;            /* vDegree */
    header[5] = nvknots-2;    /* vCP */
    header[6] = nvknots;      /* vKnots */

    ndata = 0;

    for (i = 0; i < nuknots; i++) {
        data[ndata++] = uknots[i];
    }
    for (j = 0; j < nvknots; j++) {
        data[ndata++] = vknots[j];
    }

    for (j = 0; j < nvknots-2; j++) {
        for (i = 0; i < nuknots-2; i++) {
            fraci = uknots[i+1] / uknots[nuknots-1];
            fracj = vknots[j+1] / vknots[nvknots-1];
            data[ndata++] = (1-fraci) * (1-fracj) * CORNERS(numUdp,0,0)
                          + (  fraci) * (1-fracj) * CORNERS(numUdp,0,1)
                          + (1-fraci) * (  fracj) * CORNERS(numUdp,0,2)
                          + (  fraci) * (  fracj) * CORNERS(numUdp,0,3);
            data[ndata++] = (1-fraci) * (1-fracj) * CORNERS(numUdp,1,0)
                          + (  fraci) * (1-fracj) * CORNERS(numUdp,1,1)
                          + (1-fraci) * (  fracj) * CORNERS(numUdp,1,2)
                          + (  fraci) * (  fracj) * CORNERS(numUdp,1,3);
            data[ndata++] = (1-fraci) * (1-fracj) * CORNERS(numUdp,2,0)
                          + (  fraci) * (1-fracj) * CORNERS(numUdp,2,1)
                          + (1-fraci) * (  fracj) * CORNERS(numUdp,2,2)
                          + (  fraci) * (  fracj) * CORNERS(numUdp,2,3);
        }
    }

    status = EG_makeGeometry(context, SURFACE, BSPLINE,NULL, header, data, &esurf);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_makeFace(esurf, SFORWARD, bounds, &(efaces[4]));
    CHECK_STATUS(EG_makeFace);

    iford = 5;
    status = EG_attributeAdd(efaces[4], "_iford", ATTRINT, 1, &iford, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* face 5: zmax (4,5,6,7) */
    header[0] = 0;            /* bitflag */
    header[1] = 1;            /* uDegree */
    header[2] = nuknots-2;    /* uCP */
    header[3] = nuknots;      /* uKnots */
    header[4] = 1;            /* vDegree */
    header[5] = nvknots-2;    /* vCP */
    header[6] = nvknots;      /* vKnots */

    ndata = 0;

    for (i = 0; i < nuknots; i++) {
        data[ndata++] = uknots[i];
    }
    for (j = 0; j < nvknots; j++) {
        data[ndata++] = vknots[j];
    }

    for (j = 0; j < nvknots-2; j++) {
        for (i = 0; i < nuknots-2; i++) {
            fraci = uknots[i+1] / uknots[nuknots-1];
            fracj = vknots[j+1] / vknots[nvknots-1];
            data[ndata++] = (1-fraci) * (1-fracj) * CORNERS(numUdp,0,4)
                          + (  fraci) * (1-fracj) * CORNERS(numUdp,0,5)
                          + (1-fraci) * (  fracj) * CORNERS(numUdp,0,6)
                          + (  fraci) * (  fracj) * CORNERS(numUdp,0,7);
            data[ndata++] = (1-fraci) * (1-fracj) * CORNERS(numUdp,1,4)
                          + (  fraci) * (1-fracj) * CORNERS(numUdp,1,5)
                          + (1-fraci) * (  fracj) * CORNERS(numUdp,1,6)
                          + (  fraci) * (  fracj) * CORNERS(numUdp,1,7);
            data[ndata++] = (1-fraci) * (1-fracj) * CORNERS(numUdp,2,4)
                          + (  fraci) * (1-fracj) * CORNERS(numUdp,2,5)
                          + (1-fraci) * (  fracj) * CORNERS(numUdp,2,6)
                          + (  fraci) * (  fracj) * CORNERS(numUdp,2,7);
        }
    }

    status = EG_makeGeometry(context, SURFACE, BSPLINE,NULL, header, data, &esurf);
    CHECK_STATUS(EG_makeGeometry);

    status = EG_makeFace(esurf, SREVERSE, bounds, &(efaces[5]));
    CHECK_STATUS(EG_makeFace);

    iford = 6;
    status = EG_attributeAdd(efaces[5], "_iford", ATTRINT, 1, &iford, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* sew the Faces together */
    status = EG_sewFaces(6, efaces, 0, 0, &emodel);
    CHECK_STATUS(EG_sewFaces);

    status = EG_getTopology(emodel, &eref, &oclass, &mtype, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    *ebody = echilds[0];

    /* set the output value(s) */
    status = EG_getMassProperties(*ebody, data);
    CHECK_STATUS(EG_getMassProperties);

    AREA(0)   = data[1];
    VOLUME(0) = data[0];

    /* remember this model (Body) */
    udps[numUdp].ebody = *ebody;
    goto cleanup;

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    FREE(uknots);
    FREE(vknots);
    FREE(wknots);
    FREE(data  );

    if (status != EGADS_SUCCESS) {
        *string = udpErrorStr(status);
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   udpSensitivity - return sensitivity derivatives for the "real" argument *
 *                                                                      *
 ************************************************************************
 */

int
udpSensitivity(ego    ebody,            /* (in)  Body pointer */
   /*@unused@*/int    npnt,             /* (in)  number of points */
   /*@unused@*/int    entType,          /* (in)  OCSM entity type */
   /*@unused@*/int    entIndex,         /* (in)  OCSM entity index (bias-1) */
   /*@unused@*/double uvs[],            /* (in)  parametric coordinates for evaluation */
   /*@unused@*/double vels[])           /* (out) velocities */
{
    int    iudp, judp;

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        return EGADS_NOTMODEL;
    }

    /* this routine is not written yet */
    return EGADS_NOLOAD;
}
