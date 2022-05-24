/*
 ************************************************************************
 *                                                                      *
 * udpFreeform -- udp file to generate a freeform brick                 *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2022  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS 8
#include "udpUtilities.h"

#ifdef GRAFIC
    #include "grafic.h"
#endif

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)    ((char   *) (udps[IUDP].arg[0].val))
#define IMAX(    IUDP)    ((int    *) (udps[IUDP].arg[1].val))[0]
#define JMAX(    IUDP)    ((int    *) (udps[IUDP].arg[2].val))[0]
#define KMAX(    IUDP)    ((int    *) (udps[IUDP].arg[3].val))[0]
#define XYZ(     IUDP,I)  ((double *) (udps[IUDP].arg[4].val))[I]
#define X(       IUDP,I)  ((double *) (udps[IUDP].arg[5].val))[I]
#define Y(       IUDP,I)  ((double *) (udps[IUDP].arg[6].val))[I]
#define Z(       IUDP,I)  ((double *) (udps[IUDP].arg[7].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename", "imax",  "jmax",  "kmax",  "xyz",    "x", "y", "z", };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,   ATTRINT, ATTRINT, ATTRINT, ATTRREAL, 0,   0,   0,   };
static int    argIdefs[NUMUDPARGS] = {0,          1,       1,       1,       0,        0,   0,   0,   };
static double argDdefs[NUMUDPARGS] = {0.,         0.,      0.,      0.,      0.,       0.,  0.,  0.,  };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* routines defined below */
static int spline1d(ego context, int imax,           double *x, double *y, double *z, ego *ecurve);
static int spline2d(ego context, int imax, int jmax, double *x, double *y, double *z, ego *esurf);

#ifdef GRAFIC
static void plotData(int*, void*, void*, void*, void*, void*,
                           void*, void*, void*, void*, void*, float*, char*, int);
#endif

static void *realloc_temp=NULL;              /* used by RALLOC macro */


/*
 ************************************************************************
 *                                                                      *
 *   macros for udpExecute                                              *
 *                                                                      *
 ************************************************************************
 */

#define CREATE_NODE(INODE, I, J, K)                                     \
    if (outLevel >= 1) printf("        creating Node %3d\n", INODE);    \
                                                                        \
    /* create data */                                                   \
    data[0] = X(numUdp,(I)+imax*((J)+jmax*(K)));                        \
    data[1] = Y(numUdp,(I)+imax*((J)+jmax*(K)));                        \
    data[2] = Z(numUdp,(I)+imax*((J)+jmax*(K)));                        \
                                                                        \
    /* make Node */                                                     \
    status = EG_makeTopology(context, NULL, NODE, 0,                    \
                             data, 0, NULL, NULL, &(enodes[INODE]));    \
    CHECK_STATUS(EG_makeTopology);

#define CREATE_EDGE(IEDGE, IBEG, IEND, IMAX)                            \
    if (outLevel >= 1) printf("        creating Edge %3d\n", IEDGE);    \
                                                                        \
    /* make spline Curve */                                             \
    status = spline1d(context, IMAX, x2d, y2d, z2d, &(ecurvs[IEDGE]));  \
    if (status != EGADS_SUCCESS) goto cleanup;                          \
                                                                        \
    /* make Edge */                                                     \
    etemp[0] = enodes[IBEG];                                            \
    etemp[1] = enodes[IEND];                                            \
    tdata[0] = 0;                                                       \
    tdata[1] = IMAX-1;                                                  \
                                                                        \
    status = EG_makeTopology(context, ecurvs[IEDGE], EDGE, TWONODE,     \
                             tdata, 2, etemp, NULL, &(eedges[IEDGE]));  \
    CHECK_STATUS(EG_makeTopology);

#define CREATE_FACE(IFACE, IS, IE, IN, IW, JMAX, KMAX)                  \
    if (outLevel >= 1) printf("        creating Face %3d\n", IFACE);    \
                                                                        \
    /* make spline Surface */                                           \
    status = spline2d(context, KMAX, JMAX, x2d, y2d, z2d,               \
                      &(esurfs[IFACE]));                                \
    CHECK_STATUS(spline2d);                                             \
                                                                        \
    tdata[0] = 0; tdata[1] = 0;                                         \
    (void) EG_evaluate(esurfs[IFACE], tdata, data);                     \
                                                                        \
    tdata[0] = KMAX-1; tdata[1] = 0;                                    \
    (void) EG_evaluate(esurfs[IFACE], tdata, data);                     \
                                                                        \
    tdata[0] = 0; tdata[1] = JMAX-1;                                    \
    (void) EG_evaluate(esurfs[IFACE], tdata, data);                     \
                                                                        \
    tdata[0] = KMAX-1; tdata[1] = JMAX-1;                               \
    (void) EG_evaluate(esurfs[IFACE], tdata, data);                     \
                                                                        \
    /* remember Edges that surround the Face */                         \
    etemp[0] = eedges[IS];                                              \
    etemp[1] = eedges[IE];                                              \
    etemp[2] = eedges[IN];                                              \
    etemp[3] = eedges[IW];                                              \
                                                                        \
    /* make PCurves */                                                  \
    data[0] = 0;         data[1] = 0;                                   \
    data[2] = KMAX-1;    data[3] = 0;                                   \
    status = EG_makeGeometry(context, PCURVE, LINE, NULL,               \
                             NULL, data, &(etemp[4]));                  \
    CHECK_STATUS(EG_makeGeometry);                                      \
                                                                        \
    data[0] = KMAX-1;    data[1] = 0;                                   \
    data[2] = 0;         data[3] = JMAX-1;                              \
    status = EG_makeGeometry(context, PCURVE, LINE, NULL,               \
                             NULL, data, &(etemp[5]));                  \
    CHECK_STATUS(EG_makeGeometry);                                      \
                                                                        \
    data[0] = 0;         data[1] = JMAX-1;                              \
    data[2] = KMAX-1;    data[3] = 0;                                   \
    status = EG_makeGeometry(context, PCURVE, LINE, NULL,               \
                             NULL, data, &(etemp[6]));                  \
    CHECK_STATUS(EG_makeGeometry);                                      \
                                                                        \
    data[0] = 0;         data[1] = 0;                                   \
    data[2] = 0;         data[3] = JMAX-1;                              \
    status = EG_makeGeometry(context, PCURVE, LINE, NULL,               \
                             NULL, data , &(etemp[7]));                 \
    CHECK_STATUS(EG_makeGeometry);                                      \
                                                                        \
    /* make the Loop */                                                 \
    senses[0] = SFORWARD; senses[1] = SFORWARD;                         \
    senses[2] = SREVERSE; senses[3] = SREVERSE;                         \
    senses[4] = SFORWARD; senses[5] = SFORWARD;                         \
    senses[6] = SREVERSE; senses[7] = SREVERSE;                         \
                                                                        \
    status = EG_makeTopology(context, esurfs[IFACE], LOOP, CLOSED,      \
                             NULL, 4, etemp, senses, &eloop);           \
    CHECK_STATUS(EG_makeTopology);                                      \
                                                                        \
    /* make the Face */                                                 \
    status = EG_makeTopology(context, esurfs[IFACE], FACE, SFORWARD,    \
                             NULL, 1, &eloop, senses, &(efaces[IFACE])); \
    CHECK_STATUS(EG_makeTopology);


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

    int     imax, jmax, kmax, i, j, k, ijk, senses[8], periodic, count, outLevel, ieof;
    double  *x2d=NULL, *y2d=NULL, *z2d=NULL;
    double  xtemp, ytemp, ztemp, data[18], tdata[2];
    char    *message=NULL, *filename, *token;
    FILE    *fp=NULL;
    ego                ecurvs[12], esurfs[6];
    ego     enodes[8], eedges[12], efaces[6], etemp[8], eloop, eshell;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename(0) = %s\n", FILENAME(0));
    printf("imax(0)     = %d\n", IMAX(    0));
    printf("jmax(0)     = %d\n", JMAX(    0));
    printf("kmax(0)     = %d\n", KMAX(    0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */

    /* get the outLevel from OpenCSM */
    outLevel = ocsmSetOutLevel(-1);

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
    printf("imax(%d)     = %d\n", numUdp, IMAX(    numUdp));
    printf("jmax(%d)     = %d\n", numUdp, JMAX(    numUdp));
    printf("kmax(%d)     = %d\n", numUdp, KMAX(    numUdp));
#endif

    /* data is in a file */
    if (strlen(FILENAME(0)) > 0) {

        /* open the file or stream */
        filename = FILENAME(numUdp);

        if (strncmp(filename, "<<\n", 3) == 0) {
            token = strtok(filename, " \t\n");
            if (token == NULL) {
                ieof = 1;
            } else {
                ieof = 0;
            }
        } else {
            fp = fopen(filename, "r");
            if (fp == NULL) {
                snprintf(message, 100, "could not open file \"%s\"", filename);
                status = EGADS_NOTFOUND;
                goto cleanup;
            }

            ieof = feof(fp);
        }
        if (ieof == 1) {
            snprintf(message, 100, "premature enf-of-file found");
            status = EGADS_NODATA;
            goto cleanup;
        }

        /* read the size of the configuration */
        if (fp != NULL) {
            count = fscanf(fp, "%d %d %d", &imax, &jmax, &kmax);
            if (count != 3) {
                status = EGADS_NODATA;
                goto cleanup;
            }
        } else {
            token = strtok(NULL, " \t\n");
            if (token == NULL) {
                snprintf(message, 100, "error while reading imax");
                status = EGADS_NODATA;
                goto cleanup;
            }
            imax = strtol(token, NULL, 10);

            token = strtok(NULL, " \t\n");
            if (token == NULL) {
                snprintf(message, 100, "error while reading jmax");
                status = EGADS_NODATA;
                goto cleanup;
            }
            jmax = strtol(token, NULL, 10);

            token = strtok(NULL, " \t\n");
            if (token == NULL) {
                snprintf(message, 100, "error while reading kmax");
                status = EGADS_NODATA;
                goto cleanup;
            }
            kmax = strtol(token, NULL, 10);
        }

        /* save the array size from the file */
        IMAX(numUdp) = imax;
        JMAX(numUdp) = jmax;
        KMAX(numUdp) = kmax;

        RALLOC(udps[numUdp].arg[5].val, double, imax*jmax*kmax);
        RALLOC(udps[numUdp].arg[6].val, double, imax*jmax*kmax);
        RALLOC(udps[numUdp].arg[7].val, double, imax*jmax*kmax);

        /* read the data.  if 3d, only the outside points are read */
        for (k = 0; k < kmax; k++) {
            for (j = 0; j < jmax; j++) {
                for (i = 0; i < imax; i++) {
                    if (i == 0 || i == imax-1 ||
                        j == 0 || j == jmax-1 ||
                        k == 0 || k == kmax-1   ) {
                        ijk = (i) + imax * ((j) + jmax * (k));

                        if (fp != NULL) {
                            count = fscanf(fp, "%lf", &xtemp);
                            if (count != 1) {
                                snprintf(message, 100, "error while reading x[%d,%d,%d]", i, j, k);
                                status = EGADS_NODATA;
                                goto cleanup;
                            }
                        } else {
                            token = strtok(NULL, " \t\n");
                            if (token == NULL) {
                                snprintf(message, 100, "error while reading x[%d,%d,%d]", i, j, k);
                                status = EGADS_NODATA;
                                goto cleanup;
                            }
                            xtemp = strtod(token, NULL);
                        }
                        X(numUdp,ijk) = xtemp;

                        if (fp != NULL) {
                            count = fscanf(fp, "%lf", &ytemp);
                            if (count != 1) {
                                snprintf(message, 100, "error while reading y[%d,%d,%d]", i, j, k);
                                status = EGADS_NODATA;
                                goto cleanup;
                            }
                        } else {
                            token = strtok(NULL, " \t\n");
                            if (token == NULL) {
                                snprintf(message, 100, "error while reading y[%d,%d,%d]", i, j, k);
                                status = EGADS_NODATA;
                                goto cleanup;
                            }
                            ytemp = strtod(token, NULL);
                        }
                        Y(numUdp,ijk) = ytemp;

                        if (fp != NULL) {
                            count = fscanf(fp, "%lf", &ztemp);
                            if (count != 1) {
                                snprintf(message, 100, "error while reading z[%d,%d,%d]", i, j, k);
                                status = EGADS_NODATA;
                                goto cleanup;
                            }
                        } else {
                            token = strtok(NULL, " \t\n");
                            if (token == NULL) {
                                snprintf(message, 100, "error while reading z[%d,%d,%d]", i, j, k);
                                status = EGADS_NODATA;
                                goto cleanup;
                            }
                            ztemp = strtod(token, NULL);
                        }
                        Z(numUdp,ijk) = ztemp;
                    }
                }
            }
        }

        /* close the file */
        if (fp != NULL) {
            fclose(fp);
        }

    /* data came in XYZ */
    } else if (udps[0].arg[4].size > 4) {
        imax = IMAX(0);
        jmax = JMAX(0);
        kmax = KMAX(0);

        IMAX(numUdp) = imax;
        JMAX(numUdp) = jmax;
        KMAX(numUdp) = kmax;

        RALLOC(udps[numUdp].arg[5].val, double, imax*jmax*kmax);
        RALLOC(udps[numUdp].arg[6].val, double, imax*jmax*kmax);
        RALLOC(udps[numUdp].arg[7].val, double, imax*jmax*kmax);

        /* copy the data */
        for (ijk = 0; ijk < imax*jmax*kmax; ijk++) {
            X(numUdp,ijk) = XYZ(numUdp,3*ijk  );
            Y(numUdp,ijk) = XYZ(numUdp,3*ijk+1);
            Z(numUdp,ijk) = XYZ(numUdp,3*ijk+2);
        }

    /* no data specified */
    } else {
        snprintf(message, 100, "filename and xyz both null\n");
        return EGADS_NODATA;
    }

    /* plot the data */
#ifdef GRAFIC
    {
        int  io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
        char pltitl[80];

        sprintf(pltitl, "~x~y~ imax=%d  jmax=%d  kmax=%d", imax, jmax, kmax);

        grinit_(&io_kbd, &io_scr, "udpFreeform", strlen("udpFreeform"));
        grctrl_(plotData, &indgr, pltitl,
                (void*)(&imax),
                (void*)(&jmax),
                (void*)(&kmax),
                (void*)(udps[numUdp].arg[5].val),
                (void*)(udps[numUdp].arg[6].val),
                (void*)(udps[numUdp].arg[7].val),
                (void*)NULL,
                (void*)NULL,
                (void*)NULL,
                (void*)NULL,
                strlen(pltitl));
    }
    #endif

    /* allocate necessary (oversized) 2D arrays */
    MALLOC(x2d, double, imax*jmax*kmax);
    MALLOC(y2d, double, imax*jmax*kmax);
    MALLOC(z2d, double, imax*jmax*kmax);

    /* create WireBody (since jmax<=1) */
    if (jmax <= 1) {
        if (outLevel >= 1) printf("    WireBody: (%d)\n", imax);

        /* create the 2 Nodes */
        CREATE_NODE(0, 0,      0, 0);
        CREATE_NODE(1, imax-1, 0, 0);

        /* create the Curve and Edge */
        for (i = 0; i < imax; i++) {
            x2d[i] = X(numUdp,(i)+imax*((0)+jmax*(0)));
            y2d[i] = Y(numUdp,(i)+imax*((0)+jmax*(0)));
            z2d[i] = Z(numUdp,(i)+imax*((0)+jmax*(0)));
        }
       CREATE_EDGE(0, 0, 1, imax);

        /* make a Loop */
        senses[0] = SFORWARD;
        status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                 NULL, 1, eedges, senses, &eloop);
        CHECK_STATUS(EG_makeTopology);

        /* make a WireBody */
        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* set the output value(s) */

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

    /* create FaceBody (since kmax<=1) */
    } else if (kmax <= 1) {
        if (outLevel >= 1) printf("    FaceBody: (%d*%d)\n", imax, jmax);

        /* make the cubic spline */
        status = spline2d(context, imax, jmax, &(X(numUdp,0)), &(Y(numUdp,0)), &(Z(numUdp,0)), &(esurfs[0]));
        CHECK_STATUS(spline2d);

        /* make a Face */
        status = EG_getRange(esurfs[0], data, &periodic);
        CHECK_STATUS(EG_getRange);

        status = EG_makeFace(esurfs[0], SFORWARD, data, &(efaces[0]));
        CHECK_STATUS(EG_makeFace);

        /* make a FaceBody */
        senses[0] = SFORWARD;
        status = EG_makeTopology(context, NULL, BODY, FACEBODY,
                                 NULL, 1, efaces, senses, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* inform the user that there is 1 surface meshes */
        *nMesh = 1;

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

    /* create a SolidBody */
    } else {
        if (outLevel >= 1) printf("    SolidBody: (%d*%d*%d)\n", imax, jmax, kmax);

        /* create the 8 Nodes */
        CREATE_NODE(0, 0,      0,      0     );
        CREATE_NODE(1, imax-1, 0,      0     );
        CREATE_NODE(2, 0,      jmax-1, 0     );
        CREATE_NODE(3, imax-1, jmax-1, 0     );
        CREATE_NODE(4, 0,      0,      kmax-1);
        CREATE_NODE(5, imax-1, 0,      kmax-1);
        CREATE_NODE(6, 0,      jmax-1, kmax-1);
        CREATE_NODE(7, imax-1, jmax-1, kmax-1);

        /* create the 12 Curves and Edges */
        for (i = 0; i < imax; i++) {
            x2d[i] = X(numUdp,(i     )+imax*((0     )+jmax*(0     )));
            y2d[i] = Y(numUdp,(i     )+imax*((0     )+jmax*(0     )));
            z2d[i] = Z(numUdp,(i     )+imax*((0     )+jmax*(0     )));
        }
        CREATE_EDGE( 0, 0, 1, imax);

        for (i = 0; i < imax; i++) {
            x2d[i] = X(numUdp,(i     )+imax*((jmax-1)+jmax*(0     )));
            y2d[i] = Y(numUdp,(i     )+imax*((jmax-1)+jmax*(0     )));
            z2d[i] = Z(numUdp,(i     )+imax*((jmax-1)+jmax*(0     )));
        }
        CREATE_EDGE( 1, 2, 3, imax);

        for (i = 0; i < imax; i++) {
            x2d[i] = X(numUdp,(i     )+imax*((0     )+jmax*(kmax-1)));
            y2d[i] = Y(numUdp,(i     )+imax*((0     )+jmax*(kmax-1)));
            z2d[i] = Z(numUdp,(i     )+imax*((0     )+jmax*(kmax-1)));
        }
        CREATE_EDGE( 2, 4, 5, imax);

        for (i = 0; i < imax; i++) {
            x2d[i] = X(numUdp,(i     )+imax*((jmax-1)+jmax*(kmax-1)));
            y2d[i] = Y(numUdp,(i     )+imax*((jmax-1)+jmax*(kmax-1)));
            z2d[i] = Z(numUdp,(i     )+imax*((jmax-1)+jmax*(kmax-1)));
        }
        CREATE_EDGE( 3, 6, 7, imax);

        for (j = 0; j < jmax; j++) {
            x2d[j] = X(numUdp,(0     )+imax*((j     )+jmax*(0     )));
            y2d[j] = Y(numUdp,(0     )+imax*((j     )+jmax*(0     )));
            z2d[j] = Z(numUdp,(0     )+imax*((j     )+jmax*(0     )));
        }
        CREATE_EDGE( 4, 0, 2, jmax);

        for (j = 0; j < jmax; j++) {
            x2d[j] = X(numUdp,(0     )+imax*((j     )+jmax*(kmax-1)));
            y2d[j] = Y(numUdp,(0     )+imax*((j     )+jmax*(kmax-1)));
            z2d[j] = Z(numUdp,(0     )+imax*((j     )+jmax*(kmax-1)));
        }
        CREATE_EDGE( 5, 4, 6, jmax);

        for (j = 0; j < jmax; j++) {
            x2d[j] = X(numUdp,(imax-1)+imax*((j     )+jmax*(0     )));
            y2d[j] = Y(numUdp,(imax-1)+imax*((j     )+jmax*(0     )));
            z2d[j] = Z(numUdp,(imax-1)+imax*((j     )+jmax*(0     )));
        }
        CREATE_EDGE( 6, 1, 3, jmax);

        for (j = 0; j < jmax; j++) {
            x2d[j] = X(numUdp,(imax-1)+imax*((j     )+jmax*(kmax-1)));
            y2d[j] = Y(numUdp,(imax-1)+imax*((j     )+jmax*(kmax-1)));
            z2d[j] = Z(numUdp,(imax-1)+imax*((j     )+jmax*(kmax-1)));
        }
        CREATE_EDGE( 7, 5, 7, jmax);

        for (k = 0; k < kmax; k++) {
            x2d[k] = X(numUdp,(0     )+imax*((0     )+jmax*(k     )));
            y2d[k] = Y(numUdp,(0     )+imax*((0     )+jmax*(k     )));
            z2d[k] = Z(numUdp,(0     )+imax*((0     )+jmax*(k     )));
        }
        CREATE_EDGE( 8, 0, 4, kmax);

        for (k = 0; k < kmax; k++) {
            x2d[k] = X(numUdp,(imax-1)+imax*((0     )+jmax*(k     )));
            y2d[k] = Y(numUdp,(imax-1)+imax*((0     )+jmax*(k     )));
            z2d[k] = Z(numUdp,(imax-1)+imax*((0     )+jmax*(k     )));
        }
        CREATE_EDGE( 9, 1, 5, kmax);

        for (k = 0; k < kmax; k++) {
            x2d[k] = X(numUdp,(0     )+imax*((jmax-1)+jmax*(k     )));
            y2d[k] = Y(numUdp,(0     )+imax*((jmax-1)+jmax*(k     )));
            z2d[k] = Z(numUdp,(0     )+imax*((jmax-1)+jmax*(k     )));
        }
        CREATE_EDGE(10, 2, 6, kmax);

        for (k = 0; k < kmax; k++) {
            x2d[k] = X(numUdp,(imax-1)+imax*((jmax-1)+jmax*(k     )));
            y2d[k] = Y(numUdp,(imax-1)+imax*((jmax-1)+jmax*(k     )));
            z2d[k] = Z(numUdp,(imax-1)+imax*((jmax-1)+jmax*(k     )));
        }
        CREATE_EDGE(11, 3, 7, kmax);

        /* create the 6 Surfaces, Pcurves, Loops, and Edges */
        for (j = 0; j < jmax; j++) {
            for (k = 0; k < kmax; k++) {
                x2d[k+j*kmax] = X(numUdp,(0     )+imax*((j     )+jmax*(k     )));
                y2d[k+j*kmax] = Y(numUdp,(0     )+imax*((j     )+jmax*(k     )));
                z2d[k+j*kmax] = Z(numUdp,(0     )+imax*((j     )+jmax*(k     )));
            }
        }
        CREATE_FACE(0,  8,  5, 10,  4, jmax, kmax);

        for (k = 0; k < kmax; k++) {
            for (j = 0; j < jmax; j++) {
                x2d[j+k*jmax] = X(numUdp,(imax-1)+imax*((j     )+jmax*(k     )));
                y2d[j+k*jmax] = Y(numUdp,(imax-1)+imax*((j     )+jmax*(k     )));
                z2d[j+k*jmax] = Z(numUdp,(imax-1)+imax*((j     )+jmax*(k     )));
            }
        }
        CREATE_FACE(1,  6, 11,  7,  9, kmax, jmax);

        for (k = 0; k < kmax; k++) {
            for (i = 0; i < imax; i++) {
                x2d[i+k*imax] = X(numUdp,(i     )+imax*((0     )+jmax*(k     )));
                y2d[i+k*imax] = Y(numUdp,(i     )+imax*((0     )+jmax*(k     )));
                z2d[i+k*imax] = Z(numUdp,(i     )+imax*((0     )+jmax*(k     )));
            }
        }
        CREATE_FACE(2,  0,  9,  2,  8, kmax, imax);

        for (i = 0; i < imax; i++) {
            for (k = 0; k < kmax; k++) {
                x2d[k+i*kmax] = X(numUdp,(i     )+imax*((jmax-1)+jmax*(k     )));
                y2d[k+i*kmax] = Y(numUdp,(i     )+imax*((jmax-1)+jmax*(k     )));
                z2d[k+i*kmax] = Z(numUdp,(i     )+imax*((jmax-1)+jmax*(k     )));
            }
        }
        CREATE_FACE(3, 10,  3, 11,  1, imax, kmax);

        for (i = 0; i < imax; i++) {
            for (j = 0; j < jmax; j++) {
                x2d[j+i*jmax] = X(numUdp,(i     )+imax*((j     )+jmax*(0     )));
                y2d[j+i*jmax] = Y(numUdp,(i     )+imax*((j     )+jmax*(0     )));
                z2d[j+i*jmax] = Z(numUdp,(i     )+imax*((j     )+jmax*(0     )));
            }
        }
        CREATE_FACE(4,  4,  1,  6,  0, imax, jmax);

        for (j = 0; j < jmax; j++) {
            for (i = 0; i < imax; i++) {
                x2d[i+j*imax] = X(numUdp,(i     )+imax*((j     )+jmax*(kmax-1)));
                y2d[i+j*imax] = Y(numUdp,(i     )+imax*((j     )+jmax*(kmax-1)));
                z2d[i+j*imax] = Z(numUdp,(i     )+imax*((j     )+jmax*(kmax-1)));
            }
        }
        CREATE_FACE(5,  2,  7,  3,  5, jmax, imax);

        /* make the Shell and SolidBody */
        if (outLevel >= 1) printf("        creating Shell\n");
        status = EG_makeTopology(context, NULL, SHELL, CLOSED,
                                 NULL, 6, efaces, NULL, &eshell);
        CHECK_STATUS(EG_makeTopology);

        if (outLevel >= 1) printf("        creating SolidBody\n");
        status = EG_makeTopology(context, NULL, BODY, SOLIDBODY,
                                 NULL, 1, &eshell, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* inform the user that there are 6 surface meshes */
        *nMesh = 6;

        /* set the output value(s) */

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
    }

cleanup:
    if (z2d != NULL) EG_free(z2d);
    if (y2d != NULL) EG_free(y2d);
    if (x2d != NULL) EG_free(x2d);

    if (strlen(message) > 0) {
        *string = message;
        printf("%s\n", message);
    } else if (status != EGADS_SUCCESS) {
        FREE(message);
        *string = udpErrorStr(status);
    } else {
        FREE(message);
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
    int iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

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


/*
 ************************************************************************
 *                                                                      *
 *   spline1d - create 1d cubic spline (uniform spacing, fixed ends)    *
 *                                                                      *
 ************************************************************************
 */

static int
spline1d(ego    context,
         int    imax,
         double *x,
         double *y,
         double *z,
         ego    *ecurv)
{
    int    status = EGADS_SUCCESS;

    int    i, kk, iknot, icp, iter, niter, header[4];
    double *cp=NULL, du, dx, dy, dz, data[18], dxyzmax;
    double *knotu=NULL, *CP=NULL;
    double dxyztol = 1.0e-7, relax = 0.10;

    ROUTINE(spline1d);

    icp   = imax + 2;
    iknot = imax + 6;

    /* indices associated with various arrays:
             x,y,z      knot       cp
               0         0*         0
                         1*
                                    1       (set by bc)
                         2*
                         3*
               1         4          2
               2         5          3

                        ...

             imax-1    imax       imax-2
             imax-2    imax+1     imax-1
                       imax+2*
                       imax+3*
                                  imax      (set by bc)
                       imax+4*
             imax-1    imax+5*    imax+1

        note: there are 4 repeateed knots at beginning and
              4 repeated knots at end */

    MALLOC(cp, double, (iknot+3*icp));

    knotu = &(cp[0    ]);
    CP    = &(cp[iknot]);

    /* create spline curve */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;

    kk = 0;

    /* knots (equally spaced) */
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;

    for (i = 1; i < imax; i++) {
        cp[kk++] = i;
    }

    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;

    /* initial control point */
    cp[kk++] = x[0];
    cp[kk++] = y[0];
    cp[kk++] = z[0];

    /* initial interior control point (for slope) */
    cp[kk++] = (3 * x[0] + x[1]) / 4;
    cp[kk++] = (3 * y[0] + y[1]) / 4;
    cp[kk++] = (3 * z[0] + z[1]) / 4;

    /* interior control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = x[i];
        cp[kk++] = y[i];
        cp[kk++] = z[i];
    }

    /* penultimate interior control point (for slope) */
    cp[kk++] = (3 * x[imax-1] + x[imax-2]) / 4;
    cp[kk++] = (3 * y[imax-1] + y[imax-2]) / 4;
    cp[kk++] = (3 * z[imax-1] + z[imax-2]) / 4;

    /* final control point */
    cp[kk++] = x[imax-1];
    cp[kk++] = y[imax-1];
    cp[kk++] = z[imax-1];

    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                             header, cp, ecurv);
    CHECK_STATUS(EG_makeGeometry);

    /* iterate to have knot evaluations match data points */
    niter = 1000;
    for (iter = 0; iter < niter; iter++) {
        dxyzmax = 0;

        /* beginning point is fixed */

        /* match finite-differenced slope d/du at beginning */
        i = 1;
        status = EG_evaluate(*ecurv, &(knotu[3]), data);
        CHECK_STATUS(EG_evaluate);

        du = knotu[4] - knotu[3];
        dx = x[1] - x[0] - du * data[3];
        dy = y[1] - y[0] - du * data[4];
        dz = z[1] - z[0] - du * data[5];

        /* natural end condition */
//$$$        dx = data[6] * du * du;
//$$$        dy = data[7] * du * du;
//$$$        dz = data[8] * du * du;

        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

        CP[3*(i)  ] += relax * dx;
        CP[3*(i)+1] += relax * dy;
        CP[3*(i)+2] += relax * dz;

        /* match interior points */
        for (i = 2; i < imax; i++) {
            status = EG_evaluate(*ecurv, &(knotu[i+2]), data);
            CHECK_STATUS(EG_evaluate);

            dx = x[i-1] - data[0];
            dy = y[i-1] - data[1];
            dz = z[i-1] - data[2];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*(i)  ] += dx;
            CP[3*(i)+1] += dy;
            CP[3*(i)+2] += dz;
        }

        /* match finite-differenced slope d/du at end */
        i = imax;
        status = EG_evaluate(*ecurv, &(knotu[imax+2]), data);
        CHECK_STATUS(EG_evaluate);

        du = knotu[imax+2] - knotu[imax+1];
        dx = x[imax-1] - x[imax-2] - du * data[3];
        dy = y[imax-1] - y[imax-2] - du * data[4];
        dz = z[imax-1] - z[imax-2] - du * data[5];

        /* natural end condition */
//$$$        dx = data[6] * du * du;
//$$$        dy = data[7] * du * du;
//$$$        dz = data[8] * du * du;

        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

        CP[3*(i)  ] -= relax * dx;
        CP[3*(i)+1] -= relax * dy;
        CP[3*(i)+2] -= relax * dz;

        /* convergence check */
        if (dxyzmax < dxyztol) break;

#ifdef GRAFIC
        {
            int    itype  = 0;
            int    io_kbd = 5;
            int    io_scr = 6;
            int    indgr  = 1 + 4 + 8 + 16 + 64;
            int    izero  = 0;
            int    nline, ilin[2], isym[2], nper[2];
            float  tplot[2001], xplot[2001], yplot[2001], zplot[2001];
            float  xmin, xmax, ymin, ymax;
            double t;
            char   pltitl[80];

            nline   =  1;
            ilin[0] = +1;
            isym[0] = -1;
            nper[0] = 1001;

            for (i = 0; i < 1001; i++) {
                t        = (imax-1) * (double)(i) / (double)(1000);
                status = EG_evaluate(*ecurv, &t, data);
                tplot[i] = (float)(t);
                xplot[i] = (float)(data[0]);
                yplot[i] = (float)(data[1]);
                zplot[i] = (float)(data[2]);
            }

            nline   =  2;
            ilin[1] = -2;
            isym[1] = +2;
            nper[1] = imax;

            for (i = 0; i < imax; i++) {
                tplot[1001+i] = (float)(  i );
                xplot[1001+i] = (float)(x[i]);
                yplot[1001+i] = (float)(y[i]);
                zplot[1001+i] = (float)(z[i]);
            }

            if (iter == 0) {
                grinit_(&io_kbd, &io_scr, "udpFreeform", strlen("udpFreeform"));
                grinpi_("0 for x, 1 for y, 2 for z", &itype, strlen("0 for x, 1 for y, 2 for z"));
            } else {
                indgr = 4 + 8 + 16 + 64;
                grsset_(&xmin, &xmax, &ymin, &ymax);
            }

            if        (itype == 0) {
                sprintf(pltitl, "~t~x~ iter=%d, dxyzmax=%12.3e", iter, dxyzmax);
                grline_(ilin, isym, &nline, pltitl, &indgr, tplot, xplot, nper, strlen(pltitl));
            } else if (itype == 1) {
                sprintf(pltitl, "~t~y~ iter=%d, dxyzmax=%12.3e", iter, dxyzmax);
                grline_(ilin, isym, &nline, pltitl, &indgr, tplot, yplot, nper, strlen(pltitl));
            } else {
                sprintf(pltitl, "~t~z~ iter=%d, dxyzmax=%12.3e", iter, dxyzmax);
                grline_(ilin, isym, &nline, pltitl, &indgr, tplot, zplot, nper, strlen(pltitl));
            }

            if (iter == 0) {
                grvalu_("XMINGR", &izero, &xmin, " ", strlen("XMINGR"), strlen(" "));
                grvalu_("XMAXGR", &izero, &xmax, " ", strlen("XMAXGR"), strlen(" "));
                grvalu_("YMINGR", &izero, &ymin, " ", strlen("YMINGR"), strlen(" "));
                grvalu_("YMAXGR", &izero, &ymax, " ", strlen("YMAXGR"), strlen(" "));
            }
        }
#endif

        /* make the new curve (after deleting old one) */
        status = EG_deleteObject(*ecurv);
        if (status != EGADS_SUCCESS) goto cleanup;

        status = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                                 header, cp, ecurv);
        if (status != EGADS_SUCCESS) goto cleanup;
    }

cleanup:
    if (cp != NULL) EG_free(cp);

    if (status != EGADS_SUCCESS) {
        if (ecurv != NULL) {
            EG_free(ecurv);
            *ecurv = NULL;
        }
    }

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   spline2d - create 2d cubic spline (uniform spacing, fixed ends)    *
 *                                                                      *
 ************************************************************************
 */

static int
spline2d(ego    context,
         int    imax,
         int    jmax,
         double *x,
         double *y,
         double *z,
         ego    *esurf)
{
    int    status = EGADS_SUCCESS;

    int    i, j, kk, iknot, jknot, icp, jcp, iter, niter, header[7];
    double *cp=NULL, du, dv, dx, dy, dz, data[18], dxyzmax, parms[2];
    double *knotu=NULL, *knotv=NULL, *CP=NULL;
    double dxyztol = 1.0e-7, relax = 0.10;

    ROUTINE(spline2d);

    icp   = imax + 2;
    iknot = imax + 6;
    jcp   = jmax + 2;
    jknot = jmax + 6;

    MALLOC(cp, double, (iknot+jknot+3*icp*jcp));

    knotu = &(cp[0          ]);
    knotv = &(cp[iknot      ]);
    CP    = &(cp[iknot+jknot]);

    /* create spline surface */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;
    header[4] = 3;
    header[5] = jcp;
    header[6] = jknot;

    kk = 0;

    /* knots in i-direction (equally spaced) */
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;

    for (i = 1; i < imax; i++) {
        cp[kk++] = i;
    }

    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;

    /* knots in j-direction (equally spaced) */
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;
    cp[kk++] = 0;

    for (j = 1; j < jmax; j++) {
        cp[kk++] = j;
    }

    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;
    cp[kk] = cp[kk-1]; kk++;

    /* map of point ID for imax=9 and jmax=5 (used in comments below)

             4   nw O  n  n  n  n  n  n  n  P ne
                 J  K  L  L  L  L  L  L  L  M  N
             3   w  H  *  *  *  *  *  *  *  I  e
             2   w  H  *  *  *  *  *  *  *  I  e
             1   w  H  *  *  *  *  *  *  *  I  e
                 C  D  E  E  E  E  E  E  E  F  G
             0   sw A  s  s  s  s  s  s  s  B se

                 0     1  2  3  4  5  6  7     8
    */

    /* southwest control point */
    cp[kk++] = x[(0)+(0)*imax];
    cp[kk++] = y[(0)+(0)*imax];
    cp[kk++] = z[(0)+(0)*imax];

    /* point A */
    cp[kk++] = (3 * x[(0)+(0)*imax] + x[(1)+(0)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(0)*imax] + y[(1)+(0)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(0)*imax] + z[(1)+(0)*imax]) / 4;

    /* south control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = x[(i)+(0)*imax];
        cp[kk++] = y[(i)+(0)*imax];
        cp[kk++] = z[(i)+(0)*imax];
    }

    /* point B */
    cp[kk++] = (3 * x[(imax-1)+(0)*imax] + x[(imax-2)+(0)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(0)*imax] + y[(imax-2)+(0)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(0)*imax] + z[(imax-2)+(0)*imax]) / 4;

    /* southeast control point */
    cp[kk++] = x[(imax-1)+(0)*imax];
    cp[kk++] = y[(imax-1)+(0)*imax];
    cp[kk++] = z[(imax-1)+(0)*imax];

    /* point C */
    cp[kk++] = (3 * x[(0)+(0)*imax] + x[(0)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(0)*imax] + y[(0)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(0)*imax] + z[(0)+(1)*imax]) / 4;

    /* point D */
    cp[kk++] = (3 * x[(0)+(0)*imax] + x[(1)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(0)*imax] + y[(1)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(0)*imax] + z[(1)+(1)*imax]) / 4;

    /* points E */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = (3 * x[(i)+(0)*imax] + x[(i)+(1)*imax]) / 4;
        cp[kk++] = (3 * y[(i)+(0)*imax] + y[(i)+(1)*imax]) / 4;
        cp[kk++] = (3 * z[(i)+(0)*imax] + z[(i)+(1)*imax]) / 4;
    }

    /* point F */
    cp[kk++] = (3 * x[(imax-1)+(0)*imax] + x[(imax-2)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(0)*imax] + y[(imax-2)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(0)*imax] + z[(imax-2)+(1)*imax]) / 4;

    /* point G */
    cp[kk++] = (3 * x[(imax-1)+(0)*imax] + x[(imax-1)+(1)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(0)*imax] + y[(imax-1)+(1)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(0)*imax] + z[(imax-1)+(1)*imax]) / 4;

    /* loop through interior j lines */
    for (j = 1; j < jmax-1; j++) {

        /* west control point */
        cp[kk++] = x[(0)+(j)*imax];
        cp[kk++] = y[(0)+(j)*imax];
        cp[kk++] = z[(0)+(j)*imax];

        /* point H */
        cp[kk++] = (3 * x[(0)+(j)*imax] + x[(1)+(j)*imax]) / 4;
        cp[kk++] = (3 * y[(0)+(j)*imax] + y[(1)+(j)*imax]) / 4;
        cp[kk++] = (3 * z[(0)+(j)*imax] + z[(1)+(j)*imax]) / 4;

        /* interior points */
        for (i = 1; i < imax-1; i++) {
            cp[kk++] = x[(i)+(j)*imax];
            cp[kk++] = y[(i)+(j)*imax];
            cp[kk++] = z[(i)+(j)*imax];
        }

        /* point I */
        cp[kk++] = (3 * x[(imax-1)+(j)*imax] + x[(imax-2)+(j)*imax]) / 4;
        cp[kk++] = (3 * y[(imax-1)+(j)*imax] + y[(imax-2)+(j)*imax]) / 4;
        cp[kk++] = (3 * z[(imax-1)+(j)*imax] + z[(imax-2)+(j)*imax]) / 4;

        /* east control point */
        cp[kk++] = x[(imax-1)+(j)*imax];
        cp[kk++] = y[(imax-1)+(j)*imax];
        cp[kk++] = z[(imax-1)+(j)*imax];
    }

    /* point J */
    cp[kk++] = (3 * x[(0)+(jmax-1)*imax] + x[(0)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(jmax-1)*imax] + y[(0)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(jmax-1)*imax] + z[(0)+(jmax-2)*imax]) / 4;

    /* point K */
    cp[kk++] = (3 * x[(0)+(jmax-1)*imax] + x[(1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(jmax-1)*imax] + y[(1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(jmax-1)*imax] + z[(1)+(jmax-2)*imax]) / 4;

    /* points L */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = (3 * x[(i)+(jmax-1)*imax] + x[(i)+(jmax-2)*imax]) / 4;
        cp[kk++] = (3 * y[(i)+(jmax-1)*imax] + y[(i)+(jmax-2)*imax]) / 4;
        cp[kk++] = (3 * z[(i)+(jmax-1)*imax] + z[(i)+(jmax-2)*imax]) / 4;
    }

    /* point M */
    cp[kk++] = (3 * x[(imax-1)+(jmax-1)*imax] + x[(imax-2)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(jmax-1)*imax] + y[(imax-2)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(jmax-1)*imax] + z[(imax-2)+(jmax-2)*imax]) / 4;

    /* point N */
    cp[kk++] = (3 * x[(imax-1)+(jmax-1)*imax] + x[(imax-1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(jmax-1)*imax] + y[(imax-1)+(jmax-2)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(jmax-1)*imax] + z[(imax-1)+(jmax-2)*imax]) / 4;

    /* northwest control point */
    cp[kk++] = x[(0)+(jmax-1)*imax];
    cp[kk++] = y[(0)+(jmax-1)*imax];
    cp[kk++] = z[(0)+(jmax-1)*imax];

    /* point O */
    cp[kk++] = (3 * x[(0)+(jmax-1)*imax] + x[(1)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * y[(0)+(jmax-1)*imax] + y[(1)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * z[(0)+(jmax-1)*imax] + z[(1)+(jmax-1)*imax]) / 4;

    /* north control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = x[(i)+(jmax-1)*imax];
        cp[kk++] = y[(i)+(jmax-1)*imax];
        cp[kk++] = z[(i)+(jmax-1)*imax];
    }

    /* point P */
    cp[kk++] = (3 * x[(imax-1)+(jmax-1)*imax] + x[(imax-2)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * y[(imax-1)+(jmax-1)*imax] + y[(imax-2)+(jmax-1)*imax]) / 4;
    cp[kk++] = (3 * z[(imax-1)+(jmax-1)*imax] + z[(imax-2)+(jmax-1)*imax]) / 4;

    /* northeast control point */
    cp[kk++] = x[(imax-1)+(jmax-1)*imax];
    cp[kk++] = y[(imax-1)+(jmax-1)*imax];
    cp[kk++] = z[(imax-1)+(jmax-1)*imax];

    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                             header, cp, esurf);
    CHECK_STATUS(EG_makeGeometry);

    /* iterate to have knot evaluations match data points */
    niter = 1000;
    for (iter = 0; iter < niter; iter++) {
        dxyzmax = 0;

        /* south boundary */
        {
            j = 0;

            /* sw corner is fixed */

            /* point A to match FD d/du */
            i = 1;
            parms[0] = knotu[3];
            parms[1] = knotv[3];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[4] - knotu[3];
            dx = x[(1)+(0)*imax] - x[(0)+(0)*imax] - du * data[3];
            dy = y[(1)+(0)*imax] - y[(0)+(0)*imax] - du * data[4];
            dz = z[(1)+(0)*imax] - z[(0)+(0)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* interior points along south boundary */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i+2];
                parms[1] = knotv[  3];
                status = EG_evaluate(*esurf, parms, data);
                CHECK_STATUS(EG_evaluate);

                dx = x[(i-1)+(0)*imax] - data[0];
                dy = y[(i-1)+(0)*imax] - data[1];
                dz = z[(i-1)+(0)*imax] - data[2];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += dx;
                CP[3*((i)+(j)*(imax+2))+1] += dy;
                CP[3*((i)+(j)*(imax+2))+2] += dz;
            }

            /* point B to match FD d/du */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[     3];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[imax+2] - knotu[imax+1];
            dx = x[(imax-1)+(0)*imax] - x[(imax-2)+(0)*imax] - du * data[3];
            dy = y[(imax-1)+(0)*imax] - y[(imax-2)+(0)*imax] - du * data[4];
            dz = z[(imax-1)+(0)*imax] - z[(imax-2)+(0)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* se corner is fixed */

        }

        /* line just above south boundary */
        {
            j = 1;

            /* point C to match FD d/dv */
            i = 0;
            parms[0] = knotu[3];
            parms[1] = knotv[3];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            dv = knotv[4] - knotv[3];
            dx = x[(0)+(1)*imax] - x[(0)+(0)*imax] - dv * data[6];
            dy = y[(0)+(1)*imax] - y[(0)+(0)*imax] - dv * data[7];
            dz = z[(0)+(1)*imax] - z[(0)+(0)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* point D to anhihilate d2/du/dv */
            i = 1;
            parms[0] = knotu[3];
            parms[1] = knotv[3];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[4] - knotu[3];
            dv = knotv[4] - knotv[3];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* points E to match FD d/dv */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i+2];
                parms[1] = knotv[  3];
                status = EG_evaluate(*esurf, parms, data);
                CHECK_STATUS(EG_evaluate);

                dv = knotv[4] - knotv[3];
                dx = x[(i-1)+(1)*imax] - x[(i-1)+(0)*imax] - dv * data[6];
                dy = y[(i-1)+(1)*imax] - y[(i-1)+(0)*imax] - dv * data[7];
                dz = z[(i-1)+(1)*imax] - z[(i-1)+(0)*imax] - dv * data[8];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
                CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
                CP[3*((i)+(j)*(imax+2))+2] += relax * dz;
            }

            /* point F to annihilate d2/du/dv */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[     3];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[imax+2] - knotu[imax+1];
            dv = knotv[     4] - knotv[     3];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* point G to match FD d/dv */
            i = imax + 1;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[     3];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            dv = knotv[4] - knotv[3];
            dx = x[(imax-1)+(1)*imax] - x[(imax-1)+(0)*imax] - dv * data[6];
            dy = y[(imax-1)+(1)*imax] - y[(imax-1)+(0)*imax] - dv * data[7];
            dz = z[(imax-1)+(1)*imax] - z[(imax-1)+(0)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;
        }

        /* interior j lines */
        for (j = 2; j < jmax; j++) {

            /* interior point along west boundary */
            i = 0;
            parms[0] = knotu[  3];
            parms[1] = knotv[j+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            dx = x[(0)+(j-1)*imax] - data[0];
            dy = y[(0)+(j-1)*imax] - data[1];
            dz = z[(0)+(j-1)*imax] - data[2];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += dx;
            CP[3*((i)+(j)*(imax+2))+1] += dy;
            CP[3*((i)+(j)*(imax+2))+2] += dz;

            /* point H to match FD d/du */
            i = 1;
            parms[0] = knotu[  3];
            parms[1] = knotv[j+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[4] - knotu[3];
            dx = x[(1)+(j-1)*imax] - x[(0)+(j-1)*imax] - du * data[3];
            dy = y[(1)+(j-1)*imax] - y[(0)+(j-1)*imax] - du * data[4];
            dz = z[(1)+(j-1)*imax] - z[(0)+(j-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* interior points */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i+2];
                parms[1] = knotv[j+2];
                status = EG_evaluate(*esurf, parms, data);
                CHECK_STATUS(EG_evaluate);

                dx = x[(i-1)+(j-1)*imax] - data[0];
                dy = y[(i-1)+(j-1)*imax] - data[1];
                dz = z[(i-1)+(j-1)*imax] - data[2];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += dx;
                CP[3*((i)+(j)*(imax+2))+1] += dy;
                CP[3*((i)+(j)*(imax+2))+2] += dz;
            }

            /* point I to match FD d/du */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[j   +2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[imax+2] - knotu[imax+1];
            dx = x[(imax-1)+(j-1)*imax] - x[(imax-2)+(j-1)*imax] - du * data[3];
            dy = y[(imax-1)+(j-1)*imax] - y[(imax-2)+(j-1)*imax] - du * data[4];
            dz = z[(imax-1)+(j-1)*imax] - z[(imax-2)+(j-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* interior point along east boundary */
            i = imax + 1;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[j   +2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            dx = x[(imax-1)+(j-1)*imax] - data[0];
            dy = y[(imax-1)+(j-1)*imax] - data[1];
            dz = z[(imax-1)+(j-1)*imax] - data[2];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += dx;
            CP[3*((i)+(j)*(imax+2))+1] += dy;
            CP[3*((i)+(j)*(imax+2))+2] += dz;
        }

        /* line just below north boundary */
        {
            j = jmax;

            /* point J to match FD d/dv */
            i = 0;
            parms[0] = knotu[     3];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = x[(0)+(jmax-1)*imax] - x[(0)+(jmax-2)*imax] - dv * data[6];
            dy = y[(0)+(jmax-1)*imax] - y[(0)+(jmax-2)*imax] - dv * data[7];
            dz = z[(0)+(jmax-1)*imax] - z[(0)+(jmax-2)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* point K to annihilate d2/du/dv */
            i = 1;
            parms[0] = knotu[     3];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[     4] - knotu[     3];
            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* points L to match FD d/dv */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i   +2];
                parms[1] = knotv[jmax+2];
                status = EG_evaluate(*esurf, parms, data);
                CHECK_STATUS(EG_evaluate);

                dv = knotv[jmax+2] - knotv[jmax+1];
                dx = x[(i-1)+(jmax-1)*imax] - x[(i-1)+(jmax-2)*imax] - dv * data[6];
                dy = y[(i-1)+(jmax-1)*imax] - y[(i-1)+(jmax-2)*imax] - dv * data[7];
                dz = z[(i-1)+(jmax-1)*imax] - z[(i-1)+(jmax-2)*imax] - dv * data[8];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
                CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
                CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;
            }

            /* point M to annihilate d2/du/dv */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[imax+2] - knotu[imax+1];
            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = du * dv * data[12];
            dy = du * dv * data[13];
            dz = du * dv * data[14];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* point N to match FD d/dv */
            i = imax + 1;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            dv = knotv[jmax+2] - knotv[jmax+1];
            dx = x[(imax-1)+(jmax-1)*imax] - x[(imax-1)+(jmax-2)*imax] - dv * data[6];
            dy = y[(imax-1)+(jmax-1)*imax] - y[(imax-1)+(jmax-2)*imax] - dv * data[7];
            dz = z[(imax-1)+(jmax-1)*imax] - z[(imax-1)+(jmax-2)*imax] - dv * data[8];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;
        }

        /* north boundary */
        {
            j = jmax + 1;

            /* nw corner is fixed */

            /* point O to match FD d/du */
            i = 1;
            parms[0] = knotu[     3];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[4] - knotu[3];
            dx = x[(1)+(jmax-1)*imax] - x[(0)+(jmax-1)*imax] - du * data[3];
            dy = y[(1)+(jmax-1)*imax] - y[(0)+(jmax-1)*imax] - du * data[4];
            dz = z[(1)+(jmax-1)*imax] - z[(0)+(jmax-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] += relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] += relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] += relax * dz;

            /* interior points along north boundary */
            for (i = 2; i < imax; i++) {
                parms[0] = knotu[i   +2];
                parms[1] = knotv[jmax+2];
                status = EG_evaluate(*esurf, parms, data);
                CHECK_STATUS(EG_evaluate);

                dx = x[(i-1)+(jmax-1)*imax] - data[0];
                dy = y[(i-1)+(jmax-1)*imax] - data[1];
                dz = z[(i-1)+(jmax-1)*imax] - data[2];

                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

                CP[3*((i)+(j)*(imax+2))  ] += dx;
                CP[3*((i)+(j)*(imax+2))+1] += dy;
                CP[3*((i)+(j)*(imax+2))+2] += dz;
            }

            /* point P to match FD d/du */
            i = imax;
            parms[0] = knotu[imax+2];
            parms[1] = knotv[jmax+2];
            status = EG_evaluate(*esurf, parms, data);
            CHECK_STATUS(EG_evaluate);

            du = knotu[imax+2] - knotu[imax+1];
            dx = x[(imax-1)+(jmax-1)*imax] - x[(imax-2)+(jmax-1)*imax] - du * data[3];
            dy = y[(imax-1)+(jmax-1)*imax] - y[(imax-2)+(jmax-1)*imax] - du * data[4];
            dz = z[(imax-1)+(jmax-1)*imax] - z[(imax-2)+(jmax-1)*imax] - du * data[5];

            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);

            CP[3*((i)+(j)*(imax+2))  ] -= relax * dx;
            CP[3*((i)+(j)*(imax+2))+1] -= relax * dy;
            CP[3*((i)+(j)*(imax+2))+2] -= relax * dz;

            /* se corner is fixed */

        }

        /* convergence check */
        if (dxyzmax < dxyztol) break;

        /* make the new curve (after deleting old one) */
        status = EG_deleteObject(*esurf);
        CHECK_STATUS(EG_deleteObject);

        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                                 header, cp, esurf);
        CHECK_STATUS(EG_makeGeometry);
    }

cleanup:
    if (cp != NULL) EG_free(cp);

    if (status != EGADS_SUCCESS) {
        if (esurf != NULL) {
            EG_free(esurf);
            *esurf = NULL;
        }
    }
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   plotData - level 3 GRAFIC routine for plotting raw data            *
 *                                                                      *
 ************************************************************************
 */

#ifdef GRAFIC
static void
plotData(int    *ifunct,                /* (in)  GRAFIC function indicator */
         void   *imaxP,                 /* (in)  points in i direction */
         void   *jmaxP,                 /* (in)  points in j direction */
         void   *kmaxP,                 /* (in)  points in k direction */
         void   *XP,                    /* (in)  x coordinates */
         void   *YP,                    /* (in)  y coordinates */
         void   *ZP,                    /* (in)  z coordinates */
         void   *a6,                    /* (in)  dummy GRAFC argument */
         void   *a7,                    /* (in)  dummy GRAFC argument */
         void   *a8,                    /* (in)  dummy GRAFC argument */
         void   *a9,                    /* (in)  dummy GRAFC argument */
         float  *scale,                 /* (out) array of scales */
         char   *text,                  /* (out) help text */
         int    textlen)                /* (in)  length of text */
{
    int     *imax   = (int    *) imaxP;
    int     *jmax   = (int    *) jmaxP;
    int     *kmax   = (int    *) kmaxP;
    double  *X      = (double *) XP;
    double  *Y      = (double *) YP;
    double  *Z      = (double *) ZP;

    int     i, j, k, ij, ijk;
    float   x4, y4, z4;
    double  xmin, xmax, ymin, ymax;

    int     iblack  = GR_BLACK;
    int     ired    = GR_RED;
    int     igreen  = GR_GREEN;

    ROUTINE(plotData);

    /* return scales */
    if (*ifunct == 0) {
        xmin = X[0];
        xmax = X[0];
        ymin = Y[0];
        ymax = Y[0];

        for (k = 0; k < *kmax; k++) {
            for (j = 0; j < *jmax; j++) {
                for (i = 0; i < *imax; i++) {
                    ijk = i + j * (*imax + k * (*jmax));
                    if (X[ijk] < xmin) xmin = X[ijk];
                    if (X[ijk] > xmax) xmax = X[ijk];
                    if (Y[ijk] < ymin) ymin = Y[ijk];
                    if (Y[ijk] > ymax) ymax = Y[ijk];
                }
            }
        }

        scale[0] = xmin;
        scale[1] = xmax;
        scale[2] = ymin;
        scale[3] = ymax;

        text[0] = '\0';

    /* plot image (only if *kmax <= 1) */
    } else if (*ifunct == 1 && *kmax <= 1) {

        /* i lines */
        grcolr_(&ired);
        for (j = 0; j < *jmax; j++) {
            i = 0;
            ij = i + j * (*imax);
            x4 = X[ij];
            y4 = Y[ij];
            z4 = Z[ij];
            grmov3_(&x4, &y4, &z4);

            for (i = 1; i < *imax; i++) {
                ij = i + j * (*imax);
                x4 = X[ij];
                y4 = Y[ij];
                z4 = Z[ij];
                grdrw3_(&x4, &y4, &z4);
            }
        }

        /* j lines */
        if (*jmax > 1) {
            grcolr_(&igreen);
            for (i = 0; i < *imax; i++) {
                j = 0;
                ij = i + j * (*imax);
                x4 = X[ij];
                y4 = Y[ij];
                z4 = Z[ij];
                grmov3_(&x4, &y4, &z4);

                for (j = 1; j < *jmax; j++) {
                    ij = i + j * (*imax);
                    x4 = X[ij];
                    y4 = Y[ij];
                    z4 = Z[ij];
                    grdrw3_(&x4, &y4, &z4);
                }
            }
        }

        grcolr_(&iblack);
    } else {
        printf("Illegal option selected\n");
    }

//cleanup:
    return;
}
#endif /* GRAFIC */
