/*
 ************************************************************************
 *                                                                      *
 * udpBezier -- udp file to generate a bezier wire, sheet, or solid body*
 *                                                                      *
 *            Written by Bridget Dixon @ Syracuse University            *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                      John Dannenhoffer @ Syracuse University         *
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

#define NUMUDPARGS 5
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)    ((char   *) (udps[IUDP].arg[0].val))
#define DEBUGIT( IUDP)    ((int    *) (udps[IUDP].arg[1].val))[0]
#define IMAX(    IUDP)    ((int    *) (udps[IUDP].arg[2].val))[0]
#define JMAX(    IUDP)    ((int    *) (udps[IUDP].arg[3].val))[0]
#define CP(      IUDP,I)  ((double *) (udps[IUDP].arg[4].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename", "debug", "imax",   "jmax",   "cp", };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,   ATTRINT, -ATTRINT, -ATTRINT, 0,    };
static int    argIdefs[NUMUDPARGS] = {0,          0,       0,        0,        0,    };
static double argDdefs[NUMUDPARGS] = {0.,         0.,      0.,       0.,       0.,   };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#define CHECK_STATUS2(X)                                           \
    if (status < EGADS_SUCCESS) {                                 \
        printf( "ERROR:: BAD STATUS = %d from %s\n", status, #X); \
        goto cleanup;                                             \
    } else if (DEBUGIT(0) != 0) {                                 \
        printf("%s -> status=%d\n", #X, status);                  \
    }

static void *realloc_temp=NULL;              /* used by RALLOC macro */


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

    int     i, j, imax, jmax, ipat, nipat, jpat, njpat, header[5], periodic;
    int     i0, i1, j0, j1, idum;
    int     inode, iedge, jedge, iface, icp, icp0, icp1, *senses=NULL, *esense=NULL;
    int     oclass, mtype, nchild, *senses2, count, newnode, newedge;
    double  xtemp, ytemp, ztemp, cp[48], trange[4], data[18], tol=1.0e-7;
    char    *message=NULL, *filename, *token;
    FILE    *fp=NULL;
    ego     ecurve, esurf, eloop, eshell, etemp[8], eref, *echild;
    ego     *enodes=NULL, *eedges=NULL, *efaces=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename(0) = %s\n", FILENAME(0));
    printf("debug(0)    = %d\n", DEBUGIT( 0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (STRLEN(FILENAME(0)) == 0) {
        snprintf(message, 100, "filename must not be null\n");
        status  = EGADS_NODATA;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
    printf("debug(%d)    = %d\n", numUdp, DEBUGIT( numUdp));
#endif

    /* open the file or stream */
    filename = FILENAME(numUdp);

    if (strncmp(filename, "<<\n", 3) == 0) {
        token = strtok(filename, " \t\n");
        if (token == NULL) {
            snprintf(message, 100, "premature enf-of-file found");
            status = EGADS_NODATA;
            goto cleanup;
        }
    } else {
        fp = fopen(filename, "r");
        if (fp == NULL) {
            snprintf(message, 100, "could not open file \"%s\"", filename);
            status = EGADS_NOTFOUND;
            goto cleanup;
        }
    }

    /* read the size of the bezier */
    if (fp != NULL) {
        count = fscanf(fp, "%d %d", &imax, &jmax);
        if (count != 2) {
            snprintf(message, 100, "error reading header");
            status = EGADS_NODATA;
            goto cleanup;
        }
    } else {
        token = strtok(NULL, " \t\n");
        if (token == NULL) {
            snprintf(message, 100, "error reading imax");
            status  = EGADS_NODATA;
            goto cleanup;
        } else {
            imax = strtol(token, NULL, 10);
        }

        token = strtok(NULL, " \t\n");
        if (token == NULL) {
            snprintf(message, 100, "error reading jmax");
            status  = EGADS_NODATA;
            goto cleanup;
        } else {
            jmax = strtol(token, NULL, 10);
        }
    }

    if        (imax < 4 || (imax-1)%3 != 0) {
        snprintf(message, 100, "bad value for imax=%d (from file)\n", imax);
        status  = EGADS_NODATA;
        goto cleanup;

    } else if (jmax < 1 || (jmax-1)%3 != 0) {
        snprintf(message, 100, "bad value for jmax=%d (from file)\n", jmax);
        status  = EGADS_NODATA;
        goto cleanup;
    }

    /* save the array size from the file */
    IMAX(numUdp) = imax;
    JMAX(numUdp) = jmax;

    RALLOC(udps[numUdp].arg[4].val, double, 3*imax*jmax);

    /* read the data */
    for (j = 0; j < jmax; j++) {
        for (i = 0; i < imax; i++) {

            if (fp != NULL) {
                count = fscanf(fp, "%lf %lf %lf", &xtemp, &ytemp, &ztemp);

                if (count == 3) {
                    CP(numUdp,3*(i+j*imax)  ) = xtemp;
                    CP(numUdp,3*(i+j*imax)+1) = ytemp;
                    CP(numUdp,3*(i+j*imax)+2) = ztemp;
                } else if (count <= 0) {
                    i--;
                } else {
                    snprintf(message, 100, "error reading point[%d,%d]", i, j);
                    status  = EGADS_NODATA;
                    goto cleanup;
                }
            } else {
                token = strtok(NULL, " \t\n");
                if (token == NULL) {
                    snprintf(message, 100, "error reading x[%d,%d]", i, j);
                    status  = EGADS_NODATA;
                    goto cleanup;
                } else {
                    xtemp = strtod(token, NULL);
                }

                token = strtok(NULL, " \t\n");
                if (token == NULL) {
                    snprintf(message, 100, "error reading x[%d,%d]", i, j);
                    status  = EGADS_NODATA;
                    goto cleanup;
                } else {
                    ytemp = strtod(token, NULL);
                }

                token = strtok(NULL, " \t\n");
                if (token == NULL) {
                    snprintf(message, 100, "error reading x[%d,%d]", i, j);
                    status  = EGADS_NODATA;
                    goto cleanup;
                } else {
                    ztemp = strtod(token, NULL);
                }
            }

            CP(numUdp,3*(i+j*imax)  ) = xtemp;
            CP(numUdp,3*(i+j*imax)+1) = ytemp;
            CP(numUdp,3*(i+j*imax)+2) = ztemp;
        }
    }

    /* close the file */
    if (fp != NULL) {
        fclose(fp);
    }

    /* create WireBody */
    if (JMAX(numUdp) == 1) {
        nipat = (imax - 1) / 3;

        if (DEBUGIT(0)) {
            printf("nipat=%d\n", nipat);
        }

        periodic = 1;
        i0 = 0;
        i1 = imax - 1;
        if (fabs(CP(numUdp,3*i0  ) - CP(numUdp,3*i1  )) > tol ||
            fabs(CP(numUdp,3*i0+1) - CP(numUdp,3*i1+1)) > tol ||
            fabs(CP(numUdp,3*i0+2) - CP(numUdp,3*i1+2)) > tol   ) {
            periodic = 0;
        }

        /* temporary storage used during construction */
        MALLOC(enodes, ego, (nipat+1));
        MALLOC(eedges, ego, (nipat  ));
        MALLOC(senses, int, (nipat  ));

        /* make the Nodes */
        for (ipat = 0; ipat <= nipat; ipat++) {
            if (ipat == nipat && periodic == 1) {
                enodes[ipat] = enodes[0];
                CHECK_STATUS2(EG_copyTopology(NODE));
            } else {
                status = EG_makeTopology(context, NULL, NODE, 0,
                                         &(CP(numUdp,9*ipat)), 0, NULL, NULL, &(enodes[ipat]));
                CHECK_STATUS2(EG_makeTopoloogy(NODE));
            }
        }

        /* loop through the patches of the bezier splines */
        for (ipat = 0; ipat < nipat; ipat++) {

            /* make the (cubic) bezier curve */
            header[0] = 0;
            header[1] = 3;
            header[2] = 4;

            status = EG_makeGeometry(context, CURVE, BEZIER, NULL,
                                     header, &(CP(numUdp,9*ipat)), &ecurve);
            CHECK_STATUS2(EG_makeGeometry(CURVE));

            status = EG_getRange(ecurve, trange, &idum);
            CHECK_STATUS2(EG_getRange);

            /* make the Edge */
            status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                     trange, 2, &(enodes[ipat]), NULL, &(eedges[ipat]));
            CHECK_STATUS2(EG_makeTopology(EDGE));

            senses[ipat] = SFORWARD;
        }

        /* make a Loop */
        if (periodic == 1) {
            status = EG_makeTopology(context, NULL, LOOP, CLOSED,
                                     NULL, nipat, eedges, senses, &eloop);
            CHECK_STATUS2(EG_makeTopology(LOOP-CLOSED));
        } else {
            status = EG_makeTopology(context, NULL, LOOP, OPEN,
                                     NULL, nipat, eedges, senses, &eloop);
            CHECK_STATUS2(EG_makeTopology(LOOP-OPEN));
        }

        /* make a WireBody */
        status = EG_makeTopology(context, NULL, BODY, WIREBODY,
                                 NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS2(EG_makeTopology(WIREBODY));

        /* set the output value(s) */
        IMAX(0) = imax;
        JMAX(0) = jmax;

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;

    /* create SolidBody (if closed) or SheetBody (if open) */
    } else {
        nipat = (imax - 1) / 3;
        njpat = (jmax - 1) / 3;

        if (DEBUGIT(0)) {
            printf("nipat=%d,  njpat=%d\n", nipat, njpat);
        }

        /* temporary storage used during construction */
        MALLOC(enodes, ego,   (nipat+1)*(njpat+1));
        MALLOC(eedges, ego, 2*(nipat+1)*(njpat+1));    /* actually bigger than needed */
        MALLOC(esense, int, 2*(nipat+1)*(njpat+1));    /* actually bigger than needed */
        MALLOC(efaces, ego,   (nipat  )*(njpat  ));
        MALLOC(senses, int,            4         );

        /* make (or copy) the Nodes */
        inode = 0;
        for (j0 = 0; j0 <= njpat; j0++) {
            for (i0 = 0; i0 <= nipat; i0++) {
                icp0    = 3 * (i0 + j0 * imax);
                newnode = 1;

                /* look to see if we have a Node at this location already */
                for (j1 = 0; j1 <= njpat; j1++) {
                    for (i1 = 0; i1 <= nipat; i1++) {
                        icp1 = 3 * (i1 + j1 * imax);

                        if (i1 == i0 && j1 == j0) {
                            break;
                        } else if (fabs(CP(numUdp,3*icp0  )-CP(numUdp,3*icp1  )) < tol &&
                                   fabs(CP(numUdp,3*icp0+1)-CP(numUdp,3*icp1+1)) < tol &&
                                   fabs(CP(numUdp,3*icp0+2)-CP(numUdp,3*icp1+2)) < tol   ) {
                            enodes[inode] = enodes[i1+j1*(nipat+1)];
                            CHECK_STATUS2(EG_copyTopology(NODE));
                            newnode = 0;
                            break;
                        }
                    }
                    if (newnode == 0 || (i1 == i0 && j1 == j0)) {
                        break;
                    }
                }

                if (newnode > 0) {
                    status = EG_makeTopology(context, NULL, NODE, 0,
                                             &(CP(numUdp,3*icp0)), 0, NULL, NULL, &(enodes[inode]));
                    CHECK_STATUS2(EG_makeTopology(NODE));
                }

                if (DEBUGIT(0)) {
                    printf("inode=%d, i0=%d, j0=%d, enodes=%llx\n",
                           inode, i0, j0, (long long)(enodes[inode]));
                }

                inode++;
            }
        }

        /* make (or copy) the j=constant Edges */
        if (DEBUGIT(0)) {
            printf("Making j=constant Edges...\n");
        }

        iedge = 0;
        for (jpat = 0; jpat <= njpat; jpat++) {
            for (ipat = 0; ipat < nipat; ipat++) {
                newedge  = 1;
                etemp[0] = enodes[ipat  +jpat*(nipat+1)];
                etemp[1] = enodes[ipat+1+jpat*(nipat+1)];

                jedge = 0;
                for (j1 = 0; j1 <= njpat; j1++) {
                    for (i1 = 0; i1 < nipat; i1++) {
                        etemp[2] = enodes[i1  +j1*(nipat+1)];
                        etemp[3] = enodes[i1+1+j1*(nipat+1)];

                        if (jedge >= iedge) {
                            break;

                        /* periodic Edge? */
                        } else if (etemp[0] != etemp[1] && etemp[0] == etemp[2] && etemp[1] == etemp[3]) {
                            eedges[iedge] =  eedges[jedge];
                            esense[iedge] = +esense[jedge];
                            CHECK_STATUS2(EG_copyGeometry(EDGE));
                            newedge = 0;
                            break;

                        /* anti-periodic edge? */
                        } else if (etemp[0] != etemp[1] && etemp[0] == etemp[3] && etemp[1] == etemp[2]) {
                            eedges[iedge] =  eedges[jedge];
                            esense[iedge] = -esense[jedge];
                            CHECK_STATUS2(EG_copyGeometry(EDGE));
                            newedge = 0;
                            break;
                        }
                        jedge++;
                    }
                    if (newedge == 0 || jedge >= iedge) break;
                }

                if (newedge == 0) {
                    /* copied above */

                } else if (etemp[0] == etemp[1]) {
                    trange[0] = 0;
                    trange[1] = 1;

                    status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                                             trange, 1, etemp, NULL, &(eedges[iedge]));
                    CHECK_STATUS2(EG_makeTopology(EDGE));
                    esense[iedge] = SFORWARD;
                } else {
                    /* make the (cubic) bezier curve */
                    header[0] = 0;
                    header[1] = 3;
                    header[2] = 4;

                    icp = 3 * (ipat + jpat * imax);
                    cp[ 0] = CP(numUdp,3*icp  );
                    cp[ 1] = CP(numUdp,3*icp+1);
                    cp[ 2] = CP(numUdp,3*icp+2); icp++;

                    cp[ 3] = CP(numUdp,3*icp  );
                    cp[ 4] = CP(numUdp,3*icp+1);
                    cp[ 5] = CP(numUdp,3*icp+2); icp++;

                    cp[ 6] = CP(numUdp,3*icp  );
                    cp[ 7] = CP(numUdp,3*icp+1);
                    cp[ 8] = CP(numUdp,3*icp+2); icp++;

                    cp[ 9] = CP(numUdp,3*icp  );
                    cp[10] = CP(numUdp,3*icp+1);
                    cp[11] = CP(numUdp,3*icp+2);

                    status = EG_makeGeometry(context, CURVE, BEZIER, NULL,
                                             header, cp, &ecurve);
                    CHECK_STATUS2(EG_makeGeometry(CURVE));

                    status = EG_getRange(ecurve, trange, &periodic);
                    CHECK_STATUS2(EG_getRange);


                    /* make the Edge */
                    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                             trange, 2, etemp, NULL, &(eedges[iedge]));
                    CHECK_STATUS2(EG_makeTopology(EDGE));
                    esense[iedge] = SFORWARD;
                }

                if (DEBUGIT(0)) {
                    printf("iedge=%d, ipat=%d, jpat=%d, eedges=%llx, esense=%d\n",
                           iedge, ipat, jpat, (long long)(eedges[iedge]), esense[iedge]);
                }

                iedge++;
            }
        }

        /* make (or copy) the i=constant Edges */
        if (DEBUGIT(0)) {
            printf("Making i=constant Edges...\n");
        }

        for (jpat = 0; jpat < njpat; jpat++) {
            for (ipat = 0; ipat <= nipat; ipat++) {
                newedge  = 1;
                etemp[0] = enodes[ipat+(jpat  )*(nipat+1)];
                etemp[1] = enodes[ipat+(jpat+1)*(nipat+1)];

                jedge = nipat * (njpat + 1);
                for (j1 = 0; j1 < njpat; j1++) {
                    for (i1 = 0; i1 <= nipat; i1++) {
                        etemp[2] = enodes[i1+(j1  )*(nipat+1)];
                        etemp[3] = enodes[i1+(j1+1)*(nipat+1)];

                        if (jedge >= iedge) {
                            break;

                        /* periodic Edge? */
                        } else if (etemp[0] != etemp[1] && etemp[0] == etemp[2] && etemp[1] == etemp[3]) {
                            eedges[iedge] =  eedges[jedge];
                            esense[iedge] = +esense[jedge];
                            CHECK_STATUS2(EG_copyGeometry(EDGE));
                            newedge = 0;
                            break;

                        /* anti-periodic edge? */
                        } else if (etemp[0] != etemp[1] && etemp[0] == etemp[3] && etemp[1] == etemp[2]) {
                            eedges[iedge] =  eedges[jedge];
                            esense[iedge] = -esense[jedge];
                            CHECK_STATUS2(EG_copyGeometry(EDGE));
                            newedge = 0;
                            break;
                        }
                        jedge++;
                    }
                    if (newedge == 0 || jedge >= iedge) break;
                }

                if (newedge == 0) {
                    /* copied above */

                } else if (etemp[0] == etemp[1]) {
                    trange[0] = 0;
                    trange[1] = 1;

                    status = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                                             trange, 1, etemp, NULL, &(eedges[iedge]));
                    CHECK_STATUS2(EG_makeTopology(EDGE));
                    esense[iedge] = SFORWARD;
                } else {
                    /* make the (cubic) bezier curve */
                    header[0] = 0;
                    header[1] = 3;
                    header[2] = 4;

                    icp = 3 * (ipat + jpat * imax);
                    cp[ 0] = CP(numUdp,3*icp  );
                    cp[ 1] = CP(numUdp,3*icp+1);
                    cp[ 2] = CP(numUdp,3*icp+2); icp+=imax;

                    cp[ 3] = CP(numUdp,3*icp  );
                    cp[ 4] = CP(numUdp,3*icp+1);
                    cp[ 5] = CP(numUdp,3*icp+2); icp+=imax;

                    cp[ 6] = CP(numUdp,3*icp  );
                    cp[ 7] = CP(numUdp,3*icp+1);
                    cp[ 8] = CP(numUdp,3*icp+2); icp+=imax;

                    cp[ 9] = CP(numUdp,3*icp  );
                    cp[10] = CP(numUdp,3*icp+1);
                    cp[11] = CP(numUdp,3*icp+2);

                    status = EG_makeGeometry(context, CURVE, BEZIER, NULL,
                                             header, cp, &ecurve);
                    CHECK_STATUS2(EG_makeGeometry(CURVE));

                    status = EG_getRange(ecurve, trange, &periodic);
                    CHECK_STATUS2(EG_getRange);

                    /* make the Edge */
                    status = EG_makeTopology(context, ecurve, EDGE, TWONODE,
                                             trange, 2, etemp, NULL, &(eedges[iedge]));
                    CHECK_STATUS2(EG_makeTopology(EDGE));
                    esense[iedge] = SFORWARD;
                }

                if (DEBUGIT(0)) {
                    printf("iedge=%d, ipat=%d, jpat=%d, eedges=%llx, esense=%d\n",
                           iedge, ipat, jpat, (long long)(eedges[iedge]), esense[iedge]);
                }

                iedge++;
            }
        }

        /* make each of the (bicubic) bezier surface patches */
        if (DEBUGIT(0)) {
            printf("Making Faces...\n");
        }

        iface = 0;
        for (jpat = 0; jpat < njpat; jpat++) {
            for (ipat = 0; ipat < nipat; ipat++) {

                /* make the (cubic) bezier curve */
                header[0] = 0;
                header[1] = 3;
                header[2] = 4;
                header[3] = 3;
                header[4] = 4;

                icp = 3 * (ipat + jpat * imax);
                cp[ 0] = CP(numUdp,3*icp  );
                cp[ 1] = CP(numUdp,3*icp+1);
                cp[ 2] = CP(numUdp,3*icp+2); icp++;

                cp[ 3] = CP(numUdp,3*icp  );
                cp[ 4] = CP(numUdp,3*icp+1);
                cp[ 5] = CP(numUdp,3*icp+2); icp++;

                cp[ 6] = CP(numUdp,3*icp  );
                cp[ 7] = CP(numUdp,3*icp+1);
                cp[ 8] = CP(numUdp,3*icp+2); icp++;

                cp[ 9] = CP(numUdp,3*icp  );
                cp[10] = CP(numUdp,3*icp+1);
                cp[11] = CP(numUdp,3*icp+2); icp+=imax-3;

                cp[12] = CP(numUdp,3*icp  );
                cp[13] = CP(numUdp,3*icp+1);
                cp[14] = CP(numUdp,3*icp+2); icp++;

                cp[15] = CP(numUdp,3*icp  );
                cp[16] = CP(numUdp,3*icp+1);
                cp[17] = CP(numUdp,3*icp+2); icp++;

                cp[18] = CP(numUdp,3*icp  );
                cp[19] = CP(numUdp,3*icp+1);
                cp[20] = CP(numUdp,3*icp+2); icp++;

                cp[21] = CP(numUdp,3*icp  );
                cp[22] = CP(numUdp,3*icp+1);
                cp[23] = CP(numUdp,3*icp+2); icp+=imax-3;

                cp[24] = CP(numUdp,3*icp  );
                cp[25] = CP(numUdp,3*icp+1);
                cp[26] = CP(numUdp,3*icp+2); icp++;

                cp[27] = CP(numUdp,3*icp  );
                cp[28] = CP(numUdp,3*icp+1);
                cp[29] = CP(numUdp,3*icp+2); icp++;

                cp[30] = CP(numUdp,3*icp  );
                cp[31] = CP(numUdp,3*icp+1);
                cp[32] = CP(numUdp,3*icp+2); icp++;

                cp[33] = CP(numUdp,3*icp  );
                cp[34] = CP(numUdp,3*icp+1);
                cp[35] = CP(numUdp,3*icp+2); icp+=imax-3;

                cp[36] = CP(numUdp,3*icp  );
                cp[37] = CP(numUdp,3*icp+1);
                cp[38] = CP(numUdp,3*icp+2); icp++;

                cp[39] = CP(numUdp,3*icp  );
                cp[40] = CP(numUdp,3*icp+1);
                cp[41] = CP(numUdp,3*icp+2); icp++;

                cp[42] = CP(numUdp,3*icp  );
                cp[43] = CP(numUdp,3*icp+1);
                cp[44] = CP(numUdp,3*icp+2); icp++;

                cp[45] = CP(numUdp,3*icp  );
                cp[46] = CP(numUdp,3*icp+1);
                cp[47] = CP(numUdp,3*icp+2);

                status = EG_makeGeometry(context, SURFACE, BEZIER, NULL,
                                         header, cp, &esurf);
                CHECK_STATUS2(EG_makeGeometry(SURFACE));

                /* make the Loop */
                etemp[ 0] =  eedges[ipat  + jpat * nipat];
                senses[0] =  esense[ipat  + jpat * nipat];

                etemp[ 1] =  eedges[ipat+1 + jpat * (nipat+1) + nipat * (njpat+1)];
                senses[1] =  esense[ipat+1 + jpat * (nipat+1) + nipat * (njpat+1)];

                etemp[ 2] =  eedges[ipat + (jpat+1) * nipat];
                senses[2] = -esense[ipat + (jpat+1) * nipat];

                etemp[ 3] =  eedges[ipat + jpat * (nipat+1) + nipat * (njpat+1)];
                senses[3] = -esense[ipat + jpat * (nipat+1) + nipat * (njpat+1)];

                if (senses[0] > 0) {
                    data[0] = 0;  data[1] = 0;  data[2] = +1;  data[3] = 0;
                } else {
                    data[0] = 1;  data[1] = 0;  data[2] = -1;  data[3] = 0;
                }
                status = EG_makeGeometry(context, PCURVE, LINE, esurf,
                                         NULL, data, &(etemp[4]));
                CHECK_STATUS2(EG_makeGeometry(PCURVE));

                if (senses[1] > 0) {
                    data[0] = 1;  data[1] = 0;  data[2] = 0;  data[3] = +1;
                } else {
                    data[0] = 1;  data[1] = 1;  data[2] = 0;  data[3] = -1;
                }
                status = EG_makeGeometry(context, PCURVE, LINE, esurf,
                                         NULL, data, &(etemp[5]));
                CHECK_STATUS2(EG_makeGeometry(PCURVE));

                if (senses[2] < 0) {
                    data[0] = 0;  data[1] = 1;  data[2] = +1;  data[3] = 0;
                } else {
                    data[0] = 1;  data[1] = 1;  data[2] = -1;  data[3] = 0;
                }
                status = EG_makeGeometry(context, PCURVE, LINE, esurf,
                                         NULL, data, &(etemp[6]));
                CHECK_STATUS2(EG_makeGeometry(PCURVE));

                if (senses[3] < 0) {
                    data[0] = 0;  data[1] = 0;  data[2] = 0;  data[3] = +1;
                } else {
                    data[0] = 0;  data[1] = 1;  data[2] = 0;  data[3] = -1;
                }
                status = EG_makeGeometry(context, PCURVE, LINE, esurf,
                                         NULL, data, &(etemp[7]));
                CHECK_STATUS2(EG_makeGeometry(PCURVE));

                status = EG_makeTopology(context, esurf, LOOP, CLOSED,
                                         NULL, 4, etemp, senses, &eloop);
                CHECK_STATUS2(EG_makeTopology(LOOP));

                /* make the Face */
                senses[0] = SFORWARD;
                status = EG_makeTopology(context, esurf, FACE, SFORWARD,
                                         NULL, 1, &eloop, senses, &(efaces[iface]));
                CHECK_STATUS2(EG_makeTopology);

                if (DEBUGIT(0)) {
                    printf("iface=%d, ipat=%d, jpat=%d, efaces=%llx\n", iface, ipat, jpat, (long long)(efaces[iface]));
                }

                iface++;
            }
        }

        /* make a Shell (assume CLOSED, but fixed below) */
        status = EG_makeTopology(context, NULL, SHELL, CLOSED,
                                 NULL, iface, efaces, NULL, &eshell);
        CHECK_STATUS2(EG_makeTopology(SHELL));

        /* determine if Shell is open or closed */
        status = EG_getTopology(eshell, &eref, &oclass, &mtype,
                                trange, &nchild, &echild, &senses2);
        CHECK_STATUS2(EG_getTopology);

        if (mtype == CLOSED) {
            /* make a SolidBody */
            status = EG_makeTopology(context, NULL, BODY, SOLIDBODY,
                                     NULL, 1, &eshell, NULL, ebody);
            CHECK_STATUS2(EG_makeTopology(SOLIDBODY));
        } else {
            /* make a SheetBody */
            status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                                     NULL, 1, &eshell, NULL, ebody);
            CHECK_STATUS2(EG_makeTopology(SHEETBODY));
        }

        /* set the output value(s) */
        IMAX(0) = imax;
        JMAX(0) = jmax;

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
    }


cleanup:
    if (enodes != NULL) EG_free(enodes);
    if (eedges != NULL) EG_free(eedges);
    if (esense != NULL) EG_free(esense);
    if (efaces != NULL) EG_free(efaces);
    if (senses != NULL) EG_free(senses);

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
