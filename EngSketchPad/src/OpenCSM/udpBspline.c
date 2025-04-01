/*
 ************************************************************************
 *                                                                      *
 * udpBspline -- udp file to generate a bsplibne WireBody or SheetBody  *
 *                                                                      *
 *            Written by Bridget Dixon @ Syracuse University            *
 *            Patterned after code written by Bob Haimes  @ MIT         *
 *                      John Dannenhoffer @ Syracuse University         *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2013/2025  John F. Dannenhoffer, III (Syracuse University)
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

/* shorthands for accessing argument values and velocities */
#define BITFLAG( IUDP  )  ((int    *) (udps[IUDP].arg[0].val))[0]
#define UKNOTS(  IUDP,I)  ((double *) (udps[IUDP].arg[1].val))[I]
#define VKNOTS(  IUDP,I)  ((double *) (udps[IUDP].arg[2].val))[I]
#define CPS(     IUDP,I)  ((double *) (udps[IUDP].arg[3].val))[I]
#define WEIGHTS( IUDP,I)  ((double *) (udps[IUDP].arg[4].val))[I]
#define UDEGREE( IUDP  )  ((int    *) (udps[IUDP].arg[5].val))[0]
#define VDEGREE( IUDP  )  ((int    *) (udps[IUDP].arg[6].val))[0]
#define FILENAME(IUDP  )  ((char   *) (udps[IUDP].arg[7].val))

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"bitflag", "uknots", "vknots", "cps",    "weights", "udegree", "vdegree", "filename", };
static int    argTypes[NUMUDPARGS] = {ATTRINT,   ATTRREAL, ATTRREAL, ATTRREAL, ATTRREAL,  ATTRINT,   ATTRINT,   ATTRFILE,   };
static int    argIdefs[NUMUDPARGS] = {0,         0,        0,        0,        0,         3,         3,         0,          };
static double argDdefs[NUMUDPARGS] = {0.,        0.,       0.,       0.,       0.,        0.,        0.,        0.,         };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"


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

    int     i, header[7], ndata, narg, ieof, icpbeg, icpend;
    double  xyz[18], range[4], *rdata=NULL;
    char    *message=NULL, *filename, *token;
    FILE    *fp=NULL;
    ego     ecurve, enodes[2], eedge, eloop, esurface, eface, eshell;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("bitflag( 0) = %d\n", BITFLAG(0));
    printf("uknots(  0) =");
    for (i = 0; i < udps[0].arg[1].size; i++) {
        printf(" %f\n", UKNOTS(0,i));
    }
    printf("\n");
    printf("vknots(  0) =");
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf(" %f\n", VKNOTS(0,i));
    }
    printf("\n");
    printf("cps(     0) =");
    for (i = 0; i < udps[0].arg[3].size; i++) {
        printf(" %f\n", CPS(0,i));
    }
    printf("\n");
    printf("weights( 0) =");
    for (i = 0; i < udps[0].arg[4].size; i++) {
        printf(" %f\n", WEIGHTS(0,i));
    }
    printf("\n");
    printf("udegree( 0) = %d\n", UDEGREE(0));
    printf("vdegree( 0) = %d\n", VDEGREE(0));
    printf("filename(0) = %s\n", FILENAME(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (strlen(FILENAME(0)) == 0) {
        if (udps[0].arg[0].size > 1) {
            snprintf(message, 100, "\"bitflag\" should be a scalar");
            status  = EGADS_RANGERR;
            goto cleanup;
        }

        for (i = 1; i < udps[0].arg[1].size; i++) {
            if (UKNOTS(0,i) < UKNOTS(0,i-1)) {
                snprintf(message, 100, "\"uknots[%d]\" < \"uknots[%d]\"", i, i-1);
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }

        for (i = 1; i < udps[0].arg[2].size; i++) {
            if (VKNOTS(0,i) < VKNOTS(0,i-1)) {
                snprintf(message, 100, "\"vknots[%d]\" < \"vknots[%d]\"", i, i-1);
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("bitflag( %d) = %d\n", numUdp, BITFLAG(numUdp));
    printf("uknots(  %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[1].size; i++) {
        printf(" %f\n", UKNOTS(numUdp,i));
    }
    printf("\n");
    printf("vknots(  %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[2].size; i++) {
        printf(" %f\n", VKNOTS(numUdp,i));
    }
    printf("\n");
    printf("cps(     %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[3].size; i++) {
        printf(" %f\n", CPS(numUdp,i));
    }
    printf("\n");
    printf("weights( %d) =", numUdp);
    for (i = 0; i < udps[numUdp].arg[4].size; i++) {
        printf(" %f\n", WEIGHTS(numUdp,i));
    }
    printf("\n");
    printf("udegree( %d) = %d\n", numUdp, UDEGREE(numUdp));
    printf("vdegree( %d) = %d\n", numUdp, VDEGREE(numUdp));
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
#endif

    /* special handling if FILENAME is given */
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

        /* read header */
        if (fp != NULL) {
            for (i = 0; i < 7; i++) {
                narg = fscanf(fp, "%d", &header[i]);
                if (narg < 1) {
                    fclose(fp);

                    snprintf(message, 100, "trouble reading header[%d]", i);
                    status = EGADS_NODATA;
                    goto cleanup;
                }
            }
        } else {
            for (i = 0; i < 7; i++) {
                token = strtok(NULL, " \t\n");
                if (token == NULL) {
                    snprintf(message, 100, "error while reading imax");
                    status = EGADS_NODATA;
                    goto cleanup;
                }
                header[i] = strtol(token, NULL, 10);
            }
        }

        /* read rdata */
        ndata = header[3] + header[6] + 3 * header[2] * header[5];
        if ((header[0] & 0x2) != 0) {
            ndata += header[2] * header[5];
        }

        MALLOC(rdata, double, ndata);

        if (fp != NULL) {
            for (i = 0; i < ndata; i++) {
                narg = fscanf(fp, "%lf", &rdata[i]);
                if (narg < 1) {
                    fclose(fp);

                    snprintf(message, 100, "trouble reading rdata[%d]", i);
                    status = EGADS_NODATA;
                    goto cleanup;
                }
            }
        } else {
            for (i = 0; i < ndata; i++) {
                token = strtok(NULL, " \t\n");
                if (token == NULL) {
                    snprintf(message, 100, "error while reading rdata[%d]", i);
                    status = EGADS_NODATA;
                    goto cleanup;
                }
                rdata[i] = strtod(token, NULL);
            }
        }

        /* close the file */
        if (fp != NULL) {
            fclose(fp);
        }

        /* WireBody */
        if (header[5] == 1 && header[6] == 0) {

            status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rdata, &ecurve);
            CHECK_STATUS(EG_makeGeometry);

            range[0] = 0;
            range[1] = 1;

            status = EG_evaluate(ecurve, &range[0], xyz);
            CHECK_STATUS(EG_evaluate);

            status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, &enodes[0]);
            CHECK_STATUS(EG_makeTopology);

            status = EG_evaluate(ecurve, &range[1], xyz);
            CHECK_STATUS(EG_evaluate);

            icpbeg = header[3] + header[6];
            icpend = header[3] + header[6] + 3 * header[2] * header[5] - 3;

            if (header[0] & 0x4) {
                status = EG_makeTopology(context, ecurve, EDGE, ONENODE, range, 1, enodes, NULL, &eedge);
                CHECK_STATUS(EG_makeTopology);
            } else if (fabs(rdata[icpbeg  ]-rdata[icpend  ]) > EPS06 ||
                       fabs(rdata[icpbeg+1]-rdata[icpend+1]) > EPS06 ||
                       fabs(rdata[icpbeg+2]-rdata[icpend+2]) > EPS06   ) {
                status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, &enodes[1]);
                CHECK_STATUS(EG_makeTopology);

                status = EG_makeTopology(context, ecurve, EDGE, TWONODE, range, 2, enodes, NULL, &eedge);
                CHECK_STATUS(EG_makeTopology);
            } else {
                status = EG_makeTopology(context, ecurve, EDGE, ONENODE, range, 1, enodes, NULL, &eedge);
                CHECK_STATUS(EG_makeTopology);
            }

            status = EG_makeLoop(1, &eedge, NULL, 0, &eloop);
            CHECK_STATUS(EG_makeLoop);

            status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, ebody);
            CHECK_STATUS(EG_makeTopology);

            /* remember this model (body) */
            udps[numUdp].ebody = *ebody;
            goto cleanup;

        /* SheetBody */
        } else {
            status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rdata, &esurface);
            CHECK_STATUS(EG_makeGeometry);

            range[0] = 0;
            range[1] = 1;
            range[2] = 0;
            range[3] = 1;

            status = EG_makeFace(esurface, SFORWARD, range, &eface);
            CHECK_STATUS(EG_makeFace);

            status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL, 1, &eface, NULL, &eshell);
            CHECK_STATUS(EG_makeTopology);

            status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1, &eshell, NULL, ebody);
            CHECK_STATUS(EG_makeTopology);
        }

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;
    }

    /* make a WireBody */
    if (udps[0].arg[2].size == 1) {

        /* build header and real data */
        header[0] = BITFLAG(0);                             // bitflag

        header[1] = UDEGREE(0);                             // degree
        header[3] = udps[0].arg[1].size;                    // nknot
        header[2] = header[3] - header[1] - 1;              // ncp

        /* check that data sizes are consistent */
        if (header[2] != udps[0].arg[3].size/3) {
            snprintf(message, 100, "degree, knots, and cps do not agree");
            status = EGADS_RANGERR;
            goto cleanup;
        }

        if (header[0] & 0x2) {
            if (udps[0].arg[4].size != udps[0].arg[3].size/3) {
                snprintf(message, 100, "cps and knots do not agree");
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }

        MALLOC(rdata, double, udps[0].arg[1].size+udps[0].arg[3].size+udps[0].arg[4].size);

        ndata = 0;
        for (i = 0; i < udps[0].arg[1].size; i++) {         // knots
            rdata[ndata++] = UKNOTS(0,i);
        }
        for (i = 0; i < udps[0].arg[3].size; i++) {         // cps
            rdata[ndata++] = CPS(0,i);
        }
        if (header[0] & 0x2) {
            for (i = 0; i < udps[0].arg[4].size; i++) {     // weights
                rdata[ndata++] = WEIGHTS(0,i);
            }
        }

        /* build Curve, Nodes, Edge, Loop, and WireBody */
        status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rdata, &ecurve);
        CHECK_STATUS(EG_makeGeometry);

        range[0] = 0;
        range[1] = 1;

        status = EG_evaluate(ecurve, &range[0], xyz);
        CHECK_STATUS(EG_evaluate);

        status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, &enodes[0]);
        CHECK_STATUS(EG_makeTopology);

        status = EG_evaluate(ecurve, &range[1], xyz);
        CHECK_STATUS(EG_evaluate);

        if (header[0] & 0x4) {
            status = EG_makeTopology(context, ecurve, EDGE, ONENODE, range, 1, enodes, NULL, &eedge);
            CHECK_STATUS(EG_makeTopology);
        } else if (fabs(CPS(0,0)-CPS(0,3*header[2]-3)) > EPS06 ||
                   fabs(CPS(0,1)-CPS(0,3*header[2]-2)) > EPS06 ||
                   fabs(CPS(0,2)-CPS(0,3*header[2]-1)) > EPS06   ) {
            status = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL, &enodes[1]);
            CHECK_STATUS(EG_makeTopology);

            status = EG_makeTopology(context, ecurve, EDGE, TWONODE, range, 2, enodes, NULL, &eedge);
            CHECK_STATUS(EG_makeTopology);
        } else {
            status = EG_makeTopology(context, ecurve, EDGE, ONENODE, range, 1, enodes, NULL, &eedge);
            CHECK_STATUS(EG_makeTopology);
        }

        status = EG_makeLoop(1, &eedge, NULL, 0, &eloop);
        CHECK_STATUS(EG_makeLoop);

        status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;

    /* make a SheetBody */
    } else {

        /* build header and real data */
        header[0] = BITFLAG(0);                             // bitflag

        header[1] = UDEGREE(0);                             // udegree
        header[3] = udps[0].arg[1].size;                    // nuknot
        header[2] = header[3] - header[1] - 1;              // nucp

        header[4] = VDEGREE(0);                             // vdegree
        header[6] = udps[0].arg[2].size;                    // nvknot
        header[5] = header[6] - header[4] - 1;              // nvcp

        /* check that data sizes are consistent */
        if (header[2]*header[5] != udps[0].arg[3].size/3) {
            snprintf(message, 100, "degree, knots, and cps do not agree");
            status = EGADS_RANGERR;
            goto cleanup;
        }

        if (header[0] & 0x2) {
            if (udps[0].arg[4].size != udps[0].arg[3].size/3) {
                snprintf(message, 100, "cps and knots do not agree");
                status = EGADS_RANGERR;
                goto cleanup;
            }
        }

        MALLOC(rdata, double, udps[0].arg[1].size+udps[0].arg[2].size+udps[0].arg[3].size+udps[0].arg[4].size);

        ndata = 0;
        for (i = 0; i < udps[0].arg[1].size; i++) {         // uknots
            rdata[ndata++] = UKNOTS(0,i);
        }
        for (i = 0; i < udps[0].arg[2].size; i++) {         // vknots
            rdata[ndata++] = VKNOTS(0,i);
        }
        for (i = 0; i < udps[0].arg[3].size; i++) {         // cps
            rdata[ndata++] = CPS(0,i);
        }
        if (header[0] & 0x2) {
            for (i = 0; i < udps[0].arg[4].size; i++) {     // weight
                rdata[ndata++] = WEIGHTS(0,i);
            }
        }

        /* build Surface, Face, Shell, and SheetBody */
        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rdata, &esurface);
        CHECK_STATUS(EG_makeGeometry);

        range[0] = 0;
        range[1] = 1;
        range[2] = 0;
        range[3] = 1;

        status = EG_makeFace(esurface, SFORWARD, range, &eface);
        CHECK_STATUS(EG_makeFace);

        status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL, 1, &eface, NULL, &eshell);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1, &eshell, NULL, ebody);
        CHECK_STATUS(EG_makeTopology);

        /* remember this model (body) */
        udps[numUdp].ebody = *ebody;
        goto cleanup;

    }


#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
    FREE(rdata);

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
    int    status = EGADS_SUCCESS;

    int    iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    if (uvs != NULL) {
        printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f)\n",
               (long long)ebody, npnt, entType, entIndex, uvs[0]);
    } else {
        printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=NULL)\n",
               (long long)ebody, npnt, entType, entIndex);
    }
#endif

    /* check that ebody matches one of the ebodys */
    iudp = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        if (ebody == udps[judp].ebody) {
            iudp = judp;
            break;
        }
    }
    if (iudp <= 0) {
        status = EGADS_NOTMODEL;
        goto cleanup;
    }

    /* this routine is not written yet */
    status = EGADS_NOLOAD;

cleanup:
    return status;
}
