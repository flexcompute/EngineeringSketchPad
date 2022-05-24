/*
 ************************************************************************
 *                                                                      *
 * udfSlices -- make slices of a Body                                   *
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

#define NUMUDPARGS       2
#define NUMUDPINPUTBODYS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define NSLICE(IUDP)  ((int    *) (udps[IUDP].arg[0].val))[0]
#define DIRN(  IUDP)  ((char   *) (udps[IUDP].arg[1].val))

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"nslice", "dirn",     };
static int    argTypes[NUMUDPARGS] = {ATTRINT,  ATTRSTRING, };
static int    argIdefs[NUMUDPARGS] = {0,        0,          };
static double argDdefs[NUMUDPARGS] = {0.,       0.,         };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototype for function defined below */

#include "OpenCSM.h"

static void *realloc_temp=NULL;              /* used by RALLOC macro */


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  emodel,                 /* (in)  Model containing Body */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, mtype2, nchild, nchild2, ichild, *senses;
    int     idir, nedge, nslice, islice, one=1, haveTparams, attrType, attrLen;
    CINT    *tempIlist;
    double  bbox[6], xslice=0, yslice=0, zslice=0, data[18];
    CDOUBLE *tempRlist;
    char    *message=NULL;
    CCHAR   *tempClist;
    ego     context, eref, *ebodys, *echilds, esurface, eface, eshell, esheet, emodel2, *eslices=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(emodel=%llx)\n", (long long)emodel);
    printf("nslice(0) = %d\n", NSLICE(0));
    printf("dirn(  0) = %s\n", DIRN(  0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check that Model was input that contains one Body */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass != MODEL) {
        snprintf(message, 100, "expecting a Model\n");
        status = EGADS_NOTMODEL;
        goto cleanup;
    } else if (nchild != 1) {
        snprintf(message, 10, "expecting Model to contain one Body (not %d)\n", nchild);
        status = EGADS_NOTBODY;
        goto cleanup;
    }

    status = EG_getContext(emodel, &context);
    CHECK_STATUS(EG_getContext);

    /* check arguments */
    if (udps[0].arg[0].size > 1) {
        snprintf(message, 100, "nslice should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (NSLICE(0) <= 0) {
        snprintf(message, 100, "nslice = %d <= 0", NSLICE(0));
        status  =  EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(emodel);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("nslice(%d) = %d\n", numUdp, NSLICE(numUdp));
    printf("dirn(  %d) = %s\n", numUdp, DIRN(  numUdp));
#endif

    if        (strncmp(DIRN(0), "x", 1) == 0 || strncmp(DIRN(0), "X", 1) == 0) {
        idir = 1;
    } else if (strncmp(DIRN(0), "y", 1) == 0 || strncmp(DIRN(0), "Y", 1) == 0) {
        idir = 2;
    } else if (strncmp(DIRN(0), "z", 1) == 0 || strncmp(DIRN(0), "Z", 1) == 0) {
        idir = 3;
    } else {
        idir = 1;
    }

    /* if there is a .tParams on ebodys[0], remember it and put it on all slices */
    status = EG_attributeRet(ebodys[0], ".tParams", &attrType, &attrLen,
                             &tempIlist, &tempRlist, &tempClist);
    if (status == SUCCESS) {
        haveTparams = 1;
    } else {
        haveTparams = 0;
    }

    /* get an array to store the slices */
    MALLOC(eslices, ego, 1);

    /* get the bounding box if the input Body */
    status = EG_getBoundingBox(ebodys[0], bbox);
    CHECK_STATUS(EG_getBoundingBox);

    /* get the type of the input Body */
    status = EG_getTopology(ebodys[0], &eref, &oclass, &mtype2, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    /* make the slices */
    nslice = 0 ;
    for (islice = 0; islice < NSLICE(0); islice++) {
        if (idir == 1) {
            xslice = bbox[0] + (bbox[3] - bbox[0]) * (double)(islice + 1) / (double)(NSLICE(0) + 1);

            data[0] = xslice;   data[1] = 0;   data[2] = 0;
            data[3] = 0;        data[4] = 1;   data[5] = 0;
            data[6] = 0;        data[7] = 0;   data[8] = 1;
            status = EG_makeGeometry(context, SURFACE, PLANE,
                                     NULL, NULL, data, &esurface);
            CHECK_STATUS(EG_makeGeometry);

            data[0] = bbox[1]-1;   data[1] = bbox[4]+1;
            data[2] = bbox[2]-1;   data[3] = bbox[5]+1;
            status = EG_makeFace(esurface, SFORWARD, data, &eface);
            CHECK_STATUS(EG_makeFace);
        } else if (idir == 2) {
            yslice = bbox[1] + (bbox[4] - bbox[1]) * (double)(islice + 1) / (double)(NSLICE(0) + 1);

            data[0] = 0;   data[1] = yslice;   data[2] = 0;
            data[3] = 0;   data[4] = 0;        data[5] = 1;
            data[6] = 1;   data[7] = 0;        data[8] = 0;
            status = EG_makeGeometry(context, SURFACE, PLANE,
                                     NULL, NULL, data, &esurface);
            CHECK_STATUS(EG_makeGeometry);

            data[0] = bbox[2]-1;   data[1] = bbox[5]+1;
            data[2] = bbox[0]-1;   data[3] = bbox[3]+1;
            status = EG_makeFace(esurface, SFORWARD, data, &eface);
            CHECK_STATUS(EG_makeFace);
        } else {
            zslice = bbox[2] + (bbox[5] - bbox[2]) * (double)(islice + 1) / (double)(NSLICE(0) + 1);

            data[0] = 0;   data[1] = 0;   data[2] = zslice;
            data[3] = 1;   data[4] = 0;   data[5] = 0;
            data[6] = 0;   data[7] = 1;   data[8] = 0;
            status = EG_makeGeometry(context, SURFACE, PLANE,
                                     NULL, NULL, data, &esurface);
            CHECK_STATUS(EG_makeGeometry);

            data[0] = bbox[0]-1;   data[1] = bbox[3]+1;
            data[2] = bbox[1]-1;   data[3] = bbox[4]+1;
            status = EG_makeFace(esurface, SFORWARD, data, &eface);
            CHECK_STATUS(EG_makeFace);
        }

        status = EG_makeTopology(context, NULL, SHELL, OPEN,
                                 NULL, 1, &eface, NULL, &eshell);
        CHECK_STATUS(EG_makeTopology);

        status = EG_makeTopology(context, NULL, BODY, SHEETBODY,
                                 NULL, 1, &eshell, NULL, &esheet);
        CHECK_STATUS(EG_makeTopology);

        if (mtype2 == SOLIDBODY) {
            status = EG_generalBoolean(ebodys[0], esheet, INTERSECTION, 0, &emodel2);
            CHECK_STATUS(EG_generalBoolean);
        } else if (mtype2 == SHEETBODY || mtype2 == FACEBODY) {
            status = EG_intersection(ebodys[0], esheet, &nedge, NULL, &emodel2);
            CHECK_STATUS(EG_intersection);
        } else {
            snprintf(message, 100, "input Body must be SolidBody or SheetBody");
            status = EGADS_CONSTERR;
            goto cleanup;
        }

        status = EG_deleteObject(esheet);
        CHECK_STATUS(EG_deleteObject);

        status = EG_getTopology(emodel2, &eref, &oclass, &mtype, data, &nchild2, &echilds, &senses);
        if (status == EGADS_SUCCESS) {
            RALLOC(eslices, ego, nslice+nchild2);

            for (ichild = 0; ichild < nchild2; ichild++) {
                status = EG_copyObject(echilds[ichild], NULL, &eslices[nslice]);
                CHECK_STATUS(EG_copyObject);

                if (haveTparams == 1) {
                    status = EG_attributeAdd(eslices[nslice], ".tParams", attrType, attrLen,
                                             tempIlist, tempRlist, tempClist);
                    CHECK_STATUS(EG_attributeAdd);
                }

                /* tell OpenCSM that the Faces do not have a _body attribute */
                status = EG_attributeAdd(eslices[nslice], "__markFaces__", ATTRINT, 1,
                                         &one, NULL, NULL);
                CHECK_STATUS(EG_attributeAdd);

                nslice++;
            }
        } else if (idir == 1) {
            printf(" udpExecute: slice at %10.5f was not generated\n", xslice);
        } else if (idir == 2) {
            printf(" udpExecute: slice at %10.5f was not generated\n", yslice);
        } else {
            printf(" udpExecute: slice at %10.5f was not generated\n", zslice);
        }

        status = EG_deleteObject(emodel2);
        CHECK_STATUS(EG_deleteObject);
    }

    /* make a Model of the slices */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL,
                             nslice, eslices, NULL, ebody);
    CHECK_STATUS(EG_makeTopology);

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(eslices);
    
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
