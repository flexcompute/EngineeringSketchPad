/*
 ************************************************************************
 *                                                                      *
 * udfDumpPmtrs -- udf file to dump all Parameters to a file as json dict *
 *                                                                      *
 *            Written by John Dannenhoffer                              *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2025  John F. Dannenhoffer, III (Syracuse University)
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

#define NUMUDPARGS        1
#define NUMUDPINPUTBODYS -1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(IUDP)  ((char   *) (udps[IUDP].arg[0].val))

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename", };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, };
static int    argIdefs[NUMUDPARGS] = {0,          };
static double argDdefs[NUMUDPARGS] = {1.,         };


/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototype for function defined below */


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  emodel,                 /* (in)  input Body */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, nchild, *senses, nrow, ncol, ipmtr, irow, icol, ival, icomma;
    double  data[4];
    char    *message=NULL;
    void    *modl;
    modl_T  *MODL;
    FILE    *fp;
    ego     context, eref, prev, next, *ebodys;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("FILENAME(0) = %s\n", FILENAME(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 1024);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size <= 1) {
        snprintf(message, 1024, "FILENAME must be given");
        status = EGADS_RANGERR;
        goto cleanup;
    }

    status = EG_getInfo(emodel, &oclass, &mtype, &eref, &prev, &next);
    CHECK_STATUS(EG_getInfo);

    if (oclass == CONTXT) {
        context = emodel;

    } else if (oclass == MODEL) {

        status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                                data, &nchild, &ebodys, &senses);
        CHECK_STATUS(EG_getTopology);

        if (nchild != 1) {
            snprintf(message, 1024, "expecting Model to contain one Body (not %d)\n", nchild);
            status = EGADS_NOTBODY;
            goto cleanup;
        }

        status = EG_getContext(emodel, &context);
        CHECK_STATUS(EG_getContext);

        /* make a copy of the Body (so that it does not get removed
           when OpenCSM deletes emodel) */
        status = EG_copyObject(ebodys[0], NULL, ebody);
        CHECK_STATUS(EG_copyObject);

    } else {
        snprintf(message, 1024, "emodel is neither a Contxt nor a Model");
        status = EGADS_NOTMODEL;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("FILENAME(%d) = %s\n", numUdp, FILENAME(numUdp));
#endif

    /* open up the file */
    fp = fopen(FILENAME(0), "w");
    if (fp == NULL) {
        snprintf(message, 1024, "File \"%s\" could not be opened for writing", FILENAME(0));
        status = EGADS_WRITERR;
        goto cleanup;
    }

    /* get a pointer to ocsm's MODL */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *) modl;

    /* write prolog */
    icomma = 0;
    fprintf(fp, "{");

    /* write all Parameters (of all types) */
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].name[0] == '@') continue;

        if (icomma > 0) {
            fprintf(fp, ", ");
        }

        fprintf(fp, "\'%s\' : ", MODL->pmtr[ipmtr].name);
        icomma = 1;

        nrow = MODL->pmtr[ipmtr].nrow;
        ncol = MODL->pmtr[ipmtr].ncol;
        ival = 0;

        /* string */
        if (nrow == 0) {
            fprintf(fp, "\'%s\'", MODL->pmtr[ipmtr].str);

        } else if (nrow == 1) {
            /* scalar */
            if (ncol == 1) {
                fprintf(fp, " %.9f", MODL->pmtr[ipmtr].value[ival++]);

            /* row vector */
            } else {
                fprintf(fp, "[");
                for (icol = 1; icol < ncol; icol++) {
                    fprintf(fp, " %.9f,", MODL->pmtr[ipmtr].value[ival++]);
                }
                fprintf(fp, " %.9f]", MODL->pmtr[ipmtr].value[ival++]);
            }

        /* column vector */
        } else if (MODL->pmtr[ipmtr].ncol == 1) {
            fprintf(fp, "[");
            for (irow = 1; irow < nrow; irow++) {
                fprintf(fp, " %.9f", MODL->pmtr[ipmtr].value[ival++]);
            }
            fprintf(fp, " %.9f]", MODL->pmtr[ipmtr].value[ival++]);

        /* matrix */
        } else {
            fprintf(fp, "[");
            for (irow = 1; irow <= nrow; irow++) {
                fprintf(fp, "[");
                for (icol = 1; icol < ncol; icol++) {
                    fprintf(fp, " %.9f,", MODL->pmtr[ipmtr].value[ival++]);
                }
                fprintf(fp, " %.9f]", MODL->pmtr[ipmtr].value[ival++]);
                if (irow < nrow) {
                    fprintf(fp, ", ");
                } else {
                    fprintf(fp, "]");
                }
            }
        }
    }

    fprintf(fp, "}");
    fclose(fp);

    /* set the output value(s) */

    /* remember this model (Body) */
    if (oclass == MODEL) {
        udps[numUdp].ebody = *ebody;
    } else {
        udps[numUdp].ebody = NULL;
    }

#ifdef DEBUG
    printf("udpExecute -> *ebody=%llx\n", (long long)(*ebody));
#endif

cleanup:
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
