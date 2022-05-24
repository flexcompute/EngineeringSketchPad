 /*
 ************************************************************************
 *                                                                      *
 * udpNurbbody -- udp file to generate a Body from NURBS                *
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

/* this can act as a UDF and write a nurbs.txt file from a Body */
#define NUMUDPINPUTBODYS 0
#define NUMUDPARGS 1
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(      IUDP)   ((char   *) (udps[IUDP].arg[ 0].val))

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"filename", };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,   };
static int    argIdefs[NUMUDPARGS] = {0,          };
static double argDdefs[NUMUDPARGS] = {0.,         };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"
#include <sys/stat.h>
#include <stdio.h>
#include <assert.h>

#ifdef WIN32
   #define SLASH '\\'
#else
   #include <unistd.h>
   #include <libgen.h>
   #define SLASH '/'
#endif

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

#define           HUGEQ           99999999.0
#define           PIo2            1.5707963267948965579989817
#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

#ifdef WIN32
// Returns filename portion of the given path
// Returns empty string if path is directory
const char *basename(const char *path)
{
  const char *filename = strrchr(path, '\\');
  if (filename == NULL)
    filename = path;
  else
    filename++;
  return filename;
}
#endif

static void *realloc_temp=NULL;              /* used by RALLOC macro */


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  context_in,             /* (in)  Egads context (or model) */
           ego  *ebody,                 /* (out) Body pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     oclass, mtype, nchild, *senses;
    int     nface, iface, i, n, nument, header[7], ieof;
    int     j, nuknots, nvknots, ncpu, ncpv, *tempI=NULL;
    double  data[18], *rdata=NULL, *tempR=NULL;
    char    *message=NULL, *filename, *token;
    FILE    *fp=NULL, *fp_output;
    ego     eref, esurf, *ebodys, *efaces, *echilds, context, emodel, emodel2, *efaces2=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context_in=%llx)\n", (long long)context_in);
    printf("filename(0) = %s\n", FILENAME(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(%d) = %s\n", numUdp, FILENAME(numUdp));
#endif

    /* special mode to write nurbs.txt if running as UDF */
    if (NUMUDPINPUTBODYS > 0) {

        /* store the context into emodel in case it is run as a UDF */
        emodel = context_in;

        /* check that Model was input that contains one Body */
        status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                                data, &nchild, &ebodys, &senses);
        CHECK_STATUS(EG_getTopology);

        if (oclass != MODEL) {
            snprintf(message, 100, "expecting a Model\n");
            status = EGADS_NOTMODEL;
            goto cleanup;
        } else if (nchild != 1) {
            snprintf(message, 100, "expecting Model to contain one Body (not %d)\n", nchild);
            status = EGADS_NOTBODY;
            goto cleanup;
        }

        status = EG_getContext(emodel, &context);
        CHECK_STATUS(EG_getContext);

        status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
        CHECK_STATUS(EG_getBodyTopos);

        fp_output = fopen("nurbs.txt", "w");
        if (fp_output == NULL) {
            snprintf(message, 100, "could not open file \"nurbs.txt\"\n");
            status = EGADS_NOTFOUND;
            goto cleanup;
        }

        for (iface = 0; iface < nface; iface++) {
            status = EG_getTopology(efaces[iface], &esurf, &oclass, &mtype,
                                    data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &tempI, &tempR);
            CHECK_STATUS(EG_getGeometry);
            if (tempI == NULL) goto cleanup;    // needed for splint
            if (tempR == NULL) goto cleanup;    // needed for splint

            if (oclass != SURFACE || mtype != BSPLINE) {
                snprintf(message, 100, "not a bspline surface\n");
                status = EGADS_NOTGEOM;
                goto cleanup;
            }

            fprintf(fp_output, "%5d %5d %5d %5d %5d %5d %5d\n",
                    tempI[0], tempI[1], tempI[2], tempI[3], tempI[4], tempI[5], tempI[6]);

            ncpu    = tempI[2];
            nuknots = tempI[3];
            ncpv    = tempI[5];
            nvknots = tempI[6];

            j = 0;
            for (i = 0; i < nuknots; i++) {
                fprintf(fp_output, "%20.13e\n", tempR[j++]);
            }
            for (i = 0; i < nvknots; i++) {
                fprintf(fp_output, "%20.13e\n", tempR[j++]);
            }
            for (i = 0; i < ncpu*ncpv; i++) {
                fprintf(fp_output, "%20.13e %20.13e %20.13e\n", tempR[j], tempR[j+1], tempR[j+2]);
                j += 3;
            }
            if (tempI[0]%2 == 1) {
                for (i = 0; i < ncpu*ncpv; i++) {
                    fprintf(fp_output, "%20.13e\n", tempR[j++]);
                }
            }

            EG_free(tempI);
            EG_free(tempR);
        }

        fclose(fp_output);

        EG_free(efaces);

        /* make a copy of the Body (so that it does not get removed
           when OpenCSM deletes emodel) */
        status = EG_copyObject(ebodys[0], NULL, ebody);
        CHECK_STATUS(EG_copyObject);

    /* normal processing */
    } else {

        context = context_in;

        /* open the file or stream */
        filename = FILENAME(numUdp);

        if (filename == NULL) {
            snprintf(message, 100, "NULL filename");
            status = EGADS_NOTFOUND;
            goto cleanup;

        } else if (strncmp(filename, "<<\n", 3) == 0) {
            token = strtok(filename, " \t\n");
            if (token == NULL) {
                ieof = 1;
            } else {
                ieof = 0;
            }
        } else {
            fp = fopen(filename, "rb");
            if (fp == NULL) {
                snprintf(message, 100, "could not open file \"%s\"", filename);
                status = EGADS_NOTFOUND;
                goto cleanup;
            }

            ieof = feof(fp);
        }

        nface = 0;
        MALLOC(efaces2, ego, 1);

        while (ieof == 0) {

            /* read the header */
            for (i = 0; i < 7; i++) {
                if (fp != NULL) {
                    nument = fscanf(fp, "%d", &header[i]);
                    if (nument != 1) {
                        ieof = 1;
                        break;
                    }
                } else {
                    token = strtok(NULL, " \t\n");
                    if (token == NULL) {
                        ieof = 1;
                        break;
                    }
                    header[i] = strtol(token, NULL, 10);
                }
            }
            if (ieof == 1) break;

            RALLOC(efaces2, ego, nface+1);

            /* read the (double) data */
            n = header[3] + header[6] + (3 + header[0]%2) * header[2] * header[5];
            MALLOC(rdata, double, n);

            for (i = 0; i < n; i++) {
                if (fp != NULL) {
                    nument = fscanf(fp, "%lf", &rdata[i]);
                    if (nument != 1) {
                        snprintf(message, 100, "error while reading rdata[%d]", i);
                        status = EGADS_NODATA;
                        goto cleanup;
                    }
                } else {
                    token = strtok(NULL, " \t\n");
                    if (token == NULL) {
                        snprintf(message, 100, "error while reading rdata[%d]", i);
                        status = EGADS_NODATA;
                        goto cleanup;
                    }
                    rdata[i] = strtod(token, NULL);
                }
            }

            /* find the minimum and maximum knot values */
            data[0] = rdata[                    0];
            data[1] = rdata[header[3]          -1];
            data[2] = rdata[header[3]            ];
            data[3] = rdata[header[3]+header[6]-1];

            /* make the BSPLINE Surface */
            status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                                     header, rdata, &esurf);
            CHECK_STATUS(EG_makeGeometry);

            FREE(rdata);

            /* make the Face */
            status = EG_makeFace(esurf, SFORWARD, data, &efaces2[nface]);
            CHECK_STATUS(EG_makeFace);

            nface++;
        }

        if (fp != NULL) {
            fclose(fp);
        }

        /* sew the Faces into a new Model */
        status = EG_sewFaces(nface, efaces2, 0, 0, &emodel2);
        CHECK_STATUS(EG_sewFaces);

        FREE(efaces2);

        /* get the Body out of the new Model */
        status = EG_getTopology(emodel2, &eref, &oclass, &mtype, data,
                                &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);

        /* make a copy of the Body so that it does not get destroyed
           when the Model is deleted */
        status = EG_copyObject(echilds[0], NULL, ebody);
        CHECK_STATUS(EG_copyObject);

        status = EG_deleteObject(emodel2);
        CHECK_STATUS(EG_deleteObject);
    }

    /* return the copied Body */
    udps[numUdp].ebody = *ebody;

cleanup:
    FREE(efaces2);
    FREE(rdata  );

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
    int    iudp, judp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d, uvs=%f %f)\n",
           (long long)ebody, npnt, entType, entIndex, uvs[0], uvs[1]);
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
        return EGADS_NOTMODEL;
    }

    /* this routine is not written yet */
    return EGADS_NOLOAD;
}
