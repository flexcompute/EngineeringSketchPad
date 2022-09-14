/*
 ************************************************************************
 *                                                                      *
 * udpImport -- udp file to import from a STEP file                     *
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

#define NUMUDPARGS 3
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define FILENAME(  IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define BODYNUMBER(IUDP)  ((int    *) (udps[IUDP].arg[1].val))[0]
#define NUMBODIES( IUDP)  ((int    *) (udps[IUDP].arg[2].val))[0]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"filename",  "bodynumber", "numbodies", };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING,  ATTRINT,      -ATTRINT,    };
static int    argIdefs[NUMUDPARGS] = {0,           1,            0,           };
static double argDdefs[NUMUDPARGS] = {0.,          0.,           0.,          };

/* declarations needed to get a file's timestamp */
#include <time.h>
#include <sys/stat.h>

#ifdef WIN32
    #define TIMELONG unsigned long long
#else
    #define TIMELONG unsigned long
#endif

/* routine to remove private data */
#define FREEUDPDATA(A) freePrivateData(A)
static int freePrivateData(void *data);

static ego       emodel   = NULL;       /* both are cleaned up at the end */
static char     *filename = NULL;
static TIMELONG  datetime = 0;

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

    int      oclass, mtype, mtype2, *senses, tessState, npts;
    int      nbody, ibody, jbody, nface, iface, nedge, iedge, nnode, inode, attrtype, attrlen, nremove;
    CINT     *tempIlist;
    char     *message=NULL;
    ego      geom, *ebodys, *efaces, *eedges, *enodes, topRef, prev, next, myBody;
    TIMELONG dt;

#ifdef WIN32
    struct _stat64i32 buf;
#else
    struct  stat      buf;
#endif

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename(0)   = %s\n", FILENAME(  0));
    printf("bodynumber(0) = %d\n", BODYNUMBER(0));
#endif

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[1].size > 1) {
        snprintf(message, 100, "bodynumber should be a scalar");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (BODYNUMBER(0) == -1) {

    } else if (BODYNUMBER(0) <= 0) {
        snprintf(message, 100, "BodyNumber = %d < 0", BODYNUMBER(0));
        status  = EGADS_RANGERR;
        goto cleanup;
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("filename(  %d) = %s\n", numUdp, FILENAME(  numUdp));
    printf("bodynumber(%d) = %d\n", numUdp, BODYNUMBER(numUdp));
#endif


#ifdef WIN32
    status = _stat(FILENAME(0), &buf);
#else
    status =  stat(FILENAME(0), &buf);
#endif
    if (status < 0) {
        status = EGADS_NOTFOUND;
        goto cleanup;
    }
    dt = buf.st_mtime;

    /* if emodel already exists, but the old FILENAME and/or datetime does not
       agree with previous call, delete the old emodel and load it again (below) */
    if (emodel != NULL) {
        if ((strcmp(filename,FILENAME(0)) != 0) || (dt != datetime)) {
            EG_deleteObject(emodel);
            emodel = NULL;

            EG_free(filename);
            filename = NULL;

            datetime = 0;
        }
    }

    /* if we have not previously loaded the model (or just removed it above), load it now */
    if (emodel == NULL) {
        status = EG_loadModel(context, 0, FILENAME(0), &emodel);
        if ((status != EGADS_SUCCESS) || (emodel == NULL)) {
            goto cleanup;
        }

        /* remember the emodel and the file timestamp */
        udps[numUdp].data = emodel;
        filename = EG_strdup(FILENAME(0));
        if (filename == NULL) {
            EG_deleteObject(emodel);
            emodel = NULL;
            status = EGADS_MALLOC;
            goto cleanup;
        }
        datetime = dt;
    }

    status = EG_getTopology(emodel, &geom, &oclass, &mtype, NULL, &nbody,
                            &ebodys, &senses);
    CHECK_STATUS(EG_getTopology);

    nremove = 0;

    /* we are asking for entire Model */
    if (BODYNUMBER(0) == -1) {
        for (ibody = 0; ibody < nbody; ibody++) {
            /* remove the _hist and __trace__ attributes (if they exist) */
            status = EG_getBodyTopos(ebodys[ibody], NULL, FACE, &nface, &efaces);
            if (status == EGADS_SUCCESS) {
                for (iface = 0; iface < nface; iface++) {
                    status = EG_attributeRet(efaces[iface], "_hist",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing _hist attribute from Face %d\n", iface+1);
#endif
                        status = EG_attributeDel(efaces[iface], "_hist");
                        CHECK_STATUS(EG_attributeDel);
                    }

                    status = EG_attributeRet(efaces[iface], "__trace__",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing __trace__ attribute from Face %d\n", iface+1);
#endif
                        status = EG_attributeDel(efaces[iface], "__trace__");
                        CHECK_STATUS(EG_attributeDel);
                    }
                }

                EG_free(efaces);
            }

            status = EG_getBodyTopos(ebodys[ibody], NULL, EDGE, &nedge, &eedges);
            if (status == EGADS_SUCCESS) {
                for (iedge = 0; iedge < nedge; iedge++) {
                    status = EG_attributeRet(eedges[iedge], "_hist",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing _hist attribute from Edge %d\n", iedge+1);
#endif
                        status = EG_attributeDel(eedges[iedge], "_hist");
                        CHECK_STATUS(EG_attributeDel);
                    }

                    status = EG_attributeRet(eedges[iedge], "__trace__",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing __trace__ attribute from Edge %d\n", iedge+1);
#endif
                        status = EG_attributeDel(eedges[iedge], "__trace__");
                        CHECK_STATUS(EG_attributeDel);
                    }
                }

                EG_free(eedges);
            }

            status = EG_getBodyTopos(ebodys[ibody], NULL, NODE, &nnode, &enodes);
            if (status == EGADS_SUCCESS) {
                for (inode = 0; inode < nnode; inode++) {
                    status = EG_attributeRet(enodes[inode], "_hist",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing _hist attribute from Node %d\n", inode+1);
#endif
                        status = EG_attributeDel(enodes[inode], "_hist");
                        CHECK_STATUS(EG_attributeDel);
                    }

                    status = EG_attributeRet(enodes[inode], "__trace__",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing __trace__ attribute from Node %d\n", inode+1);
#endif
                        status = EG_attributeDel(enodes[inode], "__trace__");
                        CHECK_STATUS(EG_attributeDel);
                    }
                }

                EG_free(enodes);
            }
        }

        /* look at the entities in the Model that are not Bodys */
        for (ibody = nbody; ibody < mtype; ibody++) {

            /* process if a tessellation object */
            status = EG_getInfo(ebodys[ibody], &oclass, &mtype2, &topRef, &prev, &next);
            CHECK_STATUS(EG_getTopology);

            if (oclass != TESSELLATION) continue;

            /* find matching Body and put an attribute on it */
            status = EG_statusTessBody(ebodys[ibody], &myBody, &tessState, &npts);
            CHECK_STATUS(EG_statusTessBody);

            for (jbody = 0; jbody < nbody; jbody++) {
                if (myBody == ebodys[jbody]) {
                    status = EG_attributeAdd(ebodys[jbody], "__hasTess__", ATTRINT, 1,
                                             &ibody, NULL, NULL);
                    CHECK_STATUS(EG_attributeAdd);
                }
            }
        }

        status = EGADS_SUCCESS;

        /* return the whole Model */
        *ebody = emodel;

    /* extract the ebody */
    } else {
        if (BODYNUMBER(0) > nbody) {
            status = EGADS_RANGERR;
        } else {
            status = EG_copyObject(ebodys[BODYNUMBER(0)-1], NULL, ebody);
            CHECK_STATUS(EG_copyObject);
            if (*ebody == NULL) goto cleanup;     // needed for splint

            /* remove the _hist and __trace__ attributes (if they exist) */
            status = EG_getBodyTopos(*ebody, NULL, FACE, &nface, &efaces);
            if (status == EGADS_SUCCESS) {
                for (iface = 0; iface < nface; iface++) {
                    status = EG_attributeRet(efaces[iface], "_hist",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing _hist attribute from Face %d\n", iface+1);
#endif
                        status = EG_attributeDel(efaces[iface], "_hist");
                        CHECK_STATUS(EG_attributeDel);
                    }

                    status = EG_attributeRet(efaces[iface], "__trace__",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing __trace__ attribute from Face %d\n", iface+1);
#endif
                        status = EG_attributeDel(efaces[iface], "__trace__");
                        CHECK_STATUS(EG_attributeDel);
                    }
                }

                EG_free(efaces);
            }

            status = EG_getBodyTopos(*ebody, NULL, EDGE, &nedge, &eedges);
            if (status == EGADS_SUCCESS) {
                for (iedge = 0; iedge < nedge; iedge++) {
                    status = EG_attributeRet(eedges[iedge], "_hist",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing _hist attribute from Edge %d\n", iedge+1);
#endif
                        status = EG_attributeDel(eedges[iedge], "_hist");
                        CHECK_STATUS(EG_attributeDel);
                    }

                    status = EG_attributeRet(eedges[iedge], "__trace__",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing __trace__ attribute from Edge %d\n", iedge+1);
#endif
                        status = EG_attributeDel(eedges[iedge], "__trace__");
                        CHECK_STATUS(EG_attributeDel);
                    }
                }

                EG_free(eedges);
            }

            status = EG_getBodyTopos(*ebody, NULL, NODE, &nnode, &enodes);
            if (status == EGADS_SUCCESS) {
                for (inode = 0; inode < nnode; inode++) {
                    status = EG_attributeRet(enodes[inode], "_hist",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing _hist attribute from Node %d\n", inode+1);
#endif
                        status = EG_attributeDel(enodes[inode], "_hist");
                        CHECK_STATUS(EG_attributeDel);
                    }

                    status = EG_attributeRet(enodes[inode], "__trace__",
                                             &attrtype, &attrlen, &tempIlist, NULL, NULL);
                    if (status == EGADS_SUCCESS) {
                        nremove++;
#ifdef DEBUG
                        printf(" udpExecute: removing __trace__ attribute from Node %d\n", inode+1);
#endif
                        status = EG_attributeDel(enodes[inode], "__trace__");
                        CHECK_STATUS(EG_attributeDel);
                    }
                }

                EG_free(enodes);
            }

            status = EGADS_SUCCESS;
        }

        if (status != EGADS_SUCCESS) {
            ebody = NULL;
            goto cleanup;
        }
    }

    /* we should not remove the emodel, since it will be deleted in freePrivateData */
//$$$
//$$$    emodel = NULL;

    if (nremove > 0) {
        printf("WARNING:: %d _hist and/or __trace__ attributes removed\n", nremove);
    }

    /* set the output value(s) */
    NUMBODIES(0) = nbody;

    /* remember this model (body) */
    udps[numUdp].ebody = *ebody;

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


/*
 ************************************************************************
 *                                                                      *
 *   freePrivateData - free private data                                *
 *                                                                      *
 ************************************************************************
 */

static int
freePrivateData(
    /*@unused@*/void  *data)            /* (in)  pointer to private data */
{

    if (emodel != NULL) {
        EG_deleteObject(emodel);
        emodel = NULL;
    }

    if (filename != NULL) {
        EG_free(filename);
        filename = NULL;
    }

    datetime = 0;

    return EGADS_SUCCESS;
}
