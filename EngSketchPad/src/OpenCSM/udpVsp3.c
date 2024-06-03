/*
 ************************************************************************
 *                                                                      *
 * udpVsp3 -- read a .vsp3 file (or a .stp file exported from vsp)      *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *                                                                      *
 ************************************************************************
 */

/*
 * Copyright (C) 2011/2024  John F. Dannenhoffer, III (Syracuse University)
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
#define FILENAME( IUDP)  ((char   *) (udps[IUDP].arg[0].val))
#define KEEPTEMPS(IUDP)  ((int    *) (udps[IUDP].arg[1].val))[0]

/* data about possible arguments */
static char*  argNames[NUMUDPARGS] = {"filename", "keeptemps", "recycle",   };
static int    argTypes[NUMUDPARGS] = {ATTRFILE,   ATTRINT,     ATTRRECYCLE, };
static int    argIdefs[NUMUDPARGS] = {0,          0,           0,           };
static double argDdefs[NUMUDPARGS] = {0.,         0.,          0.,          };

/* routines for handling private data */
#define FREEUDPDATA(A) freePrivateData(A)
static int freePrivateData(void *data);

#define COPYUDPDATA(SRC,TGT) copyPrivateData(SRC, TGT)
static int copyPrivateData(/*@null@*/void *src, void **tgt);

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

#include "egads_dot.h"

#ifdef GRAFIC
   #include "grafic.h"
#endif

#ifdef WIN32
    #define  SLASH       '\\'
    #define  SLEEP(msec)  Sleep(msec)
#else
    #include <unistd.h>
    #define  SLASH       '/'
    #define  SLEEP(msec)  usleep(1000*msec)
#endif

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

#define           EPS06           1.0e-06
#define           EPS12           1.0e-12
#define           EPS20           1.0e-20
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))

/* message (to be shared by all functions) */
static char message[1024];

static int runVsp3(modl_T *MODL, char filename[], ego  *emodel, int *NumUdp, udp_T *udps);
static int processStepFile(ego context, char filename[], ego *ebody);

typedef struct {
    int     ipmtr;            /* Parameter index (1:npmtr) or -1 for end */
    int     irow;             /* row       index (1:nrow) */
    int     icol;             /* column    index (1:ncol) */
    double  value;            /* value */
    double  dot;              /* velocity */
} despmtr_T;

typedef struct {
    int       magic;          /* magic number */
    int       ifirst;         /* first Body in MODL created by this call */
    int       ilast;          /* last  Body in MODL created by this call */
    int       npmtr;          /* number of DESPMTRs */
    despmtr_T *despmtrs;      /* array  of DESPMTRs */
    ego       emodel;         /* MODEL of SolidBodys created in this call */
} privdata_T;


/*
 ************************************************************************
 *                                                                      *
 *   udpExecute - execute the primitive                                 *
 *                                                                      *
 ************************************************************************
 */

int
udpExecute(ego  context,                /* (in)  EGADS context */
           ego  *emodel,                /* (out) MODEL pointer */
           int  *nMesh,                 /* (out) number of associated meshes */
           char *string[])              /* (out) error message */
{
    int     status = EGADS_SUCCESS;

    int     lenname, ipmtr, ptype, nrow, ncol, i, iudp, judp, nchange;
    int     oclass, mtype, nchild, *senses, attrType, attrLen;
    double  data[18], value, dot;
    char    filename[MAX_LINE_LEN+1], pname[MAX_NAME_LEN];
    CCHAR   *tempClist;
    ego     eref, *echilds;
    void    *modl;
    modl_T  *MODL;
    udp_T   *udps = *Udps;

    int        npmtr, irow, icol;
    privdata_T *privdata=NULL;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("filename( 0) = %s\n", FILENAME( 0));
    printf("keeptemps(0) = %d\n", KEEPTEMPS(0));
#endif

    /* default return values */
    *emodel = NULL;
    *nMesh  = 0;
    *string = NULL;

    message[0] = '\0';

    /* check arguments */
    if (udps[0].arg[0].size == 0) {
        snprintf(message, 1023, "\"filename\" must be given");
        status  = EGADS_RANGERR;
        goto cleanup;

    } else if (udps[0].arg[1].size > 1) {
        snprintf(message, 1023, "\"keeptemps\" must be a scalar");
        status = EGADS_RANGERR;
        goto cleanup;

    }

    /* get pointer to OpenCSM MODL */
    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)modl;
    if (MODL->perturb != NULL) {
        MODL = MODL->perturb;
    }

    /* if there is private data from a previous call, see if the current
       DESPMTRs match */
    iudp    = 0;
    nchange = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        privdata = (privdata_T *)(udps[judp].data);
        if (privdata        == NULL   ) continue;
        if (privdata->magic != 4433341) continue;

        status = EG_attributeRet(privdata->emodel, "vsp3Filename", &attrType, &attrLen,
                                 NULL, NULL, &tempClist);
        CHECK_STATUS(EG_attributeRet);

        if (strcmp(FILENAME(0), tempClist) != 0) continue;

        iudp    = judp;
        nchange = 0;
        for (ipmtr = 0; ipmtr < privdata->npmtr; ipmtr++) {
            status = ocsmGetValu(MODL, privdata->despmtrs[ipmtr].ipmtr,
                                       privdata->despmtrs[ipmtr].irow,
                                       privdata->despmtrs[ipmtr].icol, &value, &dot);
            CHECK_STATUS(ocsmGetValu);

            if (fabs(privdata->despmtrs[ipmtr].value-value) > EPS12) {
                nchange++;
                break;
            }
        }
        break;
    }

    if (iudp > 0 && nchange <=  0) {
#ifdef DEBUG
        printf("we can use previously generated MODEL\n");
#endif

        privdata = (privdata_T *)(udps[iudp].data);

        *emodel = privdata->emodel;
    } else {
#ifdef DEBUG
        printf("we need to generate a new MODEL (iudp=%d, nchange=%d)\n",
               iudp, nchange);
#endif
        privdata = NULL;

        /* cache copy of arguments for future use */
        status = cacheUdp(NULL);
        CHECK_STATUS(cacheUdp);

#ifdef DEBUG
        printf("udpExecute(context=%llx)\n", (long long)context);
        printf("filename( %d) = %s\n", numUdp, FILENAME( numUdp));
        printf("keeptemps(%d) = %d\n", numUdp, KEEPTEMPS(numUdp));
#endif

        /* process based upon type file filename given */
        strncpy(filename, FILENAME(numUdp), MAX_LINE_LEN);
        lenname = strlen(filename);

        /* since OpenCSM converts all forward slashes in filenames to backslashes
           on Windows (before we enter the UDP), we need to convert them back
           since anglescript (in vspscript) only allows forward slashes in filenames */
        for (i = 0; i < lenname; i++) {
            if (filename[i] == SLASH) filename[i] = '/';
        }

        /* filename is a .stp file */
        if        (lenname > 4 && strcmp(&filename[lenname-4], ".stp" ) == 0) {
            status = processStepFile(context, filename, emodel);
            CHECK_STATUS(processStepFile);

        /* filename is a .vsp3 file */
        } else if (lenname > 5 && strcmp(&filename[lenname-5], ".vsp3") == 0) {

            status = runVsp3(MODL, filename, emodel, NumUdp, udps);
            CHECK_STATUS(runVsp3);

        /* unknown filename type */
        } else {
            snprintf(message, 1023, "\"%s\" is not a .stp or .vsp3 file", filename);
            status = EGADS_RANGERR;
            goto cleanup;
        }

        SPLINT_CHECK_FOR_NULL(*emodel);

        /* put a filename attribute on the Model so that we can determine
           the correct Model to use during sensitivities */
        status = EG_attributeAdd(*emodel, "vsp3Filename", ATTRSTRING, strlen(filename),
                                 NULL, NULL, filename);
        CHECK_STATUS(EG_attributeAdd);

        /* store a table of the DESPMTRs in the private data */
        npmtr = 0;
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            status = ocsmGetPmtr(MODL, ipmtr, &ptype, &nrow, &ncol, pname);
            CHECK_STATUS(ocsmGetPmtr);

            if (ptype == OCSM_DESPMTR) {
                npmtr += nrow * ncol;
            }
        }

        MALLOC(privdata, privdata_T, 1);
        privdata->magic = 4433341;

        status = EG_getTopology(*emodel, &eref, &oclass, &mtype,
                                data, &nchild, &echilds, &senses);
        CHECK_STATUS(EG_getTopology);
        privdata->ifirst   = MODL->nbody + 1;
        privdata->ilast    = MODL->nbody + nchild;
        privdata->npmtr    = npmtr;
        privdata->despmtrs = NULL;
        privdata->emodel   = *emodel;

        MALLOC(privdata->despmtrs, despmtr_T, npmtr);

        npmtr = 0;
        for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
            status = ocsmGetPmtr(MODL, ipmtr, &ptype, &nrow, &ncol, pname);
            CHECK_STATUS(ocsmGetPmtr);

            if (ptype == OCSM_DESPMTR) {
                for (irow = 1; irow <= nrow; irow++) {
                    for (icol = 1; icol <= ncol; icol++) {
                        status = ocsmGetValu(MODL, ipmtr, irow, icol, &value, &dot);
                        CHECK_STATUS(ocsmGetValu);

                        privdata->despmtrs[npmtr].ipmtr  = ipmtr;
                        privdata->despmtrs[npmtr].irow   = irow;
                        privdata->despmtrs[npmtr].icol   = icol;
                        privdata->despmtrs[npmtr].value  = value;
                        privdata->despmtrs[npmtr].dot    = dot;

                        npmtr++;
                    }
                }
            }
        }

        udps[numUdp].data = (void *)(privdata);
    }

    /* set the output value(s) */

    /* remember this model (body) */
    /*@-kepttrans@*/
    udps[numUdp].ebody = *emodel;
    /*@+kepttrans@*/

cleanup:
    if (strlen(message) > 0) {
        MALLOC(*string, char, 1024);
        strncpy(*string, message, 1023);
    } else if (status != EGADS_SUCCESS) {
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
               int    npnt,             /* (in)  number of points */
               int    entType,          /* (in)  OCSM entity type */
               int    entIndex,         /* (in)  OCSM entity index (bias-1) */
               double uvs[],            /* (in)  parametric coordinates for evaluation */
               double vels[])           /* (out) velocities */
{
    int    status = EGADS_SUCCESS;

    int    iudp, judp, ipnt;

    int        npmtr, ipmtr;
    privdata_T *privdata0=NULL, *privdata1=NULL;

    int     ptype, nrow, irow, ncol, icol, builtTo, nbody;
    int     ibody, jbody, iprim, i, periodic, attrType, attrLen;
    int     oclass, oclass0, oclass1, mtype, mtype0, mtype1, nchild, nchild0, nchild1, nbody0, nbody1;
    int     nnode0, nnode1, inode, jnode, knode;
    int     nedge0, nedge1, iedge, jedge, kedge;
    int     nface0, nface1, iface, jface, kface;
    int     ichild, *senses, ntemp, iknot, icp;
    int     *header0, *header1;
    double  value, dot, data[18], uvs_dot[2], data_dot[18], *rdata0=NULL, *rdata1=NULL, *rdata_dot=NULL;
    double  trange[2], trange_dot[2];
    double  dtime = 0.00001;
    char    pname[MAX_NAME_LEN];
    CCHAR   *filename0, *filename1;
    ego     context, eref;
    ego     *echilds, *echilds0, *echilds1, *ebodys0, *ebodys1, *efaces0, *efaces1, *eedges0, *eedges1, *enodes0, *enodes1;
    ego     esurf, esurf0, esurf1, ecurve, ecurve0, ecurve1, *etemps;
    void    *modl, *ptrb;
    modl_T  *MODL=NULL, *PTRB=NULL;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpSensitivity(ebody=%llx, npnt=%d, entType=%d, entIndex=%d)\n", (long long)ebody, npnt, entType, entIndex);
#endif

    /* get the MODL and Body index */
    status = EG_getContext(ebody, &context);
    CHECK_STATUS(EG_getContext);

    status = EG_getUserPointer(context, (void**)(&(modl)));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)(modl);

    /* determine the iudp by comparing ebody with those
       between ifirst and ilast */
    iudp  = 0;
    ibody = 0;
    for (judp = 1; judp <= numUdp; judp++) {
        privdata0 = (privdata_T *)(udps[judp].data);

        for (jbody = privdata0->ifirst; jbody <= privdata0->ilast; jbody++) {
            if (ebody == MODL->body[jbody].ebody) {
                iudp  = judp;
                ibody = jbody;
                break;
            }
        }
        if (iudp > 0) break;
    }

    SPLINT_CHECK_FOR_NULL(privdata0);

    if (iudp == 0 || ibody == 0) {
        printf("ebody cannot be found\n");
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* if we already have _dots, we can skip directly to evaluation */
    status = EG_hasGeometry_dot(ebody);

    if (status == EGADS_SUCCESS && udps[iudp].ndotchg == 0) {
#ifdef DEBUG
        printf("we already have sensitivity info\n");
#endif
        goto evaluate;
    }

#ifdef DEBUG
    printf("we need to set up sensitivity (status=%d, ndotchg[%d]=%d)\n",
           status, iudp, udps[iudp].ndotchg);
#endif

    /* we are taking care of sensitivities, so reset the ndotchg flag */
    for (judp = 1; judp <= numUdp; judp++) {
        udps[judp].ndotchg = 0;
    }

    /* save the latest DESPMTRs */
    npmtr = 0;
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        status = ocsmGetPmtr(MODL, ipmtr, &ptype, &nrow, &ncol, pname);
        CHECK_STATUS(ocsmGetPmtr);

        if (ptype == OCSM_DESPMTR) {
            for (irow = 1; irow <= nrow; irow++) {
                for (icol = 1; icol <= ncol; icol++) {
                    status = ocsmGetValu(MODL, ipmtr, irow, icol, &value, &dot);
                    CHECK_STATUS(ocsmGetValu);

                    privdata0->despmtrs[npmtr].ipmtr = ipmtr;
                    privdata0->despmtrs[npmtr].irow  = irow;
                    privdata0->despmtrs[npmtr].icol  = icol;
                    privdata0->despmtrs[npmtr].value = value;
                    privdata0->despmtrs[npmtr].dot   = dot;

                    npmtr++;
                }
            }
        }
    }

    /* create a perturbation */
#ifdef DEBUG
    printf("creating perturbation\n");
#endif
    status = ocsmCopy(MODL, &ptrb);
    CHECK_STATUS(ocsmCopy);

    PTRB = (modl_T *)(ptrb);

    /* remove reference to udpVsp3 private data in the perturbation
       (since we do not want the emodel associated with MODL
       to get removed when we ocsmFree(PTRB) at the end) */
    for (iprim = 0; iprim < MAXPRIM; iprim++) {
        for (i = PTRB->NumUdp[iprim]; i > 0; i--) {
            privdata1 = (privdata_T *)(PTRB->Udps[iprim][i].data);
            if (privdata1        == NULL   ) break;
            if (privdata1->magic != 4433341) break;

            PTRB->Udps[iprim][i].data = NULL;
            (PTRB->NumUdp[iprim])--;
        }
    }

    /* make the perturbed DESPMTRs in PTRB */
    for (ipmtr = 0; ipmtr < npmtr; ipmtr++) {
        if (privdata0->despmtrs[ipmtr].dot != 0) {
            status = ocsmSetValuD(PTRB, privdata0->despmtrs[ipmtr].ipmtr,
                                  privdata0->despmtrs[ipmtr].irow,
                                  privdata0->despmtrs[ipmtr].icol,
                                  privdata0->despmtrs[ipmtr].value+dtime*privdata0->despmtrs[ipmtr].dot);
            CHECK_STATUS(ocsmSetValuD);
        }
    }

    /* build the perturbation */
    nbody = 0;
    status = ocsmBuild(PTRB, 0, &builtTo, &nbody, NULL);
    CHECK_STATUS(ocsmBuild);

    /* find the emodel associated with the perturbation */
    status = EG_attributeRet(privdata0->emodel, "vsp3Filename",
                             &attrType, &attrLen, NULL, NULL, &filename0);
    CHECK_STATUS(EG_attributeRet);

    for (iprim = 0; iprim < MAXPRIM; iprim++) {

        /* start at end so that you get the one created by the
           perturbation (and not the one left over from MODL) */
        for (i = PTRB->NumUdp[iprim]; i > 0; i--) {
            privdata1 = (privdata_T *)(PTRB->Udps[iprim][i].data);
            if (privdata1 != NULL) {
                if (privdata1->magic == 4433341) {
                    status = EG_attributeRet(privdata1->emodel, "vsp3Filename",
                                             &attrType, &attrLen, NULL, NULL, &filename1);
                    CHECK_STATUS(EG_attributeRet);

                    if (strcmp(filename0, filename1) == 0) break;
                }

                privdata1 = NULL;
            }
        }
        if (privdata1 != NULL) break;
    }
    if (privdata1 == NULL) {
        printf("could not find PTRB->privdata\n");
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* make sure that the base and perturbed MODELS match */
    status = EG_getTopology(privdata0->emodel, &eref, &oclass0, &mtype0,
                            data, &nbody0, &ebodys0, &senses);
    CHECK_STATUS(EG_getTopology);

    status = EG_getTopology(privdata1->emodel, &eref, &oclass1, &mtype1,
                            data, &nbody1, &ebodys1, &senses);
    CHECK_STATUS(EG_getTopology);

    if (oclass0 != oclass1) {
        printf("oclass0=%d, oclass1=%d\n",
               oclass0, oclass1);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    } else if (mtype0 != mtype1) {
        printf("mtype0=%d, mtype1=%d\n",
               mtype0, mtype1);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    } else if (nbody0 != nbody1) {
        printf("nbody0=%d, nbody1=%d\n",
               nbody0, nbody1);
        status = OCSM_INTERNAL_ERROR;
        goto cleanup;
    }

    /* compute the sensitivities of the control point locations and
           attach them to the _dots of MODL ... */
    for (ichild = 0; ichild < nbody0; ichild++) {
        ibody = privdata0->ifirst + ichild;

        /* ... Nodes first */
        status = EG_getBodyTopos(ebodys0[ichild], NULL, NODE, &nnode0, &enodes0);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_getBodyTopos(ebodys1[ichild], NULL, NODE, &nnode1, &enodes1);
        CHECK_STATUS(EG_getBodyTopos);

        if (nnode0 != nnode1) {
            printf("ichild=%d, nnode0=%d, nnode1=%d\n",
                   ichild, nnode0, nnode1);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }

        for (inode = 0; inode < nnode0; inode++) {
            jnode = 0;
            for (knode = 1; knode <= MODL->body[ibody].nnode; knode++) {
                status = EG_isEquivalent(enodes0[inode], MODL->body[ibody].node[knode].enode);
                if (status == EGADS_SUCCESS) {
                    jnode = knode;
                    break;
                }
            }

            MALLOC(rdata0,    double, 3);
            MALLOC(rdata1,    double, 3);
            MALLOC(rdata_dot, double, 3);

            status = EG_getTopology(enodes0[inode], &eref, &oclass, &mtype,
                                    rdata0, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_getTopology(enodes1[inode], &eref, &oclass, &mtype,
                                    rdata1, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            for (i = 0; i < 3; i++) {
                rdata_dot[i] = (rdata1[i] - rdata0[i]) / dtime;
            }

            status = EG_setGeometry_dot(MODL->body[ibody].node[jnode].enode, NODE, 0, NULL,
                                        rdata0, rdata_dot);
            CHECK_STATUS(EG_setGeometry_dot);

            FREE(rdata0   );
            FREE(rdata1   );
            FREE(rdata_dot);
        }

        EG_free(enodes0);
        EG_free(enodes1);

        /* ... Edges second */
        status = EG_getBodyTopos(ebodys0[ichild], NULL, EDGE, &nedge0, &eedges0);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_getBodyTopos(ebodys1[ichild], NULL, EDGE, &nedge1, &eedges1);
        CHECK_STATUS(EG_getBodyTopos);

        if (nedge0 != nedge1) {
            printf("ichild=%d, nedge0=%d, nedge1=%d\n",
                   ichild, nedge0, nedge1);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }

        /* ...... trange (even for degenerate Edges) */
        for (jedge = 1; jedge <= MODL->body[ibody].nedge; jedge++) {
            status = EG_getRange(MODL->body[ibody].edge[jedge].eedge, trange, &periodic);
            CHECK_STATUS(EG_getRange);

            trange_dot[0] = 0;
            trange_dot[1] = 0;

            status = EG_setRange_dot(MODL->body[ibody].edge[jedge].eedge, EDGE, trange, trange_dot);
            CHECK_STATUS(EG_setRange_dot);
        }

        /* ...... curves associated with Edges */
        for (iedge = 0; iedge < nedge0; iedge++) {
            jedge = 0;
            for (kedge = 1; kedge <= MODL->body[ibody].nedge; kedge++) {
                status = EG_isEquivalent(eedges0[iedge], MODL->body[ibody].edge[kedge].eedge);
                if (status == EGADS_SUCCESS) {
                    jedge = kedge;
                    break;
                }
            }

            status = EG_getTopology(eedges0[iedge], &ecurve0, &oclass, &mtype,
                                    data, &ntemp, &etemps, &senses);
            CHECK_STATUS(EG_getTopology);

            if (mtype == DEGENERATE) continue;

            status = EG_getTopology(eedges1[iedge], &ecurve1, &oclass, &mtype,
                                    data, &ntemp, &etemps, &senses);
            CHECK_STATUS(EG_getTopology);

            if (mtype == DEGENERATE) continue;

            header0 = NULL;
            header1 = NULL;

            status = EG_getGeometry(ecurve0, &oclass0, &mtype0, &eref, &header0, &rdata0);
            CHECK_STATUS(EG_getGeometry);

            status = EG_getGeometry(ecurve1, &oclass1, &mtype1, &eref, &header1, &rdata1);
            CHECK_STATUS(EG_getGeometry);

            SPLINT_CHECK_FOR_NULL(header0);
            SPLINT_CHECK_FOR_NULL(header1);
            SPLINT_CHECK_FOR_NULL(rdata0 );
            SPLINT_CHECK_FOR_NULL(rdata1 );

            if (oclass0 != CURVE || oclass1 != CURVE) {
                printf("ichild=%d, iedge=%d, oclass0=%d, oclass1=%d\n",
                       ichild, iedge, oclass0, oclass1);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            } else if (mtype0 != BSPLINE || mtype1 != BSPLINE) {
                printf("ichild=%d, iedge=%d, mtype0=%d, mtype1=%d\n",
                       ichild, iedge, mtype0, mtype1);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            } else {
                for (i = 0; i < 4; i++) {
                    if (header0[i] != header1[i]) {
                        printf("ichild=%d, iedge=%d, i=%d, header0=%d, header1=%d\n",
                               ichild, iedge, i, header0[i], header1[i]);
                        status = OCSM_INTERNAL_ERROR;
                        goto cleanup;
                    }
                }
            }

            /* store the _dots on the curve */
            MALLOC(rdata_dot, double, header0[3]+3*header0[2]);

            for (iknot = 0; iknot < header0[3]; iknot++) {
                if (rdata0[iknot] != rdata1[iknot]) {
                    printf("ichild=%d, iedge=%d, iknot=%d, rdata0=%f, rdata1=%f\n",
                           ichild, iedge, iknot, rdata0[iknot], rdata1[iknot]);
                    status = OCSM_INTERNAL_ERROR;
                    goto cleanup;
                }
                rdata_dot[iknot] = 0;
            }

            for (icp = header0[3]; icp < header0[3]+3*header0[2]; icp++) {
                rdata_dot[icp] = (rdata1[icp] - rdata0[icp]) / dtime;
            }

            status = EG_getTopology(MODL->body[ibody].edge[jedge].eedge, &ecurve, &oclass, &mtype,
                                    data, &ntemp, &etemps, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_setGeometry_dot(ecurve, CURVE, BSPLINE, header0, rdata0, rdata_dot);
            CHECK_STATUS(EG_setGeometry_dot);

            FREE(rdata_dot);

            EG_free(header0);    header0 = NULL;
            EG_free(header1);    header1 = NULL;
            EG_free(rdata0 );    rdata0  = NULL;
            EG_free(rdata1 );    rdata1  = NULL;
        }

        EG_free(eedges0);
        EG_free(eedges1);

        /* ... Faces (and their PCurves) third */
        status = EG_getBodyTopos(ebodys0[ichild], NULL, FACE, &nface0, &efaces0);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_getBodyTopos(ebodys1[ichild], NULL, FACE, &nface1, &efaces1);
        CHECK_STATUS(EG_getBodyTopos);

        if (nface0 != nface1) {
            printf("ichild=%d, nface0=%d, nface1=%d\n",
                   ichild, nface0, nface1);
            status = OCSM_INTERNAL_ERROR;
            goto cleanup;
        }

        for (iface = 0; iface < nface0; iface++) {
            jface = 0;
            for (kface = 1; kface <= MODL->body[ibody].nface; kface++) {
                status = EG_isEquivalent(efaces0[iface], MODL->body[ibody].face[kface].eface);
                if (status == EGADS_SUCCESS) {
                    jface = kface;
                    break;
                }
            }

            status = EG_getTopology(efaces0[iface], &esurf0, &oclass, &mtype,
                                    data, &nchild0, &echilds0, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_getTopology(efaces1[iface], &esurf1, &oclass, &mtype,
                                    data, &nchild1, &echilds1, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_getGeometry(esurf0, &oclass0, &mtype0, &eref, &header0, &rdata0);
            CHECK_STATUS(EG_getGeometry);

            status = EG_getGeometry(esurf1, &oclass1, &mtype1, &eref, &header1, &rdata1);
            CHECK_STATUS(EG_getGeometry);

            SPLINT_CHECK_FOR_NULL(header0);
            SPLINT_CHECK_FOR_NULL(header1);
            SPLINT_CHECK_FOR_NULL(rdata0 );
            SPLINT_CHECK_FOR_NULL(rdata1 );

            if (nchild0 != 1 || nchild1 != 1) {
                printf("ichild=%d, iface=%d, nchild0=%d, nchild1=%d\n",
                       ichild, iface, nchild0, nchild1);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            } else if (oclass0 != SURFACE || oclass1 != SURFACE) {
                printf("ichild=%d, iface=%d, oclass0=%d, oclass1=%d\n",
                       ichild, iface, oclass0, oclass1);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            } else if (mtype0 != BSPLINE || mtype1 != BSPLINE) {
                printf("ichild=%d, iface=%d, mtype0=%d, mtype1=%d\n",
                       ichild, iface, mtype0, mtype1);
                status = OCSM_INTERNAL_ERROR;
                goto cleanup;
            } else {
                for (i = 0; i < 7; i++) {
                    if (header0[i] != header1[i]) {
                        printf("ichild=%d, iface=%d, i=%d, header0=%d, header1=%d\n",
                               ichild, iface, i, header0[i], header1[i]);
                        status = OCSM_INTERNAL_ERROR;
                        goto cleanup;
                    }
                }
            }

            /* store the _dots on the surface */
            MALLOC(rdata_dot, double, header0[3]+header0[6]+3*header0[2]*header0[5]);

            for (iknot = 0; iknot < header0[3]+header0[6]; iknot++) {
                if (rdata0[iknot] != rdata1[iknot]) {
                    printf("ichild=%d, iface=%d, iknot=%d, rdata0=%f, rdata1=%f\n",
                           ichild, iface, iknot, rdata0[iknot], rdata1[iknot]);
                    status = OCSM_INTERNAL_ERROR;
                    goto cleanup;
                }
                rdata_dot[iknot] = 0;
            }

            for (icp = header0[3]+header0[6]; icp < header0[3]+header0[6]+3*header0[2]*header0[5]; icp++) {
                rdata_dot[icp] = (rdata1[icp] - rdata0[icp]) / dtime;
            }

            status = EG_getTopology(MODL->body[ibody].face[jface].eface, &esurf, &oclass, &mtype,
                                    data, &ntemp, &etemps, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_setGeometry_dot(esurf, SURFACE, BSPLINE, header0, rdata0, rdata_dot);
            CHECK_STATUS(EG_setGeometry_dot);

            FREE(rdata_dot);

            EG_free(header0);   header0 = NULL;
            EG_free(header1);   header1 = NULL;
            EG_free(rdata0 );   rdata0  = NULL;
            EG_free(rdata1 );   rdata1  = NULL;

            /* process the PCurves associated with this Face. note that
               echilds*[0] is the only loop in the Faces */

        }

        EG_free(efaces0);
        EG_free(efaces1);

        /* verify that the Body has dots */
        status = EG_hasGeometry_dot(MODL->body[ibody].ebody);
        CHECK_STATUS(EG_hasGeometry_dot);
    }

     /* remove the perturbed Body */
    status = ocsmFree(PTRB);
    CHECK_STATUS(ocsmFree);

    /* reset the user pointer back to MODL (from PTRB) */
    status = EG_setUserPointer(context, (void *)MODL);
    CHECK_STATUS(EG_setUserPointer);

    /* verify that the Body has dots */
    status = EG_hasGeometry_dot(ebody);
    CHECK_STATUS(EG_hasGeometry_dot);

evaluate:
    /* evaluate the velocities from the _dots stored in MODL */
    if (entType == OCSM_NODE) {
        inode = entIndex;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            status = EG_evaluate_dot(MODL->body[ibody].node[inode].enode, NULL, NULL, data, data_dot);
            CHECK_STATUS(EG_evaluate_dot);

            vels[3*ipnt  ] = data_dot[0];
            vels[3*ipnt+1] = data_dot[1];
            vels[3*ipnt+2] = data_dot[2];
        }
    } else if (entType == OCSM_EDGE) {
        iedge = entIndex;

        uvs_dot[0] = 0;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            status = EG_evaluate_dot(MODL->body[ibody].edge[iedge].eedge, &uvs[ipnt], uvs_dot, data, data_dot);
            CHECK_STATUS(EG_evaluate_dot);

            vels[3*ipnt  ] = data_dot[0];
            vels[3*ipnt+1] = data_dot[1];
            vels[3*ipnt+2] = data_dot[2];
        }
    } else {
        iface = entIndex;

        uvs_dot[0] = 0;
        uvs_dot[1] = 0;

        for (ipnt = 0; ipnt < npnt; ipnt++) {
            status = EG_evaluate_dot(MODL->body[ibody].face[iface].eface, &uvs[2*ipnt], uvs_dot, data, data_dot);
            CHECK_STATUS(EG_evaluate_dot);

            vels[3*ipnt  ] = data_dot[0];
            vels[3*ipnt+1] = data_dot[1];
            vels[3*ipnt+2] = data_dot[2];
        }
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   runVsp3 - run vsp and create a MODEL                               *
 *                                                                      *
 ************************************************************************
 */

static int
runVsp3(modl_T  *MODL,
        char    filename[],
        ego     *emodel,
/*@unused@*/int *NumUdp,
        udp_T   *udps)
{
    int     status = EGADS_SUCCESS;

    int     ipmtr, ptype, nrow, ncol, i;
    double  value, dot;
    char    *vsp3_root, command[1024], pname[MAX_NAME_LEN];
    FILE    *fp_vspscript;

    ROUTINE(runVsp3);

    /* --------------------------------------------------------------- */

    /* create the TeMpVsP3.vspscript file */
    fp_vspscript = fopen("TeMpVsP3.vspscript", "w");
    if (fp_vspscript == NULL) {
        snprintf(message, 1023, "could not create \"TeMpVsP3.vspscript\"");
        status = EGADS_NOTFOUND;
        goto cleanup;
    }

    /* .vspscript prolog */
    fprintf(fp_vspscript, "void main()\n");
    fprintf(fp_vspscript, "{\n");
    fprintf(fp_vspscript, "    string vspname = \"%s\";\n", filename);
    fprintf(fp_vspscript, "    ReadVSPFile(vspname);\n\n");

    /* update vsp UserParms:ESP_Group:* from DESPMTRs */
    fprintf(fp_vspscript, "    string user_ctr = FindContainer(\"UserParms\", 0);\n");
    fprintf(fp_vspscript, "    SilenceErrors();\n");
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        status = ocsmGetPmtr(MODL, ipmtr, &ptype, &nrow, &ncol, pname);
        CHECK_STATUS(ocsmGetPmtr);

        if (ptype == OCSM_DESPMTR && nrow == 1 && ncol == 1) {
            status = ocsmGetValu(MODL, ipmtr, 1, 1, &value, &dot);
            CHECK_STATUS(ocsmGetValu);

            for (i = 0; i < STRLEN(pname); i++) {
                if (pname[i] == ':') pname[i] = '.';    // OpenCSM uses ':' but VSP uses '.'
            }

            fprintf(fp_vspscript, "    SetParmVal(FindParm(user_ctr, \"%s\", \"ESP_Group\"), %20.14e);\n", pname, value);
        }
    }
    fprintf(fp_vspscript, "    PrintOnErrors();\n\n");

    /* .vspscript epilog */
    fprintf(fp_vspscript, "    string veh_id = GetVehicleID();\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"SplitSurfs\", \"STEPSettings\"), 1);\n\n");

    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelID\", \"STEPSettings\"), 1);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelName\", \"STEPSettings\"), 1);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelSurfNo\", \"STEPSettings\"), 1);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"LabelDelim\", \"STEPSettings\"), DELIM_COMMA);\n\n");

    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"SplitSubSurfs\", \"STEPSettings\"), 0);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"MergePoints\", \"STEPSettings\"), 0);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"ToCubic\", \"STEPSettings\"), 0);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"TrimTE\", \"STEPSettings\"), 0);\n");
    fprintf(fp_vspscript, "    SetParmVal(FindParm(veh_id, \"ExportPropMainSurf\", \"STEPSettings\"), 0);\n\n");

    fprintf(fp_vspscript, "    string stpname = \"TeMpVsP3.stp\";\n");
    fprintf(fp_vspscript, "    ExportFile(stpname, SET_ALL, EXPORT_STEP);\n");
    fprintf(fp_vspscript, "}\n");

    fclose(fp_vspscript);

    /* exeucute vspscript */
    vsp3_root = getenv("VSP3_ROOT");
    if (vsp3_root != NULL) {
        snprintf(command, 1023, "%s%cvspscript -script TeMpVsP3.vspscript", vsp3_root, SLASH);
    } else {
        snprintf(command, 1023, "vspscript -script TeMpVsP3.vspscript");
    }

    printf("\n====================\nRunning: %s\n", command);
    system(command);
    SLEEP(1000);
    printf("vspscript has completed\n====================\n\n");

    /* process the .stp file */
    status = processStepFile(MODL->context, "TeMpVsP3.stp", emodel);
    CHECK_STATUS(processStepFile);

    /* clean up temporary files */
    if (KEEPTEMPS(0) == 0) {
        remove("TeMpVsP3.stp");
        remove("TeMpVsP3.vspscript");
    }

cleanup:
    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   processStepFile - read .stp file and create Model with the Bodys   *
 *                                                                      *
 ************************************************************************
 */

static int
processStepFile(ego    context,         /* (in)  EGADS context */
                char   filename[],      /* (in)  filename */
                ego    *emodel)         /* (out) Model */
{
    int     status = EGADS_SUCCESS;

    int     attrType, attrLen, oclass, mtype, nchild, ichild, *senses;
    int     i, j, igeomID, igeomName, isurfNum, nbody=0, nface=0, nnn, *sss;
    int     nchild2, nface2, periodic, *header, nknot, imax, jmax;
    CINT    *tempIlist;
    double  data[18], uvrange[4], bbox[6], *rdata, mprops[14];
    double  xmin=+HUGEQ, xmax=-HUGEQ, ymin=+HUGEQ, ymax=-HUGEQ, zmin=+HUGEQ, zmax=-HUGEQ, eps06=EPS06;
    CDOUBLE *tempRlist;
    char    oldGeomID[80], geomID[80], oldGeomName[80], geomName[80], oldSurfNum[80], surfNum[80], bodyName[80];
    CCHAR   *tempClist;
    ego     emodel1, eref, *echilds, *efaces=NULL, *newBodys=NULL, newModel, *eee, eshell;
    ego     *efaces2, esurf, *echilds2;

    ROUTINE(processStepFile);

    /* --------------------------------------------------------------- */

    /* read the .stp file (created by exporting from OpenVSP) */
    status = EG_loadModel(context, 0, filename, &emodel1);
    CHECK_STATUS(EG_loadModel);

    /* get the Bodys contained in emodel1 */
    status = EG_getTopology(emodel1, &eref, &oclass, &mtype, data,
                            &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    printf("There are %d children in emodel1\n", nchild);

    MALLOC(efaces,   ego, nchild);
    MALLOC(newBodys, ego, nchild);

    nbody          = 0;
    oldGeomID[0]   = '\0';
    oldGeomName[0] = '\0';
    oldSurfNum[0]  = '\0';

    /* get the Name of each Body */
    for (ichild = 0; ichild < nchild; ichild++) {

        /* if the Body does not have a Name, skip it */
        status = EG_attributeRet(echilds[ichild], "Name", &attrType, &attrLen,
                                 &tempIlist, &tempRlist, &tempClist);
        if (status != EGADS_SUCCESS || attrType != ATTRSTRING) {
            printf("Skipping   Child %3d (does not have Name)\n", ichild+1);
            continue;
        } else {
            /* if the Body has no area, skip it */
            status = EG_getMassProperties(echilds[ichild], mprops);
            CHECK_STATUS(EG_getMassProperties);

            if (mprops[1] < EPS12) {
                printf("Skipping   Child %3d (area=%12.6ef)\n", ichild+1, mprops[1]);
                continue;
            }

            printf("Processing Child %3d (%s), nface=%d\n", ichild+1, tempClist, nface);
        }


        /* extract the geomID from the Name */
        igeomID   = 0;
        geomID[0] = '\0';

        for (i = 0; i < attrLen; i++) {
            if (tempClist[i] == ',') {
                i++;
                break;
            }

            geomID[igeomID  ] = tempClist[i];
            geomID[igeomID+1] = '\0';
            igeomID++;
        }

        /* extract the geomName from the Name */
        igeomName   = 0;
        geomName[0] = '\0';

        for ( ; i < attrLen; i++) {
            if (tempClist[i] == ',') {
                i++;
                break;
            }

            if (tempClist[i] == ' ') {
                if (igeomName == 0) {
                    continue;
                } else {
                    geomName[igeomName  ] = '_';
                    geomName[igeomName+1] = '\0';
                    igeomName++;
                }
            } else {
                geomName[igeomName  ] = tempClist[i];
                geomName[igeomName+1] = '\0';
                igeomName++;
            }
        }

        /* extract the surfNum from the Name */
        isurfNum   = 0;
        surfNum[0] = '\0';

        for ( ; i < attrLen; i++) {
            if (tempClist[i] == ',') {
                i++;
                break;
            }

            if (tempClist[i] == ' ') {
                if (isurfNum == 0) {
                    continue;
                } else {
                    surfNum[isurfNum  ] = '_';
                    surfNum[isurfNum+1] = '\0';
                    isurfNum++;
                }
            } else {
                surfNum[isurfNum  ] = tempClist[i];
                surfNum[isurfNum+1] = '\0';
                isurfNum++;
            }
        }

        /* if this has a different component name than the previous Body,
           make the SolidBody (or SheetBody if it is open) from the FaceBodys
           that have been previously processed (if any) */
        if (strcmp(geomID, oldGeomID) != 0 || strcmp(surfNum, oldSurfNum) != 0) {
            if (nface > 0) {
                eps06 = EPS06 * MAX(MAX(xmax-xmin,ymax-ymin),MAX(zmax-zmin,100));

                if (nface > 1) {
                    status = EG_sewFaces(nface, efaces, eps06, 1, &newModel);
                    CHECK_STATUS(EG_sewFaces);

                    status = EG_getTopology(newModel, &eref, &oclass, &mtype, data,
                                            &nnn, &eee, &sss);
                    CHECK_STATUS(EG_getTopology);

                    if (nnn != 1) {
                        printf("EG_sewFaces(nface=%d) generated %d Bodys.  Only using first Body\n", nface, nnn);
                        ocsmPrintEgo(newModel);
                    }

                    status = EG_copyObject(eee[0], NULL, &(newBodys[nbody++]));
                    CHECK_STATUS(EG_copyObject);

                    status = EG_deleteObject(newModel);
                    CHECK_STATUS(EG_deleteObject);
                } else {
                    status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL,
                                             1, efaces, NULL, &eshell);
                    CHECK_STATUS(EG_makeTopology);

                    status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL,
                                             1, &eshell, NULL, &(newBodys[nbody++]));
                    CHECK_STATUS(EG_makeTopology);
                }

                snprintf(bodyName, 79, "%s.%s:%d", oldGeomName, oldSurfNum, nbody);

                status = EG_attributeAdd(newBodys[nbody-1], "_name", ATTRSTRING, STRLEN(bodyName),
                                         NULL, NULL, bodyName);
                CHECK_STATUS(EG_attributeAdd);

                status = EG_attributeAdd(newBodys[nbody-1], "_vspBody", ATTRINT, 1, &nbody, NULL, NULL);
                CHECK_STATUS(EG_attributeAdd);

                printf("   Made Body %3d (%s) with %d Faces\n", nbody-1, bodyName, nface);

                nface = 0;
                xmin  = +HUGEQ;
                xmax  = -HUGEQ;
                ymin  = +HUGEQ;
                ymax  = -HUGEQ;
                zmin  = +HUGEQ;
                zmax  = -HUGEQ;
            }
            strcpy(oldGeomID,   geomID  );
            strcpy(oldGeomName, geomName);
            strcpy(oldSurfNum,  surfNum );
        }

        /* since older versions of OpenVSP can create knot vectors with
           jumps near the end, we will extract the Surface from the Face,
           adjust the knot vector, and remake the Surface and Face */
        status = EG_getBodyTopos(echilds[ichild], NULL, FACE, &nface2, &efaces2);
        CHECK_STATUS(EG_getBodyTopos);

        status = EG_getTopology(efaces2[0], &esurf, &oclass, &mtype, data,
                                &nchild2, &echilds2, &senses);
        CHECK_STATUS(EG_getTopology);

        EG_free(efaces2);

        status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &header, &rdata);
        CHECK_STATUS(EG_getGeometry);

        /* fix jumps in knot vectors (fixes error in OpenVSP) */
        for (i = 1; i < header[3]; i++) {
            if (rdata[i]-rdata[i-1] > 1.01) {
                for (j = i; j < header[3]; j++) {
                    rdata[j] -= 1;
                }
                i--;
            }
        }

        for (i = 1; i < header[6]; i++) {
            if (rdata[header[3]+i]-rdata[header[3]+i-1] > 1.01) {
                for (j = i; j < header[6]; j++) {
                    rdata[header[3]+j] -= 1;
                }
                i--;
            }
        }

        /* OpenVSP can generate degeneracies at the boundaries of a Face.  when
           these are detected, adjust the internal control point to the average
           of its neighbors ... */
        nknot = header[3] + header[6];
        imax  = header[2];
        jmax  = header[5];

#define XCP(I,J)  rdata[nknot+3*((I)+imax*(J))  ]
#define YCP(I,J)  rdata[nknot+3*((I)+imax*(J))+1]
#define ZCP(I,J)  rdata[nknot+3*((I)+imax*(J))+2]

        /* ... imin boundary */
        if (imax > 2) {
            i = 1;
            for (j = 0; j < jmax; j++) {
                if (fabs(XCP(i,j)-XCP(i-1,j)) < EPS12 &&
                    fabs(YCP(i,j)-YCP(i-1,j)) < EPS12 &&
                    fabs(ZCP(i,j)-ZCP(i-1,j)) < EPS12   ) {
                    XCP(i,j) = (XCP(i+1,j) + XCP(i-1,j)) / 2;
                    YCP(i,j) = (YCP(i+1,j) + YCP(i-1,j)) / 2;
                    ZCP(i,j) = (ZCP(i+1,j) + ZCP(i-1,j)) / 2;
                }
            }
        }

        /* ... imax boundary */
        if (imax > 2) {
            i = imax - 2;
            for (j = 0; j < jmax; j++) {
                if (fabs(XCP(i,j)-XCP(i+1,j)) < EPS12 &&
                    fabs(YCP(i,j)-YCP(i+1,j)) < EPS12 &&
                    fabs(ZCP(i,j)-ZCP(i+1,j)) < EPS12   ) {
                    XCP(i,j) = (XCP(i+1,j) + XCP(i-1,j)) / 2;
                    YCP(i,j) = (YCP(i+1,j) + YCP(i-1,j)) / 2;
                    ZCP(i,j) = (ZCP(i+1,j) + ZCP(i-1,j)) / 2;
                }
            }
        }

        /* ... jmin boundary */
        if (jmax > 2) {
            j = 1;
            for (i = 0; i < imax; i++) {
                if (fabs(XCP(i,j)-XCP(i,j-1)) < EPS12 &&
                    fabs(YCP(i,j)-YCP(i,j-1)) < EPS12 &&
                    fabs(ZCP(i,j)-ZCP(i,j-1)) < EPS12   ) {
                    XCP(i,j) = (XCP(i,j+1) + XCP(i,j-1)) / 2;
                    YCP(i,j) = (YCP(i,j+1) + YCP(i,j-1)) / 2;
                    ZCP(i,j) = (ZCP(i,j+1) + ZCP(i,j-1)) / 2;
                }
            }
        }

        /* ... jmax boundary */
        if (jmax > 2) {
            j = jmax - 2;
            for (i = 0; i < imax; i++) {
                if (fabs(XCP(i,j)-XCP(i,j+1)) < EPS12 &&
                    fabs(YCP(i,j)-YCP(i,j+1)) < EPS12 &&
                    fabs(ZCP(i,j)-ZCP(i,j+1)) < EPS12   ) {
                    XCP(i,j) = (XCP(i,j+1) + XCP(i,j-1)) / 2;
                    YCP(i,j) = (YCP(i,j+1) + YCP(i,j-1)) / 2;
                    ZCP(i,j) = (ZCP(i,j+1) + ZCP(i,j-1)) / 2;
                }
            }
        }

#undef XCP
#undef YCP
#undef ZCP
        
        status = EG_makeGeometry(context, oclass, mtype, NULL, header, rdata, &esurf);
        CHECK_STATUS(EG_makeGeometry);

        EG_free(header);
        EG_free(rdata );

        status = EG_getRange(esurf, uvrange, &periodic);
        CHECK_STATUS(EG_getRange);

        status = EG_makeFace(esurf, SFORWARD, uvrange, &(efaces[nface]));
        CHECK_STATUS(EG_makeFace);

        status = EG_attributeAdd(efaces[nface], "_vspID", ATTRSTRING, strlen(oldGeomID), NULL, NULL, oldGeomID);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_attributeAdd(efaces[nface], "_vspFace", ATTRINT, 1, &nface, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_getBoundingBox(efaces[nface], bbox);
        CHECK_STATUS(EG_getBoundingBox);

        xmin = MIN(xmin, bbox[0]);
        ymin = MIN(ymin, bbox[1]);
        zmin = MIN(zmin, bbox[2]);
        xmax = MAX(xmax, bbox[3]);
        ymax = MAX(ymax, bbox[4]);
        zmax = MAX(zmax, bbox[5]);

        nface++;
    } /* next Body from the .stp file */

    /* make Sheet/SolidBodys from last FaceBody(s) */
    if (nface > 0) {
        eps06 = EPS06 * MAX(MAX(xmax-xmin,ymax-ymin),MAX(zmax-zmin,100));

        if (nface > 1) {
            status = EG_sewFaces(nface, efaces, eps06, 1, &newModel);
            CHECK_STATUS(EG_sewFaces);

            status = EG_getTopology(newModel, &eref, &oclass, &mtype, data,
                                    &nnn, &eee, &sss);
            CHECK_STATUS(EG_getTopology);

            if (nnn != 1) {
                printf("EG_sewFaces(nface=%d) generated %d Bodys.  Only using first Body\n", nface, nnn);
                ocsmPrintEgo(newModel);
            }

            status = EG_copyObject(eee[0], NULL, &(newBodys[nbody++]));
            CHECK_STATUS(EG_copyObject);

            status = EG_deleteObject(newModel);
            CHECK_STATUS(EG_deleteObject);

        } else {
            status = EG_makeTopology(context, NULL, SHELL, OPEN, NULL,
                                     1, efaces, NULL, &eshell);
            CHECK_STATUS(EG_makeTopology);

            status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL,
                                     1, &eshell, NULL, &(newBodys[nbody++]));
            CHECK_STATUS(EG_makeTopology);
        }

        snprintf(bodyName, 79, "%s.%s:%d", oldGeomName, oldSurfNum, nbody);

        status = EG_attributeAdd(newBodys[nbody-1], "_name", ATTRSTRING, STRLEN(bodyName),
                                 NULL, NULL, bodyName);
        CHECK_STATUS(EG_attributeAdd);

        status = EG_attributeAdd(newBodys[nbody-1], "_vspBody", ATTRINT, 1, &nbody, NULL, NULL);
        CHECK_STATUS(EG_attributeAdd);

        printf("   Made Body %d (%s) with %d Faces\n", nbody-1, bodyName, nface);
    }

    /* remove the MODEL read from the .stp file and all its children */
    status = EG_deleteObject(emodel1);
    CHECK_STATUS(EG_deleteObject);

    /* make a Model of the SolidBodys to return to OpenCSM */
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, nbody, newBodys, NULL, emodel);
    CHECK_STATUS(EG_makeTopology);

cleanup:
    FREE(efaces  );
    FREE(newBodys);

    return status;
}


/*
 ************************************************************************
 *                                                                      *
 *   freePrivateData - free private data                                *
 *                                                                      *
 ************************************************************************
 */

static int
freePrivateData(void  *data)            /* (in)  pointer to private data */
{
    int    status = EGADS_SUCCESS;

    privdata_T *privdata = (privdata_T *)data;

    ROUTINE(freePrivateData);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("freePrivateData(%llx)\n", (long long)data);
#endif

    if (privdata->emodel != NULL) {
        status = EG_deleteObject(privdata->emodel);
        CHECK_STATUS(EG_deleteObject);

        privdata->emodel = NULL;
    }

    FREE(privdata->despmtrs);

    FREE(privdata);

cleanup:
    return status;
}



/*
 ************************************************************************
 *                                                                      *
 *   copyPrivateData - copy private data                                *
 *                                                                      *
 ************************************************************************
 */

static int
copyPrivateData(
      /*@null@*/void  *src,             /* (in)  pointer to source private data */
                void  **tgt)            /* (in)  pointer to target private data */
{
    int    status = EGADS_SUCCESS;

    int        ipmtr;
    despmtr_T  *despmtrs=NULL;
    privdata_T *SRC, *TGT;

    ROUTINE(copyPrivateData);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("copyPrivateData(%llx)\n", (long long)src);
#endif

    *tgt = NULL;

    if (src == NULL) goto cleanup;

    SRC = (privdata_T *)src;

    MALLOC(*tgt, privdata_T *, 1);

    TGT = (privdata_T *)(*tgt);

    TGT->magic  = SRC->magic;
    TGT->ifirst = SRC->ifirst;
    TGT->ilast  = SRC->ilast;
    TGT->npmtr  = SRC->npmtr;

    MALLOC(despmtrs, despmtr_T, SRC->npmtr);

    for (ipmtr = 0; ipmtr < SRC->npmtr; ipmtr++) {
        despmtrs[ipmtr].ipmtr = SRC->despmtrs[ipmtr].ipmtr;
        despmtrs[ipmtr].irow  = SRC->despmtrs[ipmtr].irow;
        despmtrs[ipmtr].icol  = SRC->despmtrs[ipmtr].icol;
        despmtrs[ipmtr].value = SRC->despmtrs[ipmtr].value;
        despmtrs[ipmtr].dot   = SRC->despmtrs[ipmtr].dot;
    }

    TGT->despmtrs = despmtrs;

    status = EG_copyObject(SRC->emodel, NULL, &(TGT->emodel));
    CHECK_STATUS(EG_copyObject);

cleanup:
    return status;
}
