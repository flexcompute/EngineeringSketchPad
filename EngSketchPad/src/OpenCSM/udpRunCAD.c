/*
 ************************************************************************
 *                                                                      *
 * udpRunCAD -- runs external CAD and returns Body with dots            *
 *                                                                      *
 *            Written by John Dannenhoffer @ Syracuse University        *
 *            Patterned after code written by Bob Haimes  @ MIT         *
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

#define NUMUDPARGS       3
#define NUMUDPINPUTBODYS 0
#include "udpUtilities.h"

/* shorthands for accessing argument values and velocities */
#define TOOLNAME(IUDP  )  ((char   *) (udps[IUDP].arg[0].val))
#define CASENAME(IUDP  )  ((char   *) (udps[IUDP].arg[1].val))
#define PMTRS(   IUDP,I)  ((double *) (udps[IUDP].arg[2].val))[I]

/* data about possible arguments */
static char  *argNames[NUMUDPARGS] = {"toolname", "casename", "pmtrs",     };
static int    argTypes[NUMUDPARGS] = {ATTRSTRING, ATTRFILE,   ATTRREALSEN, };
static int    argIdefs[NUMUDPARGS] = {0,          0,          0,           };
static double argDdefs[NUMUDPARGS] = {0.,         0.,         0.,          };

/* get utility routines: udpErrorStr, udpInitialize, udpReset, udpSet,
                         udpGet, udpVel, udpClean, udpMesh */
#include "udpUtilities.c"

/* prototype for function defined below */

#include "OpenCSM.h"
#include "egads_dot.h"
#include <assert.h>


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

    int     outLevel, oclass, mtype, nchild, *senses, ipmtr, nface, iface, seqnum;
    int     judp, match, i;
    double  data[18];
    char    csmfile[257], pmtrfile[257], errfile[257], stpfile[257], dotsfile[257];
    char    tempfile[257], command[257];
    char    *message=NULL;
    ego     eref, emodel, *echilds, *efaces=NULL;
    void    *modl, *temp;
    FILE    *fp_temp;
    modl_T  *MODL;
    udp_T   *udps = *Udps;

    ROUTINE(udpExecute);

    /* --------------------------------------------------------------- */

#ifdef DEBUG
    printf("udpExecute(context=%llx)\n", (long long)context);
    printf("toolname(0) = %s\n", TOOLNAME(0));
    printf("casename(0) = %s\n", CASENAME(0));
    printf("pmtrs(   0) =");
    for (i = 0; i < udps[0].arg[2].size; i++) {
        printf(" %10.5f", PMTRS(0,i));
    }
    printf("\n");
#endif

    outLevel = ocsmSetOutLevel(-1);

    SPRINT0(2, "\nEnter udpRunCAD/sensitivity\n");

    /* default return values */
    *ebody  = NULL;
    *nMesh  = 0;
    *string = NULL;

    MALLOC(message, char, 100);
    message[0] = '\0';

    /* check arguments */
    if (STRLEN(TOOLNAME(0)) <= 0) {
        snprintf(message, 100, "\"toolname\" must be specified");
        status = EGADS_RANGERR;
        goto cleanup;

    } else if (STRLEN(CASENAME(0)) <= 0) {
        snprintf(message, 100, "\"casename\" must be specified");
        status = EGADS_RANGERR;
        goto cleanup;

    }

    /* if the pmtrs agree with the pmtrs from a previous call, we can
       just re-use the previous Body and skip the rest of this routine */
    for (judp = 1; judp <= numUdp; judp++) {
        if (udps[0].arg[2].size != udps[judp].arg[2].size) continue;

        match = 1;
        for (i = 0; i < udps[0].arg[2].size; i++) {
            if (fabs(PMTRS(0,i)-PMTRS(judp,i)) > EPS12) {
                match = 0;
                break;
            }
        }
#ifdef DEBUG
        printf("judp=%3d, match=%d\n", judp, match);
#endif
        if (match == 1) {
            *ebody = udps[judp].ebody;
            goto cleanup;
        }
    }

    /* cache copy of arguments for future use */
    status = cacheUdp(NULL);
    CHECK_STATUS(cacheUdp);

#ifdef DEBUG
    printf("toolname(%d) = %s\n", numUdp, TOOLNAME(numUdp));
    printf("casename(%d) = %s\n", numUdp, CASENAME(numUdp));
#endif

    /* get a (semi-)unique sequence number */
    seqnum = time(NULL) % 1000;
    SPRINT1(2, "seqnum=%04d", seqnum);

    snprintf(csmfile,  256, "%s.csm",       CASENAME(0)        );
    snprintf(pmtrfile, 256, "%s.pmtr",      CASENAME(0)        );
    snprintf(errfile,  256, "%s.err",       CASENAME(0)        );
    snprintf(stpfile,  256, "%s_%04d.stp",  CASENAME(0), seqnum);
    snprintf(dotsfile, 256, "%s_%04d.dots", CASENAME(0), seqnum);

#ifdef DEBUG
    printf("csmfile =%s\n", csmfile );
    printf("pmtrfile=%s\n", pmtrfile);
    printf("errfile =%s\n", errfile );
    printf("stpfile =%s\n", stpfile );
    printf("dotsfile=%s\n", dotsfile);
#endif

    /* make sure the .csm file exists */
    fp_temp = fopen(csmfile, "r");
    if (fp_temp == NULL) {
        printf("\"%s\" does not exist\n", csmfile);
        status = -999;
        goto cleanup;
    } else {
        fclose(fp_temp);
    }

    /* create the .pmtr file from the current scalar DESPMTRs */
    fp_temp = fopen(pmtrfile, "w");
    if (fp_temp == NULL) {
        printf("\"%s\" could not be created\n", pmtrfile);
        status = -998;
        goto cleanup;
    }

    status = EG_getUserPointer(context, (void**)(&modl));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)modl;

    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR &&
            MODL->pmtr[ipmtr].nrow == 1            &&
            MODL->pmtr[ipmtr].ncol == 1              ) {
            fprintf(fp_temp, "%-30s %16.8f\n", MODL->pmtr[ipmtr].name, MODL->pmtr[ipmtr].value[0]);
        }
    }

    fclose(fp_temp);

    /* remove output files if they exist */
    remove(errfile );
    remove(stpfile );
    remove(dotsfile);

    /* execute the tool */
    snprintf(command, 256, "miniESP %s -outLevel %d", CASENAME(0), outLevel);
#ifdef DEBUG
    printf("command=%s=\n", command);
#endif
    system(command);

    /* check if errfile was written */
    fp_temp = fopen(errfile, "r");
    if (fp_temp != NULL) {
        temp = fgets(message, 99, fp_temp);
        if (temp == NULL) {
            snprintf(message, 100, "error checking errfile");
        }

        fclose(fp_temp);

        status = OCSM_UDP_ERROR1;
        goto cleanup;
    }

    /* make sure stpfile and dots file were written */
    snprintf(tempfile, 256, "%s.stp", CASENAME(0));
    fp_temp = fopen(tempfile, "r");
    if (fp_temp == NULL) {
        snprintf(message, 100, "\"%s\" was not created\n", stpfile);
        status = EGADS_NOTFOUND;
        goto cleanup;
    } else {
        fclose(fp_temp);

        rename(tempfile, stpfile);
    }

    snprintf(tempfile, 256, "%s.dots", CASENAME(0));
    fp_temp = fopen(tempfile, "r");
    if (fp_temp == NULL) {
        snprintf(message, 100, "\"%s\" was not created\n", dotsfile);
        status = EGADS_NOTFOUND;
        goto cleanup;
    } else {
        fclose(fp_temp);

        rename(tempfile, dotsfile);
    }

    /* load the .stp file */
    status = EG_loadModel(context, 0, stpfile, &emodel);
    CHECK_STATUS(EG_loadModel);

    /* extract the one Body from emodel */
    status = EG_getTopology(emodel, &eref, &oclass, &mtype,
                            data, &nchild, &echilds, &senses);
    CHECK_STATUS(EG_getTopology);

    if (nchild != 1) {
        snprintf(message, 100, "\"%s\" contained %d Bodys (not 1)", stpfile, nchild);
        status = EGADS_INDEXERR;
        goto cleanup;
    }

    status = EG_copyObject(echilds[0], NULL, ebody);
    CHECK_STATUS(EG_copyObject);

    status = EG_deleteObject(emodel);
    CHECK_STATUS(EG_deleteObject);

    /* put the __fromCAD__ attribute on the Body so that the velocities are
       gotten via the _dot routines */
    SPLINT_CHECK_FOR_NULL(*ebody);

    status = EG_attributeAdd(*ebody, "__fromCAD__", ATTRINT, 1, &seqnum, NULL, NULL);
    CHECK_STATUS(EG_attributeAdd);

    /* put the __trimmed__ attribute on all the Faces so that the Edge
       sensitivities get propagated into the interior of the Face */
    status = EG_getBodyTopos(*ebody, NULL, FACE, &nface, &efaces);
    CHECK_STATUS(EG_getBodyTopos);

    for (iface = 0; iface < nface; iface++) {
        SPLINT_CHECK_FOR_NULL(efaces);

        status = EG_attributeAdd(efaces[iface], "__trimmed__", ATTRSTRING, 0, NULL, NULL, "yes");
        CHECK_STATUS(EG_attributeAdd);
    }

    EG_free(efaces);  efaces = NULL;

    /* set the output value */

    /* the copy of the Body that was annotated is returned */
    udps[numUdp].ebody = *ebody;

    SPRINT0(2, "\nExit  udpRunCAD/sensitivity\n");

cleanup:
    if (efaces != NULL) EG_free(efaces);

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
    int     status = EGADS_SUCCESS;

    int iudp, judp;

    int      outLevel, ibody, jbody, count, ipmtr,jpmtr, irow, icol, seqnum;
    int      inode, jnode, nnode, iedge, jedge, nedge, iface, jface, nface;
    int      oclass, mtype, nrdata, i, attrType, attrLen, idata[7];
    int      nchild, *senses;
    CINT     *tempIlist;
    double   pmtrValue, data[18], dotmax;
    double   *rdata=NULL, *rdata_dot=NULL, trange[2], trange_dot[2];
    CDOUBLE  *tempRlist;
    char     dotsfile[257], pmtrName[80], theName[80], edgeType[80], faceType[80];
    char     *message=NULL;
    CCHAR    *tempClist;
    ego      context, ecurve, esurface, *echilds;
    void     *modl;
    modl_T   *MODL;
    FILE     *fp;

    ROUTINE(udpSensitivity);

    /* --------------------------------------------------------------- */

    MALLOC(message, char, 100);
    message[0] = '\0';

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

    outLevel = ocsmSetOutLevel(-1);

    status = EG_attributeRet(ebody, "__fromCAD__", &attrType, &attrLen,
                             &tempIlist, &tempRlist, &tempClist);
    CHECK_STATUS(EG_attributeRet);

    seqnum = tempIlist[0];

    SPRINT1(2, "\nEnter udpRunCAD/sensitivity(seqnum=%04d)\n", seqnum);

    /* open the .dots file */
    snprintf(dotsfile, 255, "%s_%04d.dots", CASENAME(iudp), seqnum);

    fp = fopen(dotsfile, "r");
    if (fp == NULL) {
        snprintf(message, 100, "could not open \"%s\"", dotsfile);
        status = -999;
        goto cleanup;
    }

    /* get the MODL associated with ebody */
    status = EG_getContext(ebody, &context);
    CHECK_STATUS(EG_getContext);

    status = EG_getUserPointer(context, (void**)(&modl));
    CHECK_STATUS(EG_getUserPointer);

    MODL = (modl_T *)modl;

    ibody = 0;
    for (jbody = MODL->nbody; jbody > 0; jbody++) {
        if (MODL->body[jbody].ebody == ebody) {
            ibody = jbody;
            break;
        }
    }

    if (ibody == 0) {
        snprintf(message, 100, "could not find ibody");
        status = -998;
        goto cleanup;
    }
    SPRINT2(1, "in udfRunCAD/udpSensitivity: MODL->body[ibody=%d].ebody=%llx", ibody, (long long)MODL->body[ibody].ebody);

    /* make sure the number of Nodes, Edges, and Faces match the dots file */
    count = fscanf(fp, "%d %d %d", &nnode, &nedge, &nface);
    assert (count == 3);

    SPRINT3(2, "nnode=%d, nedge=%d, nface=%d", nnode, nedge, nface);

    if (MODL->body[ibody].nnode != nnode ||
        MODL->body[ibody].nedge != nedge ||
        MODL->body[ibody].nface != nface   ) {
        printf("MODL->body[%d].nnode=%5d, nnode=%5d\n", ibody, MODL->body[ibody].nnode, nnode);
        printf("MODL->body[%d].nnode=%5d, nnode=%5d\n", ibody, MODL->body[ibody].nnode, nnode);
        printf("MODL->body[%d].nnode=%5d, nnode=%5d\n", ibody, MODL->body[ibody].nnode, nnode);
        snprintf(message, 100, "mismatch betwen .stp and .dots");
        status = -997;
        goto cleanup;
    }

    /* initialize the dots on the Body to zero */
    status = EG_zeroGeometry_dot(ebody);
    CHECK_STATUS(EG_zeroGeometry_dot);

    /* read each of the DESPMTR sections of the file */
    while (1) {
        count = fscanf(fp, "%s %d %d %lf", pmtrName, &irow, &icol, &pmtrValue);
        if (count != 4) break;

        SPRINT3(1, "\n**> reading %s[%d,%d] section of .dots file", pmtrName, irow, icol);

        /* determine if this section of the file matches a DESPMTR with a
           non-zero dot */
        ipmtr = 0;
        for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
            if (MODL->pmtr[jpmtr].type != OCSM_DESPMTR) continue;

            if (strcmp(MODL->pmtr[jpmtr].name, pmtrName) == 0) {
                if (MODL->pmtr[jpmtr].dot[0] != 0) {
                    ipmtr = jpmtr;
                }
                break;
            }
        }

        for (inode = 1; inode <= nnode; inode++) {
            count = fscanf(fp, "%*s %d", &jnode);
            assert (count == 1);

            snprintf(theName, 79, "Node %d", jnode);

            for (jnode = 1; jnode <= nnode; jnode++) {
                status = EG_attributeRet(MODL->body[ibody].node[jnode].enode, "Name",
                                         &attrType, &attrLen, &tempIlist, &tempRlist, &tempClist);
                CHECK_STATUS(EG_attributeRet);

                if (strcmp(theName, tempClist) != 0) continue;

                SPRINT2x(2, "inode=%3d (in dots file) maps to jnode=%3d (in ebody)", inode, jnode);

                MALLOC(rdata,     double, 3);
                MALLOC(rdata_dot, double, 3);

                count = fscanf(fp, "%lf %lf %lf", &rdata[0],     &rdata[1],     &rdata[2]);
                assert (count == 3);

                count = fscanf(fp, "%lf %lf %lf", &rdata_dot[0], &rdata_dot[1], &rdata_dot[2]);
                assert (count == 3);

                dotmax = 0;
                for (i = 0; i < 3; i++) {
                    if (fabs(rdata_dot[i]) > dotmax) dotmax = fabs(rdata_dot[i]);
                }

                if (ipmtr > 0) {
                    SPRINT1x(2, "   dotmax=%10.5f", dotmax);

                    for (i = 0; i < 3; i++) {
                        rdata_dot[i] *= MODL->pmtr[ipmtr].dot[0];
                    }

                    status = EG_setGeometry_dot(MODL->body[ibody].node[jnode].enode, NODE, 0,
                                                NULL, rdata, rdata_dot);
                    CHECK_STATUS(EG_setGeometry_dot);
                }

                FREE(rdata    );
                FREE(rdata_dot);

                SPRINT0(2, " ");

                break;
            }
        }

        for (iedge = 1; iedge <= nedge; iedge++) {
            count = fscanf(fp, "%*s %d %s", &jedge, edgeType);
            assert (count == 2);

            snprintf(theName, 79, "Edge %d", jedge);

            if (strcmp(edgeType, "DEGENERATE") == 0) {
                count = fscanf(fp, "%lf %lf", &trange[0],     &trange[1]);
                assert (count == 2);

                count = fscanf(fp, "%lf %lf", &trange_dot[0], &trange_dot[1]);
                assert (count == 2);

                if (ipmtr > 0) {
                    status = EG_setRange_dot(MODL->body[ibody].edge[jedge].eedge, EDGE, trange, trange_dot);
                    CHECK_STATUS(EG_setRange_dot);
                }
                continue;
            }

            for (jedge = 1; jedge <= nedge; jedge++) {
                status = EG_attributeRet(MODL->body[ibody].edge[jedge].eedge, "Name",
                                         &attrType, &attrLen, &tempIlist, &tempRlist, &tempClist);
                if (status != EGADS_SUCCESS) continue;     // skip degenerate Edges

                if (strcmp(theName, tempClist) != 0)  continue;

                status = EG_getTopology(MODL->body[ibody].edge[jedge].eedge, &ecurve, &oclass, &mtype,
                                        data, &nchild, &echilds, &senses);
                CHECK_STATUS(EG_getTopology);

                SPRINT2x(2, "iedge=%3d (in dots file) maps to jedge=%3d (in ebody)", iedge, jedge);

                count = fscanf(fp, "%lf %lf", &trange[0],     &trange[1]);
                assert (count == 2);

                count = fscanf(fp, "%lf %lf", &trange_dot[0], &trange_dot[1]);
                assert (count == 2);

                if        (strcmp(edgeType, "LINE") == 0) {
                    nrdata = 6;
                } else if (strcmp(edgeType, "CIRCLE") == 0) {
                    nrdata = 10;
                } else if (strcmp(edgeType, "ELLIPSE") == 0) {
                    nrdata = 11;
                } else if (strcmp(edgeType, "PARABOLA") == 0) {
                    nrdata = 10;
                } else if (strcmp(edgeType, "HYPERBOLA") == 0) {
                    nrdata = 11;
                } else if (strcmp(edgeType, "BEZIER") == 0) {
                    for (i = 0; i < 3; i++) {
                        count = fscanf(fp, "%d", &idata[i]);
                        assert (count == 1);
                    }

                    nrdata =                       3 * idata[2];
                    if ((idata[0] & 2) != 0) nrdata += idata[2];
                } else if (strcmp(edgeType, "BSPLINE") == 0) {
                    for (i = 0; i < 4; i++) {
                        count = fscanf(fp, "%d", &idata[i]);
                        assert (count == 1);
                    }

                    nrdata = idata[3]            + 3 * idata[2];
                    if ((idata[0] & 2) != 0) nrdata += idata[2];
                } else {
                    printf("unknown edgeType \"%s\"\n", edgeType);
                    nrdata = 0;
                }

                MALLOC(rdata,     double, nrdata+1);
                MALLOC(rdata_dot, double, nrdata+1);

                for (i = 0; i < nrdata; i++) {
                    count = fscanf(fp, "%lf", &rdata[i]);
                    assert (count == 1);
                }
                dotmax = 0;
                for (i = 0; i < nrdata; i++) {
                    count = fscanf(fp, "%lf", &rdata_dot[i]);
                    assert (count == 1);

                    if (fabs(rdata_dot[i]) > dotmax) dotmax = fabs(rdata_dot[i]);
                }

                if (ipmtr > 0) {
                    SPRINT1x(2, "   dotmax=%10.5f", dotmax);

                    for (i = 0; i < nrdata; i++) {
                        rdata_dot[i] *= MODL->pmtr[ipmtr].dot[0];
                    }

                    status = EG_getInfo(ecurve, &oclass, &mtype, NULL, NULL, NULL);
                    CHECK_STATUS(EG_getInfo);

                    status = EG_setGeometry_dot(ecurve, oclass, mtype, idata, rdata, rdata_dot);
                    CHECK_STATUS(EG_setGeometry_dot);

                    status = EG_setRange_dot(MODL->body[ibody].edge[jedge].eedge, EDGE, trange, trange_dot);
                    CHECK_STATUS(EG_setRange_dot);
                }

                FREE(rdata    );
                FREE(rdata_dot);

                SPRINT0(2, " ");

                break;
            }
        }

        for (iface = 1; iface <= nface; iface++) {
            count = fscanf(fp, "%*s %d %s", &jface, faceType);
            assert (count == 2);

            snprintf(theName, 79, "Face %d", jface);

            for (jface = 1; jface <= nface; jface++) {
                status = EG_attributeRet(MODL->body[ibody].face[jface].eface, "Name",
                                         &attrType, &attrLen, &tempIlist, &tempRlist, &tempClist);
                CHECK_STATUS(EG_attributeRet);

                if (strcmp(theName, tempClist) != 0) continue;

                status = EG_getTopology(MODL->body[ibody].face[jface].eface, &esurface, &oclass, &mtype,
                                        data, &nchild, &echilds, &senses);
                CHECK_STATUS(EG_getTopology);

                SPRINT2x(2, "iface=%3d (in dots file) maps to jface=%3d (in ebody)", iface, jface);

                if        (strcmp(faceType, "PLANE") == 0) {
                    nrdata = 9;
                } else if (strcmp(faceType, "SPHERICAL") == 0) {
                    nrdata = 10;
                } else if (strcmp(faceType, "CONICAL") == 0) {
                    nrdata = 14;
                } else if (strcmp(faceType, "CYLINDRICAL") == 0) {
                    nrdata = 13;
                } else if (strcmp(faceType, "EXTRUSION") == 0) {
                    nrdata = 3;
                } else if (strcmp(faceType, "TOROIDAL") == 0) {
                    nrdata = 14;
                } else if (strcmp(faceType, "REVOLUTION") == 0) {
                    nrdata = 6;
                } else if (strcmp(faceType, "BEZIER") == 0) {
                    for (i = 0; i < 5; i++) {
                        count = fscanf(fp, "%d", &idata[i]);
                        assert (count == 1);
                    }

                    nrdata =                       3 * idata[2] * idata[4];
                    if ((idata[0] & 2) != 0) nrdata += idata[2] * idata[4];
                } else if (strcmp(faceType, "BSPLINE") == 0) {
                    for (i = 0; i < 7; i++) {
                        count = fscanf(fp, "%d", &idata[i]);
                        assert (count == 1);
                    }

                    nrdata = idata[3] + idata[6] + 3 * idata[2] * idata[5];
                    if ((idata[0] & 2) != 0) nrdata += idata[2] * idata[5];
                } else {
                    printf("unknown faceType \"%s\"\n", faceType);
                    nrdata = 0;
                }

                MALLOC(rdata,     double, nrdata+1);
                MALLOC(rdata_dot, double, nrdata+1);

                for (i = 0; i < nrdata; i++) {
                    count = fscanf(fp, "%lf", &rdata[i]);
                    assert (count == 1);
                }

                dotmax = 0;
                for (i = 0; i < nrdata; i++) {
                    count = fscanf(fp, "%lf", &rdata_dot[i]);
                    assert (count == 1);

                    if (fabs(rdata_dot[i]) > dotmax) dotmax = fabs(rdata_dot[i]);
                }

                if (ipmtr > 0) {
                    SPRINT1x(2, "   dotmax=%10.5f", dotmax);

                    for (i = 0; i < nrdata; i++) {
                        rdata_dot[i] *= MODL->pmtr[ipmtr].dot[0];
                    }

                    status = EG_getInfo(esurface, &oclass, &mtype, NULL, NULL, NULL);
                    CHECK_STATUS(EG_getInfo);

                    status = EG_setGeometry_dot(esurface, oclass, mtype, idata, rdata, rdata_dot);
                    CHECK_STATUS(EG_setGeometry_dot);
                }

                FREE(rdata    );
                FREE(rdata_dot);

                SPRINT0(2, " ");

                break;
            }
        }
    }

    /* close the dots file */
    fclose(fp);

    /* make sure the Body now has dots */
    status = EG_hasGeometry_dot(ebody);
    CHECK_STATUS(EG_hasGeometry_dot);

    SPRINT0(2, "\nExit  udpRunCAD/sensitivity\n");

cleanup:
    if (strlen(message) > 0) {
        printf("message from udpRunCAD/sensitivity: %s\n", message);
    }

    FREE(message);

    return status;
}
