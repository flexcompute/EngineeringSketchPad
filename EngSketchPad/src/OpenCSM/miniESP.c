/*
 ************************************************************************
 *                                                                      *
 * miniESP.c -- generate .stp and .dots file from a .csm/udc file       *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include "OpenCSM.h"

/***********************************************************************/
/*                                                                     */
/* global variables                                                    */
/*                                                                     */
/***********************************************************************/

/* global variable holding a MODL */
static char      casename[257];        /* name of case */

static double    dtime     = 0.00001;  /* nominal dtime for perturbation */
static int       outLevel  = 1;        /* default output level */

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char
#define STRNCPY(TO, FROM, LEN) strncpy(TO, FROM, LEN); TO[LEN-1] = '\0';
#define STRNCAT(TO, FROM, LEN) strncat(TO, FROM, LEN); TO[LEN-1] = '\0';
#define EPRINT  printf

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

/* declarations for high-level routines defined below */

/* declarations for support routines defined below */

/* external functions not declared in an .h file */


/***********************************************************************/
/*                                                                     */
/*   main - main program                                               */
/*                                                                     */
/***********************************************************************/

int
main(int       argc,                    /* (in)  number of arguments */
     char      *argv[])                 /* (in)  array of arguments */
{

    int       status, status2, i, nbody, ibody, iface, iedge, inode;
    int       imajor, iminor, builtTo, showUsage=0, count;
    int       oclass, mtype, nchild, *senses, *idata_base=NULL, *idata_ptrb=NULL;
    int       ipmtr, jpmtr, irow, icol, irowcol, nidata, nrdata, degen;
    double    data_base[18], data_ptrb[18], *rdata_base=NULL, *rdata_ptrb=NULL;
    double    value, value_base, dot_base, scaled_dtime;
    char      csmfile[257], udcfile[257], pmtrfile[257], errfile[257], stpfile[257], dotsfile[257];
    char      casename2[257];
    char      name[80], entname[80], message[101];
    CCHAR     *OCC_ver;
    void      *modl, *ptrb;
    ego       context, emodel, *echilds, esurface_base, esurface_ptrb;
    ego       ecurve_base, ecurve_ptrb, eref_base, eref_ptrb, eref;
    FILE      *fp, *fp2;

    modl_T    *MODL, *PTRB;

    ROUTINE(MAIN);

    /* --------------------------------------------------------------- */

    /* get the flags and casename from the command line */
    casename[0] = '\0';

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-dtime") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%lf", &dtime);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-help") == 0 ||
                   strcmp(argv[i], "-h"   ) == 0   ) {
            showUsage = 1;
            break;
        } else if (strcmp(argv[i], "-outLevel") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &outLevel);
                if (outLevel < 0) outLevel = 0;
                if (outLevel > 3) outLevel = 3;
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "--version") == 0 ||
                   strcmp(argv[i], "-version" ) == 0 ||
                   strcmp(argv[i], "-v"       ) == 0   ) {
            (void) ocsmVersion(&imajor, &iminor);
            SPRINT2(0, "OpenCSM version: %2d.%02d", imajor, iminor);
            EG_revision(&imajor, &iminor, &OCC_ver);
            SPRINT3(0, "EGADS   version: %2d.%02d (with %s)", imajor, iminor, OCC_ver);
            exit(EXIT_SUCCESS);
        } else if (strlen(casename) == 0) {
            strncpy(casename, argv[i], 256);
        } else {
            SPRINT0(0, "two casenames given");
            showUsage = 1;
            break;
        }
    }

    (void) ocsmVersion(&imajor, &iminor);

    /* command line information */
    if (showUsage) {
        SPRINT2(0, "miniESP version %2d.%02d\n", imajor, iminor);
        SPRINT0(0, "proper usage: 'miniESP casename [options...]");
        SPRINT0(0, "   where [options...] = -dtime dtime");
        SPRINT0(0, "                        -help  -or-  -h");
        SPRINT0(0, "                        -outLevel X");
        SPRINT0(0, "STOPPING...\a");
        return EXIT_FAILURE;
    }

    if (STRLEN(casename) == 0) {
        SPRINT0(0, "A casename must be given");
        return EXIT_FAILURE;
    }

    /* welcome banner */
    SPRINT0(1, "**********************************************************");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*                    Program miniESP                     *");
    SPRINT2(1, "*                     version %2d.%02d                      *", imajor, iminor);
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*           written by John Dannenhoffer, 2025           *");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "**********************************************************\n");

    SPRINT1(1, "    casename = %s", casename  );
    SPRINT1(1, "    dtime    = %f", dtime     );
    SPRINT1(1, "    outLevel = %d", outLevel  );
    SPRINT0(1, " ");

    /* set OCSMs output level */
    (void) ocsmSetOutLevel(outLevel);

    /* set up the name of the output files */
    snprintf(csmfile,  256, "miniESP.csm"         );
    snprintf(udcfile,  256, "%s.udc",     casename);
    snprintf(pmtrfile, 256, "%s.pmtr",    casename);
    snprintf(errfile,  256, "%s.err",     casename);
    snprintf(stpfile,  256, "%s.stp",     casename);
    snprintf(dotsfile, 256, "%s.dots",    casename);

    SPRINT0(1, "Files used:");
    SPRINT1(1, "    %-30s   (input  .udc file)",            udcfile );
    SPRINT1(1, "    %-30s   (input  file with DESPMTRs)",   pmtrfile);
    SPRINT1(1, "    %-30s   (output file with error msgs)", errfile );
    SPRINT1(1, "    %-30s   (output .stp file with names)", stpfile );
    SPRINT1(1, "    %-30s   (output file with dots)\n",     dotsfile);

    /* remove the output files */
    remove(errfile );
    remove(stpfile );
    remove(dotsfile);

    /* make a temporary casename in which all backslashes are converted
       to forward slashes (since OpenCSM will do the approporiate conversions) */
    for (i = 0; i < STRLEN(casename); i++) {
        if (casename[i] != '\\') {
            casename2[i] = casename[i];
        } else {
            casename2[i] = '/';
        }
        casename2[i+1] = '\0';
    }

    /* create the .csm file (which calls the .udc) */
    fp2 = fopen(csmfile, "w");
    if (fp2 != NULL) {
        fprintf(fp2, "UDPRIM /%s\n", casename2);
        fprintf(fp2, "DUMP   miniESP.stp\n");
        fprintf(fp2, "IMPORT miniESP.stp\n");
        fprintf(fp2, "END\n");
        fclose(fp2);
    }

    fp2 = fopen(csmfile, "r");
    if (fp2 != NULL) {
        SPRINT1(1, "csmfile=%s (which was written in miniESP) contains:", csmfile);
        while (1) {
            char buffer[100];
            (void) fgets(buffer, 99, fp2);
            if (feof(fp2) != 0) break;
            SPRINT1x(1, "     %s", buffer);
        }
        fclose(fp2);
    } else {
        printf("csmfile=%s does not exist\n", csmfile);
        exit(0);
    }

    /* read the .csm file and create the MODL.  after this call, the MODL
       will contain the Design Parameters and Branches, but no Bodys */
    status   = ocsmLoad(csmfile, &modl);
    MODL = (modl_T*)modl;

    SPRINT3(1, "--> ocsmLoad(%s) -> status=%d (%s)",
            csmfile, status, ocsmGetText(status));
    CHECK_STATUS(ocsmLoad);

    /* load the pmtr file (and create any needed DESPMTRs) */
    fp = fopen(pmtrfile, "r");
    if (fp != NULL) {
        while (1) {
            count = fscanf(fp, "%s %lf", name, &value);
            if (count != 2) break;

            ipmtr = 0;
            for (jpmtr = 1; jpmtr <= MODL->npmtr; jpmtr++) {
                if (strcmp(MODL->pmtr[jpmtr].name, name) == 0) {
                    if (MODL->pmtr[jpmtr].type != OCSM_DESPMTR) {
                        snprintf(message, 100, "Parameter \"%s\" is not a DESPMTR", name);
                        status = -900;
                        goto cleanup;
                    } else if (MODL->pmtr[jpmtr].nrow != 1 ||
                               MODL->pmtr[jpmtr].ncol != 1   ) {
                        snprintf(message, 100, "Parameter \"%s\" is not a scalar", name);
                        status = -901;
                        goto cleanup;
                    } else {
                        ipmtr = jpmtr;
                        break;
                    }
                }
            }

            if (ipmtr == 0) {
                status = ocsmNewPmtr(MODL, name, OCSM_DESPMTR, 1, 1);
                CHECK_STATUS(ocsmNewPmtr);

                ipmtr = MODL->npmtr;
            }

            status = ocsmSetValuD(MODL, ipmtr, 1, 1, value);
            CHECK_STATUS(ocsmSetValuD);

            SPRINT2(1, "--> setting \"%s\" = %10.5f", name, value);
        }

        fclose(fp);
    }

    /* check that Branches are properly ordered.  this also has the
       side-effect that it sets up the Branch's left and rite parents,
       children, and indentation */
    status   = ocsmCheck(modl);
    SPRINT2(0, "--> ocsmCheck() -> status=%d (%s)",
            status, ocsmGetText(status));
    CHECK_STATUS(ocsmCheck);

    /* build the Bodys from the MODL */
    nbody  = 0;
    status = ocsmBuild(modl, 0, &builtTo, &nbody, NULL);
    SPRINT4(1, "--> ocsmBuild -> status=%d (%s), builtTo=%d, nbody=%d",
            status, ocsmGetText(status), builtTo, nbody);

    if (status < EGADS_SUCCESS) {
        strncpy(message, MODL->sigMesg, 100);
        status = -906;
        goto cleanup;
    } else {
        nbody = MODL->nbody;
    }

    /* put a name on each Node, Edge, and Face */
    for (inode = 1; inode <= MODL->body[nbody].nnode; inode++) {
        snprintf(entname, 80, "Node %d", inode);
        SPRINT2(2, "putting Name=\"%s\" on Node %d", entname, inode);

        status = EG_attributeAdd(MODL->body[nbody].node[inode].enode, "Name",
                                 ATTRSTRING, 0, NULL, NULL, entname);
        CHECK_STATUS(EG_attributeAdd);
    }

    for (iedge = 1; iedge <= MODL->body[nbody].nedge; iedge++) {
        snprintf(entname, 80, "Edge %d", iedge);
        SPRINT2(2, "putting Name=\"%s\" on Edge %d", entname, iedge);

        status = EG_attributeAdd(MODL->body[nbody].edge[iedge].eedge, "Name",
                                 ATTRSTRING, 0, NULL, NULL, entname);
        CHECK_STATUS(EG_attributeAdd);
    }

    for (iface = 1; iface <= MODL->body[nbody].nface; iface++) {
        snprintf(entname, 80, "Face %d", iface);
        SPRINT2(2, "putting Name=\"%s\" on Face %d", entname, iface);

        status = EG_attributeAdd(MODL->body[nbody].face[iface].eface, "Name",
                                 ATTRSTRING, 0, NULL, NULL, entname);
        CHECK_STATUS(EG_attributeAdd);
    }

    /* dump a .stp file */
    status = EG_copyObject(MODL->body[nbody].ebody, NULL, &eref);
    CHECK_STATUS(EG_copyObject);

    status = EG_makeTopology(MODL->context, NULL, MODEL, 0,
                             NULL, 1, &eref, NULL, &emodel);
    CHECK_STATUS(EG_makeTopology);

    status = EG_saveModel(emodel, stpfile);
    CHECK_STATUS(EG_saveModel);

    status = EG_deleteObject(emodel);
    CHECK_STATUS(EG_deleteObject);

    SPRINT1(2, "\"%s\" has been written", stpfile);

    /* open the .dots file */
    fp = fopen(dotsfile, "w");
    if (fp == NULL) {
        snprintf(message, 100, "Could not open \"%s\" for writing", dotsfile);
        status = -902;
        goto cleanup;
    }

    fprintf(fp, "%5d %5d %5d\n", MODL->body[nbody].nnode,
                                      MODL->body[nbody].nedge,
                                      MODL->body[nbody].nface);

    /* loop through all the despmtrs to compute the sensitivities */
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].type != OCSM_DESPMTR) continue;

        irowcol = 0;
        for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {
                fprintf(fp, "%-30s %5d %5d %12.6f\n",
                        MODL->pmtr[ipmtr].name, irow, icol, MODL->pmtr[ipmtr].value[irowcol]);

                SPRINT4(1, "\nWorking on %s[%d,%d] (nominal=%12.6f)",
                        MODL->pmtr[ipmtr].name, irow, icol, MODL->pmtr[ipmtr].value[irowcol]);

                /* set up a perturbation and rebuild */
                status = ocsmCopy(MODL, &ptrb);
                SPRINT2(1, "--> ocsmCopy(MODL->PTRB) -> status=%d (%s)",
                        status, ocsmGetText(status));
                CHECK_STATUS(ocsmCopy);

                PTRB = (modl_T *)ptrb;
                PTRB->tessAtEnd = 0;
                PTRB->matchSeq  = MODL;

                status = ocsmGetValu(PTRB, ipmtr, irow, icol, &value_base, &dot_base);
                CHECK_STATUS(ocsmGetValu);

                scaled_dtime = dtime * MAX(fabs(value_base), +1.0);

                status = ocsmSetValuD(PTRB, ipmtr, irow, icol, value_base+scaled_dtime);
                CHECK_STATUS(ocsmSetValuD);

                SPRINT4(1, "--> ocsmSetValuD(PTRB, %s=%10.5f) -> status=%d (%s)",
                        MODL->pmtr[ipmtr].name, value_base+scaled_dtime, status , ocsmGetText(status));

                nbody = 0;
                (void) ocsmSetOutLevel(0);
                status = ocsmBuild(ptrb, 0, &builtTo, &nbody, NULL);
                (void) ocsmSetOutLevel(outLevel);
                SPRINT2(1, "--> ocsmBuild -> status=%d (%s)",
                        status, ocsmGetText(status));
                CHECK_STATUS(ocsmBuild);

                nbody = PTRB->nbody;
                SPRINT1(2, "\n     writing data to .dots file for nbody=%d", nbody);

                /* make sure that the topologies match */
                if (MODL->body[nbody].nnode != PTRB->body[nbody].nnode ||
                    MODL->body[nbody].nedge != PTRB->body[nbody].nedge ||
                    MODL->body[nbody].nface != PTRB->body[nbody].nface   ) {
                    snprintf(message, 100, "Topology mismatch");
                    SPRINT3(0, "        MODL: nnode=%5d, nedge=%5d, nface=%5d",
                           MODL->body[nbody].nnode, MODL->body[nbody].nedge, MODL->body[nbody].nface);
                    SPRINT3(0, "        PTRB: nnode=%5d, nedge=%5d, nface=%5d",
                           PTRB->body[nbody].nnode, PTRB->body[nbody].nedge, PTRB->body[nbody].nface);
                    status = -903;
                    goto cleanup;
                }

                /* loop through all the Nodes and find the dots via finite differences */
                for (inode = 1; inode <= MODL->body[nbody].nnode; inode++) {
                    status = EG_getTopology(MODL->body[nbody].node[inode].enode, &eref,
                                            &oclass, &mtype, data_base, &nchild, &echilds, &senses);
                    CHECK_STATUS(EG_getTopology);

                    status = EG_getTopology(PTRB->body[nbody].node[inode].enode, &eref,
                                            &oclass, &mtype, data_ptrb, &nchild, &echilds, &senses);
                    CHECK_STATUS(EG_getTopology);

                    SPRINT1(2,  "   Node %5d",   inode);
                    fprintf(fp, "   Node %5d\n", inode);
                    fprintf(fp, "       %18.11e %18.11e %18.11e\n",
                            data_base[0], data_base[1], data_base[2]);
                    fprintf(fp, "       %18.11e %18.11e %18.11e\n",
                            (data_ptrb[0]-data_base[0])/scaled_dtime,
                            (data_ptrb[1]-data_base[1])/scaled_dtime,
                            (data_ptrb[2]-data_base[2])/scaled_dtime);
                }

                /* loop through all the Edges and find the dots via finite differences */
                for (iedge = 1; iedge <= MODL->body[nbody].nedge; iedge++) {
                    status = EG_getInfo(MODL->body[nbody].edge[iedge].eedge, &oclass, &mtype, NULL, NULL, NULL);
                    CHECK_STATUS(EG_getInfo);
                    if (mtype == DEGENERATE) {
                        degen = 1;
                    } else {
                        degen = 0;
                    }

                    status = EG_getTopology(MODL->body[nbody].edge[iedge].eedge, &ecurve_base,
                                            &oclass, &mtype, data_base, &nchild, &echilds, &senses);
                    CHECK_STATUS(EG_getTopology);

                    status = EG_getTopology(PTRB->body[nbody].edge[iedge].eedge, &ecurve_ptrb,
                                            &oclass, &mtype, data_ptrb, &nchild, &echilds, &senses);
                    CHECK_STATUS(EG_getTopology);

                    if (degen == 0) {
                        status = EG_getGeometry(ecurve_base, &oclass, &mtype, &eref_base, &idata_base, &rdata_base);
                        CHECK_STATUS(EG_getGeometry);

                        status = EG_getGeometry(ecurve_ptrb, &oclass, &mtype, &eref_ptrb, &idata_ptrb, &rdata_ptrb);
                        CHECK_STATUS(EG_getGeometry);
                    }

                    /* TRIMMED curves get their info from their underlying geometry */
                    if (degen == 0 && mtype == TRIMMED) {
                        EG_free(rdata_base);   rdata_base = NULL;
                        EG_free(rdata_ptrb);   rdata_ptrb = NULL;

                        status = EG_getGeometry(eref_base, &oclass, &mtype, &eref, &idata_base, &rdata_base);
                        CHECK_STATUS(EG_getGeometry);

                        status = EG_getGeometry(eref_ptrb, &oclass, &mtype, &eref, &idata_ptrb, &rdata_ptrb);
                        CHECK_STATUS(EG_getGeometry);
                    }

                    /* determine the type of Curve associate with this Edge */
                    if (degen == 1) {
                        SPRINT1(2,  "   Edge %5d  DEGENERATE",   iedge);
                        fprintf(fp, "   Edge %5d  DEGENERATE\n", iedge);
                        nidata = 0;
                        nrdata = 0;
                    } else if (mtype == LINE) {
                        SPRINT1(2,  "   Edge %5d  LINE",   iedge);
                        fprintf(fp, "   Edge %5d  LINE\n", iedge);
                        nidata = 0;
                        nrdata = 6;
                    } else if (mtype == CIRCLE) {
                        SPRINT1(2,  "   Edge %5d  CIRCLE",   iedge);
                        fprintf(fp, "   Edge %5d  CIRCLE\n", iedge);
                        nidata = 0;
                        nrdata = 10;
                    } else if (mtype == ELLIPSE) {
                        SPRINT1(2,  "   Edge %5d  ELLIPSE",   iedge);
                        fprintf(fp, "   Edge %5d  ELLIPSE\n", iedge);
                        nidata = 0;
                        nrdata = 11;
                    } else if (mtype == PARABOLA) {
                        SPRINT1(2,  "   Edge %5d  PARABOLA",   iedge);
                        fprintf(fp, "   Edge %5d  PARABOLA\n", iedge);
                        nidata = 0;
                        nrdata = 10;
                    } else if (mtype == HYPERBOLA) {
                        SPRINT1(2,  "   Edge %5d  HYPERBOLA",   iedge);
                        fprintf(fp, "   Edge %5d  HYPERBOLA\n", iedge);
                        nidata = 0;
                        nrdata = 11;
                    } else if (mtype == BEZIER) {
                        SPLINT_CHECK_FOR_NULL(idata_base);

                        SPRINT1(2,  "   Edge %5d  BEZIER",   iedge);
                        fprintf(fp, "   Edge %5d  BEZIER\n", iedge);
                        nidata = 3;
                        nrdata =                            3 * idata_base[2];
                        if ((idata_base[0] & 2) != 0) nrdata += idata_base[2];
                    } else if (mtype == BSPLINE) {
                        SPLINT_CHECK_FOR_NULL(idata_base);

                        SPRINT1(2,  "   Edge %5d  BSPLINE",   iedge);
                        fprintf(fp, "   Edge %5d  BSPLINE\n", iedge);
                        nidata = 4;
                        nrdata = idata_base[3]            + 3 * idata_base[2];
                        if ((idata_base[0] & 2) != 0) nrdata += idata_base[2];
                    } else {
                        snprintf(message, 100, "Unknown Curve mtype=%d (at this time)", mtype);
                        status = -904;
                        goto cleanup;
                    }

                    fprintf(fp, "       %18.11e %18.11e\n",
                            data_base[0], data_base[1]);
                    fprintf(fp, "       %18.11e %18.11e\n",
                            (data_ptrb[0]-data_base[0])/scaled_dtime,
                            (data_ptrb[1]-data_base[1])/scaled_dtime);

                    if (degen == 1) continue;

                    /* write out the data (and dots) associated with this Face */
                    if (nidata > 0) {
                        SPLINT_CHECK_FOR_NULL(idata_base);

                        fprintf(fp, "      ");
                        for (i = 0; i < nidata; i++) {
                            fprintf(fp, " %18d", idata_base[i]);
                        }
                        fprintf(fp, "\n");
                    }
                    if (nrdata > 0) {
                        SPLINT_CHECK_FOR_NULL(rdata_base);
                        SPLINT_CHECK_FOR_NULL(rdata_ptrb);

                        fprintf(fp, "      ");
                        for (i = 0; i < nrdata; i++) {
                            fprintf(fp, " %18.11e", rdata_base[i]);
                        }
                        fprintf(fp, "\n      ");
                        for (i = 0; i < nrdata; i++) {
                            fprintf(fp, " %18.11e", (rdata_ptrb[i] - rdata_base[i])/scaled_dtime);
                        }
                        fprintf(fp, "\n");
                    }

                    if (idata_base != NULL) {EG_free(idata_base);   idata_base = NULL;}
                    if (idata_ptrb != NULL) {EG_free(idata_ptrb);   idata_ptrb = NULL;}
                    if (rdata_base != NULL) {EG_free(rdata_base);   rdata_base = NULL;}
                    if (rdata_ptrb != NULL) {EG_free(rdata_ptrb);   rdata_ptrb = NULL;}
                }

                /* loop through all the Faces and find the dots via finite differences */
                for (iface = 1; iface <= MODL->body[nbody].nface; iface++) {
                    status = EG_getTopology(MODL->body[nbody].face[iface].eface, &esurface_base,
                                            &oclass, &mtype, data_base, &nchild, &echilds, &senses);
                    CHECK_STATUS(EG_getTopology);

                    status = EG_getTopology(PTRB->body[nbody].face[iface].eface, &esurface_ptrb,
                                            &oclass, &mtype, data_ptrb, &nchild, &echilds, &senses);
                    CHECK_STATUS(EG_getTopology);

                    status = EG_getGeometry(esurface_base, &oclass, &mtype, &eref_base, &idata_base, &rdata_base);
                    CHECK_STATUS(EG_getGeometry);

                    status = EG_getGeometry(esurface_ptrb, &oclass, &mtype, &eref_ptrb, &idata_ptrb, &rdata_ptrb);
                    CHECK_STATUS(EG_getGeometry);

                    /* TRIMMED surfaces get their info from their underlying geometry */
                    if (mtype == TRIMMED) {
                        EG_free(rdata_base);   rdata_base = NULL;
                        EG_free(rdata_ptrb);   rdata_ptrb = NULL;

                        status = EG_getGeometry(eref_base, &oclass, &mtype, &eref, &idata_base, &rdata_base);
                        CHECK_STATUS(EG_getGeometry);

                        status = EG_getGeometry(eref_ptrb, &oclass, &mtype, &eref, &idata_ptrb, &rdata_ptrb);
                        CHECK_STATUS(EG_getGeometry);
                    }

                    /* determine the type of Surface associate with this Face */
                    if        (mtype == PLANE) {
                        SPRINT1(2,  "   Face %5d  PLANE",   iface);
                        fprintf(fp, "   Face %5d  PLANE\n", iface);
                        nidata = 0;
                        nrdata = 9;
                    } else if (mtype == SPHERICAL) {
                        SPRINT1(2,  "   Face %5d  SPHERICAL",   iface);
                        fprintf(fp, "   Face %5d  SPHERICAL\n", iface);
                        nidata = 0;
                        nrdata = 10;
                    } else if (mtype == CONICAL) {
                        SPRINT1(2,  "   Face %5d  CONICAL",   iface);
                        fprintf(fp, "   Face %5d  CONICAL\n", iface);
                        nidata = 0;
                        nrdata = 14;
                    } else if (mtype == CYLINDRICAL) {
                        SPRINT1(2,  "   Face %5d  CYLINDRICAL",   iface);
                        fprintf(fp, "   Face %5d  CYLINDRICAL\n", iface);
                        nidata = 0;
                        nrdata = 13;
                    } else if (mtype == EXTRUSION) {
                        SPRINT1(2,  "   Face %5d  EXTRUSION",   iface);
                        fprintf(fp, "   Face %5d  EXTRUSION\n", iface);
                        nidata = 0;
                        nrdata = 3;
                    } else if (mtype == TOROIDAL) {
                        SPRINT1(2,  "   Face %5d  TOROIDAL",   iface);
                        fprintf(fp, "   Face %5d  TOROIDAL\n", iface);
                        nidata = 0;
                        nrdata = 14;
                    } else if (mtype == REVOLUTION) {
                        SPRINT1(2,  "   Face %5d  REVOLUTION",   iface);
                        fprintf(fp, "   Face %5d  REVOLUTION\n", iface);
                        nidata = 0;
                        nrdata = 6;
                    } else if (mtype == BEZIER) {
                        SPLINT_CHECK_FOR_NULL(idata_base);

                        SPRINT1(2,  "   Face %5d  BEZIER",   iface);
                        fprintf(fp, "   Face %5d  BEZIER\n", iface);
                        nidata = 5;
                        nrdata =                            3 * idata_base[2] * idata_base[4];
                        if ((idata_base[0] & 2) != 0) nrdata += idata_base[2] * idata_base[4];
                    } else if (mtype == BSPLINE) {
                        SPLINT_CHECK_FOR_NULL(idata_base);

                        SPRINT1(2,  "   Face %5d  BSPLINE",   iface);
                        fprintf(fp, "   Face %5d  BSPLINE\n", iface);
                        nidata = 7;
                        nrdata = idata_base[3] + idata_base[6] + 3 * idata_base[2] * idata_base[5];
                        if ((idata_base[0] & 2) != 0) nrdata +=      idata_base[2] * idata_base[5];
                    } else {
                        snprintf(message, 100, "Unknown Surface mtype=%d (at this time)", mtype);
                        status = -905;
                        goto cleanup;
                    }

                    /* write out the data (and dots) associated with this Face */
                    if (nidata > 0) {
                        SPLINT_CHECK_FOR_NULL(idata_base);

                        fprintf(fp, "      ");
                        for (i = 0; i < nidata; i++) {
                            fprintf(fp, " %18d", idata_base[i]);
                        }
                        fprintf(fp, "\n");
                    }
                    if (nrdata > 0) {
                        SPLINT_CHECK_FOR_NULL(rdata_base);
                        SPLINT_CHECK_FOR_NULL(rdata_ptrb);

                        fprintf(fp, "      ");
                        for (i = 0; i < nrdata; i++) {
                            fprintf(fp, " %18.11e", rdata_base[i]);
                        }
                        fprintf(fp, "\n      ");
                        for (i = 0; i < nrdata; i++) {
                            fprintf(fp, " %18.11e", (rdata_ptrb[i] - rdata_base[i])/scaled_dtime);
                        }
                        fprintf(fp, "\n");
                    }

                    if (idata_base != NULL) {EG_free(idata_base);   idata_base = NULL;}
                    if (idata_ptrb != NULL) {EG_free(idata_ptrb);   idata_ptrb = NULL;}
                    if (rdata_base != NULL) {EG_free(rdata_base);   rdata_base = NULL;}
                    if (rdata_ptrb != NULL) {EG_free(rdata_ptrb);   rdata_ptrb = NULL;}
                }

                /* get ready for next perturbation */
                status = ocsmFree(ptrb);
                CHECK_STATUS(ocsmFree);

                irowcol++;
            }
        }
    }

    fclose(fp);

    remove("miniESP.csm");
    remove("miniESP.stp");

    SPRINT0(1, "\n    ******************************");
    SPRINT0(1,   "==> miniESP completed successfully");
    SPRINT0(1,   "    ******************************\n");

    status = EXIT_SUCCESS;

    /* cleanup and exit */
cleanup:
    if (status < 0) {
        fp = fopen(errfile, "w");
        if (fp != NULL) {
            fprintf(fp, "ERROR %3d:: %s\n", status, message);
            fclose(fp);
        }
    }

    context = MODL->context;

    /* remove all Bodys and etess objects */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].etess != NULL) {
            status2 = EG_deleteObject(MODL->body[ibody].etess);
            SPRINT3(2, "--> EG_deleteObject(etess[%d]) -> status=%d (%s)",
                    ibody, status2, ocsmGetText(status2));

            MODL->body[ibody].etess = NULL;
        }

        if (MODL->body[ibody].ebody != NULL) {
            status2 = EG_deleteObject(MODL->body[ibody].ebody);
            SPRINT3(2, "--> EG_deleteObject(ebody[%d]) => status=%d (%s)",
                    ibody, status2, ocsmGetText(status2));

            MODL->body[ibody].ebody = NULL;
        }
    }

    /* free up the modl */
    status2 = ocsmFree(modl);
    SPRINT2(1, "--> ocsmFree() -> status=%d (%s)",
            status2, ocsmGetText(status2));

    /* clean up the udp storage */
    status2 = ocsmFree(NULL);
    SPRINT2(1, "--> ocsmFree(NULL) -> status=%d (%s)",
            status2, ocsmGetText(status2));

    /* remove the context */
    if (context != NULL) {
        status2 = EG_setOutLevel(context, 0);
        if (status2 < 0) {
            SPRINT2(0, "EG_setOutLevel -> status=%d (%s)",
                    status2, ocsmGetText(status2));
        }

        status2 = EG_close(context);
        SPRINT2(1, "--> EG_close() -> status=%d (%s)",
                status2, ocsmGetText(status2));
    }

    if (status < 0) status = EXIT_FAILURE;

    return status;
}
