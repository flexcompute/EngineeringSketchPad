/*
 ************************************************************************
 *                                                                      *
 * sensCSM.c -- test geometric and tessellation sensitivities           *
 *                                                                      *
 *              Written by John Dannenhoffer @ Syracuse University      *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "egads.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define STRNCPY(A, B, LEN) strncpy(A, B, LEN); A[LEN-1] = '\0';

#include "egads.h"
#include "common.h"
#include "OpenCSM.h"
#include "udp.h"

#if !defined(__APPLE__) && !defined(WIN32)
// floating point exceptions
   #define __USE_GNU
   #include <fenv.h>
#endif

#ifdef WIN32
    #define  SLASH '\\'
#else
    #define  SLASH '/'
#endif

#define  ACCEPTABLE_ERROR   1.0e-6
#define  ERROR_RATIO        2.0
#define  ERROR_TOLER        1e-4
#define  ERROR_REPORT       1e-4

/***********************************************************************/
/*                                                                     */
/* macros (including those that go along with common.h)                */
/*                                                                     */
/***********************************************************************/

                                               /* used by RALLOC macro */
//static void *realloc_temp=NULL;

#define  RED(COLOR)      (float)(COLOR / 0x10000        ) / (float)(255)
#define  GREEN(COLOR)    (float)(COLOR / 0x00100 % 0x100) / (float)(255)
#define  BLUE(COLOR)     (float)(COLOR           % 0x100) / (float)(255)

/***********************************************************************/
/*                                                                     */
/* global variables                                                    */
/*                                                                     */
/***********************************************************************/

/* global variable holding a MODL */
static char      casename[255];        /* name of case */
static char      pmtrname[255];        /* name of single design Parameter */
static void      *modl;                /* pointer to MODL */

static int       addVerify= 0;         /* =1 to write .gsen or .tsen file */
static int       geom     = 0;         /* =1 for geometric sensitivites */
static double    dtime    = 0;         /* dtime for perturbation */
static int       outLevel = 1;         /* default output level */
static int       tess     = 0;         /* =1 for tessellation sensitivities */
static int       showAll  = 0;         /* =1 to show all velocities */

static clock_t   geom_time;            /* total CPU time associated with ocsmGetVel */
static clock_t   tess_time;            /* total CPU time associated with ocsmGetTessVel */

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

/* declarations for high-level routines defined below */
static int checkGeomSens(int ipmtr, int irow, int icol, int *ntotal, int *nsuppress, double *errmax);
static int checkTessSens(int ipmtr, int irow, int icol, int *ntotal, int *nsuppress, double *errmax);


/***********************************************************************/
/*                                                                     */
/*   main - main program                                               */
/*                                                                     */
/***********************************************************************/

int
main(int       argc,                    /* (in)  number of arguments */
     char      *argv[])                 /* (in)  array of arguments */
{

    int       status, fileStatus, status2, i, nbody;
    int       imajor, iminor, builtTo, showUsage=0;
    int       ipmtr, irow, icol, iirow, iicol, ntotal, nsuppress=0, nerror=0;
    int       onlyrow=-1, onlycol=-1;
    double    errmax, error, errmaxGeom=0, errmaxTess=0, dtime_in=0;
    char      basename[MAX_FILENAME_LEN], dirname[MAX_FILENAME_LEN];
    char      filename[MAX_FILENAME_LEN], pname[  MAX_FILENAME_LEN];
    char      *beg, *mid, *end;
    CCHAR     *OCC_ver;
    ego       context;
    FILE      *fp_data=NULL;

    modl_T    *MODL;

    ROUTINE(MAIN);

    /* --------------------------------------------------------------- */

    fileStatus = EXIT_SUCCESS;

    /* get the flags and casename from the command line */
    casename[0] = '\0';
    pmtrname[0] = '\0';

    for (i = 1; i < argc; i++) {
        if        (strcmp(argv[i], "--") == 0) {
            /* ignore (needed for gdb) */
        } else if (strcmp(argv[i], "-addVerify") == 0) {
            addVerify = 1;
        } else if (strcmp(argv[i], "-geom") == 0) {
            geom = 1;
        } else if (strcmp(argv[i], "-despmtr") == 0) {
            if (i < argc-1) {
                strcpy(pmtrname, argv[++i]);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-dtime") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%lf", &dtime_in);
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
        } else if (strcmp(argv[i], "-showAll") == 0) {
            showAll = 1;
        } else if (strcmp(argv[i], "-tess") == 0) {
            tess = 1;
        } else if (strcmp(argv[i], "--version") == 0 ||
                   strcmp(argv[i], "-version" ) == 0 ||
                   strcmp(argv[i], "-v"       ) == 0   ) {
            (void) ocsmVersion(&imajor, &iminor);
            SPRINT2(0, "OpenCSM version: %2d.%02d", imajor, iminor);
            EG_revision(&imajor, &iminor, &OCC_ver);
            SPRINT3(0, "EGADS   version: %2d.%02d (with %s)", imajor, iminor, OCC_ver);
            exit(EXIT_SUCCESS);
        } else if (strlen(casename) == 0) {
            strcpy(casename, argv[i]);
        } else {
            SPRINT0(0, "two casenames given");
            showUsage = 1;
            break;
        }
    }

    (void) ocsmVersion(&imajor, &iminor);

    if (showUsage) {
        SPRINT2(0, "sensCSM version %2d.%02d\n", imajor, iminor);
        SPRINT0(0, "proper usage: 'sensCSM [casename[.csm]] [options...]");
        SPRINT0(0, "   where [options...] = -addVerify");
        SPRINT0(0, "                        -despmtr pmtrname");
        SPRINT0(0, "                        -dtime dtime");
        SPRINT0(0, "                        -geom");
        SPRINT0(0, "                        -help  -or-  -h");
        SPRINT0(0, "                        -outLevel X");
        SPRINT0(0, "                        -showAll");
        SPRINT0(0, "                        -tess");
        SPRINT0(0, "STOPPING...\a");
        return EXIT_FAILURE;
    }

    /* if neither geom or tess are set, raise an error */
    if (geom ==0 && tess == 0) {
        SPRINT0(0, "ERROR:: either -geom or -tess must be set");
        SPRINT0(0, "STOPPING...\a");
        return EXIT_FAILURE;
    }

    /* set the dtime as either the default or the user's specification */
    if (dtime_in > 0) {
        dtime = dtime_in;
    } else if (geom == 1) {
        dtime = 1.0e-6;
    } else if (tess == 1) {
        dtime = 1.0e-3;
    }

    /* break pmtrname into name, iirow, and iicol */
    /* if pmtrname contains a subscript, estract it now */
    SPLINT_CHECK_FOR_NULL(pmtrname);

    if (strlen(pmtrname) > 0 && strchr(pmtrname, '[') != NULL) {
        beg = strchr(pmtrname, '[');
        mid = strchr(pmtrname, ',');
        end = strchr(pmtrname, ']');

        if (beg == NULL || mid == NULL || end == NULL) {
            SPRINT0(0, "if -despmtr is given, pmtrname must be of form \"name\" or \"name[irow,icol]\"\a");
            SPRINT0(0, "STOPPING...\a");
            return EXIT_FAILURE;
        }

        *end = '\0';
        onlycol = atoi(mid+1);

        *mid = '\0';
        onlyrow = atoi(beg+1);

        *beg = '\0';
    }

    /* welcome banner */
    SPRINT0(1, "**********************************************************");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*                    Program sensCSM                     *");
    SPRINT2(1, "*                     version %2d.%02d                      *", imajor, iminor);
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*        written by John Dannenhoffer, 2010/2022         *");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "**********************************************************\n");

    SPRINT1(1, "    casename   = %s", casename  );
    SPRINT1(1, "    addVerify  = %d", addVerify );
    SPRINT1(1, "    despmtr    = %s", pmtrname  );
    SPRINT1(1, "    geom       = %d", geom      );
    SPRINT1(1, "    onlyrow    = %d", onlyrow   );
    SPRINT1(1, "    onlycol    = %d", onlycol   );
    SPRINT1(1, "    dtime      = %f", dtime     );
    SPRINT1(1, "    outLevel   = %d", outLevel  );
    SPRINT1(1, "    showAll    = %d", showAll   );
    SPRINT1(1, "    tess       = %d", tess      );
    SPRINT0(1, " ");

    /* set OCSMs output level */
    (void) ocsmSetOutLevel(outLevel);

    /* strip off .csm (which is assumed to be at the end) if present */
    if (strlen(casename) > 0) {
        strcpy(filename, casename);
        if (strstr(casename, ".csm") == NULL) {
            strcat(filename, ".csm");
        }
    } else {
        filename[0] = '\0';
    }

    /* get the OpenCASCADE version */
    EG_revision(&imajor, &iminor, &OCC_ver);

    /* get basename and dirname */
    SPLINT_CHECK_FOR_NULL(casename);

    i = strlen(casename) - 1;
    while (i >= 0) {
        if (casename[i] == '/' || casename[i] == '\\') {
            i++;
            break;
        }
        i--;
    }
    if (i == -1) {
        strcpy(dirname, ".");
        strcpy(basename, casename);
    } else {
        strcpy(basename, &(filename[i]));
        strcpy(dirname, casename);
        dirname[i-1] = '\0';
    }

    /* remove .csm or .cpc extension */
    i = strlen(basename);
    basename[i-4] = '\0';

    /* read the .csm file and create the MODL */
    status   = ocsmLoad(filename, &modl);
    MODL = (modl_T*)modl;

    SPRINT3(1, "--> ocsmLoad(%s) -> status=%d (%s)",
            filename, status, ocsmGetText(status));
    CHECK_STATUS(ocsmLoad);

    /* check that Branches are properly ordered */
    status   = ocsmCheck(modl);
    SPRINT2(0, "--> ocsmCheck() -> status=%d (%s)",
            status, ocsmGetText(status));
    CHECK_STATUS(ocsmCheck);

    /* print out the global Attributes, Parameters, and Branches */
    SPRINT0(1, "External Parameter(s):");
    if (outLevel > 0) {
        status = ocsmPrintPmtrs(modl, "");
        CHECK_STATUS(ocsmPrintPmtrs);
    }

    SPRINT0(1, "Branch(es):");
    if (outLevel > 0) {
        status = ocsmPrintBrchs(modl, "");
        CHECK_STATUS(ocsmPrintBrchs);
    }

    SPRINT0(1, "Global Attribute(s):");
    if (outLevel > 0) {
        status = ocsmPrintAttrs(modl, "");
        CHECK_STATUS(ocsmPrintAttrs);
    }

    /* build the Bodys from the MODL */
    nbody    = 0;
    status   = ocsmBuild(modl, 0, &builtTo, &nbody, NULL);
    SPRINT4(1, "--> ocsmBuild -> status=%d (%s), builtTo=%d, nbody=%d",
            status, ocsmGetText(status), builtTo, nbody);
    CHECK_STATUS(ocsmBuild);

    /* print out the Bodys */
    SPRINT0(1, "Body(s):");
    if (outLevel > 0) {
        status = ocsmPrintBodys(modl, "");
        CHECK_STATUS(ocsmPrintBodys);
    }

    /* create the full filename */
    if (geom) {
        snprintf(filename, MAX_FILENAME_LEN-1, "%s%cverify_%s%c%s.gsen", dirname, SLASH, &(OCC_ver[strlen(OCC_ver)-5]), SLASH, basename);
        if (addVerify) {
            fp_data = fopen(filename, "w");
        } else {
            fp_data = fopen(filename, "r");
        }

        if (fp_data == NULL) {
            SPRINT0(0, "ERROR:: geom error with .gsen file");
            fileStatus = EXIT_FAILURE;
        }
    } else {
        snprintf(filename, MAX_FILENAME_LEN-1, "%s%cverify_%s%c%s.tsen", dirname, SLASH, &(OCC_ver[strlen(OCC_ver)-5]), SLASH, basename);
        if (addVerify) {
            fp_data = fopen(filename, "w");
        } else {
            fp_data = fopen(filename, "r");
        }

        if (fp_data == NULL) {
            SPRINT0(0, "ERROR:: tess error with .tsen file");
            fileStatus = EXIT_FAILURE;
        }
    }

    /* check geometric and tessellation sensitivities */
    ntotal     = 0;              // total number of errors exceeding tolerance
    errmaxGeom = 0;              // maximum error (geometric)
    errmaxTess = 0;              // maximum error (tessellation)
    nerror     = 0;              // number of sensitivites whose error has increased

    /* loop through all design Parameters */
    for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
        if (MODL->pmtr[ipmtr].type != OCSM_DESPMTR) continue;

        if (strlen(pmtrname) > 0 &&
            strcmp(pmtrname, MODL->pmtr[ipmtr].name) != 0) continue;

        for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {

                if (onlyrow > 0 && onlycol > 0) {
                    irow = onlyrow;
                    icol = onlycol;
                }

                if (geom) {
                    errmax = EPS20;
                    status = checkGeomSens(ipmtr, irow, icol, &ntotal, &nsuppress, &errmax);
                    if (status != SUCCESS) {
                        SPRINT1(0, "ERROR:: geom error detected in checkGeomSens (status=%d)", status);
                        fileStatus = EXIT_FAILURE;
                    }
                    CHECK_STATUS(check_status);

                    if (addVerify) {
                        if (fp_data != NULL) {
                            fprintf(fp_data, "%-32s %5d %5d %12.5e\n", MODL->pmtr[ipmtr].name, irow, icol, errmax);
                        }

                        SPRINT4(1, "INFO:: geom error for %32s[%d,%d] is%12.5e being written to file",
                                MODL->pmtr[ipmtr].name, irow, icol, errmax);
                    } else if (fp_data != NULL) {
                        fscanf(fp_data, "%s %d %d %lf", pname, &iirow, &iicol, &error);

                        if (strcmp(MODL->pmtr[ipmtr].name, pname) != 0 || irow != iirow || icol != iicol) {
                            SPRINT0(0, "ERROR:: .gsen file does not match case");
                            fileStatus = EXIT_FAILURE;
                        } else if (errmax < ACCEPTABLE_ERROR) {
                        } else if (errmax > error*ERROR_RATIO) {
                            SPRINT5(0, "ERROR:: geom error for %32s[%d,%d] increased from %12.5e to %12.5e",
                                    MODL->pmtr[ipmtr].name, irow, icol, error, errmax);
                            nerror++;
                        } else if (errmax < error/ERROR_RATIO) {
                            SPRINT5(1, "INFO:: geom error for %32s[%d,%d] decreased from %12.5e to %12.5e",
                            MODL->pmtr[ipmtr].name, irow, icol, error, errmax);
                        }
                    }

                    if (errmax > errmaxGeom) errmaxGeom = errmax;
                }

                if (tess) {
                    errmax = EPS20;
                    status = checkTessSens(ipmtr, irow, icol, &ntotal, &nsuppress, &errmax);
                    if (status != SUCCESS) {
                        SPRINT1(0, "ERROR:: tess error detected in checkTessSens (status=%d)", status);
                        fileStatus = EXIT_FAILURE;
                    }
                    CHECK_STATUS(checkTessSens);

                    if (addVerify) {
                        if (fp_data != NULL) {
                            fprintf(fp_data, "%-32s %5d %5d %12.5e\n", MODL->pmtr[ipmtr].name, irow, icol, errmax);
                        }

                        SPRINT4(1, "INFO:: tess error for %32s[%d,%d] is%12.5e being written to file",
                                MODL->pmtr[ipmtr].name, irow, icol, errmax);
                    } else if (fp_data != NULL) {
                        fscanf(fp_data, "%s %d %d %lf", pname, &iirow, &iicol, &error);

                        if (strcmp(MODL->pmtr[ipmtr].name, pname) != 0 || irow != iirow || icol != iicol) {
                            SPRINT0(0, "ERROR:: .tsen file does not match case");
                            fileStatus = EXIT_FAILURE;
                        } else if (errmax < ACCEPTABLE_ERROR) {
                        } else if (errmax > error*ERROR_RATIO) {
                            SPRINT5(0, "ERROR:: tess error for %32s[%d,%d] increased from %12.5e to %12.5e",
                                    MODL->pmtr[ipmtr].name, irow, icol, error, errmax);
                            nerror++;
                        } else if (errmax < error/ERROR_RATIO) {
                            SPRINT5(0, "INFO:: tess error for %32s[%d,%d] decreased from %12.5e to %12.5e",
                                    MODL->pmtr[ipmtr].name, irow, icol, error, errmax);
                        }
                    }

                    if (errmax > errmaxTess) errmaxTess = errmax;
                }

                if (onlyrow > 0 && onlycol > 0) break;
            }

            if (onlyrow > 0 && onlycol > 0) break;
        }
        SPRINT0(0, " ");
    }

    if (fp_data != NULL) fclose(fp_data);

    SPRINT0(0, "==> sensCSM completed successfully");
    status = EXIT_SUCCESS;

    /* cleanup and exit */
cleanup:
    context = MODL->context;

    /* free up the modl (which removes all Bodys and etess objects) */
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

    /* report total CPU times */
    if (geom) {
        SPRINT1(0, "\nTotal CPU time in ocsmGetVel     -> %10.3f sec",
                (double)(geom_time)/(double)(CLOCKS_PER_SEC));
    }
    if (tess) {
        SPRINT1(0, "\nTotal CPU time in ocsmGetTessVel -> %10.3f sec",
                (double)(tess_time)/(double)(CLOCKS_PER_SEC));
    }

    /* report final statistics */
    if (status == SUCCESS) {
        if (geom) {
            SPRINT3(0, "\nSensitivity checks complete with %8d total errors (max geom err=%12.4e) with %d suppressions",
                    ntotal, errmaxGeom+1.0e-20, nsuppress);
        }
        if (tess) {
            SPRINT3(0, "\nSensitivity checks complete with %8d total errors (max tess err=%12.4e) with %d suppressions",
                    ntotal, errmaxTess+1.0e-20, nsuppress);
        }
    } else if (status == EXIT_FAILURE) {
        SPRINT0(0, "\nSensitivity checks not complete because \"EXIT_FAILURE\" was detected");
    } else {
        SPRINT1(0, "\nSensitivity checks not complete because error \"%s\" was detected",
                ocsmGetText(status));
    }

    /* these statements are in case we used an error return to go to cleanup */
    if (fileStatus == EXIT_FAILURE) status = EXIT_FAILURE;
    if (status     < 0            ) status = EXIT_FAILURE;
    if (nerror     > 0            ) status = EXIT_FAILURE;

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   checkGeomSens - check geometric sensitivities                     */
/*                                                                     */
/*                     this is done be comparing analytic results      */
/*                     with those computed via finite differenced      */
/*                                                                     */
/***********************************************************************/

static int
checkGeomSens(int    ipmtr,             /* (in)  Parameter index (bias-1) */
              int    irow,              /* (in)  row    index (bias-1) */
              int    icol,              /* (in)  column index (bias-1) */
              int    *ntotal,           /* (out) total number of points beyond toler */
              int    *nsuppress,        /* (out) number of suppressions */
              double *errmax)           /* (out) maximum error */
{
    int       status = SUCCESS;

    int       i, ibody, nerror, inode, iedge, iface, ipnt, ixyz, npnt_tess, ntri_tess, nrms;
    int       ntemp, builtTo, atype, alen;
    CINT      *ptype, *pindx, *tris, *tric, *tempIlist;
    double    errrms, face_errmax, edge_errmax, node_errmax, errX, errY, errZ;
    double    **face_anal=NULL, **face_fdif=NULL, face_err;
    double    **edge_anal=NULL, **edge_fdif=NULL, edge_err;
    double    **node_anal=NULL, **node_fdif=NULL, node_err;
    CDOUBLE   *xyz, *uv, *tempRlist;
    CCHAR     *tempClist;
    clock_t   old_time, new_time;

    modl_T    *MODL = (modl_T *)modl;

    ROUTINE(checkGeomSens);

    /* --------------------------------------------------------------- */

    SPRINT0(0, "\n*********************************************************");
    if (MODL->pmtr[ipmtr].nrow == 1 &&
        MODL->pmtr[ipmtr].ncol == 1   ) {
        SPRINT1(0, "Starting geometric sensitivity wrt \"%s\"",
                MODL->pmtr[ipmtr].name);
    } else {
        SPRINT3(0, "Starting geometric sensitivity wrt \"%s[%d,%d]\"",
                MODL->pmtr[ipmtr].name, irow, icol);
    }
    SPRINT0(0, "*********************************************************\n");

    SPRINT0(0, "Propagating velocities throughout feature tree");

    ibody = MODL->nbody + 1;            /* protect for warning in cleanup: */

    /* set the velocity */
    status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
    CHECK_STATUS(ocsmSetVelD);
    status = ocsmSetVelD(MODL, ipmtr, irow, icol, 1.0);
    CHECK_STATUS(ocsmSetVelD);

    /* build needed to propagate velocities to the Branch arguments */
    ntemp    = 0;
    status   = ocsmBuild(MODL, 0, &builtTo, &ntemp, NULL);
    CHECK_STATUS(ocsmBuild);

    /* perform geometric sensitivity checks for each Body on the stack */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        /* analytic geometric sensitivities (if possible) */
        SPRINT1(0, "Computing analytic sensitivities (if possible) for ibody=%d",
                ibody);

        status = ocsmSetDtime(MODL, 0);
        CHECK_STATUS(ocsmSetDtime);

        /* save analytic geometric sensitivity of each Face */
        MALLOC(face_anal, double*, MODL->body[ibody].nface+1);
        for (iface = 0; iface <= MODL->body[ibody].nface; iface++) {
            face_anal[iface] = NULL;
        }

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                    &ntri_tess, &tris, &tric);
            CHECK_STATUS(EG_getTessFace);
            if (npnt_tess <= 0) {
                SPRINT3(0, "ERROR:: EG_getTessFace -> status=%d (%s), npnt_tess=%d",
                        status, ocsmGetText(status), npnt_tess);
                status = EXIT_FAILURE;
                goto cleanup;
            }

            MALLOC(face_anal[iface], double, 3*npnt_tess);

            old_time = clock();
            status = ocsmGetVel(MODL, ibody, OCSM_FACE, iface, npnt_tess, NULL, face_anal[iface]);
            CHECK_STATUS(ocsmGetVel);
            new_time = clock();
            geom_time += (new_time - old_time);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (face_anal[iface][ixyz] != face_anal[iface][ixyz]) {
                    face_anal[iface][ixyz] = HUGEQ;
                }
            }
        }

        /* save analytic geometric sensitivity of each Edge */
        MALLOC(edge_anal, double*, MODL->body[ibody].nedge+1);
        for (iedge = 0; iedge <= MODL->body[ibody].nedge; iedge++) {
            edge_anal[iedge] = NULL;
        }

        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
            status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                    &npnt_tess, &xyz, &uv);
            CHECK_STATUS(EG_getTessEdge);
            if (npnt_tess <= 0) {
                SPRINT3(0, "ERROR:: EG_getTessEdge -> status=%d (%s), npnt_tess=%d",
                        status, ocsmGetText(status), npnt_tess);
                status = EXIT_FAILURE;
                goto cleanup;
            }

            MALLOC(edge_anal[iedge], double, 3*npnt_tess);

            old_time = clock();
            status = ocsmGetVel(MODL, ibody, OCSM_EDGE, iedge, npnt_tess, NULL, edge_anal[iedge]);
            CHECK_STATUS(ocsmGetVel);
            new_time = clock();
            geom_time += (new_time - old_time);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (edge_anal[iedge][ixyz] != edge_anal[iedge][ixyz]) {
                    edge_anal[iedge][ixyz] = HUGEQ;
                }
            }
        }

        /* save analytic geometric sensitivity of each Node */
        MALLOC(node_anal, double*, MODL->body[ibody].nnode+1);
        for (inode = 0; inode <= MODL->body[ibody].nnode; inode++) {
            node_anal[inode] = NULL;
        }

        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
            npnt_tess = 1;

            MALLOC(node_anal[inode], double, 3*npnt_tess);

            old_time = clock();
            status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, npnt_tess, NULL, node_anal[inode]);
            CHECK_STATUS(ocsmGetVel);
            new_time = clock();
            geom_time += (new_time - old_time);

            for (ixyz = 0; ixyz < 3; ixyz++) {
                if (node_anal[inode][ixyz] != node_anal[inode][ixyz]) {
                    node_anal[inode][ixyz] = HUGEQ;
                }
            }
        }

        /* if there is a perturbation, return that finite differences were used */
        if (MODL->perturb != NULL) {
            status = EXIT_SUCCESS;
            goto cleanup;
        }

        /* finite difference geometric sensitivities */
        SPRINT1(0, "Computing finite difference sensitivities for ibody=%d",
                ibody);

        status = ocsmSetDtime(MODL, dtime);
        CHECK_STATUS(ocsmSetDtime);

        /* save finite difference geometric sensitivity of each Face */
        MALLOC(face_fdif, double*, MODL->body[ibody].nface+1);
        for (iface = 0; iface <= MODL->body[ibody].nface; iface++) {
            face_fdif[iface] = NULL;
        }

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                    &ntri_tess, &tris, &tric);
            CHECK_STATUS(EG_getTessFace);

            MALLOC(face_fdif[iface], double, 3*npnt_tess);

            status = ocsmGetVel(MODL, ibody, OCSM_FACE, iface, npnt_tess, NULL, face_fdif[iface]);
            CHECK_STATUS(ocsmGetVel);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (face_fdif[iface][ixyz] != face_fdif[iface][ixyz]) {
                    face_fdif[iface][ixyz] = -HUGEQ;
                }
            }
        }

        /* save finite difference geometric sensitivity of each Edge */
        MALLOC(edge_fdif, double*, MODL->body[ibody].nedge+1);
        for (iedge = 0; iedge <= MODL->body[ibody].nedge; iedge++) {
            edge_fdif[iedge] = NULL;
        }

        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
            status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                    &npnt_tess, &xyz, &uv);
            CHECK_STATUS(EG_getTessEdge);

            MALLOC(edge_fdif[iedge], double, 3*npnt_tess);

            status = ocsmGetVel(MODL, ibody, OCSM_EDGE, iedge, npnt_tess, NULL, edge_fdif[iedge]);
            CHECK_STATUS(ocsmGetVel);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (edge_fdif[iedge][ixyz] != edge_fdif[iedge][ixyz]) {
                    edge_fdif[iedge][ixyz] = -HUGEQ;
                }
            }
        }

        /* save finite difference geometric sensitivity of each Node */
        MALLOC(node_fdif, double*, MODL->body[ibody].nnode+1);
        for (inode = 0; inode <= MODL->body[ibody].nnode; inode++) {
            node_fdif[inode] = NULL;
        }

        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
            npnt_tess = 1;

            MALLOC(node_fdif[inode], double, 3*npnt_tess);

            status = ocsmGetVel(MODL, ibody, OCSM_NODE, inode, npnt_tess, NULL, node_fdif[inode]);
            CHECK_STATUS(ocsmGetVel);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (node_fdif[inode][ixyz] != node_fdif[inode][ixyz]) {
                    node_fdif[inode][ixyz] = -HUGEQ;
                }
            }
        }

        status = ocsmSetDtime(MODL, 0);
        CHECK_STATUS(ocsmSetDtime);

        /* compare geometric sensitivities */
        node_errmax = 0;
        edge_errmax = 0;
        face_errmax = 0;

        if (MODL->pmtr[ipmtr].nrow == 1 &&
            MODL->pmtr[ipmtr].ncol == 1   ) {
            SPRINT2(0, "\nComparing geometric sensitivities wrt \"%s\" for ibody=%d",
                    MODL->pmtr[ipmtr].name, ibody);
        } else {
            SPRINT4(0, "\nComparing geometric sensitivities wrt \"%s[%d,%d]\" for ibody=%d",
                    MODL->pmtr[ipmtr].name, irow, icol, ibody);
        }

        /* print velocities for each Face, Edge, and Node */
        if (showAll == 1) {
            SPRINT0(0, "              ipnt     X_anal       X_fdif        Y_anal       Y_fdif        Z_anal       Z_fdif          error");
            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                        &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                        &ntri_tess, &tris, &tric);
                CHECK_STATUS(EG_getTessFace);

                for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                    errX = fabs(face_anal[iface][3*ipnt  ]-face_fdif[iface][3*ipnt  ]);
                    errY = fabs(face_anal[iface][3*ipnt+1]-face_fdif[iface][3*ipnt+1]);
                    errZ = fabs(face_anal[iface][3*ipnt+2]-face_fdif[iface][3*ipnt+2]);
                    SPRINT10(0, "Face %3d:%-3d %5d  %12.7f %12.7f  %12.7f %12.7f  %12.7f %12.7f  %15.10f",
                             ibody, iface, ipnt, face_anal[iface][3*ipnt  ], face_fdif[iface][3*ipnt  ],
                                                 face_anal[iface][3*ipnt+1], face_fdif[iface][3*ipnt+1],
                                                 face_anal[iface][3*ipnt+2], face_fdif[iface][3*ipnt+2],
                             MAX(MAX(errX,errY),errZ));
                }
            }
            for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                        &npnt_tess, &xyz, &uv);
                CHECK_STATUS(EG_getTessEdge);

                for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                    errX = fabs(edge_anal[iedge][3*ipnt  ]-edge_fdif[iedge][3*ipnt  ]);
                    errY = fabs(edge_anal[iedge][3*ipnt+1]-edge_fdif[iedge][3*ipnt+1]);
                    errZ = fabs(edge_anal[iedge][3*ipnt+2]-edge_fdif[iedge][3*ipnt+2]);
                    SPRINT10(0, "Edge %3d:%-3d %5d  %12.7f %12.7f  %12.7f %12.7f  %12.7f %12.7f  %15.10f",
                             ibody, iedge, ipnt, edge_anal[iedge][3*ipnt  ], edge_fdif[iedge][3*ipnt  ],
                                                 edge_anal[iedge][3*ipnt+1], edge_fdif[iedge][3*ipnt+1],
                                                 edge_anal[iedge][3*ipnt+2], edge_fdif[iedge][3*ipnt+2],
                             MAX(MAX(errX,errY),errZ));
                }
            }
            for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                ipnt = 0;
                errX = fabs(node_anal[inode][3*ipnt  ]-node_fdif[inode][3*ipnt  ]);
                errY = fabs(node_anal[inode][3*ipnt+1]-node_fdif[inode][3*ipnt+1]);
                errZ = fabs(node_anal[inode][3*ipnt+2]-node_fdif[inode][3*ipnt+2]);
                SPRINT10(0, "Node %3d:%-3d %5d  %12.7f %12.7f  %12.7f %12.7f  %12.7f %12.7f  %15.10f",
                         ibody, inode, ipnt, node_anal[inode][3*ipnt  ], node_fdif[inode][3*ipnt  ],
                                             node_anal[inode][3*ipnt+1], node_fdif[inode][3*ipnt+1],
                                             node_anal[inode][3*ipnt+2], node_fdif[inode][3*ipnt+2],
                         MAX(MAX(errX,errY),errZ));
            }
        }

        /* compare geometric sensitivities for interior points in each Face */
        nerror = 0;
        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                    &ntri_tess, &tris, &tric);
            CHECK_STATUS(EG_getTessFace);

            status = EG_attributeRet(MODL->body[ibody].face[iface].eface, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                SPRINT2(0, "Tests suppressed for ibody=%3d, iface=%3d", ibody, iface);
                FREE(face_anal[iface]);
                FREE(face_fdif[iface]);
                (*nsuppress)++;
                continue;
            } else {
                status = EGADS_SUCCESS;
            }

            errrms = 0;
            nrms   = 0;
            for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                if (ptype[ipnt] >= 0) continue;

                for (i = 0; i < 3; i++) {
                    face_err = face_anal[iface][3*ipnt+i] - face_fdif[iface][3*ipnt+i];
                    errrms  += face_err * face_err;
                    nrms++;
                }
            }
            if (nrms > 0) errrms = sqrt(errrms / nrms);

            if (errrms > face_errmax) face_errmax = errrms;
            if (errrms > ERROR_TOLER) (*ntotal)++;
            if (errrms > ERROR_REPORT) {
                SPRINT3(1, "      Face %4d:%-4d has errrms=%12.5e", ibody, iface, MAX(errrms, EPS20));
            }

            FREE(face_anal[iface]);
            FREE(face_fdif[iface]);
        }
        (*errmax) = MAX((*errmax), face_errmax);
        SPRINT3(0, "    d(Face)/d(%s) check complete with %8d total errors (errmax=%12.4e)",
                MODL->pmtr[ipmtr].name, nerror, face_errmax);

        /* compare geometric sensitivities for interior points in each Edge */
        nerror = 0;
        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
            status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                    &npnt_tess, &xyz, &uv);
            CHECK_STATUS(EG_getTessEdge);

            status = EG_attributeRet(MODL->body[ibody].edge[iedge].eedge, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                SPRINT2(0, "Tests suppressed for ibody=%3d, iedge=%3d", ibody, iedge);
                FREE(edge_anal[iedge]);
                FREE(edge_fdif[iedge]);
                (*nsuppress)++;
                continue;
            } else {
                status = EGADS_SUCCESS;
            }

            errrms = 0;
            nrms   = 0;
            for (ipnt = 1; ipnt < npnt_tess-1; ipnt++) {
                for (i = 0; i < 3; i++) {
                    edge_err = edge_anal[iedge][3*ipnt+i] - edge_fdif[iedge][3*ipnt+i];
                    errrms  += edge_err * edge_err;
                }
            }
            if (nrms > 0) errrms = sqrt(errrms / nrms);

            if (errrms > edge_errmax) edge_errmax = errrms;
            if (errrms > ERROR_TOLER) (*ntotal)++;
            if (errrms > ERROR_REPORT) {
                SPRINT3(1, "      Edge %4d:%-4d has errrms=%12.5e", ibody, iedge, MAX(errrms, EPS20));
            }

            FREE(edge_anal[iedge]);
            FREE(edge_fdif[iedge]);
        }
        (*errmax) = MAX((*errmax), edge_errmax);
        SPRINT3(0, "    d(Edge)/d(%s) check complete with %8d total errors (errmax=%12.4e)",
                MODL->pmtr[ipmtr].name, nerror, edge_errmax);

        /* compare geometric sensitivities for each Node */
        nerror = 0;
        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
            npnt_tess = 1;

            status = EG_attributeRet(MODL->body[ibody].node[inode].enode, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                SPRINT2(0, "Tests suppressed for ibody=%3d, inode=%3d", ibody, inode);
                FREE(node_anal[inode]);
                FREE(node_fdif[inode]);
                (*nsuppress)++;
                continue;
            } else {
                status = EGADS_SUCCESS;
            }

            errrms = 0;
            nrms   = 0;
            for (i = 0; i < 3; i++) {
                node_err = node_anal[inode][i] - node_fdif[inode][i];
                errrms  += node_err * node_err;
                nrms++;
            }
            if (nrms > 0) errrms = sqrt(errrms / nrms);

            if (errrms > node_errmax) node_errmax = errrms;
            if (errrms > ERROR_TOLER) (*ntotal)++;
            if (errrms > ERROR_REPORT) {
                SPRINT3(1, "      Node %4d:%-4d has errrms=%12.5e", ibody, inode, MAX(errrms, EPS20));
            }

            FREE(node_anal[inode]);
            FREE(node_fdif[inode]);
        }
        (*errmax) = MAX((*errmax), node_errmax);
        SPRINT3(0, "    d(Node)/d(%s) check complete with %8d total errors (errmax=%12.4e)",
                MODL->pmtr[ipmtr].name, nerror, node_errmax);

        FREE(face_anal);
        FREE(face_fdif);
        FREE(edge_anal);
        FREE(edge_fdif);
        FREE(node_anal);
        FREE(node_fdif);
    }

cleanup:
    /* these FREEs only get used if an error was encountered above */
    if (ibody <= MODL->nbody) {
        if (face_anal != NULL) {
            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                FREE(face_anal[iface]);
            }
        }
        if (face_fdif != NULL) {
            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                FREE(face_fdif[iface]);
            }
        }
        if (edge_anal != NULL) {
            for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                FREE(edge_anal[iedge]);
            }
        }
        if (edge_fdif != NULL) {
            for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                FREE(edge_fdif[iedge]);
            }
        }
        if (node_anal != NULL) {
            for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                FREE(node_anal[inode]);
            }
        }
        if (node_fdif != NULL) {
            for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                FREE(node_fdif[inode]);
            }
        }
    }

    FREE(face_anal);
    FREE(face_fdif);
    FREE(edge_anal);
    FREE(edge_fdif);
    FREE(node_anal);
    FREE(node_fdif);

    return status;
}


/***********************************************************************/
/*                                                                     */
/*   checkTessSens - check tessellation sensitivities                  */
/*                                                                     */
/*                     this is done by computing the distance of each  */
/*                     tessellation point plus its sensitivity from    */
/*                     a purturbed surface                             */
/*                                                                     */
/***********************************************************************/

static int
checkTessSens(int    ipmtr,             /* (in)  Parameter index (bias-1) */
              int    irow,              /* (in)  row    index (bias-1) */
              int    icol,              /* (in)  column index (bias-1) */
              int    *ntotal,           /* (out) ntotal number of points beyond toler */
              int    *nsuppress,        /* (out) number of suppressions */
              double *errmax)           /* (out) maximum error */
{
    int       status = SUCCESS;

    int       ibody, atype, alen, nbody, builtTo, nbad, nrms;
    int       nerror=0, inode, jnode, iedge, jedge, iface, jface, ipnt, itri, ixyz;
    int       npnt_tess, ntri_tess, oclass, mtype, nchild, *senses;
    CINT      *ptype, *pindx, *tris, *tric, *tempIlist, *nMap, *eMap, *fMap;
    double    data[18], old_value, old_dot, scaled_dtime, uv_clos[2], xyz_clos[18], dist, val, *Dist=NULL;
    double    errrms, face_errmax, edge_errmax, node_errmax;
    double    **face_ptrb=NULL, **edge_ptrb=NULL, **node_ptrb=NULL;
    CDOUBLE   *tempRlist, *xyz, *uv, *t, *dxyz;
    CCHAR     *tempClist;
    ego       eref, *echilds, eNewBody=NULL;
    clock_t   old_time, new_time;
    FILE      *fp_badtri=NULL, *fp_logdist=NULL;
    void      *ptrb;

    modl_T    *MODL = (modl_T *)modl;
    modl_T    *PTRB = NULL;

    ROUTINE(checkTessSens);

    /* --------------------------------------------------------------- */

    SPRINT0(0, "\n*********************************************************");
    if (MODL->pmtr[ipmtr].nrow == 1 &&
        MODL->pmtr[ipmtr].ncol == 1   ) {
        SPRINT1(0, "Starting tessellation sensitivity wrt \"%s\"",
                MODL->pmtr[ipmtr].name);
    } else {
        SPRINT3(0, "Starting tessellation sensitivity wrt \"%s[%d,%d]\"",
                MODL->pmtr[ipmtr].name, irow, icol);
    }
    SPRINT0(0, "*********************************************************\n");

    ibody = MODL->nbody + 1;            /* protect for warning in cleanup: */

    /* start the badtri file */
    fp_badtri = fopen("bad.triangles", "w");
    if (fp_badtri != NULL) {
        fprintf(fp_badtri, "xxxxx   -2 badTriangles\n");
        nbad = 0;
    }

    /* start the logdist file if we are only looking at one DESPMTR */
    if (strlen(pmtrname) > 0) {
        SPRINT0(0, "\n*****************************************");
        SPRINT0(0, "logdist.plot file is being generated");
        SPRINT0(0, "    color scheme: blue=1e-12 (or smaller)");
        SPRINT0(0, "                  red =1e-03 (or larger)");
        SPRINT0(0, "*****************************************\n");

        fp_logdist = fopen("logdist.plot", "w");
    }

    /* set the velocity */
    status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
    CHECK_STATUS(ocsmSetVelD);
    status = ocsmSetVelD(MODL, ipmtr, irow, icol, 1.0);
    CHECK_STATUS(ocsmSetVelD);

    /* build needed to propagate velocities to the Branch arguments */
    nbody  = 0;
    status = ocsmBuild(MODL, 0, &builtTo, &nbody, NULL);
    CHECK_STATUS(ocsmBuild);

    /* make a copy of the MODL */
    status = ocsmCopy(MODL, &ptrb);
    CHECK_STATUS(ocsmCopy);

    PTRB = (modl_T *)ptrb;

    PTRB->tessAtEnd = 0;

    /* perturb the configuration */
    status = ocsmGetValu(PTRB, ipmtr, irow, icol, &old_value, &old_dot);
    CHECK_STATUS(ocsmGetValu);

    scaled_dtime = dtime * MAX(fabs(old_value), +1.0);

    status = ocsmSetValuD(PTRB, ipmtr, irow, icol, old_value+scaled_dtime);
    CHECK_STATUS(ocsmSetValuD);

    SPRINT4(0, "Generating perturbed configuration with delta-%s[%d,%d]=%13.8f",
            MODL->pmtr[ipmtr].name, irow, icol, scaled_dtime);

    /* build the perturbed configuration */
    nbody  = 0;
    status = ocsmBuild(ptrb, 0, &builtTo, &nbody, NULL);
    if (status != SUCCESS) {
        SPRINT0(1, "ERROR:: tess error: perturbed configuration could not be built");
        SPRINT0(1, "        tessellation sensitivity checks skipped for all Bodys");
        nerror++;
        status = SUCCESS;
        goto cleanup;
    }

    CHECK_STATUS(ocsmBuild);

    /* initialize the maximum errors */
    face_errmax = 0;
    edge_errmax = 0;
    node_errmax = 0;

    /* perform tessellation sensitivity checks for each Body on the stack */
    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
        if (MODL->body[ibody].onstack != 1) continue;

        /* make sure perturbed and base models have the same number of Nodes, Edges, and Faces */
        if        (MODL->body[ibody].nface != PTRB->body[ibody].nface) {
            SPRINT3(1, "ERROR:: tess error: perturbed Body %d has different .nface (%d vs %d)",
                    ibody, MODL->body[ibody].nface, PTRB->body[ibody].nface);
            SPRINT0(1, "        tessellation sensitivity checks skipped for this Body");
            continue;
        } else if (MODL->body[ibody].nedge != PTRB->body[ibody].nedge) {
            SPRINT3(1, "ERROR:: tess error: perturbed Body %d has different .nedge (%d vs %d)",
                    ibody, MODL->body[ibody].nedge, PTRB->body[ibody].nedge);
            SPRINT0(1, "        tessellation sensitivity checks skipped for this Body");
            continue;
        } else if (MODL->body[ibody].nnode != PTRB->body[ibody].nnode) {
            SPRINT3(1, "ERROR:: tess error: perturbed Body %d has different .nnode (%d vs %d)",
                    ibody, MODL->body[ibody].nnode, PTRB->body[ibody].nnode);
            SPRINT0(1, "        tessellation sensitivity checks skipped for this Body");
            continue;
        }

        /* analytic tessellation sensitivities (if possible) */
        SPRINT1(0, "Computing analytic sensitivities (if possible) for ibody=%d",
                ibody);

        status = ocsmSetDtime(MODL, 0);
        CHECK_STATUS(ocsmSetDtime);

        /* determine if there is a mapping between the MODL and the PTRB */
        if (MODL->body[ibody].nface > 0) {
            status = EG_mapBody(MODL->body[ibody].ebody, PTRB->body[ibody].ebody, "_faceID", &eNewBody);
            CHECK_STATUS(EG_mapBody);

            if (eNewBody == NULL) {
                nMap = NULL;
                eMap = NULL;
                fMap = NULL;
            } else {
                status = EG_attributeRet(eNewBody, ".nMap",
                                         &atype, &alen, &nMap, &tempRlist, &tempClist);
                if (status != EGADS_SUCCESS) {
                    nMap = NULL;
                }

                status = EG_attributeRet(eNewBody, ".eMap",
                                         &atype, &alen, &eMap, &tempRlist, &tempClist);
                if (status != EGADS_SUCCESS) {
                    eMap = NULL;
                }

                status = EG_attributeRet(eNewBody, ".fMap",
                                     &atype, &alen, &fMap, &tempRlist, &tempClist);
                if (status != EGADS_SUCCESS) {
                    fMap = NULL;
                }
            }
        } else {
            nMap = NULL;
            eMap = NULL;
            fMap = NULL;
        }

        /* save perturbed points (which are the base plus their sensitivities) for each Face */
        MALLOC(face_ptrb, double*, MODL->body[ibody].nface+1);
        for (iface = 0; iface <= MODL->body[ibody].nface; iface++) {
            face_ptrb[iface] = NULL;
        }

        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                    &ntri_tess, &tris, &tric);
            CHECK_STATUS(EG_getTessFace);
            if (npnt_tess <= 0) {
                SPRINT3(0, "ERROR:: EG_getTessFace -> status=%d (%s), npnt_tess=%d",
                        status, ocsmGetText(status), npnt_tess);
                status = EXIT_FAILURE;
                goto cleanup;
            }

            old_time = clock();
            status = ocsmGetTessVel(MODL, ibody, OCSM_FACE, iface, &dxyz);
            CHECK_STATUS(ocsmGetTessVel);
            new_time = clock();
            tess_time += (new_time - old_time);

            MALLOC(face_ptrb[iface], double, 3*npnt_tess);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (dxyz[ixyz] == dxyz[ixyz]) {
                    face_ptrb[iface][ixyz] = xyz[ixyz] + scaled_dtime * dxyz[ixyz];
                } else {
                    face_ptrb[iface][ixyz] = HUGEQ;
                }
            }

            status = EG_attributeRet(MODL->body[ibody].face[iface].eface, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                SPRINT2(0, "Tests suppressed for ibody=%3d, iface=%3d", ibody, iface);
                FREE(face_ptrb[iface]);
                (*nsuppress)++;
                continue;
            }
        }

        /* save perturbed points (which are the base plus their sensitivities) for each Edge */
        MALLOC(edge_ptrb, double*, MODL->body[ibody].nedge+1);
        for (iedge = 0; iedge <= MODL->body[ibody].nedge; iedge++) {
            edge_ptrb[iedge] = NULL;
        }

        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
            status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                    &npnt_tess, &xyz, &uv);
            CHECK_STATUS(EG_getTessEdge);
            if (npnt_tess <= 0) {
                SPRINT3(0, "ERROR:: EG_getTessEdge -> status=%d (%s), npnt_tess=%d",
                        status, ocsmGetText(status), npnt_tess);
                status = EXIT_FAILURE;
                goto cleanup;
            }

            old_time = clock();
            status = ocsmGetTessVel(MODL, ibody, OCSM_EDGE, iedge, &dxyz);
            CHECK_STATUS(ocsmGetTessVel);
            new_time = clock();
            tess_time += (new_time - old_time);

            MALLOC(edge_ptrb[iedge], double, 3*npnt_tess);

            for (ixyz = 0; ixyz < 3*npnt_tess; ixyz++) {
                if (dxyz[ixyz] == dxyz[ixyz]) {
                    edge_ptrb[iedge][ixyz] = xyz[ixyz] + scaled_dtime * dxyz[ixyz];
                } else {
                    edge_ptrb[iedge][ixyz] = HUGEQ;
                }
            }

            status = EG_attributeRet(MODL->body[ibody].edge[iedge].eedge, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                SPRINT2(0, "Tests suppressed for ibody=%3d, iedge=%3d", ibody, iedge);
                FREE(edge_ptrb[iedge]);
                (*nsuppress)++;
                continue;
            }
        }

        /* save perturbed points (which are the base plus their sensitivities) for each Node */
        MALLOC(node_ptrb, double*, MODL->body[ibody].nnode+1);
        for (inode = 0; inode <= MODL->body[ibody].nnode; inode++) {
            node_ptrb[inode] = NULL;
        }

        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
            status = EG_getTopology(MODL->body[ibody].node[inode].enode, &eref,
                                    &oclass, &mtype, data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            old_time = clock();
            status = ocsmGetTessVel(MODL, ibody, OCSM_NODE, inode, &dxyz);
            CHECK_STATUS(ocsmGetTessVel);
            new_time = clock();
            tess_time += (new_time - old_time);

            MALLOC(node_ptrb[inode], double, 3);

            for (ixyz = 0; ixyz < 3; ixyz++) {
                if (dxyz[ixyz] == dxyz[ixyz]) {
                    node_ptrb[inode][ixyz] = data[ixyz] + scaled_dtime * dxyz[ixyz];
                } else {
                    node_ptrb[inode][ixyz] = HUGEQ;
                }
            }
        }

        /* if there is a perturbation, return that finite differences were used */
        if (MODL->perturb != NULL) {
            status = EXIT_SUCCESS;
            goto cleanup;
        }

        SPRINT1(0, "Computing distances of perturbed points from perturbed configuration for ibody=%d",
                ibody);

        /* check how far each perturbed Face point is from the perturbed configuration */
        nerror = 0;
        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                    &ntri_tess, &tris, &tric);
            CHECK_STATUS(EG_getTessFace);

            if (npnt_tess <= 0) {
                SPRINT3(0, "ERROR:: EG_getTessFace -> status=%d (%s), npnt_tess=%d",
                        status, ocsmGetText(status), npnt_tess);
                status = EXIT_FAILURE;
                goto cleanup;
            }

            /* find the equivalent Face in PTRB */
            if (fMap == NULL) {
                jface = iface;
            } else {
                jface = fMap[iface-1];
            }

            status = EG_attributeRet(MODL->body[ibody].face[iface].eface, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                continue;
            }

            status = EG_getTessFace(MODL->body[ibody].etess, iface,
                                    &npnt_tess, &xyz, &uv, &ptype, &pindx,
                                    &ntri_tess, &tris, &tric);
            CHECK_STATUS(EG_getTessFace);

            MALLOC(Dist, double, npnt_tess);

            errrms = 0;
            nrms   = 0;
            for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                if (ptype[ipnt] >= 0) continue;

                uv_clos[ 0] = uv[2*ipnt  ];
                uv_clos[ 1] = uv[2*ipnt+1];
                status = EG_invEvaluateGuess(PTRB->body[ibody].face[jface].eface,
                                             &(face_ptrb[iface][3*ipnt]), uv_clos, xyz_clos);
                CHECK_STATUS(EG_invEvaluateGuess);

                Dist[ipnt] = sqrt((face_ptrb[iface][3*ipnt  ] - xyz_clos[0]) * (face_ptrb[iface][3*ipnt  ] - xyz_clos[0])
                                + (face_ptrb[iface][3*ipnt+1] - xyz_clos[1]) * (face_ptrb[iface][3*ipnt+1] - xyz_clos[1])
                                + (face_ptrb[iface][3*ipnt+2] - xyz_clos[2]) * (face_ptrb[iface][3*ipnt+2] - xyz_clos[2]));
                if (Dist[ipnt] != Dist[ipnt]) {
                    printf("Dist[%d] = nan\n", ipnt);
                    Dist[ipnt] = 1e+99;
                }

                errrms += Dist[ipnt] * Dist[ipnt];
                nrms++;
            }
            if (nrms > 0) errrms = sqrt(errrms / nrms);
            FREE(face_ptrb[iface]);

            if (errrms > face_errmax) face_errmax = errrms;
            if (errrms > ERROR_TOLER) (*ntotal)++;
            if (errrms > ERROR_REPORT) {
                SPRINT3(1, "      Face %4d:%-4d has errrms=%12.5e", ibody, iface, MAX(errrms, EPS20));
            }

            /* add the distances to this Face to the logdist file */
            /*     -1 (red)   corresponds to Dist < 1e-12 */
            /*     +1 (blue)  corresponds to Dist > 1e-03 */

            if (fp_logdist != NULL) {
                fprintf(fp_logdist, "%5d%5d Body_%d_%d\n", ntri_tess, -3, ibody, iface);

                for (itri = 0; itri < ntri_tess; itri++) {
                    ipnt = tris[3*itri  ] - 1;
                    val  = MIN((7.5 + log10(MAX(Dist[ipnt], 1e-12))) / 4.5, 1);
                    fprintf(fp_logdist, "%15.8f %15.8f %15.8f %15.8f ",  xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2], val);

                    ipnt = tris[3*itri+1] - 1;
                    val  = MIN((7.5 + log10(MAX(Dist[ipnt], 1e-12))) / 4.5, 1);
                    fprintf(fp_logdist, "%15.8f %15.8f %15.8f %15.8f ",  xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2], val);

                    ipnt = tris[3*itri+2] - 1;
                    val  = MIN((7.5 + log10(MAX(Dist[ipnt], 1e-12))) / 4.5, 1);
                    fprintf(fp_logdist, "%15.8f %15.8f %15.8f %15.8f\n", xyz[3*ipnt], xyz[3*ipnt+1], xyz[3*ipnt+2], val);
                }
            }

            FREE(Dist);
        }
        (*errmax) = MAX((*errmax), face_errmax);
        SPRINT3(0, "    d(Face)/d(%s) check complete with %8d total errors (errmax=%12.4e)",
                MODL->pmtr[ipmtr].name, nerror, face_errmax);

        /* check how far each perturbed Edge point is from the perturbed configuration */
        nerror = 0;
        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {

            /* find the equivalent Edge in PTRB */
            if (eMap == NULL) {
                jedge = iedge;
            } else {
                jedge = eMap[iedge-1];
            }

            status = EG_attributeRet(MODL->body[ibody].edge[iedge].eedge, "_sensCheck",
                                     &atype, &alen, &tempIlist, &tempRlist, &tempClist);
            if (status == EGADS_SUCCESS && atype == ATTRSTRING && strcmp(tempClist, "skip") == 0) {
                FREE(edge_ptrb[iedge]);
                continue;
            }

            status = EG_getTopology(MODL->body[ibody].edge[iedge].eedge, &eref, &oclass, &mtype,
                                    data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            if (mtype == DEGENERATE) {
                FREE(edge_ptrb[iedge]);
                continue;
            }

            status = EG_getTessEdge(MODL->body[ibody].etess, iedge,
                                    &npnt_tess, &xyz, &t);
            CHECK_STATUS(EG_getTessEdge);

            errrms = 0;
            nrms   = 0;
            for (ipnt = 1; ipnt < npnt_tess-1; ipnt++) {
                uv_clos[ 0] = t[ipnt];
                status = EG_invEvaluateGuess(PTRB->body[ibody].edge[jedge].eedge,
                                             &(edge_ptrb[iedge][3*ipnt]), uv_clos, xyz_clos);
                CHECK_STATUS(EG_invEvaluateGuess);

                dist = sqrt((edge_ptrb[iedge][3*ipnt  ] - xyz_clos[0]) * (edge_ptrb[iedge][3*ipnt  ] - xyz_clos[0])
                          + (edge_ptrb[iedge][3*ipnt+1] - xyz_clos[1]) * (edge_ptrb[iedge][3*ipnt+1] - xyz_clos[1])
                          + (edge_ptrb[iedge][3*ipnt+2] - xyz_clos[2]) * (edge_ptrb[iedge][3*ipnt+2] - xyz_clos[2]));
                errrms += dist * dist;
                nrms++;
            }
            if (nrms > 0) errrms = sqrt(errrms / nrms);
            FREE(edge_ptrb[iedge]);

            if (errrms > edge_errmax) edge_errmax = errrms;
            if (errrms > ERROR_TOLER) (*ntotal)++;
            if (errrms > ERROR_REPORT) {
                SPRINT3(1, "      Edge %4d:%-4d has errrms=%12.5e", ibody, iedge, MAX(errrms, EPS20));
            }
        }
        (*errmax) = MAX((*errmax), edge_errmax);
        SPRINT3(0, "    d(Edge)/d(%s) check complete with %8d total errors (errmax=%12.4e)",
                MODL->pmtr[ipmtr].name, nerror, edge_errmax);

        /* check how far each perturbed Node point is from the perturbed configuration */
        nerror = 0;
        for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {

            /* find the equivalent Node in PTRB */
            if (nMap == NULL) {
                jnode = inode;
            } else {
                jnode = nMap[inode-1];
            }

            status = EG_getTopology(PTRB->body[ibody].node[jnode].enode, &eref,
                                    &oclass, &mtype, data, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            status = EG_getTopology(PTRB->body[ibody].node[jnode].enode, &eref,
                                    &oclass, &mtype, xyz_clos, &nchild, &echilds, &senses);
            CHECK_STATUS(EG_getTopology);

            dist = sqrt((node_ptrb[inode][0] - data[0]) * (node_ptrb[inode][0] - data[0])
                      + (node_ptrb[inode][1] - data[1]) * (node_ptrb[inode][1] - data[1])
                      + (node_ptrb[inode][2] - data[2]) * (node_ptrb[inode][2] - data[2]));
            errrms = dist;

            if (errrms > node_errmax) node_errmax = errrms;
            if (errrms > ERROR_TOLER) (*ntotal)++;
            if (errrms > ERROR_REPORT) {
                SPRINT3(1, "      Node %4d:%-4d has errrms=%12.5e", ibody, inode, MAX(errrms, EPS20));
            }

            FREE(node_ptrb[inode]);
        }
        (*errmax) = MAX((*errmax), node_errmax);
        SPRINT3(0, "    d(Node)/d(%s) check complete with %8d total errors (errmax=%12.4e)",
                MODL->pmtr[ipmtr].name, nerror, node_errmax);

        FREE(face_ptrb);
        FREE(edge_ptrb);
        FREE(node_ptrb);

        if (eNewBody != NULL) {
            status = EG_deleteObject(eNewBody);
            CHECK_STATUS(EG_deleteObject);
        }
    }

cleanup:
    if (fp_badtri != NULL) {
        fprintf(fp_badtri, "    0    0 end (of %d triangles)\n", nbad);
        fclose(fp_badtri);
    }

    if (fp_logdist != NULL) {
        fprintf(fp_logdist, "    0    0 end\n");
        fclose(fp_logdist);
    }

    /* these FREEs only get used if an error was encountered above */
    if (ibody <= MODL->nbody) {
        if (face_ptrb != NULL) {
            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                FREE(face_ptrb[iface]);
            }
        }
        if (edge_ptrb != NULL) {
            for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                FREE(edge_ptrb[iedge]);
            }
        }
        if (node_ptrb != NULL) {
            for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                FREE(node_ptrb[inode]);
            }
        }
    }

    if (PTRB != NULL) {
        (void) ocsmFree(PTRB);
    }

    FREE(face_ptrb);
    FREE(edge_ptrb);
    FREE(node_ptrb);
    FREE(Dist);

    return status;
}
