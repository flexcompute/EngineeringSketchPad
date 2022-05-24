/*
 ************************************************************************
 *                                                                      *
 * buildCSM.c -- use OpenCSM code to build Bodys                        *
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

#ifdef WIN32
    #include <windows.h>
#endif

#include "egads.h"

#define CINT    const int
#define CDOUBLE const double
#define CCHAR   const char

#define STRNCPY(A, B, LEN) strncpy(A, B, LEN); A[LEN-1] = '\0';

#include "common.h"
#include "OpenCSM.h"
#include "udp.h"

#include "gv.h"
#include "Graphics.h"

/***********************************************************************/
/*                                                                     */
/* macros (including those that go along with common.h)                */
/*                                                                     */
/***********************************************************************/

static int outLevel = 1;

                                               /* used by RALLOC macro */
//static void *realloc_temp=NULL;

#define  RED(COLOR)      (float)(COLOR / 0x10000        ) / (float)(255)
#define  GREEN(COLOR)    (float)(COLOR / 0x00100 % 0x100) / (float)(255)
#define  BLUE(COLOR)     (float)(COLOR           % 0x100) / (float)(255)

/***********************************************************************/
/*                                                                     */
/* structures                                                          */
/*                                                                     */
/***********************************************************************/


/***********************************************************************/
/*                                                                     */
/* global variables                                                    */
/*                                                                     */
/***********************************************************************/

/* global variable holding a MODL */
       char      casename[255];        /* name of case */
static void      *modl;                /* pointer to MODL */

/* global variables associated with graphical user interface (gui) */
static int        ngrobj   = 0;        /* number of GvGraphic objects */
static GvGraphic  **grobj;             /* list   of GvGraphic objects */
static int        new_data = 1;        /* =1 means image needs to be updated */
static FILE       *script  = NULL;     /* pointer to script file (if it exists) */
static int        numarg   = 0;        /* numeric argument */
static int        fly_mode = 0;        /* =0 arrow keys rotate, =1 arrow keys translate */
static int        sclrType = -1;       /* =-1 for monochrome, =0 for u or v, =1 for velocity */
static double     tuftLen  = 0;        /* length of tuft (=0 for no tuft) */
static double     bigbox[6];           /* bounding box of config */
static int        builtTo  = 0;        /* last branch built to */

/* global variables associated with paste buffer */
#define MAX_PASTE         10

static int        npaste   = 0;                         /* number of Branches in paste buffer */
static int        type_paste[MAX_PASTE];                /* type for Branches in paste buffer */
static char       name_paste[MAX_PASTE][MAX_EXPR_LEN];  /* name for Branches in paste buffer */
static char       str1_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str1 for Branches in paste buffer */
static char       str2_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str2 for Branches in paste buffer */
static char       str3_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str3 for Branches in paste buffer */
static char       str4_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str4 for Branches in paste buffer */
static char       str5_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str5 for Branches in paste buffer */
static char       str6_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str6 for Branches in paste buffer */
static char       str7_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str7 for Branches in paste buffer */
static char       str8_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str8 for Branches in paste buffer */
static char       str9_paste[MAX_PASTE][MAX_EXPR_LEN];  /* str9 for Branches in paste buffer */

/* global variables associated with bodyList */
#define MAX_BODYS        999
static int        nbody    = 0;        /* number of built Bodys */
static int        bodyList[MAX_BODYS]; /* array  of built Bodys */

/***********************************************************************/
/*                                                                     */
/* declarations                                                        */
/*                                                                     */
/***********************************************************************/

/* declarations for high-level routines defined below */

/* declarations for support routines defined below */
static int    evalRuled(modl_T *MODL, int ibody, int iface, int isketchs, int isketchn, int iedge, double uv[], double xyz[]);

/* declarations for graphic routines defined below */
       int        gvupdate();
       void       gvdata(int ngraphics, GvGraphic *graphic[]);
       int        gvscalar(int key, GvGraphic *graphic, int len, float *scalar);
       void       gvevent(int *win, int *type, int *xpix, int *ypix, int *state);
       void       transform(double xform[3][4], double *point, float *out);
static int        pickObject(int *utype);
static int        getInt(char *prompt);
static double     getDbl(char *prompt);
static void       getStr(char *prompt, char *string);

/* external functions not declared in an .h file */

/* functions and variables associated with picking */
extern GvGraphic  *gv_picked;
extern int        gv_pickon;
extern int        gv_pickmask;
extern void       PickGraphic(float x, float y, int flag);

/* functions and variables associated with background ccolor, etc. */
extern GWinAtt    gv_wAux, gv_wDial, gv_wKey, gv_w3d, gv_w2d;
extern float      gv_black[3], gv_white[3];
extern float      gv_xform[4][4];
extern void       GraphicGCSetFB(GID gc, float *fg, float *bg);
extern void       GraphicBGColor(GID win, float *bg);

/* window defines */
#define           DataBase        1
#define           TwoD            2
#define           ThreeD          3
#define           Dials           4
#define           Key             5

/* event types */
#define           KeyPress        2
#define           KeyRelease      3
#define           ButtonPress     4
#define           ButtonRelease   5
#define           Expose          12
#define           NoExpose        14
#define           ClientMessage   33

/* error codes */
/*                user by CAPRI             -1 to  -99 */
/*                used by OpenCSM         -201 to -199 */


/***********************************************************************/
/*                                                                     */
/*   CAPRI_MAIN - main program                                         */
/*                                                                     */
/***********************************************************************/

int
main(int       argc,                    /* (in)  number of arguments */
     char      *argv[])                 /* (in)  array of arguments */
{

    int       status, mtflag = -1, batch, readonly, recover, i, ibody, jbody;
    int       imajor, iminor, buildTo, nbrch, npmtr, showUsage=0;
    float     focus[4];
    double    box[6];
    char      jnlname[255], filename[255];
    CCHAR     *OCC_ver;
    clock_t   old_time, new_time;

    modl_T    *MODL;
    void      *orig_modl;

    ego       ebody, context;

    int       nkeys      = 3;
    int       keys[3]    = {'u',             'v',             'w'            };
    int       types[3]   = {GV_SURF,         GV_SURF,         GV_SURF        };
    char      titles[48] =  "u Parameter     v Parameter      velocity      " ;
    float     lims[6]    = {0.,1.,           0.,1.,           -1.,1.         };

    ROUTINE(MAIN);

    /* --------------------------------------------------------------- */

    /* get the flags and casename(s) from the command line */
    casename[0] = '\0';
    jnlname[ 0] = '\0';
    batch       = 0;
    readonly    = 0;
    recover     = 0;

    for (i = 1; i < argc; i++) {
        if        (strcmp(argv[i], "-batch") == 0) {
            batch = 1;
        } else if (strcmp(argv[i], "-jnl") == 0) {
            if (i < argc-1) {
                strcpy(jnlname, argv[++i]);
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-readonly") == 0) {
            readonly = 1;
        } else if (strcmp(argv[i], "-recover") == 0) {
            recover  = 1;
        } else if (strcmp(argv[i], "-outLevel") == 0) {
            if (i < argc-1) {
                sscanf(argv[++i], "%d", &outLevel);
                if (outLevel < 0) outLevel = 0;
                if (outLevel > 3) outLevel = 3;
            } else {
                showUsage = 1;
                break;
            }
        } else if (strcmp(argv[i], "-version") == 0 ||
                   strcmp(argv[i], "-v"      ) == 0   ) {
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

    if (showUsage) {
        SPRINT0(0, "proper usage: 'buildCSM [casename[.csm]] [options...]'");
        SPRINT0(0, "   where [options...] = -batch");
        SPRINT0(0, "                        -jnl jnlname");
        SPRINT0(0, "                        -readonly");
        SPRINT0(0, "                        -recover");
        SPRINT0(0, "                        -outLevel X");
        SPRINT0(0, "                        -version  -or-  -v");
        SPRINT0(0, "STOPPING...\a");
        exit(0);
    }

    /* welcome banner */
    (void) ocsmVersion(&imajor, &iminor);

    SPRINT0(1, "**********************************************************");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*                    Program buildCSM                    *");
    SPRINT2(1, "*                     version %2d.%02d                      *", imajor, iminor);
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "*        written by John Dannenhoffer, 2010/2022         *");
    SPRINT0(1, "*                                                        *");
    SPRINT0(1, "**********************************************************");

    /* set OCSMs output level */
    (void) ocsmSetOutLevel(outLevel);

    /* initialize a script (if given) */
    if (strlen(jnlname) > 0) {
        script = fopen(jnlname, "r");
        SPRINT1(1, "Opening script file \"%s\" ...", jnlname);

        if (script == NULL) {
            SPRINT0(0, "ERROR opening script file");
            exit(0);
        }
    }

    /* strip off .csm (which is assumed to be at the end) if present */
    if (strlen(casename) > 0) {
        strcpy(filename, casename);
        if (strstr(casename, ".csm") == NULL) {
            strcat(filename, ".csm");
        }
    } else {
        filename[0] = '\0';
    }

    /* read the .csm file and create the MODL */
    old_time = clock();
    status   = ocsmLoad(filename, &orig_modl);
    new_time = clock();
    SPRINT3(1, "--> ocsmLoad(%s) -> status=%d (%s)", filename, status, ocsmGetText(status));
    SPRINT1(1, "==> ocsmLoad CPUtime=%9.3f sec",
            (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

    if (status < 0) {
        SPRINT0(0, "ERROR:: problem in ocsmLoad");
        goto cleanup;
    }

    /* make a copy of the MODL */
    old_time = clock();
    status   = ocsmCopy(orig_modl, &modl);
    new_time = clock();
    SPRINT2(1, "--> ocsmCopy() -> status=%d (%s)", status, ocsmGetText(status));
    SPRINT1(1, "==> ocsmCopy CPUtime=%9.3f sec",
            (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

    if (status < 0) {
        SPRINT0(0, "ERROR:: problem in ocsmCopy");
        goto cleanup;
    }

    /* delete the original MODL */
    old_time = clock();
    status   = ocsmFree(orig_modl);
    new_time = clock();
    SPRINT2(1, "--> ocsmFree() -> status=%d (%s)", status, ocsmGetText(status));
    SPRINT1(1, "==> ocsmFree CPUtime=%9.3f sec",
            (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

    if (status < 0) {
        SPRINT0(0, "ERROR:: problem in ocsmFree");
        goto cleanup;
    }

    /* check that Branches are properly ordered */
    old_time = clock();
    status   = ocsmCheck(modl);
    new_time = clock();
    SPRINT2(0, "--> ocsmCheck() -> status=%d (%s)",
            status, ocsmGetText(status));
    SPRINT1(0, "==> ocsmCheck CPUtime=%10.3f sec",
            (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

    if (status < 0) {
        SPRINT0(0, "ERROR: problem in ocsmCheck");
        goto cleanup;
    }

    MODL = (modl_T*)modl;

    /* print out the global Attributes, Parameters, and Branches */
    SPRINT0(1, "External Parameter(s):");
    if (outLevel > 0) {
        status = ocsmPrintPmtrs(modl, stdout);
        if (status != SUCCESS) {
            SPRINT1(0, "ocsmPrintPmtrs -> status=%d", status);
        }
    }

    SPRINT0(1, "Branch(es):");
    if (outLevel > 0) {
        status = ocsmPrintBrchs(modl, stdout);
        if (status != SUCCESS) {
            SPRINT1(0, "ocsmPrintBrchs -> status=%d", status);
        }
    }

    SPRINT0(1, "Global Attribute(s):");
    if (outLevel > 0) {
        status = ocsmPrintAttrs(modl, stdout);
        if (status != SUCCESS) {
            SPRINT1(0, "ocsmPrintAttrs -> status=%d", status);
        }
    }

    (void)ocsmInfo(modl, &nbrch, &npmtr, &nbody);

    /* skip ocsmBuild if readonly==1 */
    if (readonly == 1 || nbrch == 0) {
        SPRINT0(0, "WARNING:: ocsmBuild skipped");

    /* build the Bodys from the MODL */
    } else {
        buildTo = 0; /* all */

        /* build the Bodys */
        nbody    = MAX_BODYS;
        old_time = clock();
        status   = ocsmBuild(modl, buildTo, &builtTo, &nbody, bodyList);
        new_time = clock();
        SPRINT5(1, "--> ocsmBuild(buildTo=%d) -> status=%d (%s), builtTo=%d, nbody=%d",
                buildTo, status, ocsmGetText(status), builtTo, nbody);
        SPRINT1(1, "==> ocsmBuild CPUtime=%9.3f sec",
                (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

        if (status < 0) {
            if (recover == 0) {
                SPRINT0(0, "ERROR:: -recover not specified\a");
                goto cleanup;
            } else {
                SPRINT1(0, "WARNING:: error caused build to terminate after Branch %d", abs(builtTo));
                SPRINT0(0, "WARNING:: configuration is shown at it was at time of error\a");
            }
        }

        /* print out the Bodys */
        SPRINT0(1, "Body(s):");
        if (outLevel > 0) {
            status = ocsmPrintBodys(modl, stdout);
            if (status != SUCCESS) {
                SPRINT1(0, "ocsmPrintBodys -> status=%d", status);
            }
        }
    }

    /* start GV */
    if (batch == 0) {
        if (nbody > 0) {
            bigbox[0] = +HUGEQ;
            bigbox[1] = +HUGEQ;
            bigbox[2] = +HUGEQ;
            bigbox[3] = -HUGEQ;
            bigbox[4] = -HUGEQ;
            bigbox[5] = -HUGEQ;

            for (jbody = 0; jbody < nbody; jbody++) {
                ibody  = bodyList[jbody];
                ebody  = MODL->body[ibody].ebody;
                status = EG_getBoundingBox(ebody, box);
                if (status != SUCCESS) {
                    SPRINT1(0, "EG_getBoundingBox -> status=%d", status);
                }

                bigbox[0] = MIN(bigbox[0], box[0]);
                bigbox[1] = MIN(bigbox[1], box[1]);
                bigbox[2] = MIN(bigbox[2], box[2]);
                bigbox[3] = MAX(bigbox[3], box[3]);
                bigbox[4] = MAX(bigbox[4], box[4]);
                bigbox[5] = MAX(bigbox[5], box[5]);
            }
        } else {
            bigbox[0] = -1;
            bigbox[1] = -1;
            bigbox[2] = -1;
            bigbox[3] = +1;
            bigbox[4] = +1;
            bigbox[5] = +1;
        }

        focus[0] = (bigbox[0] + bigbox[3]) / 2;
        focus[1] = (bigbox[1] + bigbox[4]) / 2;
        focus[2] = (bigbox[2] + bigbox[5]) / 2;
        focus[3] = sqrt(SQR(bigbox[3]-bigbox[0]) + SQR(bigbox[4]-bigbox[1]) + SQR(bigbox[5]-bigbox[2]));

        gv_black[0] = gv_black[1] = gv_black[2] = 1.0;
        gv_white[0] = gv_white[1] = gv_white[2] = 0.0;

        status = gv_init("                buildCSM     ", mtflag,
                         nkeys, keys, types, lims,
                         titles, focus);
        SPRINT1(1, "--> gv_init() -> status=%d", status);
    }

    /* cleanup and exit */
    context = MODL->context;

    /* remove all Bodys and etess objects */
    for (jbody = 0; jbody < nbody; jbody++) {
        ibody = bodyList[jbody];

        if (MODL->body[ibody].etess != NULL) {
            status = EG_deleteObject(MODL->body[ibody].etess);
            SPRINT2(2, "--> EG_deleteObject(etess[%d]) -> status=%d", ibody, status);

            MODL->body[ibody].etess = NULL;
        }

        if (MODL->body[ibody].ebody != NULL) {
            status = EG_deleteObject(MODL->body[ibody].ebody);
            SPRINT2(2, "--> EG_deleteObject(ebody[%d]) => status=%d", ibody, status);

            MODL->body[ibody].ebody = NULL;
        }
    }

    /* clean up graphics objects */
    for (i = 0; i < ngrobj; i++) {
        gv_free(grobj[i], 2);
    }
    ngrobj = 0;

    /* free up the modl */
    status = ocsmFree(modl);
    SPRINT2(1, "--> ocsmFree() -> status=%d (%s)", status, ocsmGetText(status));

    /* clean up the udp storage */
    status = ocsmFree(NULL);
    SPRINT2(1, "--> ocsmFree(NULL) -> status=%d (%s)", status, ocsmGetText(status));

    /* remove the context */
    if (context != NULL) {
        status = EG_setOutLevel(context, 0);
        if (status < 0) {
            SPRINT1(0, "EG_setOutLevel -> status=%d", status);
        }

        status = EG_close(context);
        SPRINT1(1, "--> EG_close() -> status=%d", status);
    }

cleanup:
    /* if builtTo is non-positive, then report that an error was found */
    if (builtTo <= 0) {
        SPRINT0(0, "ERROR:: build not completed because an error was detected");
    } else {
        SPRINT0(1, "==> buildCSM completed successfully");
    }

    return SUCCESS;
}


/***********************************************************************/
/*                                                                     */
/*   gvupdate - called by gv to allow the changing of data             */
/*                                                                     */
/***********************************************************************/

int
gvupdate( )
{
    int       nobj;                     /* =0   no data has changed */
                                        /* >0  new number of graphic objects */

    int        i, ibody, jbody;
    static int init = 0;
    char       *family_name;

    modl_T    *MODL = (modl_T*)modl;

//    ROUTINE(gvupdate);

    /* --------------------------------------------------------------- */

//$$$    SPRINT0(2, "enter gvupdate()");

    /* one-time initialization */
    if (init == 0) {
        GraphicGCSetFB(gv_wAux.GCs,  gv_white, gv_black);
        GraphicBGColor(gv_wAux.wid,  gv_black);

        GraphicGCSetFB(gv_wDial.GCs, gv_white, gv_black);
        GraphicBGColor(gv_wDial.wid, gv_black);

        init++;
    }

    /* simply return if no new data */
    if (new_data == 0) {
        return 0;
    }

    /* remove any previous families */
    for (i = gv_numfamily()-1; i >= 0; i--) {
        (void) gv_returnfamily(0, &family_name);
        (void) gv_freefamily(family_name);
    }

    /* remove any previous graphic objects */
    for (i = 0; i < ngrobj; i++) {
        gv_free(grobj[i], 2);
    }
    ngrobj = 0;

    /* get the new number of Edges and Faces */
    nobj = 3;                 // axes
    for (jbody = 0; jbody < nbody; jbody++) {
        ibody  = bodyList[jbody];

        nobj += (MODL->body[ibody].nedge + MODL->body[ibody].nface);
    }
    if (tuftLen > 0) {
        nobj += nbody;
    }

    /* return the number of graphic objects */
    new_data = 0;

//cleanup:
    return nobj;
}


/***********************************************************************/
/*                                                                     */
/*   gvdata - called by gv to (re)set the graphic objects              */
/*                                                                     */
/***********************************************************************/

void
gvdata(int       ngraphics,             /* (in)  number of objects to fill */
       GvGraphic *graphic[])            /* (out) graphic objects to fill */
{
    GvColor   color;
    GvObject  *object;
    int       i, j, npnt, ntri, ntuft, ibrch, iattr, showBody;
    int       utype, iedge, iface, mask, ibody, jbody, nedge, nface, attr;
    CINT      *tris, *tric, *ptype, *pindx;
    CDOUBLE   *xyz, *uv;
    char      title[16], bodyName[16];

    int       ipnt, inode;
    double    *dxyz;
    ego       etess;

    modl_T    *MODL = (modl_T*)modl;

//    ROUTINE(gvdata);

    /* --------------------------------------------------------------- */

//$$$    SPRINT1(2, "enter gvdata(ngraphics=%d)",
//$$$           ngraphics);

    /* create the graphic objects */
    grobj  = graphic;
    ngrobj = ngraphics;
    i      = 0;

    /* if the family does not exist, create it */
    if (gv_getfamily(  "Axes", 1, &attr) == -1) {
        gv_allocfamily("Axes");
    }

    /* x-axis */
    mask = GV_FOREGROUND | GV_FORWARD;

    color.red   = 1.0;
    color.green = 0.0;
    color.blue  = 0.0;

    sprintf(title, "X axis");
    utype = 999;
    graphic[i] = gv_alloc(GV_NONINDEXED, GV_DISJOINTLINES, mask, color, title, utype, 0);
    if (graphic[i] != NULL) {
        graphic[i]->number    = 1;
        graphic[i]->lineWidth = 3;
        graphic[i]->ddata = (double*)  malloc(18*sizeof(double));

        graphic[i]->ddata[ 0] = 0.0;
        graphic[i]->ddata[ 1] = 0.0;
        graphic[i]->ddata[ 2] = 0.0;
        graphic[i]->ddata[ 3] = 1.0;
        graphic[i]->ddata[ 4] = 0.0;
        graphic[i]->ddata[ 5] = 0.0;

        graphic[i]->ddata[ 6] = 1.1;
        graphic[i]->ddata[ 7] = 0.1;
        graphic[i]->ddata[ 8] = 0.0;
        graphic[i]->ddata[ 9] = 1.3;
        graphic[i]->ddata[10] = -.1;
        graphic[i]->ddata[11] = 0.0;

        graphic[i]->ddata[12] = 1.3;
        graphic[i]->ddata[13] = 0.1;
        graphic[i]->ddata[14] = 0.0;
        graphic[i]->ddata[15] = 1.1;
        graphic[i]->ddata[16] = -.1;
        graphic[i]->ddata[17] = 0.0;
        graphic[i]->object->length = 3;
        graphic[i]->object->type.plines.len = (int*)   malloc(3*sizeof(int));
        graphic[i]->object->type.plines.len[0] = 2;
        graphic[i]->object->type.plines.len[1] = 2;
        graphic[i]->object->type.plines.len[2] = 2;

        gv_adopt("Axes", graphic[i]);

        i++;
    }

    /* y-axis */
    mask = GV_FOREGROUND | GV_FORWARD;

    color.red   = 0.0;
    color.green = 1.0;
    color.blue  = 0.0;

    sprintf(title, "Y axis");
    utype = 999;
    graphic[i] = gv_alloc(GV_NONINDEXED, GV_DISJOINTLINES, mask, color, title, utype, 0);
    if (graphic[i] != NULL) {
        graphic[i]->number    = 1;
        graphic[i]->lineWidth = 3;
        graphic[i]->ddata = (double*)  malloc(24*sizeof(double));

        graphic[i]->ddata[ 0] = 0.0;
        graphic[i]->ddata[ 1] = 0.0;
        graphic[i]->ddata[ 2] = 0.0;
        graphic[i]->ddata[ 3] = 0.0;
        graphic[i]->ddata[ 4] = 1.0;
        graphic[i]->ddata[ 5] = 0.0;

        graphic[i]->ddata[ 6] = 0.0;
        graphic[i]->ddata[ 7] = 1.1;
        graphic[i]->ddata[ 8] = 0.1;
        graphic[i]->ddata[ 9] = 0.0;
        graphic[i]->ddata[10] = 1.2;
        graphic[i]->ddata[11] = 0.0;

        graphic[i]->ddata[12] = 0.0;
        graphic[i]->ddata[13] = 1.3;
        graphic[i]->ddata[14] = 0.1;
        graphic[i]->ddata[15] = 0.0;
        graphic[i]->ddata[16] = 1.2;
        graphic[i]->ddata[17] = 0.0;

        graphic[i]->ddata[18] = 0.0;
        graphic[i]->ddata[19] = 1.2;
        graphic[i]->ddata[20] = 0.0;
        graphic[i]->ddata[21] = 0.0;
        graphic[i]->ddata[22] = 1.2;
        graphic[i]->ddata[23] = -.1;
        graphic[i]->object->length = 4;
        graphic[i]->object->type.plines.len = (int*)   malloc(4*sizeof(int));
        graphic[i]->object->type.plines.len[0] = 2;
        graphic[i]->object->type.plines.len[1] = 2;
        graphic[i]->object->type.plines.len[2] = 2;
        graphic[i]->object->type.plines.len[3] = 2;

        gv_adopt("Axes", graphic[i]);

        i++;
    }

    /* z-axis */
    mask = GV_FOREGROUND | GV_FORWARD;

    color.red   = 0.0;
    color.green = 0.0;
    color.blue  = 1.0;

    sprintf(title, "Z axis");
    utype = 999;
    graphic[i] = gv_alloc(GV_NONINDEXED, GV_DISJOINTLINES, mask, color, title, utype, 0);
    if (graphic[i] != NULL) {
        graphic[i]->number    = 1;
        graphic[i]->lineWidth = 3;
        graphic[i]->ddata = (double*)  malloc(24*sizeof(double));

        graphic[i]->ddata[ 0] = 0.0;
        graphic[i]->ddata[ 1] = 0.0;
        graphic[i]->ddata[ 2] = 0.0;
        graphic[i]->ddata[ 3] = 0.0;
        graphic[i]->ddata[ 4] = 0.0;
        graphic[i]->ddata[ 5] = 1.0;

        graphic[i]->ddata[ 6] = 0.1;
        graphic[i]->ddata[ 7] = 0.0;
        graphic[i]->ddata[ 8] = 1.1;
        graphic[i]->ddata[ 9] = 0.1;
        graphic[i]->ddata[10] = 0.0;
        graphic[i]->ddata[11] = 1.3;

        graphic[i]->ddata[12] = 0.1;
        graphic[i]->ddata[13] = 0.0;
        graphic[i]->ddata[14] = 1.3;
        graphic[i]->ddata[15] = -.1;
        graphic[i]->ddata[16] = 0.0;
        graphic[i]->ddata[17] = 1.1;

        graphic[i]->ddata[18] = -.1;
        graphic[i]->ddata[19] = 0.0;
        graphic[i]->ddata[20] = 1.1;
        graphic[i]->ddata[21] = -.1;
        graphic[i]->ddata[22] = 0.0;
        graphic[i]->ddata[23] = 1.3;
        graphic[i]->object->length = 4;
        graphic[i]->object->type.plines.len = (int*)   malloc(4*sizeof(int));
        graphic[i]->object->type.plines.len[0] = 2;
        graphic[i]->object->type.plines.len[1] = 2;
        graphic[i]->object->type.plines.len[2] = 2;
        graphic[i]->object->type.plines.len[3] = 2;

        gv_adopt("Axes", graphic[i]);

        i++;
    }

    /* Bodys */
    for (jbody = 0; jbody < nbody; jbody++) {
        ibody = bodyList[jbody];
        nedge = MODL->body[ibody].nedge;
        nface = MODL->body[ibody].nface;

        showBody = 1;
        ibrch    = MODL->body[ibody].ibrch;
        for (iattr = 0; iattr < MODL->brch[ibrch].nattr; iattr++) {
            if (strcmp(MODL->brch[ibrch].attr[iattr].name, "ShowBody") == 0) {
                if (MODL->brch[ibrch].attr[iattr].defn[0] == 0) {
                    showBody = 0;
                }
            }
        }

        /* if the family does not exist, create it */
        sprintf(bodyName, "Body %d", ibody);

        if (gv_getfamily( bodyName, 1, &attr) == -1){
            gv_allocfamily(bodyName);
        }

        /* create a graphic object for each Edge (in ivol) */
        for (iedge = 1; iedge <= nedge; iedge++) {

            /* get the Edge info */
            etess = MODL->body[ibody].etess;
            (void) EG_getTessEdge(etess, iedge, &npnt, &xyz, &uv);

            /* set up the new graphic object */
            mask = MODL->body[ibody].edge[iedge].gratt.render;
            if (showBody == 0 && (mask & GV_FOREGROUND)) {
                mask -= GV_FOREGROUND;
            }

            color.red   = RED(  MODL->body[ibody].edge[iedge].gratt.color);
            color.green = GREEN(MODL->body[ibody].edge[iedge].gratt.color);
            color.blue  = BLUE( MODL->body[ibody].edge[iedge].gratt.color);

            sprintf(title, "Edge %d:%d", ibody, iedge);
            utype = 1 + 10 * ibody;
            graphic[i] = gv_alloc(GV_NONINDEXED, GV_POLYLINES, mask, color, title, utype, iedge);

            if (graphic[i] != NULL) {
                graphic[i]->number     = 1;
                graphic[i]->lineWidth  = MODL->body[ibody].edge[iedge].gratt.lwidth;
                graphic[i]->pointSize  = 3;
                graphic[i]->mesh.red   = 0;
                graphic[i]->mesh.green = 0;
                graphic[i]->mesh.blue  = 0;

                /* load the data */
                graphic[i]->ddata = (double*)  malloc(3*npnt*sizeof(double));
                for (j = 0; j < npnt; j++) {
                    graphic[i]->ddata[3*j  ] = xyz[3*j  ];
                    graphic[i]->ddata[3*j+1] = xyz[3*j+1];
                    graphic[i]->ddata[3*j+2] = xyz[3*j+2];
                }

                object = graphic[i]->object;
                object->length = 1;
                object->type.plines.len = (int *)   malloc(sizeof(int));
                object->type.plines.len[0] = npnt;

                gv_adopt(bodyName, graphic[i]);
            }
            i++;
        }

        /* create a graphic object for each Face (in ivol) */
        for (iface = 1; iface <= nface; iface++) {

            /* get the Face info */
            etess = MODL->body[ibody].etess;
            if (etess == NULL) continue;
            (void) EG_getTessFace(etess, iface, &npnt, &xyz, &uv, &ptype, &pindx,
                                  &ntri, &tris, &tric);

            /* set up new grafic object */
            mask = MODL->body[ibody].face[iface].gratt.render;
            if (sclrType >= 0) {
                mask |= GV_SCALAR;
            }
            if (showBody == 0 && (mask & GV_FOREGROUND)) {
                mask -= GV_FOREGROUND;
            }
            color.red   = RED(  MODL->body[ibody].face[iface].gratt.color);
            color.green = GREEN(MODL->body[ibody].face[iface].gratt.color);
            color.blue  = BLUE( MODL->body[ibody].face[iface].gratt.color);

            sprintf(title, "Face %d:%d", ibody, iface);
            utype = 2 + 10 * ibody;
            graphic[i] = gv_alloc(GV_INDEXED, GV_DISJOINTTRIANGLES, mask, color, title, utype, iface);

            if (graphic[i] != NULL) {
                graphic[i]->back.red   = RED(  MODL->body[ibody].face[iface].gratt.bcolor);
                graphic[i]->back.green = GREEN(MODL->body[ibody].face[iface].gratt.bcolor);
                graphic[i]->back.blue  = BLUE( MODL->body[ibody].face[iface].gratt.bcolor);

                graphic[i]->mesh.red   = RED(  MODL->body[ibody].face[iface].gratt.mcolor);
                graphic[i]->mesh.green = GREEN(MODL->body[ibody].face[iface].gratt.mcolor);
                graphic[i]->mesh.blue  = BLUE( MODL->body[ibody].face[iface].gratt.mcolor);

                graphic[i]->number    = 1;
                graphic[i]->lineWidth = MODL->body[ibody].face[iface].gratt.lwidth;

                /* load the data */
                graphic[i]->ddata = (double*)  malloc(3*npnt*sizeof(double));
                for (j = 0; j < npnt; j++) {
                    graphic[i]->ddata[3*j  ] = xyz[3*j  ];
                    graphic[i]->ddata[3*j+1] = xyz[3*j+1];
                    graphic[i]->ddata[3*j+2] = xyz[3*j+2];
                }


                object = graphic[i]->object;
                object->length = ntri;
                object->type.distris.index = (int *)   malloc(3*ntri*sizeof(int));

                for (j = 0; j < 3*ntri; j++) {
                    object->type.distris.index[j] = tris[j]-1;
                }

                gv_adopt(bodyName, graphic[i]);
            }
            i++;
        }

        /* create a graphic object for the tufts */
        if (tuftLen > 0) {

            /* count the number of point in all Nodes and Edges */
            ntuft = MODL->body[ibody].nnode;

            for (iedge = 1; iedge <= nedge; iedge++) {

                /* get the Edge info */
                etess = MODL->body[ibody].etess;
                (void) EG_getTessEdge(etess, iedge, &npnt, &xyz, &uv);

                ntuft += npnt;
            }

            /* set up new grafic object */
            mask = GV_SCALAR;
            color.red   = 0;
            color.green = 0;
            color.blue  = 0;

            sprintf(title, "Tufts %d", ibody);
            utype = 3 + 10 * ibody;
            graphic[i] = gv_alloc(GV_NONINDEXED, GV_DISJOINTLINES, mask, color, title, utype, ibody);

            if (graphic[i] != NULL) {
                graphic[i]->number    = 1;
                graphic[i]->lineWidth = 3;

                /* load the data */
                graphic[i]->ddata = (double*) malloc(6*ntuft*sizeof(double));
                j = 0;
                for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                    dxyz = (double*) malloc(3*sizeof(double));
                    (void) ocsmGetVel(MODL, ibody, OCSM_NODE, inode, 1, NULL, dxyz);

                    graphic[i]->ddata[6*j  ] = MODL->body[ibody].node[inode].x;
                    graphic[i]->ddata[6*j+1] = MODL->body[ibody].node[inode].y;
                    graphic[i]->ddata[6*j+2] = MODL->body[ibody].node[inode].z;
                    graphic[i]->ddata[6*j+3] = MODL->body[ibody].node[inode].x + tuftLen * dxyz[0];
                    graphic[i]->ddata[6*j+4] = MODL->body[ibody].node[inode].y + tuftLen * dxyz[1];
                    graphic[i]->ddata[6*j+5] = MODL->body[ibody].node[inode].z + tuftLen * dxyz[2];
                    j++;

                    free(dxyz);
                }
                for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                    etess = MODL->body[ibody].etess;
                    (void) EG_getTessEdge(etess, iedge, &npnt, &xyz, &uv);

                    dxyz = (double*) malloc(3*npnt*sizeof(double));
                    (void) ocsmGetVel(MODL, ibody, OCSM_EDGE, iedge, npnt, NULL, dxyz);

                    for (ipnt = 0; ipnt < npnt; ipnt++) {
                        graphic[i]->ddata[6*j  ] = xyz[3*ipnt  ];
                        graphic[i]->ddata[6*j+1] = xyz[3*ipnt+1];
                        graphic[i]->ddata[6*j+2] = xyz[3*ipnt+2];
                        graphic[i]->ddata[6*j+3] = xyz[3*ipnt  ] + tuftLen * dxyz[3*ipnt  ];
                        graphic[i]->ddata[6*j+4] = xyz[3*ipnt+1] + tuftLen * dxyz[3*ipnt+1];
                        graphic[i]->ddata[6*j+5] = xyz[3*ipnt+2] + tuftLen * dxyz[3*ipnt+2];
                        j++;
                    }

                    free(dxyz);
                }

                object = graphic[i]->object;
                object->length = ntuft;
            }
            i++;
        }
    }

//cleanup:
    return;
}


/***********************************************************************/
/*                                                                     */
/*   gvscalar - called by gv for color rendering of graphic objects    */
/*                                                                     */
/***********************************************************************/

int
gvscalar(int       key,                 /* (in)  scalar index (from gv_init) */
         GvGraphic *graphic,            /* (in)  GvGraphic structure for scalar fill */
/*@unused@*/int    len,                 /* (in)  length of scalar to be filled */
         float     *scalar)             /* (out) scalar to be filled */
{
    int       answer;                   /* return value */

    int       utype, iface, i, ibody;
    int       ntri, npnt;
    CINT      *tris, *tric, *ptype, *pindx;
    double    umin, umax, vmin, vmax, *vel;
    CDOUBLE   *xyz, *uv;

    ego       etess;

    modl_T    *MODL = (modl_T*)modl;

//    ROUTINE(gvscalar);

    /* --------------------------------------------------------------- */

//$$$    SPRINT3(2, "enter gvscalar(key=%d, len=%d, scalar=%f",
//$$$           key, len, *scalar);

    utype  = graphic->utype;
    iface  = graphic->uindex;
    ibody  = utype / 10;

    etess = MODL->body[ibody].etess;

    /* Face */
    if (utype%10 == 2) {
         (void) EG_getTessFace(etess, iface, &npnt, &xyz, &uv, &ptype, &pindx,
                               &ntri, &tris, &tric);

         if (key == 0) {
             umin = uv[0];
             umax = uv[0];
             for (i = 0; i < npnt; i++) {
                 umin = MIN(umin, uv[2*i]);
                 umax = MAX(umax, uv[2*i]);
             }

             for (i = 0; i < npnt; i++) {
                 scalar[i] = (uv[2*i] - umin) / (umax - umin);
             }
         } else if (key == 1) {
             vmin = uv[1];
             vmax = uv[1];
             for (i = 0; i < npnt; i++) {
                 vmin = MIN(vmin, uv[2*i+1]);
                 vmax = MAX(vmax, uv[2*i+1]);
             }

             for (i = 0; i < npnt; i++) {
                 scalar[i] = (uv[2*i+1] - vmin) / (vmax - vmin);
             }
         } else if (key == 2) {
             vel = (double*) malloc(3*npnt*sizeof(double));

             (void) ocsmGetVel(MODL, ibody, OCSM_FACE, iface, npnt, NULL, vel);

             for (i = 0; i < npnt; i++) {
                 scalar[i] = (float)(sqrt(  vel[3*i  ]*vel[3*i  ]
                                          + vel[3*i+1]*vel[3*i+1]
                                          + vel[3*i+2]*vel[3*i+2]));
             }

             free(vel);
         } else {
             for (i = 0; i < npnt; i++) {
                 scalar[i] = 0.;
             }
         }
    }

//cleanup:
    answer = 1;
    return answer;
}


/***********************************************************************/
/*                                                                     */
/*   gvevent - called by gv calls to process callbacks                 */
/*                                                                     */
/***********************************************************************/

void
gvevent(int       *win,                 /* (in)  window of event */
        int       *type,                /* (in)  type of event */
/*@unused@*/int   *xpix,                /* (in)  x-pixel location of event */
/*@unused@*/int   *ypix,                /* (in)  y-pixel location of event */
        int       *state)               /* (in)  aditional event info */
{
    FILE      *fp=NULL;

    int       uindex, utype, inode, iedge, iface, ibrch, ipmtr, irow, icol, nrow, ncol;
    int       itype, status, idum, ibody, ipaste, dum1, dum2;
    int       i, j, iclass, iactv, ichld, ileft, irite, narg, iarg, nattr, iattr;
    int       nlist, ngood, jbody;
    CINT      *tempIlist;
    double    value, dot, dum, uv[2], xyz[3], vel[1000];
    CDOUBLE   *tempRlist;
    char      jnlName[255], fileName[255], tempName[255];
    char      pmtrName[MAX_NAME_LEN], defn[255], brchName[255];
    char      aName[255], aValue[255];
    CCHAR     *tempClist;

    static int utype_save = 0, uindex_save = 0;

    int     ntemp;
    double  massprops[18];
    CCHAR   *attrName;

    modl_T    *MODL = (modl_T*)modl;

    ROUTINE(gvevent);

    /* --------------------------------------------------------------- */

//$$$    SPRINT5(2, "enter gvevent(*win=%d, *type=%d, *state=%d (%c), script=%d)",
//$$$            *win, *type, *state, *state, (int)script);

    /* repeat as long as we are reading a script (or once if
       not reading a script) */
    do {

        /* get the next script line if we are reading a script (and insert
           a '$' if we have detected an EOF) */

        if (script != NULL) {
            if (fscanf(script, "%1s", (char*)state) != 1) {
                *state = '$';
            }
            *win  = ThreeD;
            *type = KeyPress;
        }

        if ((*win == ThreeD) && (*type == KeyPress)) {
            if (*state == '\0') {

                /* these calls should never be made */
                idum = getInt("Dummy call to use getInt");
                dum  = getDbl("Dummy call to use getDbl");
                SPRINT2(0, "idum=%d   dum=%f", idum, dum);

            /* 'a' - add Parameter */
            } else if (*state == 'a') {
                SPRINT0(0, "--> Option 'a' chosen (add Parameter)");

                getStr("Enter Parameter name: ", pmtrName);
                nrow = getInt("Enter number of rows: ");
                ncol = getInt("Enter number of cols: ");

                status = ocsmNewPmtr(modl, pmtrName, OCSM_DESPMTR, nrow, ncol);
                SPRINT5(0, "--> ocsmNewPmtr(name=%s, nrow=%d, ncol=%d) -> status=%d (%s)",
                        pmtrName, nrow, ncol, status, ocsmGetText(status));

                status = ocsmInfo(modl, &dum1, &ipmtr, &dum2);
                if (status != SUCCESS) {
                    SPRINT1(0, "ocsmInfo -> status=%d", status);
                }

                for (icol = 1; icol <= ncol; icol++) {
                    for (irow = 1; irow <= nrow; irow++) {
                        SPRINT1x(0, "Enter value for %s", pmtrName);
                        SPRINT1x(0, "[%d,",               irow    );
                        SPRINT1x(0, "%d]",                icol    );
                        getStr(": ", defn);

                        status = ocsmSetValu(modl, ipmtr, irow, icol, defn);
                        SPRINT5(0, "--> ocsmSetValu(irow=%d, icol=%d, defn=%s) -> status=%d (%s)",
                                irow, icol, defn, status, ocsmGetText(status));
                    }
                }

            /* 'A' - add Branch */
            } else if (*state == 'A') {
                char str1[257], str2[257], str3[257], str4[257], str5[257];
                char str6[257], str7[257], str8[257], str9[257];
                SPRINT0(0, "--> Option 'A' chosen (add Branch)");

                SPRINT0(0, "1 box        11 extrude    31 intersect  51 set   ");
                SPRINT0(0, "2 sphere     12 loft       32 subtract   52 macbeg");
                SPRINT0(0, "3 cone       13 revolve    33 union      53 macend");
                SPRINT0(0, "4 cylinder                               54 recall");
                SPRINT0(0, "5 torus      21 fillet     41 translate  55 patbeg");
                SPRINT0(0, "6 import     22 chamfer    42 rotatex    56 patend");
                SPRINT0(0, "7 udprim                   43 rotatey    57 mark  ");
                SPRINT0(0, "                           44 rotatez    58 dump  ");
                SPRINT0(0, "                           45 scale               ");
                itype = getInt("Enter type to add: ");

                strcpy(str1, "$");
                strcpy(str2, "$");
                strcpy(str3, "$");
                strcpy(str4, "$");
                strcpy(str5, "$");
                strcpy(str6, "$");
                strcpy(str7, "$");
                strcpy(str8, "$");
                strcpy(str9, "$");

                if        (itype ==  1) {
                    itype = OCSM_BOX;
                    getStr("Enter xbase : ",  str1   );
                    getStr("Enter ybase : ",  str2   );
                    getStr("Enter zbase : ",  str3   );
                    getStr("Enter dx    : ",  str4   );
                    getStr("Enter dy    : ",  str5   );
                    getStr("Enter dz    : ",  str6   );
                } else if (itype ==  2) {
                    itype = OCSM_SPHERE;
                    getStr("Enter xcent : ",  str1   );
                    getStr("Enter ycent : ",  str2   );
                    getStr("Enter zcent : ",  str3   );
                    getStr("Enter radius: ",  str4   );
                } else if (itype ==  3) {
                    itype = OCSM_CONE;
                    getStr("Enter xvrtx : ",  str1   );
                    getStr("Enter yvrtx : ",  str2   );
                    getStr("Enter zvrtx : ",  str3   );
                    getStr("Enter xbase : ",  str4   );
                    getStr("Enter ybase : ",  str5   );
                    getStr("Enter zbase : ",  str6   );
                    getStr("Enter radius: ",  str7   );
                } else if (itype ==  4) {
                    itype = OCSM_CYLINDER;
                    getStr("Enter xbeg  : ",  str1   );
                    getStr("Enter ybeg  : ",  str2   );
                    getStr("Enter zbeg  : ",  str3   );
                    getStr("Enter xend  : ",  str4   );
                    getStr("Enter yend  : ",  str5   );
                    getStr("Enter zend  : ",  str6   );
                    getStr("Enter radius: ",  str7   );
                } else if (itype ==  5) {
                    itype = OCSM_TORUS;
                    getStr("Enter xcent : ",  str1   );
                    getStr("Enter ycent : ",  str2   );
                    getStr("Enter zcent : ",  str3   );
                    getStr("Enter dxaxis: ",  str4   );
                    getStr("Enter dyaxis: ",  str5   );
                    getStr("Enter dzaxis: ",  str6   );
                    getStr("Enter majrad: ",  str7   );
                    getStr("Enter minrad: ",  str8   );
                } else if (itype ==  6) {
                    itype = OCSM_IMPORT;
                    getStr("Enter filNam: ", &str1[1]);
                } else if (itype ==  7) {
                    itype = OCSM_UDPRIM;
                    getStr("Enter ptype : ", &str1[1]);
                    getStr("Enter name1 : ", &str2[1]);
                    getStr("Enter value1: ", &str3[1]);
                    getStr("Enter name2 : ", &str4[1]);
                    getStr("Enter value2: ", &str5[1]);
                    getStr("Enter name3 : ", &str6[1]);
                    getStr("Enter value3: ", &str7[1]);
                    getStr("Enter name4 : ", &str8[1]);
                    getStr("Enter value4: ", &str9[1]);
                } else if (itype == 11) {
                    itype = OCSM_EXTRUDE;
                    getStr("Enter dx    : ",  str1   );
                    getStr("Enter dy    : ",  str2   );
                    getStr("Enter dz    : ",  str3   );
                } else if (itype == 12) {
                    itype = OCSM_LOFT;
                    getStr("Enter smooth: ",  str1   );
                } else if (itype == 13) {
                    itype = OCSM_REVOLVE;
                    getStr("Enter xorig : ",  str1   );
                    getStr("Enter yorig : ",  str2   );
                    getStr("Enter zorig : ",  str3   );
                    getStr("Enter dxaxis: ",  str4   );
                    getStr("Enter dyaxis: ",  str5   );
                    getStr("Enter dzaxis: ",  str6   );
                    getStr("Enter angDeg: ",  str7   );
                } else if (itype == 21) {
                    itype = OCSM_FILLET;
                    getStr("Enter radius: ",  str1   );
                    getStr("Enter iford1: ",  str2   );
                    getStr("Enter iford2: ",  str3   );
                } else if (itype == 22) {
                    itype = OCSM_CHAMFER;
                    getStr("Enter radius: ",  str1   );
                    getStr("Enter iford1: ",  str2   );
                    getStr("Enter iford2: ",  str3   );
                } else if (itype == 31) {
                    itype = OCSM_INTERSECT;
                    getStr("Enter order : ", &str1[1]);
                    getStr("Enter index : ",  str2   );
                    strcpy(str3, "0");
                } else if (itype == 32) {
                    itype = OCSM_SUBTRACT;
                    getStr("Enter order : ", &str1[1]);
                    getStr("Enter index : ",  str2   );
                    strcpy(str3, "0");
                } else if (itype == 33) {
                    itype = OCSM_UNION;
                    getStr("Enter tomark: ",  str1   );
                    getStr("Enter trmLst: ", &str2[1]);
                    strcpy(str3, "0");
                } else if (itype == 41) {
                    itype = OCSM_TRANSLATE;
                    getStr("Enter dx    : ",  str1   );
                    getStr("Enter dy    : ",  str2   );
                    getStr("Enter dz    : ",  str3   );
                } else if (itype == 42) {
                    itype = OCSM_ROTATEX;
                    getStr("Enter angDeg: ",  str1   );
                    getStr("Enter yaxis : ",  str2   );
                    getStr("Enter zaxis : ",  str3   );
                } else if (itype == 43) {
                    itype = OCSM_ROTATEY;
                    getStr("Enter angDeg: ",  str1   );
                    getStr("Enter zaxis : ",  str2   );
                    getStr("Enter xaxis : ",  str3   );
                } else if (itype == 44) {
                    itype = OCSM_ROTATEZ;
                    getStr("Enter angDeg: ",  str1   );
                    getStr("Enter xaxis : ",  str2   );
                    getStr("Enter yaxish: ",  str3   );
                } else if (itype == 45) {
                    itype = OCSM_SCALE;
                    getStr("Enter fact  : ",  str1   );
                } else if (itype == 51) {
                    itype = OCSM_SET;
                    getStr("Enter pname : ", &str1[1]);
                    getStr("Enter defn  : ",  str2   );
                } else if (itype == 52) {
                    itype = OCSM_MACBEG;
                    getStr("Enter istore: ",  str1   );
                } else if (itype == 53) {
                    itype = OCSM_MACEND;
                } else if (itype == 54) {
                    itype = OCSM_RECALL;
                    getStr("Enter istore: ",  str1   );
                } else if (itype == 55) {
                    itype = OCSM_PATBEG;
                    getStr("Enter pname : ", &str1[1]);
                    getStr("Enter ncopy : ",  str2   );
                } else if (itype == 56) {
                    itype = OCSM_PATEND;
                } else if (itype == 57) {
                    itype = OCSM_MARK;
                } else if (itype == 58) {
                    itype = OCSM_DUMP;
                    getStr("Enter filNam: ", &str1[1]);
                    getStr("Enter remove: ",  str2   );
                } else {
                    SPRINT1(0, "Illegal type (%d)", itype);
                    goto next_option;
                }

                status = ocsmNewBrch(MODL, MODL->nbrch, itype, "<ESP>", -1,
                                     str1, str2, str3, str4, str5,
                                     str6, str7, str8, str9);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsnNewBrch(ibrch=%d) -> status=%d (%s)",
                            MODL->nbrch, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT1(0, "Branch %d has been added", MODL->nbrch);
                SPRINT0(0, "Use 'B' to rebuild");

            /* 'b' - undefined option */
            } else if (*state == 'b') {
                SPRINT0(0, "--> Option 'b' (undefined)");

                for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                    if (MODL->body[ibody].onstack) {
                        int periodic;
                        double range[4];

                        for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                            EG_getRange(MODL->body[ibody].edge[iedge].eedge, range, &periodic);

                            SPRINT4(0, "ibody=%2d,  iedge=%3d,  range=%10.4f %10.4f",
                                    ibody, iedge, range[0], range[1]);
                        }

                        for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                            EG_getRange(MODL->body[ibody].face[iface].eface, range, &periodic);

                            SPRINT6(0, "ibody=%2d,  iface=%3d,  range=%10.4f %10.4f %10.4f %10.4f",
                                    ibody, iface, range[0], range[1], range[2], range[3]);
                        }
                    }
                }

            /* 'B' - build to Branch */
            } else if (*state == 'B') {
                int     buildTo;
                clock_t old_time, new_time;
                SPRINT0(0, "--> Option 'B' chosen (build to Branch)");

                if (numarg > 0) {
                    buildTo = numarg;
                    numarg  = 0;
                } else {
                    buildTo = 0;
                }

                old_time = clock();
                status   = ocsmCheck(modl);
                new_time = clock();
                SPRINT2(0, "--> ocsmCheck() -> status=%d (%s)",
                        status, ocsmGetText(status));
                SPRINT1(0, "==> ocsmCheck CPUtime=%10.3f sec",
                        (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

                if (status < SUCCESS) goto next_option;

                nbody    = MAX_BODYS;
                old_time = clock();
                status   = ocsmBuild(modl, buildTo, &builtTo, &nbody, bodyList);
                new_time = clock();
                SPRINT5(0, "--> ocsmBuild(buildTo=%d) -> status=%d (%s), builtTo=%d, nbody=%d",
                        buildTo, status, ocsmGetText(status), builtTo, nbody);
                SPRINT1(0, "==> ocsmBuild CPUtime=%10.3f sec",
                        (double)(new_time-old_time) / (double)(CLOCKS_PER_SEC));

                if (status < SUCCESS) goto next_option;

                if (status >= 0) {
                    new_data = 1;
                }

            /* 'c' - undefined option */
            } else if (*state == 'c') {
                int     npnt, ntri, ipnt;
                int     nsketch=0, isketch[100], kbody, isketchs, isketchn;
                CINT    *ptype, *pindx, *tris, *tric;
                double  uvruled[2], xyzruled[3], err, errmax;
                CDOUBLE *xyz1, *uv1;

                SPRINT0(0, "--> Option 'c' (test RULEd surface)");

                for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                    if (MODL->body[ibody].onstack != 1) continue;

                    if (MODL->body[ibody].brtype == OCSM_RULE) {
                        for (kbody = MODL->body[ibody].irite; kbody >= MODL->body[ibody].ileft; kbody--) {
                            if (MODL->body[kbody].ichld == ibody) {
                                isketch[nsketch++] = kbody;
                            }
                        }

                        for (iface = 1; iface <= MODL->body[ibody].nface-2; iface++) {
                            isketchs = isketch[    (iface - 1) % (nsketch-1)];
                            isketchn = isketch[1 + (iface - 1) % (nsketch-1)];
                            iedge    =         1 + (iface - 1) / (nsketch-1);

                            (void) EG_getTessFace(MODL->body[ibody].etess, iface, &npnt, &xyz1, &uv1, &ptype, &pindx,
                                                  &ntri, &tris, &tric);

                            errmax = 0;
                            for (ipnt = 0; ipnt < npnt; ipnt++) {
                                uvruled[0] = uv1[2*ipnt  ];
                                uvruled[1] = uv1[2*ipnt+1];

                                status = evalRuled(MODL, ibody, iface, isketchs, isketchn, iedge, uvruled, xyzruled);
                                if (status != SUCCESS) {
                                    SPRINT1(0, "evalRuled -> status=%d", status);
                                }

                                err = sqrt((xyz1[3*ipnt  ] - xyzruled[0]) * (xyz1[3*ipnt  ] - xyzruled[0])
                                          +(xyz1[3*ipnt+1] - xyzruled[1]) * (xyz1[3*ipnt+1] - xyzruled[1])
                                          +(xyz1[3*ipnt+2] - xyzruled[2]) * (xyz1[3*ipnt+2] - xyzruled[2]));
                                if (err > errmax) errmax = err;
                            }
                        }
                    }
                }

            /* 'C' - write STL file */
            } else if (*state == 'C') {
                int     npnt, ntri, ntri_tot, ip0, ip1, ip2, itri, icolr;
                CINT    *pindx, *ptype, *tris, *tric;
                float   norm[3], vert[3], area;
                CDOUBLE *xyz2, *uv2;
                char    stl_filename[255], header[255];
                FILE    *stl_fp;
                ego     etess;

                #define UINT32 unsigned int
                #define UINT16 unsigned short int
                #define REAL32 float

                SPRINT0(0, "--> Option 'C' (write STL file)");

                getStr("Enter .stl filename: ", stl_filename);
                if (strstr(stl_filename, ".stl") == NULL) {
                    strcat(stl_filename, ".csm");
                }

                /* open the output file */
                stl_fp = fopen(stl_filename, "wb");
                if (stl_fp == NULL) {
                    SPRINT1(0, "ERROR:: problem opening file \"%s\"", stl_filename);
                } else {

                    /* create header */
                    sprintf(header, "written by buildCSM, ncolr=%d", 0);
                    (void) fwrite(header, sizeof(char), 80, stl_fp);

                    /* number of Triangles in Bodys on the stack */
                    ntri_tot = 0;
                    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                        if (MODL->body[ibody].onstack == 1) {
                            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                                etess = MODL->body[ibody].etess;
                                if (etess == NULL) continue;
                                (void) EG_getTessFace(etess, iface, &npnt, &xyz2, &uv2, &ptype, &pindx,
                                                      &ntri, &tris, &tric);
                                ntri_tot += ntri;
                            }
                        }
                    }
                    (void) fwrite(&ntri_tot, sizeof(UINT32), 1, stl_fp);

                    SPRINT1(0, "--> ntri_tot=%d", ntri_tot);

                    /* write the Triangles */
                    for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                        if (MODL->body[ibody].onstack == 1) {
                            for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                                etess = MODL->body[ibody].etess;
                                if (etess == NULL) continue;
                                (void) EG_getTessFace(etess, iface, &npnt, &xyz2, &uv2, &ptype, &pindx,
                                                      &ntri, &tris, &tric);

                                for (itri = 0; itri < ntri; itri++) {
                                    ip0 = tris[3*itri  ] - 1;
                                    ip1 = tris[3*itri+1] - 1;
                                    ip2 = tris[3*itri+2] - 1;

                                    norm[0] = (xyz2[3*ip1+1] - xyz2[3*ip0+1])
                                            * (xyz2[3*ip2+2] - xyz2[3*ip0+2])
                                            - (xyz2[3*ip2+1] - xyz2[3*ip0+1])
                                            * (xyz2[3*ip1+2] - xyz2[3*ip0+2]);
                                    norm[1] = (xyz2[3*ip1+2] - xyz2[3*ip0+2])
                                            * (xyz2[3*ip2  ] - xyz2[3*ip0  ])
                                            - (xyz2[3*ip2+2] - xyz2[3*ip0+2])
                                            * (xyz2[3*ip1  ] - xyz2[3*ip0  ]);
                                    norm[2] = (xyz2[3*ip1  ] - xyz2[3*ip0  ])
                                            * (xyz2[3*ip2+1] - xyz2[3*ip0+1])
                                            - (xyz2[3*ip2  ] - xyz2[3*ip0  ])
                                            * (xyz2[3*ip1+1] - xyz2[3*ip0+1]);

                                    area = sqrt(SQR(norm[0]) + SQR(norm[1]) + SQR(norm[2]));

                                    norm[0] /= area;
                                    norm[1] /= area;
                                    norm[2] /= area;
                                    (void) fwrite(norm, sizeof(REAL32), 3, stl_fp);

                                    vert[0] = xyz2[3*ip0  ];
                                    vert[1] = xyz2[3*ip0+1];
                                    vert[2] = xyz2[3*ip0+2];
                                    (void) fwrite(vert, sizeof(REAL32), 3, stl_fp);

                                    vert[0] = xyz2[3*ip1  ];
                                    vert[1] = xyz2[3*ip1+1];
                                    vert[2] = xyz2[3*ip1+2];
                                    (void) fwrite(vert, sizeof(REAL32), 3, stl_fp);

                                    vert[0] = xyz2[3*ip2  ];
                                    vert[1] = xyz2[3*ip2+1];
                                    vert[2] = xyz2[3*ip2+2];
                                    (void) fwrite(vert, sizeof(REAL32), 3, stl_fp);

                                    icolr = 0;
                                    (void) fwrite(&icolr, sizeof(UINT16), 1, stl_fp);
                                }
                            }
                        }
                    }

                    fclose(stl_fp);

                    SPRINT1(0, "--> \"%s\" has been written", stl_filename);
                }

            /* 'd' - derivative of Parameter */
            } else if (*state == 'd') {
                SPRINT0(0, "--> Option 'd' (derivative of Parameter)");

                if (numarg > 0) {
                    ipmtr  = numarg;
                    numarg = 0;

                    status = ocsmSetVel(modl, 0, 0, 0, "0");
                    if (status != SUCCESS) {
                        SPRINT1(0, "ocsmSetVel -> status=%d", status);
                    }

                    status = ocsmSetVel(modl, ipmtr, 1, 1, "1");
                    if (status != SUCCESS) {
                        SPRINT1(0, "ocsmSetVel -> status=%d", status);
                    }
                } else {
                    while (1) {
                        /* list the external Parameters */
                        status = ocsmPrintPmtrs(MODL, stdout);
                        if (status != SUCCESS) {
                            SPRINT1(0, "ocsmPrintPmtrs -> status=%d", status);
                        }

                        /* get name of Parameter to set derivative */
                        ipmtr = getInt("Enter Parameter index: ");

                        if (ipmtr < 1 || ipmtr > MODL->npmtr) {
                            break;
                        }

                        status = ocsmGetPmtr(modl, ipmtr, &itype, &nrow, &ncol, pmtrName);
                        if (status != SUCCESS) {
                            SPRINT1(0, "ocsmGetPmtr -> status=%d", status);
                        }

                        if (nrow > 1) {
                            irow = getInt("Enter row number:      ");
                        } else {
                            irow = 1;
                        }
                        if (ncol > 1) {
                            icol = getInt("Enter col number:      ");
                        } else {
                            icol = 1;
                        }
                        getStr("Enter new derivative:  ",   defn);

                        status = ocsmSetVel(modl, ipmtr, irow, icol, defn);
                        if (status != SUCCESS) {
                            SPRINT4(0, "**> ocsmSetVel(ipmtr=%d, defn=%s) -> status=%d (%s)",
                                    ipmtr, defn, status, ocsmGetText(status));
                            goto next_option;
                        }
                    }
                }

                sclrType = -1;
                tuftLen  =  0;
                new_data =  1;

                SPRINT0(0, "Use 'B' and 'w' option to redisplay velocities");

            /* 'D' - delete Branch */
            } else if (*state == 'D') {
                SPRINT0(0, "--> Option 'D' (delete Branch)");

                status = ocsmDelBrch(MODL, MODL->nbrch);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsmDelBrch(ibrch=%d) -> status=%d (%s)",
                            MODL->nbrch, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT0(0, "Branch deleted");
                SPRINT0(0, "Use 'B' to rebuild");

            /* 'e' - edit Parameter */
            } else if (*state == 'e') {
                SPRINT0(0, "--> Option 'e' chosen (edit Parameter)");

                while (1) {
                    /* list the external Parameters */
                    status = ocsmPrintPmtrs(MODL, stdout);
                    if (status != SUCCESS) {
                        SPRINT1(0, "ocsmPrintPmtrs -> status=%d", status);
                    }

                    /* get name of Parameter to set derivative */
                    ipmtr = getInt("Enter Parameter index: ");

                    if (ipmtr < 1 || ipmtr > MODL->npmtr) {
                        break;
                    }

                    status = ocsmGetPmtr(modl, ipmtr, &itype, &nrow, &ncol, pmtrName);
                    if (status != SUCCESS) {
                        SPRINT1(0, "ocsmGetPmtr -> status=%d", status);
                    }

                    if (nrow > 1) {
                        irow = getInt("Enter row number:      ");
                    } else {
                        irow = 1;
                    }
                    if (ncol > 1) {
                        icol = getInt("Enter col number:      ");
                    } else {
                        icol = 1;
                    }
                    getStr("Enter new value:       ",   defn);

                    status = ocsmSetValu(modl, ipmtr, irow, icol, defn);
                    if (status != SUCCESS) {
                        SPRINT4(0, "**> ocsmSetValu(ipmtr=%d, defn=%s) -> status=%d (%s)",
                                ipmtr, defn, status, ocsmGetText(status));
                        goto next_option;
                    }
                }

                sclrType = -1;
                tuftLen  =  0;
                new_data =  1;

                SPRINT0(0, "Use 'B' to rebuild");

            /* 'E' - edit Branch */
            } else if (*state == 'E') {
                SPRINT0(0, "--> Option 'E' chosen (edit Branch)");

                if (numarg > 0) {
                    ibrch  = numarg;
                    numarg = 0;
                } else {
                    ibrch = getInt("Enter Branch to edit: ");
                }

                if (ibrch < 1 || ibrch > MODL->nbrch) {
                    SPRINT2(0, "Illegal ibrch=%d (should be between 1 and %d)", ibrch, MODL->nbrch);
                    goto next_option;
                }

                status = ocsmGetBrch(modl, ibrch, &itype, &iclass, &iactv,
                                     &ichld, &ileft, &irite, &narg, &nattr);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsmGetBranch(ibrch=%d) -> status=%d (%s)",
                            ibrch, status, ocsmGetText(status));
                    goto next_option;
                }

                for (iarg = 1; iarg <= narg; iarg++) {
                    status = ocsmGetArg(modl, ibrch, iarg, defn, &value, &dot);
                    if (status != SUCCESS) {
                        SPRINT4(0, "**> ocsmGetArg(ibrch=%d, iarg=%d) -> status=%d (%s)",
                                ibrch, iarg, status, ocsmGetText(status));
                        goto next_option;
                    }

                    SPRINT2(0, "Old       definition for arg %d: %s", iarg, defn);
                    getStr( "Enter new definition ('.' to unchange): ",  defn);

                    if (strcmp(defn, ".") == 0) {
                        SPRINT0(0, "Definition unchanged");
                    } else {
                        status = ocsmSetArg(modl, ibrch, iarg, defn);
                        if (status != SUCCESS) {
                            SPRINT5(0, "**> ocsmSetArg(ibrch=%d, iarg=%d, defn=%s) -> status=%d (%s)",
                                    ibrch, iarg, defn, status, ocsmGetText(status));
                            goto next_option;
                        }

                        SPRINT2(0, "New       definition for arg %d: %s", iarg,  defn);
                    }
                }

                SPRINT0(0, "Use 'B' to rebuild");

            /* 'f' - change tuft length */
            } else if (*state == 'f') {
                SPRINT0(0, "--> Option 'f' (change tuft length)");

                SPRINT1(0, "Old       tuft length: %f", tuftLen);

                if (numarg > 0) {
                    tuftLen = numarg;
                    numarg  = 0;
                } else {
                    tuftLen = getDbl("Enter new tuft length: ");
                }

                new_data =  1;

            /* 'F' - undefined option */
            } else if (*state == 'F') {
                SPRINT0(0, "--> Option 'F' (undefined)");

            /* 'g' - undefined option */
            } else if (*state == 'g') {
                SPRINT0(0, "--> Option 'g' (undefined)");

            /* 'G' - undefined option */
            } else if (*state == 'G') {
                SPRINT0(0, "--> Option 'G' (undefined)");

            /* 'h' - hide Edge or Face at cursor */
            } else if (*state == 'h') {
                uindex = pickObject(&utype);

                if (utype%10 == 1) {
                    ibody = utype / 10;
                    iedge = uindex;

                    MODL->body[ibody].edge[iedge].gratt.render = 0;
                    SPRINT2(0, "Hiding Edge %d (body %d)", iedge, ibody);

                    new_data    = 1;
                    utype_save  = utype;
                    uindex_save = uindex;
                } else if (utype%10 == 2) {
                    ibody = utype / 10;
                    iface = uindex;

                    MODL->body[ibody].face[iface].gratt.render = 0;
                    SPRINT2(0, "Hiding Face %d (body %d)", iface, ibody);

                    new_data    = 1;
                    utype_save  = utype;
                    uindex_save = uindex;
                } else {
                    SPRINT0(0, "nothing to hide");
                }

            /* 'H' - undefined option */
            } else if (*state == 'H') {
                SPRINT0(0, "--> Option 'H' (undefined)");

            /* 'i' - undefined option */
            } else if (*state == 'i') {
                SPRINT0(0, "--> Option 'i' (undefined)");

            /* 'I' - undefined option */
            } else if (*state == 'I') {
                SPRINT0(0, "--> Option 'I' (undefined)");

            /* 'j' - undefined option */
            } else if (*state == 'j') {
                SPRINT0(0, "--> Option 'j' (undefined)");

            /* 'J' - undefined option */
            } else if (*state == 'J') {
                SPRINT0(0, "--> Option 'J' (undefined)");

            /* 'k' - undefined option */
            } else if (*state == 'k') {
                SPRINT0(0, "--> Option 'k' (undefined)");

            /* 'K' - undefined option */
            } else if (*state == 'K') {
                SPRINT0(0, "--> Option 'K' (undefined)");

            /* 'l' - list Parameters */
            } else if (*state == 'l') {
                SPRINT0(0, "--> Option 'l' chosen (list Parameters)");

                status = ocsmPrintPmtrs(modl, stdout);
                SPRINT2(0, "--> ocsmPrintPmtrs() -> status=%d (%s)",
                        status, ocsmGetText(status));

            /* 'L' - list Branches */
            } else if (*state == 'L') {
                SPRINT0(0, "--> Option 'L' chosen (list Branches)");

                status = ocsmPrintBrchs(modl, stdout);
                SPRINT2(0, "--> ocsmPrintBrchs() -> status=%d (%s)",
                        status, ocsmGetText(status));

            /* 'm' - view in monochrome */
            } else if (*state == 'm') {
                SPRINT0(0, "--> Option 'm' (view in monochrome)");

                sclrType = -1;
                tuftLen  =  0;
                new_data =  1;

            /* 'M' - find unmatched Edges */
            } else if (*state == 'M') {
                SPRINT0(0, "--> Option 'M' (find unmatched Edges)");

                ngood = 0;
                for (jbody = 0; jbody < nbody; jbody++) {
                    ibody = bodyList[jbody];

                    for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                        if (MODL->body[ibody].edge[iedge].ileft <= 0 ||
                            MODL->body[ibody].edge[iedge].irite <= 0   ) {
                            EG_getMassProperties(MODL->body[ibody].edge[iedge].eedge, massprops);
                            SPRINT3(0, "Edge %d:%d is unmathced  (length=%f)",
                                    ibody, iedge, massprops[1]);
                        } else {
                            ngood++;
                        }
                    }
                }

                SPRINT1(0, "there are %d good Edges", ngood);

            /* 'n' - print sensitivity on Nodes and Edges */
            } else if (*state == 'n') {
                int     buildTo;

                SPRINT0(0, "--> Option 'n' chosen (compute sensitivity on Nodes and Edges)");

                if (numarg > 0) {
                    ipmtr = numarg;
                    numarg = 0;
                } else {
                    status = ocsmPrintPmtrs(MODL, stdout);
                    if (status != SUCCESS) {
                        SPRINT1(0, "ocsmPrintPmtrs -> status=%d", status);
                    }

                    ipmtr = getInt("Enter Pmtr for sensitivity: ");
                }

                if (ipmtr >= 1 && ipmtr <= MODL->npmtr) {
                    if (MODL->pmtr[ipmtr].type == OCSM_DESPMTR) {
                        status = ocsmSetVelD(modl, 0, 0, 0, 0.0);
                        if (status != SUCCESS) {
                            SPRINT1(0, "ocsmSetVelD -> status=%d", status);
                        }

                        status = ocsmSetVelD(modl, ipmtr, 0, 0, 1.0);
                        if (status != SUCCESS) {
                            SPRINT1(0, "ocsmSetVelD -> status=%d", status);
                        }

                        /* build needed to propagate sensitivities */
                        buildTo  = 0;
                        ntemp    = 0;
                        status   = ocsmBuild(modl, buildTo, &builtTo, &ntemp, NULL);
                        if (status != SUCCESS) {
                            SPRINT1(0, "ocsmBuild -> status=%d", status);
                        }

                        for (jbody = 1; jbody <= MODL->nbody; jbody++) {
                            if (MODL->body[jbody].onstack != 1) continue;

                            for (inode = 1; inode <= MODL->body[jbody].nnode; inode++) {
                                status = ocsmGetVel(MODL, jbody, OCSM_NODE, inode, 1, NULL, vel);
                                if (status != SUCCESS) {
                                    SPRINT1(0, "ocsmGetVel -> status=%d", status);
                                }

                                SPRINT6(0, "Node %3d:%-3d       %10.4f %10.4f %10.4f   %10.4f", jbody, inode,
                                        vel[0], vel[1], vel[2],
                                        sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]));
                            }
                        }

                        sclrType = 1;
                        if (tuftLen == 0) tuftLen  = 0.1;
                        new_data = 1;

                        SPRINT0(0, "Use 'w' option to redisplay velocities");
                    } else {
                        SPRINT1(0, "ERROR:: ipmtr=%d is not an external Pmtr\a", ipmtr);
                    }
                } else {
                    SPRINT2(0, "ERROR:: ipmtr=%d is not between 1 and %d\a", ipmtr, MODL->npmtr);
                }

            /* 'N' - name Branch */
            } else if (*state == 'N') {
                SPRINT0(0, "--> Option 'N' chosen (name Branch)");

                if (numarg > 0) {
                    ibrch  = numarg;
                    numarg = 0;
                } else {
                    ibrch = getInt("Enter Branch to rename: ");
                }

                if (ibrch < 1 || ibrch > MODL->nbrch) {
                    SPRINT2(0, "Illegal ibrch=%d (should be between 1 and %d)", ibrch, MODL->nbrch);
                    goto next_option;
                }

                status = ocsmGetName(modl, ibrch, brchName);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsmGetName(ibrch=%d) -> status=%d (%s)",
                            ibrch, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT2(0, "--> Name of Branch %d is %s", ibrch, brchName);

                getStr("Enter new Branch name (. for none): ", brchName);
                if (strcmp(brchName, ".") == 0) {
                    SPRINT1(0, "Branch %4d has not been renamed", ibrch);
                    goto next_option;
                }

                status = ocsmSetName(modl, ibrch, brchName);
                if (status != SUCCESS) {
                    SPRINT4(0, "**> ocsmSetName(ibrch=%d, name=%s) -> status=%d (%s)",
                            ibrch, brchName, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT1(0, "Branch %4d has been renamed", ibrch);

            /* 'o' - undefined option */
            } else if (*state == 'o') {
                SPRINT0(0, "--> Option 'o' (undefined)");

            /* 'O' - undefined option */
            } else if (*state == 'O') {
                SPRINT0(0, "--> Option 'O' (undefined)");

            /* 'p' - get parametric coordinates */
            } else if (*state == 'p') {
                SPRINT0(0, "--> Option 'p' (get parametric coordinates)");

                ibody  = getInt("Enter ibody: ");
                iface  = getInt("Enter iface: ");
                xyz[0] = getDbl("Enter x:     ");
                xyz[1] = getDbl("Enter y:     ");
                xyz[2] = getDbl("Enter z:     ");

                status = ocsmGetUV(MODL, ibody, OCSM_FACE, iface, 1, xyz, uv);
                if (status == SUCCESS) {
                    SPRINT1(0, "u = %12.5e", uv[0]);
                    SPRINT1(0, "v = %12.5e", uv[1]);
                } else {
                    SPRINT4(0, "**> ocsmGetUV(ibody=%d, iface=%d) -> status=%d (%s)",
                            ibody, iface, status, ocsmGetText(status));
                }

            /* 'P' - get physical coordinates */
            } else if (*state == 'P') {
                SPRINT0(0, "--> Option 'P' (get physical coordinates)");

                ibody  = getInt("Enter ibody: ");
                iface  = getInt("Enter iface: ");
                uv[0]  = getDbl("Enter u:     ");
                uv[1]  = getDbl("Enter v:     ");

                status = ocsmGetXYZ(MODL, ibody, OCSM_FACE, iface, 1, uv, xyz);
                if (status == SUCCESS) {
                    SPRINT1(0, "x   = %12.5e", xyz[0]);
                    SPRINT1(0, "y   = %12.5e", xyz[1]);
                    SPRINT1(0, "z   = %12.5e", xyz[2]);
                } else {
                    SPRINT4(0, "**> ocsmGetXYZ(ibody=%d, iface=%d) -> status=%d (%s)",
                            ibody, iface, status, ocsmGetText(status));
                }

                status = ocsmGetVel(MODL, ibody, OCSM_FACE, iface, 1, uv, vel);
                if (status == SUCCESS) {
                    SPRINT3(0, "vel = %12.5e %12.5e %12.5e", vel[0], vel[1], vel[2]);
                } else {
                    SPRINT4(0, "**> ocsmGetVel(ibody=%d, iface=%d) -> status=%d (%s)",
                            ibody, iface, status, ocsmGetText(status));
                }

            /* 'q' - query Edge/Face at cursor */
            } else if (*state == 'q') {
                ego  eedge, eface;
                SPRINT0(0, "--> Option q chosen (query Edge/Face at cursor) ");

                uindex = pickObject(&utype);

                /* Edge found */
                if        (utype%10 == 1) {
                    ibody = utype / 10;
                    iedge = uindex;

                    SPRINT2(0, "Body %4d Edge %4d:", ibody, iedge);

                    eedge  = MODL->body[ibody].edge[iedge].eedge;
                    status = EG_attributeNum(eedge, &nattr);
                    if (status != SUCCESS) {
                        SPRINT1(0, "EG_attributeNum -> status=%d", status);
                    }

                    for (iattr = 1; iattr <= nattr; iattr++) {
                        status = EG_attributeGet(eedge, iattr, &attrName, &itype, &nlist,
                                                 &tempIlist, &tempRlist, &tempClist);
                        if (status != SUCCESS) {
                            SPRINT1(0, "EG_attributeGet -> status=%d", status);
                        }
                        SPRINT1x(0, "                     %-20s =", attrName);

                        if        (itype == ATTRINT) {
                            for (i = 0; i < nlist ; i++) {
                                SPRINT1x(0, "%5d ", tempIlist[i]);
                            }
                            SPRINT0(0, " ");
                        } else if (itype == ATTRREAL) {
                            for (i = 0; i < nlist ; i++) {
                                SPRINT1x(0, "%11.5f ", tempRlist[i]);
                            }
                            SPRINT0(0, " ");
                        } else if (itype == ATTRSTRING) {
                            SPRINT1(0, "%s", tempClist);
                        }
                    }

                /* Face found */
                } else if (utype%10 == 2) {
                    ibody = utype / 10;
                    iface = uindex;

                    SPRINT2(0, "Body %4d Face %4d:", ibody, iface);

                    eface  = MODL->body[ibody].face[iface].eface;
                    status = EG_attributeNum(eface, &nattr);
                    if (status != SUCCESS) {
                        SPRINT1(0, "EG_attributeNum -> status=%d", status);
                    }

                    for (iattr = 1; iattr <= nattr; iattr++) {
                        status = EG_attributeGet(eface, iattr, &attrName, &itype, &nlist,
                                                 &tempIlist, &tempRlist, &tempClist);
                        if (status != SUCCESS) {
                            SPRINT1(0, "EG_attributeGet -> status=%d", status);
                        }

                        SPRINT1x(0, "                     %-20s =", attrName);

                        if        (itype == ATTRINT) {
                            for (i = 0; i < nlist ; i++) {
                                SPRINT1x(0, "%5d ", tempIlist[i]);
                            }
                            SPRINT0(0, " ");
                        } else if (itype == ATTRREAL) {
                            for (i = 0; i < nlist ; i++) {
                                SPRINT1x(0, "%11.5f ", tempRlist[i]);
                            }
                            SPRINT0(0, " ");
                        } else if (itype == ATTRSTRING) {
                            SPRINT1(0, "%s", tempClist);
                        }
                    }

                /* other (Axis) found */
                } else {
                    SPRINT0(0, "Nothing found");
                }

                numarg = 0;

            /* 'Q' - query all attributes */
            } else if (*state == 'Q') {
                SPRINT0(0, "--> Option 'Q' chosen (quary all attributes)");

                for (ibody = 1; ibody <= MODL->nbody; ibody++) {
                    if (MODL->body[ibody].onstack != 1) continue;

                    SPRINT1(0, "ibody     %5d", ibody);
                    status = EG_attributeNum(MODL->body[ibody].ebody, &nattr);
                    if (status < SUCCESS) {
                        SPRINT1(0, "EG_attributeNum -> status=%d", status);
                    }

                    for (iattr = 1; iattr <= nattr; iattr++) {
                        status = EG_attributeGet(MODL->body[ibody].ebody, iattr,
                                                 &attrName, &itype, &nlist,
                                                 &tempIlist, &tempRlist, &tempClist);
                        if (status < SUCCESS) {
                            SPRINT1(0, "EG_attributeGet -> status=%d", status);
                        }

                        SPRINT1x(0, "                     %-20s =", attrName);

                        if        (itype == ATTRINT) {
                            for (i = 0; i < nlist ; i++) {
                                SPRINT1x(0, "%5d ", tempIlist[i]);
                            }
                            SPRINT0(0, " ");
                        } else if (itype == ATTRREAL) {
                            for (i = 0; i < nlist ; i++) {
                                SPRINT1x(0, "%11.5f ", tempRlist[i]);
                            }
                            SPRINT0(0, " ");
                        } else if (itype == ATTRSTRING) {
                            SPRINT1(0, " %s", tempClist);
                        }
                    }

                    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                        SPRINT1(0, "    iface %5d", iface);

                        status = EG_attributeNum(MODL->body[ibody].face[iface].eface, &nattr);
                        if (status < SUCCESS) {
                            SPRINT1(0, "EG_attributeNum -> status=%d", status);
                        }

                        for (iattr = 1; iattr <= nattr; iattr++) {
                            status = EG_attributeGet(MODL->body[ibody].face[iface].eface, iattr,
                                                     &attrName, &itype, &nlist,
                                                     &tempIlist, &tempRlist, &tempClist);
                            if (status < SUCCESS) {
                                SPRINT1(0, "EG_attributeGet -> status=%d", status);
                            }

                            SPRINT1x(0, "                     %-20s =", attrName);

                            if        (itype == ATTRINT) {
                                for (i = 0; i < nlist ; i++) {
                                    SPRINT1x(0, "%5d ", tempIlist[i]);
                                }
                                SPRINT0(0, " ");
                            } else if (itype == ATTRREAL) {
                                for (i = 0; i < nlist ; i++) {
                                    SPRINT1x(0, "%11.5f ", tempRlist[i]);
                                }
                                SPRINT0(0, " ");
                            } else if (itype == ATTRSTRING) {
                                SPRINT1(0, " %s", tempClist);
                            }
                        }
                    }

                    for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                        SPRINT1(0, "    iedge %5d", iedge);

                        status = EG_attributeNum(MODL->body[ibody].edge[iedge].eedge, &nattr);
                        if (status < SUCCESS) {
                            SPRINT1(0, "EG_attributeNum -> status=%d", status);
                        }

                        for (iattr = 1; iattr <= nattr; iattr++) {
                            status = EG_attributeGet(MODL->body[ibody].edge[iedge].eedge, iattr,
                                                     &attrName, &itype, &nlist,
                                                     &tempIlist, &tempRlist, &tempClist);
                            if (status < SUCCESS) {
                                SPRINT1(0, "EG_attributeGet -> status=%d", status);
                            }

                            SPRINT1x(0, "                     %-20s =", attrName);

                            if        (itype == ATTRINT) {
                                for (i = 0; i < nlist ; i++) {
                                    SPRINT1x(0, "%5d ", tempIlist[i]);
                                }
                                SPRINT0(0, " ");
                            } else if (itype == ATTRREAL) {
                                for (i = 0; i < nlist ; i++) {
                                    SPRINT1x(0, "%11.5f ", tempRlist[i]);
                                }
                                SPRINT0(0, " ");
                            } else if (itype == ATTRSTRING) {
                                SPRINT1(0, " %s", tempClist);
                            }
                        }
                    }
                }

            /* 'r' - undefined option */
            } else if (*state == 'r') {
                SPRINT0(0, "--> Option 'r' (undefined)");

            /* 'R' - resume a Branch */
            } else if (*state == 'R') {
                SPRINT0(0, "--> Option 'R' chosen (resume a Branch)");

                if (numarg > 0) {
                    ibrch  = numarg;
                    numarg = 0;
                } else {
                    ibrch = getInt("Enter Branch to resume (9999 for all): ");
                }

                if (ibrch == 9999) {
                    for (ibrch = 1; ibrch <= MODL->nbrch; ibrch++) {
                        status = ocsmSetBrch(modl, ibrch, OCSM_ACTIVE);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmSetBrch(ibrch=%d) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    }

                    SPRINT0(0, "All Branches have been resumed");
                    goto next_option;
                }

                if (ibrch < 1 || ibrch > MODL->nbrch) {
                    SPRINT2(0, "Illegal ibrch=%d (should be between 1 and %d)", ibrch, MODL->nbrch);
                    goto next_option;
                }

                status = ocsmSetBrch(modl, ibrch, OCSM_ACTIVE);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsmSetBrch(ibrch=%d) -> status=%d (%s)",
                            ibrch, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT1(0, "Brch %4d has been resumed", ibrch);
                SPRINT0(0, "Use 'B' to rebuild");

            /* 's' - compare sensitivities on Faces, Edges, and Nodes */
            } else if (*state == 's') {
                int     buildTo, npnt_tess, ntri_tess, ipnt, nerror, ntotal;
                CINT    *ptype, *pindx, *tris, *tric;
                double  **face_vel_anal, **face_vel_fd, face_vel_err;
                double  **edge_vel_anal, **edge_vel_fd, edge_vel_err;
                double  **node_vel_anal, **node_vel_fd, node_vel_err;
                double  errmax, face_errmax, edge_errmax, node_errmax, EPS05=1.0e-5;
                CDOUBLE *xyz2, *uv2;

                SPRINT0(0, "--> Option 's' chosen (compare sensitivities on Faces, Edges, and Nodes)");

                ntotal = 0;
                errmax = 0;
                for (ipmtr = 1; ipmtr <= MODL->npmtr; ipmtr++) {
                    if (MODL->pmtr[ipmtr].type != OCSM_DESPMTR) continue;
                    if (numarg > 0 && ipmtr != numarg) continue;

                    for (irow = 1; irow <= MODL->pmtr[ipmtr].nrow; irow++) {
                        for (icol = 1; icol <= MODL->pmtr[ipmtr].ncol; icol++) {

                            status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
                            if (status < 0) {
                                SPRINT1(0, "ERROR:: ocsmSetVelD -> status=%d", status);
                                goto next_option;
                            }
                            status = ocsmSetVelD(MODL, ipmtr, irow, icol, 1.0);
                            if (status < 0) {
                                SPRINT1(0, "ERROR:: ocsmSetVelD -> status=%d", status);
                                goto next_option;
                            }

                            /* build needed to propagate sensitivities */
                            SPRINT3(0, "Propagating sensitivities of parameters for \"%s[%d,%d]\"",
                                    MODL->pmtr[ipmtr].name, irow, icol);
                            buildTo  = 0;
                            ntemp    = 0;
                            status   = ocsmBuild(MODL, buildTo, &builtTo, &ntemp, NULL);
                            if (status < 0) {
                                SPRINT1(0, "ERROR:: ocsmBuild -> status=%d", status);
                                goto next_option;
                            }

                            /* perform checks for each Body on the stack */
                            for (jbody = 1; jbody <= MODL->nbody; jbody++) {
                                if (MODL->body[jbody].onstack != 1) continue;
                                nerror = 0;

                                /* analytic sensitivities (if possible) */
                                SPRINT1(0, "Computing analytic sensitivities (if possible) for jbody=%d", jbody);
                                status = ocsmSetDtime(MODL, 0);
                                if (status < 0) {
                                    SPRINT1(0, "ERROR:: ocsmSetDtime -> status=%d", status);
                                    goto next_option;
                                }

                                face_vel_anal = (double**) malloc((MODL->body[jbody].nface+1)*sizeof(double*));
                                for (iface = 1; iface <= MODL->body[jbody].nface; iface++) {
                                    status = EG_getTessFace(MODL->body[jbody].etess, iface,
                                                            &npnt_tess, &xyz2, &uv2, &ptype, &pindx,
                                                            &ntri_tess, &tris, &tric);
                                    if (status < 0 || npnt_tess <= 0) {
                                        SPRINT2(0, "ERROR:: EG_getTessFace -> status=%d, npnt_tess=%d", status, npnt_tess);
                                        goto next_option;
                                    }

                                    face_vel_anal[iface] = (double*) malloc(3*npnt_tess*sizeof(double));
                                    status = ocsmGetVel(MODL, jbody, OCSM_FACE, iface, npnt_tess, NULL, face_vel_anal[iface]);
                                    if (status < 0) {
                                        SPRINT1(0, "ERROR:: ocsmGetVel -> status=%d", status);
                                        goto next_option;
                                    }
                                    if (status != SUCCESS) nerror++;
                                }

                                edge_vel_anal = (double**) malloc((MODL->body[jbody].nedge+1)*sizeof(double*));
                                for (iedge = 1; iedge <= MODL->body[jbody].nedge; iedge++) {
                                    status = EG_getTessEdge(MODL->body[jbody].etess, iedge,
                                                            &npnt_tess, &xyz2, &uv2);
                                    if (status < 0 || npnt_tess <= 0) {
                                        SPRINT2(0, "ERROR:: EG_getTessEdge -> status=%d, npnt_tess=%d", status, npnt_tess);
                                        goto next_option;
                                    }

                                    edge_vel_anal[iedge] = (double*) malloc(3*npnt_tess*sizeof(double));
                                    status = ocsmGetVel(MODL, jbody, OCSM_EDGE, iedge, npnt_tess, NULL, edge_vel_anal[iedge]);
                                    if (status < 0) {
                                        SPRINT1(0, "ERROR:: ocsmGetVel -> status=%d", status);
                                        goto next_option;
                                    }
                                    if (status != SUCCESS) nerror++;
                                }

                                node_vel_anal = (double**) malloc((MODL->body[jbody].nnode+1)*sizeof(double*));
                                for (inode = 1; inode <= MODL->body[jbody].nnode; inode++) {
                                    npnt_tess = 1;

                                    node_vel_anal[inode] = (double*) malloc(3*npnt_tess*sizeof(double));
                                    status = ocsmGetVel(MODL, jbody, OCSM_NODE, inode, npnt_tess, NULL, node_vel_anal[inode]);
                                    if (status < 0) {
                                        SPRINT1(0, "ERROR:: ocsmGetVel -> status=%d", status);
                                        goto next_option;
                                    }
                                    if (status != SUCCESS) nerror++;
                                }

                                /* finite difference sensitivities */
                                SPRINT1(0, "Computing finite difference sensitivities for jbody=%d", jbody);
                                status = ocsmSetDtime(MODL, 0.001);
                                if (status != SUCCESS) {
                                    SPRINT1(0, "ocsmSetDtime -> status=%d", status);
                                }

                                face_vel_fd   = (double**) malloc((MODL->body[jbody].nface+1)*sizeof(double*));
                                for (iface = 1; iface <= MODL->body[jbody].nface; iface++) {
                                    status = EG_getTessFace(MODL->body[jbody].etess, iface,
                                                            &npnt_tess, &xyz2, &uv2, &ptype, &pindx,
                                                            &ntri_tess, &tris, &tric);
                                    if (status != SUCCESS) {
                                        SPRINT1(0, "EG_getTessFace -> status=%d", status);
                                    }

                                    face_vel_fd[iface] = (double*) malloc(3*npnt_tess*sizeof(double));
                                    status = ocsmGetVel(MODL, jbody, OCSM_FACE, iface, npnt_tess, NULL, face_vel_fd[iface]);
                                    if (status != SUCCESS) nerror++;
                                }

                                edge_vel_fd   = (double**) malloc((MODL->body[jbody].nedge+1)*sizeof(double*));
                                for (iedge = 1; iedge <= MODL->body[jbody].nedge; iedge++) {
                                    status = EG_getTessEdge(MODL->body[jbody].etess, iedge,
                                                            &npnt_tess, &xyz2, &uv2);
                                    if (status != SUCCESS) {
                                        SPRINT1(0, "EG_getTessEdge -> status=%d", status);
                                    }

                                    edge_vel_fd[iedge] = (double*) malloc(3*npnt_tess*sizeof(double));
                                    status = ocsmGetVel(MODL, jbody, OCSM_EDGE, iedge, npnt_tess, NULL, edge_vel_fd[iedge]);
                                    if (status != SUCCESS) nerror++;
                                }

                                node_vel_fd   = (double**) malloc((MODL->body[jbody].nnode+1)*sizeof(double*));
                                for (inode = 1; inode <= MODL->body[jbody].nnode; inode++) {
                                    npnt_tess = 1;

                                    node_vel_fd[inode] = (double*) malloc(3*npnt_tess*sizeof(double));
                                    status = ocsmGetVel(MODL, jbody, OCSM_NODE, inode, npnt_tess, NULL, node_vel_fd[inode]);
                                    if (status != SUCCESS) nerror++;
                                }

                                status = ocsmSetDtime(MODL, 0);
                                if (status != SUCCESS) {
                                    SPRINT1(0, "ocsmSetDtime -> status=%d", status);
                                }

                                SPRINT0(0, "Removing perturbation");

                                if (nerror > 0) {
                                    SPRINT1(0, "WARNING:: Sensitivities not compared since %d errors were detected during setup", nerror);
                                    ntotal++;
                                    break;
                                }

                                /* compare sensitivities */
                                node_errmax = 0;
                                edge_errmax = 0;
                                face_errmax = 0;

                                SPRINT4(0, "Comparing sensitivities for \"%s[%d,%d]\" for jbody=%d",
                                        MODL->pmtr[ipmtr].name, irow, icol, jbody);
                                nerror = 0;
                                for (iface = 1; iface <= MODL->body[jbody].nface; iface++) {
                                    status = EG_getTessFace(MODL->body[jbody].etess, iface,
                                                            &npnt_tess, &xyz2, &uv2, &ptype, &pindx,
                                                            &ntri_tess, &tris, &tric);
                                    if (status != SUCCESS) {
                                        SPRINT1(0, "EG_getTessFace -> status=%d", status);
                                    }

                                    for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                                        for (i = 0; i < 3; i++) {
                                            face_vel_err = face_vel_anal[iface][3*ipnt+i] - face_vel_fd[iface][3*ipnt+i];
                                            if (fabs(face_vel_err) > face_errmax) face_errmax = fabs(face_vel_err);
                                            if (fabs(face_vel_err) > EPS05) {
                                                if (nerror < 20 || fabs(face_vel_err) >= face_errmax) {
                                                    SPRINT8(0, "iface=%4d,  ipnt=%4d,    anal=%16.8f,  fd=%16.8f,  err=%16.8f (at %10.4f %10.4f %10.4f)",
                                                            iface, ipnt, face_vel_anal[iface][3*ipnt+i], face_vel_fd[iface][3*ipnt+i], face_vel_err,
                                                            xyz2[3*ipnt], xyz2[3*ipnt+1], xyz2[3*ipnt+2]);
                                                }
                                                nerror++;
                                                ntotal++;
                                            }
                                        }
                                    }

                                    free(face_vel_anal[iface]);
                                    free(face_vel_fd[  iface]);
                                }
                                errmax = MAX(errmax, face_errmax);
                                SPRINT3(0, "    d(Face)/d(%s) check complete with %d total errors (errmax=%12.4e)",
                                        MODL->pmtr[ipmtr].name, nerror, face_errmax);

                                nerror = 0;
                                for (iedge = 1; iedge <= MODL->body[jbody].nedge; iedge++) {
                                    status = EG_getTessEdge(MODL->body[jbody].etess, iedge,
                                                            &npnt_tess, &xyz2, &uv2);
                                    if (status != SUCCESS) {
                                        SPRINT1(0, "EG_getTessEdge -> status=%d", status);
                                    }

                                    for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                                        for (i = 0; i < 3; i++) {
                                            edge_vel_err = edge_vel_anal[iedge][3*ipnt+i] - edge_vel_fd[iedge][3*ipnt+i];
                                            if (fabs(edge_vel_err) > edge_errmax) edge_errmax = fabs(edge_vel_err);
                                            if (fabs(edge_vel_err) > EPS05) {
                                                if (nerror < 20 || fabs(edge_vel_err) >= edge_errmax) {
                                                    SPRINT9(0, "iedge=%4d,  ipnt=%4d,  %1d anal=%16.8f,  fd=%16.8f,  err=%16.8f (at %10.4f %10.4f %10.4f)",
                                                            iedge, ipnt, i, edge_vel_anal[iedge][3*ipnt+i],
                                                            edge_vel_fd[  iedge][3*ipnt+i], edge_vel_err,
                                                            xyz2[3*ipnt], xyz2[3*ipnt+1], xyz2[3*ipnt+2]);
                                                }
                                                nerror++;
                                                ntotal++;
                                            }
                                        }
                                    }

                                    free(edge_vel_anal[iedge]);
                                    free(edge_vel_fd[  iedge]);
                                }
                                errmax = MAX(errmax, edge_errmax);
                                SPRINT3(0, "    d(Edge)/d(%s) check complete with %d total errors (errmax=%12.4e)",
                                        MODL->pmtr[ipmtr].name, nerror, edge_errmax);

                                nerror = 0;
                                for (inode = 1; inode <= MODL->body[jbody].nnode; inode++) {
                                    npnt_tess = 1;

                                    for (ipnt = 0; ipnt < npnt_tess; ipnt++) {
                                        for (i = 0; i < 3; i++) {
                                            node_vel_err = node_vel_anal[inode][3*ipnt+i] - node_vel_fd[inode][3*ipnt+i];
                                            if (fabs(node_vel_err) > node_errmax) node_errmax = fabs(node_vel_err);
                                            if (fabs(node_vel_err) > EPS05) {
                                                if (nerror < 20 || fabs(node_vel_err) >= node_errmax) {
                                                    SPRINT9(0, "inode=%4d,  ipnt=%4d,  %1d anal=%16.8f,  fd=%16.8f,  err=%16.8f (at %10.4f %10.4f %10.4f)",
                                                            inode, ipnt, i, node_vel_anal[inode][3*ipnt+i], node_vel_fd[inode][3*ipnt+i], node_vel_err,
                                                            MODL->body[jbody].node[inode].x,
                                                            MODL->body[jbody].node[inode].y,
                                                            MODL->body[jbody].node[inode].z);
                                                }
                                                nerror++;
                                                ntotal++;
                                            }
                                        }
                                    }

                                    free(node_vel_anal[inode]);
                                    free(node_vel_fd[  inode]);
                                }
                                errmax = MAX(errmax, node_errmax);
                                SPRINT3(0, "    d(Node)/d(%s) check complete with %d total errors (errmax=%12.4e)",
                                        MODL->pmtr[ipmtr].name, nerror, node_errmax);

                                free(face_vel_anal);
                                free(face_vel_fd  );
                                free(edge_vel_anal);
                                free(edge_vel_fd  );
                                free(node_vel_anal);
                                free(node_vel_fd  );
                            }

                            status = ocsmSetVelD(MODL, 0, 0, 0, 0.0);
                            if (status != SUCCESS) {
                                SPRINT1(0, "ocsmSetVelD -> status=%d", status);
                            }
                        }
                    }
                }

                if (errmax < 1e-20) {
                    SPRINT1(0, "\nSensitivity checks complete with %8d total errors (errmax=            )",
                            ntotal);
                } else {
                    SPRINT2(0, "\nSensitivity checks complete with %8d total errors (errmax=%12.4e)",
                            ntotal, errmax);
                }

                numarg   =  0;
                sclrType =  1;
                tuftLen  =  0;
                new_data =  1;

            } else if (*state == 'S') {
                SPRINT0(0, "--> Option 'S' chosen (suppress a Branch)");

                if (numarg > 0) {
                    ibrch  = numarg;
                    numarg = 0;
                } else {
                    ibrch = getInt("Enter Branch to suppress: ");
                }

                if (ibrch < 1 || ibrch > MODL->nbrch) {
                    SPRINT2(0, "Illegal ibrch=%d (should be between 1 and %d)", ibrch, MODL->nbrch);
                    goto next_option;
                }

                status = ocsmSetBrch(modl, ibrch, OCSM_SUPPRESSED);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsmSetBrch(ibrch=%d) -> status=%d (%s)",
                            ibrch, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT1(0, "Branch %4d has been suppressed", ibrch);
                SPRINT0(0, "Use 'B' to rebuild");

            /* 't' - write .topo file */
            } else if (*state == 't') {
                char tempname[255];
                FILE *fp2;
                SPRINT0(0, "--> Option 't' chosen (write .topo file)");

                strcpy(tempname, casename);
                strcat(tempname, ".topo");
                fp2 = fopen(tempname, "w");

                for (jbody = 0; jbody < nbody; jbody++) {
                    ibody = bodyList[jbody];
                    fprintf(fp2, "Body %d\n", ibody);

                    fprintf(fp2, "inode nedge      x           y           z\n");
                    for (inode = 1; inode <= MODL->body[ibody].nnode; inode++) {
                        fprintf(fp2, "%5d %5d %11.4f %11.4f %11.4f\n", inode,
                                MODL->body[ibody].node[inode].nedge,
                                MODL->body[ibody].node[inode].x,
                                MODL->body[ibody].node[inode].y,
                                MODL->body[ibody].node[inode].z);
                    }

                    fprintf(fp2, "iedge  ibeg  iend ileft irite nface ibody iford\n");
                    for (iedge = 1; iedge <= MODL->body[ibody].nedge; iedge++) {
                        fprintf(fp2, "%5d %5d %5d %5d %5d %5d %5d %5d\n", iedge,
                                MODL->body[ibody].edge[iedge].ibeg,
                                MODL->body[ibody].edge[iedge].iend,
                                MODL->body[ibody].edge[iedge].ileft,
                                MODL->body[ibody].edge[iedge].irite,
                                MODL->body[ibody].edge[iedge].nface,
                                MODL->body[ibody].edge[iedge].ibody,
                                MODL->body[ibody].edge[iedge].iford);
                    }

                    fprintf(fp2, "iface ibody iford\n");
                    for (iface = 1; iface <= MODL->body[ibody].nface; iface++) {
                        fprintf(fp2, "%5d %5d %5d\n", iface,
                                MODL->body[ibody].face[iface].ibody,
                                MODL->body[ibody].face[iface].iford);
                    }
                }

                fclose(fp2);

                SPRINT1(0, "--> Option 't' (\"%s\" has been written)", tempname);

            /* 'T' - attribute Branch */
            } else if (*state == 'T') {
                SPRINT0(0, "--> Option 'T' (attribute Branch)");

                if (numarg > 0) {
                    ibrch  = numarg;
                    numarg = 0;
                } else {
                    ibrch = getInt("Enter Branch to attribute: ");
                }

                if (ibrch < 1 || ibrch > MODL->nbrch) {
                    SPRINT2(0, "Illegal ibrch=%d (should be between 1 and %d)", ibrch, MODL->nbrch);
                    goto next_option;
                }

                status = ocsmGetBrch(modl, ibrch, &itype, &iclass, &iactv,
                                     &ichld, &ileft, &irite, &narg, &nattr);
                if (status != SUCCESS) {
                    SPRINT3(0, "**> ocsmGetBrch(ibrch=%d) -> status=%d (%s)",
                            ibrch, status, ocsmGetText(status));
                    goto next_option;
                }

                for (iattr = 1; iattr <= nattr; iattr++) {
                    status = ocsmRetAttr(modl, ibrch, iattr, aName, aValue);
                    if (status != SUCCESS) {
                        SPRINT4(0, "**> ocsmRetAttr(ibrch=%d, iattr=%d) -> status=%d (%s)",
                                ibrch, iattr, status, ocsmGetText(status));
                        goto next_option;
                    }

                    SPRINT2(0, "   %-24s=%s", aName, aValue);
                }

                getStr("Enter Attribute name (. for none): ", aName );

                if (strcmp(aName, ".") == 0) {
                    SPRINT0(0, "Attribute has not been saved");
                    goto next_option;
                }

                getStr("Enter Attribute value            : ", aValue);
                status = ocsmSetAttr(modl, ibrch, aName, aValue);
                if (status != SUCCESS) {
                    SPRINT4(0, "**> ocsmSetAttr(ibrch=%d, aName=%s) -> status=%d (%s)",
                            ibrch, aName, status, ocsmGetText(status));
                    goto next_option;
                }

                SPRINT1(0, "Attribute '%s' has been saved", aName);

            /* 'u' - gv scalar  */
            } else if (*state == 'u') {
                SPRINT0(0, "--> Option 'u' (color by u parameter)");

                sclrType = 0;
                tuftLen  = 0;
                new_data = 1;

            /* 'U' - unhide last hidden */
            } else if (*state == 'U') {

                if        (utype_save == 0) {
                    SPRINT0(0, "nothing to unhide");
                } else if (utype_save%10 == 1) {
                    ibody = utype_save / 10;
                    iedge = uindex_save;

                    MODL->body[ibody].edge[iedge].gratt.render = 2 + 64;
                    SPRINT2(0, "Unhiding Edge %d (body %d)", iedge, ibody);

                    new_data = 1;
                } else if (utype_save%10 == 2) {
                    ibody = utype_save / 10;
                    iface = uindex_save;

                    MODL->body[ibody].face[iface].gratt.render = 2 + 4 + 64;
                    SPRINT2(0, "Unhiding Face %d (body %d)", iface, ibody);

                    new_data = 1;
                } else {
                    SPRINT0(0, "nothing to unhide");
                }

                utype_save = 0;
                utype_save = 0;

            /* 'v' - undefined option */
            } else if (*state == 'v') {
                SPRINT0(0, "--> Option 'v' (color by v parameter)");

                sclrType = 0;
                tuftLen  = 0;
                new_data = 1;

            /* 'V' - paste Branches */
            } else if (*state == 'V') {
                SPRINT0(0, "--> Option 'V' (paste Branches)");

                if (npaste <= 0) {
                    SPRINT0(0, "Nothing to paste");
                    goto next_option;
                }

                for (ipaste = npaste-1; ipaste >= 0; ipaste--) {
                    status = ocsmNewBrch(MODL, MODL->nbrch,
                                               type_paste[ipaste],
                                               "<ESP>", -1,
                                               str1_paste[ipaste],
                                               str2_paste[ipaste],
                                               str3_paste[ipaste],
                                               str4_paste[ipaste],
                                               str5_paste[ipaste],
                                               str6_paste[ipaste],
                                               str7_paste[ipaste],
                                               str8_paste[ipaste],
                                               str9_paste[ipaste]);
                    if (status != SUCCESS) {
                        SPRINT3(0, "**> ocsmNewBrch(type=%d) -> status=%d (%s)",
                                type_paste[ipaste], status, ocsmGetText(status));
                        goto next_option;
                    }

                    if (strncmp(name_paste[ipaste], "Brch_", 5) != 0) {
                        status = ocsmSetName(MODL, MODL->nbrch, name_paste[ipaste]);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmSetName(ibrch=%d) -> status=%d (%s)",
                                    MODL->nbrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    }

                    SPRINT1(0, "New Branch (%d) added from paste buffer", type_paste[ipaste]);
                }

                SPRINT0(0, "Use 'B' to rebuild");

            /* 'w' - gv scalar */
            } else if (*state == 'w') {
                SPRINT0(0, "--> Option 'w' (color by design velocity)");

                sclrType = 1;
                tuftLen  = 0;
                new_data = 1;

            /* 'W' - write .csm file */
            } else if (*state == 'W') {
                SPRINT0(0, "--> Option 'W' chosen (write .csm file)");

                getStr("Enter filename: ", fileName);
                if (strstr(fileName, ".csm") == NULL) {
                    strcat(fileName, ".csm");
                }

                status = ocsmSave(MODL, fileName);
                SPRINT3(0, "--> ocsmSave(%s) -> status=%d (%s)",
                        fileName, status, ocsmGetText(status));

            /* 'x' - undefined option */
            } else if (*state == 'x') {
                SPRINT0(0, "--> Option 'x' (undefined)");

            /* 'X' - cut Branches */
            } else if (*state == 'X') {
                SPRINT0(0, "--> Option 'X' (cut Branches)");

                /* remove previous contents from paste buffer */
                npaste = 0;

                if (numarg > 0) {
                    npaste = numarg;
                    numarg = 0;
                } else {
                    npaste = getInt("Enter number of Branches to cut: ");
                }

                if (npaste > MAX_PASTE) {
                    SPRINT2(0, "Illegal npaste=%d (should be between 1 and %d)", npaste, MAX_PASTE);
                    npaste = 0;
                    goto next_option;
                }

                if (npaste < 1 || npaste > MODL->nbrch) {
                    SPRINT2(0, "Illegal npaste=%d (should be between 1 and %d)", npaste, MODL->nbrch);
                    npaste = 0;
                    goto next_option;
                }

                for (ipaste = 0; ipaste < npaste; ipaste++) {
                    ibrch = MODL->nbrch;

                    status = ocsmGetBrch(MODL, ibrch, &type_paste[ipaste], &iclass, &iactv,
                                         &ichld, &ileft, &irite, &narg, &nattr);
                    if (status != SUCCESS) {
                        SPRINT3(0, "**> ocsmGetBrch(ibrch=%d) -> status=%d (%s)",
                                ibrch, status, ocsmGetText(status));
                        goto next_option;
                    }

                    status = ocsmGetName(MODL, ibrch, name_paste[ipaste]);
                    if (status != SUCCESS) {
                        SPRINT3(0, "**> ocsmGetName(ibrch=%d) => status=%d (%s)",
                                ibrch, status, ocsmGetText(status));
                        goto next_option;
                    }

                    if (narg >= 1) {
                        status = ocsmGetArg(modl, ibrch, 1, str1_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=1) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str1_paste[ipaste], "");
                    }
                    if (narg >= 2) {
                        status = ocsmGetArg(modl, ibrch, 2, str2_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=2) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str2_paste[ipaste], "");
                    }
                    if (narg >= 3) {
                        status = ocsmGetArg(modl, ibrch, 3, str3_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=3) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str3_paste[ipaste], "");
                    }
                    if (narg >= 4) {
                        status = ocsmGetArg(modl, ibrch, 4, str4_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=4) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str4_paste[ipaste], "");
                    }
                    if (narg >= 5) {
                        status = ocsmGetArg(modl, ibrch, 5, str5_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=5) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str5_paste[ipaste], "");
                    }
                    if (narg >= 6) {
                        status = ocsmGetArg(modl, ibrch, 6, str6_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=6) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str6_paste[ipaste], "");
                    }
                    if (narg >= 7) {
                        status = ocsmGetArg(modl, ibrch, 7, str7_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=7) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str7_paste[ipaste], "");
                    }
                    if (narg >= 8) {
                        status = ocsmGetArg(modl, ibrch, 8, str8_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=8) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str8_paste[ipaste], "");
                    }
                    if (narg >= 9) {
                        status = ocsmGetArg(modl, ibrch, 9, str9_paste[ipaste], &value, &dot);
                        if (status != SUCCESS) {
                            SPRINT3(0, "**> ocsmGetArg(ibrch=%d, iarg=9) -> status=%d (%s)",
                                    ibrch, status, ocsmGetText(status));
                            goto next_option;
                        }
                    } else {
                        strcpy(str9_paste[ipaste], "");
                    }

                    status = ocsmDelBrch(MODL, MODL->nbrch);
                    if (status != SUCCESS) {
                        SPRINT3(0, "**> ocsmDelBrch(ibrch=%d) -> status=%d (%s)",
                                MODL->nbrch, status, ocsmGetText(status));
                        goto next_option;
                    }

                    SPRINT1(0, "Old Branch (%d) deleted", type_paste[ipaste]);
                }

                SPRINT0(0, "Use 'B' to rebuild");

            /* 'y' - undefined option */
            } else if (*state == 'y') {
                SPRINT0(0, "--> Option 'y' (undefined)");

            /* 'Y' - undefined option */
            } else if (*state == 'Y') {
                SPRINT0(0, "--> Option 'Y' (undefined)");

            /* 'z' - undefined option */
            } else if (*state == 'z') {
                SPRINT0(0, "--> Option 'z' (undefined)");

            /* 'Z' - undefined option */
            } else if (*state == 'Z') {
                SPRINT0(0, "--> Option 'Z' (undefined)");

            /* '0' - append "0" to numarg */
            } else if (*state == '0') {
                numarg = 0 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '1' - append "1" to numarg */
            } else if (*state == '1') {
                numarg = 1 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '2' - append "2" to numarg */
            } else if (*state == '2') {
                numarg = 2 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '3' - append "3" to numarg */
            } else if (*state == '3') {
                numarg = 3 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '4' - append "4" to numarg */
            } else if (*state == '4') {
                numarg = 4 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '5' - append "5" to numarg */
            } else if (*state == '5') {
                numarg = 5 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '6' - append "6" to numarg */
            } else if (*state == '6') {
                numarg = 6 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '7' - append "7" to numarg */
            } else if (*state == '7') {
                numarg = 7 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '8' - append "8" to numarg */
            } else if (*state == '8') {
                numarg = 8 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '9' - append "9" to numarg */
            } else if (*state == '9') {
                numarg = 9 + numarg * 10;
                SPRINT1(0, "numarg = %d", numarg);

           /* 'bksp' - erase last digit of numarg */
            } else if (*state == 65288) {
                numarg = numarg / 10;
                SPRINT1(0, "numarg = %d", numarg);

            /* '>' - write viewpoint */
            } else if (*state == '>') {
                sprintf(tempName, "ViewMatrix%d.dat", numarg);
                fp = fopen(tempName, "w");
                fprintf(fp, "%f %f %f %f\n", gv_xform[0][0], gv_xform[1][0],
                                             gv_xform[2][0], gv_xform[3][0]);
                fprintf(fp, "%f %f %f %f\n", gv_xform[0][1], gv_xform[1][1],
                                             gv_xform[2][1], gv_xform[3][1]);
                fprintf(fp, "%f %f %f %f\n", gv_xform[0][2], gv_xform[1][2],
                                             gv_xform[2][2], gv_xform[3][2]);
                fprintf(fp, "%f %f %f %f\n", gv_xform[0][3], gv_xform[1][3],
                                             gv_xform[2][3], gv_xform[3][3]);
                fclose(fp);

                SPRINT1(0, "%s has been saved", tempName);

                numarg = 0;

            /* '<' - read viewpoint */
            } else if (*state == '<') {
                int     count = 0;
                double  size;

                sprintf(tempName, "ViewMatrix%d.dat", numarg);
                fp = fopen(tempName, "r");
                if (fp != NULL) {
                    SPRINT1(0, "resetting to %s", tempName);

                    count += fscanf(fp, "%f%f%f%f", &(gv_xform[0][0]), &(gv_xform[1][0]),
                                                    &(gv_xform[2][0]), &(gv_xform[3][0]));
                    count += fscanf(fp, "%f%f%f%f", &(gv_xform[0][1]), &(gv_xform[1][1]),
                                                    &(gv_xform[2][1]), &(gv_xform[3][1]));
                    count += fscanf(fp, "%f%f%f%f", &(gv_xform[0][2]), &(gv_xform[1][2]),
                                                    &(gv_xform[2][2]), &(gv_xform[3][2]));
                    count += fscanf(fp, "%f%f%f%f", &(gv_xform[0][3]), &(gv_xform[1][3]),
                                                    &(gv_xform[2][3]), &(gv_xform[3][3]));
                    fclose(fp);

                    if (count != 16) {
                        size = 0.5 * sqrt(SQR(bigbox[3]-bigbox[0]) + SQR(bigbox[4]-bigbox[1]) + SQR(bigbox[5]-bigbox[2]));

                        gv_xform[0][0] = +1 / size;
                        gv_xform[1][0] =  0;
                        gv_xform[2][0] =  0;
                        gv_xform[3][0] = -(bigbox[0] + bigbox[3]) / 2 / size;
                        gv_xform[0][1] =  0;
                        gv_xform[1][1] = +1 / size;
                        gv_xform[2][1] =  0;
                        gv_xform[3][1] = -(bigbox[1] + bigbox[4]) / 2 / size;
                        gv_xform[0][2] =  0;
                        gv_xform[1][2] =  0;
                        gv_xform[2][2] = +1 / size;
                        gv_xform[3][2] = -(bigbox[2] + bigbox[5]) / 2 / size;
                        gv_xform[0][3] =  0;
                        gv_xform[1][3] =  0;
                        gv_xform[2][3] =  0;
                        gv_xform[3][3] =  1;
                    }
                } else {
                    SPRINT1(0, "%s does not exist", tempName);
                }

                numarg   = 0;

            /* '$' - read journal file */
            } else if (*state == '$') {
                int    count;

                SPRINT0(0, "--> Option $ chosen (read journal file)");

                if (script == NULL) {
                    SPRINT0(0, "Enter journal filename: ");
                    count = scanf("%s", jnlName);

                    if (count == 1) {
                        SPRINT1x(0, "Opening journal file \"%s\" ...", jnlName);

                        script = fopen(jnlName, "r");
                        if (script != NULL) {
                            SPRINT0(0, "okay");
                        } else {
                            SPRINT0(0, "ERROR detected");
                        }
                    }
                } else {
                    fclose(script);
                    SPRINT0(0, "Closing journal file");

                    script = NULL;
                    *win   =    0;
                }

            /* <home> - original viewpoint */
            } else if (*state == 65360) {
                double  size;

                size = 0.5 * sqrt(SQR(bigbox[3]-bigbox[0]) + SQR(bigbox[4]-bigbox[1]) + SQR(bigbox[5]-bigbox[2]));

                gv_xform[0][0] = +1 / size;
                gv_xform[1][0] =  0;
                gv_xform[2][0] =  0;
                gv_xform[3][0] = -(bigbox[0] + bigbox[3]) / 2 / size;
                gv_xform[0][1] =  0;
                gv_xform[1][1] = +1 / size;
                gv_xform[2][1] =  0;
                gv_xform[3][1] = -(bigbox[1] + bigbox[4]) / 2 / size;
                gv_xform[0][2] =  0;
                gv_xform[1][2] =  0;
                gv_xform[2][2] = +1 / size;
                gv_xform[3][2] = -(bigbox[2] + bigbox[5]) / 2 / size;
                gv_xform[0][3] =  0;
                gv_xform[1][3] =  0;
                gv_xform[2][3] =  0;
                gv_xform[3][3] =  1;

            /* & - toggle flying mode */
            } else if (*state == '&') {
                if (fly_mode == 0) {
                    SPRINT0(0, "--> turning fly mode on");
                    fly_mode = 1;
                } else {
                    SPRINT0(0, "--> turning fly mode off");
                    fly_mode = 0;
                }

            /* <left> - rotate left 30 deg or translate object left */
            } else if (*state == 65361) {
                if (fly_mode == 0) {
                    double  cosrot=cos(PI/6), sinrot=sin(PI/6), temp0, temp2;

                    for (i = 0; i < 4; i++) {
                        temp0 = gv_xform[i][0];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][0] = cosrot * temp0 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp0 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][0] -= 0.5;
                }

            /* <up> - rotate up 30 deg or translate object  up */
            } else if (*state == 65362) {
                if (fly_mode == 0) {
                    double   cosrot=cos(-PI/6), sinrot=sin(-PI/6), temp1, temp2;

                    for (i = 0; i < 4; i++) {
                        temp1 = gv_xform[i][1];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][1] = cosrot * temp1 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp1 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][1] += 0.5;
                }

            /* <rite> - rotate right 30 deg or translate object right */
            } else if (*state == 65363) {
                if (fly_mode == 0) {
                    double   cosrot=cos(-PI/6), sinrot=sin(-PI/6), temp0, temp2;

                    for (i = 0; i < 4; i++) {
                        temp0 = gv_xform[i][0];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][0] = cosrot * temp0 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp0 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][0] += 0.5;
                }

            /* <down> - rotate down 30 deg or translate object down */
            } else if (*state == 65364) {
                if (fly_mode == 0) {
                    double   cosrot=cos(PI/6), sinrot=sin(PI/6), temp1, temp2;

                    for (i = 0; i < 4; i++) {
                        temp1 = gv_xform[i][1];
                        temp2 = gv_xform[i][2];

                        gv_xform[i][1] = cosrot * temp1 - sinrot * temp2;
                        gv_xform[i][2] = sinrot * temp1 + cosrot * temp2;
                    }
                } else {
                    gv_xform[3][1] -= 0.5;
                }

            /* <PgUp> - zoom in */
            } else if (*state == 65365) {
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 3; j++) {
                        gv_xform[i][j] *= 2;
                    }
                }

            /* <PgDn> - zoom out */
            } else if (*state == 65366) {
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 3; j++) {
                        gv_xform[i][j] /= 2;
                    }
                }

            /* '?' - help */
            } else if (*state == '?') {
                SPRINT0(0, "===========================   ===========================   ===========================");
                SPRINT0(0, "                              3D Window - special options                              ");
                SPRINT0(0, "===========================   ===========================   ===========================");
                SPRINT0(0, "L list     Branches           l list Parameters           0-9 build numeric arg (#)    ");
                SPRINT0(0, "E edit     Branch (#)         e edit Parameter           Bksp edit  numeric arg (#)    ");
                SPRINT0(0, "A add      Branch             a add  Parameter                                         ");
                SPRINT0(0, "N name     Branch             d deriv of Parameter (#)  Arrow rot 30 deg or xlate obj  ");
                SPRINT0(0, "T attrib.  Branch             f change tuft length          & toggle lfy mode          ");
                SPRINT0(0, "S suppress Branch (#)         h hide Edge/Face at curs   Home original viewpoint       ");
                SPRINT0(0, "R resume   Branch (#)         U unhide last hidden       PgUp zoom in                  ");
                SPRINT0(0, "D delete   Branch             q query Edge/Face at curs  PgDn zoom out                 ");
                SPRINT0(0, "X cut      Branches (#)       p get para coords             W write .csm file          ");
                SPRINT0(0, "V paste    Branches           P get physical coords         C write .stl file          ");
                SPRINT0(0, "                              M find unmatched Edges        t write .topo file         ");
                SPRINT0(0, "B build to Branch (#)         m view in monochrome          $ read  journal file       ");
                SPRINT0(0, "                              u view u-parameter            < read  viewpoint (#)      ");
                SPRINT0(0, "n print sens at Nodes/Edges   v view v-parameter            > write viewpoint (#)      ");
                SPRINT0(0, "s compare sensitivities       w view surface velcity (#)  ESC exit                     ");
//$$$                SPRINT0(0, "unused options: b* c* C F g G H i I j J k K M* o O Q r x y Y z Z                       ");

            /* 'ESC' - exit program */
            } else if (*state == 65307 || *state == 1) {
                SPRINT0(1, "--> Exiting buildCSM");
                *state = 65307;
                break;
            }

        next_option:
            continue;
        }

        /* repeat as long as we are in a script */
    } while (script != NULL);

//cleanup:
    return;
}


/***********************************************************************/
/*                                                                     */
/*   transform - perform graphic transformation                        */
/*                                                                     */
/***********************************************************************/

void
transform(double    xform[3][4],        /* (in)  transformation matrix */
          double    *point,             /* (in)  input Point */
          float     *out)               /* (out) output Point */
{
    out[0] = xform[0][0]*point[0] + xform[0][1]*point[1]
           + xform[0][2]*point[2] + xform[0][3];
    out[1] = xform[1][0]*point[0] + xform[1][1]*point[1]
           + xform[1][2]*point[2] + xform[1][3];
    out[2] = xform[2][0]*point[0] + xform[2][1]*point[1]
           + xform[2][2]*point[2] + xform[2][3];
}


/***********************************************************************/
/*                                                                     */
/*  pickObject - return the objected pointed to by user                */
/*                                                                     */
/***********************************************************************/

int
pickObject(int       *utype)            /* (out) type of object picked */
{
    int       uindex;                   /* index of object picked */

    int       xpix, ypix, saved_pickmask;
    float     xc, yc;

    (void) GraphicCurrentPointer(&xpix, &ypix);

    xc   =  (2.0 * xpix) / (gv_w3d.xsize - 1.0) - 1.0;
    yc   =  (2.0 * ypix) / (gv_w3d.ysize - 1.0) - 1.0;

    saved_pickmask = gv_pickmask;
    gv_pickmask = -1;
    PickGraphic(xc, -yc, 0);
    gv_pickmask = saved_pickmask;
    if (gv_picked == NULL) {
        *utype = 0;
        uindex = 0;
    } else {
        *utype = gv_picked->utype;
        uindex = gv_picked->uindex;
    }

    return uindex;
}


/***********************************************************************/
/*                                                                     */
/*   getInt - get an int from the user or from a script                */
/*                                                                     */
/***********************************************************************/

static int
getInt(char      *prompt)               /* (in)  prompt string */
{
    int       answer;                   /* integer to return */

    int       count;

    /* integer from screen */
    if (script == NULL) {
        SPRINT1x(0, "%s", prompt); fflush(0);
        count = scanf("%d", &answer);
        if (count != 1) {
            answer = -99999;
        }

    /* integer from script (journal file) */
    } else {
        count = fscanf(script, "%d", &answer);
        if (count != 1) {
            answer = -99999;
        }
        SPRINT2(0, "%s %d", prompt, answer);
    }

    return answer;
}


/***********************************************************************/
/*                                                                     */
/*   getDbl - get a double from the user or from a script              */
/*                                                                     */
/***********************************************************************/

static double
getDbl(char      *prompt)               /* (in)  prompt string */
{
    double    answer;                   /* double to return */

    int       count;

    /* double from screen */
    if (script == NULL) {
        SPRINT1x(0, "%s", prompt); fflush(0);
        count = scanf("%lf", &answer);
        if (count != 1) {
            answer = -99999;
        }

    /* double from script (journal file) */
    } else {
        count = fscanf(script, "%lf", &answer);
        if (count != 1) {
            answer = -99999;
        }
        SPRINT2(0, "%s %f", prompt, answer);
    }

    return answer;
}


/***********************************************************************/
/*                                                                     */
/*   getStr - get a string from the user or from a script              */
/*                                                                     */
/***********************************************************************/

static void
getStr(char      *prompt,               /* (in)  prompt string */
       char      *string)               /* (out) string */
{
    int       count;

    /* double from screen */
    if (script == NULL) {
        SPRINT1x(0, "%s", prompt); fflush(0);
        count = scanf("%254s", string);
        if (count != 1) {
            string[0] = '\0';
        }

    /* double from script (journal file) */
    } else {
        count = fscanf(script, "%s", string);
        if (count != 1) {
            string[0] = '\0';
        }
        SPRINT2(0, "%s %s", prompt, string);
    }
}


/*
 ************************************************************************
 *                                                                      *
 *   evalRuled - evaluate location on a ruled surface                   *
 *                                                                      *
 ************************************************************************
 */

static int
evalRuled(modl_T   *MODL,
          int      ibody,
          int      iface,
          int      isketchs,
          int      isketchn,
          int      iedge,
          double   uv[],
          double   xyz[])
{
    int status = SUCCESS;

    int      periodic;
    double   uvlimits[4], ulimits[4], ulimitn[4], ubar, vbar, t, datas[18], datan[18];

    ROUTINE(evalRuled);

    /* --------------------------------------------------------------- */

    /* default return */
    xyz[0] = 0;
    xyz[1] = 0;
    xyz[2] = 0;

    /* get face and edge limits */
    status = EG_getRange(MODL->body[ibody].face[iface].eface, uvlimits, &periodic);
    CHECK_STATUS(EG_getRange);

    status = EG_getRange(MODL->body[isketchs].edge[iedge].eedge, ulimits, &periodic);
    CHECK_STATUS(EG_getRange);

    status = EG_getRange(MODL->body[isketchn].edge[iedge].eedge, ulimitn, &periodic);
    CHECK_STATUS(EG_getRange);

    /* convert given uv to ubar and vbar */
    ubar = (uv[0] - uvlimits[0]) / (uvlimits[1] - uvlimits[0]);
    vbar = (uv[1] - uvlimits[2]) / (uvlimits[3] - uvlimits[2]);

    /* evaluate coordinates on south and north boundaries */
    t = ulimits[0] + (ulimits[1] - ulimits[0]) * ubar;
    status = EG_evaluate(MODL->body[isketchs].edge[iedge].eedge, &t, datas);
    CHECK_STATUS(EG_evaluate);

    t = ulimitn[0] + (ulimitn[1] - ulimitn[0]) * ubar;
    status = EG_evaluate(MODL->body[isketchn].edge[iedge].eedge, &t, datan);
    CHECK_STATUS(EG_evaluate);

    /* linearly interpolate between south and north */
    xyz[0] = (1-vbar) * datas[0] + vbar * datan[0];
    xyz[1] = (1-vbar) * datas[1] + vbar * datan[1];
    xyz[2] = (1-vbar) * datas[2] + vbar * datan[2];

cleanup:
    return status;
}
